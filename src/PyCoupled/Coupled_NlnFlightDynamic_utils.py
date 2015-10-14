'''@package PyCoupled.Coupled_NlnFlightDynamic
@brief      NonlinearDynamic Beam + Rigid-body dynamics + UVLM.
@author     S. Maraniello
@contact    salvatore.maraniello10@imperial.ac.uk
@version    0.0
@date       14/10/2015
@pre        None
@summary:   I/O required by Coupled_NlnFlightDynamic are wrapped in this module
@warning    MOVE THESE METHODS ELSEWHERE


'''
#----------------------------------------------------------------------- Packages
import sys
import Main.SharPySettings as Settings
import DerivedTypes
import BeamIO
#import BeamLib
#import BeamInit
import numpy as np
import ctypes as ct
#from PyFSI.Beam2UVLM import InitSection
#from PyFSI.Beam2UVLM import CoincidentGrid
#from PyAero.UVLM.Utils import UVLMLib
#from PyFSI.Beam2UVLM import CoincidentGridForce
#from PyAero.UVLM.Utils import DerivedTypesAero
from PyAero.UVLM.Utils import PostProcess
#from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels
#from PyAero.UVLM.Solver.VLM import InitSteadyWake
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
#import PyCoupled.Coupled_NlnStatic as Static
#import PyBeam.Utils.XbeamLib as xbl
#from PyCoupled.Coupled_NlnStatic import AddGravityLoads
from DerivedTypesAero import ControlSurf
from collections import OrderedDict
import re
from math import pow
from PyBeam.Utils.XbeamLib import Skew
from PyAero.UVLM.Utils.DerivedTypesAero import Gust
import h5py, time                                          # sm added packages
import PyLibs.io.save


        

def panellingFromFreq(freq,c=1.0,Umag=1.0):
    """@brief Calculate adequate spatial/temporal resolution on UVLM grid
    based on a frequency of interest.
    @param freq Frequency of interest.
    @param c chord of wing.
    @param Umag mean reltive free-stream velocity magnitude.
    @returns M Number of chordwise panels for wing.
    @returns DelTime Suggested timestep [seconds.]
    """
    k = freq*c/(2*Umag) #get expected reduced freq
    M = int(50.0*k/np.pi) #get adequate panelling
    DelTime = c/(Umag*M) #get resulting DelTime
    return M, DelTime


def saveh5(SaveDict, AELAOPTS, VMINPUT, VMOPTS, XBOPTS, XBINPUT, XBOUT ):
    
    h5filename=SaveDict['OutputDir'] + SaveDict['OutputFileRoot'] + '.h5'
    hdfile=h5py.File(h5filename,'w')
    #print ('created %s'%h5filename)
    PyLibs.io.save.add_class_as_grp(AELAOPTS,hdfile)
    PyLibs.io.save.add_class_as_grp(VMINPUT,hdfile)
    PyLibs.io.save.add_class_as_grp(VMOPTS,hdfile)
    PyLibs.io.save.add_class_as_grp(XBOPTS,hdfile)
    PyLibs.io.save.add_class_as_grp(XBINPUT,hdfile)
    PyLibs.io.save.add_class_as_grp(XBOUT,hdfile)
    hdfile.close()
    return None


def write_SOL912_def(XBOPTS,XBINPUT,XBELEM,NumNodes_tot,PosDefor,PsiDefor,SaveDict):
        
    # Write deformed configuration to file. TODO: tidy this away inside function.
    ofile = SaveDict['OutputDir'] + SaveDict['OutputFileRoot'] + '_SOL912_def.dat'
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Writing file %s ... ' %(ofile))
    fp = open(ofile,'w')
    fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
    fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
    fp.close()
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('done\n')
    WriteMode = 'a'
    # Write
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM,
                       PosDefor, PsiDefor, ofile, WriteMode)
    return None
    

def write_SOL912_final(Time, PosPsiTime, NumNodes_tot, DynOut, Vrel, VrelDot, SaveDict):

    # Write _dyn file.
    ofile = SaveDict['OutputDir'] + SaveDict['OutputFileRoot'] + '_SOL912_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
    #Write _shape file
    ofile = SaveDict['OutputDir'] + SaveDict['OutputFileRoot'] + '_SOL912_shape.dat'
    fp = open(ofile,'w')
    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
    fp.close()
    
    #Write rigid-body velocities
    ofile = SaveDict['OutputDir'] + SaveDict['OutputFileRoot'] + '_SOL912_rigid.dat'
    fp = open(ofile,'w')
    BeamIO.Write_rigid_File(fp, Time, Vrel, VrelDot)
    fp.close()
    
    return None


def write_SOL912_out(time, PosDefor, PsiDefor, PosIni, PsiIni, XBELEM, 
                     writeDict, SaveDict , **kwargs):

    '''
    Writes structural output. Function called at each time-step to update the 
    solution. 
    If a file object is not passed in input, this is created from scratch. 
    '''
    
    if 'FileObject' in kwargs:
        fp = kwargs['FileObject']
    else:
        # create file
        ofile = SaveDict['OutputDir'] +SaveDict['OutputFileRoot'] + '_SOL912_out.dat'
        fp = open(ofile,'w')
        fp.write("{:<14}".format("Time"))
        for output in writeDict.keys():
            fp.write("{:<14}".format(output))
        fp.write("\n")
        fp.flush()
     
        
    # Write initial outputs to file.
    locForces = None # Stops recalculation of forces
    fp.write("{:<14,e}".format(time))
    for myStr in writeDict.keys():
        
        # Search for Forces
        if re.search(r'^R_.',myStr):
            
            if re.search(r'^R_._.', myStr): index = int(myStr[4])
            elif re.search(r'root', myStr): index = 0
            elif re.search(r'tip', myStr): index = -1
            else:
                raise IOError("Node index not recognised.")
            
            if myStr[2] == 'x':  component = 0
            elif myStr[2] == 'y': component = 1
            elif myStr[2] == 'z': component = 2
            else: raise IOError("Displacement component not recognised.")
            
            fp.write("{:<14,e}".format(PosDefor[index,component]))
        
        # Search for Moments    
        elif re.search(r'^M_.',myStr):
            
            if re.search(r'^M_._.', myStr): index = int(myStr[4])
            elif re.search(r'root', myStr): index =  0
            elif re.search(r'tip', myStr): index = -1
            else: raise IOError("Node index not recognised.")
            
            if   myStr[2] == 'x': component = 0
            elif myStr[2] == 'y': component = 1
            elif myStr[2] == 'z': component = 2
            else: raise IOError("Moment component not recognised.")
            
            if locForces == None:
                locForces = BeamIO.localElasticForces(PosDefor, PsiDefor,
                                                      PosIni, PsiIni,
                                                      XBELEM,
                                                      [index])
            
            fp.write("{:<14,e}".format(locForces[0,3+component]))
        else:
            raise IOError("writeDict key not recognised.")
    # END for myStr
    fp.write("\n")
    fp.flush()
 
    return fp
    

def write_TecPlot(Zeta, ZetaStar, Gamma, GammaStar, NumTimeSteps, iStep, time, SaveDict, **kwargs):
    '''
    Writes TecPlot file. Function called during the time-stepping to update the
    solution. 
    If a file object is not passed in input, this is created from scratch. 
    '''
    
    if 'FileObject' in kwargs: 
        FileObject=kwargs['FileObject']
    else:
        FileName = SaveDict['OutputDir']+SaveDict['OutputFileRoot']+'AeroGrid.dat'
        Variables = ['X', 'Y', 'Z','Gamma']    
        FileObject = PostProcess.WriteAeroTecHeader(FileName,'Default',Variables)
    
    # Plot results of static analysis
    PostProcess.WriteUVLMtoTec(    FileObject,
                                   Zeta, ZetaStar,
                                   Gamma, GammaStar,
                                   TimeStep = iStep,
                                   NumTimeSteps = NumTimeSteps,#XBOPTS.NumLoadSteps.value,
                                   Time = time,
                                   Text = True)
    
    return FileObject



if __name__ == '__main__':
    # Beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 OutInaframe = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(1),
                                 Solution = ct.c_int(912),
                                 MinDelta = ct.c_double(1e-5),
                                 NewmarkDamp = ct.c_double(5e-3))
    # beam inputs.
    XBINPUT = DerivedTypes.Xbinput(2,4)
    XBINPUT.BeamLength = 6.096
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 0.99e+06
    XBINPUT.BeamStiffness[4,4] = 9.77e+06
    XBINPUT.BeamStiffness[5,5] = 1.0e+09
    XBINPUT.BeamStiffness[:,:] = 1.0e0*XBINPUT.BeamStiffness[:,:]
    XBINPUT.BeamMass[0,0] = 35.71
    XBINPUT.BeamMass[1,1] = 35.71
    XBINPUT.BeamMass[2,2] = 35.71
    XBINPUT.BeamMass[3,3] = 8.64*10
    XBINPUT.BeamMass[4,4] = 1.0e+9
    XBINPUT.BeamMass[5,5] = 1.0e+9
    # Off diagonal terms (in Theodorsen sectional coordinates)
    ElasticAxis = -0.34
    InertialAxis = -7.0/50.0
    x_alpha = -(InertialAxis - ElasticAxis)
    # pitch-plunge coupling term (b-frame coordinates)
    c = 1.8288
    cgLoc = 0.5*c*np.array([0.0, x_alpha, 0.0])
    cgSkew = Skew(cgLoc)
#     mOff = x_alpha*(c/2)*XBINPUT.BeamMass[0,0]
#     XBINPUT.BeamMass[2,3] = -mOff
#     XBINPUT.BeamMass[0,5] = mOff
#     XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    XBINPUT.BeamMass[:3,3:] = XBINPUT.BeamMass[0,0] * cgSkew
    XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    XBINPUT.BeamMass[4,4] += XBINPUT.BeamMass[0,0]*pow(cgLoc[2],2.0)
    XBINPUT.BeamMass[5,5] += XBINPUT.BeamMass[0,0]*pow(cgLoc[1],2.0)
    
    # Get suggested panelling.
    Umag = 50.0
    AOA  = 5.0*np.pi/180.0
    M = 4
    delTime = c/(Umag*M)
    # Unsteady parameters.
    XBINPUT.dt = delTime
    XBINPUT.t0 = 0.0
    XBINPUT.tfin = 2
    
    # Set motion of wing.
    NumSteps = np.ceil( (XBINPUT.tfin + XBINPUT.dt - XBINPUT.t0) / XBINPUT.dt)
    XBINPUT.ForcedVel = np.zeros((NumSteps,6),ct.c_double,'F')
    for i in range(XBINPUT.ForcedVel.shape[0]):
        XBINPUT.ForcedVel[i,:] = [0.0, Umag*np.cos(AOA), -Umag*np.sin(AOA), 0.0, 0.0, 0.0]
#         XBINPUT.ForcedVel[i,:] = [0.0, Umag, 0.0, 0.0, 0.0, 0.0]
    XBINPUT.ForcedVelDot = np.zeros((NumSteps,6),ct.c_double,'F')
     
    # aero params.
    WakeLength = 30.0*c
    Mstar = int(WakeLength/(delTime*Umag))
    # aero options.
    N = XBINPUT.NumNodesTot - 1
    VMOPTS = DerivedTypesAero.VMopts(M = M,
                                     N = N,
                                     ImageMethod = False,
                                     Mstar = Mstar,
                                     Steady = False,
                                     KJMeth = False,
                                     NewAIC = True,
                                     DelTime = delTime,
                                     NumCores = 4)
    # Aero inputs.
    iMin = M - M/4
    iMax = M
    jMin = N - N/4
    jMax = N
    typeMotion = 'sin'
    betaBar = 0.0*np.pi/180.0
    omega = 30.0
    ctrlSurf = ControlSurf(iMin,
                           iMax,
                           jMin,
                           jMax,
                           typeMotion,
                           betaBar,
                           omega)
    
    # Gust inputs
    gust = Gust(uMag = 0.0*Umag,
                h = 10.0,
                r = 0.0)
    
    VMINPUT = DerivedTypesAero.VMinput(c = c,
                                       b = XBINPUT.BeamLength,
                                       U_mag = 0.0,
                                       alpha = 0.0*np.pi/180.0,
                                       theta = 0.0,
                                       WakeLength = WakeLength,
                                       ctrlSurf = ctrlSurf,
                                       gust = None)
    
    # Aerolastic simulation results.
    AELAOPTS = AeroelasticOps(ElasticAxis = ElasticAxis,
                              InertialAxis = InertialAxis,
                              AirDensity = 1.02,
                              Tight = False,
                              ImpStart = False)
    
    # Live output options.
    writeDict = OrderedDict()
    writeDict['R_z (tip)'] = 0
    writeDict['M_x (root)'] = 0
    writeDict['M_y (root)'] = 0
    writeDict['M_z (root)'] = 0
    
    SaveDict=Settings.SaveDict
    SaveDict['OutputDir'] = Settings.SharPyProjectDir + "output/temp/"
    SaveDict['OutputFileRoot'] = "FlyingWing"
    
    # Solve nonlinear dynamic simulation.
    Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS,writeDict =  writeDict,SaveDict=SaveDict)