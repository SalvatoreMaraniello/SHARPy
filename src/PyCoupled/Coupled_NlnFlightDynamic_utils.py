'''@package PyCoupled.Coupled_NlnFlightDynamic
@brief      NonlinearDynamic Beam + Rigid-body dynamics + UVLM.
@author     S. Maraniello
@contact    salvatore.maraniello10@imperial.ac.uk
@version    0.0
@date       14/10/2015
@pre        None
@summary    I/O required by Coupled_NlnFlightDynamic are wrapped in this module

            !!! All io routines moved to PiLibs.io.dat: remove after testing !!!

'''
#----------------------------------------------------------------------- Packages
import sys
import Main.SharPySettings as Settings
#import DerivedTypes
import BeamIO
#import BeamLib
#import BeamInit
import numpy as np
#import ctypes as ct
#from PyFSI.Beam2UVLM import InitSection
#from PyFSI.Beam2UVLM import CoincidentGrid
#from PyAero.UVLM.Utils import UVLMLib
#from PyFSI.Beam2UVLM import CoincidentGridForce
#from PyAero.UVLM.Utils import DerivedTypesAero
from PyAero.UVLM.Utils import PostProcess
#from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels
#from PyAero.UVLM.Solver.VLM import InitSteadyWake
#from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
#import PyCoupled.Coupled_NlnStatic as Static
#import PyBeam.Utils.XbeamLib as xbl
#from PyCoupled.Coupled_NlnStatic import AddGravityLoads
#from DerivedTypesAero import ControlSurf
#from collections import OrderedDict
import re
#from math import pow
#from PyBeam.Utils.XbeamLib import Skew
#from PyAero.UVLM.Utils.DerivedTypesAero import Gust
#import h5py, time                                          # sm added packages
#import PyLibs.io.save




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




        
        
    
    