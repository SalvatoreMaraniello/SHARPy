'''@package PyCoupled.Coupled_NlnFlightDynamic
@brief      NonlinearDynamic Beam + Rigid-body dynamics + UVLM.
@author     Rob Simpson & Henrik Hesse
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       07/11/2013
@pre        None
@warning    None

Modified S. Maraniello, 25 Sep 2015

'''
#----------------------------------------------------------------------- Packages
import sys, os, h5py
import Main.SharPySettings as Settings
import DerivedTypes
#import BeamIO
import BeamLib
import BeamInit
import numpy as np
import ctypes as ct
from PyFSI.Beam2UVLM import InitSection
from PyFSI.Beam2UVLM import CoincidentGrid
from PyFSI.Beam2UVLM import CoincidentGridForce
from PyAero.UVLM.Utils import UVLMLib
#from PyAero.UVLM.Utils import DerivedTypesAero
from PyAero.UVLM.Utils import PostProcess
from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels
from PyAero.UVLM.Solver.VLM import InitSteadyWake
#from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
import PyCoupled.Coupled_NlnStatic as Static
import PyBeam.Utils.XbeamLib as xbl
from PyCoupled.Coupled_NlnStatic import AddGravityLoads
#from DerivedTypesAero import ControlSurf
#from collections import OrderedDict
#import re
#from math import pow
#from PyBeam.Utils.XbeamLib import Skew
#from PyAero.UVLM.Utils.DerivedTypesAero import Gust
import time                                          # sm added packages
import PyLibs.io.save

from PyCoupled.Coupled_NlnFlightDynamic_utils import write_SOL912_def, write_SOL912_final, \
                                                     write_SOL912_out, write_TecPlot



def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS,**kwords):
    """@brief Nonlinear dynamic solver using Python to solve aeroelastic
    equation.
    @details Assembly of structural matrices is carried out with 
    Fortran subroutines. Aerodynamics solved using PyAero\.UVLM.
    @param XBINPUT Beam inputs (for initialization in Python).
    @param XBOPTS Beam solver options (for Fortran).
    @param VMOPTS UVLM solver options (for C/C++).
    @param VMINPUT UVLM solver inputs (for initialization in Python).
    @param VMUNST Unsteady input information for aero solver.
    @param AELAOPTS Options relevant to coupled aeroelastic simulations.
    @param writeDict OrderedDict of 'name':tuple of outputs to write.
    """
        
    # Check correct solution code.
    assert XBOPTS.Solution.value == 912, ('NonlinearFlightDynamic requested' +
                                          ' with wrong solution code')
    
    # I/O management
    XBOUT=DerivedTypes.Xboutput()  
    SaveDict=Settings.SaveDict
    if 'SaveDict' in kwords: SaveDict=kwords['SaveDict']
    if SaveDict['Format']=='h5':
        Settings.WriteOut=False
        Settings.PlotTec=False
        OutList=[AELAOPTS, VMINPUT, VMOPTS, XBOPTS, XBINPUT, XBOUT]
        #if VMINPUT.ctrlSurf!=None:
        #    for cc in range(len(VMINPUT.ctrlSurf)):
        #        OutList.append(VMINPUT.ctrlSurf[cc])
        if SaveDict['SaveWake'] is True:
            dirwake=SaveDict['OutputDir']+'wake'+SaveDict['OutputFileRoot']+'/'
            os.system('mkdir -p %s' %dirwake)
            XBOUT.dirwake=dirwake  
            
    XBOUT.cputime.append(time.clock()) # time.processor_time more appropriate but equivalent
    # for debugging
    XBOUT.ForceDofList=[]
    XBOUT.ForceRigidList=[]

    # Initialise static beam data.
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    # special BCs
    SphFlag=False
    if XBNODE.Sflag.any(): SphFlag=True
    
    # Debugging Flags

    if 'HardCodeAero' in kwords: HardCodeAero=kwords['HardCodeAero']
    SaveExtraVariables = False  

    
    #------------------------- Initial Displacement: ImpStart vs Static Solution
    
    # Calculate initial displacements.
    if AELAOPTS.ImpStart == False:
        XBOPTS.Solution.value = 112 # Modify options.
        VMOPTS.Steady = ct.c_bool(True)
        Rollup = VMOPTS.Rollup.value
        VMOPTS.Rollup.value = False
        # Solve Static Aeroelastic.
        PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, Force = \
                    Static.Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS)
        XBOPTS.Solution.value = 912 # Reset options.
        VMOPTS.Steady = ct.c_bool(False)
        VMOPTS.Rollup.value = Rollup
    elif AELAOPTS.ImpStart == True:
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        Force = np.zeros((XBINPUT.NumNodesTot,6),ct.c_double,'F')


    if SaveDict['Format']=='dat': 
        write_SOL912_def(XBOPTS,XBINPUT,XBELEM,NumNodes_tot,PosDefor,PsiDefor,SaveDict)
    
    
    
    #------------------------------------------------ Initialise Dynamic Problem
    
    # Initialise structural variables for dynamic analysis.
    Time, NumSteps, ForceTime, Vrel, VrelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    # Delete unused variables.
    del OutGrids, VelocTime
    
    # sm I/O
    ### why forced velocity with Sol912 ???
    ### If forced velocities are prescribed, then is Sol312        
    XBOUT.PosDefor=PosDefor                # ...SOL912_def.dat
    XBOUT.PsiDefor=PsiDefor   
    XBOUT.ForceTime_force=ForceTime        # ...SOL912_force.dat
    XBOUT.Vrel_force=Vrel
    XBOUT.VrelDot_force=VrelDot   
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    
    
    #------------------------------------------------------ Initialise Variables
    #Initialise structural system tensors
    MssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    CssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    KssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    FstrucFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    
    ms = ct.c_int()
    cs = ct.c_int()
    ks = ct.c_int()
    fs = ct.c_int()
    
    Msr = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Csr = np.zeros((NumDof.value,6), ct.c_double, 'F')
    
    X     = np.zeros(NumDof.value, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    Qstruc = np.zeros(NumDof.value, ct.c_double, 'F')
    
    #Initialise rigid-body system tensors
    MrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F')
    CrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F') 
    KrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F') 
    FrigidFull = np.zeros((6,NumDof.value+6), ct.c_double, 'F')
    
    mr = ct.c_int()
    cr = ct.c_int()
    kr = ct.c_int()
    fr = ct.c_int()
    
    Mrr = np.zeros((6,6), ct.c_double, 'F')
    Crr = np.zeros((6,6), ct.c_double, 'F')
        
    Qrigid = np.zeros(6, ct.c_double, 'F')
    
    #Initialise full system tensors
    Q     = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    DQ    = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    dQdt  = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    dQddt = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    Force_All = np.zeros(NumDof.value+6, ct.c_double, 'F')

    Msys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F')
    Csys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F') 
    Ksys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F') 
    Asys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F')
    
    Qsys = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    

    
    #---------------------------------------------------- Start Dynamic Solution
    
    #Initialise rotation operators. TODO: include initial AOA here
    currVrel=Vrel[0,:].copy('F')
    
    
    # Initialise attitude:
    Quat =  xbl.psi2quat(XBINPUT.PsiA_G)
    #Quat=   XBINPUT.quat0
    
    #### sm debug
    XBOUT.Quat0=Quat
    XBOUT.currVel0=currVrel
    
    Cao  = xbl.Rot(Quat)
    ACoa = np.zeros((6,6), ct.c_double, 'F')
    ACoa[:3,:3] = np.transpose(Cao)
    ACoa[3:,3:] = np.transpose(Cao)
    Cqr = np.zeros((4,6), ct.c_double, 'F')
    Cqq = np.zeros((4,4), ct.c_double, 'F')
        
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4):
        Unit4[i,i] = 1.0
    
    # Extract initial displacements and velocities.
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,
                                  PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                                  X, dXdt)
    
    # Approximate initial accelerations.
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                             ct.c_double, 'F')
    
    #Populate state vector
    Q[:NumDof.value]=X.copy('F')
    dQdt[:NumDof.value]=dXdt.copy('F')
    dQdt[NumDof.value:NumDof.value+6] = Vrel[0,:].copy('F')
    dQdt[NumDof.value+6:]= Quat.copy('F')
    
    #Force at the first time-step
    Force += (XBINPUT.ForceDyn*ForceTime[0]).copy('F')
    

    #Assemble matrices and loads for structural dynamic analysis
    currVrel=Vrel[0,:].copy('F')
    tmpQuat=Quat.copy('F')
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                         PosIni, PsiIni, PosDefor, PsiDefor,
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                         Force, currVrel, 0*currVrel,
                         NumDof, Settings.DimMat,
                         ms, MssFull, Msr,
                         cs, CssFull, Csr,
                         ks, KssFull, fs, FstrucFull,
                         Qstruc, XBOPTS, Cao)
       
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),
                              ct.byref(ct.c_int(6)),
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                              ct.byref(NumDof),
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )

    Qstruc -= np.dot(FstrucFull, Force_Dof)
    
    
    #Assemble matrices for rigid-body dynamic analysis
    BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                         PosIni, PsiIni, PosDefor, PsiDefor,
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                         currVrel, 0*currVrel, tmpQuat,
                         NumDof, Settings.DimMat,
                         mr, MrsFull, Mrr,
                         cr, CrsFull, Crr, Cqr, Cqq,
                         kr, KrsFull, fr, FrigidFull,
                         Qrigid, XBOPTS, Cao)
    
    BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),
                               ct.byref(ct.c_int(6)),
                               Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                               ct.byref(ct.c_int(NumDof.value+6)),
                               Force_All.ctypes.data_as(ct.POINTER(ct.c_double)) )

    Qrigid -= np.dot(FrigidFull, Force_All)
         
#     #Separate assembly of follower and dead loads   
#     tmpForceTime=ForceTime[0].copy('F') 
#     tmpQforces,Dummy,tmpQrigid = xbl.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, 
#                                  PosIni, PsiIni, PosDefor, PsiDefor, 
#                                  (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), 
#                                  (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), 
#                                   Cao,1)
#                            
#     Qstruc -= tmpQforces      
#     Qrigid -= tmpQrigid
    
    #Assemble system matrices
    Msys[:NumDof.value,:NumDof.value] = MssFull.copy('F')
    Msys[:NumDof.value,NumDof.value:NumDof.value+6] = Msr.copy('F')
    Msys[NumDof.value:NumDof.value+6,:NumDof.value] = MrsFull.copy('F')
    Msys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Mrr.copy('F')
    Msys[NumDof.value+6:,NumDof.value+6:] = Unit4.copy('F')
       
    Qsys[:NumDof.value] = Qstruc
    Qsys[NumDof.value:NumDof.value+6] = Qrigid
    Qsys[NumDof.value+6:] = np.dot(Cqq,dQdt[NumDof.value+6:])
     
      
    # -------------------------------------------------------------------   
    # special BCs
    iiblock=[]
    
    # block translations
    if SphFlag:
        iiblock = [ ii for ii in range(NumDof.value,NumDof.value+3) ]
    
    # block translations (redundant:)
    for ii in range(3):
        if XBINPUT.EnforceTraVel_FoRA[ii] is True:
            iiblock.append(NumDof.value+ii)
    
    # block rotations
    iirotfree=[] # free rotational dof 
    for ii in range(3):
        if XBINPUT.EnforceAngVel_FoRA[ii] is True:
            iiblock.append(NumDof.value+3+ii)
        else:
            iirotfree.append(NumDof.value+3+ii)
    
    # modify matrices
    if len(iiblock)>0:
        # block dof
        Msys[iiblock,:] = 0.0
        Msys[iiblock,iiblock] = 1.0
        Qsys[iiblock] = 0.0
    # ------------------------------------------------------------------
    
        # add damp at the spherical joints
        if XBINPUT.sph_joint_damping is not None:
            Qsys[iirotfree]+= XBINPUT.sph_joint_damping*dQdt[iirotfree]
    
            
    # add structural damping term
    if XBINPUT.str_damping_model is not None:
        Cdamp = XBINPUT.str_damping_param['alpha'] * MssFull + \
                XBINPUT.str_damping_param['beta']  * KssFull
        Qsys[:NumDof.value] += np.dot( Cdamp, dQdt[:NumDof.value] )                  
        pass
        
    #store initial matrices for eigenvalues analysis
    XBOUT.MssFull0 = MssFull.copy()
    XBOUT.CssFull0 = CssFull.copy()
    XBOUT.KssFull0 = KssFull.copy()

    # Initial Accel.
    ###dQddt[:] = np.dot(np.linalg.inv(Msys), -Qsys)
    dQddt[:] = np.linalg.solve(Msys,-Qsys)
    
    XBOUT.dQddt0=dQddt.copy()
    
    #Record position of all grid points in global FoR at initial time step
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    #Position/rotation of the selected node in initial deformed configuration
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    
    #Get gamma and beta for Newmark scheme
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*(gamma + 0.5)**2
    
    
    #---------------------------------------------- Initialise Aerodynamic Force
    # Initialise Aero       
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
    
    # Declare memory for Aero variables.
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    K = VMOPTS.M.value*VMOPTS.N.value
    AIC = np.zeros((K,K),ct.c_double,'C')
    BIC = np.zeros((K,K),ct.c_double,'C')
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Initialise A-frame location and orientation to be zero.
    OriginA_a = np.zeros(3,ct.c_double,'C')
    PsiA_G = XBINPUT.PsiA_G.copy() #xbl.quat2psi(Quat) # CRV at iStep    
    
    # Init external velocities.  
    Ufree = InitSteadyExternalVels(VMOPTS,VMINPUT)
    if AELAOPTS.ImpStart == True:
        Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')             
        Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
        # Generate surface, wake and gamma matrices.
        CoincidentGrid(PosDefor, PsiDefor, Section, currVrel[:3], 
                       currVrel[3:], PosDotDef, PsiDotDef, XBINPUT,
                       Zeta, ZetaDot, OriginA_a, PsiA_G,
                       VMINPUT.ctrlSurf)
        # init wake grid and gamma matrix.
        ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta,currVrel[:3])
        
    # sm save
    #XBOUT.Zeta0 = Zeta.copy('C')
    #XBOUT.ZetaStar0 = ZetaStar.copy('C')
    #XBOUT.ZetaStarList.append(np.float32( ZetaStar.copy('C') ))
    
    # Define TecPlot stuff
    if Settings.PlotTec==True:
        FileObject=write_TecPlot(Zeta, ZetaStar, Gamma, GammaStar, NumSteps.value, 0, Time[0], SaveDict)
    
    if 'writeDict' in kwords and Settings.WriteOut == True:
        fp= write_SOL912_out(Time[0], PosDefor, PsiDefor, PosIni, PsiIni, XBELEM, kwords['writeDict'], SaveDict)
            
    
    # sm write class
    XBOUT.QuatList.append(Quat.copy())
    XBOUT.CRVList.append(PsiA_G)
    XBOUT.PosIni=PosIni
    XBOUT.PsiIni=PsiIni    
    
    XBOUT.AsysListStart=[]
    XBOUT.AsysListEnd=[]
    XBOUT.MsysList=[]
    XBOUT.CsysList=[]
    XBOUT.KsysList=[]
    

    #---------------------------------------------------------------- Time loop
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        #calculate dt
        dt = Time[iStep+1] - Time[iStep]
        
        # Set dt for aero force calcs.
        VMOPTS.DelTime = ct.c_double(dt)
        
        #Predictor step
        Q       += dt*dQdt + (0.5-beta)*dQddt*np.power(dt,2.0)
        dQdt    += (1.0-gamma)*dQddt*dt
        dQddt[:] = 0.0
        
        # Quaternion update for orientation.
        Quat = dQdt[NumDof.value+6:].copy('F')
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = xbl.Rot(Quat)
        
        #nodal displacements and velocities from state vector
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F'); 
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT,NumNodes_tot,XBELEM,XBNODE,
                                       PosIni,PsiIni,NumDof,X,dXdt,
                                       PosDefor,PsiDefor,PosDotDef,PsiDotDef)
        

                      
          
        #Reset convergence parameters
        Iter = 0
        ResLog10 = 1.0
        
        
        #Newton-Raphson loop      
        while ( (ResLog10 > XBOPTS.MinDelta.value) & (Iter < XBOPTS.MaxIterations.value) ):
             
                        
            #------------------------------------------------------- Aero Force Loop
            # Force at current time-step. TODO: Check communication flow. 
            if (iStep > 0 and AELAOPTS.Tight == False)  and (ResLog10>1.0 or Iter==0):
                # - Force not computed at first time-step 
                # - If iStep>0: 
                #    - the aerodynamic force is always computed at least once (Iter==0)
                #    - the aerodynamic force is update until the residual does not fall 
                #      below a prescribed factor
                
                # zero aero forces.
                AeroForces[:,:,:] = 0.0
                
                # Update CRV.
                PsiA_G = xbl.quat2psi(Quat) # CRV at iStep
                
                # Update origin.
                # sm: origin position projected in FoR A
                currVrel=Vrel[iStep-1,:].copy('F')
                OriginA_a[:] = OriginA_a[:] + currVrel[:3]*dt #sm: OriginA_a initialised to zero
                
                # Update control surface deflection.
                if VMINPUT.ctrlSurf != None:
                    # open-loop control
                    for cc in range(len(VMINPUT.ctrlSurf)):
                        VMINPUT.ctrlSurf[cc].update(Time[iStep],iStep=iStep)
                
                # Generate surface grid.
                currVrel=Vrel[iStep,:].copy('F')
                CoincidentGrid(PosDefor, PsiDefor, Section, currVrel[:3], 
                               currVrel[3:], PosDotDef, PsiDotDef, XBINPUT,
                               Zeta, ZetaDot, OriginA_a, PsiA_G,
                               VMINPUT.ctrlSurf)
                
                # Update wake geom       
                #'roll' data.
                ZetaStar = np.roll(ZetaStar,1,axis = 0)
                GammaStar = np.roll(GammaStar,1,axis = 0)
                #overwrite grid points with TE.
                ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
                # overwrite Gamma with TE value from previous timestep.
                GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
                
                # Apply gust velocity.
                if VMINPUT.gust != None:
                    Utot = Ufree + VMINPUT.gust.Vels(Zeta)
                else:
                    Utot = Ufree
                
                # Solve for AeroForces
                UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Utot, ZetaStar, VMOPTS, 
                               AeroForces, Gamma, GammaStar, AIC, BIC)
                
                # Apply density scaling
                AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
                
                if Settings.PlotTec==True:
                    FileObject=write_TecPlot(Zeta, ZetaStar, Gamma, GammaStar, NumSteps.value, iStep, 
                                             Time[iStep], SaveDict,FileObject=FileObject)
    
                # map AeroForces to beam.
                CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                                    Force, PsiA_G)
    
                ForceAero = Force.copy('C')  
                
                # Add gravity loads.
                AddGravityLoads(Force, XBINPUT, XBELEM, AELAOPTS,
                                PsiDefor, VMINPUT.c)
                
                # Add thrust and other point loads
                Force += (XBINPUT.ForceStatic + 
                          XBINPUT.ForceDyn*ForceTime[iStep+1]).copy('F')
                
                
                #END if iStep > 0

                                    
            #set tensors to zero 
            MssFull[:,:] = 0.0; CssFull[:,:] = 0.0
            KssFull[:,:] = 0.0; FstrucFull[:,:] = 0.0
            Msr[:,:] = 0.0; Csr[:,:] = 0.0
            Qstruc[:] = 0.0
            
            MrsFull[:,:] = 0.0; CrsFull[:,:] = 0.0
            KrsFull[:,:] = 0.0; FrigidFull[:,:] = 0.0
            Mrr[:,:] = 0.0; Crr[:,:] = 0.0
            Qrigid[:] = 0.0
    
            Msys[:,:] = 0.0; Csys[:,:] = 0.0
            Ksys[:,:] = 0.0; Asys[:,:] = 0.0
            Qsys[:] = 0.0
            
            # Update counter.
            Iter += 1
            if XBOPTS.PrintInfo.value==True: sys.stdout.write('   %-7d ' %(Iter))
                        
            # Nodal displacements and velocities from state vector
            X=Q[:NumDof.value].copy('F') 
            dXdt=dQdt[:NumDof.value].copy('F'); 
            BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot,
                                           XBELEM, XBNODE,
                                           PosIni, PsiIni,
                                           NumDof, X, dXdt,
                                           PosDefor, PsiDefor,
                                           PosDotDef, PsiDotDef)


            #rigid-body velocities and orientation from state vector
            Vrel[iStep+1,:]    = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep+1,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
            Quat = dQdt[NumDof.value+6:].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            Cao  = xbl.Rot(Quat)


            #Update matrices and loads for structural dynamic analysis
            tmpVrel=Vrel[iStep+1,:].copy('F')
            tmpQuat=Quat.copy('F')
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                                 PosIni, PsiIni, PosDefor, PsiDefor,
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                                 Force, tmpVrel, 0*tmpVrel,
                                 NumDof, Settings.DimMat,
                                 ms, MssFull, Msr,
                                 cs, CssFull, Csr,
                                 ks, KssFull, fs, FstrucFull,
                                 Qstruc, XBOPTS, Cao)
            
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),
                              ct.byref(ct.c_int(6)),
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                              ct.byref(NumDof),
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
                    
            
            #Update matrices for rigid-body dynamic analysis
            BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                                 PosIni, PsiIni, PosDefor, PsiDefor,
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                                 tmpVrel, 0*tmpVrel, tmpQuat,
                                 NumDof, Settings.DimMat,
                                 mr, MrsFull, Mrr,
                                 cr, CrsFull, Crr, Cqr, Cqq,
                                 kr, KrsFull, fs, FrigidFull,
                                 Qrigid, XBOPTS, Cao)
    
            BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),
                                       ct.byref(ct.c_int(6)),
                                       Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                                       ct.byref(ct.c_int(NumDof.value+6)),
                                       Force_All.ctypes.data_as(ct.POINTER(ct.c_double)) )
        
        
            #Residual at first iteration
            if(Iter == 1):
                Res0_Qglobal = max(max(abs(Qsys)),1)
                Res0_DeltaX  = max(max(abs(DQ)),1)
              
            
            #Assemble discrete system matrices with linearised quaternion equations          
            Msys[:NumDof.value,:NumDof.value] = MssFull.copy('F')
            Msys[:NumDof.value,NumDof.value:NumDof.value+6] = Msr.copy('F')
            Msys[NumDof.value:NumDof.value+6,:NumDof.value] = MrsFull.copy('F')
            Msys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Mrr.copy('F')
            Msys[NumDof.value+6:,NumDof.value+6:] = Unit4.copy('F')
            
            Csys[:NumDof.value,:NumDof.value] = CssFull.copy('F')
            Csys[:NumDof.value,NumDof.value:NumDof.value+6] = Csr.copy('F')
            Csys[NumDof.value:NumDof.value+6,:NumDof.value] = CrsFull.copy('F')
            Csys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Crr.copy('F')
            
            Csys[NumDof.value+6:,NumDof.value:NumDof.value+6] = Cqr.copy('F')
            Csys[NumDof.value+6:,NumDof.value+6:] = Cqq.copy('F')
            
            Ksys[:NumDof.value,:NumDof.value] = KssFull.copy('F')
            Ksys[NumDof.value:NumDof.value+6,:NumDof.value] = KrsFull.copy('F')
                     
#             #Separate assembly of follower and dead loads   
#             tmpForceTime=ForceTime[iStep+1].copy('F') 
#             tmpQforces,Dummy,tmpQrigid = xbl.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
#                                             PosIni, PsiIni, PosDefor, PsiDefor, \
#                                             (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
#                                             (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
#                                             Cao,1)
#                                    
#             Qstruc -= tmpQforces      
#             Qrigid -= tmpQrigid

            #Compute residual to solve update vector
            Qstruc += -np.dot(FstrucFull, Force_Dof)
            Qrigid += -np.dot(FrigidFull, Force_All)
            
            Qsys[:NumDof.value] = Qstruc
            Qsys[NumDof.value:NumDof.value+6] = Qrigid
            Qsys[NumDof.value+6:] = np.dot(Cqq,dQdt[NumDof.value+6:])

            Qsys += np.dot(Msys,dQddt)
            
            # include damping
            if XBINPUT.str_damping_model == 'prop':
                Cdamp = XBINPUT.str_damping_param['alpha'] * MssFull + \
                        XBINPUT.str_damping_param['beta']  * KssFull
                Csys[:NumDof.value,:NumDof.value] += Cdamp
                Qsys[:NumDof.value] += np.dot(Cdamp, dQdt[:NumDof.value])
                                

            # 
            # special BCs
            # if  SphFlag:
            if len(iiblock)>0: # allow to enforce only attitude while keeping velocity free
                Msys[iiblock,:] = 0.0
                Msys[iiblock,iiblock] = 1.0
                Csys[iiblock,:] = 0.0
                Ksys[iiblock,:] = 0.0
                Qsys[iiblock]   = 0.0
                if XBINPUT.sph_joint_damping is not None:
                    Csys[iirotfree,iirotfree] += XBINPUT.sph_joint_damping
                    Qsys[iirotfree] += XBINPUT.sph_joint_damping*dQdt[iirotfree]
          
            #Calculate system matrix for update calculation
            Asys = Ksys + Csys*gamma/(beta*dt) + Msys/(beta*dt**2)
            
            #Compute correction
            
            ###DQ[:] = np.dot(np.linalg.inv(Asys), -Qsys)
            DQ[:] = np.linalg.solve(Asys,-Qsys)

            Q += DQ
            dQdt += DQ*gamma/(beta*dt)
            dQddt += DQ/(beta*dt**2)
            
            
            #Update convergence criteria
            if XBOPTS.PrintInfo.value==True:                 
                sys.stdout.write('%-10.4e ' %(max(abs(Qsys))))
            
            Res_Qglobal = max(abs(Qsys))
            Res_DeltaX  = max(abs(DQ))
            
            ResLog10 = max(Res_Qglobal/Res0_Qglobal,Res_DeltaX/Res0_DeltaX)
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('%-10.4e %8.4f\n' %(max(abs(DQ)),ResLog10))

            if SaveExtraVariables is True:
                if Iter == 1:
                    XBOUT.AsysListStart.append(Asys.copy())
                if ( (ResLog10 < XBOPTS.MinDelta.value) or (Iter >= XBOPTS.MaxIterations.value) ):
                    XBOUT.AsysListEnd.append(Asys.copy())
                    XBOUT.MsysList.append(Msys.copy())
                    XBOUT.CsysList.append(Csys.copy())
                    XBOUT.KsysList.append(Ksys.copy())
        # END Netwon-Raphson.
        
        
        # sm: here to avoid crash at first time-step
        if iStep > 0 and AELAOPTS.Tight == False:
            XBOUT.ForceAeroList.append(ForceAero.copy('C')) 
  
        # sm: save aero data
        if ( SaveDict['SaveWake']==True           and 
             iStep%SaveDict['SaveWakeFreq'] == 0  and 
             iStep==2                                 ):
            nfile=iStep//SaveDict['SaveWakeFreq']
            hdwake=h5py.File(dirwake+'%.4d.h5'%nfile,'w')
            hdwake['iStep']=iStep
            hdwake['Zeta']= np.float32(Zeta.copy('C'))
            hdwake['ZetaStar']= np.float32(ZetaStar.copy('C'))            
            hdwake.close()

        
        # sm debug:
        #XBOUT.ForceDofList.append( np.dot(FstrucFull, Force_Dof).copy() )
        XBOUT.ForceRigidList.append( np.dot(FrigidFull, Force_All).copy() )
              
        #update to converged nodal displacements and velocities
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F'); 
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
        
        PosPsiTime[iStep+1,:3] = PosDefor[-1,:]
        PosPsiTime[iStep+1,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        #Position of all grid points in global FoR
        i1 = (iStep+1)*NumNodes_tot.value
        i2 = (iStep+2)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
        
        #Export rigid-body velocities/accelerations
        if XBOPTS.OutInaframe.value==True:
            Vrel[iStep,:] = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
        else:
            Quat = dQdt[NumDof.value+6:].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            Cao  = xbl.Rot(Quat)
            ACoa[:3,:3] = np.transpose(Cao)
            ACoa[3:,3:] = np.transpose(Cao)
            
            Vrel[iStep,:] = np.dot(ACoa,dQdt[NumDof.value:NumDof.value+6].copy('F'))
            VrelDot[iStep,:] = np.dot(ACoa,dQddt[NumDof.value:NumDof.value+6].copy('F'))
        
        
        if 'writeDict' in kwords and Settings.WriteOut == True:
            fp= write_SOL912_out(Time[iStep+1], PosDefor, PsiDefor, PosIni, PsiIni, XBELEM, 
                                 kwords['writeDict'], SaveDict, FileObject=fp)
                    
        # 'Rollup' due to external velocities. TODO: Must add gusts here!
        ZetaStar[:,:] = ZetaStar[:,:] + VMINPUT.U_infty*dt
        
        # sm: append outputs
        XBOUT.QuatList.append(Quat.copy())
        XBOUT.CRVList.append(PsiA_G.copy())
 
        # sm I/O: FoR A velocities/accelerations
        XBOUT.Time=Time                     # ...dyn.dat
        #XBOUT.PosPsiTime = PosPsiTime       
        
        XBOUT.DynOut=DynOut                 # ...shape.dat
        
        XBOUT.Vrel=Vrel                     # ...rigid.dat
        XBOUT.VrelDot=VrelDot
        #XBOUT.PosPsiTime=PosPsiTime          
        
        XBOUT.PsiList.append(PsiDefor.copy())   
        
        XBOUT.cputime.append( time.clock() - XBOUT.cputime[0] )
        
        if SaveDict['SaveProgress']:
            iisave=np.arange(1,NumSteps.value,np.ceil(NumSteps.value/SaveDict['NumSavePoints']))
            if any(iisave==iStep):
                PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5', 
                                      *OutList)
        
    # END Time loop
    
    
    if SaveDict['Format'] == 'dat': 
        write_SOL912_final(Time, PosPsiTime, NumNodes_tot, DynOut, Vrel, VrelDot, SaveDict) 
        
    # Close output file if it exists.
    if 'writeDict' in kwords and Settings.WriteOut == True:
        fp.close()
        
    # Close TecPlot ASCII FileObject.
    if Settings.PlotTec==True:
        PostProcess.CloseAeroTecFile(FileObject)
        
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')
        
    # For interactive analysis at end of simulation set breakpoint.
    pass

    # sm I/O: FoR A velocities/accelerations

    XBOUT.Time=Time                     # ...dyn.dat
    XBOUT.PosPsiTime = PosPsiTime       
    
    XBOUT.DynOut=DynOut                 # ...sgape.dat
    
    XBOUT.Vrel=Vrel                     # ...rigid.dat
    XBOUT.VrelDot=VrelDot
    XBOUT.PosPsiTime=PosPsiTime   
    
    if  SaveDict['SaveWake'] is True:
        #XBOUT.dirwake=dirwake  
        XBOUT.NTwake=NumSteps.value//SaveDict['SaveWakeFreq']
    
    #saveh5(SaveDict, AELAOPTS, VMINPUT, VMOPTS, XBOPTS, XBINPUT, XBOUT )
    PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5', *OutList)
    
    return XBOUT
        
