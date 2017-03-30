'''@package PyCoupled.Coupled_NlnFlightDynamic
@brief      NonlinearDynamic Beam + Rigid-body dynamics + UVLM.
@author     Rob Simpson & Henrik Hesse, S. Maraniello
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       07/11/2013
@pre        None
@warning    Follower and dead forces are only modelled as external forces that follow the FoR A or 
            not. In the second case, forces are defined in  FoR G components. No structural 
            deformation is followed!

'''
#----------------------------------------------------------------------- Packages

import numpy as np
import ctypes as ct
import sys, os, h5py, time

import Main.SharPySettings as Settings
import DerivedTypes
import BeamLib
import BeamInit
import PyBeam.Utils.XbeamLib as xbl

from PyFSI.Beam2UVLM import InitSection
from PyFSI.Beam2UVLM import CoincidentGrid

from PyAero.UVLM.Utils import PostProcess
from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels
from PyAero.UVLM.Solver.VLM import InitSteadyWake

import PyCoupled.Coupled_NlnStatic as Static
from PyBeam.Utils.XbeamLib import AddGravityLoads

import PyLibs.io.save
from PyLibs.io.dat import write_SOL912_def, write_SOL912_final, \
                          write_SOL912_out, write_TecPlot

from PyCoupled.Coupled_NlnFlightDynamic_support_aug import update_UVLMsol, init_GEBM_tensors, \
                                                       assemble_GEBM_tensors, set_res_tight, \
                                                       lagsolver


def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS,**kwords):
    """
    @brief Nonlinear dynamic solver using Python to solve aeroelastic equation.
    @details Assembly of structural matrices is carried out with  Fortran subroutines. Aerodynamics 
             solved using PyAero\.UVLM.
    @param XBINPUT Beam inputs (for initialization in Python).
    @param XBOPTS Beam solver options (for Fortran).
    @param VMOPTS UVLM solver options (for C/C++).
    @param VMINPUT UVLM solver inputs (for initialization in Python).
    @param VMUNST Unsteady input information for aero solver.
    @param AELAOPTS Options relevant to coupled aeroelastic simulations.
    @param writeDict OrderedDict of 'name':tuple of outputs to write.
    
    @todo: Add list of variables with description
    
    """

    # Check correct solution code.
    assert XBOPTS.Solution.value == 912, ('NonlinearFlightDynamic requested' +
                                          ' with wrong solution code')
    # Check loads options
    assert XBOPTS.FollowerForceRig.value == True, ('For NonlinearFlightDynamic '
        'solution XBOPTS.FollowerForceRig = ct.c_bool(True) is required. \n'
        'Note that gravity is treated as a dead load by default.')


    # I/O options
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
            
    XBOUT.cputime.append(time.clock()) # time.processor_time more appropriate 
    XBOUT.ForceDofList=[]
    XBOUT.ForceRigidList=[]


    #------------------ Initial Displacement/Forces: ImpStart vs Static Solution
    
    # Initialise static beam data.
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    
    # Calculate initial displacements.
    if AELAOPTS.ImpStart == True:

        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        # Extract forces
        XBOUT.ForceStaticTotal = XBINPUT.ForceStatic.copy('F')
        AddGravityLoads(XBOUT.ForceStaticTotal,XBINPUT,XBELEM,
                        AELAOPTS=None, # allows defining inertial/elastic axis
                        PsiDefor=PsiDefor,
                        chord = 0.0, # used to define inertial/elastic axis
                        PsiA_G=XBINPUT.PsiA_G,
                        FollowerForceRig=True)

        ForceAero = 0.0*XBOUT.ForceStaticTotal.copy('C')

    else:
        XBOPTS.Solution.value = 112 # Modify options.
        VMOPTS.Steady = ct.c_bool(True)
        Rollup = VMOPTS.Rollup.value
        VMOPTS.Rollup.value = False

        if VMINPUT.ctrlSurf != None:
            # open-loop control
            for cc in range(len(VMINPUT.ctrlSurf)):
                VMINPUT.ctrlSurf[cc].update(0.0,iStep=0)
        # Solve Static Aeroelastic.
        # Note: output force includes gravity loads
        #PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, Force = \
        #    Static.Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS)
        XBSTA=Static.Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS)
        # Extract forces
        XBOUT.ForceStaticTotal=XBSTA.ForceStaticTotal.copy(order='F') # gravity/aero/applied
        ForceAero = XBSTA.ForceAeroStatic.copy(order='C')
        # Define Pos/Psi as fortran array (crash otherwise)
        PosDefor=XBSTA.PosDeforStatic.copy(order='F')
        PsiDefor=XBSTA.PsiDeforStatic.copy(order='F')
        XBOUT.PosDeforStatic = XBSTA.PosDeforStatic
        XBOUT.PsiDeforStatic = XBSTA.PsiDeforStatic
        # Wake variables
        Zeta=XBSTA.ZetaStatic
        ZetaStar=XBSTA.ZetaStarStatic
        Gamma=XBSTA.GammaStatic
        GammaStar=XBSTA.GammaStarStatic
        del XBSTA

        XBOPTS.Solution.value = 912 # Reset options.
        VMOPTS.Steady = ct.c_bool(False)
        VMOPTS.Rollup.value = Rollup


    # Save output   
    if SaveDict['Format']=='dat': 
        write_SOL912_def(XBOPTS,XBINPUT,XBELEM,NumNodes_tot,PosDefor,PsiDefor,SaveDict)
    else: # HDF5
        XBOUT.drop( PosIni=PosIni, PsiIni=PsiIni,PosDeforStatic=PosDefor.copy(), 
                    PsiDeforStatic=PsiDefor.copy() )
        XBOUT.ForceAeroList.append( ForceAero )
        PyLibs.io.save.h5file(SaveDict['OutputDir'], 
                                     SaveDict['OutputFileRoot']+'.h5', *OutList)

    
    #------------------------------------------- Initialise Structural Variables    

    # Build arrays for dynamic analysis
    #   1. initial accelerations  {Pos/Psi}DotDotDef approximated to zero
    #   2. {Pos/Psi}DotDotDef remain zero
    ( Time, NumSteps, ForceTime, Vrel, VrelDot, PosDotDef, PsiDotDef, 
      OutGrids, PosPsiTime, VelocTime, DynOut ) = BeamInit.Dynamic( XBINPUT, XBOPTS)
    del OutGrids, VelocTime
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),ct.c_double, 'F')
    

    # Allocate memory for GEBM tensors:
    #   1. all set to zero
    #   2. all arrays as ct_c_double with Fortran ordering
    if 'S' in XBINPUT.BConds: Nconstr=2
    else: Nconstr=0
    #Nconstr=0
    ( MssFull, CssFull, KssFull, FstrucFull, Qstruc,
      MrsFull, CrsFull, KrsFull, FrigidFull, Qrigid,
      Mrr, Crr, Msr, Csr, 
      Force_Dof, Force_All,
      X, dXdt, Q, DQ, dQdt, dQddt,
      Msys, Csys, Ksys, Asys, Qsys ) = init_GEBM_tensors(NumDof,Nconstr)
    
    
    # Extract initial displacements (and velocities ?).
    #   1. also for non impulsive start, this step only populates X...
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,
                                 PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                                  X, dXdt)
    

    # Assemble Tensors for Dynamic Solution
    #  1. extract variables used twice or more (Quat, currVel)
    #  2. update Force
    #  3. Populate state vectors (Q, dQdt)
    if XBINPUT.PsiA_G_dyn != None: PsiA_G = XBINPUT.PsiA_G_dyn.copy()
    else: PsiA_G = XBINPUT.PsiA_G.copy()
    Quat=xbl.psi2quat(PsiA_G)
    currVrel=Vrel[0,:].copy('F')   # also shared by aero initialisation
    
    # ForceStatic will include aero only in the non-impulsive case
    Force = XBOUT.ForceStaticTotal.copy('F') + XBINPUT.ForceDyn[0,:,:].copy('F')  
    
    Q[:NumDof.value]=X.copy('F')
    dQdt[:NumDof.value]=dXdt.copy('F')
    dQdt[NumDof.value:NumDof.value+6] = currVrel.copy('F')
    dQdt[NumDof.value+6:NumDof.value+10]= Quat.copy('F')
    
    (Msys, Csys, Ksys, Qsys) = assemble_GEBM_tensors( 
                                       0, # iStep
                                       NumDof, Nconstr, NumNodes_tot,
                                       XBELEM, XBNODE,
                                       XBINPUT, VMINPUT, XBOPTS, 
                                       AELAOPTS, VMOPTS,  
                                       currVrel,     
                                       PosIni, PsiIni,          
                                       PosDefor, PsiDefor,  
                                       PosDotDef, PsiDotDef,
                                       PosDotDotDef, PsiDotDotDef,
                                       Force,
                                       MssFull, CssFull, KssFull, FstrucFull,
                                       Msr, Csr,
                                       Force_Dof, Qstruc,
                                       MrsFull, CrsFull, KrsFull, FrigidFull,
                                       Mrr, Crr, Qrigid,
                                       Q, dQdt, dQddt, 
                                       Force_All,
                                       Msys, Csys, Ksys, Asys, Qsys )


    # I/O
    DynOut[0:NumNodes_tot.value,:] = PosDefor                # nodal position
    PosPsiTime[0,:3] = PosDefor[-1,:]                        # nodal position
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]  # nodal CRV
    XBOUT.drop( Vrel=Vrel, VrelDot=VrelDot, ForceTime=ForceTime, 
                AsysListStart =[], AsysListEnd=[], 
                MsysList=[], CsysList=[], KsysList=[])
    XBOUT.QuatList.append(Quat.copy())
    XBOUT.CRVList.append( PsiA_G.copy())
    XBOUT.drop( Msys0 = Msys.copy(), Csys0 = Csys.copy(), 
                                       Ksys0 = Ksys.copy(), Qsys0 = Qsys.copy())
    

    #--------------------------------------------------------------------- Initialise UVLM variables
    
    # 1. Section properties
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
    
    # Declare memory for Aero variables.
    K = VMOPTS.M.value*VMOPTS.N.value
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    AIC = np.zeros((K,K),ct.c_double,'C')
    BIC = np.zeros((K,K),ct.c_double,'C')
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Initialise A-frame location and orientation to be zero.
    OriginA_a = np.zeros(3,ct.c_double,'C')
    #PsiA_G=XBINPUT.PsiA_G.copy()
    
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
    
    # Define TecPlot stuff
    if Settings.PlotTec==True:
        FileObject=write_TecPlot( Zeta, ZetaStar, Gamma, GammaStar, NumSteps.value, 0, Time[0], 
                                  SaveDict)
    if 'writeDict' in kwords and Settings.WriteOut == True:
        fp= write_SOL912_out( Time[0], PosDefor, PsiDefor, PosIni, PsiIni, XBELEM, 
                              kwords['writeDict'], SaveDict)
            
    
    
    #------------------------------------------------------------------------------- Time-Loop Start
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    
    #Get gamma and beta for Newmark scheme
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*(gamma + 0.5)**2 
    
    
    # Initial Accelerations:
    #   1. These can be non-zero (e.g Impulsive start or unbalanced start: in these cases, the 
    # structure has an initially undeformed/statically-deformed configuration, with zero velocities
    # but accelerations are non-zero).
    #   2. lagged solution compensates for Msys mal-conditioned
    if XBOPTS.RigidDynamics is False and AELAOPTS.ImpStart is False:
        lagsol=True
        if lagsol==True: dQddt[:] = lagsolver(Msys,-Qsys, NumDof.value, MaxIter=1)
        else: dQddt[:] = np.linalg.solve(Msys,-Qsys)
    else:
        # solve only rigid body dynamics
        dQddt[:NumDof.value]=0.0
        dQddt[NumDof.value:]=np.linalg.solve( Msys[NumDof.value:,NumDof.value:],
                                             -Qsys[NumDof.value:] ) 
    XBOUT.dQddt0=dQddt.copy()


    for iStep in range(1,NumSteps.value+1): 
        # Note len(Time)=NumSteps+1
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Time: %-10.4e\n' %(Time[iStep]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        #calculate dt
        dt = Time[iStep]-Time[iStep-1]
        VMOPTS.DelTime = ct.c_double(dt)
        
        #Predictor step
        Q    += dt*dQdt + (0.5-beta)*dQddt*np.power(dt,2.0)
        dQdt += (1.0-gamma)*dQddt*dt
        
        ### Corrector
        #  1. assume initial guess for the acceleration dQddt(iStep) = dQddt(iStep-1)
        Q += beta*dQddt*dt**2
        dQdt += gamma*dQddt*dt
        
        #nodal displacements and velocities from state vector
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F'); 
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT,NumNodes_tot,XBELEM,XBNODE,
                                       PosIni,PsiIni,NumDof,X,dXdt,
                                       PosDefor,PsiDefor,PosDotDef,PsiDotDef)
        
        #Reset convergence parameters
        Iter = 0
        ResLog10 = 1.0
        AELAMinRes=0.0
        
        # Residual at previous time-step
        # warning: the convergence check is not well implemented: these variables are going to be 
        #  equal to 1 always.
        if XBOPTS.RigidDynamics==False: iicheck=[ii for ii in range(NumDof.value+10+Nconstr)]  # to check
        else: iicheck=[ii for ii in range(NumDof.value, NumDof.value+10+Nconstr)]
        Res0_Qglobal = max(max(abs(Qsys[iicheck])),1)
        Res0_DeltaX  = max(max(abs(DQ[iicheck])),1)
        
        #---------------------------------------------------------------------------- Newton-Raphson 
        while ( (ResLog10 > XBOPTS.MinDelta.value) & (Iter < XBOPTS.MaxIterations.value) ):
 
            # Unpack state variables
            Vrel[iStep,:]    = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
            Quat = dQdt[NumDof.value+6:NumDof.value+10].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            X=Q[:NumDof.value].copy('F') 
            dXdt=dQdt[:NumDof.value].copy('F')
            Cao  = xbl.Rot(Quat)
            PsiA_G = xbl.quat2psi(Quat)

            # Aero Force 
            #   The aerodynamic force is update until the residual does not fall below a prescribed 
            # factor, AELAMinRes. This is evaluated after Iter=0 to access the value of the residual 
            # ResLog10. According to the setting AELAOPTS.MinRes, AELAOPTS.MaxRes, AELAOPTS.LimRes 
            # the threshold can then be lowered or raised.
            
            if Iter==1: AELAMinRes = set_res_tight( ResLog10, AELAOPTS.MinRes, 
                                                    AELAOPTS.MaxRes,AELAOPTS.LimRes)
            
            if (AELAOPTS.Tight == False)  and (ResLog10>AELAMinRes or Iter==0):
                Force, Zeta, ZetaStar, Gamma, GammaStar = update_UVLMsol( iStep, dt, Time, Force,
                                            XBINPUT, VMINPUT, XBELEM, AELAOPTS, VMOPTS, 
                                            PsiA_G, Vrel[iStep,:].copy('F'), OriginA_a, Section,
                                            PosDefor, PsiDefor, PosDotDef, PsiDotDef, 
                                            Ufree, AIC, BIC,
                                            AeroForces, Zeta, ZetaDot, ZetaStar, Gamma, GammaStar )
                ForceAero = Force.copy()        
            else: Force[:,:] = ForceAero 

            # Add gravity loads.
            AddGravityLoads(Force, XBINPUT, XBELEM, AELAOPTS, PsiDefor, 
                                       VMINPUT.c, PsiA_G, FollowerForceRig=True)

            # Add prescribed loads
            Force += (XBINPUT.ForceStatic + XBINPUT.ForceDyn[iStep,:,:]).copy('F')
            # Dead forces are defined in FoR G, follower in FoR A.
            # None follows the structure as it deflects
                             
            #Update matrices and loads for structural dynamic analysis
            BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot,
                                           XBELEM, XBNODE,
                                           PosIni, PsiIni,
                                           NumDof, X, dXdt,
                                           PosDefor, PsiDefor,
                                           PosDotDef, PsiDotDef)           
            
            (Msys, Csys, Ksys, Qsys) =  assemble_GEBM_tensors(
                                             iStep, 
                                             NumDof, Nconstr, NumNodes_tot,
                                             XBELEM, XBNODE,
                                             XBINPUT, VMINPUT, XBOPTS, 
                                             AELAOPTS, VMOPTS,  
                                             Vrel[iStep,:].copy('F')  ,                   
                                             PosIni, PsiIni,
                                             PosDefor, PsiDefor,
                                             PosDotDef, PsiDotDef,
                                             PosDotDotDef, PsiDotDotDef,
                                             Force,
                                             MssFull, CssFull, KssFull, FstrucFull,
                                             Msr, Csr,
                                             Force_Dof, Qstruc,
                                             MrsFull, CrsFull, KrsFull, FrigidFull,
                                             Mrr, Crr, Qrigid,
                                             Q, dQdt, dQddt,
                                             Force_All,
                                             Msys, Csys, Ksys, Asys, Qsys )  
            
            #Calculate system matrix for update calculation
            Asys = Ksys + Csys*gamma/(beta*dt) + Msys/(beta*dt**2)
       
            #Compute correction
            if XBOPTS.RigidDynamics==False:
                DQ[:]  = np.linalg.solve(Asys,-Qsys)
            else:    
                DQ[:NumDof.value]=0.0
                ### correction is zero, no need of this
                #Qsys[NumDof.value:NumDof.value+6]+=np.dot(
                #                        Ksys[NumDof.value:NumDof.value+6,:NumDof.value],
                #                        DQ[:NumDof.value] )
                DQ[NumDof.value:]=np.linalg.solve( Asys[NumDof.value:,NumDof.value:],
                                                  -Qsys[NumDof.value:])  
            Q     += DQ
            dQdt  += DQ*gamma/(beta*dt)
            dQddt += DQ/(beta*dt**2)

            #Update convergence criteria
            Res_Qglobal = max(abs(Qsys[iicheck]))
            Res_DeltaX  = max(abs(DQ[iicheck]))
            ResLog10 = max(Res_Qglobal/Res0_Qglobal, Res_DeltaX/Res0_DeltaX)
            
            # Update counter.
            Iter += 1
            if XBOPTS.PrintInfo.value==True: 
                sys.stdout.write( '   %-7d %-10.4e %-10.4e %8.4f\n' 
                                  %(Iter, max(abs(Qsys)), max(abs(DQ)), ResLog10) )             
        # END Netwon-Raphson.
        

        
        #----------------------------------------------------------------------- Terminate time-Loop
        
        # Unpack state variables
        Vrel[iStep,:]    = dQdt[NumDof.value:NumDof.value+6].copy('F')
        VrelDot[iStep,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
        Quat = dQdt[NumDof.value+6:NumDof.value+10].copy('F')
        Quat = Quat/np.linalg.norm(Quat)
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F')
        Cao  = xbl.Rot(Quat)
        PsiA_G = xbl.quat2psi(Quat)
        

        # sm: here to avoid crash at first time-step
        if AELAOPTS.Tight == False:
            XBOUT.ForceAeroList.append(ForceAero.copy('C'))

        # sm: save aero data
        if ( SaveDict['SaveWake']==True           and 
             iStep%SaveDict['SaveWakeFreq'] == 0  ):
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
        
        PosPsiTime[iStep,:3] = PosDefor[-1,:]
        PosPsiTime[iStep,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        #Position of all grid points in global FoR
        i1 = (iStep)*NumNodes_tot.value
        i2 = (iStep+1)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
        
        #Export rigid-body velocities/accelerations
        if XBOPTS.OutInaframe.value==False:
            ACoa = np.zeros((6,6), ct.c_double, 'F')
            ACoa[:3,:3] = np.transpose(Cao)
            ACoa[3:,3:] = np.transpose(Cao)
            Vrel[iStep,:] = np.dot(ACoa,dQdt[NumDof.value:NumDof.value+6].copy('F'))
            VrelDot[iStep,:] = np.dot(ACoa,dQddt[NumDof.value:NumDof.value+6].copy('F'))
        
        
        if 'writeDict' in kwords and Settings.WriteOut == True:
            fp= write_SOL912_out(Time[iStep], PosDefor, PsiDefor, PosIni, PsiIni, XBELEM, 
                                 kwords['writeDict'], SaveDict, FileObject=fp)
                    
        # 'Rollup' due to external velocities. TODO: Must add gusts here!
        ZetaStar[:,:] = ZetaStar[:,:] + VMINPUT.U_infty*dt
        
        # sm: append outputs

 
        # sm I/O: FoR A velocities/accelerations
        XBOUT.drop( Time=Time, DynOut=DynOut, Vrel=Vrel, VrelDot=VrelDot  )
        XBOUT.PsiList.append(PsiDefor.copy())   
        XBOUT.QuatList.append(Quat.copy())
        XBOUT.CRVList.append(PsiA_G.copy())
        
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
    
    # sm I/O: FoR A velocities/accelerations
    XBOUT.drop( Time=Time, PosPsiTime=PosPsiTime, DynOut=DynOut, Vrel=Vrel, VrelDot=VrelDot  )
    if  SaveDict['SaveWake'] is True:
        XBOUT.NTwake=NumSteps.value//SaveDict['SaveWakeFreq']
    PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5', *OutList)
    
    return XBOUT
   
   
 
