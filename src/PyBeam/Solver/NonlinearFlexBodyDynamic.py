'''@package PyBeam.Solver.NonlinearDynamic
@brief      Nonlinear dynamic solvers for unconstrained problems.
@author     Rob Simpson, Henrik Hesse, S. Maraniello
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       24/10/2013 
@pre        None
@warning    Memory usage not optimised (see ForceDyn, ForceDyn_dead, 
            ForceDyn_doll variables)
'''

import sys
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import BeamLib
import BeamInit
import numpy as np
import ctypes as ct
import XbeamLib
from DerivedTypes import dump
from PyBeam.Solver.NonlinearStatic import Solve_Py as Solve_Py_Static
import PyBeam.Utils.XbeamLib as xbl
import h5py, time  
import PyLibs.io.save, PyLibs.io.dat
from IPython import embed
from XbeamLib import AddGravityLoads


def Solve_F90(XBINPUT,XBOPTS,**kwords):

    """@brief Nonlinear dynamic structural solver using f90 solve routine."""
    
    #Check correct solution code
    assert XBOPTS.Solution.value == 912, ('NonlinearDynamic (F90) requested' +\
                                              ' with wrong solution code')  

    # I/O management
    XBOUT=DerivedTypes.Xboutput()
    XBOUT.cputime.append(time.clock())  
    OutList=[XBOPTS, XBINPUT, XBOUT]  
    if 'SaveDict' in kwords:
        SaveDict=kwords['SaveDict']
        Settings.OutputDir     =SaveDict['OutputDir']
        Settings.OutputFileRoot=SaveDict['OutputFileRoot']     

    #Initialise variables for static analysis
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT, XBOPTS)
    
    #Change solution code to NonlinearStatic
    XBOPTS.Solution.value = 112
    
    #Set initial conditions as undef config
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    # Add Gravity - record applied ForceStatic for I/O
    if XBOPTS.FollowerForceRig is True and np.linalg.abs(XBINPUT.g)>1e-14:
        raise NameError('For automatically include gravity loads, '
                  'XBOPTS.FollowerForceRig=False is required.\n'
                  'Use PyCoupled.NonlinearFlexBodyDynamic.Solve_Py instead!')
    ForceStaticCopy = XBINPUT.ForceStatic.copy('F')
    AddGravityLoads(XBINPUT.ForceStatic,XBINPUT,XBELEM,
                    AELAOPTS=None, # allows defining inertial/elastic axis
                    PsiDefor=PsiDefor,
                    chord = 0.0, # used to define inertial/elastic axis
                    PsiA_G=None) # no G/A frame in this solution

    #Solve static:
    if XBOPTS.ImpStart == False:
        BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                            PosIni, PsiIni, XBNODE, NumDof,\
                            PosDefor, PsiDefor)
        #Change solution code back to NonlinearDynamic
        XBOPTS.Solution.value = 912

    XBOUT.drop( PosIni=PosIni, PsiIni=PsiIni, PosDeforStatic=PosDefor.copy(), 
                PsiDeforStatic=PsiDefor.copy() )

    #Write deformed configuration to file
    if SaveDict['Format']=='dat':
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_def.dat'
        if XBOPTS.PrintInfo==True:
            sys.stdout.write('Writing file %s ... ' %(ofile))
        fp = open(ofile,'w')
        fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
        fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
        fp.close()
        if XBOPTS.PrintInfo==True:
            sys.stdout.write('done\n')
        WriteMode = 'a'
        
        BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                           PosDefor, PsiDefor, ofile, WriteMode)
    else:
        PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5',
                              *OutList)

        
    #Initialise variables for dynamic analysis
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
    = BeamInit.Dynamic(XBINPUT,XBOPTS)

    # Delete unused variables.
    del PosPsiTime, VelocTime
    
    ###Initialise quaternions for rigid-body dynamics
    #Quat = np.zeros(4, ct.c_double, 'F')
    #Quat[0] = 1.0
    if XBINPUT.PsiA_G_dyn != None: PsiA_G = XBINPUT.PsiA_G_dyn.copy()
    else: PsiA_G = XBINPUT.PsiA_G.copy()
    Quat=XbeamLib.psi2quat(PsiA_G)
    
    
    #Write _force file
    if SaveDict['Format']=='dat':
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_force.dat'
        fp = open(ofile,'w')
        BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
        fp.close() 

    # sm write class
    XBOUT.QuatList.append(Quat)
    XBOUT.Time=Time
    #XBOUT.ForcedVel=ForcedVel        # ...SOL912_force.dat
    #XBOUT.ForcedVelDot=ForcedVelDot
    XBOUT.ForceTime=ForceTime
        
    #Solve dynamic using f90 solver
    BeamLib.Xbeam_Solv_FreeNonlinDynamic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
            PosIni, PsiIni, XBNODE, NumDof,\
            PosDefor, PsiDefor, Quat,\
            NumSteps, Time, ForceTime, ForcedVel, ForcedVelDot,\
            PosDotDef, PsiDotDef, DynOut, OutGrids)
    
    # Divide gravity and applied load for I/O
    XBOUT.ForceTotal = XBINPUT.ForceStatic.copy('F')
    XBINPUT.ForceStatic = ForceStaticCopy

    # collect output for saving
    XBOUT.DynOut=DynOut
    XBOUT.ForcedVel=ForcedVel
    XBOUT.ForcedVelDot=ForcedVelDot

    if SaveDict['Format']=='dat':
        #Write _shape file
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_shape.dat'
        fp = open(ofile,'w')
        BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
        fp.close()
    
        #Write rigid-body velocities
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_rigid.dat'
        fp = open(ofile,'w')
        BeamIO.Write_rigid_File(fp, Time, ForcedVel, ForcedVelDot)
        fp.close()
    else:
        PyLibs.io.save.h5file(SaveDict['OutputDir'], 
                                     SaveDict['OutputFileRoot']+'.h5', *OutList)


    return XBOUT     


def Solve_Py(XBINPUT,XBOPTS,SaveDict=Settings.SaveDict):

    """Nonlinear dynamic structural solver in Python. Assembly of matrices 
    is carried out with Fortran subroutines."""
    
    #Check correct solution code
    assert XBOPTS.Solution.value == 912, ('NonlinearDynamic requested' +\
                                              ' with wrong solution code')
    
    # I/O management
    Settings.SaveDict=SaveDict # overwrite for compatibility with 'dat' format output
    XBOUT=DerivedTypes.Xboutput()  
    XBOUT.cputime.append(time.clock())
    if SaveDict['Format']=='h5':
        Settings.WriteOut=False
        Settings.PlotTec=False
        OutList=[XBOPTS, XBINPUT, XBOUT]
            
    #Initialise beam
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
     
    # Solve static
    # Note: gravity loads added according to XBINPUT.PsiA_G into Solve_Py_Static
    if XBOPTS.ImpStart==True:
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        XBOUT.ForceStaticTotal = XBINPUT.ForceStatic.copy('F')
        AddGravityLoads(XBOUT.ForceStaticTotal,XBINPUT,XBELEM,
                        AELAOPTS=None, # allows defining inertial/elastic axis
                        PsiDefor=PsiDefor,
                        chord=0.0, # used to define inertial/elastic axis
                        PsiA_G=XBINPUT.PsiA_G,
                        FollowerForceRig=XBOPTS.FollowerForceRig.value)
    else:
        if XBNODE.Sflag.any():
            raise NameError('Static solution with spherical joint not '
                                                 'implemented and/or possible!')
        XBOPTS.Solution.value = 112
        XBSTA = Solve_Py_Static(XBINPUT, XBOPTS, SaveDict=SaveDict)
        XBOPTS.Solution.value = 912
        PosDefor=XBSTA.PosDeforStatic.copy(order='F')
        PsiDefor=XBSTA.PsiDeforStatic.copy(order='F')
        XBOUT.PosDeforStatic  =XBSTA.PosDeforStatic
        XBOUT.PsiDeforStatic  =XBSTA.PsiDeforStatic
        XBOUT.ForceStaticTotal=XBSTA.ForceStaticTotal # includes gravity
        del XBSTA


    if SaveDict['Format']=='dat': 
        PyLibs.io.dat.write_SOL912_def(XBOPTS,XBINPUT,XBELEM,
                                NumNodes_tot,PosDefor,PsiDefor,SaveDict)
    
    ## sm I/O
    #XBOUT.PosDefor=PosDefor
    #XBOUT.PsiDefor=PsiDefor

    #Initialise variables for dynamic analysis
    Time, NumSteps, ForceTime, Vrel, VrelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut = BeamInit.Dynamic(XBINPUT,XBOPTS)
    # Delete unused variables.
    del OutGrids, VelocTime
    
    #Write _force file
    if SaveDict['Format']=='dat': 
        PyLibs.io.dat.write_SOL912_force(XBOPTS,XBINPUT,XBELEM,
                                         Time, ForceTime, Vrel, VrelDot)
    
    # sm I/O
    ### why forced velocity with Sol912 ???
    ### If forced velocities are prescribed, then is Sol312
    XBOUT.Time_force=Time               # ...SOL912_force.dat
    XBOUT.ForceTime_force=ForceTime
    XBOUT.Vrel_force=Vrel
    XBOUT.VrelDot_force=VrelDot  
    # debugging 
    XBOUT.KsysList=[]
    XBOUT.MsysList=[]
    
    #---------------------------------------------------- Start Dynamic Solution

    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
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
    
    
    #Initialise rotation operators 
    Quat =  xbl.psi2quat(XBINPUT.PsiA_G)
    Cao  = XbeamLib.Rot(Quat)
    ACoa = np.zeros((6,6), ct.c_double, 'F')
    Cqr = np.zeros((4,6), ct.c_double, 'F')
    Cqq = np.zeros((4,4), ct.c_double, 'F')
        
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4): Unit4[i,i] = 1.0
        
    #Extract initial displacements and velocities
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,
                          PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                          X, dXdt)
          
    #Initialise accelerations as zero arrays
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                            ct.c_double, 'F')
    
    #Populate state vector
    Q[:NumDof.value]=X.copy('F')
    dQdt[:NumDof.value]=dXdt.copy('F')
    dQdt[NumDof.value:NumDof.value+6] = Vrel[0,:].copy('F')
    dQdt[NumDof.value+6:]= Quat.copy('F')
    
    # Force at the first time-step
    # Gravity needs including for correctly computing the accelerations at the 
    # 1st time-step
    # Note: rank of Force increased after removing ForceTime
    #Force = (XBINPUT.ForceStatic + XBINPUT.ForceDyn*ForceTime[0]).copy('F')
    #Force = (XBINPUT.ForceStatic + XBINPUT.ForceDyn[0,:,:]).copy('F')
    Force = (XBOUT.ForceStaticTotal + XBINPUT.ForceDyn[0,:,:]).copy('F')


    #Assemble matrices and loads for structural dynamic analysis
    tmpVrel=Vrel[0,:].copy('F')
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
    
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )

    Qstruc -= np.dot(FstrucFull, Force_Dof)
    
    
    #Assemble matrices for rigid-body dynamic analysis
    BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         tmpVrel, 0*tmpVrel, tmpQuat,\
                         NumDof, Settings.DimMat,\
                         mr, MrsFull, Mrr,\
                         cr, CrsFull, Crr, Cqr, Cqq,\
                         kr, KrsFull, fr, FrigidFull,\
                         Qrigid, XBOPTS, Cao)
    
    BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                               ct.byref(ct.c_int(6)),\
                               Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                               ct.byref(ct.c_int(NumDof.value+6)),\
                               Force_All.ctypes.data_as(ct.POINTER(ct.c_double)) )

    Qrigid -= np.dot(FrigidFull, Force_All)
    
          
    #Separate assembly of follower and dead loads   
    tmpForceTime=ForceTime[0].copy('F') 
    tmpQforces,Dummy,tmpQrigid = XbeamLib.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                                    PosIni, PsiIni, PosDefor, PsiDefor, \
                                    ### sm increased rank of ForceDyn_*
                                    #(XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
                                    #(XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
                                    (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll[0,:,:]), \
                                    (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead[0,:,:]), \
                                    Cao,1)
                           
    Qstruc -= tmpQforces      
    Qrigid -= tmpQrigid
    
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
    if XBNODE.Sflag.any():
        # block translations (redundant:)
        for ii in range(3):
            if XBINPUT.EnforceTraVel_FoRA[ii] == True: iiblock.append(NumDof.value+ii)
        # block rotations
        iirotfree=[] # free rotational dof 
        for ii in range(3):
            if XBINPUT.EnforceAngVel_FoRA[ii] is True: iiblock.append(NumDof.value+3+ii)
            else: iirotfree.append(NumDof.value+3+ii)

    if len(iiblock)>0:
        Msys[iiblock,:] = 0.0
        Msys[:,iiblock] = 0.0 
        Msys[iiblock,iiblock] = 1.0
        Csys[iiblock,:] = 0.0
        Ksys[iiblock,:] = 0.0
        Qsys[iiblock]   = 0.0
        if XBINPUT.sph_joint_damping is not None:
            Csys[iirotfree,iirotfree] += XBINPUT.sph_joint_damping
            Qsys[iirotfree] += XBINPUT.sph_joint_damping*dQdt[iirotfree]
    #from IPython import embed; embed()     
    # add structural damping term
    if XBINPUT.str_damping_model is not None:
        Cdamp = XBINPUT.str_damping_param['alpha'] * MssFull + \
                XBINPUT.str_damping_param['beta']  * KssFull
        Qsys[:NumDof.value] += np.dot( Cdamp, dQdt[:NumDof.value] )                  
        pass
     
    # -------------------------------------------------------------------  
    
    #store initial matrices for eigenvalues analysis
    XBOUT.Msys0 = Msys.copy()
    XBOUT.Csys0 = Csys.copy()
    XBOUT.Ksys0 = Ksys.copy()
    XBOUT.Qsys0 = Qsys.copy()

    #Initial Accel
    #dQddt[:] = np.dot(np.linalg.inv(Msys), -Qsys)
    #dQddt[:] = np.linalg.solve(Msys,-Qsys) # most correct but inaccurate / Msys ill-conditioned
    # solve only rigid body dynamics
    dQddt[:NumDof.value]=0.0
    dQddt[NumDof.value:]=np.linalg.solve( Msys[NumDof.value:,NumDof.value:],
                                                          -Qsys[NumDof.value:] ) 


    XBOUT.dQddt0=dQddt.copy()
    
    #Record position of all grid points in global FoR at initial time step
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    #Position/rotation of the selected node in initial deformed configuration
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    
    #Get gamma and beta for Newmark scheme
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*(gamma + 0.5)**2
    
    # sm write class
    XBOUT.QuatList.append(Quat)
    XBOUT.PosIni=PosIni
    XBOUT.PsiIni=PsiIni 
   
    XBOUT.AsysListStart=[]
    XBOUT.AsysListEnd=[]
    XBOUT.MsysList=[]
    XBOUT.CsysList=[]
    XBOUT.KsysList=[]
    
    #Time loop
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('  Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        
        #calculate dt
        dt = Time[iStep+1] - Time[iStep]
        
        ###Predictor step
        Q += dt*dQdt + (0.5-beta)*dQddt*dt**2
        dQdt += (1.0-gamma)*dQddt*dt
        ### Corrector
        #dQddt[:] = 0.0 # initial guess for acceleration at next time-step is zero
        Q += beta*dQddt*dt**2
        dQdt += gamma*dQddt*dt

        #Reset convergence parameters
        Iter = 0
        ResLog10 = 1.0

        #Newton-Raphson loop      
        while ( (ResLog10 > XBOPTS.MinDelta.value) \
                & (Iter < XBOPTS.MaxIterations.value) ):
            
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
            Ksys[:,:] = 0.0; Asys[:,:] = 0.0;
            Qsys[:] = 0.0
            
            
            #Update counter
            Iter += 1
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('   %-7d ' %(Iter))

            #nodal displacements and velocities from state vector
            X=Q[:NumDof.value].copy('F') 
            dXdt=dQdt[:NumDof.value].copy('F'); 
            BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
            
            
            #rigid-body velocities and orientation from state vector
            Vrel[iStep+1,:] = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep+1,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
            Quat = dQdt[NumDof.value+6:].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            Cao  = XbeamLib.Rot(Quat)


            ### sm: removed ForceTime and increased rank of ForceDyn 
            # Note: do not add gravity load here!
            #Force at current time-step
            #Force = (XBINPUT.ForceStatic+XBINPUT.ForceDyn*ForceTime[iStep+1]).copy('F')
            Force = (XBINPUT.ForceStatic+XBINPUT.ForceDyn[iStep+1,:,:]).copy('F')

            # Add gravity loads (accounting for new orientation)
            AddGravityLoads(Force, XBINPUT, XBELEM,
                        AELAOPTS=None, # allows defining inertial/elastic axis
                        PsiDefor=PsiDefor,
                        chord = 0.0, # used to define inertial/elastic axis
                        PsiA_G=xbl.quat2psi(Quat),
                        FollowerForceRig=XBOPTS.FollowerForceRig.value)

            #Assemble matrices and loads for structural dynamic analysis
            tmpVrel=Vrel[iStep+1,:].copy('F')
            tmpQuat=Quat.copy('F')
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                                 PosIni, PsiIni, PosDefor, PsiDefor,\
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                                 Force, tmpVrel, 0*tmpVrel,\
                                 NumDof, Settings.DimMat,\
                                 ms, MssFull, Msr,\
                                 cs, CssFull, Csr,\
                                 ks, KssFull, fs, FstrucFull,\
                                 Qstruc, XBOPTS, Cao)
            
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
                    
            
            #Assemble matrices for rigid-body dynamic analysis
            BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                                 PosIni, PsiIni, PosDefor, PsiDefor,\
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                                 tmpVrel, 0*tmpVrel, tmpQuat,\
                                 NumDof, Settings.DimMat,\
                                 mr, MrsFull, Mrr,\
                                 cr, CrsFull, Crr, Cqr, Cqq,\
                                 kr, KrsFull, fs, FrigidFull,\
                                 Qrigid, XBOPTS, Cao)
    
            BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                                       ct.byref(ct.c_int(6)),\
                                       Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                       ct.byref(ct.c_int(NumDof.value+6)),\
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
            
          
            #Separate assembly of follower and dead loads   
            #tmpForceTime=ForceTime[iStep+1].copy('F') 
            tmpQforces,Dummy,tmpQrigid = XbeamLib.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                                            PosIni, PsiIni, PosDefor, PsiDefor, \
                                            ### sm: increased rank of ForceDyn_*
                                            #(XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
                                            #(XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
                                            (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll[iStep+1,:,:] ), \
                                            (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead[iStep+1,:,:] ), \
                                            Cao,1)
                                   
            Qstruc -= tmpQforces      
            Qrigid -= tmpQrigid
    
            #Compute residual
            Qstruc += -np.dot(FstrucFull, Force_Dof)
            Qrigid += -np.dot(FrigidFull, Force_All)
            
            # final of last iter
            XBOUT.Qstruc=Qstruc.copy()
            XBOUT.Qrigid=Qrigid.copy()
            # final of initial iter
            if iStep==0:
                XBOUT.Qstruc0=Qstruc.copy()
                XBOUT.Qrigid0=Qrigid.copy()
  
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

                #XBOUT.Msys=Msys.copy('F')
                #XBOUT.Qsys=Qsys.copy('F')
                #XBOUT.Csys=Csys.copy('F')
                #XBOUT.Ksys=Ksys.copy('F')
  
            #Calculate system matrix for update calculation
            Asys = Ksys + Csys*gamma/(beta*dt) + Msys/(beta*dt**2)

            #Compute correction
            DQ[:] = np.dot(np.linalg.inv(Asys), -Qsys)

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
                
        #END Newton-Raphson
        
        # debug
        #XBOUT.ForceRigidList.append( np.dot(FrigidFull, Force_All).copy() )
        
        #update to converged nodal displacements and velocities
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F'); 
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
        
        PosPsiTime[iStep+1,:3] = PosDefor[(NumNodes_tot.value-1)/2+1,:]
        PosPsiTime[iStep+1,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        #Position of all grid points in global FoR
        i1 = (iStep+1)*NumNodes_tot.value
        i2 = (iStep+2)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
        
        #Export rigid-body velocities/accelerations
        if XBOPTS.OutInaframe.value==True:
            Vrel[iStep+1,:] = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep+1,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
        else:
            Quat = dQdt[NumDof.value+6:].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            Cao  = XbeamLib.Rot(Quat)
            ACoa[:3,:3] = np.transpose(Cao)
            ACoa[3:,3:] = np.transpose(Cao)
            
            Vrel[iStep+1,:] = np.dot(ACoa,dQdt[NumDof.value:NumDof.value+6].copy('F'))
            VrelDot[iStep+1,:] = np.dot(ACoa,dQddt[NumDof.value:NumDof.value+6].copy('F'))
    
        # sm: append outputs
        XBOUT.QuatList.append(Quat.copy())
 
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
    
    #END Time loop

    if SaveDict['Format'] == 'dat': 
        PyLibs.io.dat.write_SOL912_final(Time, PosPsiTime, 
                                         NumNodes_tot, DynOut, Vrel, VrelDot, SaveDict) 
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')

    # sm I/O: FoR A velocities/accelerations
    XBOUT.Time=Time                     # ...dyn.dat
    XBOUT.PosPsiTime = PosPsiTime       
    
    XBOUT.DynOut=DynOut                 # ...shape.dat
    
    XBOUT.Vrel=Vrel                     # ...rigid.dat
    XBOUT.VrelDot=VrelDot
    XBOUT.PosPsiTime=PosPsiTime    
    
    # save h5 
    XBINPUT.ForceDyn = XBINPUT.ForceDyn + XBINPUT.ForceDyn_foll + XBINPUT.ForceDyn_dead
    del(XBINPUT.ForceDyn_dead)
    del(XBINPUT.ForceDyn_foll)
    PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5', *OutList)
     
    
    return XBOUT    
