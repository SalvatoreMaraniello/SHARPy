'''@package PyBeam.Solver.NonlinearDynamic
@brief      Nonlinear dynamic solvers for unconstrained problems.
@author     Rob Simpson, Henrik Hesse, S. Maraniello
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       24/10/2013 
@pre        None
@warning    Large memory usage due to the division between ForceDyn, 
            ForceDyn_dead, ForceDyn_doll

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



def Solve_F90(XBINPUT,XBOPTS,**kwords):
    """@brief Nonlinear dynamic structural solver using f90 solve routine."""
    
    #Check correct solution code
    assert XBOPTS.Solution.value == 912, ('NonlinearDynamic (F90) requested' +\
                                              ' with wrong solution code')  
    
    # I/O management
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
    
    #Solve static
    # sm: to be defined
    print('static solution commented out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    solve_static=False
    if solve_static:
        BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                            PosIni, PsiIni, XBNODE, NumDof,\
                            PosDefor, PsiDefor)
    else:
        pass
    
    #Change solution code back to NonlinearDynamic
    XBOPTS.Solution.value = 912
    
    
    #Write deformed configuration to file
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
        
    
    #Initialise variables for dynamic analysis
    #### sm: previous code
    #Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    #PosDotDef, PsiDotDef,\
    #OutGrids, PosPsiTime, VelocTime, DynOut\
    #    = BeamInit.Dynamic(XBINPUT,XBOPTS, NumNodes_tot)
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
    = BeamInit.Dynamic(XBINPUT,XBOPTS)
    
    # Delete unused variables.
    del PosPsiTime, VelocTime
    
    #Initialise quaternions for rigid-body dynamics
    Quat = np.zeros(4, ct.c_double, 'F')
    Quat[0] = 1.0
    
    
    #Write _force file
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_force.dat'
    fp = open(ofile,'w')
    BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
    fp.close() 

    # sm write class
    XBOUT=DerivedTypes.Xboutput()
    XBOUT.QuatList.append(Quat)
    XBOUT.Time=Time
    XBOUT.PosIni=PosIni
    XBOUT.PsiIni=PsiIni  
    #XBOUT.ForcedVel=ForcedVel        # ...SOL912_force.dat
    #XBOUT.ForcedVelDot=ForcedVelDot
    XBOUT.ForceTime=ForceTime
        
    #Solve dynamic using f90 solver
    BeamLib.Xbeam_Solv_FreeNonlinDynamic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
            PosIni, PsiIni, XBNODE, NumDof,\
            PosDefor, PsiDefor, Quat,\
            NumSteps, Time, ForceTime, ForcedVel, ForcedVelDot,\
            PosDotDef, PsiDotDef, DynOut, OutGrids)
    
    
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
    
    # collect output for saving
    XBOUT.DynOut=DynOut
    XBOUT.ForcedVel=ForcedVel
    XBOUT.ForcedVelDot=ForcedVelDot
    
    try:
        h5filename=Settings.OutputDir + Settings.OutputFileRoot + '.h5'
        hdfile=h5py.File(h5filename,'w')
        print ('created %s'%h5filename)
        PyLibs.io.save.add_class_as_grp(XBOPTS,hdfile)
        PyLibs.io.save.add_class_as_grp(XBINPUT,hdfile)
        PyLibs.io.save.add_class_as_grp(XBOUT,hdfile)
        hdfile.close()
    except:
        pass  
    
    return XBOUT     


def Solve_Py(XBINPUT,XBOPTS,**kwords):
    """Nonlinear dynamic structural solver in Python. Assembly of matrices 
    is carried out with Fortran subroutines."""
    
    #Check correct solution code
    assert XBOPTS.Solution.value == 912, ('NonlinearDynamic requested' +\
                                              ' with wrong solution code')
    
    # I/O management
    XBOUT=DerivedTypes.Xboutput()  
    SaveDict=Settings.SaveDict
    if 'SaveDict' in kwords: SaveDict=kwords['SaveDict']
    if SaveDict['Format']=='h5':
        Settings.WriteOut=False
        Settings.PlotTec=False
        OutList=[XBOPTS, XBINPUT, XBOUT]
            
    XBOUT.cputime.append(time.clock()) # or time.processor_time
    XBOUT.ForceDofList=[]
    XBOUT.ForceRigidList=[]
            
    #Initialise beam
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    # special BCs
    SphFlag=False
    if XBNODE.Sflag.any(): SphFlag=True
     
    #Solve static
    ### sm commented
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    print('STATIC SOLUTION COMMENTED OUT!!!')
    solve_static=False
    if solve_static:
        XBOPTS.Solution.value = 112 
        PosDefor, PsiDefor = Solve_Py_Static(XBINPUT,XBOPTS)
    
    XBOPTS.Solution.value = 912 

    if SaveDict['Format']=='dat': 
        PyLibs.io.dat.write_SOL912_def(XBOPTS,XBINPUT,XBELEM,
                                NumNodes_tot,PosDefor,PsiDefor,SaveDict)
    
    # sm I/O
    XBOUT.PosDefor=PosDefor
    XBOUT.PsiDefor=PsiDefor
    
    
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
    
    # sm: consistency amongst solvers
    #Quat = np.zeros(4, ct.c_double, 'F'); Quat[0] = 1.0
    Quat =  xbl.psi2quat(XBINPUT.PsiA_G)
    Cao  = XbeamLib.Rot(Quat)
    ACoa = np.zeros((6,6), ct.c_double, 'F')
    Cqr = np.zeros((4,6), ct.c_double, 'F')
    Cqq = np.zeros((4,4), ct.c_double, 'F')
        
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4):
        Unit4[i,i] = 1.0
        
    
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
    
    ### sm: removed ForceTime and increased rank of ForceDyn
    # Force at the first time-step
    #Force = (XBINPUT.ForceStatic + XBINPUT.ForceDyn*ForceTime[0]).copy('F')
    Force = (XBINPUT.ForceStatic + XBINPUT.ForceDyn[0,:,:]).copy('F')
    

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
    
    if SphFlag:

        # block translations (redundant:)
        for ii in range(3):
            if XBINPUT.EnforceTraVel_FoRA[ii] == True:
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

        # add damp at the spherical joints
        if XBINPUT.sph_joint_damping is not None:
            Qsys[iirotfree]+= XBINPUT.sph_joint_damping*dQdt[iirotfree]
         
    # add structural damping term
    if XBINPUT.str_damping_model is not None:
        Cdamp = XBINPUT.str_damping_param['alpha'] * MssFull + \
                XBINPUT.str_damping_param['beta']  * KssFull
        Qsys[:NumDof.value] += np.dot( Cdamp, dQdt[:NumDof.value] )                  
        pass
     
    # -------------------------------------------------------------------  


    #Initial Accel
    #dQddt[:] = np.dot(np.linalg.inv(Msys), -Qsys)
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
        

        #Predictor step
        Q += dt*dQdt + (0.5-beta)*dQddt*dt**2
        dQdt += (1.0-gamma)*dQddt*dt
        dQddt[:] = 0.0
        
        ### sm: removed ForceTime and increased rank of ForceDyn 
        #Force at current time-step
        #Force = (XBINPUT.ForceStatic+XBINPUT.ForceDyn*ForceTime[iStep+1]).copy('F')
        Force = (XBINPUT.ForceStatic+XBINPUT.ForceDyn[iStep+1,:,:]).copy('F')
        
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
            tmpForceTime=ForceTime[iStep+1].copy('F') 
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
        
        # sm debug
        XBOUT.ForceRigidList.append( np.dot(FrigidFull, Force_All).copy() )
        
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
    XBINPUT.ForceDyn = XBINPUT.ForceDyn + XBINPUT.ForceDyn_foll +  XBINPUT.ForceDyn_dead
    del(XBINPUT.ForceDyn_dead)
    del(XBINPUT.ForceDyn_foll)
    PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5', *OutList)
     
    
    return XBOUT    
    
    
