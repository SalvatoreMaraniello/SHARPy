'''@package PyBeam.Solver.NonlinearDynamic
@brief      Nonlinear dynamic solvers.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       10/12/2012
@pre        None
@warning    None
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
import PyLibs.io.save
import time
from XbeamLib import AddGravityLoads, psi2quat
from IPython import embed


def Solve_F90(XBINPUT,XBOPTS,SaveDict=Settings.SaveDict):
    """@brief Nonlinear dynamic structural solver using f90 solve routine."""
    
    #"Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic (F90) requested' +\
                                              ' with wrong solution code')
    assert XBINPUT.ForceDynType=='modal', ('Fortran version of Sol 312 only sup'
        'ports modal dynamic forces input: modify ForceDynType or use Solve_Py')
    
    # I/O management
    Settings.SaveDict=SaveDict # overwrite for compatibility with 'dat' format output methods
    XBOUT=DerivedTypes.Xboutput()
    XBOUT.cputime.append(time.clock())  
    OutList=[XBOPTS, XBINPUT, XBOUT]  

    #"Initialise variables for static analysis"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT, XBOPTS)
    
    
    #"Change solution code to NonlinearStatic"
    XBOPTS.Solution.value = 112
    
    #"Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    # Add Gravity - record applied ForceStatic for I/O
    ForceStaticCopy = XBINPUT.ForceStatic.copy('F')
    AddGravityLoads(XBINPUT.ForceStatic,XBINPUT,XBELEM,
                    AELAOPTS=None, # allows defining inertial/elastic axis
                    PsiDefor=PsiDefor,
                    chord = 0.0, # used to define inertial/elastic axis
                    PsiA_G=None) # no G/A frame in this solution

    #"Solve static"
    BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                            PosIni, PsiIni, XBNODE, NumDof,\
                            PosDefor, PsiDefor)
    
    #"Change solution code back to NonlinearDynamic"
    XBOPTS.Solution.value = 312
    
    
    #"Write deformed configuration to file"
    if SaveDict['Format']=='dat':
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_def.dat'
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
        XBOUT.drop( PosIni=PosIni,PsiIni=PsiIni,PosDeforStatic=PosDefor.copy(), 
                    PsiDeforStatic=PsiDefor.copy(), QuatList=np.array([[1,0,0,0]]))
        PyLibs.io.save.h5file(SaveDict['OutputDir'],
                                     SaveDict['OutputFileRoot']+'.h5', *OutList)
    
    
    #"Change solution code to NonlinearDynamic"
    XBOPTS.Solution.value = 312
    
    
    #"Initialise variables for dynamic analysis"
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    
    XBOUT.Time=Time
    XBOUT.Vrel=ForcedVel
    XBOUT.VrelDot=ForcedVelDot
    XBOUT.ForceTime=ForceTime


    #"Write _force file"
    if SaveDict['Format']=='dat':
        ofile = Settings.OutputDir+Settings.OutputFileRoot+'_SOL312_force.dat'
        fp = open(ofile,'w')
        BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
        fp.close() 

    
    #"Solve dynamic using f90 solver"
    BeamLib.Cbeam3_Solv_NonlinearDynamic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
            PosIni, PsiIni, XBNODE, NumDof,\
            PosDefor, PsiDefor,\
            NumSteps, Time, ForceTime, ForcedVel, ForcedVelDot,\
            PosDotDef, PsiDotDef,\
            PosPsiTime, VelocTime, DynOut, OutGrids)

    XBOUT.DynOut=DynOut

    # Divide gravity and applied load for I/O
    XBOUT.ForceStaticTotal = XBINPUT.ForceStatic.copy('F')
    XBINPUT.ForceStatic = ForceStaticCopy

    if SaveDict['Format']=='dat':
        #"Write _dyn file"
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_dyn.dat'
        fp = open(ofile,'w')
        BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
        fp.close()
        #"Write _shape file"
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_shape.dat'
        fp = open(ofile,'w')
        BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
        fp.close()
    else:
        PyLibs.io.save.h5file(SaveDict['OutputDir'], 
                                     SaveDict['OutputFileRoot']+'.h5', *OutList)

    return XBOUT
    
    
def Solve_F90_steps(XBINPUT,XBOPTS):
    """@brief Nonlinear dynamic structural solver with time-steps controlled
    in Python but solved using F90 routines at each time-step."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic (F90) requested' +\
                                              ' with wrong solution code')
    
    
    "Initialise variables for static analysis"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT, XBOPTS)
    
    
    "Change solution code to NonlinearStatic"
    XBOPTS.Solution.value = 112
    
    
    "Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    
    "Solve static"
    BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                            PosIni, PsiIni, XBNODE, NumDof,\
                            PosDefor, PsiDefor)
    
    
    "Write deformed configuration to file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_def.dat'
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
    
    
    "Change solution code to NonlinearDynamic"
    XBOPTS.Solution.value = 312
    
    
    "Initialise variables for dynamic analysis"
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    
    
    "Write _force file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_force.dat'
    fp = open(ofile,'w')
    BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
    fp.close() 
    
    
    "Write _vel file"   
    #TODO: write _vel file
    
    
    "Write .mrb file"
    #TODO: write .mrb file
    
    
    "Save time-dependant variables with only one time step (2 time-steps)"
    TimeElements = Time[0:2].copy('F')
    Ftime = ForceTime[0:2].copy('F')
    Fvel = ForcedVel[0:2,:].copy('F')
    FvelDot = ForcedVelDot[0:2,:].copy('F')
    
    
    "Create Output vectors, one for each time-step"
    PosPsiTimeStep = PosPsiTime[0:2,:].copy('F')
    VelocTimeStep = VelocTime[0:2,:].copy('F')
    DynOutStep = DynOut[0:2*NumNodes_tot.value , : ].copy('F')
    
    
    "Initialise and approximate initial accelerations as zero arrays"
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                           ct.c_double, 'F')
    
    
    "Start Time-loop"
    for TimeStep in range(0,NumSteps.value):
        
        "Update 2-step time arrays"
        TimeElements = Time[TimeStep:TimeStep+2]
        Ftime = ForceTime[TimeStep:TimeStep+2]
        Fvel[:,:] = ForcedVel[TimeStep:TimeStep+2,:]
        FvelDot[:,:] = ForcedVelDot[TimeStep:TimeStep+2,:]
        
        
        "Solve dynamic step using f90 solver"
        BeamLib.Cbeam3_Solv_NonlinearDynamicAccel(XBINPUT, XBOPTS, NumNodes_tot,\
                XBELEM,PosIni, PsiIni, XBNODE, NumDof,\
                PosDefor, PsiDefor,\
                ct.c_int(1), TimeElements, Ftime, Fvel, FvelDot,\
                PosDotDef, PsiDotDef,\
                PosDotDotDef, PsiDotDotDef,\
                PosPsiTimeStep, VelocTimeStep, DynOutStep, OutGrids)
        
        
        "Append output information"
        PosPsiTime[TimeStep+1,:] = PosPsiTimeStep[1,:]
        VelocTime[TimeStep+1,:] = VelocTimeStep[1,:]
        DynOut[(TimeStep)*NumNodes_tot.value:(TimeStep+1)*NumNodes_tot.value-1,\
                : ] = DynOutStep[NumNodes_tot.value:2*NumNodes_tot.value-1]
    
    "END Time-loop"
    
    
    "Write _dyn file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
    
    "Write _shape file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_shape.dat'
    fp = open(ofile,'w')
    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
    fp.close()
    

def Solve_Py(XBINPUT,XBOPTS, moduleName = None, SaveDict=Settings.SaveDict):
    """Nonlinear dynamic structural solver in Python. Assembly of matrices 
    is carried out with Fortran subroutines."""
    
    #"Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic requested' +\
                                              ' with wrong solution code')
    
    # I/O management
    Settings.SaveDict=SaveDict # overwrite for compatibility with 'dat' format output
    XBOUT=DerivedTypes.Xboutput()
    XBOUT.cputime.append(time.clock())  
    if SaveDict['Format']=='h5':
        Settings.WriteOut=False
        Settings.PlotTec=False
        OutList=[XBOPTS, XBINPUT, XBOUT] 

    #"Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS, moduleName)
    
    # Solve static solution
    if XBOPTS.ImpStart==True:
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        XBOUT.ForceStaticTotal = XBINPUT.ForceStatic.copy('F')
        AddGravityLoads(XBOUT.ForceStaticTotal,XBINPUT,XBELEM,
                        AELAOPTS=None, # allows defining inertial/elastic axis
                        PsiDefor=PsiDefor,
                        chord = 0.0, # used to define inertial/elastic axis
                        PsiA_G=XBINPUT.PsiA_G,
                        FollowerForceRig=XBOPTS.FollowerForceRig.value)
    else:
        XBOPTS.Solution.value = 112 
        XBSTA = Solve_Py_Static(XBINPUT, XBOPTS, moduleName, SaveDict=SaveDict)
        XBOPTS.Solution.value = 312
        # Define Pos/Psi as fortran array (crash otherwise)
        PosDefor=XBSTA.PosDeforStatic.copy(order='F')
        PsiDefor=XBSTA.PsiDeforStatic.copy(order='F')
        XBOUT.PosDeforStatic=XBSTA.PosDeforStatic
        XBOUT.PsiDeforStatic=XBSTA.PsiDeforStatic
        XBOUT.ForceStaticTotal=XBSTA.ForceStaticTotal # includes gravity loads
        del XBSTA

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
        XBOUT.drop( PosIni=PosIni, PsiIni=PsiIni, PosDeforStatic=PosDefor.copy(), 
                    PsiDeforStatic=PsiDefor.copy())
        PyLibs.io.save.h5file(SaveDict['OutputDir'],SaveDict['OutputFileRoot']+'.h5',
                              *OutList)
    
    
    #"Initialise variables for dynamic analysis"
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS, moduleName)
        
    # Delete unused vars
    del OutGrids, VelocTime
    
    #"Write _force file"
    if SaveDict['Format']=='dat':
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_force.dat'
        fp = open(ofile,'w')
        BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
        fp.close() 

    # sm write class
    XBOUT.Time=Time
    XBOUT.Vrel=ForcedVel
    XBOUT.VrelDot=ForcedVelDot

    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    
    #"Initialise structural system tensors"
    MglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    CglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    KglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    FglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    Asys = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    
    ms = ct.c_int()
    cs = ct.c_int()
    ks = ct.c_int()
    fs = ct.c_int()
    
    Mvel = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Cvel = np.zeros((NumDof.value,6), ct.c_double, 'F')
    
    X     = np.zeros(NumDof.value, ct.c_double, 'F')
    DX    = np.zeros(NumDof.value, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof.value, ct.c_double, 'F')
    dXddt = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
    
    
    #"Initialise rotation operators"
    Unit = np.zeros((3,3), ct.c_double, 'F')
    for i in range(3):
        Unit[i,i] = 1.0
    
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4):
        Unit4[i,i] = 1.0
        
    Cao = Unit.copy('F')
    Temp = Unit4.copy('F')
    
    #Quat = np.zeros(4, ct.c_double, 'F')
    #Quat[0] = 1.0
    Quat =  psi2quat(XBINPUT.PsiA_G)
    XBOUT.QuatList.append(Quat)
    
    #"Extract initial displacements and velocities"
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,\
                          PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                          X, dXdt)
    
    
    #"Approximate initial accelerations"
    
    #"Initialise accelerations as zero arrays"
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                           ct.c_double, 'F')
    
    #"Force at the first time-step"
    # Note: rank of Force increased after removing ForceTime
    Force = (XBOUT.ForceStaticTotal + XBINPUT.ForceDyn[0,:,:]).copy('F')

    #"Assemble matrices for dynamic analysis"
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         Force, ForcedVel[0,:], ForcedVelDot[0,:],\
                         NumDof, Settings.DimMat,\
                         ms, MglobalFull, Mvel,\
                         cs, CglobalFull, Cvel,\
                         ks, KglobalFull, fs, FglobalFull,\
                         Qglobal, XBOPTS, Cao)


    #store initial matrices for eigenvalues analysis
    XBOUT.M0 = MglobalFull.copy()
    XBOUT.C0 = CglobalFull.copy()
    XBOUT.K0 = KglobalFull.copy()
       
    #"Get force vector for unconstrained nodes (Force_Dof)"
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    
    
    #"Get RHS at initial condition"
    Qglobal += -np.dot(FglobalFull, Force_Dof)
    
    
    #Separate assembly of follower and dead loads"   
    tmpForceTime=ForceTime[0].copy('F')  
    Qforces = XbeamLib.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                                    PosIni, PsiIni, PosDefor, PsiDefor, \
                                    ### sm increased rank of ForceDyn_*
                                    #(XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
                                    #(XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
                                    (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll[0,:,:]), \
                                    (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead[0,:,:]), \
                                    Cao, 1)[0]
                         
    Qglobal -= Qforces

    #"Initial Accel"
    dXddt[:] = np.dot(np.linalg.inv(MglobalFull), -Qglobal)
    #embed()
    #dXddt[iiblock]=0.0
    
    #"Record position of all grid points in global FoR at initial time step"
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    #"Position/rotation of the selected node in initial deformed configuration"
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    
    #"Get gamma and beta for Newmark scheme"
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*np.power((gamma + 0.5),2.0)
    
    
    #"Time loop"
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('  Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        
        #"calculate dt"
        dt = Time[iStep+1] - Time[iStep]
        
        
        #"Update transformation matrix for given angular velocity"
        Temp = np.linalg.inv(Unit4 + 0.25*XbeamLib.QuadSkew(ForcedVel[iStep+1,3:])*dt)
        Quat = np.dot(Temp, np.dot(Unit4 - 0.25*XbeamLib.QuadSkew(ForcedVel[iStep,3:])*dt, Quat))
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = XbeamLib.Rot(Quat)
        
        
        #"Predictor step"
        X += dt*dXdt + (0.5-beta)*dXddt*dt**2
        dXdt += (1.0-gamma)*dXddt*dt
        dXddt[:] = 0.0
        
        # Force calculation: moved in while as depending on PsiDefor
        #Force=(XBINPUT.ForceStatic+XBINPUT.ForceDyn*ForceTime[iStep+1]).copy('F')
        #Force=(XBINPUT.ForceStatic+XBINPUT.ForceDyn[iStep+1,:,:]).copy('F')


        #"Reset convergence parameters"
        Iter = 0
        ResLog10 = 1.0
        
        
        #"Newton-Raphson loop"        
        while ( (ResLog10 > XBOPTS.MinDelta.value) \
                & (Iter < XBOPTS.MaxIterations.value) ):
            
            #"set tensors to zero"
            Qglobal[:] = 0.0 
            Mvel[:,:] = 0.0
            Cvel[:,:] = 0.0
            MglobalFull[:,:] = 0.0
            CglobalFull[:,:] = 0.0
            KglobalFull[:,:] = 0.0
            FglobalFull[:,:] = 0.0
            
            #"Update counter"
            Iter += 1
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('   %-7d ' %(Iter))
                
            
            #"nodal diplacements and velocities from state vector"
            BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
    
            # update force

            Force = (XBINPUT.ForceStatic+XBINPUT.ForceDyn[iStep+1,:,:]).copy('F')
            AddGravityLoads(Force, XBINPUT, XBELEM,
                        AELAOPTS=None, # allows defining inertial/elastic axis
                        PsiDefor=PsiDefor,
                        chord = 0.0, # used to define inertial/elastic axis
                        PsiA_G=XbeamLib.quat2psi(Quat),
                        FollowerForceRig=XBOPTS.FollowerForceRig.value)  

            #"update matrices"
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         Force, ForcedVel[iStep+1,:], ForcedVelDot[iStep+1,:],\
                         NumDof, Settings.DimMat,\
                         ms, MglobalFull, Mvel,\
                         cs, CglobalFull, Cvel,\
                         ks, KglobalFull, fs, FglobalFull,\
                         Qglobal, XBOPTS, Cao)
            
            #"Get force vector for unconstrained nodes (Force_Dof)"
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
            

            Q1 =  np.dot(MglobalFull, dXddt)
            Q2 = np.dot(Mvel,ForcedVelDot[iStep+1,:])
            Q3 = - np.dot(FglobalFull, Force_Dof)

            #"Solve for update vector"
            #"Residual"
            Qglobal += np.dot(MglobalFull, dXddt) \
                     + np.dot(Mvel,ForcedVelDot[iStep+1,:]) \
                     - np.dot(FglobalFull, Force_Dof)
                              
            #Separate assembly of follower and dead loads  
            tmpForceTime=ForceTime[iStep+1].copy('F') 
            Qforces = XbeamLib.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                                            PosIni, PsiIni, PosDefor, PsiDefor, \
                                            (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
                                            (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
                                            Cao,1)[0]
                           
            Qglobal -= Qforces
            
            if XBOPTS.PrintInfo.value==True:                 
                sys.stdout.write('%-10.4e ' %(max(abs(Qglobal))))
            
            
            #"Calculate system matrix for update calculation"
            Asys = KglobalFull + \
                      CglobalFull*gamma/(beta*dt) + \
                      MglobalFull/(beta*np.power(dt,2.0))
            

            #"Solve for update"
            DX[:] = np.dot(np.linalg.inv(Asys), -Qglobal)
            
            #"Corrector step"
            X += DX
            dXdt += DX*gamma/(beta*dt)
            dXddt += DX/(beta*dt**2)


            #"Residual at first iteration"
            if(Iter == 1):
                Res0_Qglobal = max(max(abs(Qglobal)),1)
                Res0_DeltaX  = max(max(abs(DX)),1)
            
            #"Update residual and compute log10"
            Res_Qglobal = max(abs(Qglobal))
            Res_DeltaX  = max(abs(DX))
            
            ResLog10 = max(Res_Qglobal/Res0_Qglobal,Res_DeltaX/Res0_DeltaX)
            
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('%-10.4e %8.4f\n' %(max(abs(DX)),ResLog10))
                
        #"END Netwon-Raphson"
        
        
        #"update to converged nodal displacements and velocities"
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
        
        
        PosPsiTime[iStep+1,:3] = PosDefor[(NumNodes_tot.value-1)/2+1,:]
        PosPsiTime[iStep+1,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        #"Position of all grid points in global FoR"
        i1 = (iStep+1)*NumNodes_tot.value
        i2 = (iStep+2)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
    
        # sm I/O: FoR A velocities/accelerations
        XBOUT.QuatList.append(Quat)
        XBOUT.DynOut=DynOut
        XBOUT.PsiList.append(PsiDefor.copy())   
        XBOUT.cputime.append( time.clock() - XBOUT.cputime[0] )

        if SaveDict['SaveProgress'] and SaveDict['Format']=='h5':
            iisave=np.arange(1,NumSteps.value,np.ceil(
                                      NumSteps.value/SaveDict['NumSavePoints']))
            if any(iisave==iStep):
                PyLibs.io.save.h5file(SaveDict['OutputDir'], 
                                      SaveDict['OutputFileRoot']+'.h5',*OutList)
    
    #"END Time loop"


    if SaveDict['Format'] == 'dat': 
        #"Write _dyn file"
        ofile = Settings.OutputDir + Settings.OutputFileRoot+'_SOL312_dyn.dat'
        fp = open(ofile,'w')
        BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
        fp.close()
        #"Write _shape file"
        ofile = Settings.OutputDir+Settings.OutputFileRoot+'_SOL312_shape.dat'
        fp = open(ofile,'w')
        BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
        fp.close()
    else:
        XBINPUT.ForceDyn = XBINPUT.ForceDyn + XBINPUT.ForceDyn_foll +  XBINPUT.ForceDyn_dead
        PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5', 
                              *OutList)

    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')


    return XBOUT
        
    


if __name__ == '__main__':
    """Main"""
    
    """Set up Xbopts for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1."""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 312 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-05
    XBOPTS.FollowerForce.value = False
    XBOPTS.PrintInfo.value = True
    XBOPTS.MaxIterations.value = 100
         
    """Set up Xbinput for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1."""
    XBINPUT = DerivedTypes.Xbinput(2,32)
#     XBINPUT.BeamLength = 5.0
#     XBINPUT.BeamStiffness[0,0] = 4.8e+08
#     XBINPUT.BeamStiffness[1,1] = 3.231e+08
#     XBINPUT.BeamStiffness[2,2] = 3.231e+08
#     XBINPUT.BeamStiffness[3,3] = 1.0e+06
#     XBINPUT.BeamStiffness[4,4] = 9.346e+06
#     XBINPUT.BeamStiffness[5,5] = 9.346e+06
#     XBINPUT.BeamMass[0,0] = 100
#     XBINPUT.BeamMass[1,1] = 100
#     XBINPUT.BeamMass[2,2] = 100
#     XBINPUT.BeamMass[3,3] = 10
#     XBINPUT.BeamMass[4,4] = 0.001 #Must be non-zero
#     XBINPUT.BeamMass[5,5] = 0.001 #Must be non-zero
#     XBINPUT.ForceStatic[-1,2] = 6.0e+05
#     XBINPUT.tfin = 0.1
#     XBINPUT.dt = 0.001
#     
#     XBINPUT.ForceDyn[-1,2] = 7.0e+05
    
    #Solve_F90(XBINPUT,XBOPTS)
    #Solve_F90_steps(XBINPUT,XBOPTS)
    Solve_Py(XBINPUT,XBOPTS,'Input_DANC')