'''@package PyBeam.Solver.NonlinearStatic
@brief      Nonlinear static solvers.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       25/10/2012
@pre        None
@warning    None
'''

import sys
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import BeamLib
import XbeamLib
import BeamInit
import numpy as np
import ctypes as ct
import PyLibs.io.save
from XbeamLib import AddGravityLoads, Rot, psi2quat
from IPython import embed

def Solve_F90(XBINPUT,XBOPTS,SaveDict=Settings.SaveDict):
    """@brief Nonlinear static structural solver using f90 solve routine."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 112, ('NonlinearStatic (F90) requested' +\
                                              ' with wrong solution code')
    
    # I/O management
    Settings.SaveDict=SaveDict # overwrite for compatibility with 'dat' format output methods
    XBOUT=DerivedTypes.Xboutput()
    OutList=[XBOPTS, XBINPUT, XBOUT]  

    "Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    
    "Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    # Add Gravity
    ForceStaticCopy = XBINPUT.ForceStatic.copy('F')
    AddGravityLoads(XBINPUT.ForceStatic,XBINPUT,XBELEM,
                    AELAOPTS=None,# allows defining inertial/elastic axis
                    PsiDefor=PsiDefor,
                    chord=0.0,# used to define inertial/elastic axis
                    PsiA_G=np.zeros((3,)),# no G/A frame in this solution
                    FollowerForceRig=XBOPTS.FollowerForceRig.value)


    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear static case (using .f90 routines) ... \n')
    
    BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                                PosIni, PsiIni, XBNODE, NumDof,\
                                PosDefor, PsiDefor)
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')
    
    # Divide gravity and applied load for I/O
    #ForceTotal = XBINPUT.ForceStatic.copy('F')
    #XBINPUT.ForceStatic = ForceStaticCopy

    # Prepare output
    XBOUT.drop( PosIni=PosIni,PsiIni=PsiIni,PosDeforStatic=PosDefor.copy(), 
                PsiDeforStatic=PsiDefor.copy(),
                ForceStaticTotal=XBINPUT.ForceStatic.copy() )
    XBINPUT.ForceStatic=ForceStaticCopy

    #"Write deformed configuration to file"
    if SaveDict['Format']=='dat':
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL112_def.dat'
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Writing file %s ... ' %(ofile))
        fp = open(ofile,'w')
        fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
        fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
        fp.close()
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('done\n')
        WriteMode = 'a'
    
        BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                       PosDefor, PsiDefor, ofile, WriteMode)
    else:
        PyLibs.io.save.h5file( SaveDict['OutputDir'], 
                               SaveDict['OutputFileRoot']+'.h5',
                               *OutList)

    
    "Print deformed configuration"
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('--------------------------------------\n')
        sys.stdout.write('NONLINEAR STATIC SOLUTION\n')
        sys.stdout.write('%10s %10s %10s\n' %('X','Y','Z'))
        for inodi in range(NumNodes_tot.value):
            sys.stdout.write(' ')
            for inodj in range(3):
                sys.stdout.write('%12.5e' %(PosDefor[inodi,inodj]))
            sys.stdout.write('\n')
        sys.stdout.write('--------------------------------------\n')


    "Return solution as optional output argument"
    return XBOUT #PosDefor, PsiDefor


def Solve_F90_steps(XBINPUT,XBOPTS):
    """@brief Nonlinear static structural solver using f90 solve routine called
    once per load-step."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 112, ('NonlinearStatic (F90) requested' +\
                                              ' with wrong solution code')
    
    "Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    
    "Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear static case (using .f90 routine at' +\
                         ' each load-step) ... \n')
    
    
    
    "Initialise load increments and set F90 load-step to 1"
    LoadSteps = XBOPTS.NumLoadSteps.value
    LoadIncrement = XBINPUT.ForceStatic/LoadSteps
    XBOPTS.NumLoadSteps.value = 1
    
    
    "Start loading loop"
    for step in range(1,LoadSteps+1):
        "Current load to be applied"
        XBINPUT.ForceStatic = step*LoadIncrement
        
        "Print load step"
        if XBOPTS.PrintInfo.value == True:
            print('     Python-based outer load step %d' %(step))
          
        "Solve with one load step"
        BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                                PosIni, PsiIni, XBNODE, NumDof,\
                                PosDefor, PsiDefor)
        
        
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')
    
     
    "Write deformed configuration to file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL112_def.dat'
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Writing file %s ... ' %(ofile))
    fp = open(ofile,'w')
    fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
    fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
    fp.close()
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('done\n')
    WriteMode = 'a'
    
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                       PosDefor, PsiDefor, ofile, WriteMode)
    
    
    "Print deformed configuration"
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('--------------------------------------\n')
        sys.stdout.write('NONLINEAR STATIC SOLUTION\n')
        sys.stdout.write('%10s %10s %10s\n' %('X','Y','Z'))
        for inodi in range(NumNodes_tot.value):
            sys.stdout.write(' ')
            for inodj in range(3):
                sys.stdout.write('%12.5e' %(PosDefor[inodi,inodj]))
            sys.stdout.write('\n')
        sys.stdout.write('--------------------------------------\n')
        
    
    "Return solution as optional output argument"
    return PosDefor, PsiDefor


def Solve_Py(XBINPUT,XBOPTS, moduleName = None, SaveDict=Settings.SaveDict):
    """Nonlinear static solver using Python to solve residual
    equation. Assembly of matrices is carried out with Fortran subroutines."""
    
    #"Check correct solution code"
    assert XBOPTS.Solution.value == 112, ('NonlinearStatic requested' +\
                                              ' with wrong solution code')
    
    # I/O management
    Settings.SaveDict=SaveDict # overwrite for compatibility with 'dat' format output
    XBOUT=DerivedTypes.Xboutput()
    OutList=[XBOPTS, XBINPUT, XBOUT]

    #"Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS, moduleName)

    # debug
    XBOUT.xbnode=XBNODE
    XBOUT.xbelem=XBELEM
    
    #"Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')

    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear static case in Python ... \n')

    # Determine orientatin of FoR A w.r.t. G (for gravity loads)
    Quat = psi2quat(XBINPUT.PsiA_G)
    Cao = Rot(Quat)
    
    #"Initialise structural eqn tensors"
    KglobalFull = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); ks = ct.c_int()
    KglobalFull_foll = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); 
    FglobalFull = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); fs = ct.c_int()    
    
    DeltaS  = np.zeros(NumDof.value, ct.c_double, 'F')
    Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
    x       = np.zeros(NumDof.value, ct.c_double, 'F')
    dxdt    = np.zeros(NumDof.value, ct.c_double, 'F')
    
    # Add Gravity
    XBOUT.ForceStaticTotal = XBINPUT.ForceStatic.copy('F')
    AddGravityLoads(XBOUT.ForceStaticTotal,XBINPUT,XBELEM,
                    AELAOPTS=None, # allows defining inertial/elastic axis
                    PsiDefor=PsiDefor,
                    chord = 0.0, # used to define inertial/elastic axis
                    PsiA_G=XBINPUT.PsiA_G,
                    FollowerForceRig=XBOPTS.FollowerForceRig.value)

    #"Load Step tensors"
    iForceStep     = np.zeros((NumNodes_tot.value,6), ct.c_double, 'F')
    iForceStep_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    rpar=1.1 # amplification factor
    IncrPercent=np.array([ rpar**(ii) for ii in range(XBOPTS.NumLoadSteps.value) ])
    step00=1./np.sum(IncrPercent)
    IncrPercent=IncrPercent*step00
    assert np.sum(IncrPercent)-1.0<1e-10, 'Algorithm for determining step-'\
       'size of boundary displacements failed. Try reducing the number of steps'
        

    #"Start Load Loop"
    # NumLoadsSteps doubled if forced displacements are enforced in the solution:
    # - the first NumLoadsSteps increase the applied load
    # - the last NumLoadsSteps increase the applied displacements at the boundaries
    if len(XBINPUT.ForcedDisp)==0: 
        Nsteps = XBOPTS.NumLoadSteps.value
    else: 
        #NumDispSteps=XBOPTS.NumLoadSteps.value
        #Nsteps = NumDispSteps+XBOPTS.NumLoadSteps.value
        Nsteps = 2* XBOPTS.NumLoadSteps.value
        # At each iLoadStep>=XBOPTS.NumLoadSteps.value, an amount of 
        # displacement, determined through the IncrPercent factor, is enforced
        # at the boundary. IncrPercent can be arbitrary but should also be such
        # that the sum of each IncrPercent=1. 
        # The initial steps should be very small, so as to avoid the solution
        # to jump to a different branch.
        # Different models are possible:
        ## equally distributed: requires a lot of loads steps
        #IncrPercent=1./NumDispSteps*np.ones((NumDispSteps,))
        ## exponential distribution: guarantees very small initial steps but 
        # reduces the total number of steps
        #rpar=1.3 # amplification factor
        #IncrPercent=np.array([ rpar**(ii) for ii in range(NumDispSteps) ])
        #step00=1./np.sum(IncrPercent)
        #IncrPercent=IncrPercent*step00

    XBOUT.CondK=[]
    #for iLoadStep in range(XBOPTS.NumLoadSteps.value):
    for iLoadStep in range(Nsteps):
        
        #"Reset convergence parameters"
        Iter = 0
        ResLog10 = 1.0
        
        #"General load case"
        #iForceStep = XBOUT.ForceStaticTotal*float( (iLoadStep+1) / \
        #                                              XBOPTS.NumLoadSteps.value)
        #ForceScaling = np.min([ float((iLoadStep+1)/XBOPTS.NumLoadSteps.value),
        #                                                                   1.0])
        #iForceStep = XBOUT.ForceStaticTotal*ForceScaling

        # debug:
        #if iLoadStep>: 1/0
        
        if iLoadStep<XBOPTS.NumLoadSteps.value:
            # Gradually increase the static load
            IncrHere=IncrPercent[iLoadStep]
            iForceStep=iForceStep+IncrHere*XBOUT.ForceStaticTotal

            # prepare initial guess
            if iLoadStep<1:
                PosDefor_0 = PosIni.copy('F')
                PsiDefor_0 = PsiIni.copy('F')
                PosDefor_1 = PosIni.copy('F')
                PsiDefor_1 = PsiIni.copy('F')
            else:
                PosDefor_0 = PosDefor_1.copy('F')
                PsiDefor_0 = PsiDefor_1.copy('F')
                PosDefor_1 = PosDefor.copy('F')
                PsiDefor_1 = PsiDefor.copy('F')
            # guess for initial Newton iteration
            # damping
            damp=0.0#5#1#1.0
            PosDefor = PosDefor + damp*(PosDefor_1-PosDefor_0)*\
                                               IncrHere/IncrPercent[iLoadStep-1]
            PsiDefor = PsiDefor + damp*(PsiDefor_1-PsiDefor_0)*\
                                               IncrHere/IncrPercent[iLoadStep-1]
            # mid node
            #embed()
            midNode=int((XBINPUT.NumNodesTot-1)/2)
            print('R0\tR1\tdR\tRguess\tIncOld\tIncHere')
            print(6*'%.4f\t'\
                     %( PosDefor_0[midNode,2], 
                        PosDefor_1[midNode,2], 
                       (PosDefor_1-PosDefor_0)[midNode,2],
                        PosDefor[midNode,2],
                        IncrPercent[iLoadStep-1],
                        IncrHere))

        else: #if iLoadStep>=XBOPTS.NumLoadSteps.value:
            # Enforce displacement at the boundary
            IncrHere=IncrPercent[iLoadStep-XBOPTS.NumLoadSteps.value]
            # update displacements (also initial guess)
            print('IncrPercent (displacements are cumulative)=%f' %IncrHere)
            PosDefor = ForceDisplacement(PosIni,PosDefor,XBINPUT.ForcedDisp,
                                               XBINPUT.PsiA_G,factor=IncrHere)
        # debug forced displacements
        #import matplotlib.pyplot as plt 
        #plt.plot(PosIni[:,0],PosIni[:,2],'k',marker='s')
        #plt.plot(PosDefor[:,0],PosDefor[:,2],'r',marker='o')
        #plt.show()
        if XBOPTS.PrintInfo.value == True:
            sys.stdout.write('  iLoad: %-10d\n' %(iLoadStep+1))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')

        #"Newton Iteration"
        while( (ResLog10 > XBOPTS.MinDelta.value) \
             & (Iter < XBOPTS.MaxIterations.value) ):
            
            #"Increment iteration counter"
            Iter += 1
            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write('   %-7d ' %(Iter))

            #"Set structural eqn tensors to zero"
            KglobalFull[:,:] = 0.0; ks = ct.c_int()
            KglobalFull_foll[:,:] = 0.0; 
            FglobalFull[:,:] = 0.0; fs = ct.c_int()
            Qglobal[:] = 0.0

            #"Assemble matrices for static problem"
            BeamLib.Cbeam3_Asbly_Static(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                        PosIni, PsiIni, PosDefor, PsiDefor,\
                        iForceStep, NumDof,\
                        ks, KglobalFull, fs, FglobalFull, Qglobal,\
                        XBOPTS, Cao)

            #"Get forces on unconstrained nodes"
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              iForceStep.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              iForceStep_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
            
            #"Separate assembly of follower and dead loads"
            Qforces, KglobalFull_foll = \
                XbeamLib.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                                      PosIni, PsiIni, PosDefor, PsiDefor, \
                                      XBINPUT.ForceStatic_foll,XBINPUT.ForceStatic_dead, \
                                      Cao, float((iLoadStep+1)/XBOPTS.NumLoadSteps.value))[:2]
                        
            KglobalFull += KglobalFull_foll
            Qglobal = Qglobal -Qforces -np.dot(FglobalFull,iForceStep_Dof)

            #import matplotlib.pyplot as plt
            #plt.spy(KglobalFull)
            #embed()
            if iLoadStep==0 and Iter==1:
                XBOUT.K0=KglobalFull.copy()
                XBOUT.Q0=Qglobal.copy()
            #XBOUT.CondK[-1].append(np.linalg.cond(KglobalFull))

            #"Calculate \Delta State Vector"
            DeltaS = - np.linalg.solve(KglobalFull, Qglobal)
            # @warning: on very stiff problem, using a matrix inversion is more 
            # robust then using np.linalg.solve. A least-squares solver can also
            # be used.
            #if np.linalg.cond(KglobalFull)>1e20:
            #    Sol=np.linalg.lstsq(KglobalFull, -Qglobal)
            #    DeltaS=Sol[0]
            #    Krank=Sol[2]
            #    sing_val=Sol[3]
            #    if Krank<NumDof.value:
            #        raise NameError('Singular matrix')
            #else:
            #    #DeltaS = - np.dot(np.linalg.inv(KglobalFull), Qglobal)
            #    DeltaS = - np.linalg.solve(KglobalFull, Qglobal)

            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write(2*'%-10.4e '%(max(abs(Qglobal)),
                                                              max(abs(DeltaS))))
            #"Residual at first iteration"
            if(Iter == 1):
                Res0_Qglobal = max(max(abs(Qglobal)),1)
                Res0_DeltaX  = max(max(abs(DeltaS)),1)
            #"Update residual and compute log10"
            Res_Qglobal = max(abs(Qglobal))
            Res_DeltaX  = max(abs(DeltaS))

            ResLog10 = max(Res_Qglobal/Res0_Qglobal, Res_DeltaX/Res0_DeltaX)
            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write('%8.4f\n' %(ResLog10))

            # Update Solution (PosDefor/PsiDefor) for next iter
            BeamLib.Cbeam3_Solv_Update_Static(XBINPUT, NumNodes_tot, XBELEM,\
                                              XBNODE, NumDof, DeltaS,\
                                              PosIni, PsiIni, PosDefor,PsiDefor)

            #"Stop the solution"
            if(ResLog10 > 1.e10):
                sys.stderr.write(' STOP\n')
                sys.stderr.write(' The max residual is %e\n' %(ResLog10))
                return XBOUT
                exit(1)
            
        #"END Newton Loop"
    #"END Load Loop"

    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')
   
    # Prepare output
    XBOUT.drop( PosIni=PosIni,PsiIni=PsiIni,PosDeforStatic=PosDefor.copy(), 
                PsiDeforStatic=PsiDefor.copy(), K=KglobalFull.copy())

    #"Print deformed configuration"
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('--------------------------------------\n')
        sys.stdout.write('NONLINEAR STATIC SOLUTION\n')
        sys.stdout.write('%10s %10s %10s\n' %('X','Y','Z'))
        for inodi in range(NumNodes_tot.value):
            sys.stdout.write(' ')
            for inodj in range(3):
                sys.stdout.write('%12.5e' %(PosDefor[inodi,inodj]))
            sys.stdout.write('\n')
        sys.stdout.write('--------------------------------------\n')

    if SaveDict['Format']=='dat':     
        #"Write deformed configuration to file"
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL112_def.dat'
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Writing file %s ... ' %(ofile))
        fp = open(ofile,'w')
        fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
        fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
        fp.close()
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('done\n')
        WriteMode = 'a'
    
        BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                           PosDefor, PsiDefor, ofile, WriteMode)
    else:
        PyLibs.io.save.h5file(SaveDict['OutputDir'], SaveDict['OutputFileRoot']+'.h5',
                              *OutList)
    
    #"Return solution as optional output argument"
    return XBOUT #PosDefor, PsiDefor




def ForceDisplacement(PosIni,PosDefor,ForcedDisp,PsiA_G=np.zeros((3,)),
                                                                    factor=1.0):
    ''' Method to deform the Position of the grid so as to enforce 
    displacements at the boundaries as defined in ForcedDisp.
    The amount of  displacement is proportional to the distance between the 
    node and the node at which the forced displacement is enforced.
    @param factor: used to scale the displacements of the boundary nodes
    @warning: algorithm will may fail if node=0 and node=NumNodesTot are not
    the most distant nodes of the model.
    assume node=0 and node=-1 are the most distant nodes in the model
    '''
    NumNodesTot=PosIni.shape[0]
    DeltaMax=np.linalg.norm(PosIni[-1,:]-PosIni[0,:])

    for ff in range(len(ForcedDisp)):
        node_disp=ForcedDisp[ff]['node']

        # check FoR and convert position in FoR A coordinates
        if ForcedDisp[ff]['FoR']=='A':
            pos_disp=ForcedDisp[ff]['pos']
        elif ForcedDisp[ff]['FoR']=='G':
            Coa=XbeamLib.RotCRV(PsiA_G).transpose()
            pos_disp=np.dot(Coa,ForcedDisp[ff]['pos'])
        else:
            raise NameError('Wrong FoR in Xbinput.ForcedDisp!')

        # displacement enforced at node=node_disp
        DispEnf=factor*(pos_disp-PosIni[node_disp,:])
        # displace nodes
        for nn in range(NumNodesTot):
            # @note: delta needs to be computed on PosIni so as to ensure 
            # that at nodes 0 and -1 the scaling is equal to 0 or 1
            delta=PosIni[nn,:]-PosIni[node_disp,:]
            scaling = 1.-np.linalg.norm(delta)/DeltaMax
            #print('Node %d scaling %.3f'%(nn,scaling))
            PosDefor[nn,:]=PosDefor[nn,:]+scaling*DispEnf  

    
    # # debugging algorithm for deforming the beam when forced displacement are 
    # # enforced
    # PosIni=np.array([ 
    #     [0., 0., 0.],
    #     [3., 0., 0.],
    #     [6., 0., 0.],
    #     [6., 4., 0.],
    #     [6., 8., 0.],
    #     ])
    # XBINPUT.NumNodesTot=PosIni.shape[0]
    # XBINPUT.ForcedDisp=[]
    # XBINPUT.addForcedDisp(node=-1,  pos=np.array([4.,6.,0]) )
    # XBINPUT.addForcedDisp(node=0,  pos=np.array([-1.,2.,0]) )

    #PosDefor = ForceDisplacement(PosIni,PosDefor,XBINPUT.ForcedDisp,XBINPUT.PsiA_G,factor=1.0)
    # debug forced displacements
    #import matplotlib.pyplot as plt 
    #plt.plot(PosIni[:,0],PosIni[:,1],'k',marker='s')
    #plt.plot(PosDefor[:,0],PosDefor[:,1],'r',marker='o')
    #plt.show()

    return PosDefor  






if __name__ == '__main__':
    """Set up Xbopts for nonlinear static analysis defined in input_rob.f90
    TPY0 test case"""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 112 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-04
    XBOPTS.PrintInfo.value = True      
    """Set up Xbinput for nonlinear static analysis defined in input_rob.f90
    TPY0 test case"""
    XBINPUT = DerivedTypes.Xbinput(2,8)
    XBINPUT.BeamLength = 16.0
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 4.0e+06
    XBINPUT.BeamMass[0,0] = 0.75
    XBINPUT.BeamMass[1,1] = 0.75
    XBINPUT.BeamMass[2,2] = 0.75
    XBINPUT.BeamMass[3,3] = 0.1
    XBINPUT.BeamMass[4,4] = 0.001
    XBINPUT.BeamMass[5,5] = 0.001
    XBINPUT.ForceStatic[-1,2] = 800

    #Solve_F90(XBINPUT,XBOPTS)
    #Solve_F90_steps(XBINPUT,XBOPTS)
    Solve_Py(XBINPUT,XBOPTS)