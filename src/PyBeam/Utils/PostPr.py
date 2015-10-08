'''
Created on 19 Sep 2014

@author: sm6110

PostProcessing tools.

'''

import ctypes as ct
import numpy as np
from warnings import warn
import copy

import sys
sys.path.append('../..')

import PyBeam.Utils.XbeamLib
import PyLibs.numerics.diff, PyLibs.numerics.integr
import BeamIO, BeamInit
import PyBeam.Utils.DerivedTypes




def reshape_DynOut(DynOut,NumSteps):
    ''' 
    Reshape output form dynamic simulation (xbcomponent.DynOut) into a 3D array
    format having:
        THPos[NS,:,:] = PosDef at time-step NS
    where:
        Pos[nn,ii] = Position of ii-th coordinate of nn-th node
    note that:
        THPos[0,:,:] = PosIni
        THPos[-1,:,:]= PosDef
    '''
    
    NumNodes = (DynOut.shape[0])/(NumSteps+1)
    THPos = DynOut.reshape((NumSteps+1,NumNodes,3))
    
    #remark:
    # this gives the wrong output
    #THPos = DynOut.reshape((NumSteps+1,NumNodes,3),order='F')

    return THPos


def compute_velocity_Iord(THPos,Time):
    '''
    Compute velocities at each node from the position array THPos (output from
    reshape_DynOut) using discrete finite differences
    '''
    
    warn('Function First Order Accuracy!')
    
    shTH = THPos.shape
    NumSteps = shTH[0]-1
    NumNodes = shTH[1]
    
    inv_2dt = 1.0/(Time[2]-Time[0])
    
    THVel = np.zeros((NumSteps,NumNodes,3), dtype=float, order='F')
    
    # timestep 1: backward 
    THVel[0,:,:] = (THPos[1,:,:]-THPos[0,:,:]) / (Time[1]-Time[0])
    
    # timestep 2... NumSteps-1
    for jj in range(2,NumSteps-1,1):
        THVel[jj,:,:] = inv_2dt * ( THPos[jj+1,:,:] - THPos[jj-1,:,:] )
    
    # timestep NumSteps: fwd
    THVel[-1,:,:] = (THPos[-1,:,:]-THPos[-2,:,:]) / (Time[-1]-Time[-2])
    
    return THVel


def compute_velocity(THPos,Time):
    '''
    Compute velocities at each node from the position array THPos.
    
    The function assumes equally spaced array and by default uses a II order
    approximation. 

    Remark: the method for derivatives is based on polynomial interpolation 
    (diff.dwgt). It is thus general and derivatives of higher accuracy can be 
    obtained.
    
    The code is, however, specifically set for II order accuracy 1st derivatives.
    To improve it:
        a. the boundary treatment should be modified.
        b. sarr variable is now specific for II order accuracy
    '''
    

    shTH = THPos.shape
    NumSteps = shTH[0]-1
    NumNodes = shTH[1]
    
    THVel = np.zeros((NumSteps+1,NumNodes,3), dtype=float, order='F')
    
    # get weights for II order polynomial derivative
    # remark: if a general order of derivation is required, special attention
    # should be put in evaluating the derivatives at the boundary 
    kder=2
    KK=kder+1
    w0   = PyLibs.numerics.diff.dwgt( Time[:KK],  Time[0], k=1)
    wend = PyLibs.numerics.diff.dwgt(Time[-KK:], Time[-1], k=1)
    wmid = PyLibs.numerics.diff.dwgt( Time[:KK], Time[ 1], k=1)
    
    # timestep 1 and final:
    for jj in range(3):
        for nn in range(NumNodes):
            THVel[ 0,nn,jj] = np.dot(  w0,THPos[:KK,nn,jj])
            THVel[-1,nn,jj] = np.dot(wend,THPos[-KK:,nn,jj])
            
    # timestep 2... NumSteps-1
    # sarr specific for II order derivative
    sarr = np.array([-1,0,1]) 
    for tt in range(1,NumSteps,1):
        ttind = tt+sarr
        for jj in range(3):
            for nn in range(NumNodes):
                THVel[tt,nn,jj] = np.dot(wmid,THPos[ttind,nn,jj])
     
    return THVel


def compute_envelope(THres):
    '''
    Defines Envelope of THres, where Thres is the output from the routines
        reshape_DynOut (or)
        compute_velocities 
    '''
     
    
    EnvMin = np.min(THres,axis=0)
    EnvMax = np.max(THres,axis=0)
    
    return EnvMin, EnvMax
           
    
def vector_time_history(Xi,vec0=np.array([1,0,0])):
    ''' 
    Given a history of quaternions (see 'lib.integr.rotations'), and initial vector
    vec0, the function computes at each time-step tt the rotation matrix associated
    with the quaternion Xi[tt,:] and rotates vec0.
    The resulting vectors are stored into the matrix V
    
    Note that, if the rotation at tt=0 is non-zero, V[0,:] ~= vec0
    '''  
    
    NumSteps = len(Xi[:,0])-1
    V = np.zeros( (NumSteps+1,3), dtype=float, order='F' ) 
    
    for tt in range(NumSteps+1):
        
        Cga = PyBeam.Utils.XbeamLib.Rot(Xi[tt,:]).transpose()
        V[tt,:] = np.dot(Cga,vec0)
        
    return V


def THPosDefGlobal(DynOut,Time,RefVel,set_origin='a',**kwargs): #xi0=np.array([1,0,0,0]),
    ''' 
    Given a time simulation, the function changes the FoR in which the Position
    of each node of the beam, at each time-step of the simulation, from local (a)
    (default output in xbeam code - PosDef variable) to ground (G).
    
    According to the value of origin the output will be:
    set_origin='a': the position vector is unchanged, i.e. its origin is in 
        the local (a) FoR.
    set_origin='G': the origin of the position vector is set in the ground 
        (G) FoR.

    The output array PosDefGlobal has shape:
        THPosDefGlobal[NS,nn,ii] 
    with:
        NS: time-step
        nn: beam node
        ii: coordinate x, y, z
    
    The orientation of the FoR A can be passed in input (Xi) or can be integrated
    from the FoR A velocities.
    
    The RefVel is intended in G FoR!!!
    '''

    # quaternion key words
    if 'Xi' in  kwargs:
        Xi=kwargs['Xi']
    else:
        if 'xi0' in kwargs:
            xi0 = kwargs['xi0']
        else:
            xi0=np.array([1,0,0,0])
        Xi = PyLibs.numerics.integr.rotations(RefVel[:,3:],Time,xi0=xi0)
        
    NumSteps=len(Time)-1
    NumNodes = int( (DynOut.shape[0])/(NumSteps+1) )
    THPosDefGlobal=np.zeros((NumSteps+1,NumNodes,3))  
     
    THPosDefLocal = DynOut.reshape((NumSteps+1,NumNodes,3)) # output in a frame

    if set_origin=='a':
        aOrigin = np.zeros((NumSteps+1,3))
    elif set_origin=='G':
        aOrigin = PyLibs.numerics.integr.function(RefVel[:,:3],Time)
    else:
        raise NameError("Set a valid origin ('a' or 'G') for the positions vectors!")
        
    #### compute quaternions associated with a frame rotation
    #if not 'Xi' in kwargs:
    #    Xi = PyLibs.numerics.integr.rotations(RefVel[:,3:],Time,xi0=xi0)
    
    # apply rotation to the position vector of each node
    # nb: two loops are requires as the position vector changes at each time-step
    for tt in range(NumSteps+1):
        Cga = PyBeam.Utils.XbeamLib.Rot(Xi[tt,:]).transpose() 
        for nn in range(NumNodes):     
            THPosDefGlobal[tt,nn,:]=np.dot(Cga,THPosDefLocal[tt,nn,:]) +aOrigin[tt,:]
    
    return THPosDefGlobal
        
  
def ElasticForces(XBINPUT,XBOPTS,DynOut,PsiList): 
    '''
    Given the results of a simulation computes the stresses along
    the beam at each time step.
    
    - PsiList is a list of CRV at each step of the simulation. 
    - DynOut is an array containing all the displacements.
    - XBINPUT & XBOPTS are either:
        1. objects of the Xbinput/Xbopts classes (see PyBeam.Utils.DerivedTypes)
           with attributes defined as c_types
        2. objects from a class having the same attributes (but not defined 
           as c_types)
    
    The elastic forces are returned in the format:
        S = [ time_step, Node, component ]
        
           
    Remark:
    - PsiList and DynOut are equivalent (one for displ one for rotation)
      but the format is inconsistent!
    '''
    
    # Prepare classes
    if not(isinstance(XBOPTS,PyBeam.Utils.DerivedTypes.Xbopts)):
        H=copy.deepcopy(XBOPTS)
        XBOPTS=PyBeam.Utils.DerivedTypes.Xbopts()
        XBOPTS.init_from_class(H)
    
    # Initialise beam
    XBINPUT, XBOPTS, NumNodes_tot, \
    XBELEM, PosIni, PsiIni, XBNODE, NumDof = BeamInit.Static(XBINPUT,XBOPTS)
    
    # Extract Elements Type
    if isinstance(XBINPUT,PyBeam.Utils.DerivedTypes.Xbinput):
        NumNodesElem=XBINPUT.NumNodesElem.value 
    else:
        NumNodesElem=XBINPUT.NumNodesElem
    
    # Prepare Node List
    if NumNodesElem == 2:
        NodeList=[ii for ii in range(int(NumNodes_tot.value))] 
    elif NumNodesElem == 3:
        # midpoint nodes not supported
        NodeList=[2*ii for ii in range(int(0.5*NumNodes_tot.value))]
    else:
        raise NameError('Incorrect Number of Nodes per Element in XBINPUT!')
    
    # get number of steps
    NumSteps = int( DynOut.shape[0]/NumNodes_tot.value ) - 1
    
    print('dynout shape', DynOut.shape       )
    print('NumNodes'    , NumNodes_tot.value )
    print('NumSteps'    , NumSteps           )
    
    PosDef=PyBeam.Utils.PostPr.reshape_DynOut(DynOut,NumSteps)
    
    # compute stresses
    Smat=np.empty(( NumSteps, len(NodeList), 6 ))
    for tt in range(NumSteps):
        # strains=BeamIO.localStrains(PosDef[tt,:,:], xbout.PsiList[tt], 
        #                             xbout.PosIni, xbout.PsiIni, XBELEM, NodeList,SO3 = False)
        Smat[tt,:,:]=BeamIO.localElasticForces(PosDef[tt,:,:], PsiList[tt], 
                                                 PosIni,  PsiIni, 
                                                 XBELEM, NodeList)
    
    return Smat, NodeList
        

''' ------------------------------------------------------------------------ '''



if __name__=='__main__':
    
    import matplotlib.pyplot as plt

    test='tran'#'tran' # 'rot' # 'rotconv'


    if test=='rotconv': 
        print('Check Convergence rate of Rotations Integration')
        print('run testing.integration module')
        

    if test=='tran':
        print('Test functions Integration')
        
        expmax=6
        eevec=np.arange(expmax)[1:]

        method='1tay'
        # create matrices to store the error
        Er = np.zeros((expmax,3))
        
        AnalSol=np.array([2, 0, np.exp(np.pi)-1.0]) 
        
        print('Analytical ', AnalSol)
        for ee in eevec:

            # set time-steps and discretised funciton
            NT=10**ee
            Time = np.linspace(0,np.pi,NT)
            
            F = np.zeros((NT,3))
            F[:,0]=np.sin(Time)
            F[:,1]=np.cos(Time)
            F[:,2]=np.exp(Time)
            
            dF = np.zeros((NT,3))
            dF[:,0] =np.cos(Time)
            dF[:,1] =-np.sin(Time)
            dF[:,2] =np.exp(Time)       
            
            RefPos = PyLibs.numerics.integr.function(F,Time,method=method,dF=dF)
            print('Computed: ', RefPos[-1,:])
            Er[ee,:] = RefPos[-1,:]-AnalSol            
            
        AbsEr=abs(Er)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(np.arange(expmax),AbsEr[:,0],'r')
        ax.plot(np.arange(expmax),AbsEr[:,1],'b')
        ax.plot(np.arange(expmax),AbsEr[:,2],'k')
        ax.set_yscale('log')
        plt.legend(('sin','cos','exp'))
        plt.title('Convergence Rate for "%s" method' %(method))
        fig.savefig('convrate_'+method+'.png')
        plt.show()
               
        
        ## uncomment to show result
        #plt.plot(Time,RefPos1[:,0],'k')
        #plt.plot(Time,RefPos1[:,1],'b')
        #plt.show()
        
        
    if test=='rot': 
        print('Test Rotations integration ')
        NT=10
        Time = np.linspace(0,2,NT)   
        
        print('constant angular velocity about the y axis')
        # in this test, the angular velocity is only about the y axis.
        # It doesn't matter if the y axis is that of the earth or body FoR: this
        # is invariant
        # The velocity is constant pi rad/s, so in 2sec the FoR is expected to 
        # perform a full loop.
        # The final quaternion is expected to be:
        # xi = [cos(2pi/2) sin(2pi/2)*j] = [-1+ 0 0- 0]
        # where -1+ is -1 approaching from up
        # 0+ is 0 approaching from up
        Omega = np.zeros((NT,3))
        Omega[:,1]=np.ones(NT)*np.pi # pi rad/sec
        
        Xi = PyLibs.numerics.integr.rotations(Omega, Time)
        print('Xi[-3:,:]',Xi[-3:,:])
        
        Vx = vector_time_history(Xi)                   # rotation of x unit vector
        Vy = vector_time_history(Xi,np.array([0,1,0])) # should stay constant
        V45=vector_time_history(Xi,np.array([0,1,1]))  # rotation of non-unit 45 deg vector
        
        print('Vx[-3:,:]', Vx[-3:,:])
        print('Vy[-3:,:]', Vy[-3:,:]) 
        print('V45[-3:,:]', V45[-3:,:])
        ###lib.plot.rig.axis_traj(Vx[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()
        ###lib.plot.rig.axis_traj(Vy[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()    
        ###lib.plot.rig.axis_traj(V45[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()
        
        print('constant angular velocity about a general axis')
        # Even here, the axis of rotation is invariant (in the inertial FoR this is
        # obvious ad the inertial fFoR is not rotationg, in the body FoR also - imagine 
        # a rigid body rotating: relative positions do not change.
        # Even here we expect to return to the initial rotation, i.e.
        # Xi[-1,:] = -1+. 0+ 0+ 0], with components 1 and 2 being identical (as
        # associated to the 45 deg axis n) 
        nrot = np.array( [np.cos(np.pi/4.0), np.sin(np.pi/4.0), 0] )
        wvec = np.pi * nrot
        Omega = np.ones((NT,3))*wvec
        
        Xi = PyLibs.numerics.integr.rotations(Omega, Time)
        print('Xi[-3:,:]',Xi[-3:,:])   
        
        Vx = vector_time_history(Xi)                   # rotation of x unit vector
        Vy = vector_time_history(Xi,np.array([0,1,0])) # should stay constant
        V45= vector_time_history(Xi,np.array([0,1,1])) # rotation of non-unit 45 deg vector
        print('Vx[-3:,:]', Vx[-3:,:])
        print('Vy[-3:,:]', Vy[-3:,:]) 
        print('V45[-3:,:]', V45[-3:,:])
        
        ###lib.plot.rig.axis_traj(Vx[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()
        ###lib.plot.rig.axis_traj(Vy[:-5,:]) # not all cycle to visualise the direction of rotation   
        plt.show() 
        ###lib.plot.rig.axis_traj(V45[:-5,:]) # not all cycle to visualise the direction of rotation 
        plt.show()     
    
