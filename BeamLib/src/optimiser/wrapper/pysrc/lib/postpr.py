'''
Created on 19 Sep 2014

@author: sm6110


PostProcessing tools.


Ref.:
- integration:
    http://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html

'''

import numpy as np
import scipy.integrate

import sys
sys.path.append('../..')
import shared
import PyBeam.Utils.XbeamLib


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



def compute_velocity(THPos,Time):
    '''
    Compute velocities at each node from the position array THPos (output from
    reshape_DynOut) using discrete finite differences
    '''
    
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



def compute_envelope(THres):
    '''
    Defines Envelope of THres, where Thres is the output from the routines
        reshape_DynOut (or)
        compute_velocities 
    '''
     
    
    EnvMin = np.min(THres,axis=0)
    EnvMax = np.max(THres,axis=0)
    
    return EnvMin, EnvMax
        
    

def integrate_function(F,x,method='trap',dF=0.0):
    
    '''
    Given the matrix F such that: 
        F[:,jj] = [ f_jj(x) ]
    i.e. such that each column is a vector valued function over the array x, 
    the method returns the matrix of integral function
        I[:,jj] = [i_jj(x)]
    
    Methods: 
    - 1tay: I order Taylor forward
    - trap: II order, very accurate for zero sums
    - 2tay: II order Taylor, with Forward evaluated derivatives
    - 2taycnt: II order Taylor with I derivative computed as average (as trap).
               the method is I order.
    - trap2tay: fusion of trap + II order term evaluated forward. method I and
                II order, depending on the function.
    
    Remarks:

    
    Warning:
    This method is not suitable to integrate rotations!!!
    '''
    
     
    # set-up variables
    sh = F.shape
    NumSteps=sh[0]-1
    I = np.zeros(sh,dtype=float,order='F')    
    dx = x[1]-x[0]
    dx2 = dx**2
    
    # ---------------------------------------------------------------- integrate
    # methods organised in order of accuracy
    
    if method=='1tay':
        # Taylor expansion I order with Forward formula
        # I[ii+1,:]=I[ii,:] + dx*F[ii,:]           
        # backward version also possible:
        # I[ii+1,:]=I[ii,:] + dx*F[ii+1,:]          
        for ii in range(NumSteps):
            I[ii+1,:]=I[ii,:] + dx*F[ii ,:]    
    
    elif method=='trap':
        # use trapezoidal rule
        # equivalent to I order Taylor expansion with average on the I derivative:
        # I[ii+1,:]=I[ii,:] + 0.5*dx*(F[ii,:]+F[ii+1,:])
        for jj in range(sh[1]):
            I[1:,jj]=scipy.integrate.cumtrapz(F[:,jj],x=x)
                    
    elif method=='2tay':
        # Taylor expansion II order with Forward formula
        # backward version also possible:
        # I[ii+1,:]=I[ii,:] + dx*F[ii+1,:] +dx2/2 * dF[ii+1,:] 
        for ii in range(NumSteps):
            I[ii+1,:]=I[ii,:] + dx*F[ii  ,:] + 0.5*dx2*dF[ii,:]        
           
    elif method=='2taycnt': 
        # Taylor II order expansion with derivatives evaluated as average between
        # current and forward time-step 
        for ii in range(NumSteps):
            dI  = 0.5*dx * (F[ii,:]+F[ii+1,:])
            ddI = 0.25*dx2 * (dF[ii,:]+dF[ii+1,:])
            I[ii+1,:]=I[ii,:] + dI + ddI
    
    elif method=='trap2tay':
        # fusion of trapezoidal rule (I derivative) plus II order term from 
        # taylor expansion evaluated forward
        # NO IMPROVEMENT
         for ii in range(NumSteps):
            dI  = 0.5*dx * (F[ii,:]+F[ii+1,:])
            ddI = 0.5*dx2 * (dF[ii,:])
            I[ii+1,:]=I[ii,:] + dI + ddI    
    
    elif method=='mytrap':
        # hand implementation of trapezoidal rule to check the correctness of
        # the formula
        for ii in range(NumSteps):
            dI  = 0.5*dx * (F[ii,:]+F[ii+1,:])
            I[ii+1,:]=I[ii,:] + dI               
        
    else:
        raise NameError('Use a valid integration method!')
                
    return I



def integrate_rotations(Omega, Time, xi0=np.array([1,0,0,0]),method='1tay',dOmega=0.0):
    ''' 
    The routine integrates the rotations starting from the angular velocity 
    vector.
    
    Input:
    - Omega: matrix such that Omega[tt,:] is the angular velocity at the 
    time-step tt
    - time: time vector
    - xi0: initial quaternion describing the initial condition. If not passed in
    input, xi0 is automatically set to represent a zero rotation, i.e.:
        xi0 = ( cos(fi/2) ; sin(fi/2)*n  )^T    with fi=0 (n arbitrary)
        xi0 = (1 0 0 0)^T
     
    Methods: 
    - 1tay: I order Taylor forward
    - trap: II order, very accurate for zero sums
    - 2tay: II order Taylor, with Forward evaluated derivatives 
        
    Ref.:  Aircraft Control and Simulation, pag. 46, by Stevens, Lewis.
    
    Remarks: tests performed with integrate_function (see below) showed II order
        accuracy with trapezoidal rule. This does not require evaluation of 
        II order terms (i.e. dOmegadt) so has an easier implementation. 
     
    '''
    
    dt = Time[1]-Time[0]
    dt2 = dt**2
    NumSteps = len(Time)-1
    
    Xi = np.zeros((NumSteps+1,4))
    Xi[0,:]=xi0
    
    
    if method=='1tay':
        # Taylor expansion I order with Forward formula
        # I[ii+1,:]=I[ii,:] + dx*F[ii,:]           
        # backward version also possible:
        # I[ii+1,:]=I[ii,:] + dx*F[ii+1,:] 
        for tt in range(NumSteps):
            # build incremental update operator
            QuadSkew = PyBeam.Utils.XbeamLib.QuadSkew(Omega[tt,:])   
            # get quaternion derivative
            dxi = -0.5* np.dot(QuadSkew,Xi[tt,:])  
            # use simple integration
            Xi[tt+1,:] = Xi[tt,:] + dt*dxi
     
  
    elif method=='trap':
        # use trapezoidal rule
        # equivalent to I order Taylor expansion with average on the I derivative:
        # I[ii+1,:]=I[ii,:] + 0.5*dx*(F[ii,:]+F[ii+1,:])
        raise NameError('Method not implemented as required iterations! ' +
                        'Use Runge-Kutta or other explicit method!')
     
           
    elif method=='2tay':
        # Taylor expansion II order with Forward formula
        # backward version also possible:
        # I[ii+1,:]=I[ii,:] + dx*F[ii+1,:] +dx2/2 * dF[ii+1,:] 
        for tt in range(NumSteps):
            
            # build incremental update operator
            QuadSkew = PyBeam.Utils.XbeamLib.QuadSkew(Omega[tt,:])   
            # get quaternion derivative
            dxi = -0.5* np.dot(QuadSkew,Xi[tt,:])  
            
            # Build II order term
            dQuadSkew = PyBeam.Utils.XbeamLib.QuadSkew(dOmega[tt,:])
            ddxi = -0.5* ( np.dot(QuadSkew,dxi) +
                           np.dot(dQuadSkew,Xi[tt,:]) )
            
            # use simple integration
            Xi[tt+1,:] = Xi[tt,:] + dt*dxi + 0.5*dt2*ddxi  
           

    else:
        raise NameError('Use a valid integration method!')     
        
    return Xi
         
        
    
def vector_time_history(Xi,vec0=np.array([1,0,0])):
    ''' 
    Given a history of quaternions (see 'integrate_rotations'), and initial vector
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



def THPosDefGlobal(DynOut,RefVel,Time,set_origin='a'):
    ''' 
    Given a time simulation, the function changes the FoR in which the Position
    pf each node of the beam, at each time-step of the simulation, from local (a)
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
    
    '''

    NumSteps=len(Time)-1
    NumNodes = (DynOut.shape[0])/(NumSteps+1)
    THPosDefGlobal=np.zeros((NumSteps+1,NumNodes,3))  
     
    THPosDefLocal = DynOut.reshape((NumSteps+1,NumNodes,3)) # output in a frame

    if set_origin=='a':
        aOrigin = np.zeros((NumSteps+1,3)) # a bit memory consuming but who cares...
    elif set_origin=='G':
        aOrigin = integrate_function(RefVel[:,:3],Time)
    else:
        raise NameError("Set a valid origin ('a' or 'G') for the positions vectors!")
        
    # compute quaternions associated with a frame rotation
    Xi = integrate_rotations(RefVel[:,3:],Time)
    
    # apply rotation to the position vector of each node
    # nb: two loops are requires as the position vector changes at each time-step
    for tt in range(NumSteps+1):
        Cga = PyBeam.Utils.XbeamLib.Rot(Xi[tt,:]).transpose() 
        for nn in range(NumNodes):     
            THPosDefGlobal[tt,nn,:]=np.dot(Cga,THPosDefLocal[tt,nn,:]) +aOrigin[tt,:]
    
    return THPosDefGlobal
        
    
 
  

''' ------------------------------------------------------------------------ '''



if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    import lib.plot.dyn

    test='tran'#'tran' # 'rot' # 'rotconv'



    if test=='rotconv': 
        print '--------------- Check Convergence rate of Rotations Integration '
        print 'run testing.integration module'
        

    if test=='tran':
        print '---------------------------------------- Test functions Integration '
        
        expmax=6
        eevec=np.arange(expmax)[1:]
        print eevec
        method='1tay'
        # create matrices to store the error
        Er = np.zeros((expmax,3))
        
        AnalSol=np.array([2, 0, np.exp(np.pi)-1.0]) 
        
        print 'Analytical ', AnalSol
        for ee in eevec:

            # set time-steps and discretised funciton
            NT=10**ee
            Time = np.linspace(0,np.pi,NT)
            print NT
            F = np.zeros((NT,3))
            F[:,0]=np.sin(Time)
            F[:,1]=np.cos(Time)
            F[:,2]=np.exp(Time)
            
            dF = np.zeros((NT,3))
            dF[:,0] =np.cos(Time)
            dF[:,1] =-np.sin(Time)
            dF[:,2] =np.exp(Time)       
            
            RefPos = integrate_function(F,Time,method=method,dF=dF)
            print 'Computed: ', RefPos[-1,:]
            Er[ee,:] = RefPos[-1,:]-AnalSol            
            
        print Er
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
        print '---------------------------------------- Test Rotations integration '
        NT=10
        Time = np.linspace(0,2,NT)   
        
        print 'constant angular velocity about the y axis'
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
        
        Xi = integrate_rotations(Omega, Time)
        print 'Xi[-3:,:]'
        print Xi[-3:,:]
        
        Vx = vector_time_history(Xi)                   # rotation of x unit vector
        Vy = vector_time_history(Xi,np.array([0,1,0])) # should stay constant
        V45=vector_time_history(Xi,np.array([0,1,1]))  # rotation of non-unit 45 deg vector
        
        print 'Vx[-3:,:]'
        print Vx[-3:,:]
        print 'Vy[-3:,:]'
        print Vy[-3:,:]    
        print 'V45[-3:,:]'
        print V45[-3:,:]
        lib.plot.rig.axis_traj(Vx[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()
        lib.plot.rig.axis_traj(Vy[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()    
        lib.plot.rig.axis_traj(V45[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()
        
        print 'constant angular velocity about a general axis'
        # Even here, the axis of rotation is invariant (in the inertial FoR this is
        # obvious ad the inertial fFoR is not rotationg, in the body FoR also - imagine 
        # a rigid body rotating: relative positions do not change.
        # Even here we expect to return to the initial rotation, i.e.
        # Xi[-1,:] = -1+. 0+ 0+ 0], with components 1 and 2 being identical (as
        # associated to the 45 deg axis n) 
        nrot = np.array( [np.cos(np.pi/4.0), np.sin(np.pi/4.0), 0] )
        wvec = np.pi * nrot
        Omega = np.ones((NT,3))*wvec
        
        Xi = integrate_rotations(Omega, Time)
        print 'Xi[-3:,:]'
        print Xi[-3:,:]    
        
        Vx = vector_time_history(Xi)                   # rotation of x unit vector
        Vy = vector_time_history(Xi,np.array([0,1,0])) # should stay constant
        V45= vector_time_history(Xi,np.array([0,1,1])) # rotation of non-unit 45 deg vector
        print 'Vx[-3:,:]'
        print Vx[-3:,:]
        print 'Vy[-3:,:]'
        print Vy[-3:,:]    
        print 'V45[-3:,:]'
        print V45[-3:,:]
        
        lib.plot.rig.axis_traj(Vx[:-5,:]) # not all cycle to visualise the direction of rotation
        plt.show()
        lib.plot.rig.axis_traj(Vy[:-5,:]) # not all cycle to visualise the direction of rotation   
        plt.show() 
        lib.plot.rig.axis_traj(V45[:-5,:]) # not all cycle to visualise the direction of rotation 
        plt.show()     
    
