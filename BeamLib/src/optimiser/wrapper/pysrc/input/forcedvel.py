'''
Created on 17 Sep 2014

@author: sm6110


The module contains methods to define the velocity at the beam support.

The methods are equivalent to those implemented in the subroutine:

    module: input
    input_forcedvel (NumNodes,Time,ForcedVel,ForcedVelDot)
    
of the original Fortran beam solver. The output variables are:

    ForcedVel   (NumSteps+1,6) ! Forced velocities at the support.
    ForcedVelDot(NumSteps+1,6) ! Derivatives of the forced velocities at the support.

these describe velocities and accelerations of the support at NumSteps, evenly 
spaced, times-steps.

'''

import numpy as np

#-------------------------------------------------- Methods from input_forcedved

def set_zero(Time):
    # VelAmp=0.d0
    
    NumSteps = len(Time)-1
    
    ForcedVel    = np.zeros( (NumSteps+1,6), dtype='float', order='F' )
    ForcedVelDot = np.zeros( (NumSteps+1,6), dtype='float', order='F' )
    
    return ForcedVel, ForcedVelDot 




def set_sin(Time, Omega, TrAmpl, RotAmpl ):
    ''' 
    Sinusoidal oscillation of frequency Omega and amplitudes TrAmpl
    (translational dof) and MomAmpl (rotational dof), such that the final support
    motion will be:
    
    Tr/Rot = {TrAmpl/RotAmpl}*sin(Omega Time)
    
    '''
 
    NumSteps = len(Time)-1    
    
    ForcedVel    = np.zeros( (NumSteps+1,6), dtype='float', order='F' )
    ForcedVelDot = np.zeros( (NumSteps+1,6), dtype='float', order='F' )

    Time = Omega*Time

    ForcedVel[:,:3] = np.sin(Time).reshape((NumSteps+1,1)) * TrAmpl
    ForcedVel[:,3:] = np.sin(Time).reshape((NumSteps+1,1)) * RotAmpl 
       
    ForcedVelDot[:,:3]= Omega * np.cos(Time).reshape((NumSteps+1,1)) * TrAmpl
    ForcedVelDot[:,3:]= Omega * np.cos(Time).reshape((NumSteps+1,1)) * RotAmpl
        
    return ForcedVel, ForcedVelDot 



def set_ramp( Time, Tfull, TrAmp, RotAmp ):
    '''
    Tfull: time-step at which the full load amplitude is reached.
    '''
 
    if Tfull > Time[-1]:
        raise NameError('Tfull has to be less-equal then Time[-1]') 
    
    NumSteps = len(Time)-1

    ForcedVel    = np.empty( (NumSteps+1,6), dtype='float', order='F' )
    ForcedVelDot = np.empty( (NumSteps+1,6), dtype='float', order='F' )

    vel_shape = np.ones( (NumSteps+1,1), dtype=float, order='F' )
    acc_shape = np.zeros( (NumSteps+1,1), dtype=float, order='F' )
    
    NSfull = np.sum(Time < Tfull)                # 1st time-step with full load
    for ii in range(NSfull):                     # NSfull is excluded
        vel_shape[ii]= Time[ii]/Tfull
        acc_shape[ii]=      1.0/Tfull
    
    ForcedVel[:,:3] = vel_shape * TrAmp
    ForcedVel[:,3:] = vel_shape * RotAmp
    
    ForcedVelDot[:,:3] = acc_shape * TrAmp
    ForcedVelDot[:,3:] = acc_shape * RotAmp 
    
    return ForcedVel, ForcedVelDot


def set_rampsin( Time, Tfull, Omega, TrAmpl, RotAmpl ):
    '''
    Tfull: time-step at which the full load amplitude is reached.
    omega: sine frequency
    '''
    
    v1 = np.ones(3)
    
    vel_sin, acc_sin = set_sin(Time,Omega, TrAmpl, RotAmpl )
    vel_ramp1, acc_ramp1 = set_ramp( Time, Tfull, v1, v1 )
    
    ForcedVel = vel_sin * vel_ramp1
    ForcedVelDot = acc_sin*vel_ramp1 + vel_sin*acc_ramp1
    
    return ForcedVel, ForcedVelDot    





if __name__=='__main__':
    
    import matplotlib.pyplot as plt
        
    Time = np.linspace(0,5,500) # sec
    Omega = 2.0*np.pi*2.0 # 2Hz 
    TrAmpl =np.array([1.,10.,100.])
    RotAmpl=np.array([2.,20.,200.]) 
    Tfull = 3.0
    

    #ForcedVel, ForcedVelDot = set_sin(Time, Omega, TrAmpl, RotAmpl )
    #ForcedVel, ForcedVelDot = set_ramp(Time, Tfull, TrAmpl, RotAmpl )
    ForcedVel, ForcedVelDot = set_rampsin(Time, Tfull, Omega, TrAmpl, RotAmpl )
    
    
    for ii in range(3):
        plt.plot(Time,ForcedVel[:,ii],'r')
        plt.plot(Time,ForcedVel[:,ii+3],'k')
        plt.title('Velocities of Tr and Rot no. %d' %(ii))
        plt.show()
        
        plt.plot(Time,ForcedVelDot[:,ii],'r')
        plt.plot(Time,ForcedVelDot[:,ii+3],'k')
        plt.title('Accelerations of Tr and Rot no. %d' %(ii))
        plt.show()       
    
    
        
        
        





