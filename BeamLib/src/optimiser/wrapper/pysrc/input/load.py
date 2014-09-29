'''
Created on 17 Sep 2014

@author: sm6110


The module contains methods to define static and dynamic loads by determining:

    a. the force distribution on the beam (set_*_Force functions)    
    b. the time variation (set_*_ForceTime functions)


Remark:

The methods are equivalent to those implemented in the subroutine:

    module: input
    input_dynforce (NumNodes,Time,ForceStatic,ForceDynAmp,ForceTime)
    
of the original Fortran beam solver. The methods here implemented can be used
to allocate the following variables of the original Fortran code:

    ForceStatic (NUmNodes,6) ! Internal force/moments at nodes.
    ForceDynAmp (NumNodes,6) ! Amplitude of the applied dynamic nodal forces.
    ForceTime   (NumSteps+1) ! Time history of the dynamic nodal forces.
    
Note that the original Fortran code only allows for loads of the kind:

    Load = ForceDynAmp(position) * ForceTime(time)

'''

import numpy as np



''' ----------------------------------------------------------- Dictionaries '''
# not implemented
# consider the idea of implementing a dictionary for the functions and a 
# a method that, given all the input, creates the specific input for each mehtod



''' ---------------------------------------------------- Forces distribution '''


def set_nodal_ForceDistr(NumNodes,NodeForce,
          ForceVec  = np.zeros((3)) ,
          MomentVec = np.zeros((3)) ):
    '''
    Applies a load at one node of the model
    
    This method is useful to define root/tip applied forces:
        root: NodeForce =  0
        tip:  NodeForce = -1
    '''
    
    ForceDistr = np.zeros( (NumNodes,6), dtype='float', order='F')

    ForceDistr[NodeForce,:3]=ForceVec
    ForceDistr[NodeForce,3:]=MomentVec
    
    ForceDistr = np.array( ForceDistr, dtype='float', order='F')
        
    return ForceDistr
    

def set_unif_ForceDistr(NumNodes,
          ForceVec  = np.zeros((3)) ,
          MomentVec = np.zeros((3)) ):
    '''
    Applies the same load to all the nodes of the model. The loads can be scaled
    during the simulation using one of the set_***_ForceTime methods but the
    "shape" of the loads will remain the same
    ''' 
    
    ForceDistr = np.zeros( (NumNodes,6), dtype='float', order='F')

    for ii in range(NumNodes):
        ForceDistr[ii,:3]=ForceVec
        ForceDistr[ii,3:]=MomentVec
        
    return ForceDistr




''' -------------------------------------------------------------- ForceTime '''

def set_const_ForceTime( Time ):
    '''    
    Remark: this is not a step load as at the initial time the beam is already 
    deformed (i.e. the load is not zero)!
    '''
    
    NumSteps = len(Time)-1
    ForceTime = np.ones( (NumSteps+1), dtype='float', order='F' )
        
    return ForceTime
    
    
def set_impulse_ForceTime( Time ):
    '''
    "initial_load in Fortran code"
    
    If ForceVec and MomentVec are not passed, a null load is generated.
    
    Remark: this is not a step load as at the initial time the beam is already 
    deformed (i.e. the load is not zero)!
    '''  
    
    NumSteps = len(Time)-1
    ForceTime = np.zeros( (NumSteps+1), dtype='float', order='F' )
    ForceTime[0]=1.0
        
    return ForceTime   


def set_ramp_ForceTime( Time, Tfull ):
    '''
    Tfull: time-step at which the full load amplitude is reached.
    '''
    
    NumSteps = len(Time)-1
    ForceTime = np.ones( (NumSteps+1), dtype='float', order='F' )
    
    if Tfull > Time[-1]:
        raise NameError('Tfull has to be less-equal then Time[-1]') 
    
    NSfull = np.sum(Time < Tfull) # 1st time-step with full load
    
    for ii in range(NSfull): # NSfull is excluded
        ForceTime[ii]= Time[ii]/Tfull
    
    return ForceTime


def set_sin_ForceTime(Time, Omega):
    
    ForceTime = np.sin( Omega * np.array(Time, dtype='float', order='F') )
    
    return ForceTime


def set_cos_ForceTime(Time, Omega):
    
    ForceTime = np.cos( Omega * np.array(Time, dtype='float', order='F') )
    
    return ForceTime


def set_omcos_ForceTime(Time, Omega):
    
    ForceTime = np.array(1.0 - set_cos_ForceTime(Time, Omega),
                         dtype='float', order='F')
    
    return ForceTime
    

def set_rampsin_ForceTime( Time, Tfull, Omega ):
    '''
    Tfull: time-step at which the full load amplitude is reached.
    omega: sine frequency
    '''
    
    ForceTime= set_ramp_ForceTime( Time, Tfull) * set_sin_ForceTime(Time,Omega)
                
    return ForceTime    



if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    
    Time = np.linspace(0,3,400)
    NumNodes = 3
    F = np.array([1,2,3])
    M = np.array([4,5,6]) 
    Omega = 2.0 * (2.0*np.pi)/Time[-1] # 2 periods in plot
    Tfull = 1.2
    
    print 'Test set_unif_ForceDistr'
    print set_unif_ForceDistr(NumNodes,F,M)
    print set_unif_ForceDistr(NumNodes,F)
    print set_unif_ForceDistr(NumNodes,MomentVec=M)
    print set_unif_ForceDistr(NumNodes)
    
    print 'Test set_nodal_ForceDistr'
    print 'Tip'
    print set_nodal_ForceDistr(NumNodes,-1,F,M)
    print 'Root'
    print set_nodal_ForceDistr(NumNodes, 0,F,M)
    
    print '--------------------------------------------------------------------'
    print 'Time'
    print Time
    
    print 'Test set_const_ForceTime'
    print set_const_ForceTime(Time)
    
    print 'Test set_impulse_ForceTime'
    print set_impulse_ForceTime(Time)    
    
    print 'Test set_ramp_ForceTime'
    FTr = set_ramp_ForceTime(Time,Tfull)
    print FTr
    plt.plot(Time,FTr,'b')
    plt.plot(Time,-FTr,'b')

    print 'Test set_rampsin_ForceTime (and set_sin_ForceTime)'
    FTsr=set_rampsin_ForceTime( Time, Tfull, Omega )
    print FTsr
    plt.plot(Time,FTsr,'k')
    plt.plot(Tfull*np.ones((3)),np.linspace(-1,1,3),'ro-')
    plt.title('Sine Ramp with full load at %f' %(Tfull))
    plt.show()
    
    print 'Test set_omcos_ForceTime (and set_cos_ForceTime)'
    FTomc=set_omcos_ForceTime( Time, Omega )
    print FTomc
    plt.plot(Time,FTomc,'k')
    plt.title('1 - Cosine')
    plt.show()
    
    