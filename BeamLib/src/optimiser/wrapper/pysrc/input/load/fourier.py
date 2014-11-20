'''
Created on 3 Nov 2014

@author: sm6110

Global reconstruction of signal using Fourier coefficients

'''

import numpy as np



def glsin_nodal(NumNodes, Time, NodeForce, A, F ):
    '''
    Given a set of amplitudes coefficient and related frequencies, the function
    reconstructs a signal over the time-domain Time creates overlapping 
    sinusoidal waves. The signal is used as a loading component acting on 
    the node of number 'NodeForce'. 
    
    The input A and F are arrays of size (cfmax,6), where cfmax is the total
    number of frequencies used. If A,F are 1D arrays, they are reshaped, if possible, 
    into a cfmax,6 dimension array. 
      
    Input:
    NumNodes: total number of nodes in the structure.
    A, F: For each of the ii-th component of the force at the node NodeForce, 
    the array A(:,ii) returns the amplitude associated to the sinusoidal waves of
    frequency F(:,ii)
    nodeForce: node where to apply the force. The nodes numbering is such that:
            root: NodeForce =  0
            tip:  NodeForce = -1 or NumNodes-1
    '''
    
    NumSteps = len(Time)-1
    if NodeForce > NumNodes-1: raise NameError("NodeForce can't be higher then NumNodes-1!")
    
    Fdyn = np.zeros( (NumNodes,6,NumSteps+1), dtype=float, order='F')
    
    # reshape if 1D array (allows A,F to be 1D) 
    Acp = reshape_cfs(A)    
    Fcp = reshape_cfs(F)
    
    # build F-time matrix
    for comp in range(6):
        Fdyn[NodeForce,comp,:]=sin_series_vec(Time,Acp[:,comp],Fcp[:,comp])
    
    return Fdyn



def sin_series_vec(Time, av, fv):
    '''
    Given a time vector Time, and a set of coefficients av and fv, the function
    builds the sines series signal:
    
        y(ii) = sin( 2*pi*Time[ii]*fv[ii] )
    '''

    NumSteps = len(Time)-1
    cfmax = len(fv)
    if cfmax != len(av): raise NameError('av and fv should have same length!')
    
    y = np.zeros(NumSteps+1)

    for jj in range(cfmax):
        y = y + av[jj]*np.sin(2.0*np.pi*fv[jj]*Time)
    
    return y


def reshape_cfs(A):
    ''' Given a 1D array of coefficients for the Fourier/Sin representation of
    loads, the function reshapes them into a (cfmax,6) array '''
    # reshape if 1D array (allows A,F to be 1D) 

    if len(A.shape) == 1:
        print 'reshaping!!!'
        B=np.array(A,order='F')
        cfmax,rem=len(B)/6,len(B)%6
        if rem!=0:
            raise NameError('If A,F are 1D arrays, their length has to be a multiple of 6')
        A=B.reshape((cfmax,6))
    return A


if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    
    NodeForce=0
    NumNodes=3
    cfmax=4
    
    A, F = np.zeros((cfmax,6)), np.ones((cfmax,6))
    Time=np.linspace(0,2,100)
    
    # build frequencies
    fo = 0.5/Time[-1]
    fv = np.zeros(cfmax)
    for ii in range(cfmax):
        fv[ii]=(ii+1)*fo
        F[ii,:]=np.ones(6)*fv[ii]    
    print 'frequencies [Hz]: ', fv
    print 'F matrix [hz]', F
    
    A[0,0]=1.0 # Fx - single freq
    A[1,1]=1.0 # Fy - single freq
    A[2,2]=1.0 # Fz - single freq
    A[3,3]=1.0 # Mx - single freq
    A[:,4]=np.array([0.2, 0.2, 0.2, 0.4])
    A[:,5]=np.array([0.8, 0.1, 0.1, 0.1])
    print 'A matrix: ', A
    
    # check reshape capability if A is 1D
    A2 = np.zeros((4*6,))
    A2[0]=1.0 # Fx - single freq
    A2[7]=1.0 # Fy - single freq
    A2[14]=1.0 # Fz - single freq
    A2[21]=1.0 # Mx - single freq
    A2[4],A2[10],A2[16],A2[22]=0.2,0.2,0.2,0.4
    A2[5],A2[11],A2[17],A2[23]=0.8,0.1,0.1,0.1
    
    
    Fdyn = glsin_nodal(NumNodes, Time, NodeForce, A2, F )
    
    
    print '----------------------------------- Check against sin_series_vec ---'
    
    for comp in range(6):
        y = sin_series_vec(Time,A[:,comp],F[:,comp])
        fig=plt.figure('check sin_series_vec')
        ax = fig.add_subplot(111)
        hman,=ax.plot(Time,y,'r+')
        hauto,=ax.plot(Time,Fdyn[NodeForce,comp,:],'k')
        plt.legend([hman, hauto],('manual','auto'))
        #fig.savefig('fig_sin_series_vec_%d.pdf' %(comp))
        plt.show()
    
    # check that other nodes are not allocated    
    for ii in range(NumNodes):
        M = np.max(np.max( np.abs(Fdyn[ii,:,:]) ))
        print 'Max abs value for node %d is %f' %(ii,M)
    
    
    
    
    
    
    
        
        
    
    
    
        
    
    
    