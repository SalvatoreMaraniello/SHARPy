'''
Created on 13 Apr 2016

@author: sm6110

Methods for smoothing

'''

import scipy.linalg
import numpy as np


def laplacian(tv,yv,alpha):
    '''
    Smooth a signal based on Sobolev gradient method:
        p(t) - alpha ddp/dt2 = P
    where p and P are, respectively, smothed and original signal
    
    The method can be used to smooth continuous gradients ot
    time histories.
    
    The solution is sought solving:
    
        (I - alpha L) p = P
    
    where L and I are, respectively, the laplacian and identity matrices,
    plus boundary conditions 
        p[0]=P[0]
        p[T]=P[T]
    
    Note: this is similar to the heat equation:
        dp/dt - alpha ddp/dx2 = f
    where typically one time-stepping for p_new after time discretisation. If the
    time stepping is implicit, the same algorithm is retrieved.
    
    Warning:  that the matrix ddp/dt2 is ill conditioned for small time-steps.
    Having p - alpha Laplacian reduces the condition number, but special care
    should be used if alpha becomes large.

    '''
  
    assert len(tv)==len(yv), 'Size of tv and yv not matching!'

    dt=tv[1]-tv[0]
    cflap=-alpha/dt**2
    cvec = np.zeros((len(yv),))
    cvec[0]=-2.0*cflap+1.0
    cvec[-1]=1.0*cflap
    cvec[ 1]=1.0*cflap
    D = scipy.linalg.circulant(cvec)
    
    # BC enforcement
    D[0,:], D[-1,:]=0.0, 0.0
    D[0,0], D[-1,-1]=1.0, 1.0
    # smooth signal
    ysmooth = np.linalg.solve(D,yv)    

    return ysmooth
    
    
def avg(yv,cicles=1):
    '''
    Smooth a function yv using: 
        yv_filt[ii] = 0.5* (yv[ii-1] + yv[ii+1])
    The process is iterated cicles times
    
    Warning: use laplacian method
    '''
    
    N=len(yv)
    yfilt=yv.copy()
    for ii in range(cicles):
        for ii in range(1,N-1):
            yfilt[ii]=0.5*(yv[ii-1]+yv[ii+1]) 
        yv=yfilt.copy() 
    return yfilt