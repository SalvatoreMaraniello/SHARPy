'''
Created on 10 Mar 2016

@author: sm6110

Includes common methods for treating constraints

'''

import numpy as np


def KSsum(gvec,k=3.0,alpha=1.0,**kwargs):
    ''' 
    Given a a vector gvec, containing positive and negative numbers, the
    function return an approximation for the maximum value of the vector. 
    
    For discrete sums (alpha=1), the KS functional is ALWAYS conservative.
    For time integration, alpha = 1/dt
    
    alpha: scaling parameter. 
    k: as k-> inf, the KS -> max(gvec).
    
    key words arguments:
    - m: can be arbitrary and does not, mathematically, influence the result. It
    is used, however, for numerical stability. By default, it is set as the max.
    value in gvec to avoid overflow during the exponential calculations.
    
    '''
    
    assert k<1e3, ('k value in KS function too high!!! The problem may be numerically '+
                   'infeasible: the simulation is stopped.')

    if 'm' in kwargs: m=kwargs['m']
    else:  m=np.max(gvec)
            
    gk = (gvec-m)*k
    
    return m + 1.0/k * np.log( np.sum(np.exp(gk)) / alpha )



def transcription(gv,dt,eps):
    '''
    
    Given the constraint gv(t,x)>=0
    The function creates a smooth approximation of fmin = min(gv,0) and returns
    
        int[fmin] dt
    
    Ref. The Control Parameterization Method for Nonlinear Optimal Control: A Survey
    
    @warning: method not verified

    '''
    
    Ng = len(gv)
    fmin = np.zeros((Ng,))
    for ii in range(Ng):
        gg = gv[ii]
        if   gg <-eps: fmin[ii]=gg[ii]
        elif gg > eps: fmin[ii]=0.0
        else:  fmin[ii] = -(gg - eps)**2/(4.0*eps)
    
    return np.trapz(fmin)*dt
        
    
    






if __name__=='__main__':
    
    print('Test KSeval: ')
    gvec = -np.array([5., 40., -7., 10., -8., -10.0])
    gvec = -np.array([0.1, 0.2, 0.22, 0.2, 0.18])
    gvec =  np.array([0.1, 0.2, 0.22, 0.2, 0.18])
    gvec =  0*np.pi/180.0*np.array([0.1, 0.02, 0.0, -0.08, -0.18])
    print('input values:', gvec)
    
    print('KS functional for different k values:')
    mval=np.max(gvec)
    a=1.0
    print('k=0.1:', KSsum(gvec,k=0.1,m=mval,alpha=a) )
    print('k=1:'  , KSsum(gvec,k=1.0,m=mval,alpha=a) )
    print('k=3:' , KSsum(gvec,k=3.0,m=mval,alpha=a)  )
    print('k=10:', KSsum(gvec,k=10.0,m=mval,alpha=a) )
    print('k=100:', KSsum(gvec,k=100.0,m=mval,alpha=a) )
    print('k=999:', KSsum(gvec,k=999.0,m=mval,alpha=a) )
    
    