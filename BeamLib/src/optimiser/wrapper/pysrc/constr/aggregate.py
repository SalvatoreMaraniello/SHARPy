'''
Created on 8 Jun 2015

@author: sm6110

Module containing aggregation techniques 

'''

import numpy as np


  
    
def KSsum(gvec,k=3.0,alpha=1.0,**kwargs):
    ''' 
    Given a a vector gvec, containing positive and negative numbers, the
    function return an approximation for the maximum value of the vector. 
    
    For discrete sums (alpha=1), the KS functional is ALWAYS conservative.
    
    alpha: scaling parameter. 
    k: as k-> inf, the KS -> max(gvec).
    
    key words argoments:
    - m: can be arbitrary and does not, mathematically, influence the result. It
    is used, however, for numerical stability. By default, it is set as the max.
    value in gvec to avoid overflow during the exponential calculations.
    '''
    
    if (kwargs.has_key('m')):
        m=kwargs['m']
    else:
        m=np.max(gvec)
        
    
    
    if k>1e3:
        raise NameError('k value in KS function too high!!! The problem may be numerically unfeasible: the simulation is stopped.')
    
    gk = (gvec-m)*k
    
    return m + 1.0/k * np.log( np.sum(np.exp(gk)) / alpha )






if __name__=='__main__':
    
    print 'Test KSeval: '
    gvec = -np.array([5., 40., -7., 10., -8., -10.0])
    print 'input values:', gvec
    
    print 'KS functional for different k values:'
    mval=np.max(gvec)
    a=1.0
    print 'k=0.1:', KSsum(gvec,k=0.1,m=mval,alpha=a)
    print 'k=1:'  , KSsum(gvec,k=1.0,m=mval,alpha=a)
    print 'k=3:' , KSsum(gvec,k=3.0,m=mval,alpha=a)
    print 'k=10:', KSsum(gvec,k=10.0,m=mval,alpha=a)
    
    