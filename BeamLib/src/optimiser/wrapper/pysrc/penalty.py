'''
Created on 31 Oct 2014

@author: sm6110

Methods for constraints

'''

import numpy as np
import scipy.integrate
###import ctypes as ct
import sys
sys.path.append('..')

###import shared
import lib.diff

# Load Dynamic Library
###xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)


def qiqi_force_reg(Force,Time):
    ''' 
    Constraint "Regularisation" term to limit the applied force based on 
    Qiqi & Yu (2012)
    '''
    
    # differentiate force
    dF = lib.diff.difffun(Force,Time)
    
    # build term to integrate
    I = Force**2 + 1e-2 * dF**2
    
    # and integrate
    J = 0.5*scipy.integrate.simps(I, Time) 
    
    return J



def quad_push(pval,pmax,plim,alpha):
    ''' 
    Push term to add to cost function to speed up convergence. 
    - pval is the value of a constraint that should be pval<pmax. 
    - plim is the limit after which the gradient dp remains constant to alpha
    
    The value returned is: 
    - push=0          if f < pmax
    - push'=alpha     if f > plim
    - push quadratic  if pmax < f < plim
     '''
    
    df = pval - pmax
    
    if pval<=pmax:
        push = 0.0
    elif pval>=plim:
        push = 0.5*alpha*(plim-pmax) + alpha*(pval-plim)
    else:
        push = 0.5*alpha/(plim-pmax) * (pval-pmax)**2
    
    return push




if __name__=='__main__':
    
    '''
    A=2.3
    Time=np.linspace(0,2.0*np.pi,2000)
    Force = A*np.sin(Time)
    
    J = qiqi_force_reg(Force,Time)
    Jexp = 0.5*A**2*(np.pi + 1e-2*np.pi)
    print 'Jnum: %f' %(J)
    print 'Jexp: %f' %(Jexp)
    '''
    
    
    import matplotlib.pyplot as plt
    
    fv = np.linspace(0,10,1000)
    
    alpha = 1.0
    pmax=6.0
    plim=6.5
    
    pushv = np.array([ quad_push(ff,pmax,plim,alpha) for ff in fv])
    dpv = np.diff(pushv)/(fv[1]-fv[0])
    
    plt.plot(fv,pushv)
    plt.show()
    
    plt.plot(fv[1:],dpv)
    plt.show()
    
    