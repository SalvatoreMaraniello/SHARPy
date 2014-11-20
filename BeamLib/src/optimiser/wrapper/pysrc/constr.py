'''
Created on 31 Oct 2014

@author: sm6110

Methods for constraints

'''



import numpy as np
import scipy.integrate
import ctypes as ct
import sys
sys.path.append('..')

import shared
import lib.diff

# Load Dynamic Library
xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)



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



if __name__=='__main__':
    
    A=2.3
    Time=np.linspace(0,2.0*np.pi,2000)
    Force = A*np.sin(Time)
    
    J = qiqi_force_reg(Force,Time)
    Jexp = 0.5*A**2*(np.pi + 1e-2*np.pi)
    print 'Jnum: %f' %(J)
    print 'Jexp: %f' %(Jexp)