'''
Created on 31 Oct 2014

@author: sm6110

Methods for constraints

'''



import numpy as np

import sys
sys.path.append('..')
import shared
import lib.postpr, lib.diff


# Load Dynamic Library
xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)



def qiqi_force_reg(Time,Force):
    ''' 
    Constraint "Regularization" term to limit the applied force based on 
    Qiqi & Yu (2012)
    '''
    
    # differentiate force
    dF = lib.diff.difffun(Time,Force)
    
    # build term to integrate
    I = Force**2 + 1e-2 * dF**2
    
    # and integrate
    J = lib.integr.function(I, Time, method='trap') 
    
    return J
    
    