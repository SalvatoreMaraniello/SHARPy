'''
Created on 17 Sep 2014

@author: sm6110


The module contains methods to define .

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

def zero(Time):
    # VelAmp=0.d0
    
    NumSteps = len(Time)-1
    
    ForcedVel    = np.zeros( (NumSteps+1,6), dtype='float', order='F' )
    ForcedVelDot = np.zeros( (NumSteps+1,6), dtype='float', order='F' )
    
    return ForcedVel, ForcedVelDot 


