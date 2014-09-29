'''
Created on 19 Sep 2014

@author: sm6110


PostProcessing tools

'''

import numpy as np


def reshape_DynOut(DynOut,NumSteps):
    ''' 
    Reshape output form dynamic simulation (xbcomponent.DynOut) into a 3D array
    format having:
        THPos[NS,:,:] = PosDef at time-step NS
    where:
        Pos[nn,ii] = Position of ii-th coordinate of nn-th node
    note that:
        THPos[0,:,:] = PosIni
        THPos[-1,:,:]= PosDef
    '''
    
    NumNodes = (DynOut.shape[0])/(NumSteps+1)
    THPos = DynOut.reshape((NumSteps+1,NumNodes,3))
    
    #remark:
    # this gives the wrong output
    #THPos = DynOut.reshape((NumSteps+1,NumNodes,3),order='F')

    return THPos



def compute_velocity(THPos,Time):
    '''
    Compute velocities at each node from the position array THPos (output from
    reshape_DynOut) using discrete finite differences
    '''
    
    shTH = THPos.shape
    NumSteps = shTH[0]-1
    NumNodes = shTH[1]
    
    inv_2dt = 1.0/(Time[2]-Time[0])
    
    THVel = np.zeros((NumSteps,NumNodes,3), dtype=float, order='F')
    
    # timestep 1: backward 
    THVel[0,:,:] = (THPos[1,:,:]-THPos[0,:,:]) / (Time[1]-Time[0])
    
    # timestep 2... NumSteps-1
    for jj in range(2,NumSteps-1,1):
        THVel[jj,:,:] = inv_2dt * ( THPos[jj+1,:,:] - THPos[jj-1,:,:] )
    
    # timestep NumSteps: fwd
    THVel[-1,:,:] = (THPos[-1,:,:]-THPos[-2,:,:]) / (Time[-1]-Time[-2])
    
    return THVel



def compute_envelope(THres):
    '''
    Defines Envelope of THres, where Thres is the output from the routines
        reshape_DynOut (or)
        compute_velocities 
    '''
     
    
    EnvMin = np.min(THres,axis=0)
    EnvMax = np.max(THres,axis=0)
    
    return EnvMin, EnvMax
        
    
    
    
    
    
    