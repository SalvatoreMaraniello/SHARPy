'''
# Salvatore Maraniello 5 Aug 2014

This module contains methods to allocate the variables for static and dynamic problems.

The assumption of 1 beam only model is used.

The pre-allocation is required to add robustness and (possibly) to speed-up the code.

'''

import numpy as np
import ctypes as ct

# Parameters
MaxElNod = 3

def total_nodes(NumElems,NumNodesElem):
    ''' Determines the total number of nodes in the model.
        A 1 beam model is assumed.

        Input:
        Numelems: number of elements in the model
        NumNodesElem [ctypes.c_int]: nodes per elem (2: linear, 3: quadratic)
        
        Output:
        NumNodes: total number of nodes in the element
    '''

    if (NumNodesElem is 2) or (NumNodesElem is 3):
        NumNodes =   (NumNodesElem-1)*NumElems + 1
    else:
        raise NameError('Wrong Input: NumNodesElem must be equal to 2 or 3. Current value: ' %(NumNodesElem))

    return NumNodes


def fwd_static(NumNodes,NumElems):
    ''' Allocates output variables for forward execution of static problem 

        Input:
        NumNodes: total number of nodes in the element
        NumElems: total number of elements
        
        Output:
        PosIni, PosDef: initial and deformed nodal position
        PsiIni, PsiDef: initial and deformend nodal rotation  
        InternalForces: internal forces      
    '''

    PosIni        =np.empty((         NumNodes,3),dtype=float,order='F')    
    PosDef        =np.empty((         NumNodes,3),dtype=float,order='F')    
    PsiIni        =np.empty((NumElems,MaxElNod,3),dtype=float,order='F')    
    PsiDef        =np.empty((NumElems,MaxElNod,3),dtype=float,order='F')  
    InternalForces=np.empty((         NumNodes,6),dtype=float,order='F')   

    return PosIni, PosDef, PsiIni, PsiDef, InternalForces

