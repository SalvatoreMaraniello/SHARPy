'''@package PyBeam.Utils.Input
@brief      Input functions to populate DerivedTypes for PyBeam simulations.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       28/11/2012
@pre        None
@warning    None

Modified: S. Maraniello, 28 Sep 2015
    - Node: New BCs allowed
'''

import sys
import DerivedTypes
import SharPySettings as Settings
import numpy as np
import ctypes as ct

def Setup(XBINPUT, XBOPTS):
    """@brief Check the beam input and options for inconsistencies.
    @details Initial functionality adapted from FxInput_PyFx.f90."""
    
    # Check number of nodes per element and Gauss points.
    if XBINPUT.NumNodesElem != 2 and XBINPUT.NumNodesElem != 3:
        sys.stderr.write('Invalid number of nodes per element\n')
    elif XBINPUT.NumNodesElem == 2 and XBOPTS.NumGauss.value != 1:
        sys.stderr.write('Number of nodes per element (2) inconsistent with' +\
                         ' number of gauss points (%d)'\
                         %(XBOPTS.NumGauss.value) +\
                         ' Gauss points set to (1)\n')
        XBOPTS.NumGauss.value = 1
    elif XBINPUT.NumNodesElem == 3 and XBOPTS.NumGauss.value != 2:
        sys.stderr.write('Number of nodes per element (3) inconsistent with' +\
                         ' number of gauss points(%d)'\
                         %(XBOPTS.NumGauss.value)  +\
                         ' Gauss points set to (2)\n')
        XBOPTS.NumGauss.value = 2    
    sys.stderr.flush()    
    return XBINPUT, XBOPTS


def Elem(XBINPUT, XBOPTS, XBELEM):
    """@brief Set up Xbelem derived type.
    @details Initial functionality adapted from 
    FxInput_PyFx.f90."""
    
    # Check number of nodes per element and set default connectivity.
    assert (XBINPUT.NumNodesElem == 2 or XBINPUT.NumNodesElem == 3),\
            'Invalid number of nodes per element'
    if XBINPUT.NumNodesElem == 2:
        for ElemNo in range(0,XBINPUT.NumElems):
            i = ElemNo*Settings.MaxElNod
            XBELEM.Conn[i]=ElemNo+1
            XBELEM.Conn[i+1]=ElemNo+2
            XBELEM.NumNodes[ElemNo] = 2
            
        NumNodes_tot = ct.c_int(XBINPUT.NumElems + 1) 
            
    elif XBINPUT.NumNodesElem == 3:
        for ElemNo in range(0,XBINPUT.NumElems):
            i = ElemNo*Settings.MaxElNod
            XBELEM.Conn[i]=2*ElemNo+1
            XBELEM.Conn[i+1]=2*ElemNo+3
            XBELEM.Conn[i+2]=2*ElemNo+2
            XBELEM.NumNodes[ElemNo] = 3
        
        NumNodes_tot = ct.c_int(2*XBINPUT.NumElems + 1)
    
    # Calculate Beam Inverse Stiffness Matrix.        
    InverseStiffness = np.linalg.inv(XBINPUT.BeamStiffness)
    
    # Populate XBELEM
    for ElemNo in range(0,XBINPUT.NumElems):
        i_start = ElemNo*6
        i_end = ElemNo*6+6
        XBELEM.Mass[i_start:i_end,:] = XBINPUT.BeamMass
        XBELEM.Stiff[i_start:i_end,:] = XBINPUT.BeamStiffness
        XBELEM.InvStiff[i_start:i_end,:] = InverseStiffness
    
    # Element orientation and member number.
    for ElemNo in range(0,XBINPUT.NumElems):
        i = ElemNo*3
        XBELEM.Vector[i+1] = 1.0
        XBELEM.MemNo[ElemNo] = 0
    
    return NumNodes_tot, XBELEM


def Node(XBINPUT, XBOPTS, NumNodes_tot, XBELEM):
    """@brief Define nodal properties.
    
    Note: the beam is always assumed to lie along the x axis. If spherical joint
    are found, this is shifted.
    
    """
    
    # Declare and populate boundary conditions.
    BoundConds = np.zeros(NumNodes_tot.value,dtype=ct.c_int,order='F')
    if XBINPUT.BConds[0] is 'M':
        BoundConds[0                     ] = -1
        BoundConds[NumNodes_tot.value - 1] = -1 
        MidNode = (NumNodes_tot.value - 1)/2
        if XBINPUT.BConds[1] is 'C': BoundConds[ MidNode ] = 1
        if XBINPUT.BConds[1] is 'S': BoundConds[ MidNode ] = 2
        x_shift=-XBINPUT.BeamLength*(float(MidNode)/float(NumNodes_tot.value - 1))
    else:    
        if XBINPUT.BConds[0] is 'C': BoundConds[                     0] =  1
        if XBINPUT.BConds[1] is 'C': BoundConds[NumNodes_tot.value - 1] =  1
        if XBINPUT.BConds[0] is 'F': BoundConds[                     0] = -1
        if XBINPUT.BConds[1] is 'F': BoundConds[NumNodes_tot.value - 1] = -1
        if XBINPUT.BConds[0] is 'S': BoundConds[                     0] =  2
        x_shift=0.0
        if XBINPUT.BConds[1] is 'S': 
            BoundConds[NumNodes_tot.value - 1] =  2
            x_shift=-XBINPUT.BeamLength 
    if np.max(BoundConds)==0 or XBINPUT.BConds=='SS':
        raise Exception('Invalid boundary conditions string.')

    
    # Default straight beam based on BeamLength
    PosIni = np.zeros((NumNodes_tot.value,3), dtype=ct.c_double, order='F')
    
    for NodeNo in range(0,NumNodes_tot.value):
        PosIni[NodeNo,0] = x_shift + XBINPUT.BeamLength*\
                            (float(NodeNo)/float(NumNodes_tot.value - 1))
        
    # Define pre-twist.
    PhiNodes = np.zeros(NumNodes_tot.value, dtype=ct.c_double, order='F')    

    return PosIni, PhiNodes, BoundConds
    

if __name__ == '__main__':
    XBINPUT = DerivedTypes.Xbinput()
    XBINPUT.NumNodesElem = 3
    XBOPTS = DerivedTypes.Xbopts()
    XBINPUT, XBOPTS = Setup(XBINPUT, XBOPTS)
    XBINPUT, XBOPTS = Setup(XBINPUT, XBOPTS)
    
    XBINPUT.NumNodesElem = 3
    XBINPUT.NumElems = 8
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 4.0e+06
    XBELEM = DerivedTypes.Xbelem(XBINPUT.NumElems, Settings.MaxElNod)
    NumNodes_tot, XBELEM = Elem(XBINPUT, XBOPTS, XBELEM)
    
    XBINPUT.ForceStatic[2] = 12345.6
    PosIni, PhiNodes, ForceStatic, BoundConds = Node(XBINPUT, XBOPTS,\
                                                      NumNodes_tot, XBELEM)
    
    