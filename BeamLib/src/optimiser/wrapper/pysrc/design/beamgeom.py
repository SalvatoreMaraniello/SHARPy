'''
Created on 13 Oct 2014

@author: sm6110

This module defines the beam geometry
'''

import numpy as np





def straight( PosIni, BeamLength=1.0, axis_vector=np.array([1,0,0]), shift=0.0 ):
    ''' 
    Defines the geometry of a straight beam.
    
    By default the beam has unitary length and goes along the x axis and has always
    origin in (0,0,0)
    
    PosIni contains the initial nodes coordinates and can be initialised 
    using beamvar.fwd_static
    
    WARNING: SHIFT NOT WORKING! The shift was introduced when attempting to
    move the spherical joint position by simply defining a beam that started
    from a negative position. [9 sep 2015] 
    Issues arise in the fortran code, as matrices allocated as sparse do not
    have enough memory to store negative numbers (arising from the negative)
    positions
    '''
    
    # check No. of nodes
    NumNodes = len(PosIni[:,0])
    
    # normalise unit vector
    axis_vector = BeamLength*axis_vector/np.linalg.norm(axis_vector)
    
    ### sm 9 sep 2015: shift introduced to allow for hinge to be placed at
    ### half beam without changing the solver.
    ### Fortran issues as the sparse matrices do not have enough memory when
    ### numers become negative!!!

    #shift=-0.5*BeamLength
    shiftperc=shift/BeamLength
    #PosIni[:,0]=np.linspace(0,axis_vector[0],NumNodes)
    #PosIni[:,1]=np.linspace(0,axis_vector[1],NumNodes)
    #PosIni[:,2]=np.linspace(0,axis_vector[2],NumNodes)
    low=shiftperc*axis_vector
    up=(1.0+shiftperc)*axis_vector
    PosIni[:,0]=np.linspace(low[0],up[0],NumNodes)
    PosIni[:,1]=np.linspace(low[1],up[1],NumNodes)
    PosIni[:,2]=np.linspace(low[2],up[2],NumNodes)
    
    return PosIni


if __name__=='__main__':
    
    #import lib.plot
    import matplotlib.pyplot as plt
    
    PosIni=np.zeros((21,3),dtype=float,order='F')
    
    ## check default values
    ##PosIni=straight(PosIni)
    
    # check general 
    PosIni=straight(PosIni,10.0,np.array([1,0,0]))
    
    #lib.plot.sta_unif(PosIni, PosIni,hold=True)

    plt.plot(PosIni[:,0],PosIni[:,1])
    plt.xlabel('x')
    plt.ylabel('y')    
    plt.show()
