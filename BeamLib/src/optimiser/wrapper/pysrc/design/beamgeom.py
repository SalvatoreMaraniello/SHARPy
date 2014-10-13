'''
Created on 13 Oct 2014

@author: sm6110

This module defines the beam geometry
'''

import numpy as np

def straight( PosIni, BeamLength=1.0, axis_vector=np.array([1,0,0]) ):
    ''' 
    Defines the geometry of a straight beam.
    
    By default the beam has unitary length and goes along the x axis and has always
    origin in (0,0,0)
    
    PosIni contains the initial nodes coordinates and can be initialised 
    using beamvar.fwd_static
    '''
    
    # check No. of nodes
    NumNodes = len(PosIni[:,0])
    
    # normalise unit vector
    axis_vector = BeamLength*axis_vector/np.linalg.norm(axis_vector)
    
    # allocate coordinates
    PosIni[:,0]=np.linspace(0,axis_vector[0],NumNodes)
    PosIni[:,1]=np.linspace(0,axis_vector[1],NumNodes)
    PosIni[:,2]=np.linspace(0,axis_vector[2],NumNodes)
    
    return PosIni


if __name__=='__main__':
    
    import lib.plot
    import matplotlib.pyplot as plt
    
    PosIni=np.zeros((21,3),dtype=float,order='F')
    
    ## check default values
    ##PosIni=straight(PosIni)
    
    # check general 
    PosIni=straight(PosIni,10.0,np.array([0,1,1]))
    
    lib.plot.sta_unif(PosIni, PosIni,hold=True)
    plt.xlabel('x')
    plt.ylabel('y')

    
    plt.show()
