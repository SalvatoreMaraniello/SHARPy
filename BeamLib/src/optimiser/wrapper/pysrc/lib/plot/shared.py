'''
Created on 15 Oct 2014

@author: sm6110

Shared settings for plotting routines

'''


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from warnings import warn

import scipy.interpolate


#---------------------------------------------------------------- Common setting

color_cycle_list=['k', 'b', 'r', 'g', 'y', 'c']
mpl.axes.set_default_color_cycle(color_cycle_list)

grayshade='0.6'
grayshade_minor='0.4'    
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}





def set_3d_limits(XYZ,method='sym'):
    
    ''' 
    Given a matrix XYZ, whose columns are the x, y, z values of a trajectory or
    set of position vectors, the methods sets the axis for the 3d plot.
    
    method:
        - 'sym': find the max absolute coordinate in the x, y, and z direction
    and returns xlim, ylim, zlim as positive scalars
    '''
    
    if method=='sym':
        # set axis limits
        xlim = 1.2*np.max(np.abs(XYZ[:,0]))
        ylim = 1.2*np.max(np.abs(XYZ[:,1]))
        zlim = 1.2*np.max(np.abs(XYZ[:,2]))
            
        xyzmax = 1.2*np.max(np.max(np.abs(XYZ)))
        tol=1e-2*xyzmax
        if xlim<tol:
            xlim=xyzmax
        if ylim<tol:
            ylim=xyzmax
        if zlim<tol:
            zlim=xyzmax
        
    
    return xlim, ylim, zlim
 
 
 


  
