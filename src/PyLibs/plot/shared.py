'''
Created on 15 Oct 2014

@author: sm6110

Shared settings for plotting routines.

Parameters change from a matplotlib version to another. Thus a record of 
different dictionaries is kept. The default is set to params.

From matplotlib  > v.1.3 
matplotlib.RcParams.find_all()
e.g.
>>> print(matplotlib.rcParams.find_all('\.size'))
RcParams({'font.size': 12,
          'xtick.major.size': 4,
          'xtick.minor.size': 2,
          'ytick.major.size': 4,
          'ytick.minor.size': 2})

'''


import numpy as np
import scipy.interpolate

import matplotlib as mpl
import matplotlib.pyplot as plt

from warnings import warn






#---------------------------------------------------------------- Common setting
### incompatible with matplotlib versions higher then 1.1.1
#color_cycle_list=['k', 'b', 'r', 'g', 'y', 'c']
#mpl.axes.set_default_color_cycle(color_cycle_list)
grayshade='0.6'
grayshade_minor='0.4'    
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}



# -----------------------------------------------------------------------------
fontlabel = 20 # if you want it different
fontlabel_pres = 24



# ----------------------------------------------------- set  parameters articles
# params adjusted to be plot on a scale 0.4 



### v1.1.1rc
params111rc = {'legend.fontsize': 18,
               'legend.linewidth': 2,
               'font.size': 16,       # for all, e.g. the ticks
               'legend.numpoints': 1} # for plotting the marker only once

### v 1.4.2
params142   = {'legend.fontsize': fontlabel,
                'font.size': fontlabel,       # for all, e.g. the ticks
                'xtick.labelsize': fontlabel-2,
                'ytick.labelsize': fontlabel-2, 
               'figure.autolayout': True,
               'legend.numpoints': 1} # for plotting the marker only once


### v 1.4.2
params142_pres = {'legend.fontsize': fontlabel_pres,
                  'font.size': fontlabel_pres,       # for all, e.g. the ticks
                  'xtick.labelsize': fontlabel_pres-2,
                  'ytick.labelsize': fontlabel_pres-2, 
                  'figure.autolayout': True,
                  'legend.numpoints': 1} # for plotting the marker only once

params = params142
params_pres = params142_pres


#-------------------------------------------------------------------------------------------- Colors

cdict = {
'dark-blue'    : '#003366',
'royal-blue'   : '#4169E1', 
'blue-violet'  : '#6666FF', 
'clear-red'    : '#CC3333',  
'dark-green'   : '#336633',
'orange'       : '#FF6600'
         }






def update_by_font(params,fontsize):
    ''' Scale all parameters based on fontsize
        pltext is an instance for matplotlib.pyplot '''
    params['legend.fontsize'] = fontsize
    params['font.size'] = fontsize       # for all, e.g. the ticks
    params['xtick.labelsize'] = fontsize-2
    params['ytick.labelsize'] = fontsize-2 
    params['legend.numpoints'] = 1 # for plotting the marker only once
    return params


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
 
 
def make_grid_from_array(x,y,z):

    # checks
    N = len(x)
    if len(y)!=N or len(z)!=N:
        warn('The inputs need to be 1D arrays of the same length')
    
    # set up a dictionary:
    raw_data = {}
    Xset,Yset=set(),set()

    # Prepare Interpolation
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')


    for cc in range(N):
        raw_data[( x[cc] , y[cc] )] = z[cc] 
        Yset.add(y[cc])
        Xset.add(x[cc])
    
    Xset =sorted(list(Xset))
    Yset=sorted(list(Yset))
    
    Xmat,Ymat = np.meshgrid(Xset,Yset)    
    Zmat = []
    
    for xx in Xset:
        row = []
        for yy in Yset:
            try:
                row.append(raw_data[(xx,yy)]) # produces transposed plots
            except:
                warn('Data for (%f,%f) required interpolation!' %(xx,yy) )
                zz = rbf(xx, yy)
                row.append(zz)                
        Zmat.append(row)
    Zmat = np.array(Zmat)
    Zmat = Zmat.transpose() # fix added cause all the plots were the transpose of what expected
    
    return Xmat,Ymat,Zmat
 
 


  
