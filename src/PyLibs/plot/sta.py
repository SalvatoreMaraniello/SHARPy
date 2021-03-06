'''
Created on 19 Aug 2014

@author: sm6110

plotting tools

ref. http://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
     http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
     http://matplotlib.org/api/pyplot_api.html

Remarks: Following functions merged from 'monoopt.myplots.threedim' module:
    - make_grid_from_array
    - surface
    - contour
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import shared




#---------------------------------------------------------------- Common setting

grayshade='0.6'
grayshade_minor='0.4'    
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}


def sta_unif(PosIni,PosDef,equal=False,hold=False):
    ''' Tool to visualise initial and deformed beam configuration '''
    
    fig = plt.figure()
    #ax = fig.gca(projection='3d',aspect='equal')
    #ax.plot(PosIni[:,0], PosIni[:,1], PosIni[:,2], zdir='z', label='Ini.',color='k',aspect='equal')
    #ax.plot(PosDef[:,0], PosDef[:,1], PosDef[:,2], zdir='z', label='Def.',color='b',aspect='equal')    

    ax = fig.add_subplot(111,projection='3d') # creates a Axes3D object
    ax.plot(PosIni[:,0], PosIni[:,1], PosIni[:,2], zdir='z', label='Ini.',color='k')
    ax.plot(PosDef[:,0], PosDef[:,1], PosDef[:,2], zdir='z', label='Def.',color='b')
    
    if equal==True:
        ax.set_aspect('equal')
    
    ax.legend()
    #ax.set_xlim3d(0, 5)
    #ax.set_ylim3d(0, 5)
    #ax.set_zlim3d(0, 5)
    
    if hold is False:
        plt.show()
    
      
def linlin2Dplot(x,y,xlab,ylab,labels,savename,**kwargs):
    ''' plot with linear axis and standard format 
        x,y: array, these can be a list of arrays or single arrays 
        savename: full path to destination file without format
        fig_format = a list with any supported format:
        http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.savefig 
        hold: if True, the figures are not saved but returned to prompt for
        further manipulations    
    '''

    # Optional Input
    if not(kwargs.has_key('yrotation')):
        kwargs['yrotation']=0.0
    if not(kwargs.has_key('format_list')):
        kwargs['format_list']=['pdf','eps','ps']    
    if not(kwargs.has_key('hold')):
        kwargs['hold']=False
    if not(kwargs.has_key('style_dict')):
        kwargs['style_dict']={}
    if not(kwargs.has_key('legend_location')):
        kwargs['legend_location']=1
    
    if (type(x) is list) and (type(y) is list):
        N=len(x)
        if len(y) != N:
            raise NameError('x and y list must have same length')
    else:
        N=0
        x,y=[x],[y]

   
    # plot
    for ii in range(N):
        plt.plot(x[ii],y[ii],**_produce_dictionary(kwargs['style_dict'],ii))
        plt.hold('True')

    lim=[np.min(x),np.max(x),np.min(y),np.max(y)]
    xrange, yrange = lim[1]-lim[0], lim[3]-lim[2]
    plt.xlabel(xlab,fontdict=font)
    plt.ylabel(ylab,rotation=kwargs['yrotation'],fontdict=font)
    plt.axis([ lim[0]-0.1*xrange, lim[1]+0.1*xrange, lim[2]-0.1*yrange, lim[3]+0.1*yrange ]) 
    # set grid
    plt.grid(b=True, which='major', color=grayshade      , linestyle='-')
    plt.grid(b=True, which='minor', color=grayshade_minor, linestyle='--')    
    
    plt.legend(labels,loc=kwargs['legend_location'])
    
    plt.show()
    if kwargs['hold'] == False:
        
        for ff in kwargs['format_list']:
            print ff
            plt.savefig(savename+'.'+ff)
        plt.close()

    ### Useful Options
    #plt.legend(['leg a','leg b'],loc='upper right')
    #plt.xticks([i*0.1 for i in range(0,11)])
    #plt.yticks([i*0.2 for i in range(-7,3)])
    #plt.axis([xmin,xmax,ymin,ymax])    
    #ax = plt.gca()
    #ax.invert_yaxis()

    return 


 
def _produce_dictionary(kwargs_list,ii):
    ''' kwargs_list is either an empty dictionary or a dictionary whose items
    values are list. 
    From kwargs_list, a new dictionary, with same keys but ii-th values ''' 
    
    kwargs = {}
    if len(kwargs_list)>0:
        for kk in kwargs_list.iteritems():
            kwargs[kk[0]] = kk[1][ii]  
            
    print 'dictionary generated: ', kwargs
    return kwargs

      
def surface(x,y,z):

    Xmat,Ymat,Zmat = shared.make_grid_from_array(x,y,z)
    
    fig=plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(Xmat,Ymat,Zmat,rstride=1,cstride=1,cmap=mpl.cm.jet,linewidth=0)

    return fig,ax


def contour(x,y,z,**kwargs):
    ''' 
    Simplifies the use of matplotlib.pyplot.contour/contourf function
    
    remarks:
    - For colormaps look at matplotlib.cm
    - Contour usage tips:
        - CS.cl_cvalues: returns an array with all the level lines values   
        - cset1 = contourf(X, Y, Z, levels,cmap=mpl.cm.get_cmap(cmap, len(levels)-1),
                        norm=norm)  
    
    '''      
 
    # kwargs 
    if not(kwargs.has_key('colorbar_flag')):
        kwargs['colorbar_flag']=False 
    if not(kwargs.has_key('contourf')):
        kwargs['contourf']=False

    # build grids
    Xmat,Ymat,Zmat = shared.make_grid_from_array(x,y,z)
 
    if not(kwargs.has_key('Nlevels')):    
        args_in = (Xmat,Ymat,Zmat)
    else:
        args_in = (Xmat,Ymat,Zmat,kwargs['Nlevels'])
    
    if kwargs['contourf']==True:
        CS=plt.contourf(*args_in,**kwargs)
    else:
        CS = plt.contour(*args_in,**kwargs)
        plt.clabel(CS, fontsize=10, **kwargs)
    if kwargs['colorbar_flag'] is True:
        CB = plt.colorbar(CS,extend='both')
    
    return
  
 