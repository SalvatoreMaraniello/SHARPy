'''
Created on 15 Oct 2014

@author: sm6110

Collects plot methods for time dependent simulations (structural dynamics and
rigid body dynamics)

'''

# packages and common variables imported in lib.plot.shared
from PyLibs.plot.shared import *

from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D


def traj(V):
    ''' 
    Given a time-history of vectors V[tt,:], where tt is the time-step, the
    function plots the trajectory described by the vector V[tt,:].
    
    '''
    
    if len(V[0,:]) is 2:
        
        plt.plot(V[0,0],V[0,1],'xr')
        plt.plot(V[:,0],V[:,1],color=grayshade)  # @UndefinedVariable

                 
    elif len(V[0,:]) is 3: 

        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d') # creates a Axes3D object
        ax.plot(V[:,0], V[:,1], V[:,2], zdir='z',color=grayshade)  # @UndefinedVariable
        
        # trick to get initial position
        V0 = np.ones((2,3))*V[0,:]
        ax.plot(V0[:,0], V0[:,1], V0[:,2], 'ro', zdir='z')  # @UndefinedVariabl                      ]
        #ax.plot(V[0,0], V[0,1], V[0,2], 'ro', zdir='z')  # @UndefinedVariable  
        #ax.plot(V[0,0], V[0,1], V[0,2], 'ro')  # @UndefinedVariable  
        #plt.plot(V[0,:],'ro') 
        
        # set axis limits
        xlim, ylim, zlim = set_3d_limits(V)
        
        #ax.set_aspect('equal')
        ax.set_xlim3d(-xlim, xlim)
        ax.set_ylim3d(-ylim, ylim)
        ax.set_zlim3d(-zlim, zlim)     

    
        
    return None
        

def axis_traj(V,Origin=None):
    ''' 
    Given a time-history of vectors V[tt,:], where tt is the time-step, the
    function plots all the vectors V[tt,:].
    
    The method is similar to traj, but:
        a. the whole position vector is plotted, from its origin to the tip
        b. if the origin of the vectors V[tt,:] changes, this can be accounted for
        in the plotting. In the traj method, this cannot be visualised
    '''
    
    def build_vec_for_plot(vec,orig):  
        Vplot = np.zeros((2,len(vec)))
        Vplot[ 1,:]=orig+vec
        Vplot[ 0,:]=orig
        return Vplot
        
    NT = len(V[:,0])
    
    if Origin is None:
        Origin=np.zeros(V.shape)
    
    if len(V[0,:]) is 2:
             
        for tt in range(NT):
            Vplot=build_vec_for_plot(V[tt,:],Origin[tt,:])
            plt.plot(Vplot[:,0],Vplot[:,1],color=grayshade)  # @UndefinedVariable
            
        Vplot=build_vec_for_plot(V[0,:],Origin[0,:])
        plt.plot(Vplot[:,0],Vplot[:,1],'r')    
                 
    elif len(V[0,:]) is 3: 

        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d') # creates a Axes3D object
        
        for tt in range(NT):
            Vplot=build_vec_for_plot(V[tt,:],Origin[tt,:])
            ax.plot(Vplot[:,0],Vplot[:,1],Vplot[:,2],zdir='z',color=grayshade)  # @UndefinedVariable    
        
        
        Vplot=build_vec_for_plot(V[0,:],Origin[0,:])
        ax.plot(Vplot[:,0],Vplot[:,1],Vplot[:,2],'r',zdir='z')        
        # set axis limits
        xlim, ylim, zlim = set_3d_limits(V+Origin)
        
        #ax.set_aspect('equal')
        ax.set_xlim3d(-xlim, xlim)
        ax.set_ylim3d(-ylim, ylim)
        ax.set_zlim3d(-zlim, zlim)     
       
    return None


def beam_snapshots(THPosDef,asp_ratio='equal'):
    '''
    Plots snapshots of the beam during a dynamic simulation. 
    
    THPosDef is output from one of these methods:
        lib.postpr.THPosDefGlobal
        lib.postpr.reshape_DynOut
     
    Output: fig, ax
    
    use:
        ax.set_xlabel('xlabel')
        ax.set_aspect(2)
        fig.savefig('equal.png')
        ax.set_aspect('auto')
        fig.savefig('auto.png')
        forceAspect(ax,aspect=1)
        fig.savefig('force.png')        
        
    Note: as opposed to axis_traj, there is no need to pass in input the origin
    of the FoR time history. If the beam is moving in space, the translation
    can be applied in lib.postpr.THPosDefGlobal
    '''
    
    Nsnap = len(THPosDef[:,0,0])
    
    fig = plt.figure()
    
    if len(THPosDef[0,0,:])==2:
        ax = fig.add_subplot(111)
        for tt in range(Nsnap):
            ax.plot(THPosDef[tt,:,0],THPosDef[tt,:,1],color=grayshade)
        ax.plot(THPosDef[0,:,0],THPosDef[0,:,1],color='r') 
              
    elif len(THPosDef[0,0,:])==3:
        ax = fig.add_subplot(111,projection='3d') 
        for tt in range(Nsnap):
            ax.plot(THPosDef[tt,:,0],THPosDef[tt,:,1],THPosDef[tt,:,2],color=grayshade)        
        ax.plot(THPosDef[0,:,0],THPosDef[0,:,1],THPosDef[0,:,2],color='r')  
    else:
        raise NameError('THPosDef is not 2D or 3D!!!')
    
    if asp_ratio is not None:
        ax.set_aspect(asp_ratio)
    
    return fig, ax
    


def polytrend(xv,tv,Fmat,fig,**kwargs):
    '''
    The function visualises the quantity Fmat, continuously defined over the 
    domain xv, at the discrete times tv.
    Fmat has size (len(xv),len(tv))
    '''

    #------------------------------------------------ Extract optional arguments
    #if 'figure' is kwargs:
    #    fig = kwargs['figure']
    #else: 
    #    print('creating a brand new figure')
    #    fig = plt.figure('Poly trend',(10,8))
        
    if 'cmap' in kwargs:
        cmap=kwargs['cmap']
    else:
        cmap=mpl.cm.winter     
     
    
    ax = fig.gca(projection='3d')
    
    verts = []
    for ii in range(len(tv)):
        
        fvec = np.abs(Fmat[:,ii])
    
        #####create vertices for polygon
        # to the set of points given in input, additional edges need to be 
        # included to 'close' the polygons
        XS=np.concatenate([[xv[0]],xv,[xv[-1]]])
        YS=np.concatenate([[0],fvec,[0]])
        verts.append(list(zip(XS, YS)))    
    
    # create a list of different colours from a colour map for each facet
    ticknumber=np.linspace(0,1,len(tv))
    clist=[]
    for nn in ticknumber:
        clist.append(cmap(nn))
    
    poly = PolyCollection(verts, facecolors = clist, edgecolors='0.6' )#, edgecolors=clist)
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=tv, zdir='y')
    
    ax.grid(False)
    
    return fig, ax   
    
    
    
    


if __name__=='__main__':
    
    NT = 30
    V = np.zeros((NT,3))
    time = np.linspace(0,2*np.pi,NT)
    V[:,0]=np.cos(time)
    V[:,2]=np.sin(time)
    V[:,1]=np.ones(NT)
    
    axis_traj(V[:-5,:])
    plt.show()
    
    O = np.zeros((NT,3))
    O[:,1]=np.linspace(0,2,NT)
    
    axis_traj(V,O)
    plt.show()
    
    
    
