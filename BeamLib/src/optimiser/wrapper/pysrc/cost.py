'''
# Salvatore Maraniello: 7 Aug 2014 --> 20 Aug 2014

This module contains methods to evaluate the cost functions related to the 
beam solver. Some of these have been wrapped from fortran code (f_total_mass &
f_node_disp), while the others are written in python and tested in this module.


Remark: The fortran routines computing the global cost and the constraints cost 
are not extracted. These can be rebuilt in python environment.

'''

import ctypes as ct
import numpy as np

import sys
sys.path.append('..')
import shared
import beamvar


# Load Dynamic Library
xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)

# dictionary of directions
DirDict={'x':0, 'y':1, 'z':2}

# ------------------------------------------------------------ Fortran routines
# The following routine are not extracted - not worth to built an interface
#fcost  =xb.__opt_cost_MOD_cost_global
#fconstr=xb.__opt_cost_MOD_cost_constraints

# Total structural mass
f_total_mass = xb.__opt_cost_MOD_cost_total_mass_wrap

# Nodal Displacement (absolute value)
f_node_disp  = xb.__opt_cost_MOD_cost_node_disp_wrap


#--------------------------------------------------------------- Python routines
def f_xyz_disp(PosIni,PosDef,**kwargs):
    ''' 
    Returns displacement along one direction with sign.
        Input:
        - PosIni,PosDef: initial and deformed nodal position
        - **kwargs:
            - direction: 'x','y','z'
            - NNode: node at which to evaluate the displacement
    '''
    
    dir   = kwargs.get('dir'      , 'z'                )
    NNode = kwargs.get('NNode'    , len(PosIni[:,0])-1 )
    
    try:
        jj=DirDict[dir]
    except KeyError:
        raise NameError('Specify a direction:x x, y, or z')

    # print 'extract: (%d,%d)' %(NNode,jj)
    return PosDef[NNode,jj] - PosIni[NNode,jj]



def f_xyz_vel(DynOut,RefVel,Time,NNode=-1,TTime=-1,dir='x',FoR='G'):
    '''
    Returns the 'dir' velocity (in the frame of reference FoR)  of the node 
    'NNode' at the time step TTime. The velocity is computed via a II order 
    polynomial interpolation
    
    Remark: tested in 
    /home/sm6110/git/SHARPy_studies/GEBM/20141024_qiqi_validation/postpr_hinge_testing.py
    
    '''
    import lib.postpr

    # computed Deformed Beam at each time-step in FoR
    THPosDef = lib.postpr.THPosDefGlobal(DynOut,RefVel,Time,set_origin=FoR)
    # and get velocities
    THVel=lib.postpr.compute_velocity(THPosDef,Time)    
    
    # get direction
    try:
        kk=DirDict[dir]
    except KeyError:
        raise NameError('Specify a direction:x x, y, or z')    
    
    vel = THVel[TTime,NNode,kk]

    return vel
    






def return_args(C,AttrList):
    '''Given a class C and a list of attributes to the class, the function
    creates a tuple with the class attributes '''
    
    T=[]
    for attr in AttrList:
        print attr
        T.append(getattr(C,attr))
    
    return tuple(T)
    


''' ----------------------------------------------------------- testing  --- '''
if __name__ == '__main__':
    
    
    
    #-------------------------------------------------------- f_xyz_disp testing
    PosIni=np.array([ [0.0, 10.0, 0.0 ],
                      [1.0, 11.0, 5.0 ] ])
    PosDef=np.array([ [1.0, 90.0, 4.0 ],
                      [5.0,  9.0, 2.0 ] ])
    
    print f_xyz_disp(PosIni,PosDef) # check defaults
    print f_xyz_disp(PosIni,PosDef,dir='x',NNode=0)
    print f_xyz_disp(PosIni,PosDef,dir='y',NNode=1)
    print f_xyz_disp(PosIni,PosDef,dir='z')
    
    
    #----------------------------------------------------------- DirDict testing
    print DirDict['x'],DirDict['y'],DirDict['z']
    #print DirDict['h']
    
    
    #------------------------------------------------- Test Velocity Calculation
    
    
    
    
    #----------------------------------------------------- Test return arguments
    class C():
        a=3
        b='lala'
        c=7
    print C.a, C.b, C.c
    
    t=return_args(C,['a','b','c'])
    print t
    

    
    
    
     











