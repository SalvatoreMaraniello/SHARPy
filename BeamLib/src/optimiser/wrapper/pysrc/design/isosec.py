'''

Salvatore Maraniello. 12/Aug/2014

Module to compute stiffness and mass matrices of isotropic beam cross sections.

The module wraps the functions contained into lib_isosec.f90 to allow design
definition outside the Fortran environment.

Input and Output of these routines are Python variables only (ctype interface
contained in the methods below).

The list of subroutine implemented is:

 - getprop: returns properties of material
      - alluminium
      - default  (E  = 1.0d7, G  = 1.0d5, rho= 1.0_8)
 - rect: returns mass and stiffness matrices of a beam element with full
   rectangular cross-section.
 - hollowrect:returns mass and stiffness matrices of a beam element with
   hollow rectangular cross-section.
 - ellip: full isotropic elliptical cross-section
 - hollowellip: isotropic hollow elliptical cross-section
 - circ: full isotropic circular cross-section
 - hollowcirc: isotropic hollow circular cross-section

# where the wrapped functions are:
getprop     = xb.__lib_isosec_MOD_getprop
circ        = xb.__lib_isosec_MOD_isocirc
ellip       = xb.__lib_isosec_MOD_isoellip
hollowcirc  = xb.__lib_isosec_MOD_isohollowcirc
hollowellip = xb.__lib_isosec_MOD_isohollowellip
hollowrect  = xb.__lib_isosec_MOD_isohollowrect
rect        = xb.__lib_isosec_MOD_isorect


->Reference:
 - http://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
 
             !!! Important: For future development   !!! 
             !!! ensure the numerical input that     !!!
             !!! can serve for design are always the !!!
             !!! first ones to be passed - see also  !!!
             !!! defbeam.py                          !!!

'''

import sys
sys.path.append('../..')

import ctypes as ct
import numpy as np

import shared
import beamvar

# Load Dynamic Library
xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)


# ------------------------------------------------------------------------------
def rect(l2, l3, material):

    # Prepare Interface
    ct_l2 = ct.c_double( l2 )
    ct_l3 = ct.c_double( l3 )
    
    ct_material = ct.c_char_p( material.ljust(5)  )
    
    M =  np.zeros((6,6),dtype=float,order='F')
    K =  np.zeros((6,6),dtype=float,order='F')
 
    # call Fortran function
    Ffun = xb.__lib_isosec_MOD_isorect
    Ffun( ct.byref(ct_l2), ct.byref( ct_l3 ), ct_material,               # input
          M.ctypes.data_as(ct.c_void_p), K.ctypes.data_as(ct.c_void_p)) # output
   
    return M, K


#------------------------------------------------------------------------------- 
def rect_fact_torsion(l2,l3,material,factor=1e3):
    ''' 
    Rectangular cross-section with factor applied to torsional stiffness.
    The default factor value is 1e3 to make to avoid large torsional rotations
    '''
    
    M, K = rect(l2,l3,material)
    K[3,3]=factor*K[3,3]
    
    return M, K


# ------------------------------------------------------------------------------
def hollowrect(l2, l3, t2, t3, material):

    # Prepare Interface
    ct_l2 = ct.c_double( l2 )
    ct_l3 = ct.c_double( l3 )
    ct_t2 = ct.c_double( t2 )
    ct_t3 = ct.c_double( t3 ) 
    
    ct_material = ct.c_char_p( material.ljust(5)  )
    
    M =  np.zeros((6,6),dtype=float,order='F')
    K =  np.zeros((6,6),dtype=float,order='F')
    
    # call Fortran function
    Ffun = xb.__lib_isosec_MOD_isohollowrect
    Ffun( ct.byref(ct_l2), ct.byref( ct_l3 ),
          ct.byref(ct_t2), ct.byref( ct_t3 ), 
          ct_material,                                                 
          M.ctypes.data_as(ct.c_void_p), K.ctypes.data_as(ct.c_void_p)) # output
   
    return M, K


# ------------------------------------------------------------------------------
def ellip(l2, l3, material):

    # Prepare Interface
    ct_l2 = ct.c_double( l2 )
    ct_l3 = ct.c_double( l3 )
    
    ct_material = ct.c_char_p( material.ljust(5)  )
    
    M =  np.zeros((6,6),dtype=float,order='F')
    K =  np.zeros((6,6),dtype=float,order='F')
 
    # call Fortran function
    Ffun = xb.__lib_isosec_MOD_isoellip
    Ffun( ct.byref(ct_l2), ct.byref( ct_l3 ), ct_material,               # input
          M.ctypes.data_as(ct.c_void_p), K.ctypes.data_as(ct.c_void_p)) # output
   
    return M, K


# ------------------------------------------------------------------------------
def hollowellip(l2, l3, t2, t3, material):

    # Prepare Interface
    ct_l2 = ct.c_double( l2 )
    ct_l3 = ct.c_double( l3 )
    ct_t2 = ct.c_double( t2 )
    ct_t3 = ct.c_double( t3 ) 
    
    ct_material = ct.c_char_p( material.ljust(5)  )
    
    M =  np.zeros((6,6),dtype=float,order='F')
    K =  np.zeros((6,6),dtype=float,order='F')
    
    # call Fortran function
    Ffun = xb.__lib_isosec_MOD_isohollowellip
    Ffun( ct.byref(ct_l2), ct.byref( ct_l3 ),
          ct.byref(ct_t2), ct.byref( ct_t3 ), 
          ct_material,                                                 
          M.ctypes.data_as(ct.c_void_p), K.ctypes.data_as(ct.c_void_p)) # output
   
    return M, K


# ------------------------------------------------------------------------------
def circ(r, material):

    # Prepare Interface
    ct_r = ct.c_double( r )
    
    ct_material = ct.c_char_p( material.ljust(5)  )
    
    M =  np.zeros((6,6),dtype=float,order='F')
    K =  np.zeros((6,6),dtype=float,order='F')
 
    # call Fortran function
    Ffun = xb.__lib_isosec_MOD_isocirc
    Ffun( ct.byref( ct_r ), ct_material,               # input
          M.ctypes.data_as(ct.c_void_p), K.ctypes.data_as(ct.c_void_p)) # output
   
    return M, K


# ------------------------------------------------------------------------------
def hollowcirc(r, t, material):

    # Prepare Interface
    ct_r = ct.c_double( r )
    ct_t = ct.c_double( t )
 
    ct_material = ct.c_char_p( material.ljust(5)  )
    
    M =  np.zeros((6,6),dtype=float,order='F')
    K =  np.zeros((6,6),dtype=float,order='F')
    
    # call Fortran function
    Ffun = xb.__lib_isosec_MOD_isohollowcirc
    Ffun( ct.byref(ct_r), ct.byref( ct_t ),
          ct_material,                                                 
          M.ctypes.data_as(ct.c_void_p), K.ctypes.data_as(ct.c_void_p)) # output
   
    return M, K

''' -------------------------------------------------------------  Testing '''

if __name__ == '__main__':
    
    material='allum'
 
  
    print 'Test rect...'
    l2=0.2; l3=0.5;
    M,K=rect(l2, l3, material)
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)
    
    print 'Test rect_fact_torsion...'
    M,K=rect_fact_torsion(l2, l3, material)
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)
    
    print 'Test hollowrect...'
    l2=0.2; l3=0.5;
    t2=0.1; t3=0.25;
    M,K=hollowrect(l2, l3, t2, t3, material)
    print 'expected full...'
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)
    t2=0.0; t3=0.0;
    M,K=hollowrect(l2, l3, t2, t3, material)
    print 'expected zero...'
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)


    print 'Test ellip...'
    l2=0.1; l3=0.25;
    M,K=ellip(l2, l3, material)
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)


    print 'Test hollowellip...'
    l2=0.1; l3=0.25;
    t2=0.1; t3=0.25;
    M,K=hollowellip(l2, l3, t2, t3, material)
    print 'expected full...'
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)
    t2=0.0; t3=0.0;
    M,K=hollowellip(l2, l3, t2, t3, material)
    print 'expected zero...'
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)


    print 'Test circ...'
    r=0.2
    M,K=circ(r, material)
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)


    print 'Test hollowcirc...'
    r=0.2; t=r;
    M,K=hollowcirc(r, t, material)
    print 'expected full...'
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)
    t=0.0;
    M,K=hollowcirc(r, t, material)
    print 'expected zero...'
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)


    print 'Matching NCB1 test case:'
    M,K = rect(l2,l3,'lala')
    print 'Mass Diagonal', np.diag(M)
    print 'Stiffness Diagonal', np.diag(K)    



