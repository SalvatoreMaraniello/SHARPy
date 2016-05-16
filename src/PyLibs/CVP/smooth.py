'''
Created on 13 May 2016

@author: sm6110


Methods to build smooth polynomials (build_vec). To achieve a smooth reconstruction, the polynomial
is built starting from the value of its p-th derivative --- assumed to be constant over some user
defined intervals [ tc[i]; tc[i+1] ) --- and the values of its derivatives at t=0.

Routines summary:
- build_vec: builds the polynomial over a time domain tv
- heaviside_step: builds Heaviside step function. 
- singularity_fun: builds singularity function <t-a>^n
- estimate_parametrisation: given a function fv(t), the methods provides a rough estimation of the
   coefficients cfv which parametrise the curve fv(t).


'''

import numpy as np
import PyLibs.numerics.diff




def build_vec( tv, tc, cfv, ic=[0], p=1, EvalDeriv=True ):
    '''
    
    A polynomial of order p, f(t), is built starting from the values of its p-th derivative. This is 
    assumed to be a piecewise constant function with  values:
    
        d^p{f}/dt^p = cfv_i for t in [ tc_i, tc_{i+1} )
        
    where tc is an array of control points defining the intervals over which the p-th derivative is
    constant.
    
    To close the problem, the values of the 0, 1, ... {p-1}-th derivative at t=0 have to be 
    specified. These are collected into the list/array 'ic'.
    
    Input:
    - tv:  domain over which the polynomial f(t) is evaluated    
    - tc:  array of length NC containing the control points which specify the intervals over which
           the p-th derivative is constant. These do not have to be equally spaced but need to be 
           compatible with the time domain array, tv, i.e.:
               tc[0]=tv[0]    
               tc[-1]=tv[-1]    
    - cfv: values of the p-th derivative over the intervals defined by tc. Note that cfv has length
           NC-1 and that cfv[i] corresponds to the value in the interval  [ tc[i]; tc[i+1] )
    -ic:   list/array such that its n-th element is equal to the n-th derivative of the polynomial 
           at time t=0
               ic[n] =d^n{f}/d_t^n(t=0), where n = 0, 1, ... p-1
           note that len(ic)=p
    - Evalderiv: id true, all the derivatives of the polynomial up to the p-th order are returned
    
    
    Output:
    - Fv: if EvalDeriv is False, is an array containing the value of the polynomials over the time
          domain specified by tv.
          if EvalDeriv is true, is a p+1 rows matrix such that the n-the row contains the value of 
          the n-th derivative of the polynomial (and for n=0 the polynomial itself)
    
    '''
    
    NT=len(tv)
    NC=len(tc)
    
    assert NC-1==len(cfv) , 'len(cfv)==len(tc)-1 not verified!'
    assert len(ic)==p , 'len(ic)==p not verified!'
    
    Fv = np.zeros((p+1,NT))
     
    # derive difference between coefficients
    dv = np.zeros((NC-1,))
    dv[1:]=np.diff(cfv)

    # build ic vector
    F0=np.array( list(ic) + [cfv[0]] )


    for kk in range(p+1):
        # add initial condition contribution (ii=0, 1, 2 .... kk)
        for ii in range(kk+1):
            Fv[p-kk,:] += F0[p-ii]/np.math.factorial(kk-ii)*singularity_fun(tv, tc[0], kk-ii)  
        
        # evaluate p-k derivative for each singularity function
        for ii in range(NC-1):
            Fv[p-kk,:] += dv[ii]/np.math.factorial(kk)*singularity_fun(tv, tc[ii], kk)
    
    if EvalDeriv==False: Fv = Fv[0,:]        
      
    return Fv    
        



def heaviside_step(tv,a):
    '''
    Return Heaviside step function
        H(tv-a)
    '''
    
    hv = np.zeros((len(tv),))
    hv[tv>=a]=1.0
    
    return hv


def singularity_fun(tv,a,n):
    
    if n>= 0:
        sv = heaviside_step(tv,a) * (tv-a)**n    
    else:
        raise NameError('Method not implemented for n < 0')
    
    return sv



def estimate_parametrisation(tv,fv,tc,p):
    '''
    
    Given a signal to reconstruct, the function provides a quick estimation of the
    coefficients required to fit that signal with a smooth polynomial of order p
    defined over the control points tc. 
    
    '''

    icvec = np.empty((p,))

    df = fv
    for nn in range(p):
        # allocate IC
        icvec[nn]=df[nn]
        # differentiate
        df = PyLibs.numerics.diff.difffun(df,tv)
    
    
    NC=len(tc)
    cfv = np.empty((NC-1,))
    for ii in range(NC-1):
        ttvec = (tv>=tc[ii]) * (tv<tc[ii+1])
        cfv[ii] = np.average( df[ttvec] )  
    
    return cfv, icvec

        





if __name__=='__main__':
  
    import os
    import sys
    sys.path.append( os.environ["SHARPYDIR"]+'/src' )
    sys.path.append( os.environ["SHARPYDIR"]+'/src/Main' )
    import SharPySettings as Settings
    import PyLibs.numerics.diff  
    import matplotlib.pyplot as plt 
    
    T=3.0
    tv=np.linspace(0,T,1000)
    tc=np.linspace(0,T,7)
    
    cfv = np.array([10,-2,2,-2,1,-1])
    
    
    p=3
    Fv=build_vec(tv,tc,cfv,ic=[-10,1,2],p=p)
  

    
    # verify signal built
    FD = 0.0*Fv
    for ii in range(1,p+1):
        FD[ii,:] = PyLibs.numerics.diff.difffun(Fv[ii-1,:], tv)
        #plt.plot(tv, Fv[ii,:], 'r-o', markevery=10)
        #plt.plot(tv, FD[ii,:], 'b')
        #plt.legend(['analytical', 'Finite Differences'])
        #plt.show()
    
    # verify estimation  
    cfv_est, ic_est = estimate_parametrisation(tv,Fv[0,:],tc,p)
    
    #plt.plot(tv, Fv[ii,:], 'r-o', markevery=10)
    #plt.plot(tv, FD[ii,:], 'b')
    #plt.legend(['analytical', 'Finite Differences'])
    #plt.show()
    
    
    
    
















    