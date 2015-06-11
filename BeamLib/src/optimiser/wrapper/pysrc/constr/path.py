'''
Created on 8 Jun 2015

@author: sm6110


General methods to handle path constraints


'''

import numpy as np
import constr.aggregate


def sampling1D(xv,fv,fmax,NS):
    '''
    Given the function fv defined over xv, and a number of sample points NS, 
    the function will return the constraint vector:
        g = fv - fmax
    such that g_i < 0 is the constraint f_i < fmax is not violated at the sample
    points.
    
    Note that:
    - NS includes the extrema xv[0] and xv[-1]
    - if the sample points do not coincide with any value in xv, a piecewise
    linear interpolation will be used
    '''
    
    xs = np.linspace(xv[0],xv[-1],NS)
    
    fs = np.interp(xs,xv,fv)
    
    return fs - fmax



def KSfun(xv,fv,fmax,NI,k=3.0):
    '''
    Given the function fv defined over xv, and a number of intervals, the 
    function will apply on each interval the KS functional to return NI constraint
    of the kind:
        g = KS(fv) - fmax
    such that g_i < 0 is the constraint is not violated in the interval NI.
    
    The KS functional integration is based on constr.aggregate.KSsum, with alpha
    parameter chosen equal to the dt of integration, and m automatically computed.
    k is the only reasonable parameter to pass.
    
    If NI is such that there are no points at the extrema of the intervals, no
    interpolation is performed. The intervals are simply restricted to the nearest
    available point. The assumption is acceptable as all points will be included
    and, even if the maximum happens at an extrema (thus in a point with low
    integration weight associated), changes in the fv are smooth so the error 
    will be acceptable.
    
    
    '''
    # divide into equal intervals:
    xext = np.linspace(xv[0],xv[-1],NI+1)
    alpha=1.0#xv[1]-xv[0]
    
    # compute KSv
    gv = np.zeros((NI,))    
    for ii in range(NI):
        iivec = (xv>xext[ii]) * (xv<xext[ii+1])
        fint = fv[iivec]
        gv[ii] = constr.aggregate.KSsum(fint, k, alpha) - fmax
        
    return gv    
        
        
    
def max(xv,fv,fmax,NI):
    '''
    Given the function fv defined over xv, and a number of intervals, the 
    function will compute on each interval the max(fv) to return NI constraint
    of the kind:
        g = max(fv) - fmax
    such that g_i < 0 is the constraint is not violated in the interval NI.
    
    '''
    # divide into equal intervals:
    xext = np.linspace(xv[0],xv[-1],NI+1)
    alpha=1.0#xv[1]-xv[0]
    
    # compute KSv
    gv = np.zeros((NI,))    
    for ii in range(NI):
        iivec = (xv>xext[ii]) * (xv<xext[ii+1])
        fint = fv[iivec]
        gv[ii] = np.max(fint) - fmax
        
    return gv     
 
 
    


if __name__=='__main__':

    import matplotlib.pyplot as plt
        
    xv=np.linspace(0,2*np.pi,2.0/0.0002)
    fv = 2.0*np.sin(3.0*xv) -1.0*np.cos(xv)
    plt.plot(xv,fv)

    """
    # check sampling1D
    NS=20
    xs=np.linspace(xv[0],xv[-1],NS)
    gs = sampling1D(xv,fv,3.0,NS)
    plt.plot(xs,gs+3.0,'r')
    """
    
    # check KSfun
    NI=4
    xs = np.linspace(0,2*np.pi,NI+1)
    for xx in xs:
        plt.plot(xx*np.ones((2,)),np.array([-6,6]),'r')
 
    fmax=2.0
    gv = KSfun(xv,fv,fmax,NI,k=20.0)
    #gv = max(xv,fv,fmax,NI)
    
    print 'KS functional with fmax = %f' %(fmax)
    print gv
  
    xext=np.linspace(0,2.0*np.pi,NI+1)
    for ii in range(NI):
        plt.plot(xext[ii:ii+2],gv[ii]*np.ones(2)+fmax,'r')
  
    plt.show()  
    
     
    
    
    
    
