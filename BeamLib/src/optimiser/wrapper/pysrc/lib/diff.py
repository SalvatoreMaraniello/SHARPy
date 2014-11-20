'''
Created on 29 Oct 2014

@author: sm6110

Low level methods for differentiation

'''

import numpy as np



def dwgt(xv,xp,k=0):
    ''' 
    Given a stencil of coordinates xv, the function returns a row vector wv
    such that, given a vector valued function fv=f(xv), the k-th derivative of
    f at the point xp is:
    
    d^k f = np.dot(wv,fv) 
    
    If k is not passed in input, the function will perform by default a 
    polynomial interpolation.
    
    Remark: the weight wv do not depend on the function value
    
    '''
   
    # initialise and check
    Nx=len(xv) 
    if k>Nx:
        raise NameError('Derivation order not compatible with stencil dimension!')
    
    # change FoR
    xv = xv-xp
    
    # get base polynomials coefficients
    I=np.eye(Nx)
    C=np.zeros((Nx,Nx))
    for ii in range(Nx):
        C[ii,:]=np.polyfit(xv,I[ii,:],Nx-1)
        
    # get weights
    wv = C[:,Nx-k-1]*np.math.factorial(k)
         
    return wv



def difffun(f,x):
    '''
    Compute II order accuracy first derivative of a function f over the equally
    spaced vector x.
    
    The code is specifically set for II order accuracy 1st derivatives.
    To improve it:
        a. the boundary treatment should be modified.
        b. sarr variable is now specific for II order accuracy
    '''
    
    # get weights for II order polynomial derivative
    # remark: if a general order of derivation is required, special attention
    # should be put in evaluating the derivatives at the boundary 
    kder=2
    KK=kder+1
    w0   = dwgt( x[:KK], x[ 0], k=1)
    wend = dwgt(x[-KK:], x[-1], k=1)
    wmid = dwgt( x[:KK], x[ 1], k=1)
    
    # 
    NumSteps = len(f)-1
    df = np.zeros(NumSteps+1)
    
    # timestep 1 and final:
    df[ 0] = np.dot(  w0,f[ :KK])
    df[-1] = np.dot(wend,f[-KK:])
            
    # timestep 2... NumSteps-1
    # sarr specific for II order derivative
    sarr = np.array([-1,0,1]) 
    for tt in range(1,NumSteps,1):
        ttind = tt+sarr
        df[tt] = np.dot(wmid,f[ttind])
     
    return df




if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    
    
    
    
    #------------------------------------------------------------- check difffun
    x=np.linspace(0,2*np.pi,100)
    f=np.sin(x)
    
    df = difffun(f,x)
    plt.plot(x,f,'r')
    plt.plot(x,df,'b')
    plt.show()
    
    
    
    
    
    #---------------------------------------------------------------- check dwgt

    # analytical
    x = np.array([-1,0,1,2])
    f = 2.0*x**2 - x + 4.0
    df = 4.0*x - 1.0
    ddf = 4.0 + 0.0*x
    
    #numerical
    N=20
    xv=np.linspace(-2,3,N)
    fv=np.zeros(N)
    dfv=np.zeros(N)
    ddfv=np.zeros(N)
    
    for ii in range(N):
        fv[ii]=np.dot(dwgt(x,xv[ii],0),f)
        dfv[ii]=np.dot(dwgt(x,xv[ii],1),f)        
        ddfv[ii]=np.dot(dwgt(x,xv[ii],2),f)  
    
    plt.figure('interpolation')
    plt.plot(x,f,'ro')
    plt.plot(xv,fv,'k')
    plt.show()
    
    plt.figure('I der')
    plt.plot(x,df,'ro')
    plt.plot(xv,dfv,'k')
    plt.show()  
      
    plt.figure('II der')
    plt.plot(x,ddf,'ro')
    plt.plot(xv,ddfv,'k')
    plt.show()    
    
    
    
    
     
     
     
    
        
    
    
    
    
    