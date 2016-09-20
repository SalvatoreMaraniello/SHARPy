'''

@author: Salvatore Maraniello
@summary: contains a collection of constraint methods for optimal
        control problems using a Fourier series based control vector 
        parametrisation
        
@warning: untested!

'''

import numpy as np
import PyLibs.CVP.fourier


def bound_constraint_DSS(acfv,fcfv,Ymax,Nsamp=4,scaler=1.0):
    ''' 
    Returns an array with all the constraint necessary to impose a boundary 
    to a control discretised using a DS series. The constraints are all of the 
    form:
    
    scaler**2 * ( Ymax**2 - q(acfv) ) > 0
    
    where q is a quadratic form returning a positive value. 
    The method is not valid for a generic sum of sine terms, as it assumes that
    all the frequencies are of the form:
        fn = n fo, fo = 1/2T
    
    Accuracy: depends on the number of sampling points Nsamp per wave. 

    '''

    Ncf = len(acfv)
    
    # find base and max frequencies
    f0, fmax = np.min(fcfv), np.max(fcfv)    

    #-------------------------------------------------------------------- checks
    # a. if f0=0, the next smallest frequency is found
    # b. checks if fmax is multiple of f0. Note this is necessary, but not 
    #    sufficient, to state fcfv builds a DS series

    tol=1e-6
    
    if np.abs(f0-0.0)<tol:
        iimin=np.argmin(fcfv)
        # build array without minimum
        fcfvnnz=np.concatenate((fcfv[:iimin],fcfv[iimin+1:]))
        f0 = np.min(fcfvnnz)
        if np.abs(f0-0.0)<tol:
            print ('fcfv=', fcfv)
            raise NameError("f0=%f found twice in the frequency array. Check your input!" %(f0))
        
    # check fmax/fo is an integer
    Nmax=np.round(fmax/f0)
    if np.abs(fmax/f0 - Nmax)>tol:
        raise NameError("The max. frequency fmax=%f is not a multiple of the base frequency f0=%f... you sure this is Discrete Sine series???" %(f0,fmax))
        
    # create sampling domain over non dimensional time tv:
    # tv = tdim / (2T)
    # only half of the time domain needs to be sampled
    tv = np.linspace(0,0.5,0.5*Nmax*Nsamp+1)
    
    # evaluate the series at these points:
    yv =  PyLibs.CVP.fourier.sin_series_vec(tv, acfv, fcfv/f0)
    
    ##### debug
    ##import matplotlib.pyplot as plt
    ##plt.plot(tv,yv)
    ##plt.show()
    
    # return constraint value
    gv = scaler**2 * ( Ymax**2 - yv**2 )
    
    return gv
    


if __name__=='__main__':
    

    import matplotlib.pyplot as plt
    
    T=2.0
    NT=501 
    Ncf=16
    acfv=1.5*np.random.rand(Ncf)-0.5
    
    
    #f0=1.0/(2.*T)
    #fcfv=np.array([ (nn+1)*f0 for nn in range(Ncf) ])
    
    tv = np.linspace(0,T,NT)
    (yv,fcfv)=PyLibs.CVP.fourier.DSSvec(tv,acfv,EvalDeriv=False)
    
    #yv = PyLibs.CVP.fourier.sin_series_vec(tv, acfv, fcfv)
    
    # find limits of the funcitons
    yabsmax = np.max(np.abs(yv))
    yabsred = 3.3#0.90*yabsmax
    
    # evaluate cosntraints and find active ones
    gv = bound_constraint_DSS(acfv,fcfv,yabsred,Nsamp=4,scaler=1.0)
    ggvec = gv<0.0
    


    if ggvec.any()==True:
        print ('Violation detected:')
        print ('gvalue = ', gv[ggvec])
        print ('gvalue relative = ', gv[ggvec]/yabsred**2)
    
    fig = plt.figure('Function and Limits')
    ax = fig.add_subplot(111)
    ax.plot(tv,yv,'0.6',linewidth=3)
    ax.plot(tv, yabsred*np.ones(NT),'r',linewidth=1)
    ax.plot(tv,-yabsred*np.ones(NT),'r',linewidth=1)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show()
    
    