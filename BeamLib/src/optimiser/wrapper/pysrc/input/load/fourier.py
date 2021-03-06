'''
Created on 3 Nov 2014

@author: sm6110

Global reconstruction of signal using Fourier coefficients

'''

import numpy as np



def glsin_nodal(NumNodes, Time, NodeForce, A, F , optdict={}):
    '''
    Given a set of amplitudes coefficient and related frequencies, the function
    reconstructs a signal over the time-domain Time creates overlapping 
    sinusoidal waves. The signal is used as a loading component acting on 
    the node of number 'NodeForce'. 
    
    The input A and F are arrays of size (cfmax,6), where cfmax is the total
    number of frequencies used. If A,F are 1D arrays, they are reshaped, if possible, 
    into a cfmax,6 dimension array. 
      
    Input:
    NumNodes: total number of nodes in the structure.
    A, F: For each of the ii-th component of the force at the node NodeForce, 
    the array A(:,ii) returns the amplitude associated to the sinusoidal waves of
    frequency F(:,ii)
    nodeForce: node where to apply the force. The nodes numbering is such that:
            root: NodeForce =  0
            tip:  NodeForce = -1 or NumNodes-1
            
    optdict: see 'sin_series_vec' method
    '''
    
    NumSteps = len(Time)-1
    if NodeForce > NumNodes-1: raise NameError("NodeForce can't be higher then NumNodes-1!")
    
    Fdyn = np.zeros( (NumNodes,6,NumSteps+1), dtype=float, order='F')
    
    # reshape if 1D array (allows A,F to be 1D) 
    Acp = reshape_cfs(A)    
    Fcp = reshape_cfs(F)
    
    # build F-time matrix
    for comp in range(6):
        Fdyn[NodeForce,comp,:]=sin_series_vec(Time,Acp[:,comp],Fcp[:,comp], optdict=optdict)
    
    return Fdyn



def sin_series_vec(Time, av, fv, optdict={}):
    '''
    Given a time vector Time, and a set of coefficients av and fv, the function
    builds the sines series signal:
    
        y(ii) = sin( 2*pi*Time[ii]*fv[ii] )
        
    the optdict dictionary allows to specify extra properties:
        - reverse: flip left/right the signal
        - smooth: allows to apply a smoothing factor to avoid non zero signals
        at time t=0 when frequencies that are not in the group fn = 1/(2 Time[-1])
        are used with a reverse option true (see smooth_factor method)
        - Twin: specify the window where to apply the smoothing
    '''

    NumSteps = len(Time)-1
    cfmax = len(fv)
    if cfmax != len(av): raise NameError('av and fv should have same length!')
    
    y = np.zeros(NumSteps+1)

    for jj in range(cfmax):
        y = y + av[jj]*np.sin(2.0*np.pi*fv[jj]*Time)
    
    if 'reverse' in optdict:
        if optdict['reverse'] == True: y = y[::-1]
    if 'smooth' in optdict:
        if optdict['smooth'] == True:
            if 'Twin' in optdict:
                y = y*smooth_factor(Time,optdict['Twin'])
            else:
                raise NameError('Specify time window Twin in which apply the smoothing!')
    
    return y


def reshape_cfs(A):
    ''' Given a 1D array of coefficients for the Fourier/Sin representation of
    loads, the function reshapes them into a (cfmax,6) array '''
    # reshape if 1D array (allows A,F to be 1D) 

    if len(A.shape) == 1:
        print 'reshaping!!!'
        B=np.array(A,order='F')
        cfmax,rem=len(B)/6,len(B)%6
        if rem!=0:
            raise NameError('If A,F are 1D arrays, their length has to be a multiple of 6')
        A=B.reshape((cfmax,6))
    return A


def myfft(t,f,fcut=np.infty):
    '''
    Returns frequency spectrum of a function f. To work properly, the inputs t 
    and f must be such that:
        f[0]=f[-1]
        t[-1]=T
    where T is the period of simulation. The function returns the complex 
    coefficients of the complex Fourier series associated to the function f 
    multiplied by 2, such that:
        cf = ( B - A i)
    where A and B are the coeff.s associted to the sines and cosines. The 
    vector of output frequencies is also given.

    Note:
    if fcut is not specified:
        Nout = Nin/2 + 1 (both if Nin is even or odd)
    otherwise, the output size changes according to the cut frequency
    
    '''
    
    N = len(t)  # signal length
    T=t[-1]-t[0]
    fo = 1./T
    if N !=len(f):
        raise NameError('t and f must have same length')
 
    ## find complex coefficients
    cfs = 2.0*np.fft.rfft(f[:-1],N)/N 

    # and related frequencies...
    fr = fo*np.arange(0,len(cfs)) 
    # nb: this is equivalent to:
    # fr = fs * linspace(0,1,len(cfs))
    # as: fr[nn]=fs*nn/Nout=Nout/T*nn/Nout=nn/T
    # where Nout=len(cfs)
    
    ## cut freq:
    iivec = fr <= fcut

    return fr[iivec], cfs[iivec]


def get_dss(tv,yv,fcut):
    # create yv signal for FFT (with only sines in output)
    yvext = np.concatenate( (yv[:-1] ,-yv[-1::-1]) )  
    tvext   = np.concatenate( (tv[:-1]   , tv+tv[-1] ) )  
    # FFT the extended signal
    fr,cf=myfft(tvext,yvext,fcut)   
    # sine series coefficients    
    acf=-np.imag(cf)
    return fr,cf,acf


def smooth_factor(tv,Twin):
    '''
    Returns a smoothing functions having values:
        0 _> 1 for t<Twin - cubic trend, zero first derivative at 0 and Twin
        1 constant for t > Twin
    Can be used to smooth a DSS series when non zero values occur at t=0    
    '''
    
    pvec = np.zeros(4)
    pvec[0] = -2.0/Twin**3
    pvec[1] = 3.0/Twin**2
    
    fsmooth = np.ones(len(tv))
    iiwin = tv<Twin
    fsmooth[iiwin]=np.polyval(pvec,tv[iiwin])
    
    return fsmooth



if __name__=='__main__':
    
    import matplotlib.pyplot as plt
    
    acfv=np.array([1.7, .5])
    fv=np.array([0.75, 2.45])
    tv=np.linspace(0.0,2.0,1000)
    yv=sin_series_vec(tv, acfv, fv )
    yvrev=sin_series_vec(tv, acfv, fv, {'reverse':True })
    
    #smooth=smooth_factor(tv,Twin=0.05)
    #yvrevsm=smooth*yvrev
    yvrevsm=sin_series_vec(tv, acfv, fv, {'reverse':True,'smooth':True,'Twin':0.05} )
    
    plt.plot(tv,yv,'b')
    plt.plot(tv,yvrev,'r')
    plt.plot(tv,yvrevsm,'k')
    
    plt.show()
    
    
    '''
    
    NodeForce=0
    NumNodes=3
    cfmax=4
    
    A, F = np.zeros((cfmax,6)), np.ones((cfmax,6))
    Time=np.linspace(0,2,100)
    
    # build frequencies
    fo = 0.5/Time[-1]
    fv = np.zeros(cfmax)
    for ii in range(cfmax):
        fv[ii]=(ii+1)*fo
        F[ii,:]=np.ones(6)*fv[ii]    
    print 'frequencies [Hz]: ', fv
    print 'F matrix [hz]', F
    
    A[0,0]=1.0 # Fx - single freq
    A[1,1]=1.0 # Fy - single freq
    A[2,2]=1.0 # Fz - single freq
    A[3,3]=1.0 # Mx - single freq
    A[:,4]=np.array([0.2, 0.2, 0.2, 0.4])
    A[:,5]=np.array([0.8, 0.1, 0.1, 0.1])
    print 'A matrix: ', A
    
    # check reshape capability if A is 1D
    A2 = np.zeros((4*6,))
    A2[0]=1.0 # Fx - single freq
    A2[7]=1.0 # Fy - single freq
    A2[14]=1.0 # Fz - single freq
    A2[21]=1.0 # Mx - single freq
    A2[4],A2[10],A2[16],A2[22]=0.2,0.2,0.2,0.4
    A2[5],A2[11],A2[17],A2[23]=0.8,0.1,0.1,0.1
    
    
    Fdyn = glsin_nodal(NumNodes, Time, NodeForce, A2, F )
    
    
    print '----------------------------------- Check against sin_series_vec ---'
    
    for comp in range(6):
        y = sin_series_vec(Time,A[:,comp],F[:,comp])
        fig=plt.figure('check sin_series_vec')
        ax = fig.add_subplot(111)
        hman,=ax.plot(Time,y,'r+')
        hauto,=ax.plot(Time,Fdyn[NodeForce,comp,:],'k')
        plt.legend([hman, hauto],('manual','auto'))
        #fig.savefig('fig_sin_series_vec_%d.pdf' %(comp))
        plt.show()
    
    # check that other nodes are not allocated    
    for ii in range(NumNodes):
        M = np.max(np.max( np.abs(Fdyn[ii,:,:]) ))
        print 'Max abs value for node %d is %f' %(ii,M)
    
    
    '''
    
    
    
    
        
        
    
    
    
        
    
    
    