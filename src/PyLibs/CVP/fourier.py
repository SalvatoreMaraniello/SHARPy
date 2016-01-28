'''
@package    PyLibs.CVP.fourier
@brief      Contains methods for building fourier based parametrisations
@author     Salvatore Maraniello
@contact    salvatore.maraniello10@imperial.ac.uk 
@version    1.0
@date       15/10/2015 (original 3/11/2014)
@pre        Adapted from python2.7 version in BeamLib/src/optimiser/wrapper/pysrc

@details    The module contains methods for parametrising signals using
            sine and cosines
@warning    Only Discrete Sine Series implemented
                
'''

import numpy as np



def sin_nodal_force(NumNodes, Time, NodeForce, A, F ):
    '''    
    Given a set of amplitudes coefficient and related frequencies, the function
    reconstructs a signal over the time-domain Time created overlapping 
    sinusoidal waves. The signals can be used as a loading component acting on 
    the node of number 'NodeForce'. 
    
    The input A and F are arrays of size (cfmax,6), where cfmax is the total
    number of frequencies used. If A,F are 1D arrays, they are reshaped, if possible, 
    into a cfmax, 6 dimension array. 
      
    Input:
    NumNodes: total number of nodes in the structure.
    A, F: For each of the ii-th component of the force at the node NodeForce, 
    the array A(:,ii) returns the amplitude associated to the sinusoidal waves of
    frequency F(:,ii)
    nodeForce: node where to apply the force. The nodes numbering is such that:
            root: NodeForce =  0
            tip:  NodeForce = -1 or NumNodes-1
            
    Output:
    Fdyn: matrix such that Fdyn[NodeForce,comp,:] contains the loading condition
    for the 'comp' component of the node 'nodeForce'
    '''
    
    NumSteps = len(Time)-1
    if NodeForce > NumNodes-1: raise NameError("NodeForce can't be higher then NumNodes-1!")
    
    Fdyn = np.zeros( (NumNodes,6,NumSteps+1), dtype=float, order='F')
    
    # reshape if 1D array (allows A,F to be 1D) 
    Acp = reshape_cfs(A)    
    Fcp = reshape_cfs(F)
    
    # build F-time matrix
    for comp in range(6):
        Fdyn[NodeForce,comp,:]=sin_series_vec(Time,Acp[:,comp],Fcp[:,comp])
    
    return Fdyn


def DSSvec(Time,av,EvalDeriv=False):
    '''
    As per sin_series_vec, but all frequencies are automatically 
    built based on the number of coefficients in av.
    
    Note that av[0] will correspond to the lower frequency component
    of the signal (and not to a constant one)
    
    If EvalDeriv is True, the time derivative of the signal is
    returned as well
    '''
    
    # build frequencies
    Ncf=len(av)
    f0 = 0.5/Time[-1]
    fv = np.zeros(Ncf)
    for ii in range(Ncf):
        fv[ii]=(ii+1)*f0
    
    y = sin_series_vec(Time,av,fv)
    
    if EvalDeriv==True:
        bv = 2.0*np.pi*fv*av
        dydt=cos_series_vec(Time,bv,fv)
        out=(y,dydt,fv)
    else:
        out=(y,fv)
    
    return out
    

def sin_series_vec(Time, av, fv):
    '''
    Given a time vector Time, and a set of coefficients av and fv, the function
    builds the sines series signal:
    
        y(ii) = sin( 2*pi*Time[ii]*fv[ii] )
    '''

    NumSteps = len(Time)-1
    cfmax = len(fv)
    if cfmax != len(av): raise NameError('av and fv should have same length!')
    
    y = np.zeros(NumSteps+1)

    for jj in range(cfmax):
        y = y + av[jj]*np.sin(2.0*np.pi*fv[jj]*Time)
    
    return y


def cos_series_vec(Time, bv, fv):
    '''
    Given a time vector Time, and a set of coefficients av and fv, the function
    builds the sines series signal:
    
        y(ii) = cos( 2*pi*Time[ii]*fv[ii] )
    '''

    NumSteps = len(Time)-1
    cfmax = len(fv)
    if cfmax != len(bv): raise NameError('bv and fv should have same length!')
    
    y = np.zeros(NumSteps+1)

    for jj in range(cfmax):
        y = y + bv[jj]*np.cos(2.0*np.pi*fv[jj]*Time)
    
    return y


def reshape_cfs(A):
    ''' Given a 1D array of coefficients for the Fourier/Sin representation of
    loads, the function reshapes them into a (cfmax,6) array '''
    # reshape if 1D array (allows A,F to be 1D) 

    if len(A.shape) == 1:
        print('reshaping!!!')
        B=np.array(A,order='F')
        cfmax,rem=len(B)/6,len(B)%6
        if rem!=0:
            raise NameError('If A,F are 1D arrays, their length has to be a multiple of 6')
        A=B.reshape((cfmax,6))
    return A


def fft(t,yv,fcut=np.infty,periodic=True):
    '''
    Returns frequency spectrum of a function yv:
    - if periodic, the input signal is assumed periodic and the inputs t and 
      yv must be such that:
          yv[0]=yv[-1]
          t[-1]=T
      in this case, the last element of yv is ignored.
    - if not periodic, the last element of yv is not discarded.
    where T is the period of simulation. The function returns the complex 
    coefficients of the complex Fourier series associated to the function yv 
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
    if N !=len(yv):
        raise NameError('t and yv must have same length')
 
    ## find complex coefficients
    if periodic==True: cfs = 2.0*np.fft.rfft(yv[:-1],N)/N 
    else: cfs = 2.0*np.fft.rfft(yv,N)/N 

    # and related frequencies...
    fr = fo*np.arange(0,len(cfs)) 
    # nb: this is equivalent to:
    # fr = fs * linspace(0,1,len(cfs))
    # as: fr[nn]=fs*nn/Nout=Nout/T*nn/Nout=nn/T
    # where Nout=len(cfs)
    
    ## cut freq:
    iivec = fr <= fcut

    return fr[iivec], cfs[iivec]


def ifft_bkp(fr,cfs, periodic=True):
    '''
    Inverse of fft function defined in this module, i.e.:
    
    fr, cfs = fft(t,y)
    t, v = ifft(fr, cfs)
    
    '''
    
    f0 = fr[1]-fr[0]#np.average(np.diff(fr))
    T = 1./f0
    if periodic==True: N = len(fr)+1
    else: N = len(fr)
    
    tv = np.linspace(0.0,T,N)
    yv = np.fft.irfft(N/2.0*cfs, N)
    
    return tv, yv


def ifft(fr,cfs, tv, periodic=True):
    '''
    Inverse of fft function defined in this module, i.e.:
    
    fr, cfs = fft(t,y)
    t, v = ifft(fr, cfs)
    
    where tv is the original time vecotr used for fft and
    fr, cfs are the output of fft.
    
    @warning: to be verified if periodic flag is required.
    The function works for aperiodic signals, and seem to
    automatically work for periodic as well!
    
    '''
    
    f0 = fr[1]-fr[0]#np.average(np.diff(fr))
    T = 1./f0
    assert np.abs(T-tv[-1])<1e-6, ('Frequencies spacing anf total time do not match!')
    
    N=len(tv)
    yv = np.fft.irfft(N/2.0*cfs, N)
    
    return  yv


def fft_filter(tv,yv,factor=0.1,periodic=True):
    '''
    removes all harmonics whose amplituce is below the factor 
    of the peak
    '''
    
    # FFT
    fr, cfs = fft(tv,yv,fcut=np.infty,periodic=periodic)
    
    # get harmonics amplitude
    CF = np.abs(cfs)
    
    # and set to zero the small ones
    ccvec = CF < CF.max()
    cfs[ccvec]=0.0
    
    # go back in time domain
    ysmooth = ifft(fr,cfs,tv,periodic)
    
    return ysmooth
    
    

def get_dss(tv,yv,fcut):
    # create yv signal for FFT (with only sines in output)
    yvext = np.concatenate( (yv[:-1] ,-yv[-1::-1]) )  
    tvext   = np.concatenate( (tv[:-1]   , tv+tv[-1] ) )  
    # FFT the extended signal
    fr,cf=fft(tvext,yvext,fcut)   
    # sine series coefficients    
    acf=-np.imag(cf)
    return fr,cf,acf




if __name__=='__main__':
    
    import matplotlib.pyplot as plt    
    import PyLibs.numerics.diff
    
    
    av = np.array([1., 4., 2., 6.])
    tv=np.linspace(0,3,400)
    
    # build DSS with derivatives
    (y, dy) = DSSvec(tv,av,EvalDeriv=True)
    
    # check derivative against numerical
    dynum=PyLibs.numerics.diff.difffun(y,tv)
    
    # plot
    fig=plt.figure('check sin_series_vec')
    ax = fig.add_subplot(211)
    ax.plot(tv,y,'r+')
    
    ax = fig.add_subplot(212)
    ax.plot(tv,dy,'r')
    ax.plot(tv,dynum,'k--')
    
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
    print('frequencies [Hz]: ', fv)
    print('F matrix [hz]', F)
    
    A[0,0]=1.0 # Fx - single freq
    A[1,1]=1.0 # Fy - single freq
    A[2,2]=1.0 # Fz - single freq
    A[3,3]=1.0 # Mx - single freq
    A[:,4]=np.array([0.2, 0.2, 0.2, 0.4])
    A[:,5]=np.array([0.8, 0.1, 0.1, 0.1])
    print('A matrix: ', A)
    
    # check reshape capability if A is 1D
    A2 = np.zeros((4*6,))
    A2[0]=1.0 # Fx - single freq
    A2[7]=1.0 # Fy - single freq
    A2[14]=1.0 # Fz - single freq
    A2[21]=1.0 # Mx - single freq
    A2[4],A2[10],A2[16],A2[22]=0.2,0.2,0.2,0.4
    A2[5],A2[11],A2[17],A2[23]=0.8,0.1,0.1,0.1
    
    
    Fdyn = sin_nodal_force(NumNodes, Time, NodeForce, A2, F )
    
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
        print('Max abs value for node %d is %f' %(ii,M))
    '''
    
    
    
    
    
    
        
        
    
    
    
        
    
    
    