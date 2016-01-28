'''

@package    PyLibs.CVP.spline
@brief      Contains methods for building spline based parametrisations
@author     Salvatore Maraniello
@contact    salvatore.maraniello10@imperial.ac.uk 
@version    1.0
@date       15/10/2015 (original 6/02/2015)
@pre        Adapted from python2.7 version in BeamLib/src/optimiser/wrapper/pysrc

@details    The module contains methods for parametrising signals using
            B-splines
            
List of Functions:
------------------
    - zero_base: get spline base of order 0
    - p_basis: gets all the base functions of order p. Tc does not have to be uniformly 
      spaced
    - get_higher_base: using recursive formula, gets the pp order basis functions from 
      the pp-1 basis
    - normalise_uniform_base: normalises the basis by the maximum value of all of them.
        Note that by default, the cases p=2, 3 are treated assuming equally spaced
        knots, thus to use analytically evaluated maxima.
    - extended_knots_vector: given a set of knot points inside the computational
        domain, the function extend them outside the domain to compute boundary 
        basis
    

Methodology:
------------
Splines degree can range from p=0 (piecewise constant) to, theoretically, 
infinity. Note, however, that to cover a domain with NI intervals using p order 
splines, p+NI spline bases will be required. As the number of knots is NK=NI+1,
we need:
(p+NI) - (NI+1)= p-1 extra conditions, meaning that:
    p=0: we can't achieve a representation passing for all the knots points
    p=1: we can find a unique combination of basis such to pass through all the
        points
    p=2: extra conditions, e.g. on the derivatives at the extrema, need to be
        set.

Knots points and intervals will be numbered as below. Extra indices at the sides 
of the domain will need to be defined to number the spline basis having degree 
p>0.  


    -3dt    -2dt    -1dt    0    dt    2dt    3dt    4dt    5dt    T     T+dt    T+2dt   T+3dt        time
    |-------|-------|------||----|-----|------|------|------|------||------|------|-------|------------>
k=  -3      -2      -1      0    1     2      3      4      5      6       7      8       9
i=     -3       -2      -1    0     1      2     3      4      5      6       7      8
 

We also define:
kmax, imax: maximum value for k and i (9 and 8 in the example above)
kmin, imin: minimum value for k and i indices. For both will be:
                kmin = imin = -p
            where p is the order of the spline bases used. (p=-3 in the example
            above)
NK, NI: total number of indices inside the computational domain (in the example
        NK=7, NI=6). Note that:
                NK=kmax-p+1    NI=imax-p+1
NS: total number of spline functions, equal to p+NI

Each spline base is numbered starting from the knot at which becomes non-zero.
The first spline base will have index -p, while the last will have index imax.
Thus:
    NS = p + NI

         
'''

import numpy as np
from warnings import warn

import PyLibs.numerics.diff




def nodal_force(NumNodes, Time, NodeForce, TCint, Scf, p ):
    '''
    Given a set of amplitudes coefficient and spline basis, the function
    reconstructs a signal over the time-domain Time. The signal is used as a 
    loading component acting on the node of number 'NodeForce'. 
    
    The input A and F are arrays of size (NS,6), where NS = p + NI is the total
    number of spline used. Each column contains the control points internal to
    the time domain and the related splines coefficients.
      
    Input:
    NumNodes: total number of nodes in the structure.
    A, F: For each of the ii-th component of the force at the node NodeForce, 
    the array A(:,ii) returns the amplitude associated to the sinusoidal waves of
    frequency F(:,ii)
    nodeForce: node where to apply the force. The nodes numbering is such that:
            root: NodeForce =  0
            tip:  NodeForce = -1 or NumNodes-1
            
    Remark:
    - consistency checks are performed inside spline_series_vec and subfunctions
    '''
    
    NumSteps = len(Time)-1
    if NodeForce > NumNodes-1: raise NameError("NodeForce can't be higher then NumNodes-1!")
    
    Fdyn = np.zeros( (NumSteps+1,NumNodes,6), dtype=float, order='F')
     
    # build F-time matrix
    for comp in range(6):
        if np.max(np.abs(Scf[:,comp]))>0.0:
            Fdyn[:,NodeForce,comp]=build_vec(Time,TCint[:,comp],Scf[:,comp],p, True)
    
    return Fdyn


def change_series_base(tcint, p, tcint_old, Scf_old, p_old, tc_uniform_flag=True):
    '''Given the coefficients of a spline series Scf_old, of order p_old and over
    the internal knots point tcint_old, the coefficients of a new spline series
    Scf, or order p and over the internal knots point tcint, are computed.
    
    Notes:
    - uniform_flag=True: assumes uniform knots when normalising.
    '''

    # get new spline series size
    NK = len( tcint )
    NI = NK-1
    NS = NI+p    
    
    # domain of interpolation
    tv=np.linspace(tcint[0],tcint[-1],NS)
    
    '''# more robust method????
    tv0=np.linspace(tcint[0],tcint[-1],NI+1)
    tv=np.zeros((NS,)
    dtm=0.5*(tv[1]-tv[0])
    for pp in range(p-1):
        if np.mod(pp,2):
            tvextra[pp]=(2.0*pp+1.)*dtm
        else:
            tvextra[pp]=tv[-1]-( pp*dtm )    
    tv=np.concatenate((tv,tvextra))
    '''
    
    # evaluate old spline on the new domain
    fv = build_vec(tv, tcint_old, Scf_old, p_old, tc_uniform_flag)
    
    # Build new series basis over tv
    tc = extended_knots_vector(tcint,p,tv)
    Smat=p_basis(tv,tc,p)
    Smat = normalise_uniform_base(p,Smat,tc_uniform_flag)
    
    # solve liner system:
    Scf=np.linalg.solve(Smat.transpose(), fv)
    
    return Scf
    

def build_vec(tv, tcint, Scf, p, tc_uniform_flag=True, EvalDeriv=False):
    '''
    Given a time vector tv, and a set of knots points tc with value Scf, the 
    function a p-order spline. 
    
    Assumptions:
    - tc[0]=tv[0] and tc[-1]=tv[-1]
    - the function automatically extends the knot points domain using
    extended_knots_vector
    '''

    # check base number and spline order
    NT=len(tv)
    NK = len( tcint )
    NI = NK-1
    NS = NI+p
    if len(Scf) != NS:
        raise NameError('Scf must have length %d!' %(NS))

    # extend domain
    tcext = extended_knots_vector(tcint, p, tv)

    # get spline base
    S = p_basis(tv,tcext,p)
    
    # normalise base
    S = normalise_uniform_base(p,S,tc_uniform_flag)
    
    # build function
    yv = np.sum( S.transpose()*Scf ,1)

    ### debugging
    #fv0 = np.zeros((NT,))
    #for ss in range(NS):
    #    fv0=fv0+S[ss,:]*Scf[ss]   
    #fig = plt.figure('comparison built function')
    #ax = fig.add_subplot(111)
    #ax.plot(tv,fv,'0.6',linewidth=3)
    #ax.plot(tv,fv0,'r')
    #ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    #ax.set_ylim([-0.5, 1.5])
    #plt.legend(('fast', 'slow'))
    #plt.show()
    
    if EvalDeriv==True:
        dyv = PyLibs.numerics.diff.difffun(yv,tv)
        out=(yv,dyv)
    else:
        out=yv

    return out


def zero_base(tv,tc,s):
    ''' 
    Computes the zero order s-th base function on the domain tv given the
    knot points tc. The control points vector will be such that:
        tc[0]=tv[0]    tc[-1]=tv[-1]
    Note that it is assumed that:
        yv = 1 if   tc[k] <= t < tc[k+1]
    with the exception of the last base, for which 
        yv = 1 if   tc[k] <= t < tc[k+1]
    
    '''
    
    # check domain
    if tc[0]!=tv[0] or tc[-1]!=tv[-1]:
        print('tc[0]=%f, tc[-1]=%f, tv[0]=%f, tv[-1]=%f' %(tc[0],tc[-1],tv[0],tv[-1]))
        raise NameError('tc[0] and tc[-1] must be equal to tv[0] and tv[-1]')
    
    # check base number
    NK = len(tc)
    NI = NK-1
    if s>NI:
        raise NameError('Index for zero order %d-th base function has to be smaller then %d' %(s,NK))
    
    # find non-zero domain
    k=s
    #yv = np.float64( (tv<=tc[k+1]) * (tv>=tc[k]) )
    yv = np.float64( (tv<tc[k+1]) * (tv>=tc[k]) )
    
    if k==NK-2:
        yv[tv==tc[k+1]] = 1.0
    

    return yv


def p_basis(tv,tc,p):
    ''' 
    Computes the p-th order basis functions on the domain tv given the
    knot points tc. The control points vector will be such that:
        tc[0]<=tv[0]    tc[-1]>=tv[-1]
    where the number of control points verifying:
        tc[0]<tv[0]
    is equal to p, order of the spline base 
    '''

    # check domain
    if tc[p]!=tv[0] or tc[-(1+p)]!=tv[-1]:
        print('tc[p]=%f, tc[-(1+p)]=%f, tv[0]=%f, tv[-1]=%f' %(tc[p],tc[-(1+p)],tv[0],tv[-1]))
        raise NameError('tc[p] and tc[-p] must be equal to tv[0] and tv[-1]')
    
    # check base number and spline order
    NT=len(tv)
    tc0 = tc[ (tc>=tv[0]) * (tc<=tv[-1]) ]
    NK = len( tc0 )
    NI = NK-1
    NS = NI+p

    # build S0 base over tv 
    NS0=NI
    S0 = np.empty((NS0,NT))
    for ss in range(NS0):
        S0[ss,:]=zero_base(tv,tc0,ss)
    if p==0:
        return S0
    
    # build higher base:
    Sold = np.zeros((NS,NT))
    Sold[p:] = S0
    for pp in range(1,p+1):
        #print 'computing bases %d' %(pp)
        
        S = get_higher_base(tv,tc,pp,Sold)
        Sold=S.copy()

    return S


def get_higher_base(tv,tc,pp,Slow,debug=False):
    '''
    computes all the basis of the pp-th order splines on the domain tv given the 
    (pp-1)-th basis function Slow over the domain tv
    
    Note that the size of Slow and S depends on the final order p = NS-NI
    '''

    # check base number and spline order
    NS, NT = Slow.shape
    NK = len( tc[ (tc>=0.0) * (tc<=tv[-1]) ] )
    NI = NK-1
    
    #NS=p+NI # number of splines in this base
    p = NS-NI
    pp0 = p-pp
    
    S = np.zeros((NS,NT))
        
    for ss in range(pp0,NS-1):
        
        alpha_1 = ( tv - tc[ss]    ) / (tc[ss+pp]   - tc[ss]  )
        alpha_2 = (-tv + tc[ss+pp+1]) / (tc[ss+pp+1] - tc[ss+1])
        S[ss,:] = alpha_1 * Slow[ss,:] +  alpha_2 * Slow[ss+1,:]
        
        if debug==True:
            import matplotlib.pyplot as plt
            print('base %d order No. %d' %(pp,ss))
            fig = plt.figure('building %d-th base of order %d' %(ss,pp))
            ax = fig.add_subplot(111)
            ax.plot(tv,S[ss,:],'r')
            ax.plot(tv,Slow[ss,:],'0.6')
            ax.plot(tv,Slow[ss+1,:],'k')
            ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
            ax.set_ylim([-0.5, 1.5])
            plt.legend(('%d-th ord.%d' %(ss,pp),  '%d-th ord.%d' %(ss,pp-1), '%d-th ord.%d' %(ss+1,pp-1) ))
            plt.show()  
        
    # last base:
    alpha_1 = ( tv - tc[NS-1]    ) / (tc[NS-1+pp]   - tc[NS-1]  )
    S[NS-1,:] = alpha_1 * Slow[NS-1,:]
        
    return S


def normalise_uniform_base(p,S,tc_uniform_flag=True):
    ''' 
    Given a set of basis functions, the functions are normalised such that:
        - basis defined inside the domain have maximum equal to 1
        - basis at the boundary of the domain will have maximum < 1
        
    WARNING: by default, the function assumes that the basis S are obtained 
    using equally spaced knots [tc_uniform] allowing precise scaling for p=2,3
    order splines
    '''
    
    if tc_uniform_flag==True:
        if p==0 or p==1:
            pass
        if p==2:
            S = 4.0/3.0*S
        if p==3:
            S = 3.0/2.0*S
        else:
            S = (1.0/np.max(S))*S
    else:
        S = (1.0/np.max(S))*S
        
    return S
              
  
def extended_knots_vector(tcint,p,tv):
    ''' Given a set of knots points inside the computational domain and such 
     that:
         tcint[0]=tv[0] & tcint[-1]=tv[-1]
     the function returns the extended set of knot points:
         tcext = [ -p*dt_left.... tcint[0] .... tcint[-1] .... p*dt_right ]
     where dt_left and dt_right are taken to be the first and last interval length
     for tcint vector not uniformely spaced.
     
     Note: tv is required only for checking
    '''       
     
    # check
    if tcint[0]!=tv[0] or tcint[-1]!=tv[-1]:
        print('tc[0]=%f, tc[-1]=%f, tv[0]=%f, tv[-1]=%f' %(tcint[0],tcint[-1],tv[0],tv[-1]))
        raise NameError('tcint[0] and tcint[-1] must be equal to tv[0] and tv[-1]')
    
    if p==0:
        return tcint
    
    # get deltas
    dt_L = tcint[ 1]  -tcint[ 0]
    dt_R = tcint[-1] - tcint[-2]
     
    NCint = len(tcint)
    tcext = np.empty((2.0*p + NCint,))
    tcext[p:-p]=tcint
     
    for pp in range(p):
        tcext[  pp   ]=         -(p-pp)*dt_L
        tcext[-(pp+1)]=tcint[-1]+(p-pp)*dt_R

    return tcext
        
        
def reconstruct(tcint,tv,fv,p=3,method='zero'):
    '''
    Given a function (tv,fv) and a set of control points tcint, such that
        tcint[0]=tv[0]
        tcint[-1]=tv[-1]
    the function returns a spline basis and coefficients of order p that 
    reconstruct the sought signal tv,fv
    
    The method is very accurate if the tcint set is included in the tv set.
    Otherwise the accuracy will depend on the interpolation error. A linear 
    interpolation was used under the assumption that len(tv)>>len(tcint) 
    '''
    
    if p!=3:
        raise NameError('function only developed for p=3!!!')
    Nc=len(tcint)
    Ns=Nc+p-1
    
    Nv=len(tv)
    if Nv<5*Nc:
        warn('Reconstruction may be inaccurate if interpolation is used !!!')
    
    # get control points value of f
    fc=0.0*tcint
    for ii in range(Nc):
        tt=tcint[ii]
        try:
            iivec=tt==tv
            fc[ii]=fv[iivec]
        except:
            print("tcint[%1.2d]=%f couldn't be found exactly in tv!"%(ii,tcint[ii]))
            warn("The value will be found by linear interpolation:")
            fc[ii]=np.interp(tcint[ii],tv,fv)
            '''
            # find by minimum tolerance
            tol=1e-6*tv[-1]
            er=(tv-tcint[ii])**2
            minerror=er.min()
            if minerror<tol**2:
                mm=er.argmin()
                fc[ii]=fv[mm]
            '''
    

    # get the set of basis functions over the tcint domain plus extra points:
    if method=='zero':
        tcrec=np.zeros((Nc+(p-1)))
        fcrec=np.zeros((Nc+(p-1)))
        tcrec[0],fcrec[0]   = tcint[0],fc[0]
        tcrec[2:-2],fcrec[2:-2]= tcint[1:-1],fc[1:-1]
        tcrec[-1],fcrec[-1]=tcint[-1],fc[-1]
        
        # ---------------------------------------------------------------------
        # determine remaining points from tv,fv to be intermediate between:
        # tcint[ 0] and tcint[ 1]
        # tcint[-2] and tcint[-1]
            #
        # Inportant: this part tries hard to avoid interpolation (as we can 
        # are free to choose any extra node. Nodes shouldn't even be at the 
        # extrema. 
        # As the algorythm may fail for small discretisations, in this case
        # the extra points are hard-coded
         
        if Nv<5*Nc:
            tcrec[1],fcrec[1]=tv[1],fv[1]
            tcrec[-2],fcrec[-2]=tv[-2],fv[-2]
        else:                                                
            iivec=tv<tcint[1]
            iiextra=np.round(.5*np.sum(iivec)) 
            tcrec[1],fcrec[1]=tv[iiextra],fv[iiextra]
            
            iivec=tv<tcint[-2]
            iiextra=np.round(0.5*(Nv+np.sum(iivec))) - 1 # cause we start counting from zero
            tcrec[-2],fcrec[-2]=tv[iiextra],fv[iiextra]
        # ---------------------------------------------------------------------
    
        # get the set of basis functions over the final domain
        tc=extended_knots_vector(tcint,p,tcrec)
        Srec = p_basis(tcrec,tc,p)
        Srec = normalise_uniform_base(p,Srec,tc_uniform_flag=True)
        
        # solve for basis coefficients
        scfv = np.linalg.solve(Srec.transpose(),fcrec)

        # build higher definition basis function
        S = p_basis(tv,tc,p)
        S = normalise_uniform_base(p,S,tc_uniform_flag=True)       
    
    
    
    if method=='spectral':
        import PyLibs.CVP.fourier as fourier
        
        if fv[0]!=fv[-1] and fv[0]!=0.0:
            raise NameError('Function has to be symmetric')
        
        #------------------------------------ # get dss representation of signal
        dtc=tcint[1]-tcint[0]
        fmax=1.0/(2.0*dtc)
        frv,cfv,acfv=fourier.get_dss(tv, fv, fcut=1.20*fmax)
        # adjust length of frv
        frv,acfv=frv[1:Ns+1],acfv[1:Ns+1] # no sine at 0 Hz
        fcut=frv[-1]+np.finfo(frv[-1]).eps 
        
        #------------------------------------ # compute FFT of each spline basis
        # build higher definition basis function
        tc=extended_knots_vector(tcint,p,tv)
        S = p_basis(tv,tc,p)
        S = normalise_uniform_base(p,S,tc_uniform_flag=True)             
        #print S.shape
    
        AcfMat=np.zeros((Nc+p-1,len(frv)))
        for ii in range(Nc+p-1):
            frv,cfv,acfsp=fourier.get_dss(tv, S[ii,:], fcut=fcut)
            AcfMat[ii,:]=acfsp[1:]

        # solve!
        scfv = np.linalg.solve(AcfMat.transpose(),acfv)
        
    '''
    # beuild higher definition basis function
    S = p_basis(tv,tc,p)
    S = normalise_uniform_base(p,S,tc_uniform_flag=True)    
    '''
        
    # build function
    fspline = np.sum( S.transpose()*scfv ,1)    

    return fspline, scfv, S
    
    
    



   

    
    
    
