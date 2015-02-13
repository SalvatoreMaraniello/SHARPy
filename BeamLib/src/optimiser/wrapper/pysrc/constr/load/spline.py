'''
Created on 12 Feb 2015

@author: sm6110

Method to compute bound and boundary constraints for spline parametrised control

'''

import numpy as np



def boundary_constraint(Scf,p,F0,FT,scaler=1.0):
    '''
    Returns boundary constraints for the initial and final value of the control
    p=3 tested
        self.driver.add_constraint('0.25*xbeam.Scf[0,4]+xbeam.Scf[1,4]+0.25*xbeam.Scf[2,4]=0.0') # Ttq[0]=0.0 
        self.driver.add_constraint('0.25*xbeam.Scf[%d,4]+xbeam.Scf[%d,4]+0.25*xbeam.Scf[%d,4]=0.0' %(NS-3,NS-2,NS-1)) # Ttq[-1]=0.0 
        
    gv(0) = scaler(f(0) - F0)
    gv(T) = scaler(f(T) - FT)
    '''


    NS = len(Scf)
    gv = np.empty((2,))

    if p==3:
        gv[0] = scaler*(0.25*Scf[0] +Scf[1]   +0.25*Scf[2]    - F0)
        gv[1] = scaler*(0.25*Scf[NS-3]+Scf[NS-2]+0.25*Scf[NS-1] - FT)                 
    elif p==2:
        gv[0] = scaler*(0.666666667*(Scf[0] +Scf[1]) - F0)
        gv[1] = scaler*(0.666666667*(Scf[NS-2] +Scf[NS-1]) - FT)
    else:
        raise NameError('Spline order must be either p = 2 or p = 3.')
        
    return gv


def bound_constraint(Scf,p,Fmax,scaler=1.0):
    ''' returns an array with all the constraint necessary to impose a boundary 
    to the control. The constraints are all of the form:
    
    scaler**2 * (q(Scf) - Fmax**2) < 0
    
    where q is a quadratic form returning a positive value. The equations are
    valid only for equally spaced control
    
    Accuracy:
    p = 2: control within 1.02 %
    p = 3: control within 3.92 % 
    
    Remark: 
    - at the boundary some constraints may be redundant. Not accounted for 
    to keep the code simple
    '''

    NS = len(Scf)
    
    if p == 3:
        Gb = np.zeros((NS-1,3)) # splinesNS-1-th not required
        Gc = np.zeros((NS-2,1)) # splines NS-2, NS-1-th not required
        
        for ii in range(0,NS-1):                    
            # conditions b1 and b1symmetry
            Gb[ii-1,0]=     Scf[ii] + 0.25*Scf[ii+1]
            Gb[ii-1,1]=0.25*Scf[ii] +      Scf[ii+1]
            # condition b2
            Gb[ii-1,2]=0.718750*(Scf[ii]+Scf[ii+1])
        
        for ii in range(0,NS-2):
            Gc[ii-1,0]=0.25*Scf[ii] + Scf[ii+1] + 0.25*Scf[ii+2]
        
    elif p == 2:
        Gb = np.zeros((NS-2,1)) # splines NS-1-th not required
        Gc = np.zeros((NS-3,1)) # splines NS-2, NS-1-th not required

        for ii in range(0,NS-1):                    
            # condition b
            Gb[ii-1,0]=0.666666667*( Scf[ii] + Scf[ii+1] )

        for ii in range(0,NS-2):
            Gc[ii-1,0]=0.166666667*Scf[ii] + Scf[ii+1] + 0.166666667*Scf[ii+2]       
    
    else:
        raise NameError('Spline order must be either p = 2 or p = 3.')
    
    # compute constraint, scale and square
    Gb=scaler**2 * ( Gb**2 - Fmax**2)
    Gc=scaler**2 * ( Gc**2 - Fmax**2)
    
    # arrange column by column, i.e. gv = [G[:,0] G[:,1] ... ]
    gb = Gb.transpose().reshape(( Gb.shape[0]*Gb.shape[1] ,))
    gc = Gc.transpose().reshape(( Gc.shape[0]*Gc.shape[1] ,))
    
    gv = np.concatenate((gb,gc))
    
    return gv





if __name__=='__main__':
    
    import input.load.spline as sp
    import matplotlib.pyplot as plt
    
    T=1.0
    NT=101
    p=2
    NI=22
    NS=p+NI
    tv=np.linspace(0,T,NT)
    tcint = np.linspace(0,T,NI+1)  
    
    Scf= -np.random.rand(NS)
    fv = sp.spline_series_vec(tv, tcint, Scf, p, True)  
    
    vabsmax = np.max(np.abs(fv))
    vabsred = 0.95*vabsmax
    
    gv = sp.bound_constraint(Scf, p, vabsred, scaler=1.0)
    
    Ngc = NS - (5-p)
    gc = gv[-Ngc:]
    gb = gv[:(len(gv)-len(gc))] 
    
    # find constraints that detected the violation
    bbvec = gb>0.0
    ccvec = gc>0.0

    if bbvec.any()==True:
        print 'gb detected:'
        print 'gb = ', gb[bbvec]
    if ccvec.any()==True:
        print 'gc detected:'
        print 'gc = ', gc[ccvec]


    fig = plt.figure('Function and Limits')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.plot(tv, vabsred*np.ones(NT),'r',linewidth=1)
    ax.plot(tv,-vabsred*np.ones(NT),'r',linewidth=1)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show()
    
    
    # check initial and final time constraints
    print 'Boundary Costraints check:'
    F0, FT = fv[0], fv[-1]
    geq = sp.boundary_constraint(Scf, p, F0, FT, scaler=1.0)
    print 'Satisfied constraints: %f, %f' %(geq[0], geq[1])
    geq = sp.boundary_constraint(Scf, p, 0.9*F0, 1.2*FT, scaler=1.0)
    print 'Unsatisfied constraints: %f, %f' %(geq[0], geq[1])