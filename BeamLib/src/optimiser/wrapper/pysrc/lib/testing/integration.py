'''
Created on 24 Oct 2014

@author: sm6110

Test lib.postpr.integrate_rotation methods


'''
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../..')
import shared

import lib.postpr



def compute_quat(fi,nrot):
    
    Xi = np.zeros(4)
    Xi[0], Xi[1:] = np.cos(0.5*fi),  np.sin(0.5*fi)*nrot
    
    return Xi



def OmegaProfile(Time,nrot,profile):
    ''' given a time vector, the function returns an angular velocity time
        profile. The outputs are:
        
        '''
    
    NT = len(Time)-1        
    wt = np.zeros(NT+1)
    Omega = np.zeros((NT+1,3))
    dOmega = np.zeros((NT+1,3))
    
    if profile=='lin':
        wt = np.linspace(0,np.pi,NT+1)
        dwt = np.ones(NT+1)*np.pi
        
    elif profile=='sin':
        wt = np.sin(0.5*np.pi*Time)
        dwt = 0.5*np.pi*np.cos(0.5*np.pi*Time)
    
    elif profile=='exp':
        wt = np.exp(Time)
        dwt = np.exp(Time)
        
    else:
        raise NameError('Profile %s not defined!' %(profile))    
           
    for ii in range(3):
        Omega[:,ii] =nrot[ii]*wt
        dOmega[:,ii]=nrot[ii]*dwt
    
    return Omega, dOmega  
        
# ---------------------------------------------------------------------------- #


if __name__=='__main__':
    
    expmax=5
    eevec=np.arange(expmax)[1:]
    print eevec
    method='1tay'
    # create matrices to store the error
    Er = np.zeros((expmax,3))


    # set a general axis of rotation
    nrot = np.array( [np.cos(np.pi/4.0), np.sin(np.pi/4.0), 0] )

    #----------------------------------------------------------- expected values
    # quaternions after 1"
    
    # lin
    fi_lin = np.pi/2
    Xi_lin = compute_quat(fi_lin, nrot)

    # sin
    fi_sin = 2.0/np.pi
    Xi_sin = compute_quat(fi_sin, nrot)
    
    # exp
    fi_exp = np.e - 1.0
    Xi_exp = compute_quat(fi_exp, nrot)
    
    print 'expected quaternions for...'
    print 'lin:', Xi_lin
    print 'sin:', Xi_sin
    print 'exp:', Xi_exp
      
      
    #--------------------------------------------------------- Convergence Check
    for ee in eevec:

        NT=10**ee
        Time = np.linspace(0,1,NT+1)
        
        # linear case
        Omega, dOmega = OmegaProfile(Time,nrot,profile='lin')
        Xi = lib.postpr.integrate_rotations(Omega, Time, method=method, dOmega=dOmega)
        Er[ee,0]= np.linalg.norm( Xi_lin - Xi[  -1,:] )
        
        # sin case
        Omega, dOmega = OmegaProfile(Time,nrot,profile='sin')
        Xi = lib.postpr.integrate_rotations(Omega, Time, method=method, dOmega=dOmega)
        Er[ee,1]= np.linalg.norm( Xi_sin - Xi[  -1,:] )
        
        # exp case
        Omega, dOmega = OmegaProfile(Time,nrot,profile='exp')
        Xi = lib.postpr.integrate_rotations(Omega, Time, method=method, dOmega=dOmega)
        Er[ee,2]= np.linalg.norm( Xi_exp - Xi[  -1,:] )
    
    

    print Er
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(expmax),Er[:,0],'r')
    ax.plot(np.arange(expmax),Er[:,1],'b')
    ax.plot(np.arange(expmax),Er[:,2],'k')
    ax.set_yscale('log')
    plt.legend(('lin','sin','exp'))
    plt.title('Convergence Rate for "%s" method' %(method))
    fig.savefig('rotconvrate_'+method+'.png')
    plt.show()


