'''
Created on 4 Feb 2015

@author: sm6110

Understand minimum control points to reproduce a 1Hz range accurately.

NI=4*T*fmax with cubic splines seems the most robust choice.

'''

import scipy.interpolate as scint
import numpy as np
import matplotlib.pyplot as plt


# reconstructed domain and reference
T = 2  # sec
xv = np.linspace(0,T,500)
fmax=6.0
frange=np.arange(5.0,fmax+1e-3,0.125)

for ff in frange:
    
    yvref = np.sin(2.0*np.pi*ff*xv)
    
    # ----------------------------- spline approximation with 4 points per frequency
    
    # control points
    NI = 4.0*T*fmax 
    xc = np.linspace(0,T,NI+1)
    yc = np.sin(2.0*np.pi*ff*xc)
    
    # spline approx
    yv0=scint.interpolate.spline(xc,yc,xv,order=0)
    yv1=scint.interpolate.spline(xc,yc,xv,order=1)
    yv2=scint.interpolate.spline(xc,yc,xv,order=2)
    yv3=scint.interpolate.spline(xc,yc,xv,order=3)
    
    # frequency domain - numerical 
    fig = plt.figure('Reconstruction of frequency %1.1f Hz with %d control points' %(ff, NI+1))
    ax = fig.add_subplot(111)
    ax.plot(xv,yvref,'r',linewidth=2)
    ax.plot(xv,yv0,'y',linewidth=1)
    ax.plot(xv,yv1,'0.6',linewidth=1)
    ax.plot(xv,yv2,'b',linewidth=1)
    ax.plot(xv,yv3,'k',linewidth=1)
    ax.set_xlabel('s')
    plt.legend(('ref', 'ord 0', 'ord 1', 'ord 2', 'ord 3'))
    fig.savefig('./fig/NI%dsin%2.2fHz.png' %(NI,ff))
    fig.show()   

    # ------------- spline approximation with 3 points per frequency not shifted
    
    # control points
    NI = 3.0*T*fmax # 2 interval minimum per period of oscillation
    xshift=T/NI
    xc = np.linspace(-xshift,T+xshift,NI+2)
    yc = np.sin(2.0*np.pi*ff*xc)
    
    # spline approx
    yv0=scint.interpolate.spline(xc,yc,xv,order=0)
    yv1=scint.interpolate.spline(xc,yc,xv,order=1)
    yv2=scint.interpolate.spline(xc,yc,xv,order=2)
    yv3=scint.interpolate.spline(xc,yc,xv,order=3)
    
    # frequency domain - numerical 
    fig = plt.figure('Reconstruction of frequency %1.1f Hz with %d shifted control points' %(ff, NI+1))
    ax = fig.add_subplot(111)
    ax.plot(xv,yvref,'r',linewidth=2)
    ax.plot(xv,yv0,'y',linewidth=1)
    ax.plot(xv,yv1,'0.6',linewidth=1)
    ax.plot(xv,yv2,'b',linewidth=1)
    ax.plot(xv,yv3,'k',linewidth=1)
    ax.set_xlabel('s')
    ax.set_ylim([-1.5, 1.5])
    plt.legend(('ref', 'ord 0', 'ord 1', 'ord 2', 'ord 3'))
    fig.savefig('./fig/NI%dsin%2.2fHz.png' %(NI,ff))
    fig.show()  


    # ------------- spline approximation with 2 points per frequency but shifted
    
    # control points
    NI = 2.0*T*fmax # 2 interval minimum per period of oscillation
    xshift=1.0/(4.0*fmax)
    xc = np.linspace(-xshift,T+xshift,NI+2)
    yc = np.sin(2.0*np.pi*ff*xc)
    
    # spline approx
    yv0=scint.interpolate.spline(xc,yc,xv,order=0)
    yv1=scint.interpolate.spline(xc,yc,xv,order=1)
    yv2=scint.interpolate.spline(xc,yc,xv,order=2)
    yv3=scint.interpolate.spline(xc,yc,xv,order=3)
    
    # frequency domain - numerical 
    fig = plt.figure('Reconstruction of frequency %1.1f Hz with %d shifted control points' %(ff, NI+1))
    ax = fig.add_subplot(111)
    ax.plot(xv,yvref,'r',linewidth=2)
    ax.plot(xv,yv0,'y',linewidth=1)
    ax.plot(xv,yv1,'0.6',linewidth=1)
    ax.plot(xv,yv2,'b',linewidth=1)
    ax.plot(xv,yv3,'k',linewidth=1)
    ax.set_xlabel('s')
    ax.set_ylim([-1.5, 1.5])
    plt.legend(('ref', 'ord 0', 'ord 1', 'ord 2', 'ord 3'))
    fig.savefig('./fig/NI%dsin%2.2fHz.png' %(NI,ff))
    fig.show() 
    
    
   
    

