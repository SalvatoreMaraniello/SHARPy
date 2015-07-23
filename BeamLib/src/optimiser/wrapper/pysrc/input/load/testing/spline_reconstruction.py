'''
Created on 23 Jul 2015

@author: sm6110

Test spline reconstruction technique for different function

Script copied in "SHARPy_studies/other/bases/Bsplines" area

'''

import numpy as np
import input.load.spline as sp

import matplotlib.pyplot as plt



#--------------------------------------------------------------- Random Function
tv=np.linspace(0,3,80)
tc=np.linspace(0,3,12)

#fv=np.sin(3*tv)+0.3*tv**2
#fspline, scfv, S =sp.reconstruct(tc,tv,fv)

#plt.plot(tv,fv,'r')
#plt.plot(tv,fspline,'kd')
#plt.show()



#----------------------------------------------- Spline basis: can we go back???
#
# case 01: the control points do not match the domain:
#    two sources of error: those related to the linear interpolation of control
#    and extra points (the extra points location should not affect)
# case 02: the control points coincide with the points in tv. In this case 
#    machine precision accuracy is obtained since regardless the size of tv
# 
Nc=25
ErList=[]
Nv=[50,100,1000,10000]
p=3
Nv=[48,240,2400,24000]

scfv0=np.random.rand(Nc+(p-1))
tc=np.linspace(0,2,Nc)

for nn in Nv:
    tv=np.linspace(0,2,nn+1)
    
    print tv
    print tc

    fv=sp.spline_series_vec(tv, tc, scfv0, p, tc_uniform_flag=True)
    fspline, scfvrec, S =sp.reconstruct(tc,tv,fv)

    print 'Original:', scfv0
    print 'Reconstructed:', scfvrec
    ErList.append( np.linalg.norm(scfvrec-scfv0) )
    
    fig = plt.figure('Spline Nc=%.2d Reconstruction Nv=%d'%(Nc,nn) )
    ax = fig.add_subplot(111)
    ax.plot(tv,fv, color='0.6')
    ax.plot(tv,fspline,'r')



fig = plt.figure('Error')
ax = fig.add_subplot(111)
ax.loglog(np.array(Nv),np.array(ErList),'r')
plt.show()



