'''
Created on 4 Feb 2015

@author: sm6110

Test spline reconstruction with control of maximum value and maximum derivative.

See spline range for assessment of number of control points and spline order 
chosen

'''

import scipy.interpolate as scint
from scipy.interpolate import InterpolatedUnivariateSpline

import numpy as np
import matplotlib.pyplot as plt

fmax=2. # Hz
T = 2.  # s
fnow = 1.37



# reconstructed domain and reference
xv = np.linspace(0,T,500)
yvref = np.sin(2.0*np.pi*fnow*xv)

# ----------------------------- spline approximation with 4 points per frequency

# control points
NI = 4.0*T*fmax 
xc = np.linspace(0,T,NI+1)

##dxc=T/NI
##xc = np.linspace(-dxc,T+dxc,NI+3)

yc = np.sin(2.0*np.pi*fnow*xc)

# spline approx - 3rd order
yv3=scint.interpolate.spline(xc,yc,xv,order=3)

# find base and reconstruct using the base coefficients
# no smoothing is applied (see documentation)
# with order=3 we get same results as above.
ord=3
spl = InterpolatedUnivariateSpline(xc, yc,k=ord)
yspl=spl(xv)
Aspl=spl.get_coeffs()



# find the bases function and reconstruct
# we reconstruct the spline obtained above with the function for generating a
# spline representation from a base
import sys
sys.path.append('../../..')
import input.load.spline as sp
ymine, B = sp.bspleval(x=xv, knots=xc, coeffs=Aspl, order=ord, debug=False)



# frequency domain - numerical 
fig = plt.figure('Interpolation of frequency %1.1f Hz with %d control points' %(fmax, NI+1))
ax = fig.add_subplot(111)
ax.plot(xv,yvref,'r',linewidth=2)
ax.plot(xv,yv3,'k',linewidth=1)
ax.plot(xv,yspl,'b',linewidth=1)
ax.plot(xv,ymine,'0.6',linewidth=2)
ax.set_xlabel('s')
plt.legend(('ref', 'ord 3','univ spline','univ mine'))
fig.show()  






