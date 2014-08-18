'''
Salvatore Maraniello, 14 Aug 2014

Tests for computing spline derivatives in respect to control points.

Being the spline reconstruction linear, the derivative in respect to the control 
points are the Basis spline functions themselves

The ii-th base function, for any order/number of control points, can be computed
as the spline interpolation of the set
    yc[jj] = 1 (if jj=ii)
    yc[jj] = 0 (otherwise)

The approach works for both polynomial and non functions (as, for each given
set of control points, there exist a polynomial passing through them), randomly
spaced control points.

Important:
The boundaries of the domain must be included inside the control points to 
achieve a global reconstruction.

'''


import numpy as np
import scipy.interpolate as scint
import matplotlib.pyplot as mp


# define functions to fit:
def fpar(x,a,b,c):
    y = a*x**2 + b*x + c
    return y         

def randomfun(x):
    y = 2.0*np.log(x+0.1) - np.exp(x-0.5)
    return y

def get_derivative(xc,yc,y,NC):
    # numerical derivatives for each point of the spline in respect to the 
    # control point
    dy = np.zeros((100,NC))
    for ii in range(NC):
        yfd = yc.copy()
        yfd[ii] = yc[ii]*(1.0+1e-6)
        yup = scint.interpolate.spline(xc,yfd,xref,order=ord)
        dy[:,ii]=(yup-y)/(yc[ii]*(1e-6))  
    return dy  



#----------------------------------------------------------------------- execute
NC, ord = 4, 3
xref = np.linspace(0.,1.,100)           # equally spaced
xc = np.array(np.linspace(0.,1.,NC))

xc = np.array([0.0, 0.3, 0.5, 1.0])     # random spacing
NC=len(xc)


#------------------------------------------------- define reference and spline 1
a,b,c = -1., 0.5, 2.0
yref = fpar(xref,a,b,c)
yc = fpar(xc, a,b,c)
y = scint.interpolate.spline(xc,yc,xref,order=ord)
mp.plot(xref,yref,'r') 
mp.plot(xref,y,'y')
mp.plot(xc,yc,'ko')
mp.legend(('ref','spline','knots')) 
mp.title('II ord polynomial')
mp.show()

# get derivative spline 01
dy1=get_derivative(xc,yc,y,NC)

for ii in range(NC):
    mp.plot(xref,dy1[:,ii])
mp.title('Derivatives in respect to control points. I case')
mp.show()

#------------------------------------------------------------ define spline No.2
a,b,c = 0.0, 0.0, 1.0
yref = fpar(xref,a,b,c)
yc = fpar(xc, a,b,c)
y = scint.interpolate.spline(xc,yc,xref,order=ord)
mp.plot(xref,yref,'r') 
mp.plot(xref,y,'y') 
mp.plot(xc,yc,'ko')
mp.legend(('ref','spline','knots')) 
mp.title('II ord polynomial - different coeff.s')
mp.show()

# get derivative spline 02
dy2=get_derivative(xc,yc,y,NC)

for ii in range(NC):
    mp.plot(xref,dy2[:,ii])
mp.title('Derivatives in respect to control points. II case')
mp.show()


#------------------------------------------------------------ define spline No.3
yref = randomfun(xref)
yc = randomfun(xc)
y = scint.interpolate.spline(xc,yc,xref,order=ord)
mp.plot(xref,yref,'r') 
mp.plot(xref,y,'y') 
mp.plot(xc,yc,'ko')
mp.legend(('ref','spline','knots')) 
mp.title('Random non-linear function')
mp.show()

# get derivative spline 03
dy3=get_derivative(xc,yc,y,NC)

for ii in range(NC):
    mp.plot(xref,dy3[:,ii])
mp.title('Derivatives in respect to control points. III case')
mp.show()


#------------------------------------------------------------------- Final Check
e1=np.max(np.abs(dy1-dy2))
e2=np.max(np.abs(dy2-dy3))
emax = np.max([e1,e2])

print 'Max error:'
print emax



#----------------------------------------------------------- Get Basis Functions
Base = np.zeros((100,NC))
for ii in range(NC):
    yc=np.zeros((NC))
    yc[ii]=1.0
    Base[:,ii] = scint.interpolate.spline(xc,yc,xref,order=ord)

for ii in range(NC):
    mp.plot(xref,Base[:,ii])
mp.title('Analytically Computed Basis (derivatives)')
mp.show()

# The Basis Functions are the derivatives in respect to the control points.
e_analytical = np.max(np.abs(Base-dy1))
print 'Error with Analytical Calculation: ', e_analytical
    
    
    




# Analytical derivatives:


'''
doesn't work with arrays

# interp test with multiple functions
a,b,c = -1., 0.5, 2.0
yref1 = fpar(xref,a,b,c)
yc1 = fpar(xc,a,b,c)
yref2 = randomfun(xref)
yc2 = randomfun(xc)
#Yref = np.array([yref1,yref2])
YC = np.array([yc1,yc2])
XC = np.array([xc,xc])
XREF = np.array([xref,xref])

Y = scint.interpolate.spline(XC,YC,XREF,order=3)

mp.plot(xref,yref1,'k')
mp.plot(xref,yref2,'b')
for ii in range(2):
    mp.plot(xref,Yref[:,ii])
mp.title('interplate multiple function contemporary')
mp.show()


'''









