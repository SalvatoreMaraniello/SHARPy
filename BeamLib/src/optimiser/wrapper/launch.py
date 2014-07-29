
# ref.1: http://www.walkingrandomly.com/?p=85
# ref.2: http://stackoverflow.com/questions/5811949/call-functions-from-a-shared-fortran-library-in-python
# ref.3: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.ctypes.html
# 
# The routine names can be printed typing (from terminal) 
# nm ./xbeamopt.so (ref.2)
# -----------------------------------------

#from ctypes import byref, cdll, c_int, c_double, c_void_p
import ctypes as ct
import numpy as np

xb = ct.cdll.LoadLibrary('./bin/xbeamopt.so')

# extract main routine
fm=xb.__opt_routine_MOD_opt_main

# prepare input:
NumElem  = ct.c_int(14)
NumNodes = ct.c_int(3)
print 'NumElem = %d, NumNodes = %d' % (NumElem.value, NumNodes.value)

# this creates a _ctype type, which is not really ctype...
#COST =np.array([])
#x.ctypes.data
#C=(c_double*10)()
#COST=c_double(10.8) # this is better, but allocates the memory
#COST=c_double(0.8)
####### call the routine

#######
#COST=c_void_p()
#fm( byref(NumElem), byref(NumNodes), byref(COST) )

#######
#COST=np.array([],dtype=float)
#fm( byref(NumElem), byref(NumNodes), COST.ctypes.data_as(c_void_p) )

#######
#COST=np.empty([],dtype=float,order='F')
#fm( byref(NumElem), byref(NumNodes), COST.ctypes.data_as(c_void_p) )

####### this works if the size exactly matches
#a=np.ones(7)
#COST=np.empty_like(a,dtype=float,order='F')
#fm( byref(NumElem), byref(NumNodes), COST.ctypes.data_as(c_void_p) )
#print COST

####### same as above
#COST=np.empty((7),dtype=float,order='F')
#fm( byref(NumElem), byref(NumNodes), COST.ctypes.data_as(c_void_p) )

###
#COST=np.empty((7),dtype=float,order='F')
#pCOST=ct.pointer()
#fm( ct.byref(NumElem), ct.byref(NumNodes), COST.ctypes.data_as(c_void_p) )


####### same as above
# the allocated array works
### works but numbers are crap. also the size must be the same
###COST=np.empty((7),dtype=float,order='F')
W_COST=np.empty((2),dtype=float,order='F')
fm( ct.byref(NumElem), ct.byref(NumNodes), COST.ctypes.data_as(ct.c_void_p), W_COST.ctypes.data_as(ct.c_void_p) )


print 'NOPT_MAX %d' % (NumElem.value)
print 'Cost Function:', COST 
print 'Cost Weights:', W_COST 



