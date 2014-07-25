
# ref.1: http://www.walkingrandomly.com/?p=85
# ref.2: http://stackoverflow.com/questions/5811949/call-functions-from-a-shared-fortran-library-in-python
# ref.3: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.ctypes.html
# 
# The routine names can be printed typing (from terminal) 
# nm ./xbeamopt.so (ref.2)
# -----------------------------------------

from ctypes import byref, cdll, c_int

xb = cdll.LoadLibrary('./bin/xbeamopt.so')

# extract main routine
fm=xb.__opt_routine_MOD_opt_main

# prepare input:
NumElem  = c_int(14)
NumNodes = c_int(3)
NOPTMAX = c_int(0)
print 'NumElem = %d, NumNodes = %d' % (NumElem.value, NumNodes.value)


#COST =array([0,0,0,0,0,0,0])
#x.ctypes.data



# call the routine
fm( byref(NumElem), byref(NumNodes), byref(NOPTMAX) )


print 'NOPT_MAX %d' % (NOPTMAX.value)



