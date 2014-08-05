#
# Salvatore Maraniello 31 July 2014
# launcher of beam solver
#
# ref.1: http://www.walkingrandomly.com/?p=85
# ref.2: http://stackoverflow.com/questions/5811949/call-functions-from-a-shared-fortran-library-in-python
# ref.3: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.ctypes.html
# 
# ps: to print the library names, type (from terminal): 
# nm ./xbeamopt.so (ref.2)
# ---------------------------------------------------------------------------------------------------------------

import ctypes as ct
import numpy as np

# Load Dynamic Library
xb = ct.cdll.LoadLibrary('./bin/xbeamopt.so')

# extract main routine
fm=xb.__opt_routine_MOD_opt_main

# Options - Default Values:
# overwritten by input_file.py
FollowerForce   = True
FollowerForceRig= True   
PrintInfo       = True    
OutInBframe     = True   
OutInaframe     = False   
ElemProj     = 0       
MaxIterations=99         
NumLoadSteps=5         
#NumGauss=1 # not an input                
Solution=111              
DeltaCurved=1e-5         
MinDelta=1e-8            
NewmarkDamp=1.e-4      


# ------------------------------------------------------------------------- import all input
# python './input_file.py'
from input_file import *


# ------------------------------------------------------------------------- prepare input:
# arrays do not require any modification
NumElems     = ct.c_int( NumElems     )
NumNodes     = ct.c_int( NumNodes     )
NumNodesElem = ct.c_int( NumNodesElem )

ElemType     = ct.c_char_p( ElemType.ljust(4) ) # the ljust is not required
TestCase     = ct.c_char_p( TestCase.ljust(4) )
BConds       = ct.c_char_p( BConds.ljust(2)   )

BeamLength1 = ct.c_double( BeamLength1 )
BeamLength2 = ct.c_double( BeamLength2 )
ThetaRoot   = ct.c_double( ThetaRoot   )
ThetaTip    = ct.c_double( ThetaTip    )
TipMass     = ct.c_double( TipMass     )
TipMassY    = ct.c_double( TipMassY    )
TipMassZ    = ct.c_double( TipMassZ    )
SectWidth   = ct.c_double( SectWidth   )
SectHeight  = ct.c_double( SectHeight  )
Omega       = ct.c_double( Omega       )


# ------------------------------------------------------------------------- Options

FollowerForce    = ct.c_bool( FollowerForce    )
FollowerForceRig = ct.c_bool( FollowerForceRig )  
PrintInfo        = ct.c_bool( PrintInfo        )
OutInBframe      = ct.c_bool( OutInBframe      )
OutInaframe      = ct.c_bool( OutInaframe      )   

ElemProj     = ct.c_int( ElemProj      )      
MaxIterations= ct.c_int( MaxIterations )         
NumLoadSteps = ct.c_int( NumLoadSteps  )   
#NumGauss     = ct.c_int( NumGauss      )          
Solution     = ct.c_int( Solution      )  
            
DeltaCurved  = ct.c_double( DeltaCurved )       
MinDelta     = ct.c_double( MinDelta    )      
NewmarkDamp  = ct.c_double( NewmarkDamp )



# ------------------------------------------------------------------------- Options
print 'Values in input:'
print 'DeltaCurved %f' %(DeltaCurved.value)
print 'MinDelta    %f' %(MinDelta.value)
print 'NewmarkDamp %f' %(NewmarkDamp.value)


print 'ElemType: ',  ElemType.value
print 'TestCase: ',  TestCase.value
print 'BConds: '  ,  BConds.value

print 'Beamlength1: ', BeamLength1.value
print 'Beamlength2: ', BeamLength2.value

# ---------------------------------------------------------------------------------------
W_COST=np.empty((2),dtype=float,order='F')



fm ( ct.byref(NumElems), ct.byref(NumNodes), W_COST.ctypes.data_as(ct.c_void_p), 
     ct.byref(NumNodesElem), ElemType, TestCase, BConds,
     ct.byref(BeamLength1), ct.byref(BeamLength2),
     BeamStiffness.ctypes.data_as(ct.c_void_p), BeamMass.ctypes.data_as(ct.c_void_p),
     ExtForce.ctypes.data_as(ct.c_void_p), ExtMomnt.ctypes.data_as(ct.c_void_p),
     ct.byref(SectWidth), ct.byref(SectHeight),
     ct.byref(ThetaRoot), ct.byref(ThetaTip),
     ct.byref(TipMass), ct.byref(TipMassY), ct.byref(TipMassZ),
     ct.byref(Omega),
     ct.byref(FollowerForce), ct.byref(FollowerForceRig), ct.byref(PrintInfo),   
     ct.byref(OutInBframe), ct.byref(OutInaframe), ct.byref(ElemProj), ct.byref(MaxIterations),
     ct.byref(NumLoadSteps), ct.byref(Solution), ct.byref(DeltaCurved), ct.byref(MinDelta),
     ct.byref(NewmarkDamp) )


# check output
print 'NOPT_MAX %d' % (NumElems.value)
print 'Cost Weights:', W_COST 
print 'NumLoadssteps: ', NumLoadSteps.value
print 'FollowerForce: ', FollowerForce.value 






