'''
 Salvatore Maraniello 31 July 2014

 Python launcher for beam solver. The file allows to define the input in 
a python file (input_file.py) and call the beam solver.
 
 Remarks:
  - Output arrays are preallocated. The variables determining the size of the 
    output are:
        NumNodes: total number of nodes in the model
    To avoid to change complitely the structure of the fortran code, this is 
    allowed to change the values of such output. In such event, the wrapper 
    would crash. A copy of the original values is, therefore, stored, and a 
    check is performed once the execution is terminated.
  - The code does call an externally compiled dynamic library. As subroutines 
    names have a different name in the library then in the original code, to 
    print the library names, type (from terminal): 
    nm ./xbeamopt.so (ref.2)
  - To improve style, check out ref.4

 References for wrapper:
    ref.1: http://www.walkingrandomly.com/?p=85
    ref.2: http://stackoverflow.com/questions/5811949/call-functions-from-a-shared-fortran-library-in-python
    ref.3: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.ctypes.html
    ref.4: http://stackoverflow.com/questions/23641754/passing-string-to-fortran-dll-using-ctypes-and-python

 random quote: << Never try to out drink a Swede, unless you happen to be a 
                  Finn or at least a Russian. >>

'''

import os
import sys
import numpy as np
import ctypes as ct

import shared # this contains paths to modules etc.
import beamvar
import cost

# Change current directory


# Load Dynamic Library
xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)

# extract main routine
fwd_run=xb.__opt_routine_MOD_opt_main

#--------------------------------------------------- # Options - Default Values:
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


# ------------------------------------------------------------ Import all Input
# python './input_file.py'
from input_file import *


# ------------------------------------------------------- Determine Arrays Size
NumNodes = beamvar.total_nodes(NumElems,NumNodesElem)
NumNodes_copy = NumNodes

# --------------------------------------------------------------- Prepare input
# arrays do not require any modification
NumElems     = ct.c_int( NumElems     )
NumNodesElem = ct.c_int( NumNodesElem )
NumNodes     = ct.c_int( NumNodes )

ElemType     = ct.c_char_p( ElemType.ljust(4) )     # the ljust is not required
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


# ------------------------------------------------------------- Prepare Options

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


# -------------------------------------------------------------- Prepare Output

# General
DensityVector = np.zeros((NumElems.value),dtype=float,order='F')
LengthVector  = np.zeros((NumElems.value),dtype=float,order='F')

# Problem dependent
[PosIni, PosDef, PsiIni, PsiDef, InternalForces] = beamvar.fwd_static(NumNodes.value,NumElems.value)





# --------------------------------------------------------------------- Options
'''
print 'Values in input:'
print 'DeltaCurved %f' %(DeltaCurved.value)
print 'MinDelta    %f' %(MinDelta.value)
print 'NewmarkDamp %f' %(NewmarkDamp.value)
print 'ElemType: ',  ElemType.value
print 'TestCase: ',  TestCase.value
print 'BConds: '  ,  BConds.value
print 'Beamlength1: ', BeamLength1.value
print 'Beamlength2: ', BeamLength2.value
'''
# ----------------------------------------------------------------- Call Solver
#W_COST=np.empty((2),dtype=float,order='F')


fwd_run( ct.byref(NumElems), ct.byref(NumNodes),
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
     ct.byref(NumLoadSteps), ct.byref(Solution), ct.byref(DeltaCurved), ct.byref(MinDelta), ct.byref(NewmarkDamp),
     PosIni.ctypes.data_as(ct.c_void_p), PsiIni.ctypes.data_as(ct.c_void_p),                  
     PosDef.ctypes.data_as(ct.c_void_p), PsiDef.ctypes.data_as(ct.c_void_p), InternalForces.ctypes.data_as(ct.c_void_p),
     DensityVector.ctypes.data_as(ct.c_void_p), LengthVector.ctypes.data_as(ct.c_void_p)    )
    

# -----------------------------------------------------------------------
TotalMass = ct.c_double( 0.0 ) 
cost.f_total_mass( DensityVector.ctypes.data_as(ct.c_void_p), LengthVector.ctypes.data_as(ct.c_void_p), ct.byref(NumElems), ct.byref(TotalMass) )
NodeDisp = ct.c_double( 0.0 )

NNode = ct.c_int(NumNodes.value)
cost.f_node_disp( PosIni.ctypes.data_as(ct.c_void_p), PosDef.ctypes.data_as(ct.c_void_p), ct.byref(NumNodes), ct.byref(NodeDisp), ct.byref(NNode) )

# call without optional argument doesnt work
#cost.f_node_disp( PosIni.ctypes.data_as(ct.c_void_p), PosDef.ctypes.data_as(ct.c_void_p), ct.byref(NumNodes), ct.byref(NodeDisp))

print 'Cost Functions:'
print 'Total Mass', TotalMass.value
print 'Tip Displacement', NodeDisp.value

# ---------------------------------------------------------------------- Checks

# Arrays size:
if NumNodes.value != NumNodes_copy:
    raise NameError('NumNodes has been changed during the fortran code exewcution!') 

# print stuff
print 'NumNodes %d' %(NumNodes.value)
print 'NOPT_MAX %d' %(NumElems.value)
print 'NumLoadssteps: ', NumLoadSteps.value
print 'FollowerForce: ', FollowerForce.value 






