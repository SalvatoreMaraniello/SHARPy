
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



# Options and Problem Setup
NumElems = 10
NumNodes = 3

NumNodesElem = 3
ElemType='DISP'
TestCase='NCB1'
BConds='CF'

# Design
BeamLength1 = 5.0
BeamLength2 = 0.0
ThetaRoot   = 0.0
ThetaTip    = 0.0
ExtForce= np.array([ 0.0, 0.0, 600.e3 ],dtype=float,order='F')
ExtMomnt= np.array([ 0.0, 0.0, 0.0    ],dtype=float,order='F')

TipMass =0.0
TipMassY=0.0
TipMassZ=0.0
SectWidth=0.0
SectHeight=0.0
Omega=0.0

BeamStiffness=np.zeros((6,6),dtype=float,order='F')
BeamStiffness[0,0]= 4.8e8   # EA [Nm]
BeamStiffness[1,1]= 3.231e8 # GA
BeamStiffness[2,2]= 3.231e8
BeamStiffness[3,3]= 1.e6    # GJ
BeamStiffness[4,4]= 9.346e6 # EI
BeamStiffness[5,5]= 9.346e6

BeamMass=np.zeros((6,6),dtype=float,order='F')
BeamMass[0,0]=100.e0        # m [kg/m]
BeamMass[1,1]=BeamMass[0,0]
BeamMass[2,2]=BeamMass[0,0]
BeamMass[3,3]=10.e0         # J [kgm]
BeamMass[4,4]=10.e0
BeamMass[5,5]=10.e0


# Options - Default Values:
FollowerForce   = True
FollowerForceRig= True   
PrintInfo       = True    
OutInBframe     = True   
OutInaframe     = False   
ElemProj     = 0       
MaxIterations=99         
NumLoadSteps=5         
NumGauss=1                
Solution=111              
DeltaCurved=1e-5         
MinDelta=1d-8            
NewmarkDamp=1.e-4      


# NCB1 case:
FollowerForce=False
PrintInfo=False              
MaxIterations=99    
NumLoadSteps=10
Solution=112
MinDelta= 1.e-5





# ------------------------------------------------------------------------- prepare input:
# arrays do not require any modification
NumElems     = ct.c_int( NumElems     )
NumNodes     = ct.c_int( NumNodes     )
NumNodesElem = ct.c_int( NumNodesElem )

ElemType     = ct.c_char_p( ElemType )
TestCase     = ct.c_char_p( TestCase )
BConds       = ct.c_char_p( BConds   )

BeamLength1 = ct.c_float( BeamLength1 )
BeamLength2 = ct.c_float( BeamLength2 )
ThetaRoot   = ct.c_float( ThetaRoot   )
ThetaTip    = ct.c_float( ThetaTip    )
TipMass     = ct.c_float( TipMass     )
TipMassY    = ct.c_float( TipMassY    )
TipMassZ    = ct.c_float( TipMassZ    )
SectWidth   = ct.c_float( SectWidth   )
SectHeight  = ct.c_float( SectHeight  )
Omega       = ct.c_float( Omega       )


# ------------------------------------------------------------------------- Options

<---------------------- CONTINUE FROM HERE!!!
< Convert all types into ctype classes

FollowerForce   = True
FollowerForceRig= True   
PrintInfo       = True    
OutInBframe     = True   
OutInaframe     = False   
ElemProj     = 0       
MaxIterations=99         
NumLoadSteps=5         
NumGauss=1                
Solution=111              
DeltaCurved=1e-5         
MinDelta=1d-8            
NewmarkDamp=1.e-4   



# ---------------------------------------------------------------------------------------
###COST=np.empty((7),dtype=float,order='F')
W_COST=np.empty((2),dtype=float,order='F')


##### Launch Solver
#fm( ct.byref(NumElems), ct.byref(NumNodes), W_COST.ctypes.data_as(ct.c_void_p) )


# check output
print 'NOPT_MAX %d' % (NumElems.value)
print 'Cost Weights:', W_COST 



