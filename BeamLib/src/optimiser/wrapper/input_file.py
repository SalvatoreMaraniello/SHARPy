#
# Beam Solver Input File
#
# ---------------------------------------------------------------------------------------------------------------
import numpy as np
# ---------------------------------------------------------------------------------------------------------------

# Options and Problem Setup
NumElems = 10
NumNodesElem = 3

ElemType="DISP"
TestCase='NCB1'
BConds='CF'

# Design
BeamLength1 = 5.0
BeamLength2 = 0.0
ThetaRoot   = 0.0
ThetaTip    = 0.0

# Applied external forces at the Tip
ExtForce= np.array([ 0.0, 0.0, 600.e3 ],dtype=float,order='F')
ExtMomnt= np.array([ 0.0, 0.0, 0.0    ],dtype=float,order='F')

TipMass    =0.0
TipMassY   =0.0
TipMassZ   =0.0
SectWidth  =0.0
SectHeight =0.0
Omega      =0.0

BeamStiffness=np.zeros((6,6),dtype=float,order='F')
BeamStiffness[0,0] = 4.8e8   # EA [Nm]
BeamStiffness[1,1] = 3.231e8 # GA
BeamStiffness[2,2] = 3.231e8
BeamStiffness[3,3] = 1.e6    # GJ
BeamStiffness[4,4] = 9.346e6 # EI
BeamStiffness[5,5] = 9.346e6

BeamMass=np.zeros((6,6),dtype=float,order='F')
BeamMass[0,0] = 100.e0        # m [kg/m]
BeamMass[1,1] = BeamMass[0,0]
BeamMass[2,2] = BeamMass[0,0]
BeamMass[3,3] = 10.e0         # J [kgm]
BeamMass[4,4] = 10.e0
BeamMass[5,5] = 10.e0


# Options - NCB1 case:
FollowerForce=False
PrintInfo=False              
MaxIterations=99    
NumLoadSteps=10
Solution=112
MinDelta= 1.e-5
# FollowerForceRig = True 
# OutInBframe      = True   
# OutInaframe      = False  
# ElemProj         = 0                  
# DeltaCurved      = 1e-5         
# MinDelta         = 1e-8            
# NewmarkDamp      = 1.e-4   




