'''
Created on 18 Aug 2014

@author: sm6110


Test Constant component against BEND case (cardona)

Remark: shear 

'''


import sys
import numpy as np
import ctypes as ct

sys.path.append('..')
from component import XBeamSolver


#------------------------------------------------------------------ Set Testcase
testcase = 'constant'

if testcase == 'constant':
    
    xbsolver=XBeamSolver()
    
    xbsolver.TestCase='BEND'
    ElemType="DISP"
    xbsolver.BConds='CF'
    
    xbsolver.NumElems = 10
    xbsolver.NumNodesElem = 3
    
    # Design
    xbsolver.BeamLength1 = 5.0
    xbsolver.beam_shape = 'constant'
    
    
    
    
    
    
    # Options - NCB1 case:
    xbsolver.FollowerForce=True
    xbsolver.PrintInfo=False              
    xbsolver.MaxIterations=99    
    xbsolver.NumLoadSteps=10
    xbsolver.Solution=112
    xbsolver.MinDelta= 1.e-5




    # cost
    xbsolver.NNode = -1 # tip displacement

    xbsolver.ExtForce[2]=600.e3

    xbsolver.execute() 
    
    print 'Cost Functions:'
    print 'Total Mass', xbsolver.TotalMass, xbsolver.ct_TotalMass.value
    print 'Node Displacement', xbsolver.NodeDisp, xbsolver.ct_NodeDisp.value 



