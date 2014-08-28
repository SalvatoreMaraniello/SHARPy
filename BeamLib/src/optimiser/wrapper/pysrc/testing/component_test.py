'''
Created on 18 Aug 2014

@author: sm6110


Test Constant component

Remark:  

'''

import numpy as np
import ctypes as ct

import sys
sys.path.append('..')
sys.path.append('../..')
import shared

from xbcomponent import XBeamSolver

#------------------------------------------------------------------ Set Testcase
testcase = 'constant'

if testcase == 'constant':
    
    from lib.plot import sta_unif
    
    
    xbsolver=XBeamSolver()
    
    xbsolver.TestCase='test'
    ElemType="DISP"
    xbsolver.BConds='CF'
    
    xbsolver.NumElems = 10
    xbsolver.NumNodesElem = 3
    
    # Design
    xbsolver.beam_shape = 'constant'
    xbsolver.material ='alumA'
    xbsolver.cross_section_type ='isorect'    
    xbsolver.BeamLength1 = 5.0   
    
    # match NCB1 bending stiffness
    E = 7.e10
    A = 4.8e8/E
    I = 9.346e6/E
    print 'I = ', I
    xbsolver.cs_l2 = 0.2
    xbsolver.cs_l3 = np.power( 12.0*I/xbsolver.cs_l2 , 1.0/3.0)
    xbsolver.cs_l3 = 0.2
    print 'l2, l3 =', xbsolver.cs_l2, xbsolver.cs_l3

    # Loading
    xbsolver.ExtForce[2]=6e3 #600.e3

    # Options - NCB1 case:
    xbsolver.FollowerForce=True
    xbsolver.PrintInfo=False  
    xbsolver.FollowerForceRig= True
    xbsolver.OutInBframe = True
    xbsolver.OutInaframe = False    
                  
    xbsolver.MaxIterations=9999    
    xbsolver.NumLoadSteps=80
    xbsolver.Solution=112
    xbsolver.MinDelta= 1.e-5

    ## cost
    xbsolver.NNode = -1 # tip displacement

    xbsolver.execute() 
    
    print 'External Force: ', xbsolver.ExtForce
    print 'Cost Functions:'
    print 'Total Mass', xbsolver.TotalMass, xbsolver.ct_TotalMass.value
    print 'Node Abs Displacement', xbsolver.NodeDisp, xbsolver.ct_NodeDisp.value 
    print 'Node z displacement', xbsolver.ZDisp
    
    # Make a nice plot
    sta_unif(xbsolver.PosIni,xbsolver.PosDef)



