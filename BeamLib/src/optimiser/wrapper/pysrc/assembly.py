'''
Salvatore Maraniello 8 Aug 2014

OpenMDAO assembly for Beam solver

Rigid Body Dynamics and Structural Dynamics are treated as a single discipline.

The coupling is cared of inside the original Fortran code.

'''
import sys
sys.path.append('..')

import shared
from xbcomponent import XBeamSolver

from openmdao.main.api import Assembly
from openmdao.lib.drivers.slsqpdriver import SLSQPdriver
from openmdao.lib.drivers.doedriver import DOEdriver



class MonoliticOpt(Assembly):
    """ A basic Optimiser """

    def configure(self):
        
        # Create Optimiser instance:
        self.add('driver', SLSQPdriver())

        # Gradient Calculation:
        #self.driver.differentiator = FiniteDifference()  
        #self.driver.differentiator.form = 'forward'
        #self.driver.differentiator.default_stepsize = 1.0
        
        # SQLSQP Flags
        self.driver.iprint = 1
        self.driver.accuracy = 1e-6

        # Create component for structural solver
        self.add('xbeam', XBeamSolver())

        # Now you can access typing self.driver/self.ellwing

        # Add to driver's workflow
        # this tells the driver what to run. you could insert here sub-drivers
        # to perform sub-optimisations
        self.driver.workflow.add('xbeam')
        
        # Objective
        #self.driver.add_objective('xbeam.TotalMass')
        self.driver.add_objective('xbeam.ZDisp')

        # Constraint for study02a
        #self.driver.add_constraint('xbeam.NodeDisp > 1.76')
        self.driver.add_constraint('xbeam.TotalMass = 135.0')
        
        # Design Variables
        #self.driver.add_parameter('xbeam.BeamLength1', 4.5, 5.5)
        self.driver.add_parameter('xbeam.cs_l2', 0.9, 1.1)
        self.driver.add_parameter('xbeam.cs_l3', 0.9, 1.1)
        
 
 
class DOE(Assembly):
    """ A basic DOE """

    def configure(self):
        
        # Create Optimiser instance:
        self.add('driver', DOEdriver())
        self.driver.DOEgenerator = FullFactorial(10)
        
        self.driver.iprint = 0
        
        # Create component for structural solver
        self.add('xbeam', XBeamSolver())
        self.driver.workflow.add('xbeam')
        
        self.driver.case_outputs=['xbeam.TotalMass','xbeam.ZDisp','xbeam.NodeDisp']
        
        # Design Variables
        self.driver.add_parameter('xbeam.cs_l2', 0.8, 1.6)
        self.driver.add_parameter('xbeam.cs_l3', 0.8, 1.6)  
        
        self.driver.recorders = [ListCaseRecorder(),] 
 
 
       

if __name__=='__main__':
    
    import numpy as np

    opt_problem = MonoliticOpt()   
    
    opt_problem.xbeam.TestCase='test'
    ElemType="DISP"
    opt_problem.xbeam.BConds='CF'
    
    opt_problem.xbeam.NumElems = 10
    opt_problem.xbeam.NumNodesElem = 3
    
    # Design
    opt_problem.xbeam.beam_shape = 'constant'
    opt_problem.xbeam.material ='allum'
    opt_problem.xbeam.cross_section_type ='isorect'    
    opt_problem.xbeam.BeamLength1 = 5.0   
    
    # match NCB1 bending stiffness
    E = 7.e10
    A = 4.8e8/E
    I = 9.346e6/E
    print 'I = ', I
    opt_problem.xbeam.cs_l2 = 0.8
    #opt_problem.xbeam.cs_l3 = np.power( 12.0*I/opt_problem.xbeam.cs_l2 , 1.0/3.0)
    opt_problem.xbeam.cs_l3 = 0.8
    
    print 'l2, l3 =', opt_problem.xbeam.cs_l2, opt_problem.xbeam.cs_l3

    # Loading
    opt_problem.xbeam.ExtForce[2]=600.e3 /1e3

    # Options - NCB1 case:
    opt_problem.xbeam.FollowerForce=True
    opt_problem.xbeam.PrintInfo=False  
    opt_problem.xbeam.FollowerForceRig= True
    opt_problem.xbeam.OutInBframe = True
    opt_problem.xbeam.OutInaframe = False    
                  
    opt_problem.xbeam.MaxIterations=999    
    opt_problem.xbeam.NumLoadSteps=20
    opt_problem.xbeam.Solution=112
    opt_problem.xbeam.MinDelta= 1.e-5

    ## cost
    opt_problem.xbeam.NNode = -1 # tip displacement

    
    opt_problem.run()
    
    print 'Cost Functions:'
    print 'Total Mass', opt_problem.xbeam.TotalMass, opt_problem.xbeam.ct_TotalMass.value
    print 'Node Displacement', opt_problem.xbeam.NodeDisp, opt_problem.xbeam.ct_NodeDisp.value
    print 'l2, l3 =', opt_problem.xbeam.cs_l2, opt_problem.xbeam.cs_l3  
                   
