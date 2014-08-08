'''
# Salvatore Maraniello 8 Aug 2014

OpenMDAO assembly for Beam solver

Rigid Body Dynamics and Structural Dynamics are treated as a single discipline.

The coupling is cared of inside the original Fortran code.

'''
import sys
sys.path.append('..')

import shared
from component import XBeamSolver

from openmdao.main.api import Assembly
from openmdao.lib.drivers.slsqpdriver import SLSQPdriver
from component import XBeamSolver


class MonoliticModel(Assembly):
    """ A basic Optimiser """

    def configure(self):
        
        # Create Optimiser instance:
        self.add('driver', SLSQPdriver())
        #self.add('driver', CONMINdriver())
        #self.add('driver', NEWSUMTdriver())
        
        # Gradient Calculation:
        #self.driver.differentiator = FiniteDifference()  
        #self.driver.differentiator.form = 'forward'
        #self.driver.differentiator.default_stepsize = 1.0
        
        # SQLSQP Flags
        self.driver.iprint = 1
        self.driver.accuracy = 1e-4

        # Create component for structural solver
        self.add('xbeam', XBeamSolver())

        # Now you can access typing self.driver/self.ellwing

        # Add to driver's workflow
        # this tells the driver what to run. you could insert here sub-drivers
        # to perform sub-optimisations
        self.driver.workflow.add('xbeam')
        
        # Objective
        self.driver.add_objective('xbeam.TotalMass')
        
        # Constraint for study02a
        self.driver.add_constraint('xbeam.NodeDisp > 1.76')
        
        # Design Variables
        self.driver.add_parameter('xbeam.BeamLength1', 4.5, 5.5)

       
       

if __name__=='__main__':
    

    opt_problem = MonoliticModel()   
    
    opt_problem.xbeam.TestCase='NCB1'
    opt_problem.xbeam.BConds='CF'
    
    opt_problem.xbeam.NumElems = 10
    opt_problem.xbeam.NumNodesElem = 3
    
    # Design
    opt_problem.xbeam.BeamLength1 = 5.0
    
    # Options - NCB1 case:
    opt_problem.xbeam.FollowerForce=False
    opt_problem.xbeam.PrintInfo=False              
    opt_problem.xbeam.MaxIterations=99    
    opt_problem.xbeam.NumLoadSteps=10
    opt_problem.xbeam.Solution=112
    opt_problem.xbeam.MinDelta= 1.e-5

    # cost
    opt_problem.xbeam.NNode = -1 # tip displacement

    # reminder:
    opt_problem.xbeam.BeamMass[0,0] = 100.e0        # m [kg/m]
    opt_problem.xbeam.BeamMass[1,1] = opt_problem.xbeam.BeamMass[0,0]
    opt_problem.xbeam.BeamMass[2,2] = opt_problem.xbeam.BeamMass[0,0]
    opt_problem.xbeam.BeamMass[3,3] = 10.e0         # J [kgm]
    opt_problem.xbeam.BeamMass[4,4] = 10.e0
    opt_problem.xbeam.BeamMass[5,5] = 10.e0     

    opt_problem.xbeam.BeamStiffness[0,0] = 4.8e8   # EA [Nm]
    opt_problem.xbeam.BeamStiffness[1,1] = 3.231e8 # GA
    opt_problem.xbeam.BeamStiffness[2,2] = 3.231e8
    opt_problem.xbeam.BeamStiffness[3,3] = 1.e6    # GJ
    opt_problem.xbeam.BeamStiffness[4,4] = 9.346e6 # EI
    opt_problem.xbeam.BeamStiffness[5,5] = 9.346e6

    opt_problem.xbeam.ExtForce[2]=600.e3

    
    opt_problem.run()
    
    print 'Cost Functions:'
    print 'Total Mass', opt_problem.xbeam.TotalMass, opt_problem.xbeam.ct_TotalMass.value
    print 'Node Displacement', opt_problem.xbeam.NodeDisp, opt_problem.xbeam.ct_NodeDisp.value     
                   
