'''

Salvatore Maraniello 8 Aug 2014

OpenMDAO component for Beam solver

- Remarks:
   _Solution_List possible values:
       102/112: cbeam3 linear/nonlinear static.
       202/212: cbeam3 linear/nonlinear structural dynamic
       302/312: cbeam3 linear/nonlinear static + structural dynamic
       900/910:        linear/nonlinear rigid-body dynamic
       902/912: cbeam3 linear/nonlinear flexible-body dynamic
       922: cbeam3 nonlinear static + flexible-body dynamic
       952: cbeam3 linear flexible with nonlinear rigid-body dynamic

'''

import sys
import numpy as np
import ctypes as ct

sys.path.append('..')

from openmdao.main.api import Component, ComponentWithDerivatives
from openmdao.main.datatypes.api import Float, Array, Enum, Int, Str, Bool

import shared
import beamvar
import cost
        
# Load Dynamic Library
wrapso = ct.cdll.LoadLibrary(shared.wrapso_abspath)
        



class XBeamSolver(Component):
 
    #-------------------------------------------------------- OpenMDAO Interface
    
    '''Options:'''
    FollowerForce   = Bool( True, iotype='in', desc='Follower Force; default: True')  
    FollowerForceRig= Bool( True, iotype='in', desc='Follower force in the body-fixed frame; default: True')  
    PrintInfo       = Bool( True, iotype='in', desc='Print info on screen')  
    OutInBframe     = Bool( True, iotype='in', desc='print velocities in B-frame (if not, use a-frame). default: True')  
    OutInaframe     = Bool(False, iotype='in', desc='print velocities in a-frame (if not, Inertial frame). default: False') 
      
    ElemProj        = Enum( 0, (0,1,2), iotype='in', desc='Element info computed in (0) global frame, (1) fixed element frame, (2) moving element frame.')  
    _Solution_List = (102,112,202,212,302,312,900,910,902,912,922,952) # see header for possible values
    Solution        = Enum(112, _Solution_List, iotype='in', desc='Solution Process')  
       
    MaxIterations   = Int( 99, iotype='in', desc='Maximum number of iterations')          
    NumLoadSteps    = Int(  5, iotype='in', desc='Number of load increments')                      
                
    DeltaCurved     = Float( 1e-5, iotype='in', desc='Minimum angle for two unit vectors to be parallel')        
    MinDelta        = Float( 1e-8, iotype='in', desc='Convergence parameter for Newton-Raphson iterations')               
    NewmarkDamp     = Float( 1e-4, iotype='in', desc='Numerical damping in the Newmark integration scheme')       
   
   
    '''Shared Input'''
    NumElems     = Int(10, iotype='in', desc='Total number of elements')
    
    NumNodesElem = Enum(3, (2,3), iotype='in', desc='Number of nodes per element; (2) linear, (3) quadratic.')
    ElemType     = Enum('DISP',('DISP','STRN'), iotype='in')
    BConds       = Enum('CF'  ,('CF','FF')    , iotype='in', desc='Boundary conditions, (CF) clamped-free, (FF) free-free.') 
      
    TestCase = Str('test', iotype='in', desc='Test case name.')


    '''Design'''
    BeamLength1 = Float( 0.0, iotype='in', desc='')
    BeamLength2 = Float( 0.0, iotype='in', desc='')
    ThetaRoot   = Float( 0.0, iotype='in', desc='')
    ThetaTip    = Float( 0.0, iotype='in', desc='')
    TipMass     = Float( 0.0, iotype='in', desc='')
    TipMassY    = Float( 0.0, iotype='in', desc='')
    TipMassZ    = Float( 0.0, iotype='in', desc='')
    SectWidth   = Float( 0.0, iotype='in', desc='')
    SectHeight  = Float( 0.0, iotype='in', desc='')
    Omega       = Float( 0.0, iotype='in', desc='') 
 

    ''' Cost '''
    TotalMass =  Float( 0.0, iotype='out', desc='Total structural mass')
    NodeDisp  =  Float( 0.0, iotype='out', desc='Absolute value of NNode displacement')
    NNode     =    Int(  -1, iotype='in', desc='Node at which monitor the displacement. If negative, the last node is picked up')
 
    
    def __init__(self):
        
        # not to overwrite Component __init__
        super(XBeamSolver, self).__init__()
        
        # extract main routine
        ### method below fails inside a class but runs from main
        ###fwd_run=xbso.__opt_routine_MOD_opt_main 
        self.fwd_run = getattr(wrapso, "__opt_routine_MOD_opt_main")       
        
        ''' set ctypes '''
        # Applied external forces at the Tip
        self.ExtForce      = np.zeros(  (3), dtype=float, order='F')
        self.ExtMomnt      = np.zeros(  (3), dtype=float, order='F')    
        self.BeamStiffness = np.zeros((6,6), dtype=float, order='F')
        self.BeamMass      = np.zeros((6,6), dtype=float, order='F')
   
    

    def pre_solve(self):
        
        # ------------------------------------------------------- Determine Arrays Size
        self.NumNodes = beamvar.total_nodes(self.NumElems,self.NumNodesElem)
        self._NumNodes_copy = self.NumNodes  
        
        # set value for node displacement monitoring
        if self.NNode==-1:
            print 'setting NNode'
            self.NNode = self.NumNodes      



    def execute(self):
        
        # set size etc.
        self.pre_solve()
        
        # set ctypes for Fortran interface
        self.set_ctypes()
        
        # call xbeam solver
        self.fwd_run(ct.byref(self.ct_NumElems), ct.byref(self.ct_NumNodes),
                     ct.byref(self.ct_NumNodesElem), self.ct_ElemType, self.ct_TestCase, self.ct_BConds,
                     ct.byref(self.ct_BeamLength1), ct.byref(self.ct_BeamLength2),
                     self.BeamStiffness.ctypes.data_as(ct.c_void_p), self.BeamMass.ctypes.data_as(ct.c_void_p),
                     self.ExtForce.ctypes.data_as(ct.c_void_p), self.ExtMomnt.ctypes.data_as(ct.c_void_p),
                     ct.byref(self.ct_SectWidth), ct.byref(self.ct_SectHeight),
                     ct.byref(self.ct_ThetaRoot), ct.byref(self.ct_ThetaTip),
                     ct.byref(self.ct_TipMass), ct.byref(self.ct_TipMassY), ct.byref(self.ct_TipMassZ),
                     ct.byref(self.ct_Omega),
                     ct.byref(self.ct_FollowerForce), ct.byref(self.ct_FollowerForceRig), ct.byref(self.ct_PrintInfo),   
                     ct.byref(self.ct_OutInBframe), ct.byref(self.ct_OutInaframe), ct.byref(self.ct_ElemProj), ct.byref(self.ct_MaxIterations),
                     ct.byref(self.ct_NumLoadSteps), ct.byref(self.ct_Solution), ct.byref(self.ct_DeltaCurved), ct.byref(self.ct_MinDelta), ct.byref(self.ct_NewmarkDamp),
                     self.PosIni.ctypes.data_as(ct.c_void_p), self.PsiIni.ctypes.data_as(ct.c_void_p),                  
                     self.PosDef.ctypes.data_as(ct.c_void_p), self.PsiDef.ctypes.data_as(ct.c_void_p), self.InternalForces.ctypes.data_as(ct.c_void_p),
                     self.DensityVector.ctypes.data_as(ct.c_void_p), self.LengthVector.ctypes.data_as(ct.c_void_p)    )
        
        # output checks
        self.check()
        
        
        # compute cost
        cost.f_total_mass( self.DensityVector.ctypes.data_as(ct.c_void_p), self.LengthVector.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumElems), ct.byref(self.ct_TotalMass) )
        cost.f_node_disp( self.PosIni.ctypes.data_as(ct.c_void_p), self.PosDef.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumNodes), ct.byref(self.ct_NodeDisp), ct.byref(self.ct_NNode) )
        
        
        # update values
        self.update_variables()
        
              
        
    def set_ctypes(self):
        
        # --------------------------------------------------------------- Prepare input
        # arrays do not require any modification
        self.ct_NumElems     = ct.c_int( self.NumElems     )
        self.ct_NumNodesElem = ct.c_int( self.NumNodesElem )
        self.ct_NumNodes     = ct.c_int( self.NumNodes )
        
        self.ct_ElemType     = ct.c_char_p( self.ElemType.ljust(4) )     # the ljust is not required
        self.ct_TestCase     = ct.c_char_p( self.TestCase.ljust(4) )
        self.ct_BConds       = ct.c_char_p( self.BConds.ljust(2)   )
        
        self.ct_BeamLength1 = ct.c_double( self.BeamLength1 )
        self.ct_BeamLength2 = ct.c_double( self.BeamLength2 )
        self.ct_ThetaRoot   = ct.c_double( self.ThetaRoot   )
        self.ct_ThetaTip    = ct.c_double( self.ThetaTip    )
        self.ct_TipMass     = ct.c_double( self.TipMass     )
        self.ct_TipMassY    = ct.c_double( self.TipMassY    )
        self.ct_TipMassZ    = ct.c_double( self.TipMassZ    )
        self.ct_SectWidth   = ct.c_double( self.SectWidth   )
        self.ct_SectHeight  = ct.c_double( self.SectHeight  )
        self.ct_Omega       = ct.c_double( self.Omega       )
        
        # ------------------------------------------------------------- Prepare Options
        self.ct_FollowerForce    = ct.c_bool( self.FollowerForce    )
        self.ct_FollowerForceRig = ct.c_bool( self.FollowerForceRig )  
        self.ct_PrintInfo        = ct.c_bool( self.PrintInfo        )
        self.ct_OutInBframe      = ct.c_bool( self.OutInBframe      )
        self.ct_OutInaframe      = ct.c_bool( self.OutInaframe      )   
        
        self.ct_ElemProj      = ct.c_int( self.ElemProj      )      
        self.ct_MaxIterations = ct.c_int( self.MaxIterations )         
        self.ct_NumLoadSteps  = ct.c_int( self.NumLoadSteps  )   
                 
        self.ct_Solution      = ct.c_int( self.Solution      )  
        
        self.ct_DeltaCurved   = ct.c_double( self.DeltaCurved )       
        self.ct_MinDelta      = ct.c_double( self.MinDelta    )      
        self.ct_NewmarkDamp   = ct.c_double( self.NewmarkDamp )
        
        
        #-------------------------------------------------------- Prepare Output
        
        # General
        self.DensityVector = np.zeros((self.NumElems),dtype=float,order='F')
        self.LengthVector  = np.zeros((self.NumElems),dtype=float,order='F')
        
        # Problem dependent
        [self.PosIni, self.PosDef, self.PsiIni, self.PsiDef, self.InternalForces] = beamvar.fwd_static(self.NumNodes,self.NumElems)
        
        
        #---------------------------------------------------------- Prepare Cost
        self.ct_TotalMass = ct.c_double( self.TotalMass ) 
        self.ct_NodeDisp  = ct.c_double( self.NodeDisp )
        self.ct_NNode     = ct.c_int(self.NNode)

        
    def update_variables(self):

        # --------------------------------------------------------------- Prepare input
        # these values are not supposed to change.
        #self.NumElems = self.ct_NumElems.value    
        #self.NumNodesElem  = self.ct_NumNodesElem.value
        self.NumNodes = self.ct_NumNodes.value
        
        # update shared design
        self.BeamLength1 = self.ct_BeamLength1.value 
        self.BeamLength2 = self.ct_BeamLength2.value
        self.ThetaRoot = self.ct_ThetaRoot.value  
        self.ThetaTip  = self.ct_ThetaTip.value
        self.TipMass = self.ct_TipMass.value
        self.TipMassY = self.ct_TipMassY.value
        self.TipMassZ = self.ct_TipMassZ.value
        self.SectWidth = self.ct_SectWidth.value
        self.SectHeight = self.ct_SectHeight.value
        self.Omega = self.ct_Omega.value
        
        #---------------------------------------------------------- Update Cost
        self.TotalMass = self.ct_TotalMass.value
        self.NodeDisp = self.ct_NodeDisp.value
        self.NNode = self.ct_NNode.value
        
        
    def set_arrays_order(self):
        print 'lala'
        
    

    def check(self):
            
        # Arrays size:
        if self.NumNodes != self._NumNodes_copy:
            raise NameError('NumNodes has been changed during the Fortran code execution!') 




#------------------------------------------------------------------------------ 
        
    
if __name__=='__main__':
    
    xbsolver=XBeamSolver()
    
    xbsolver.TestCase='NCB1'
    ElemType="DISP"
    xbsolver.BConds='CF'
    
    xbsolver.NumElems = 10
    xbsolver.NumNodesElem = 3
    
    # Design
    xbsolver.BeamLength1 = 5.0
    
    # Options - NCB1 case:
    xbsolver.FollowerForce=False
    xbsolver.PrintInfo=False              
    xbsolver.MaxIterations=99    
    xbsolver.NumLoadSteps=10
    xbsolver.Solution=112
    xbsolver.MinDelta= 1.e-5

    # cost
    xbsolver.NNode = -1 # tip displacement

    # reminder:
    # pay attention here not to reallocate (order='F' to be kept)
    xbsolver.BeamMass[0,0] = 100.e0        # m [kg/m]
    xbsolver.BeamMass[1,1] = xbsolver.BeamMass[0,0]
    xbsolver.BeamMass[2,2] = xbsolver.BeamMass[0,0]
    xbsolver.BeamMass[3,3] = 10.e0         # J [kgm]
    xbsolver.BeamMass[4,4] = 10.e0
    xbsolver.BeamMass[5,5] = 10.e0     

    xbsolver.BeamStiffness[0,0] = 4.8e8   # EA [Nm]
    xbsolver.BeamStiffness[1,1] = 3.231e8 # GA
    xbsolver.BeamStiffness[2,2] = 3.231e8
    xbsolver.BeamStiffness[3,3] = 1.e6    # GJ
    xbsolver.BeamStiffness[4,4] = 9.346e6 # EI
    xbsolver.BeamStiffness[5,5] = 9.346e6

    xbsolver.ExtForce[2]=600.e3

    xbsolver.execute() 
    
    print 'Cost Functions:'
    print 'Total Mass', xbsolver.TotalMass, xbsolver.ct_TotalMass.value
    print 'Node Displacement', xbsolver.NodeDisp, xbsolver.ct_NodeDisp.value 
    
    

        
    
