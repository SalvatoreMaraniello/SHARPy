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

import sys, os
import numpy as np
import ctypes as ct

sys.path.append('..')

from openmdao.main.api import Component, ComponentWithDerivatives
from openmdao.main.datatypes.api import Float, Array, Enum, Int, Str, Bool

import shared
import design.beamelem, design.beamnodes, design.beamgeom
import input.load, input.forcedvel
import beamvar
import cost
import lib.save

        
# Load Dynamic Library
wrapso = ct.cdll.LoadLibrary(shared.wrapso_abspath)
        



class XBeamSolver(Component):
 
    #-------------------------------------------------------- OpenMDAO Interface

    '''Shared Input'''
    NumElems     = Int(10, iotype='in', desc='Total number of elements')
    NumNodesElem = Enum(3, (2,3), iotype='in', desc='Number of nodes per element; (2) linear, (3) quadratic.')
    ElemType     = Enum('DISP',('DISP','STRN'), iotype='in')
    
    # BCs available:
    # - CF, CC: static/dynamics/rigid-body
    # - CH, HC, HH: static only
    # - HF: non-linear coupled rigid body + structural
    BConds       = Enum('CF'  ,('CF','CC','CS','SC','SS','SF'), iotype='in', 
                        desc='Boundary conditions: C (clamped end), F (free end), S (spherical joint).')   
    TestCase = Str('test', iotype='in', desc='Test case name.')

    '''Options:'''
    FollowerForce   = Bool( True, iotype='in', desc='Follower Force; default: True')  
    FollowerForceRig= Bool( True, iotype='in', desc='Follower force in the body-fixed frame; default: True')  
    PrintInfo       = Bool( True, iotype='in', desc='Print info on screen')  
    OutInBframe     = Bool( True, iotype='in', desc='print velocities in B-frame (if not, use a-frame). default: True')  
    OutInaframe     = Bool(False, iotype='in', desc='print velocities in a-frame (if not, Inertial frame). default: False') 
      
    ElemProj        = Enum( 0, (0,1,2), iotype='in', desc='Element info computed in (0) global frame, (1) fixed element frame, (2) moving element frame.')  
    _Solution_List = (102,112,142,202,212,302,312,900,910,902,912,922,952) # see header for possible values
    Solution        = Enum(112, _Solution_List, iotype='in', desc='Solution Process')  
       
    MaxIterations   = Int( 99, iotype='in', desc='Maximum number of iterations')          
    NumLoadSteps    = Int(  5, iotype='in', desc='Number of load increments')                      
                
    DeltaCurved     = Float( 1e-5, iotype='in', desc='Minimum angle for two unit vectors to be parallel')        
    MinDelta        = Float( 1e-8, iotype='in', desc='Convergence parameter for Newton-Raphson iterations')               
    NewmarkDamp     = Float( 1e-4, iotype='in', desc='Numerical damping in the Newmark integration scheme')
   
    '''Design'''
    # general
    beam_geometry = Enum('straight',('straight'),iotype='in', desc='Defines beam geometry')
    BeamAxis = Array(np.array([1,0,0]),iotype='in',desc='Axis along which the straight beam develops')    
    BeamLength1 = Float( 0.0, iotype='in', desc='beam length')
    TipMass     = Float( 0.0, iotype='in', desc='')
    TipMassY    = Float( 0.0, iotype='in', desc='')
    TipMassZ    = Float( 0.0, iotype='in', desc='')

    # beam span reconstruction
    beam_shape = Enum('constant',('constant','spline','constant_hardcoded'),iotype='in',desc='Define how elements properties vary along the span')

 
    # isotropic material cross-sectional models(see design.isosec)
    material = Enum( 'alumA',('alumA','titnA','test'), iotype='in', desc='Material for isotropic cross section (see design.isosec)') 
    cross_section_type = Enum(('isorect','isoellip','isocirc',
                               'isohollowrect','isohollowellip','isohollowcirc',
                               'isorect_fact_torsion'),
                               iotype='in')
    cs_factor = Float(1e3, iotype='in',desc='Factor to be applied in isorec_fact_torsion. For spline method, this is constant along the span.')

    # parameters for constant beam span reconstruction 
    cs_l2   = Float(iotype='in',desc='size along the x2 axis - (see design.isosec)')
    cs_l3   = Float(iotype='in',desc='size along the x3 axis - (see design.isosec)')
    cs_t2   = Float(iotype='in',desc='thickness along the x2 axis - (see design.isosec)')
    cs_t3   = Float(iotype='in',desc='thickness along the x2 axis - (see design.isosec)')
    cs_r    = Float(iotype='in',desc='radius for circular cross-sections - (see design.isosec)')
    cs_t    = Float(iotype='in',desc='thickness for circular cross-sections - (see design.isosec)')
    cs_pretwist = Float(0.0,iotype='in',desc='beam pretwist')
    
    # parameters for spline reconstruction
    ControlElems = Array(iotype='in',dtype=np.int,desc="Control sections for beam shape deformation")
    SplineOrder = Int(3,iotype='in',desc='spline order (as per scipy.interpolate.interpolate.spline)')
    CS_l2   = Array(iotype='in',dtype=np.float,desc='size along the x2 axis at control sections - (see design.isosec)')
    CS_l3   = Array(iotype='in',dtype=np.float,desc='size along the x3 axis at control sections  - (see design.isosec)')
    CS_t2   = Array(iotype='in',dtype=np.float,desc='thickness along the x2 axis at control sections  - (see design.isosec)')
    CS_t3   = Array(iotype='in',dtype=np.float,desc='thickness along the x2 axis at control sections  - (see design.isosec)')
    CS_r    = Array(iotype='in',dtype=np.float,desc='radius for circular cross-sections at control sections  - (see design.isosec)')
    CS_t    = Array(iotype='in',dtype=np.float,desc='thickness for circular cross-sections at control sections  - (see design.isosec)')
    
    ControlNodes = Array(iotype='in',dtype=np.int,desc="Control nodes for beam shape deformation") 
    CS_pretwist = Array(iotype='in',dtype=np.float,desc='beam pre-twist at nodes')

   
    # Static Load Distribution
    AppStaLoadType = Enum('nodal', ('uniform','nodal','userdef'), 
                          desc='Static prescribed external load uniform, concentrated to a node or user defined. see load')
    NodeAppStaForce = Int( -1, iotype='in', desc='Node at which a static force is applied. -1 = tip')    
    # Dynamic Load distribution
    AppDynLoadType = Enum('nodal', ('uniform','nodal','userdef'), 
                          desc='Dynamic prescribed external load uniform or concentrated to a node or user defined. see load')
    NodeAppDynForce = Int( -1, iotype='in', desc='Node at which a dynamic force is applied. -1 = tip')    
    # Dynamic Load variation
    AppDynLoadVariation = Enum('constant', 
                               ('constant','impulse','sin','cos','ramp','omcos','rampsin','rect'), 
                               desc='Method for time variation of prescribed dynamic load')
    NumSteps = 0 # steps for dynamic simulation   
    TimeRamp = Float(0.0, iotype='in', desc ='For ramp AppDynLoadVariation or PrVelType, time required to reach full loading')
    Omega    = Float( 0.0, iotype='in', desc='Frequency for sin, cos, rampsin AppDynLoadVariation or PrVelType')   
    TimeRectStart = Float(0.0, iotype='in', desc ='For rect AppDynLoadVariation or PrVelType, time when signal starts')
    TimeRectEnd   = Float(0.0, iotype='in', desc ='For rect AppDynLoadVariation or PrVelType, time when signal ends')
 
    # Prescribed velocities at the support
    PrVelType = Enum('zero', ('zero','sin','ramp','rampsin','userdef'), desc='Prescribed velocity at beam support')
    
 
    ''' Control '''
    # Static
    ExtForce = Array(default_value=np.zeros((3),order='F'), # order='F' is not necessary
                     iotype='in', dtype=float, shape=(3,) ,
                     desc='Nodal applied force - single node or distributed according to self.AppStaLoadType')
    ExtMomnt = Array(default_value=np.zeros((3),order='F'), # order='F' is not necessary
                     iotype='in', dtype=float, shape=(3,) ,
                     desc='Nodal applied moments - single node or distributed according to self.AppStaLoadType')   
    # Dynamic
    ExtForceDyn = Array(default_value=np.zeros((3),order='F'), # order='F' is not necessary
                     iotype='in', dtype=float, shape=(3,) ,
                     desc='Nodal applied force - single node or distributed according to self.AppSDynLoadType (load.set_***_ForceDynAmp)')
    ExtMomntDyn = Array(default_value=np.zeros((3),order='F'), # order='F' is not necessary
                     iotype='in', dtype=float, shape=(3,) ,
                     desc='Nodal applied moments - single node or distributed according to self.AppDynLoadType (load.set_***_ForceDynAmp)')   
 
 
    ''' Optimisation '''
    # this allows to define costumised optimisation functions. The output value 
    # is stored in fval (scalar). The function is ffun and the arguments are 
    # in fargs. 
    # In some situations, the pre_solve method may have to be run before being
    # able to define the arguments.
    fval = Float( 0.0, iotype='out', desc='generic cost function')
    
    
    # Specific Cost/Constrains [superceeded]
    TotalMass =  Float( 0.0, iotype='out', desc='Total structural mass')
    NodeDisp  =  Float(  iotype='out', desc='Absolute value of NNode displacement')
    ZDisp     =  Float(  iotype='out', desc='z displacement of NNode')
    YDisp     =  Float(  iotype='out', desc='y displacement of NNode')
    XDisp     =  Float(  iotype='out', desc='x displacement of NNode')
    NNode     =    Int(  -1, iotype='in', desc='Node at which monitor the displacement. If negative, the last node is picked up')    
     
 
    def __init__(self):
        
        # not to overwrite Component __init__
        super(XBeamSolver, self).__init__()
        
        # extract main routine
        ### method below fails inside a class but runs from main
        ###fwd_run=xbso.__opt_routine_MOD_opt_main 
        self.fwd_run = getattr(wrapso, "__opt_routine_MOD_opt_main")       
        
        
        #--------------------------------------------------------- Input Related  
        self.BeamSpanStiffness = np.zeros( (self.NumElems,6,6), dtype=float, order='F' )
        self.BeamSpanMass      = np.zeros( (self.NumElems,6,6), dtype=float, order='F' )  
        
        # for constant_hardcoded method only
        self.Mroot = np.zeros((6,6),dtype=float,order='F')
        self.Kroot = np.zeros((6,6),dtype=float,order='F')


        # ------------------------------------------------------- External Loads
        ### commented below with #~# have been moved to head
        
        # Static        
        #~#self.ExtForce      = np.zeros(  (3), dtype=float, order='F')
        #~#self.ExtMomnt      = np.zeros(  (3), dtype=float, order='F')  
        # Dynamic Simulation
        self.Time        = np.zeros((1),dtype=float,order='F')     # simulation time
        #~#self.ExtForceDyn = np.zeros((3), dtype='float', order='F') # input for load.set_***_ForceDynAmp
        #~#self.ExtMomntDyn = np.zeros((3), dtype='float', order='F') # input for load.set_***_ForceDynAmp
        
        #------------------------------------------------- Prescribed Velocities
        self.PrTrVelAmpl  = np.zeros(  (3), dtype=float, order='F')
        self.PrRotVelAmpl = np.zeros(  (3), dtype=float, order='F')
        
        #------------------------------------------------------------- Utilities
        # save directory
        self._savedir = '.'
    
        # multiple runs
        self._counter = 0
        self._use_previous_solution = False
        
        # append code version
        self._version = shared.CodeVersion()    
        
        
        #------------------------------------------------------------------ Cost
        # dummy arguments for generic cost function
        self.ffun = cost.f_xyz_disp
        self.fargs=['PosIni','PosDef']
        self.fname = self.ffun.func_name


    def pre_solve(self):
        
        # ------------------------------------------------ Determine Arrays Size
        self.NumNodes = beamvar.total_nodes(self.NumElems,self.NumNodesElem)
        self._NumNodes_copy = self.NumNodes  


        #------------------------------------------- Take care of dynamics input
        # Set Arrays Order
        ### not required for 1D array, here for memo!
        self.Time = np.array(self.Time, dtype=float, order='F')
        self.NumSteps = len(self.Time)-1
        
        # Prescribed Velocities at support
        # only one method available here

        #----------------------------------------- Set External Prescribed Loads
        self.set_loads()
        
        
        #--------------------------------------------- set Prescribed Velocities
        self.set_prescr_vel()


        # output 
        # Current nodal position vector
        self.PosDotDef = np.zeros( (self.NumNodes,3), dtype='float', order='F' )
        # Current element orientation vectors
        self.PsiDotDef = np.zeros( (self.NumElems,3,3), dtype='float', order='F' ) # the 2nd index is the max nodes per element
        # Position vector/rotation history at beam tip.
        self.PosPsiTime= np.zeros( (self.NumSteps+1,6), dtype='float', order='F' )
        # History of Velocities
        self.VelocTime = np.zeros( (self.NumSteps+1,self.NumNodes), dtype='float', order='F' )
        # Position of all nodes wrt to global frame a for each time step
        self.DynOut    = np.zeros( ((self.NumSteps+1)*self.NumNodes,3), dtype='float', order='F' )
        

        if self.Solution == 142:
            self.ForcedVel = np.zeros( (1,6), dtype='float', order='F' )


        #-------------------------------- Set Rigid-Body + Dynamics input/output

        self.RefVel     =np.zeros( (self.NumSteps+1,6), dtype='float', order='F') 
        self.RefVelDot  =np.zeros( (self.NumSteps+1,6), dtype='float', order='F')         
        ###self.RefVel = self.ForcedVel.copy()
        ###self.RefVelDot=self.ForcedVelDot.copy()
        self.Quat =np.array([1,0,0,0],dtype=float,order='F')

        #-------------------------------------------------------- Construct Beam   
        if self.beam_shape=='constant':
            # set elements
            argslist = self.set_cross_section_arguments()
            self.BeamSpanMass, self.BeamSpanStiffness = design.beamelem.constant(self.NumElems,self.cross_section_type,argslist)
            # set nodes
            self.PhiNodes = design.beamnodes.constant(self.NumNodes, self.cs_pretwist)
        elif self.beam_shape == 'spline':
            # set elements
            ArgsList = []
            for cc in len(self.ControlElems):
                argslist = self.set_cross_section_arguments(ElementNumber=cc)       
                ArgsList.append(list(argslist))
            self.BeamSpanMass, self.BeamSpanStiffness = design.beamelem.spline(self.NumElems,self.cross_section_type,ArgsList,self.ControlElems,self.SplineOrder)
            # set nodes
            self.PhiNodes = design.beamnodes.spline(self.NumNodes, self.CS_pretwist, self.ControlNodes, self.SplineOrder)      
        elif self.beam_shape == 'constant_hardcoded':
            self.BeamSpanMass, self.BeamSpanStiffness = design.beamelem.constant_hardcoded(self.NumElems,self.Mroot,self.Kroot)
            self.PhiNodes = design.beamnodes.constant(self.NumNodes, self.cs_pretwist)

        
        #------------------------------------- Define Cost/Constrains parameters
        # set value for node displacement monitoring
        if self.NNode==-1:
            print 'setting NNode'
            self.NNode = self.NumNodes      


        #-------------------------------------------------------- Prepare Output
        # Design Related
        self.DensityVector = np.zeros((self.NumElems),dtype=float,order='F')
        self.LengthVector  = np.zeros((self.NumElems),dtype=float,order='F')
        

        ### sm 12 sep 2014:
        # Use old solution for initial guess in the optimisation. This if 
        # statement only prevents self.{Psi/Pos}Def not to be overwritten as
        # the self.{Pos/Psi}Ini are allocated here but always overwritten inside 
        # the Fortran code.
        if self._counter < 1 or self._use_previous_solution is False:
            # self.PosDef, self.PsiDef are set to zero
            # all others are empty arrays
            [self.PosIni, self.PosDef, self.PsiIni, self.PsiDef, self.InternalForces] = beamvar.fwd_static(self.NumNodes,self.NumElems) 
            # assign initial coordinates
            self.PosIni = design.beamgeom.straight(self.PosIni,self.BeamLength1,self.BeamAxis)


        # prepare saving directory:
        self._savedir=os.path.abspath(self._savedir)
        if os.path.isdir(self._savedir) is False:
            os.makedirs(self._savedir)
        

    def execute(self):
        
        # set size etc.
        self.pre_solve()

        # set ctypes for Fortran interface
        self.set_ctypes()
        
        #--------------------------------------------------------------- Testing
        #print self.Time
        #print self.ForceTime
        
        #self.Time = np.array(self.Time, dtype=float, order='F')
        #self.ForceTime = np.array(self.ForceTime, dtype=float, order='F')
        #self.ForceDynAmp = np.array(self.ForceDynAmp, dtype=float, order='F')
        #self.Time = np.zeros((5), dtype=float, order='F')
        
        #------------------------------------------------------------------------------ 
        
        ###print 'Starting Solution:'
        ###print 'Initial Tip Displacements:'
        ###if self._counter>1:
        ###    self.PosIni=self.PosDef.copy(order='F')
        ###print self.PosIni[-1,:]
        
        # save before run (in case of crash)
        counter_string =  '%.3d' % (self._counter)
        savename = self._savedir + '/' + self.TestCase + '_' + counter_string + '.h5'
        lib.save.h5comp( self, savename)
        
        # call xbeam solver

        self.fwd_run(ct.byref(self.ct_NumElems), ct.byref(self.ct_NumNodes),
                     ct.byref(self.ct_NumNodesElem), self.ct_ElemType, self.ct_TestCase, self.ct_BConds,
                     self.BeamSpanStiffness.ctypes.data_as(ct.c_void_p),        # Properties along the span
                     self.BeamSpanMass.ctypes.data_as(ct.c_void_p),                       
                     self.PhiNodes.ctypes.data_as(ct.c_void_p),                 # PreTwist angle along the span
                     self.ForceStatic.ctypes.data_as(ct.c_void_p), 
                     ct.byref(self.ct_TipMass), ct.byref(self.ct_TipMassY), ct.byref(self.ct_TipMassZ),
                     ct.byref(self.ct_FollowerForce), ct.byref(self.ct_FollowerForceRig), ct.byref(self.ct_PrintInfo),    # options
                     ct.byref(self.ct_OutInBframe), ct.byref(self.ct_OutInaframe), ct.byref(self.ct_ElemProj), ct.byref(self.ct_MaxIterations),
                     ct.byref(self.ct_NumLoadSteps), ct.byref(self.ct_Solution), ct.byref(self.ct_DeltaCurved), ct.byref(self.ct_MinDelta), ct.byref(self.ct_NewmarkDamp),
                     self.PosIni.ctypes.data_as(ct.c_void_p), self.PsiIni.ctypes.data_as(ct.c_void_p),       # output static           
                     self.PosDef.ctypes.data_as(ct.c_void_p), self.PsiDef.ctypes.data_as(ct.c_void_p), self.InternalForces.ctypes.data_as(ct.c_void_p),
                     self.DensityVector.ctypes.data_as(ct.c_void_p), self.LengthVector.ctypes.data_as(ct.c_void_p) , # output design - part 1
                     ct.byref(self.ct_NumSteps), self.Time.ctypes.data_as(ct.c_void_p),                            # dynamic input 
                     self.ForceTime.ctypes.data_as(ct.c_void_p), self.ForceDynAmp.ctypes.data_as(ct.c_void_p),     # input_dynforce
                     self.ForcedVel.ctypes.data_as(ct.c_void_p), self.ForcedVelDot.ctypes.data_as(ct.c_void_p),    # input_forcedvel
                     self.PosDotDef.ctypes.data_as(ct.c_void_p), self.PsiDotDef.ctypes.data_as(ct.c_void_p),
                     self.PosPsiTime.ctypes.data_as(ct.c_void_p), self.VelocTime.ctypes.data_as(ct.c_void_p),
                     self.DynOut.ctypes.data_as(ct.c_void_p)  ,                                              # output from sol 202, 212, 302, 312, 322
                     self.RefVel.ctypes.data_as(ct.c_void_p), self.RefVelDot.ctypes.data_as(ct.c_void_p), 
                     self.Quat.ctypes.data_as(ct.c_void_p)   )                # output rigid-body + dynamics
        
        print 'Out of Fortran!!!'
        
        # output checks
        self.check()
        
        # compute cost from Fortran functions [superceeded]
        cost.f_total_mass( self.DensityVector.ctypes.data_as(ct.c_void_p), self.LengthVector.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumElems), ct.byref(self.ct_TotalMass) )
        cost.f_node_disp( self.PosIni.ctypes.data_as(ct.c_void_p), self.PosDef.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumNodes), ct.byref(self.ct_NodeDisp), ct.byref(self.ct_NNode) )
        
        # update values
        self.update_variables()
        
        # compute cost from Python functions [superceeded]
        self.ZDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='z',NNode=self.NNode-1) # python indices starts from 0
        self.YDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='y',NNode=self.NNode-1)
        self.XDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='x',NNode=self.NNode-1)
        
        
        # compute general cost
        
        #self.fval = self.ffun(*self.fargs)
        tplargs = cost.return_args(self, self.fargs)
        self.fval = self.ffun(*tplargs)
        
        # sm: savename etc moved before solver started
        ### save
        #counter_string =  '%.3d' % (self._counter)
        #savename = self._savedir + '/' + self.TestCase + '_' + counter_string + '.h5'
        lib.save.h5comp( self, savename)
        
        # counter: for multiple executions
        self._counter = self._counter+1
        
                          
    def set_ctypes(self):
        
        # --------------------------------------------------------------- Prepare input
        # arrays do not require any modification
        self.ct_NumElems     = ct.c_int( self.NumElems     )
        self.ct_NumNodesElem = ct.c_int( self.NumNodesElem )
        self.ct_NumNodes     = ct.c_int( self.NumNodes )
        self.ct_NumSteps     = ct.c_int( self.NumSteps )
        
        self.ct_ElemType     = ct.c_char_p( self.ElemType.ljust(4) )     # the ljust is not required
        self.ct_TestCase     = ct.c_char_p( self.TestCase.ljust(4) )
        self.ct_BConds       = ct.c_char_p( self.BConds.ljust(2)   )
        
        self.ct_BeamLength1 = ct.c_double( self.BeamLength1 )
        self.ct_TipMass     = ct.c_double( self.TipMass     )
        self.ct_TipMassY    = ct.c_double( self.TipMassY    )
        self.ct_TipMassZ    = ct.c_double( self.TipMassZ    )
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
        self.TipMass = self.ct_TipMass.value
        self.TipMassY = self.ct_TipMassY.value
        self.TipMassZ = self.ct_TipMassZ.value
        self.Omega = self.ct_Omega.value
        
        #---------------------------------------------------------- Update Cost
        self.TotalMass = self.ct_TotalMass.value
        self.NodeDisp = self.ct_NodeDisp.value
        self.NNode = self.ct_NNode.value
                   

    def check(self):    
        # Arrays size:
        if self.NumNodes != self._NumNodes_copy:
            raise NameError('NumNodes has been changed during the Fortran code execution!') 


    def set_cross_section_arguments(self,**kwargs):
        ''' 
        Defines the list of input for the cross-section definition 
        
        ElementNumber: optional argument in case of spline shape beam. The
        parameters points to the elements of the self.CS_* arrays.
        
        Remark: for modified cross-sections (e.g. isosec.rect_fact_torsion), the
        stiffness factor is kept constant along the span also for spline 
        reconstruction
        
        '''
        
        cc = kwargs.get('ElementNumber', None)
        
        if cc is None:
            # get number of numerical input:
            if self.cross_section_type == 'isocirc':
                argslist=[self.cs_r, self.material]
            elif self.cross_section_type == 'isohollowcirc':
                argslist=[self.cs_r, self.cs_t, self.material]       
            elif self.cross_section_type == 'isorect' or self.cross_section_type == 'isoellip':              
                argslist=[self.cs_l2, self.cs_l3, self.material] 
            elif  self.cross_section_type == 'isorect_fact_torsion':
                argslist=[self.cs_l2, self.cs_l3, self.material, self.cs_factor ]                   
            elif self.cross_section_type == 'isohollowrect' or self.cross_section_type == 'isohollowellip':
                argslist=[self.cs_l2, self.cs_l3, self.cs_t2,self.cs_t3, self.material] 
            else:
                raise NameError('Cross Section Type "%s" not found!' %self.cross_section_type )           
        else:
            # get number of numerical input:
            if self.cross_section_type == 'isocirc':
                argslist=[self.CS_r[cc], self.material ]
            elif self.cross_section_type == 'isohollowcirc':
                argslist=[self.CS_r[cc], self.CS_t[cc], self.material ]       
            elif self.cross_section_type == 'isorect' or self.cross_section_type == 'isoellip':
                argslist=[self.CS_l2[cc], self.CS_l3[cc], self.material ] 
            elif  self.cross_section_type == 'isorect_fact_torsion':
                argslist=[self.CS_l2[cc], self.CS_l3[cc], self.material, self.cs_factor ] 
            elif self.cross_section_type == 'isohollowrect' or self.cross_section_type == 'isohollowellip':
                argslist=[self.CS_l2[cc], self.CS_l3[cc], self.CS_t2[cc],self.CS_t3[cc], self.material] 
            else:
                raise NameError('Cross Section Type "%s" not found!' %self.cross_section_type )              
        
        return argslist





#------------------------------------------------------------------------------ 
        

    def set_loads(self):
        
        #------------------------------------------------- Set Static Load input
        
        # Span Distribution of Prescribed Load           
        if self.AppStaLoadType == 'uniform':
            self.ForceStatic = input.load.set_unif_ForceDistr(
                                self.NumNodes, 
                                self.ExtForce, self.ExtMomnt)    
        elif self.AppStaLoadType == 'nodal':
            self.ForceStatic = input.load.set_nodal_ForceDistr(
                                self.NumNodes, self.NodeAppStaForce,
                                self.ExtForce, self.ExtMomnt)
        elif self.AppStaLoadType == 'userdef':
            # make sure the array order is 'F'
            self.ForceStatic = np.array(self.ForceStatic, dtype=float, order='F')
        else:
            raise NameError('AppStaLoadType not valid!')
            
        
        #--------------------------------------------- Set Dynamics input/output
        # these are always allocated (as they are passed in input to fortran
        # solver. 
        
        # Span Distribution of Prescribed Load           
        if self.AppDynLoadType == 'uniform':
            self.ForceDynAmp = input.load.set_unif_ForceDistr(
                                   self.NumNodes, self.ExtForceDyn, self.ExtMomntDyn)    
        elif self.AppDynLoadType == 'nodal':
            self.ForceDynAmp = input.load.set_nodal_ForceDistr(
                                self.NumNodes, self.NodeAppDynForce,
                                self.ExtForceDyn, self.ExtMomntDyn)
        elif self.AppDynLoadType == 'userdef':
            # make sure the array order is 'F'
            self.ForceDynAmp = np.array(self.ForceDynAmp, dtype=float, order='F')
        else:
            raise NameError('AppDynLoadType not valid!')
            
        # Time variation of Prescribed Load
        if self.AppDynLoadVariation == 'constant':
            self.ForceTime = input.load.set_const_ForceTime(self.Time)
        elif self.AppDynLoadVariation == 'impulse':
            self.ForceTime = input.load.set_impulse_ForceTime(self.Time)
        elif self.AppDynLoadVariation == 'ramp':
            self.ForceTime = input.load.set_ramp_ForceTime(self.Time,self.TimeRamp)
        elif self.AppDynLoadVariation == 'sin':
            self.ForceTime = input.load.set_sin_ForceTime(self.Time,self.Omega)
        elif self.AppDynLoadVariation == 'cos':
            self.ForceTime = input.load.set_cos_ForceTime(self.Time,self.Omega)
        elif self.AppDynLoadVariation == 'omcos':
            self.ForceTime = input.load.set_omcos_ForceTime(self.Time,self.Omega)
        elif self.AppDynLoadVariation == 'rampsin':
            self.ForceTime = input.load.set_omcos_ForceTime(self.Time,self.TimeRamp,self.Omega) 
        elif self.AppDynLoadVariation == 'rect':
            self.ForceTime = input.load.set_rect_ForceTime(self.Time,self.TimeRectStart,self.TimeRectEnd)

        else:
            raise NameError('AppDynLoadVariation not valid!')           


    def set_prescr_vel(self):
        ''' set prescribed velocities (input.forcedvel)'''
     
        if self.PrVelType == 'sin':
            self.ForcedVel, self.ForcedVelDot = input.forcedvel.set_sin(self.Time, self.Omega, self.PrTrVelAmpl, self.PrRotVelAmpl)
        elif self.PrVelType == 'ramp':
            self.ForcedVel, self.ForcedVelDot = input.forcedvel.set_ramp(self.Time, self.TimeRamp, self.PrTrVelAmpl, self.PrRotVelAmpl) 
        elif self.PrVelType == 'rampsin':
            self.ForcedVel, self.ForcedVelDot = input.forcedvel.set_rampsin(self.Time, self.TimeRamp, self.Omega, self.PrTrVelAmpl, self.PrRotVelAmpl)         
        elif self.PrVelType == 'zero':
            self.ForcedVel, self.ForcedVelDot = input.forcedvel.set_zero(self.Time)             
        elif self.PrVelType == 'userdef':
            self.ForcedVel = np.array(self.ForcedVel, dtype=float, order='F')
            self.ForcedVelDot = np.array( self.ForcedVelDot, dtype=float, order='F')
        else:
            raise NameError('PrVelType not valid!')         


''' ------------------------------------------------------------------------ ''' 
    
if __name__=='__main__':
    
    from lib.plot import sta_unif
    import lib.save, lib.read
    
    xbsolver=XBeamSolver()
    
    xbsolver.TestCase='NCB1'
    ElemType="DISP"
    xbsolver.BConds='CF'
    
    xbsolver.NumElems = 10
    xbsolver.NumNodesElem = 3
    
    # Design
    xbsolver.beam_shape = 'constant'
    xbsolver.material ='alumA'
    xbsolver.cross_section_type ='isorect'  
    xbsolver.BeamLength1 = 5.0   
    xbsolver.cs_l2 = 0.1135
    xbsolver.cs_l3 = 0.2417

    # Loading
    xbsolver.AppStaLoadType='uniform'
    xbsolver.ExtForce[2]=6e5
    

    # Options - NCB1 case:
    xbsolver.FollowerForce=False
    xbsolver.PrintInfo=False              
    xbsolver.MaxIterations=30    
    xbsolver.NumLoadSteps=20
    xbsolver.Solution=112
    xbsolver.MinDelta= 1.e-5


    
    xbsolver.FollowerForceRig= True
    xbsolver.OutInBframe = True
    xbsolver.OutInaframe = False

    
    lib.save.h5comp(xbsolver, './component_before.h5')

    ## cost
    xbsolver.NNode = -1 # tip displacement

    xbsolver.execute() 
    print 'External force', xbsolver.ExtForce
    print 'Cost Functions:'
    print 'Total Mass', xbsolver.TotalMass, xbsolver.ct_TotalMass.value
    print 'Node Abs Displacement', xbsolver.NodeDisp, xbsolver.ct_NodeDisp.value 
    print 'Node z displacement', xbsolver.ZDisp
    
    # II run test
    xbsolver._use_previous_solution=True
    xbsolver.NumLoadSteps=1
    xbsolver.execute() 
    print 'External force', xbsolver.ExtForce
    print 'Cost Functions:'
    print 'Total Mass', xbsolver.TotalMass, xbsolver.ct_TotalMass.value
    print 'Node Abs Displacement', xbsolver.NodeDisp, xbsolver.ct_NodeDisp.value 
    print 'Node z displacement', xbsolver.ZDisp    
    
    # Make a nice plot
    # sta_unif(xbsolver.PosIni,xbsolver.PosDef,equal=True)
    sta_unif(xbsolver.PosIni,xbsolver.PosDef)
    
    # save:
    lib.save.h5comp(xbsolver, './component_after.h5')
    
    #### read test
    ##XBread=lib.read.h5comp('./component_after.h5')
    ##for tt in XBread.items():
    ##    print tt
        
    # check version
    print xbsolver._version.number
    print xbsolver._version.description
    
    
    
    
    
    

        
    
