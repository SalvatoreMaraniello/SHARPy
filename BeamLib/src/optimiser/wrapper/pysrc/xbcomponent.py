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
import design.beamelem, design.beamnodes
import beamvar
import cost
import lib.save

        
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
    _Solution_List = (102,112,142,202,212,302,312,900,910,902,912,922,952) # see header for possible values
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
    # general
    BeamLength1 = Float( 0.0, iotype='in', desc='')
    BeamLength2 = Float( 0.0, iotype='in', desc='')
    TipMass     = Float( 0.0, iotype='in', desc='')
    TipMassY    = Float( 0.0, iotype='in', desc='')
    TipMassZ    = Float( 0.0, iotype='in', desc='')
    Omega       = Float( 0.0, iotype='in', desc='') 

    # beam span reconstruction
    beam_shape = Enum('constant',('constant','spline','constant_hardcoded'),iotype='in')
 
    # isotropic material cross-sectional models(see design.isosec)
    material = Enum( 'alumA',('alumA','test'), iotype='in', desc='Material for isotropic cross ection (see design.isosec)') 
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

    # Cost/Constrains
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
        
        # Applied external forces at the Tip
        self.ExtForce      = np.zeros(  (3), dtype=float, order='F')
        self.ExtMomnt      = np.zeros(  (3), dtype=float, order='F')  
          
        self.BeamSpanStiffness = np.zeros( (self.NumElems,6,6), dtype=float, order='F' )
        self.BeamSpanMass      = np.zeros( (self.NumElems,6,6), dtype=float, order='F' )  
        
        # for constant_hardcoded method only
        self.Mroot = np.zeros((6,6),dtype=float,order='F')
        self.Kroot = np.zeros((6,6),dtype=float,order='F')
    
        # save directory
        self._savedir = '.'
    
        # multiple runs
        self._counter = 0
        self._use_previous_solution = False
        

    def pre_solve(self):
        
        # ------------------------------------------------ Determine Arrays Size
        self.NumNodes = beamvar.total_nodes(self.NumElems,self.NumNodesElem)
        self._NumNodes_copy = self.NumNodes  
                       

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
        # General
        self.DensityVector = np.zeros((self.NumElems),dtype=float,order='F')
        self.LengthVector  = np.zeros((self.NumElems),dtype=float,order='F')
        
        # Problem dependent
        
        ### sm 12 sep 2014:
        # Use old solution for initial guess in the optimisation. This if 
        # statement only prevents self.{Psi/Pos}Def not to be overwritten as
        # the self.{Pos/Psi}Ini are allocated here but always overwritten inside 
        # the Fortran code.
        if self._counter < 1 or self._use_previous_solution is False:
            # self.PosDef, self.PsiDef are set to zero
            # all others are empty arrays
            [self.PosIni, self.PosDef, self.PsiIni, self.PsiDef, self.InternalForces] = beamvar.fwd_static(self.NumNodes,self.NumElems) 
             

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
        #print 'PhiNodes', self.PhiNodes
        #print 'BeamMass', self.BeamSpanMass
        #print 'BeamStiff'
        #for ii in range(self.NumElems):
        #    print np.diag(self.BeamSpanStiffness[ii,:,:])
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
                     ct.byref(self.ct_BeamLength1), ct.byref(self.ct_BeamLength2),
                     self.BeamSpanStiffness.ctypes.data_as(ct.c_void_p),        # Properties along the span
                     self.BeamSpanMass.ctypes.data_as(ct.c_void_p),                       
                     self.PhiNodes.ctypes.data_as(ct.c_void_p),                 # PreTwist angle along the span
                     self.ExtForce.ctypes.data_as(ct.c_void_p), self.ExtMomnt.ctypes.data_as(ct.c_void_p),
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
        
        # compute cost from Fortran functions
        cost.f_total_mass( self.DensityVector.ctypes.data_as(ct.c_void_p), self.LengthVector.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumElems), ct.byref(self.ct_TotalMass) )
        cost.f_node_disp( self.PosIni.ctypes.data_as(ct.c_void_p), self.PosDef.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumNodes), ct.byref(self.ct_NodeDisp), ct.byref(self.ct_NNode) )
        
        # update values
        self.update_variables()
        
        # compute cost from Python functions
        self.ZDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='z',NNode=self.NNode-1) # python indices starts from 0
        self.YDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='y',NNode=self.NNode-1)
        self.XDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='x',NNode=self.NNode-1)
        
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
        
        self.ct_ElemType     = ct.c_char_p( self.ElemType.ljust(4) )     # the ljust is not required
        self.ct_TestCase     = ct.c_char_p( self.TestCase.ljust(4) )
        self.ct_BConds       = ct.c_char_p( self.BConds.ljust(2)   )
        
        self.ct_BeamLength1 = ct.c_double( self.BeamLength1 )
        self.ct_BeamLength2 = ct.c_double( self.BeamLength2 )
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
        self.BeamLength2 = self.ct_BeamLength2.value
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
    xbsolver.cs_l2 = 0.15
    xbsolver.cs_l3 = 0.15

    # Loading
    xbsolver.ExtForce[2]=6e5

    # Options - NCB1 case:
    xbsolver.FollowerForce=False
    xbsolver.PrintInfo=False              
    xbsolver.MaxIterations=99    
    xbsolver.NumLoadSteps=20
    xbsolver.Solution=142
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
    
    # read test
    XBread=lib.read.h5comp('./component_after.h5')
    for tt in XBread.items():
        print tt
    
    
    
    
    

        
    
