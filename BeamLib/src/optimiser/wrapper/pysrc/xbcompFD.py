'''

Salvatore Maraniello 19 Feb 2015

OpenMDAO component for Beam solver with FD parallel gradient

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
import h5py 
import multiprocessing as mpr
from warnings import warn

sys.path.append('..')
import shared
from openmdao.main.api import ComponentWithDerivatives
from openmdao.main.datatypes.api import Float, Array, Enum, Int, Str, Bool


import design.beamelem, design.beamnodes, design.beamgeom
import input.load.modal, input.load.fourier, input.load.spline, input.forcedvel
import beamvar
import cost, constr.common
import lib.save
      
# Load Dynamic Library
wrapso = ct.cdll.LoadLibrary(shared.wrapso_abspath)



class XBeamSolver(ComponentWithDerivatives):
 
    #-------------------------------------------------------- OpenMDAO Interface

    '''Shared Input'''
    NumElems     = Int(10, iotype='in', desc='Total number of elements')
    NumNodesElem = Enum(3, (2,3), iotype='in', desc='Number of nodes per element; (2) linear, (3) quadratic.')
    ElemType     = Enum('DISP',('DISP','STRN'), iotype='in')
    
    # BCs available:
    # - CF, CC: static/dynamics/rigid-body
    # - CH, HC, HH: static only
    # - HF: non-linear coupled rigid body + structural
    BConds       = Enum('CF'  ,('FC','CF','CC','CS','SC','SS','SF','FS','MS','MC'), iotype='in', 
                        desc='Boundary conditions: C (clamped end), F (free end), S (spherical joint). MS, MC applies a spherical joint or clamp in the middle')   
    TestCase = Str('test', iotype='in', desc='Test case name.')

    '''Options:'''
    FollowerForce   = Bool( True, iotype='in', desc='Follower Force; default: True')  
    FollowerForceRig= Bool( True, iotype='in', desc='Follower force in the body-fixed frame; default: True')  
    PrintInfo       = Bool( True, iotype='in', desc='Print info on screen')  
    OutInBframe     = Bool( True, iotype='in', desc='print velocities in B-frame (if not, use a-frame). default: True')  
    OutInaframe     = Bool(False, iotype='in', desc='print velocities in a-frame (if not, Inertial frame). default: False') 
      
    ElemProj        = Enum( 0, (0,1,2), iotype='in', desc='Element info computed in (0) global frame, (1) fixed element frame, (2) moving element frame.')  
    _Solution_List = (102,112,142,202,212,302,312,900,910,902,912,922,932,952) # see header for possible values
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
    material = Enum( 'alumA',('alumA','titnA','test','qiqi1'), iotype='in', desc='Material for isotropic cross section (see design.isosec)') 
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
    cs_pretwist = Float(0.0,iotype='in',desc='beam pre-twist')
    
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
    DynForceParam = Enum('modal', ('modal','fourier','spline'), desc='Method to model the dynamic force (see self.set_load related modules)' )
    ExtForceDyn = Array(default_value=np.zeros((3),order='F'), # order='F' is not necessary
                     iotype='in', dtype=float, shape=(3,) ,
                     desc='Nodal applied force - single node or distributed according to self.AppSDynLoadType (load.set_***_ForceDynAmp)')
    ExtMomntDyn = Array(default_value=np.zeros((3),order='F'), # order='F' is not necessary
                     iotype='in', dtype=float, shape=(3,) ,
                     desc='Nodal applied moments - single node or distributed according to self.AppDynLoadType (load.set_***_ForceDynAmp)')   
    DynForceSpecialOptions = {'reverse':False,'smooth':False,'Twin':0.0}
    '''# spline control
    DynFrc_spline_order=Int(3, iotype='in', desc='order of B-spline control parametrisation')
    TCint = Array( iotype='in', dtype=float, desc='TCint[:,cc] contains the NC+1 knots points for generating the spline representation over the cc-th force/moment component')
    Scf = Array( iotype='in', dtype=float, desc='spline weights related to the NS = spline_order+NC')  
    '''  
    '''# fourier control
    DynForceSpecialOptions = {'reverse':False,'smooth':False,'Twin':0.0}
    Acf = Array(default_value=np.zeros((6),order='F'), iotype='in', dtype=float, desc='Sin related to Fourier coefficients (see design.load.fourier)')   
    Bcf = Array(default_value=np.zeros((6),order='F'), iotype='in', dtype=float, desc='Cin related to Fourier coefficients (see design.load.fourier)')
    Fcf = Array(default_value=np.zeros((6),order='F'), iotype='in', dtype=float, desc='Frequencies related to Fourier coefficients (see design.load.fourier)')
    '''
 
    ''' Optimisation '''
    # this allows to define costumised optimisation functions. The output value 
    # is stored in fval (scalar). The function is ffun and the arguments are 
    # in fargs. 
    # In some situations, the pre_solve method may have to be run before being
    # able to define the arguments.
    fval = Float( 0.0, iotype='out', desc='generic cost function')
    geqval=Array(iotype='out', desc='value of equality constraints')
    gdisval=Array(iotype='out', desc='value of disequality constraints')
    
    # Specific Cost/Constrains [superceeded]
    TotalMass =  Float( 0.0, iotype='out', desc='Total structural mass')
    NodeDisp  =  Float(  iotype='out', desc='Absolute value of NNode displacement')
    ZDisp     =  Float(  iotype='out', desc='z displacement of NNode')
    YDisp     =  Float(  iotype='out', desc='y displacement of NNode')
    XDisp     =  Float(  iotype='out', desc='x displacement of NNode')
    NNode     =  Int(-1,iotype='in', desc='Node at which monitor the displacement. If negative, the last node is picked up')    
     
 
     
 
    def __init__(self):
        
       
        # not to overwrite Component __init__
        super(XBeamSolver, self).__init__()
        
        self.missing_deriv_policy = 'error' 
        #self.missing_deriv_policy = 'assume_zero' 
        #self.missing_deriv_policy = 'fd' #NOTE: This will be developed later        
        
        # extract main routine
        ### method below fails inside a class but runs from main
        ###fwd_run=xbso.__opt_routine_MOD_opt_main 
        self.fwd_run = getattr(wrapso, "__opt_routine_MOD_opt_main")       
        self.SUCCESS=True # handling exceptions (sol 932 only)

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
        
        
        # initial conditions
        
        # initial a FoR orientation (0 deg.s default)
        self.Quat0 = np.array( [1, 0, 0, 0], dtype=float, order='F')
        # initial a FoR velocity in inertial FoR G
        self.RefVel0 = np.zeros((1,6), dtype=float, order='F')
        self.RefVelDot0 = np.zeros((1,6), dtype=float, order='F')
        # DynOut0: initial deformation
        # note that each NumNodes rows of DynOut are equal to PosDef at that time-step
        # because the restart is only applicable to dynamic cases, the variable to specify
        # the initial beam deformation has been labelled DynOut0 and not PosDef0.
        # Note that DynOut0 requires the knowledge of the number of nodes in the
        # beam, so assumes a restart.
        ##self.DynOut0 = np.array((1,3),dtype=float,order='F')

        # Initial nodal position 
        self.PosDef0 = np.zeros( (1,3), dtype='float', order='F' )
        # Initial element orientation vectors
        self.PsiDef0 = np.zeros( (1,3,3), dtype='float', order='F' ) # the 2nd index is the max nodes per element

        # Initial nodal position velocity
        self.PosDotDef0 = np.zeros( (1,3), dtype='float', order='F' )
        # Initial element orientation velocity vectors
        self.PsiDotDef0 = np.zeros( (1,3,3), dtype='float', order='F' ) # the 2nd index is the max nodes per element
        
        
        #------------------------------------------------------------- Utilities
        # save directory
        self._savedir = '.'
    
        # multiple runs
        self._counter = 0
        self._use_previous_solution = False
        
        # restart option
        self._restart= False
        self._restart_file=''
        
        # append code version
        self._version = shared.CodeVersion()    
        
        # Parallel FDs
        self.PROCESSORS = 1       # number of processors
        self.parallelFDs = False # FLAG to set parallel calculations
        self._gradmode=False # if true, only computes gradient and then stops
        
        #------------------------------------------------------------------ Cost
        # dummy arguments for generic cost function
        #self.ffun = cost.f_xyz_disp
        #self.fargs=['PosIni','PosDef']
        #self.fname = self.ffun.func_name
        
        # dummy arguments for constraint functions
        #self.geqfun  =[]
        #self.geqargs =[]
        #self.gdisfun =[]
        #self.gdisargs=[]


    def run_custom(self):
        ''' 
        Allows to run custom operations. 
         - connect new design variables to attributes in the XBeamSolver class, 
           e.g. if only a portion of an array is required as design variable
         - implement easy constraint, thus reducing the amount of design 
         variables (e.g. if cs_l2*cs_l3 = 0.10, cs_l3 can be computed from 
         cs_l2 each time)
        '''
        return self


    def pre_solve(self):
        
        # ------------------------------------------------ run custom operations
        self.run_custom()
        
        
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
            shift=0.0  
            if 'S' in self.BConds:
                if self.BConds=='MS':
                    shift=-0.5*self.BeamLength1
                elif self.BConds=='FS':
                    shift=-self.BeamLength1
                
            self.PosIni = design.beamgeom.straight(self.PosIni,self.BeamLength1,self.BeamAxis,shift)



        #--------------------------------------- Set output & initial conditions
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


        # Set Rigid-Body + Dynamics input/output
        # initial values can be added to displacements and velocities
        self.RefVel    = np.zeros( (self.NumSteps+1,6), dtype='float', order='F')
        self.RefVelDot = np.zeros( (self.NumSteps+1,6), dtype='float', order='F')
        
        
        #------------------------------------------------------ set general ICs:

        # a FoR 
        self.Quat =np.array(self.Quat0,dtype=float,order='F')
        self.RefVel[0,:]=self.RefVel0
        self.RefVelDot[0,:]=self.RefVelDot0
        ### Fortran initialisation discarded
        ### self.RefVel[0,:] = self.ForcedVel[0,:]          
        ### self.RefVelDot[0,:]=self.ForceVelDot[0,:] # discarded: acceleration not a valid IC for a II order problem            

        if self._restart==True:
            # initial deformation (in a FoR)
            # WARNING: this overrides self.PoseDef defined above           
            ##self.DynOut[:self.NumNodes,:]=self.DynOut0 ## <-- does nothing!
            self.PosDef=np.array(self.PosDef0, dtype='float', order='F')
            self.PsiDef=np.array(self.PsiDef0, dtype='float', order='F')
            self.PosDotDef=np.array(self.PosDotDef0, dtype='float', order='F')
            self.PsiDotDef=np.array(self.PsiDotDef0, dtype='float', order='F')


        # prepare saving directory:
        self._savedir=os.path.abspath(self._savedir)
        if os.path.isdir(self._savedir) is False:
            os.makedirs(self._savedir)
            
        return self
        

    def execute(self):
        
        # set size etc.
        self.pre_solve()

        # set ctypes for Fortran interface
        self.set_ctypes()
             
        # save before run (in case of crash)
        counter_string =  '%.3d' % (self._counter)
        savename = self._savedir + '/' + self.TestCase + '_' + counter_string + '.h5'
        
        try:  
            lib.save.h5comp( self, savename)
            self.call_fwd_run()
        except NameError:
            if self._counter<1:
                # optimiser repeats twice the initial conditions, under runs 0 and 1.
                # However, as the IC is identical, an error will be raised for counter=0
                raise NameError('xbeam solver did not converge for the IC set!')
            else:
                print 'xbeam solver crashed: solution will be restarted after decreasing the step size'
                
                # update design variables values
                red_fact=np.array([0.5,0.25,0.125,0.0625])
                
                # get value of design variables set by the optimiser
                DesignValCrash = self.get_DesignVal()
                Ndes=len(self.DesignVal)
                
                # with auto FD self.DesignValOld is not retained correctly
                # it must be re-read
                old_filename = self._savedir + '/' + self.TestCase + '_%.3d.h5'%(self._counter-1)
                hdfile = h5py.File(old_filename,'r')
                DesignValOld = hdfile['DesignVal'].value
                #val = hdfile['DesignVal'].value
                #setattr(self,'DesignValOld',val)    #<--- 
                hdfile.close()

                # and step size set by the optimiser
                #DesignStepCrash = [(DesignValCrash[ii] - self.DesignValOld[ii]) for ii in range(Ndes)]
                DesignStepCrash = [(DesignValCrash[ii] - DesignValOld[ii]) for ii in range(Ndes)]
                
                # update design variable values
                for rr in range(len(red_fact)):
                    print 'Restart with factor %1.4f' %(red_fact[rr])
                    DesignStepNew = [red_fact[rr]*DesignStepCrash[ii] for ii in range(Ndes)]
                    #DesignValNew=[(self.DesignValOld[ii]+DesignStepNew[ii]) for ii in range(Ndes)]
                    DesignValNew=[(DesignValOld[ii]+DesignStepNew[ii]) for ii in range(Ndes)]
                    self.set_DesignVal(DesignValNew)
                    
                    try:
                        # reset all the variables
                        self.pre_solve()
                        self.del_ctypes()
                        self.set_ctypes()
                        lib.save.h5comp( self, savename)
                        # rerun the problem
                        self.call_fwd_run()
                        break
                    except NameError:
                        if rr==len(red_fact)-1:
                            raise NameError('Even with reduction factor of %1.4f, the run has crashed at iteration %d!' %(red_fact))
                    finally:
                        # saves. Overrides previous reduction steps if more then one is executed
                        saverep = self._savedir + '/' + self.TestCase + '_' +counter_string + '_REDUCEDSTEPREPORT' + '.h5'
                        hdfile = h5py.File(saverep ,'w')
                        
                        DesignTpl ,FunctionalsTpl = self.list_deriv_vars()
                        DesignNames=list(DesignTpl)
                        
                        for dd in range(Ndes):
                            hdfile[DesignNames[dd]+'_old']  =DesignValOld[dd] #self.DesignValOld[dd]
                            hdfile[DesignNames[dd]+'_crash']=DesignValCrash[dd]
                            hdfile[DesignNames[dd]+'_step_crash']=DesignStepCrash[dd]
                            hdfile[DesignNames[dd]+'_new']=DesignValNew[dd]
                            hdfile[DesignNames[dd]+'_step_new']=DesignStepNew[dd]
                            
                        hdfile['red_fact']=red_fact[rr]
                        hdfile.close()               
                        
        # output checks
        self.check()
        
        # compute cost from Fortran functions [superceeded]
        cost.f_total_mass( self.DensityVector.ctypes.data_as(ct.c_void_p), self.LengthVector.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumElems), ct.byref(self.ct_TotalMass) )
        cost.f_node_disp( self.PosIni.ctypes.data_as(ct.c_void_p), self.PosDef.ctypes.data_as(ct.c_void_p), ct.byref(self.ct_NumNodes), ct.byref(self.ct_NodeDisp), ct.byref(self.ct_NNode) )
        
        # update values
        self.update_variables()
        
        # compute cost from Python functions [superceeded]
        self.ZDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='z',NNode=self.NNode-1)
        self.YDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='y',NNode=self.NNode-1)
        self.XDisp = cost.f_xyz_disp(self.PosIni, self.PosDef, dir='x',NNode=self.NNode-1)
        
        
        # evaluate cost, constraints and design
        self.eval_functionals()      
        try:
            ####self.DesignValOld = self.get_DesignVal() # not retained with FD gradient
            self.DesignVal = self.get_DesignVal()        # this will be saved in h5 file
        except AttributeError:
            warn('One or more design variables not found!!!')
            self.DesignVal=[]
                        
        # sm: savename etc moved before solver started
        ### save
        #counter_string =  '%.3d' % (self._counter)
        #savename = self._savedir + '/' + self.TestCase + '_' + counter_string + '.h5'
        lib.save.h5comp( self, savename)
        
        # counter: for multiple executions
        self._counter = self._counter+1
        # special exit (FD convergence test) 
        if self._counter > 2 and self._gradmode==True:
            print('!!! Gradient mode finished! !!!') 
            1/0 # only way to exit from OpenMDAO. Use try/except in main file
            
        

    def call_fwd_run(self):
        '''
        Simple wrap function to store the full call of the fortran code.
        
        IMPORTANT:
        though most of the variables in input to self.fwd_run are defined
        in the Fortran original code as optional, not passing them will generate
        memory errors.
        
        '''
        
        if self._gradmode == True and self._counter>1: 1/0
        
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
                     self.ForceDynamic.ctypes.data_as(ct.c_void_p),                                                # new input dynforce
                     self.ForcedVel.ctypes.data_as(ct.c_void_p), self.ForcedVelDot.ctypes.data_as(ct.c_void_p),    # input_forcedvel
                     self.PosDotDef.ctypes.data_as(ct.c_void_p), self.PsiDotDef.ctypes.data_as(ct.c_void_p),
                     self.PosPsiTime.ctypes.data_as(ct.c_void_p), self.VelocTime.ctypes.data_as(ct.c_void_p),
                     self.DynOut.ctypes.data_as(ct.c_void_p)  ,                                              # output from sol 202, 212, 302, 312, 322
                     self.RefVel.ctypes.data_as(ct.c_void_p), self.RefVelDot.ctypes.data_as(ct.c_void_p), 
                     self.Quat.ctypes.data_as(ct.c_void_p)   ,                # output rigid-body + dynamics
                     ct.byref(self.ct_SUCCESS) ) # handlying exceptions for sol 932
        
        '''
        class FwdNotConverged(Exception):
            def __init__(self, value):
                self.SUCCESS = True
            def __str__(self):
                return repr(self.value)
        '''
        if self.ct_SUCCESS.value is False:
            #raise  FwdNotConverged(self.ct_SUCCESS.value)
            raise NameError('Sol 932 did not converge')

                          
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
        
        #------------------------------------ handling exceptions (sol 932 only)
        self.ct_SUCCESS = ct.c_bool(self.SUCCESS)


    def del_ctypes(self,printinfo=False):
        
        # list attributes
        ListNames=dir(self) # lists attributes and methods
        for ll in ListNames:
            if ll[:3]=='ct_':
                delattr(self,ll)
                if printinfo==True:
                    print 'Attribute %s deleted!' %(ll)
                
        '''
        # -------------------------------------------------------- Prepare input
        # arrays do not require any modification
        self.ct_NumElems     = None
        self.ct_NumNodesElem = None
        self.ct_NumNodes     = None
        self.ct_NumSteps     = None
        
        self.ct_ElemType     = None    # the ljust is not required
        self.ct_TestCase     = None
        self.ct_BConds       = None
        
        self.ct_BeamLength1 = None
        self.ct_TipMass     = None
        self.ct_TipMassY    = None
        self.ct_TipMassZ    = None
        self.ct_Omega       = None
        
        # ------------------------------------------------------------- Prepare Options
        self.ct_FollowerForce    = None
        self.ct_FollowerForceRig = None
        self.ct_PrintInfo        = None
        self.ct_OutInBframe      = None
        self.ct_OutInaframe      = None 
        
        self.ct_ElemProj      = None
        self.ct_MaxIterations = None
        self.ct_NumLoadSteps  = None
                 
        self.ct_Solution      = None
        
        self.ct_DeltaCurved   = None
        self.ct_MinDelta      = None
        self.ct_NewmarkDamp   = None
        
        #---------------------------------------------------------- Prepare Cost
        self.ct_TotalMass = None
        self.ct_NodeDisp  = None
        self.ct_NNode     = None
        
        #------------------------------------ handling exceptions (sol 932 only)
        self.ct_SUCCESS = None    
        
        '''

        
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
        # Static Load always use shape methods from input.load.modal
        # Span Distribution of Prescribed Load
        if self.AppStaLoadType == 'uniform':
            print 'setting static uniform force'
            self.ForceStatic = input.load.modal.set_unif_ForceDistr(
                                self.NumNodes, 
                                self.ExtForce, self.ExtMomnt)    
        elif self.AppStaLoadType == 'nodal':
            self.ForceStatic = input.load.modal.set_nodal_ForceDistr(
                                self.NumNodes, self.NodeAppStaForce,
                                self.ExtForce, self.ExtMomnt)
        elif self.AppStaLoadType == 'userdef':
            # make sure the array order is 'F'
            self.ForceStatic = np.array(self.ForceStatic, dtype=float, order='F')
        else:
            raise NameError('AppStaLoadType not valid!')
        
            
        if self.DynForceParam=='modal': 
            # In this method, loads are decomposed into a:
            #   - shape: self.ForceDynAmp - only space dependent
            #   - amplitude: self.ForceTime - only time dependent 
            # such that the force at time-step tt on the node nn is:
            # F(nn,tt) = shape(nn) * amplitude(tt)
            
            # Span Distribution of Prescribed Load           
            if self.AppDynLoadType == 'uniform':
                self.ForceDynAmp = input.load.modal.set_unif_ForceDistr(
                                       self.NumNodes, self.ExtForceDyn, self.ExtMomntDyn)    
            elif self.AppDynLoadType == 'nodal':
                self.ForceDynAmp = input.load.modal.set_nodal_ForceDistr(
                                    self.NumNodes, self.NodeAppDynForce,
                                    self.ExtForceDyn, self.ExtMomntDyn)
            elif self.AppDynLoadType == 'userdef':
                self.ForceDynAmp = np.array(self.ForceDynAmp, dtype=float, order='F')
            else:
                raise NameError('AppDynLoadType not valid!')
                
            # Time variation of Prescribed Load
            if self.AppDynLoadVariation == 'constant':
                self.ForceTime = input.load.modal.set_const_ForceTime(self.Time)
            elif self.AppDynLoadVariation == 'impulse':
                self.ForceTime = input.load.modal.set_impulse_ForceTime(self.Time)
            elif self.AppDynLoadVariation == 'ramp':
                self.ForceTime = input.load.modal.set_ramp_ForceTime(self.Time,self.TimeRamp)
            elif self.AppDynLoadVariation == 'sin':
                self.ForceTime = input.load.modal.set_sin_ForceTime(self.Time,self.Omega)
            elif self.AppDynLoadVariation == 'cos':
                self.ForceTime = input.load.modal.set_cos_ForceTime(self.Time,self.Omega)
            elif self.AppDynLoadVariation == 'omcos':
                self.ForceTime = input.load.modal.set_omcos_ForceTime(self.Time,self.Omega)
            elif self.AppDynLoadVariation == 'rampsin':
                self.ForceTime = input.load.modal.set_omcos_ForceTime(self.Time,self.TimeRamp,self.Omega) 
            elif self.AppDynLoadVariation == 'rect':
                self.ForceTime = input.load.modal.set_rect_ForceTime(self.Time,self.TimeRectStart,self.TimeRectEnd)
    
            else:
                raise NameError('AppDynLoadVariation not valid!')  
        
            # Build the final vector: DO NOT ADD STATIC PART
            # in some solutions (e.g.312) the static solver is called first
            self.ForceDynamic =  input.load.modal.assembly(self.ForceDynAmp,self.ForceTime)
        

        if self.DynForceParam=='fourier': 
            # detect special options
            
            self.ForceDynAmp = np.zeros( (self.NumNodes,6), dtype='float', order='F')
            self.ForceTime = np.ones( (self.NumSteps+1), dtype='float', order='F' )
            self.ForceDynamic =  input.load.fourier.glsin_nodal(
                                self.NumNodes, self.Time, self.NodeAppDynForce, 
                                self.Acf, self.Fcf, optdict=self.DynForceSpecialOptions)
            
        if self.DynForceParam=='spline': 
            self.ForceDynAmp = np.zeros( (self.NumNodes,6), dtype='float', order='F')
            self.ForceTime = np.ones( (self.NumSteps+1), dtype='float', order='F' )
            #spline_nodal(NumNodes, Time, NodeForce, TCint, Scf, p )
            self.ForceDynamic =  input.load.spline.nodal_force(
                                self.NumNodes, self.Time, self.NodeAppDynForce, 
                                self.TCint, self.Scf, self.DynFrc_spline_order)            
            

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


    def eval_functionals(self):
        '''
        Evaluate the cost and constraintsfunctionals given in the lists
        self.ffun, self.geqfun and self.gdisfun
        '''
        # compute general cost
        if hasattr(self,'ffun'):
            #self.fval = self.ffun(*self.fargs)
            tplargs = cost.return_args(self, self.fargs)
            self.fval = self.ffun(*tplargs)
            self.fname = self.ffun.func_name
        
        # compute equality constraints
        if hasattr(self, 'geqfun'):
            #if len(self.geqfun)>0:
            geqval_list=[]
            self.geqname=[]
            Ngeq = len(self.geqfun)
            for gg in range(Ngeq):
                ggeqfun  = self.geqfun[gg]
                ggeqargs = self.geqargs[gg]
                #geqval_list.append( ggeqfun(*constr.common.return_args(self, ggeqargs)) )
                ggeqval = ggeqfun(*constr.common.return_args(self, ggeqargs))
                if np.isscalar(ggeqval): # output is a scalar
                    geqval_list.append( ggeqval )
                else: # output is an array
                    Nelems = np.prod(ggeqval.shape)
                    gvec = np.reshape( ggeqval, ( Nelems ,) )
                    for gelem in gvec:
                        geqval_list.append( gelem )                
                self.geqname.append( ggeqfun.func_name )
            self.geqval=np.array(geqval_list)
            #print 'eq constr: ', self.geqval
        
        # compute disequality constraints
        if hasattr(self, 'gdisfun'):
            #if len(self.gdisfun)>0:
            gdisval_list=[]
            self.gdisname=[]
            Ngdis = len(self.gdisfun)
            for gg in range(Ngdis):
                ggdisfun  = self.gdisfun[gg]
                ggdisargs = self.gdisargs[gg]
                #gdisval_list.append( ggdisfun(*constr.common.return_args(self, ggdisargs)) )
                ggdisval = ggdisfun(*constr.common.return_args(self, ggdisargs))
                if np.isscalar(ggdisval): # output is a scalar
                    gdisval_list.append( ggdisval )
                else: # output is an array
                    Nelems = np.prod(ggdisval.shape)
                    gvec = np.reshape( ggdisval, ( Nelems ,) )
                    for gelem in gvec:
                        gdisval_list.append( gelem )                         
                self.gdisname.append( ggdisfun.func_name ) 
            self.gdisval = np.array(gdisval_list)       
            #print 'dis constr:', self.gdisval
        
        return


    def list_deriv_vars(self):
        """
        Specified the inputs and outputs where derivatives are defined
        
        Best practice: define as list and convert to tuple to aovid cases where
        t = ('one_variable') is a string, not a tuple.
        
        """
        
        Design = ['']
        Functionals = ['']
        
        return tuple(Design), tuple(Functionals)


    def provideJ(self):
        '''Note that the function will fail if the variables are not 1d arrays or
        scalar'''
        
        # get list of variables/functionals
        DesignList ,FunctionalsList = self.list_deriv_vars()

        # extract length of each element
        def get_variables_size(VarList):
            '''
            Given a list/tuple of N variables, scalar or 1d arrays, the 
            function returns a N length array containing the length of each 
            variable
            
            Reminder: the method doesn't work with 2D arrays!
            '''
            
            # check: if VarList=('anystring'), varList will be saved as string.
            # the code will crash unless this is redifined as tuple
            if isinstance(VarList,str) is True:
                raise NameError('VarList detected as a string, not a tuple: check list_deriv_vars!')
                
            Nv = len(VarList)
            VarLen=np.ones(Nv,dtype=np.int64)
            
            for vv in range(Nv):
                # get the variable
                var = getattr(self,VarList[vv]) 
                if np.isscalar(var) is False:
                    if len(var.shape)>1:
                        raise NameError('FD Jacobian can only handle 1D arrays!')
                    else:
                        VarLen[vv]=len(var)        
        
            return VarLen
        
        DesLen = get_variables_size(DesignList)
        FunLen = get_variables_size(FunctionalsList)
        ND = len(DesLen)
        
        # allocate jacobian
        Nx = np.sum(DesLen)
        Nf = np.sum(FunLen)
        J = np.zeros((Nf,Nx))
        
        # check if fd_steps have been added
        if hasattr(self,'fd_steps') is False:
            raise NameError('self.fd_steps not found! Add it in the input file using self.get_fd_steps()!')
        else:
            dx = self.fd_steps
        
        # -------------------------------------------------------- start looping  
        
        # get global index for design variables: Svec[dd] contains starting 
        # index for the dd-th variable (if dd-th variable is an array, this is
        # the index of the first element 
        Svec = np.zeros((ND,),dtype=np.int64) 
        for dd in range(ND):
            Svec[dd]=np.sum(DesLen[:dd])

  
        if self.parallelFDs is False:
            # ------------------------------------------------------serial code:
            for dd in range(ND):
                for ii in range(DesLen[dd]):
                    jvii = self.perturb(self.copy(),dx,dd,ii,DesignList,FunctionalsList,DesLen,FunLen)
                    J[:,Svec[dd]+ii]=jvii
        
        else:
            # ---------------------------------------------------- parallel code:
            
            # create pool
            self.fwd_run=None
            self.del_ctypes(printinfo=False)
            pool = mpr.Pool(processes=self.PROCESSORS) 
            # send the jobs
            results=[]
            for dd in range(ND):
                for ii in range(DesLen[dd]):        
                    results.append( pool.apply_async(self.perturb,args=( self.copy(),dx,dd,ii,DesignList,FunctionalsList,DesLen,FunLen ) ))              
            # and retrieve the results!
            jvList = [p.get() for p in results]
            if len(jvList)!=Nx:
                raise NameError('Debugging: Number of jacobian vectors computed with Parallel FDs wrongly equal to %d!' %(len(jvList)))
            for jj in range(Nx):
                J[:,jj]=jvList[jj]
            # - 1. close the pool (memory in workers goes to zero) 
            # - 2. exit the worker processes (processes are killed)
            pool.close()
            pool.join() 
            # reset the ctypes for future executions
            self.fwd_run=getattr(wrapso, "__opt_routine_MOD_opt_main") 
            self.set_ctypes()            
            
            
        # --------------------- save jacobian data in the next forward execution
        self.Jacobian = J
        self.DesignList = DesignList
        self.FunctionalsList = FunctionalsList
        # resave the whole run 
        counter_string =  '%.3d' % (self._counter-1)
        savename = self._savedir + '/' + self.TestCase + '_' + counter_string + '.h5'
        lib.save.h5comp( self, savename)


        return J
        

    def perturb(self,xbpert,dx,dd,ii,DesignList,FunctionalsList,DesLen,FunLen):
        ''' Computes (dfdx_ii, dgdx_ii) '''
        
        # rebuild link to Fortran library
        if xbpert.parallelFDs: 
            xbpert.fwd_run=getattr(wrapso, "__opt_routine_MOD_opt_main") 
            xbpert.set_ctypes()
        
        # get unperturbed functionalsallocate unperturbed functionals
        def put_functionals_in_array(xbpert,FunctionalsList,FunLen):
            
            Nlist = len(FunctionalsList)
            Nf = np.sum(FunLen)
            jv = np.zeros((Nf,))
                
            for ll in range(Nlist):
                # functional value
                varname=FunctionalsList[ll]
                f0 = getattr(xbpert,varname)
                #get starting index
                gg = np.sum(FunLen[:ll])
                if np.isscalar(f0):
                    jv[gg]=f0
                else:
                    jv[gg:gg+FunLen[ll]]=f0
                        
            return jv
            
        jv0 = put_functionals_in_array(xbpert,FunctionalsList,FunLen)        
  
        # get variable and starting index in global ordering
        varname = DesignList[dd]
        x0 = getattr(xbpert,varname)  # array or scalar   
            
        ss = np.sum(DesLen[:dd])
            
        if np.isscalar(x0):
            if ii != 0:
                raise NameError('For scalar, dimension should be 1')
            setattr(xbpert,varname,x0 + dx[ss+ii]) # or dx[ss], is the same
        else:
            xpert=x0
            xpert[ii]=x0[ii]+dx[ss+ii]
            setattr(xbpert,varname,xpert)
        
        # reallocate fwd_method - not copied with pool, overwritten for serial
        xbpert.fwd_run=getattr(wrapso, "__opt_routine_MOD_opt_main")  
        xbpert.TestCase=xbpert.TestCase+'_FD%.3d' %(ss+ii)   
        xbpert._counter=xbpert._counter-1                  
        xbpert.execute()
            
        # Finite Difference
        jv = ( put_functionals_in_array(xbpert,FunctionalsList,FunLen) - jv0 )/dx[ss+ii]
            
        return jv


    def get_fd_steps(self,AssemblyObj):
        '''This method can only be used once the XBeamSolver is added to a
        driver workflow. The Assembly class needs to be passed in input.
        The method has to be called in the input setup
        '''
        
        self.fd_steps=AssemblyObj.driver.get_fd_steps()
        
        # check fd step
        if any(self.fd_steps==0):
            print 'FD array:'
            print self.fd_steps
            raise NameError('Detected 0 FDs step! Remind to define the FD step inside the input file!')
        if any(np.isnan(self.fd_steps)):
            print 'FD array:'
            print self.fd_steps
            raise NameError('Detected NaN in FDs step! Remind to define the FD step inside the input file!')
 
        
        return self
                
                
    def get_DesignVal(self):
        '''
        Gets current value of design vector. These are not stored as an array
        but as a list, whose ii-th element is related to the ii-th design variable
        (being this a float or array) defined using self.list_deriv_vars()
       
        The choice of storing the results in list format is due to the fact that
        this method is not used for computing the jacobian, but just for resetting
        the optimiser step in case of crash during the optimisation
        
        Remark:
        the method is analoguous to SLSQPdriver.eval_parameters
        
        ''' 
        
        # get tuples of variables/functionals
        DesignTpl ,FunctionalsTpl = self.list_deriv_vars()
        DesignList=list(DesignTpl)
        
        # extract values
        DesignVal=[ getattr(self,dd) for dd in DesignList ]

                    
        return DesignVal
    
    
    def set_DesignVal(self,DesignValNew):
        ''' 
        Complementary method for get_DesignVal. In this case, the values of the
        Design variables can be reset given the values specified into DesignValNew
        
        The ii-the element of DesignValNew is related to the ii-th design variable 
        (being it float or array) returned by the self.list_deriv_vars() method.
        
        Warning:
        the method does not check if the design variables are still in the design 
        space.
        '''
        # get tuples of variables/functionals
        DesignTpl ,FunctionalsTpl = self.list_deriv_vars()
        DesignList=list(DesignTpl)
        
        for ii in range(len(DesignList)):
            setattr(self, DesignList[ii], DesignValNew[ii])    
        
        return self    
        
        
                
        
