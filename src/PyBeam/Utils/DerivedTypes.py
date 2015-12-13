'''@package PyBeam.Utils.DerivedTypes
@brief      A collection of derived types that mirror those in 
xbeam_shared.f90 plus extras.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@author     A.Da-Ronch (Liverpool) 
@version    2.0
@date       14/11/2012
@pre        None
@warning    Stiffness parameter Sigma not implemented.

@modified   S.Maraniello
@date       08/09/2015
@comment    Sflag for spherical BCs added
            Methods to convert class to Xbopts
@note       Methods to set/remove ctypes specifications added to
            Xbopts and Xbinput. These are not required to use
            multiprocessing functionalities in python3.2
            
'''

import sys
import numpy as np #http://docs.scipy.org/doc/numpy/reference/
import ctypes as ct

class Xbopts:
    """@brief Simple struct-like class containing xbeam options from 
    xbeam_shared.f90.
    
    @param FollowerForce Follower force flag, bool.
    @param FollowerForceRigid Follower force in body-fixed frame flag, bool.
    @param PrintInfo Print information to stdout, bool.
    @param OutInBframe Print velocities in B-frame (if not, use a-frame)
    @param OutInaframe Print velocities in a-frame (if not, use inertial frame)
    @param ElemProj Element info computed in (1) global frame (2) fixed element
                    frame (3) moving element frame.
    @param MaxIterations Maximum number of iterations, int.
    @param NumLoadSteps Number of load increments, int.
    @param NumGauss Number of Gauss points in element integration
    @param Solution Solution process 
                = 102/112: cbeam3 linear/nonlinear static
                = 202/212: cbeam3 linear/nonlinear structural dynamic
                = 302/312: cbeam3 linear/nonlinear static + structural dynamic
                = 900/910:        linear/nonlinear rigid-body dynamic
                = 902/912: cbeam3 linear/nonlinear flexible-body dynamic
                =     922: cbeam3 nonlinear static + flexible-body dynamic
                =     952: cbeam3 linear flexible with nonlinear rigid-body dyn
    @param DeltaCurved Min. angle for two vectors to be parallel, double.
    @param MinDelta Convergence param for Newton-Rhaphson iterations, double.
    @param NewmarkDamp Numerical damping in Newmark integration scheme, double.
    @param EnforceAngVel_FoRA: if the n-th element of the list is True, the
            n-th component of the angular velocity of the FoR A (expressed in 
            FoR A components) is kept constant. 
    
    @warning 
        - If FollowerForce = ct.c_bool(True), beam forces must be applied
    in the local FoR, B.
        - init_from_class method unnecessary after adding set_ctypes, del_ctypes
        
    """

    def __init__(self, FollowerForce = ct.c_bool(False), \
                 FollowerForceRig = ct.c_bool(True), \
                 PrintInfo = ct.c_bool(False), OutInBframe = ct.c_bool(True), \
                 OutInaframe = ct.c_bool(False), ElemProj = ct.c_int(0), \
                 MaxIterations = ct.c_int(99), NumLoadSteps = ct.c_int(5), \
                 NumGauss = ct.c_int(1), Solution = ct.c_int(111), \
                 DeltaCurved = ct.c_double(1.0e-5), \
                 MinDelta = ct.c_double(1.0e-8), \
                 NewmarkDamp = ct.c_double(1.0e-4) ):
        """@brief Default initialisation is as defined in the original fortran
        derived type."""
        self.FollowerForce = FollowerForce
        self.FollowerForceRig = FollowerForceRig
        self.PrintInfo = PrintInfo
        self.OutInBframe = OutInBframe
        self.OutInaframe = OutInaframe
        self.ElemProj = ElemProj
        self.MaxIterations = MaxIterations
        self.NumLoadSteps = NumLoadSteps
        self.NumGauss = NumGauss
        self.Solution = Solution
        self.DeltaCurved = DeltaCurved
        self.MinDelta = MinDelta
        self.NewmarkDamp = NewmarkDamp
        
        self._ctypes_links=True  
        self._ctypes_attributes = ['FollowerForce', 'FollowerForceRig', 'PrintInfo',  
                 'OutInBframe', 'OutInaframe', 'ElemProj', 'MaxIterations', 'NumLoadSteps',
                 'NumGauss', 'Solution', 'DeltaCurved', 'MinDelta', 'NewmarkDamp' ]
        self._ctypes_conversion= [ [ ct.c_int, ct.c_bool, ct.c_double ], 
                                   [      int,      bool,       float ] ]
        
        
    def init_from_class(self,H):
        '''
        Method unnecessary!!!
        
        Given a class object H (containing all or part of the attributes 
        created in the  __init__ method but not defined as c_types), 
        overrides all the  values of  the attributes in self according to 
        the values of H attributes.
        '''
        
        #### working but not robust:
        #### If an attribute in misspelled, an error will not occur
        #
        ## extract all attributes
        attrdict=self.__dict__
        #for attrname in dict:
        #    ### get link to c type
        #    link=getattr(self,attrname)
        #    try: link.value=getattr(H,attrname)
        #    except AttributeError: print('%s will not be changed!!!')
        
        #### more robust:
        #
        # extract all attributes
        for attrname in self._ctypes_attributes:
            ### get link to c type
            link=getattr(self,attrname)
            if type(link) in self._ctypes_conversion[0]:
                link.value=getattr(H,attrname)
            else: pass # do nothing

        self._ctypes_links=True
        
        return self
                             

    def del_ctypes(self):
        '''
        converts all the ctypes attributes into normal
        python variables. The method assumes all attributes are
        ctypes
        '''

        if self._ctypes_links is False:
            raise NameError('ctypes already deleted!!!')   

        for attrname in self._ctypes_attributes:

            link = getattr(self,attrname)
            if type(link) in self._ctypes_conversion[0]:
                setattr(self,attrname,link.value)
            else: pass # nothing required

        self._ctypes_links=False
        
        return self        
        

    def set_ctypes(self):
        '''
        Convert all the attributes into ctypes. The conversion
        is according to self._ctypes_conversion.
        @warning: all attributes are converted to c_types!
        '''

        if self._ctypes_links is True:
            raise NameError('ctypes already set-up!!!')

        for attrname in self._ctypes_attributes:
            val = getattr(self,attrname)
            #print('Setting %s'%attrname)
            if type(val) in self._ctypes_conversion[1]:
                # find correct c_types function/type
                cc = self._ctypes_conversion[1].index( type(val) )
                cfun = self._ctypes_conversion[0][cc]
                # and overwrite the attribute
                setattr(self,attrname,cfun(val))

        self._ctypes_links=True
        
        return self                
        

    def print_all(self):
        for attr in self.__dict__:
            print("%s = %s" % (attr, getattr(self, attr)))
        
        
class Xbinput:
    """@brief Contains inputs for PyBeam functions.
    
    @param NumNodesElem Number of nodes in each beam element, either 2 or 3.
    @param NumElems Number of elements in the beam model.
    @param NumSteps Number of timesteps.
    @param BeamLength Length of beam.
    @param BeamStiffness 6-by-6 sectional stiffness matrix.
    @param BeamMass 6-by-6 sectional mass matrix.
    @param BConds 2-char string with boundary condition, 'CF' = Clamped-Free.
    @param Sigma Stiffness parameter. Not Implemented.
    @param iOut Output file.
    @param t0 Initial time.
    @param tfin Final time.
    @param dt Timestep.
    @param Omega Angular velocity of oscillatory motions.
    @param NumNodesTot total number of nodes in the model.
    @param ForceStatic NumNodes static force vector at nodes.
    @param ForceStatic_foll NumNodes static force vector of follower loads at nodes.
    @param ForceStatic_dead NumNodes static force vector of dead loads at nodes.
    @param ForceDyn Numnodes dynamic forces at nodes.
    @param ForceDyn_foll Numnodes dynamic forces of follower loads  at nodes.
    @param ForceDyn_dead Numnodes dynamic forces of dead loads at nodes.
    @param ForcingType Type of dynamic forcing.
    @param g acceleration due to gravity.
    @param PsiA_G CRV associated to the angle/axis that requited to rotate G over A
    
    @param str_damping_model: None/Prop apply no or proportional damping to structure
    @param str_damping_param: parameters for damping model. Default: alpha/beta for
            proportional damping
    @param sph_joint_damping if a spherical joint is in the model, this can be sued
           to add viscous damping to the rotational degrees of freedom of FoR A origin.
    @param EnforceAngVel_FoRA=3*[False]: enforce angular velocity in ii direction in FoR A
    @param EnforceTraVel_FoRA=3*[False]: enforce translational velocity in ii direction in FoR A
    @param 
    """

    def __init__(self, NumNodesElem, NumElems,
                 BeamLength = 1.0,
                 BeamStiffness = np.zeros((6,6),ct.c_double,'F'),
                 BeamMass = np.zeros((6,6),ct.c_double,'F'),
                 BConds = 'CF',
                 Sigma = 1.0,
                 iOut = 1,
                 t0 = 0.0,
                 tfin = 0.0,
                 dt = 0.0,
                 Omega = 0.0,
                 ForcingType = 'Const',
                 RampTime = 0.0,
                 g = 0.0, # Leave this alone!
                 PsiA_G = np.array([0.0,0.0,0.0]),
                 str_damping_model=None,
                 str_damping_param={'alpha': 0.0, 'beta':0.0},
                 sph_joint_damping=None
                 ):
        """@brief NumNodesElem and NumElems must be specified for initialisation
                  of force arrays.
        
        @param NumNodes Local variable storing number of nodes in model for
               array init.
        """
        self.NumNodesElem = NumNodesElem
        self.NumElems = NumElems
        self.BeamLength = BeamLength
        self.BeamStiffness = BeamStiffness
        self.BeamMass = BeamMass
        self.BConds = BConds
        self.Sigma = Sigma
        self.iOut = iOut
        self.t0 = t0
        self.tfin = tfin
        self.dt = dt
        self.Omega = Omega
        self.ForcingType = ForcingType
        self.RampTime = RampTime
        self.g = g
        self.PsiA_G = PsiA_G
        #self.quat0 = quat0
        
        self.EnforceAngVel_FoRA=3*[False]
        self.EnforceTraVel_FoRA=3*[False]
        
        # dihedral angle left/right
        self.dihedral_angle = [0.0, 0.0]
        # dihedral length in % total length
        self.dihedral_lambda = [0.0, 0.0]
        

        # Check number of nodes per element.
        if self.NumNodesElem != 2 and self.NumNodesElem != 3:
            sys.stderr.write('Invalid number of nodes per element\n')
        elif self.NumNodesElem == 2:
            NumNodesTot = NumElems + 1
        elif self.NumNodesElem == 3:
            NumNodesTot = 2*NumElems + 1
            
        self.NumNodesTot = NumNodesTot
        
        # Initialise nodal arrays.
        self.ForceStatic = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceStatic_foll = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceStatic_dead = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceDyn = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceDyn_foll = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceDyn_dead = np.zeros((NumNodesTot,6),ct.c_double,'F')
        
        self._ctypes_links=True  
        self._ctypes_attributes = [ 'BeamStiffness', 'BeamMass',
                'ForceStatic', 'ForceStatic_foll' 'ForceStatic_dead' 
                'ForceDyn', 'ForceDyn_foll', 'ForceDyn_dead' ]
        #self._ctypes_conversion= [ [ ct.c_int, ct.c_bool, ct.c_double ], 
        #                           [      int,      bool,       float ] ] 
        
        self.str_damping_model=str_damping_model
        self.str_damping_param=str_damping_param
        self.sph_joint_damping = sph_joint_damping


    def set_ctypes(self):
        '''
        Add ctypes specification to numpy.ndarray arrays
        '''

        if self._ctypes_links is True:
            raise NameError('ctypes already set-up!!!')

        for attrname in self._ctypes_attributes:
            
            val = getattr(self,attrname)
            assert( type(val) is np.ndarray)

            setattr(self,attrname, np.array( val, ct.c_double,'F' ) )

        self._ctypes_links=True
        
        return self  
    
    
    def del_ctypes(self):
        '''
        Remove ctypes specifications from numpy.ndarray
        '''

        if self._ctypes_links is False:
            raise NameError('ctypes already deleted!!!')

        for attrname in self._ctypes_attributes:
            
            val = getattr(self,attrname)
            assert( type(val) is np.ndarray)

            setattr(self,attrname, np.array( val ) )

        self._ctypes_links=False
        
        return self       
        
        
    def print_all(self):
        for attr in self.__dict__:
            print("%s = %s" % (attr, getattr(self, attr))) 


class Xbelem:
    """@brief Pythonic version of fortran arrays containing derived type 'Elem'
    data, one set of arrays for each element, as defined in xbeam_shared.f90."""
    def __init__(self, NumElems, MaxElNod = 3):
        """@brief Create 'NumElems' arrays with zero entries.
        
        @param NumNodes Number of nodes in each element.
        @param MemNo Member to which the element belongs.
        @param Conn Connectivities.
        @param Master Master node for each node j in each element:
        (j,m): Node belongs to master elem m (or 0 if current is master).
        (j,n): Node n within master element (or 0 if current is master).
        @param Length Undeformed length of each element.
        @param PreCurv Undeformed curvature of the element.
        @param Psi Undeformed rotation vector of element frame.
        @param Vector Element orientation vector - goes along the local Y axis.
        @param Mass Mass matrix (constant along the element).
        @param Stiff Stiffness matrix (constant along the element).
        @param InvStiff Inverse of stiffness matrix.
        @param RBMass Non-Structural (lumped) mass at element nodes.
        @details Memory mapping in line with f90:Elem(i)%Array(j,k,l) reference 
        requirement using existing f90:do_xbelem_var protocol in 
        Fx_Wrapper_PyFx.f90 by using the 
        np.Array((j*i,k,l),dtype=ct.*,order='F') syntax."""
        
        self.NumNodes = np.zeros(NumElems, dtype=ct.c_int, order='F')
        self.MemNo = np.zeros(NumElems, dtype=ct.c_int, order='F')
        self.Conn = np.zeros((MaxElNod*NumElems), dtype=ct.c_int, order='F')
        self.Master = np.zeros((MaxElNod*NumElems*2), dtype=ct.c_int, order='F')
        self.Length = np.zeros(NumElems, dtype=ct.c_double, order='F')
        self.PreCurv = np.zeros((3*NumElems), dtype=ct.c_double, order='F')
        self.Psi = np.zeros((3*NumElems), dtype=ct.c_double, order='F')
        self.Vector = np.zeros((3*NumElems), dtype=ct.c_double, order='f')
        self.Mass = np.zeros((6*NumElems,6), dtype=ct.c_double, order='F')
        self.Stiff = np.zeros((6*NumElems,6), dtype=ct.c_double, order='F')
        self.InvStiff = np.zeros((6*NumElems,6), dtype=ct.c_double, order='F')
        self.RBMass = np.zeros((MaxElNod*NumElems,6,6), dtype=ct.c_double, \
                                                               order='F')

        
class Xbnode:
    """@brief Pythonic nodal information as defined in xbeam_shared.f90
    
    @param Master Master node info for each node:
     (master elem, node within master elem).
    @param Vdof Degree of freedom in velocity vector.
    @param Fdof Degree of freedom in force vector.
    @param Sflag Spherical BC
    
    """
    def __init__(self, NumNodesTot):
        self.Master = np.zeros(2*NumNodesTot,dtype=ct.c_int,order='F')       
        self.Vdof = np.zeros(NumNodesTot,dtype=ct.c_int,order='F') 
        self.Fdof = np.zeros(NumNodesTot,dtype=ct.c_int,order='F')
        self.Sflag= np.zeros(NumNodesTot,dtype=ct.c_int,order='F')


class Xboutput:
    '''
    Class to store all structural related output of a simulation. Attributes can
    be preallocated.
    '''
    
    def __init__(self):
        self.QuatList=[]
        self.PsiList=[]
        self.ZetaList=[]
        self.ZetaStarList=[]
        self.ForceAeroList=[]
        self.GammaStarList=[]
        self.GammaList=[]
        self.CRVList=[]
        # for performance
        self.cputime=[]
     
    
def dump(obj):
    """@brief Prints all attributes of an object"""
    for attr in dir(obj):
        print("obj.%s = %s" % (attr, getattr(obj, attr)))



if(__name__ == '__main__'):
    print('Batch run of PyAsblyStruct.py...')
    
    print('Test: Xbopts class (see script)...')

    print('Default initialisation:')
    XBOPTS = Xbopts()
    print('FollowerForce = %r' %(XBOPTS.FollowerForce))

    print('Custom initialisation:')
    XBOPTS2 = Xbopts(False)
    print('FollowerForce = %r' %(XBOPTS2.FollowerForce))
    print('Solution = %r' %(XBOPTS.Solution))

    print('Default arguments:')
    print(Xbopts.__init__.__defaults__)
    print('FINISHED')
    
    print('Test: Xbinput class (see script)...')
    XBINPUT = Xbinput(3,8)
    print(XBINPUT.__dict__)
    print('FINISHED')
    
    print('Test Xbelem class memory')
    XBELEM=Xbelem(2,3)
    print(XBELEM.Mass)