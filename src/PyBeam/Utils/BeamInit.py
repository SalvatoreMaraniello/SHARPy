'''@package PyBeam.Utils.BeamInit
@brief      For initialising variables required for beam simulations.
@author     Rob Simpson, Salvatore Maraniello
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       10/12/2012 -> Jan 2016
@pre        None
@warning    ForceTime term in Dynamic method obsolete. To be removed
'''

import sys
import DerivedTypes
import SharPySettings as Settings
import numpy as np
import ctypes as ct
import BeamIO
import BeamLib
import PyLibs.CVP.fourier, PyLibs.CVP.spline
from IPython import embed


def Static(XBINPUT,XBOPTS, moduleName = None):
    """@brief Initialise everything needed for Static beam simulation."""

    # Declare variables that are not dependent on NumNodes_tot.
    XBELEM      = DerivedTypes.Xbelem(XBINPUT.NumElems,Settings.MaxElNod)
    NumDof      = ct.c_int()
    PsiIni      = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                           dtype=ct.c_double, order='F')
        
    # Load customised inputs
    if moduleName != None:
        Input = __import__(moduleName)
    else:
        Input = __import__('Input')
        
    # Check inputs.
    XBINPUT, XBOPTS = Input.Setup(XBINPUT, XBOPTS)
        
    # Set-up element properties
    NumNodes_tot, XBELEM = Input.Elem(XBINPUT, XBOPTS, XBELEM)
    # Set-up nodal properties.
    PosIni, PhiNodes, BoundConds = \
        Input.Node(XBINPUT, XBOPTS, NumNodes_tot, XBELEM)
    # Compute initial (undeformed) geometry.
    BeamLib.f_xbeam_undef_geom( 
        ct.byref(ct.c_int(XBINPUT.NumElems)), 
        ct.byref(NumNodes_tot), 
        XBELEM.NumNodes.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.MemNo.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.Conn.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.Master.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.Length.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.PreCurv.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Psi.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Vector.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Mass.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Stiff.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.InvStiff.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.RBMass.ctypes.data_as(ct.POINTER(ct.c_double)),
        PosIni.ctypes.data_as(ct.POINTER(ct.c_double)), 
        PhiNodes.ctypes.data_as(ct.POINTER(ct.c_double)), 
        PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.byref(XBOPTS.FollowerForce),
        ct.byref(XBOPTS.FollowerForceRig),
        ct.byref(XBOPTS.PrintInfo),
        ct.byref(XBOPTS.OutInBframe),
        ct.byref(XBOPTS.OutInaframe),
        ct.byref(XBOPTS.ElemProj),
        ct.byref(XBOPTS.MaxIterations),
        ct.byref(XBOPTS.NumLoadSteps),
        ct.byref(XBOPTS.NumGauss),
        ct.byref(XBOPTS.Solution),
        ct.byref(XBOPTS.DeltaCurved),
        ct.byref(XBOPTS.MinDelta),
        ct.byref(XBOPTS.NewmarkDamp) )
    # Write to undeformed geometry to file.
    WriteMode = 'w'
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Writing file (mode: %s) %s' %(WriteMode,
                         Settings.OutputDir + Settings.OutputFileRoot + 
                         '_SOL' + XBOPTS.Solution.value + '_und.dat\n'))
    
    BeamIO.WriteUndefGeometry(XBINPUT.NumElems,NumNodes_tot.value,XBELEM,
                              PosIni,PsiIni,
                              Settings.OutputDir + Settings.OutputFileRoot +
                              '_INIT' , WriteMode)
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Identify nodal degrees of freedom ... ')
    # Initialize nodal data derived type.
    XBNODE = DerivedTypes.Xbnode(NumNodes_tot.value)

    BeamLib.f_xbeam_undef_dofs( 
        ct.byref(ct.c_int(XBINPUT.NumElems)), 
        ct.byref(NumNodes_tot), 
        XBELEM.NumNodes.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.MemNo.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.Conn.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.Master.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBELEM.Length.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.PreCurv.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Psi.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Vector.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Mass.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.Stiff.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.InvStiff.ctypes.data_as(ct.POINTER(ct.c_double)), 
        XBELEM.RBMass.ctypes.data_as(ct.POINTER(ct.c_double)), 
        BoundConds.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)), 
        XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)), 
        ct.byref(NumDof),
        XBNODE.Sflag.ctypes.data_as(ct.POINTER(ct.c_int)) )

    if XBOPTS.PrintInfo==True:
        sys.stdout.write('done\n')
    
    # Add RBMasses
    # - this is added only to the master nodes (XBNODE.Master), otherwise
    # is double accounted for.
    Master = XBNODE.Master.reshape((len(XBNODE.Master)//2,2))
    for RBM in XBINPUT.RBMass:
        # Find elements associated to the node
        #  - node numbering starts from 1 in XBELEM.Conn
        #  - numbering of ee and ll starts from 1
        #  - len(ccvec)>1 if the node is at the edge of more then 1 element
        ccvec = np.where(XBELEM.Conn == RBM['node']+1)[0]
        for cc in ccvec:
            # Check if node is a master
            ee = cc//3 + 1  # element to which node belongs
            ll =  cc%3 + 1  # local position in element
            if ll in Master[Master[:,0]==ee,1]:
                XBELEM.RBMass[cc,:,:] = RBM['Mass']

    return XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            XBNODE, NumDof
            

def Dynamic(XBINPUT,XBOPTS, moduleName = None):
    """@brief Initialise everything for dynamic analysis.
    
    @param moduleName Full path to module containing custom PyBeam setup
            scripts. E.g "/home/rjs10/git/SHARPy/output/DancingFrame/PyBeamInput
            "
            """

    # Load customised inputs
    if moduleName != None:
        Input = __import__(moduleName)
        #Set up dynamic testcase.
        try: 
            XBINPUT, XBOPTS = Input.DynSetup(XBINPUT, XBOPTS)
        except AttributeError: 
            Input = __import__('Input')
        else:
            raise
    else:
        Input = __import__('Input')

    #Create time vector.
    Time = np.arange(XBINPUT.t0,XBINPUT.tfin  + XBINPUT.dt, XBINPUT.dt,
                     ct.c_double)
    
    #Calculate number of timesteps.
    NumSteps = ct.c_int(len(Time) - 1)
    
    # sm: set dynamic loads
    XBINPUT = DynLoads(XBINPUT)
    
    # -------------------------------------------------- sm: move into DynLoads

    #Create force-amp-in-time array.
    ForceTime = np.zeros(NumSteps.value+1,ct.c_double,'F')
    if XBINPUT.ForcingType == 'Const':
        ForceTime[:] = 1.0
        
    elif XBINPUT.ForcingType == 'Sin':
        assert XBINPUT.Omega != 0.0, ('Dynamic forcing requested with zero ' +
                                       'frequency')
        for i in np.arange(len(Time)):
            ForceTime[i] = np.sin(XBINPUT.Omega * Time[i])
    elif XBINPUT.ForcingType == 'Ramp':
        assert XBINPUT.RampTime <= XBINPUT.tfin,'Ramp time greater than final time'
        assert XBINPUT.RampTime >= XBINPUT.t0,'Ramp time less than start time'
        assert Time.__contains__(XBINPUT.RampTime),('RampTime is not equal to '+
                                                    'element in Time array')
        ForceTime[:] = 1.0
        Gradient = 1.0/(XBINPUT.RampTime - XBINPUT.t0)
        t = XBINPUT.t0
        i = 0
        while t <= XBINPUT.RampTime:
            ForceTime[i] = Gradient*i*XBINPUT.dt
            i += 1
            t += XBINPUT.dt
    
    elif XBINPUT.ForcingType == 'RampSin':
        assert XBINPUT.Omega != 0.0, ('Dynamic forcing requested with zero ' +
                                       'frequency')
        assert XBINPUT.RampTime <= XBINPUT.tfin,'Ramp time greater than final time'
        assert XBINPUT.RampTime >= XBINPUT.t0,'Ramp time less than start time'
        assert Time.__contains__(XBINPUT.RampTime),('RampTime is not equal to '+
                                                    'element in Time array')
        
        #Define sinusoid.
        SinVect = np.zeros_like(ForceTime, ct.c_double, 'F')
        for i in np.arange(len(Time)):
            SinVect[i] = np.sin(XBINPUT.Omega * Time[i])
            
        #Define ramp.
        RampVect = np.zeros_like(ForceTime, ct.c_double, 'F')
        RampVect[:] = 1.0
        Gradient = 1.0/(XBINPUT.RampTime - XBINPUT.t0)
        t = XBINPUT.t0
        i = 0
        while t <= XBINPUT.RampTime:
            RampVect[i] = Gradient*i*XBINPUT.dt
            i += 1
            t += XBINPUT.dt
            
        #Combine.
        for i in np.arange(len(ForceTime)):
            ForceTime[i] = RampVect[i] * SinVect[i]
            
    elif XBINPUT.ForcingType == 'InitialLoad':
        ForceTime[0] = 1.0
        
    elif  XBINPUT.ForcingType == '1-cos':
        GustTime = 0.01
        for iStep in range(len(ForceTime)):
            if (Time[iStep] > 0.0 and Time[iStep] < GustTime):
                ForceTime[iStep] = 1.0 - np.cos(2*np.pi*Time[iStep]/GustTime)/2.0
            #END if
        #END iStep
            
    else:
        raise ValueError('ForcingType not recognised')

    # -------------------------------------------------- sm: move into DynLoads

    
    # Create default forced velocities and accelerations if not specified.
    if hasattr(XBINPUT, 'ForcedVel'):
        assert XBINPUT.ForcedVel.shape[0] == NumSteps.value+1
        ForcedVel = XBINPUT.ForcedVel.copy('F')
        ForcedVelDot = XBINPUT.ForcedVelDot.copy('F')
    else:
        ForcedVel = np.zeros((NumSteps.value+1,6),ct.c_double,'F')
        ForcedVelDot = np.zeros((NumSteps.value+1,6),ct.c_double,'F')

    # Rates of Dof Arrays.
    PosDotDef  = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    PsiDotDef  = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                           ct.c_double, 'F')
    
    # Create Outgrids and check only 1 node specified.
    OutGrids = np.zeros(XBINPUT.NumNodesTot, ct.c_bool, 'F')
    OutGrids[-1] = True #at tip
    Nonzeros = np.nonzero(OutGrids)
    if len(Nonzeros[0]) == 0:
        raise Exception('Zero Outgrid output nodes: min is 1.')
    elif len(Nonzeros[0]) > 1:
        raise Exception('Too many Outgrid output nodes: max is 1.')
    
    # Position/rotation history at Outgrids = True element.
    PosPsiTime = np.zeros((NumSteps.value+1,6), ct.c_double, 'F')
    VelocTime  = np.zeros((NumSteps.value+1,XBINPUT.NumNodesTot),
                           ct.c_double, 'F')
    
    # Position of all nodes at all times.
    DynOut     = np.zeros(((NumSteps.value+1)*XBINPUT.NumNodesTot,3),
                           ct.c_double, 'F')
    
    return Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
            PosDotDef, PsiDotDef,\
            OutGrids, PosPsiTime, VelocTime, DynOut
    
    
def DynLoads(XBINPUT):
    '''
    Compute forcing terms for beam model based on a collection of methods
    
    @todo: move here the loads methods inside Dynamic function
    '''
    
    # pre-allocate...
    Time = np.arange(XBINPUT.t0,XBINPUT.tfin  + XBINPUT.dt, XBINPUT.dt,
                     ct.c_double)
    shape = (len(Time),XBINPUT.NumNodesTot,6)
    XBINPUT.ForceDyn_foll = np.zeros(shape,ct.c_double,'F')
    XBINPUT.ForceDyn_dead = np.zeros(shape,ct.c_double,'F')


    if XBINPUT.ForceDynType=='hard-coded':
        assert XBINPUT.ForceDyn.shape==shape, 'If using XBINPUT.ForceDynType=%s'\
            ' XBINPUT.ForceDyn needs to have shape ('+ 3*'%d, '%shape + ') !'

    else:
        # allocate
        XBINPUT.ForceDyn = np.zeros(shape,ct.c_double,'F')

        # ... and fill
        if XBINPUT.ForceDynType=='dss': 
            print('NumNodesTot is %f'%XBINPUT.NumNodesTot)
            params = XBINPUT.ForceDynDict['dss']
            #assert params['Node']<XBINPUT.NumNodesTot-1 \
            #("params[Node] can't be higher then NumNodes-1!")
            for comp in range(6):
                XBINPUT.ForceDyn[:,params['Node'],comp] = \
                    PyLibs.CVP.fourier.sin_series_vec( Time,
                                    params['Acf'][:,comp],params['Fcf'][:,comp])

        elif XBINPUT.ForceDynType=='spline': 
            params = XBINPUT.ForceDynDict['spline']
            XBINPUT.ForceDyn = PyLibs.CVP.spline.nodal_force(
                XBINPUT.NumNodesTot, Time,params['Node'],params['TCint'], 
                                                  params['Scf'],params['order'])   

        elif XBINPUT.ForceDynType=='modal': 
            params = XBINPUT.ForceDynDict['modal']
            NT=shape[0]
            assert NT==len(params['ForceTime']), 'ForceTime parameter needs to '\
                                  'have the same length as the number of time-steps'
            for tt in range(shape[0]):
                XBINPUT.ForceDyn[tt,:,:]=params['ForceTime'][tt]*params['ForceNodes']

        elif XBINPUT.ForceDynType is None:
            pass
        else:
            raise NameError('XBINPUT.ForceDynType method not recognised!')




    return XBINPUT




if __name__ == '__main__':
    XBINPUT = DerivedTypes.Xbinput(2,1)
    XBINPUT.tfin = 1.0
    XBINPUT.dt = 0.01
    
    XBOPTS = DerivedTypes.Xbopts()
    
    XBINPUT.BeamStiffness[0,0] = 1.0e+09 #otherwise inverse fails
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 4.0e+06
    
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = Static(XBINPUT,XBOPTS)
    
    XBINPUT.Omega = 4*np.pi
    XBINPUT.ForcingType = '1-cos'
    XBINPUT.RampTime = 0.5
                
    Dynamic(XBINPUT,XBOPTS)        
