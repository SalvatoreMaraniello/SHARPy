
import numpy as np
import ctypes as ct

from PyFSI.Beam2UVLM import CoincidentGrid
from PyFSI.Beam2UVLM import CoincidentGridForce
from PyAero.UVLM.Utils import UVLMLib

import Main.SharPySettings as Settings
import PyBeam.Utils.XbeamLib as xbl
import BeamLib



def assemble_GEBM_tensors(   iStep,
                             NumDof, NumNodes_tot,
                             XBELEM, XBNODE,
                             XBINPUT, VMINPUT, XBOPTS, 
                             AELAOPTS, VMOPTS,  
                             tmpVrel,     
                             PosIni, PsiIni,          
                             PosDefor, PsiDefor,  
                             PosDotDef, PsiDotDef,
                             PosDotDotDef, PsiDotDotDef,
                             Force, #?
                             MssFull, CssFull, KssFull, FstrucFull,
                             Msr, Csr,
                             Force_Dof, Qstruc,
                             MrsFull, CrsFull, KrsFull, FrigidFull,
                             Mrr, Crr, Qrigid,
                             Q, dQdt, dQddt, 
                             Force_All,
                             Msys, Csys, Ksys, Asys, Qsys ):
    
    # utilities:
    # (dummy variables to store unused output from assembly functions)
    ms, mr = ct.c_int(), ct.c_int()
    cs, cr = ct.c_int(), ct.c_int()
    ks, kr = ct.c_int(), ct.c_int()
    fs, fr = ct.c_int(), ct.c_int()
       
    #set tensors to zero 
    MssFull[:,:] = 0.0; CssFull[:,:] = 0.0
    KssFull[:,:] = 0.0; FstrucFull[:,:] = 0.0
    Msr[:,:] = 0.0; Csr[:,:] = 0.0
    Qstruc[:] = 0.0
    MrsFull[:,:] = 0.0; CrsFull[:,:] = 0.0
    KrsFull[:,:] = 0.0; FrigidFull[:,:] = 0.0
    Mrr[:,:] = 0.0; Crr[:,:] = 0.0
    Qrigid[:] = 0.0
    Msys[:,:],Csys[:,:],Ksys[:,:],Asys[:,:],Qsys[:] = 0.0,0.0,0.0,0.0,0.0
     
    #rigid-body velocities and orientation from state vector
    Quat = dQdt[NumDof.value+6:].copy('F')
    Quat = Quat/np.linalg.norm(Quat)
    Cao  = xbl.Rot(Quat)                       # input for xbeam assembly
    Cqr = np.zeros((4,6), ct.c_double, 'F')    # dummy: output for xbeam assembly
    Cqq = np.zeros((4,4), ct.c_double, 'F')    # dummy: output for xbeam assembly
    
    #Update matrices and loads for structural dynamic analysis
    tmpQuat=Quat.copy('F')
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,     # in
                         PosIni, PsiIni, PosDefor, PsiDefor,                # in
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,  # in
                         Force, tmpVrel, 0*tmpVrel,                         # in
                         NumDof, Settings.DimMat,                           # out
                         ms, MssFull, Msr,                                  # out
                         cs, CssFull, Csr,                                  # out
                         ks, KssFull, fs, FstrucFull,Qstruc,                # out
                         XBOPTS, Cao)                                       # in
    
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),
                      ct.byref(ct.c_int(6)),
                      Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                      ct.byref(NumDof),
                      Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),
                      XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
            
    #Update matrices for rigid-body dynamic analysis
    BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,      # in
                         PosIni, PsiIni, PosDefor, PsiDefor,                # in
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,  # in
                         tmpVrel, 0*tmpVrel, tmpQuat,                       # in
                         NumDof, Settings.DimMat,                           # out
                         mr, MrsFull, Mrr,                                  # out
                         cr, CrsFull, Crr, Cqr, Cqq,                        # out
                         kr, KrsFull, fr, FrigidFull, Qrigid,               # out
                         XBOPTS, Cao)                                       # in

    BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),
                               ct.byref(ct.c_int(6)),
                               Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                               ct.byref(ct.c_int(NumDof.value+6)),
                               Force_All.ctypes.data_as(ct.POINTER(ct.c_double)) )

    #Assemble discrete system matrices with linearised quaternion equations  
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4): Unit4[i,i] = 1.0     
    Msys[:NumDof.value,:NumDof.value] = MssFull.copy('F')
    Msys[:NumDof.value,NumDof.value:NumDof.value+6] = Msr.copy('F')
    Msys[NumDof.value:NumDof.value+6,:NumDof.value] = MrsFull.copy('F')
    Msys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Mrr.copy('F')
    Msys[NumDof.value+6:,NumDof.value+6:] = Unit4.copy('F')
    
    Csys[:NumDof.value,:NumDof.value] = CssFull.copy('F')
    Csys[:NumDof.value,NumDof.value:NumDof.value+6] = Csr.copy('F')
    Csys[NumDof.value:NumDof.value+6,:NumDof.value] = CrsFull.copy('F')
    Csys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Crr.copy('F')
    
    Csys[NumDof.value+6:,NumDof.value:NumDof.value+6] = Cqr.copy('F')
    Csys[NumDof.value+6:,NumDof.value+6:] = Cqq.copy('F')
    
    Ksys[:NumDof.value,:NumDof.value] = KssFull.copy('F')
    Ksys[NumDof.value:NumDof.value+6,:NumDof.value] = KrsFull.copy('F')
             
    #Separate assembly of follower and dead loads   
    #tmpForceTime=ForceTime[iStep+1].copy('F') 
    tmpQforces,Dummy,tmpQrigid = xbl.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof,
                                     PosIni, PsiIni, PosDefor, PsiDefor,
                                    (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll[iStep,:,:] ),
                                    (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead[iStep,:,:] ),
                                    Cao,1)  

    Qstruc -= tmpQforces      
    Qrigid -= tmpQrigid

    #Compute residual to solve update vector
    Qstruc -= np.dot(FstrucFull, Force_Dof)
    Qrigid -= np.dot(FrigidFull, Force_All)
    
    Qsys[:NumDof.value] = Qstruc
    Qsys[NumDof.value:NumDof.value+6] = Qrigid
    Qsys[NumDof.value+6:] = np.dot(Cqq,dQdt[NumDof.value+6:])

    Qsys += np.dot(Msys,dQddt)
    
    # include damping
    if XBINPUT.str_damping_model == 'prop':
        Cdamp = XBINPUT.str_damping_param['alpha'] * MssFull + \
                XBINPUT.str_damping_param['beta']  * KssFull
        Csys[:NumDof.value,:NumDof.value] += Cdamp
        Qsys[:NumDof.value] += np.dot(Cdamp, dQdt[:NumDof.value])


    #---------------------------------------------------------------------------------- special BCs
    
    iiblock=[]
    if XBNODE.Sflag.any():
        # block translations (redundant:)
        for ii in range(3):
            if XBINPUT.EnforceTraVel_FoRA[ii] == True: iiblock.append(NumDof.value+ii)
        # block rotations
        iirotfree=[] # free rotational dof 
        for ii in range(3):
            if XBINPUT.EnforceAngVel_FoRA[ii] is True: iiblock.append(NumDof.value+3+ii)
            else: iirotfree.append(NumDof.value+3+ii)
    
    if len(iiblock)>0:
        Msys[iiblock,:] = 0.0
        Msys[:,iiblock] = 0.0 
        Msys[iiblock,iiblock] = 1.0
        Csys[iiblock,:] = 0.0
        Ksys[iiblock,:] = 0.0
        Qsys[iiblock]   = 0.0
        if XBINPUT.sph_joint_damping is not None:
            Csys[iirotfree,iirotfree] += XBINPUT.sph_joint_damping
            Qsys[iirotfree] += XBINPUT.sph_joint_damping*dQdt[iirotfree]

    return Msys, Csys, Ksys, Qsys


def init_GEBM_tensors(NumDof):
    '''
    Initialise structural variables
    '''
    
    
    MssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    CssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    KssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    FstrucFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    
    Msr = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Csr = np.zeros((NumDof.value,6), ct.c_double, 'F')
    
    X     = np.zeros(NumDof.value, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    Qstruc = np.zeros(NumDof.value, ct.c_double, 'F')
    
    #Initialise rigid-body system tensors
    MrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F')
    CrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F') 
    KrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F') 
    FrigidFull = np.zeros((6,NumDof.value+6), ct.c_double, 'F')
    
    Mrr = np.zeros((6,6), ct.c_double, 'F')
    Crr = np.zeros((6,6), ct.c_double, 'F')   
    Qrigid = np.zeros(6, ct.c_double, 'F')
    
    #Initialise full system tensors
    Q     = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    DQ    = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    dQdt  = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    dQddt = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    Force_All = np.zeros(NumDof.value+6, ct.c_double, 'F')

    Msys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F')
    Csys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F') 
    Ksys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F') 
    Asys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F')
    
    Qsys = np.zeros(NumDof.value+6+4, ct.c_double, 'F')

    return     ( MssFull, CssFull, KssFull, FstrucFull, Qstruc,
      MrsFull, CrsFull, KrsFull, FrigidFull, Qrigid,
      Mrr, Crr, Msr, Csr, 
      Force_Dof, Force_All,
      X, dXdt, Q, DQ, dQdt, dQddt,
      Msys, Csys, Ksys, Asys, Qsys )


def update_UVLMsol(  iStep, dt, Time, Force,
                     XBINPUT, VMINPUT, XBELEM, AELAOPTS, VMOPTS, 
                     PsiA_G, currVrel, OriginA_a, Section,
                     PosDefor, PsiDefor, PosDotDef, PsiDotDef, 
                     Ufree, AIC, BIC,
                     AeroForces, Zeta, ZetaDot, ZetaStar, Gamma, GammaStar ): 
    
    '''
    Wraps operations required to update the aerodynamic force.
    
    
    @warning: while the aim of the routine is to provide the nodal values of the 
              aerodynamic forces (Force array) the wake information (Zeta*, Gamma* 
              variables) have to be returned to the main solver so as to allow:
                  a. to update at following time-steps (code crash otherwise)
                  b. to saving wake information
    '''
    
    # zero aero forces.
    Force.fill(0.0)
    AeroForces[:,:,:] = 0.0
    
    # Update origin.
    # sm: origin position projected in FoR A

    OriginA_a[:] = OriginA_a[:] + currVrel[:3]*dt #sm: OriginA_a initialised to zero
    
    # Update control surface deflection.
    if VMINPUT.ctrlSurf != None:
        # open-loop control
        for cc in range(len(VMINPUT.ctrlSurf)):
            VMINPUT.ctrlSurf[cc].update(Time[iStep],iStep=iStep)
    
    # Generate surface grid.
    CoincidentGrid(PosDefor, PsiDefor, Section, currVrel[:3], 
                   currVrel[3:], PosDotDef, PsiDotDef, XBINPUT,
                   Zeta, ZetaDot, OriginA_a, PsiA_G,
                   VMINPUT.ctrlSurf)
    
    # Update wake geom       
    #'roll' data.
    ZetaStar = np.roll(ZetaStar,1,axis = 0)
    GammaStar = np.roll(GammaStar,1,axis = 0)
    #overwrite grid points with TE.
    ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
    # overwrite Gamma with TE value from previous timestep.
    GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
    
    # Apply gust velocity.
    if VMINPUT.gust != None:   Utot = Ufree + VMINPUT.gust.Vels(Zeta)
    else:                      Utot = Ufree
    
    # Solve for AeroForces
    UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Utot, ZetaStar, VMOPTS, 
                   AeroForces, Gamma, GammaStar, AIC, BIC)
    
    # Apply density scaling
    AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
    
    #if Settings.PlotTec==True:
    #    FileObject=write_TecPlot(Zeta, ZetaStar, Gamma, GammaStar, NumSteps.value, iStep, 
    #                             Time[iStep], SaveDict,FileObject=FileObject)

    # map AeroForces to beam.
    CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                        Force, PsiA_G)

    return Force, Zeta, ZetaStar, Gamma, GammaStar   


def set_res_tight(ResLog0,ResLogMin=1e1,ResLogMax=1e3,ResLogLim=1e-1):
    '''
    Returns an adaptive form of the residual R such that, if
        ResLog < R
    the aeroelastic wake is frozen.
    
    ResLogMin: limit value below which the wake is frozen.
    ResLogMax: value above which the residual R is kept constant to ResLogLim
    '''

    LogMax=np.log10(ResLogMax)
    LogMin=np.log10(ResLogMin)
    Log0 = np.log10(ResLog0)
    LogLim=np.log10(ResLogLim)
    
    if Log0 <=LogMin: 
        R=ResLogMin
    elif Log0 > LogMax: 
        R= np.power(10,LogLim)
    else:
        L = LogMin + (Log0 - LogMin)/(LogMax - LogMin)*(LogLim - LogMin)
        R = np.power(10.0,L)
    
    return R
    

def panellingFromFreq(freq,c=1.0,Umag=1.0):
    """@brief Calculate adequate spatial/temporal resolution on UVLM grid
    based on a frequency of interest.
    @param freq Frequency of interest.
    @param c chord of wing.
    @param Umag mean relative free-stream velocity magnitude.
    @returns M Number of chordwise panels for wing.
    @returns DelTime Suggested timestep [seconds.]
    """
    k = freq*c/(2*Umag) #get expected reduced freq
    M = int(50.0*k/np.pi) #get adequate panelling
    DelTime = c/(Umag*M) #get resulting DelTime
    return M, DelTime


def lagsolver(M,Q,MinTol=1e-3, MaxIter=1):
    '''
    Lagged solution of couple rigid-flexible body dynamics
    
    solves M x = Q
    
    where M = [ [Mstr , Msr],
                [ Mrs, Mrig ] ]
    Q = [Qstr, Qrig]
    
    @warning: no criterion for checking tolerance implemented: use for
    single iteration only!
    '''
    
    assert MaxIter==1, 'No criterion for checking tolerance implemented: use for single '\
                       'iteration only!'

    NumTot = M.shape[0]
    NumDof = NumTot-10
    
    iistr=[ii for ii in range(NumDof)]
    iirig=[ii for ii in range(NumDof,NumTot)]
    
    Mstr=M[:NumDof,:NumDof]
    Mrig=M[NumDof:,NumDof:]
    
    x=np.zeros((NumTot,))
    xstr = np.zeros((NumDof,))
    xrig = np.zeros((10,))
    
    tol=1
    Iter=0
    
    while tol>MinTol and Iter<MaxIter:
        # compute inertia due to structural vibrations
        qinertia = np.dot(M[iirig,:NumDof],xstr)
        # solve rigid-body dynamics
        xrig=np.linalg.solve(Mrig, Q[iirig] - qinertia)
        # compute inertia due rigid-body dynamics
        qinertia = np.dot(M[iistr,NumDof:], xrig)
        # solve for structure
        xstr = np.linalg.solve(Mstr, Q[iistr] - qinertia)
        # update
        Iter += 1
    
    x[:NumDof]=xstr
    x[NumDof:]=xrig
    
    return x


