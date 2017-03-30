'''@package PyBeam.Utils.XbeamLib
@brief      Functions found in lib_xbeam.f90, lib_rotvect.f90
@author     Rob Simpson, Salvatore Maraniello
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       15/01/2013
@pre        None
@warning    None
'''
import sys, os
sys.path.append( os.environ["SHARPYDIR"]+'/src' )
sys.path.append( os.environ["SHARPYDIR"]+'/src/Main' )

import numpy as np
import ctypes as ct
import SharPySettings as Settings
from scipy import linalg
from PyBeam.Utils.Misc import isodd
import BeamLib
import lib_rotvect




def QuatMult(quat0, quat1):
    w0, x0, y0, z0 = quat0
    w1, x1, y1, z1 = quat1
    qtot = np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                      x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                      x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0])
    return qtot



def RotCRV(Psi):
    '''
    Calculates the rotation matrix from CRV. 
    - Psi is the CRV that defines the rotation required to move a frame A over a 
      frame B.
    - C will be the rotation matrix from A to B or, alternatively, the projection 
      matrix from B to A


      @warning: for the same rotation, this function returns the transpose of
      Rot. i.e.
      RotCRV(Psi).T = Rot(psi2quat(Psi))
      @warning: equivalent to Psi2TransMat !!
    '''
    
    C = np.zeros((3,3))
    I = np.eye(3)
    PsiSkew=Skew(Psi)
    Psival = np.linalg.norm(Psi)
    
    if Psival>1e-6:
       C = I + np.sin(Psival)/Psival * PsiSkew \
             + (1.-np.cos(Psival))/Psival**2 * np.dot(PsiSkew,PsiSkew)
    else: 
        Cadd=I
        PsiPower = np.eye(3)
        tol=1e-8
        kk=0
        while np.max(np.abs(Cadd))>tol:
            # add current
            C+=Cadd
            # new term
            kk+=1
            PsiPower = np.dot(PsiPower,PsiSkew)
            Cadd = 1./np.math.factorial(kk) * PsiPower
        
    return C


def QuadSkew(Omega):
    """@brief Calculate incremental update operator for Quaternion based on Omega.
    Quaternion ODE to compute orientation of body-fixed frame a.
    Duplicated from lib_xbeam.f90
    See Shearer and Cesnik, 'Nonlinear Flight Dynamics of Very Flexible Aircraft',
    Journal of Aircraft 2007; 44 (5): 1528-1545
    """
    
    QS = np.zeros((4,4), ct.c_double, 'F')
    
    QS[0,1], QS[3,2], QS[1,0], QS[2,3] = \
        Omega[0], Omega[0], -Omega[0], -Omega[0]
        
    QS[0,2], QS[1,3], QS[2,0], QS[3,1] = \
        Omega[1], Omega[1], -Omega[1], -Omega[1]
        
    QS[0,3], QS[2,1], QS[1,2], QS[3,0] = \
        Omega[2], Omega[2], -Omega[2], -Omega[2]
        
    return QS


def Rot(q1):
    """@brief Calculate rotation matrix based on quaternions.
    See Aircraft Control and Simulation, pag. 31, by Stevens, Lewis.
    
    Remark: if A is a FoR obtained rotating a FoR G of angle fi about an axis n 
    (remind n will be invariant during the rotation), and q is the related 
    quaternion q(fi,n), the function will return the matrix Cag such that:
        - Cag rotates A to G
        - Cag transforms the coordinates of a vector defined in G component to A 
        components.
    """
    
    q = q1.copy(order='F')
    q = q/np.linalg.norm(q)
    
    RotMat = np.zeros((3,3), ct.c_double, 'F')
    
    RotMat[0,0] = q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
    RotMat[1,1] = q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2
    RotMat[2,2] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2
    
    RotMat[0,1] = 2.*(q[1]*q[2] + q[0]*q[3])
    RotMat[1,0] = 2.*(q[1]*q[2] - q[0]*q[3])
    
    RotMat[0,2] = 2.*(q[1]*q[3] - q[0]*q[2])
    RotMat[2,0] = 2.*(q[1]*q[3] + q[0]*q[2])
    
    RotMat[1,2] = 2.*(q[2]*q[3] + q[0]*q[1])
    RotMat[2,1] = 2.*(q[2]*q[3] - q[0]*q[1])
    
    return RotMat


def Euler2Quat(phi,theta,psi):
    """
    @brief Calculate quaternions based on Euler angles.
    See Aircraft Control and Simulation, pag. 31, by Stevens, Lewis.
    @var phi: roll
    @var theta: pitch
    @var psi: yaw
    Euler angles are
    @note: Euler angles are defined as a 321 sequence of rotations (i.e. about 
    the z (yaw), y'(pitch) and x'' (roll or bank) axis (see Stevens, Lewis, sec.1.4) 
    """
    q = np.zeros(4, ct.c_double, 'F')
    
    #Quaternion components
    q[0]=np.cos(phi/2)*np.cos(theta/2)*np.cos(psi/2)+np.sin(phi/2)*np.sin(theta/2)*np.sin(psi/2);
    q[1]=np.sin(phi/2)*np.cos(theta/2)*np.cos(psi/2)-np.cos(phi/2)*np.sin(theta/2)*np.sin(psi/2);
    q[2]=np.cos(phi/2)*np.sin(theta/2)*np.cos(psi/2)+np.sin(phi/2)*np.cos(theta/2)*np.sin(psi/2);
    q[3]=np.cos(phi/2)*np.cos(theta/2)*np.sin(psi/2)-np.sin(phi/2)*np.sin(theta/2)*np.cos(psi/2);

    #Normalize
    q = q/np.linalg.norm(q)
    
    return q


def Quat2Euler(qv):
    """
    @brief calculate euler angles from quaternions. Euler angles are
    defined as a 321 seuqence of rotations (i.e. about the z (yaw), 
    y'(pitch) and x'' (roll or bank) axis (see Stevens, Lewis, sec.1.4) 
    @var qv: quaternion vector
    @var ev: euler angles vector [roll, pitch, yaw]
    """
    
    ev=np.array([  
                np.arctan2( 2.0*(qv[0]*qv[1]+qv[2]*qv[3]) , 1.0-2.0*(qv[1]**2+qv[2]**2) ),
                np.arcsin ( 2.0*(qv[0]*qv[2]-qv[3]*qv[1])                               ), 
                np.arctan2( 2.0*(qv[0]*qv[3]+qv[1]*qv[2]) , 1.0-2.0*(qv[2]**2+qv[3]**2) )
                 ])
    
    return ev


def Psi2TransMat(Psi):
    """
    @brief Calculates the transformation matrix associated with CRV Psi.
    @details This gives the transformation from B to a, or the rotation from
              a to B.
    @warning: equivalent to RotCRV."""
    TransMat = np.zeros((3,3))
    Norm = TriadNorm(Psi)
    SkewPsi=Skew(Psi)
    # Check norm of Psi for linearised OR fully-nonlinear matrix.
    if Norm < Settings.RotEpsilon:
        TransMat[:,:] = np.identity(3) + SkewPsi + 0.5*np.dot(SkewPsi,SkewPsi) 
    else:
        TransMat[:,:] = np.identity(3) \
                      + np.sin(Norm)/Norm*SkewPsi \
                      + (1 - np.cos(Norm))/Norm**2 * np.dot(SkewPsi,SkewPsi)
    
    return TransMat


def TriadNorm(Triad):
    """@brief Returns the norm of a 3x1 tensor (Triad)."""
    return np.sqrt(Triad[0]**2+Triad[1]**2+Triad[2]**2)


def Skew(Vector):
    """@brief Returns the skew-symmetric Matrix associated to Vector"""
    
    SkewMat = np.zeros((3,3))
    
    SkewMat[0,1]=-Vector[2]
    SkewMat[0,2]= Vector[1]
    SkewMat[1,0]= Vector[2]
    SkewMat[1,2]=-Vector[0]
    SkewMat[2,0]=-Vector[1]
    SkewMat[2,1]= Vector[0]
    
    return SkewMat


def VectofSkew(Skew):
    """@brief Returns the vector defining skew-symmetric matrix Skew."""
    
    Vect = np.array([Skew[2,1],Skew[0,2],Skew[1,0]])
    
    return Vect


def Tangential(Psi):
    """@brief Calculates the tangential operator related to Psi, a cartesian
     rotation vector."""
    Tang = np.zeros((3,3))
    Norm = TriadNorm(Psi)
    SkewPsi=Skew(Psi)
    # Check norm of Psi for linearised OR fully-nonlinear matrix.
    if Norm < Settings.RotEpsilon:
        Tang[:,:] = np.identity(3) - 0.5* SkewPsi + 1.0/6.0*np.dot(SkewPsi,SkewPsi) 
    else:
        Tang[:,:] = np.identity(3) \
                    - (1-np.cos(Norm))/Norm**2 * SkewPsi \
                    + (1 - np.sin(Norm)/Norm) * np.dot(SkewPsi,SkewPsi)/Norm**2
    
    return Tang


def AddGravityLoads(BeamForces,XBINPUT,XBELEM,AELAOPTS,PsiDefor,
                      chord = 1.0, PsiA_G=None, FollowerForceRig=True):
    """@brief Apply nodal gravity loads. If the solver calling this function is
           using a FollowerForceRig=True option, loads are rotated so as the 
           guarantee that gravity are always considered as dead loads aligned
           along the z axis of the FoR G.
    @param BeamForces Nodal forces to update.
    @param XBINPUT Xbeam inputs.
    @param XBELEM Xbeam element information.
    @param AELAOPTS Aeroelastic options. If none, input will be neglected
    @param PsiA_G Cartesian rotation vector describing orientation of a-frame
           with respect to earth - rotation is performed only if 
           FollForceRig=True.
    @param chord Aerofoil chord, assumed to be 1.0 if ommited. 
    @details Offset of mass centroid from elastic axis is currently calculated
             using the AELAOPTS.ElasticAxis and .InertialAxis parameters which 
             are analogous to that used by Theodorsen. It therefore assumes the 
             aerofoil section is defined on the y-axis and hence the moment arm
             points in the B-frame y-direction.
    @warning Assumes evenly spaced distribution of nodes along beam.
    @warning Not valid for twisted aerofoil sections, i.e those that are defined
             with some theta angle.
    """
    
    MassPerLength = XBINPUT.BeamMass[0,0]
    
    if XBINPUT.NumNodesElem == 2:
        NodeSpacing = XBELEM.Length[0]
    elif XBINPUT.NumNodesElem == 3:
        NodeSpacing = 0.5*XBELEM.Length[0]
        
    ForcePerNode = -NodeSpacing * MassPerLength * XBINPUT.g 
    
    # Obtain transformation from Earth to a-frame.
    #if PsiA_G == None: CGa = Psi2TransMat(XBINPUT.PsiA_G)
    #else: CGa = Psi2TransMat(PsiA_G)
    if FollowerForceRig is True:
        # rotate gravity loads
        if PsiA_G == None: 
            PsiA_G = XBINPUT.PsiA_G
    else:
        # do not rotate gravity loads
        PsiA_G = np.zeros((3,))

    CGa = Psi2TransMat(PsiA_G)
    CaG = CGa.T
    
    # Force in a-frame.
    Force_a = np.dot(CaG,np.array([0.0,0.0,ForcePerNode]))
    
    # Indices for boundary nodes.
    Root = 0
    Tip = BeamForces.shape[0]-1
    
    # Apply forces.   CHANGE HERE
    BeamForces[Root+1:Tip,:3] += Force_a
    BeamForces[Root,:3] += 0.5*Force_a
    BeamForces[Tip,:3] += 0.5*Force_a
    
    # Add RB masses contribution
    for RBM in XBINPUT.RBMass:
        ForceRBM_a = -XBINPUT.g*RBM['Mass'][0,0]
        BeamForces[RBM['node'],:3]+=np.dot(CaG,np.array([0.0,0.0,ForceRBM_a]))

    # Get number of nodes per beam element.
    NumNodesElem = XBINPUT.NumNodesElem
    
    # Loop through nodes to get moment arm at each.
    for iNode in range(XBINPUT.NumNodesTot):
        
        # Work out what element we are in (works for 2 and 3-noded).
        if iNode == 0:
            iElem = 0
        elif iNode < XBINPUT.NumNodesTot-1:
            iElem = int(iNode/(NumNodesElem-1))
        elif iNode == XBINPUT.NumNodesTot-1:
            iElem = int((iNode-1)/(NumNodesElem-1))
            
        # Work out what sub-element node (iiElem) we are in.
        if NumNodesElem == 2:
            if iNode < XBINPUT.NumNodesTot-1:
                iiElem = 0
            elif iNode == XBINPUT.NumNodesTot-1:
                iiElem = 1
                
        elif NumNodesElem == 3:
            iiElem = 0
            if iNode == XBINPUT.NumNodesTot-1:
                iiElem = 2 
            elif isodd(iNode):
                iiElem = 1
        
        # Calculate transformation matrix for each node.
        CaB = Psi2TransMat(PsiDefor[iElem,iiElem,:])
        
        # Define moment arm in B-frame
        # Moment arm to elastic axis defined using EA and IA of Theodorsen.
        ###
        # The if/else allows using this method in structural only solvers.
        # TO IMPROVE: the inertial and elastic axis properties should be in 
        # Xbinput.
        if AELAOPTS != None:
            armY = -(AELAOPTS.InertialAxis - AELAOPTS.ElasticAxis)*chord/2.0
            armY_a = np.dot(CaB,np.array([0.0, armY, 0.0]))
        else:
            armY_a=np.zeros((3,))
        
        ### Calculate moment
        BeamForces[iNode,3:] += np.cross( armY_a, BeamForces[iNode,:3] )     
        #if (iNode == Root or iNode == Tip):
        #    BeamForces[iNode,3:] += np.cross(armY_a, 0.5*Force_a)
        #else:
        #    BeamForces[iNode,3:] += np.cross(armY_a, Force_a)
    
    

    
def LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                 PosIni, PsiIni, PosDefor, PsiDefor, \
                 Force_foll, Force_dead, CAG, iFactor):
    """@brief Assemble separate follower and dead loads.
    @param XBINPUT Xbeam inputs.
    @param XBELEM Xbeam element information.
    @param XBNODE Xbeam nodal information.
    @param XBOPTS Xbeam simulation options.
    @param NumDof numbers of DoF in the problem.
    @param PosInI nodal position at reference configuration.
    @param PsiIni nodal CRV at reference configuration.
    @param PosDefor nodal position at reference configuration.
    @param PsiDefor nodal CRV at reference configuration.
    @param Force_foll Matrix of applied nodal follower forces.
    @param Force_dead Matrix of applied nodal dead forces.
    @param CAG transformation matrix from inertial (G) to body-fixed frame (A).
    @param iFactor factor to account for optional load stepping.
    """
    
    # Initialise function variables
    NumNodes_tot=ct.c_int((XBINPUT.NumNodesElem - 1)*XBINPUT.NumElems + 1)
    
    KglobalFull_foll = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); ksf = ct.c_int()
                            
    FglobalFull_foll = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); fsf = ct.c_int()
                            
    FglobalFull_dead = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); fsd = ct.c_int()
                            
    Force_foll_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_dead_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    
    # Assemble separate force coefficient matrices
    BeamLib.Cbeam3_Asbly_Fglobal(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                        PosIni, PsiIni, PosDefor, PsiDefor,\
                        Force_foll*iFactor,NumDof,\
                        ksf, KglobalFull_foll, fsf, FglobalFull_foll,\
                        fsd, FglobalFull_dead, CAG)
    
    # Get follower and dead forces on constrained nodes
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                      ct.byref(ct.c_int(6)),\
                      Force_foll.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      ct.byref(NumDof),\
                      Force_foll_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                      ct.byref(ct.c_int(6)),\
                      Force_dead.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      ct.byref(NumDof),\
                      Force_dead_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    
    # Project follower and dead forces
    Qforces = np.dot(FglobalFull_foll, Force_foll_Dof*iFactor) + np.dot(FglobalFull_dead, Force_dead_Dof*iFactor)  

    #--------------------------------------
    # Repeat in case of rigid-body dynamics
    if XBOPTS.Solution.value == 912:
        FrigidFull_foll = np.zeros((6,NumDof.value+6), ct.c_double, 'F'); frf = ct.c_int()
        FrigidFull_dead = np.zeros((6,NumDof.value+6), ct.c_double, 'F'); frd = ct.c_int()
        
        Force_foll_All = np.zeros(NumDof.value+6, ct.c_double, 'F')
        Force_dead_All = np.zeros(NumDof.value+6, ct.c_double, 'F')
        
        # Assemble separate force coefficient matrices
        BeamLib.Xbeam_Asbly_Frigid(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                                     PosIni, PsiIni, PosDefor, PsiDefor,\
                                     NumDof, frf, FrigidFull_foll,\
                                     frd, FrigidFull_dead, CAG)
    
        # Get follower and dead forces on unconstrained nodes
        BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                                   ct.byref(ct.c_int(6)),\
                                   Force_foll.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                   ct.byref(ct.c_int(NumDof.value+6)),\
                                   Force_foll_All.ctypes.data_as(ct.POINTER(ct.c_double)) )
        BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                                   ct.byref(ct.c_int(6)),\
                                   Force_dead.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                   ct.byref(ct.c_int(NumDof.value+6)),\
                                   Force_dead_All.ctypes.data_as(ct.POINTER(ct.c_double)) )
        
        # Project follower and dead forces
        Qrigid = np.dot(FrigidFull_foll,Force_foll_All*iFactor) + np.dot(FrigidFull_dead,Force_dead_All*iFactor) 
    
    else:
        Qrigid = 0.0

    #------------------
    # Return everything
    return Qforces, KglobalFull_foll, Qrigid

def quat2psi(quat):
    """@brief Convert quaternion into Cartesian rotation vector.
    
    @param quat Quaternion.
    @returns Cartesian rotation vector.
    """
    return(lib_rotvect.lib_rotvect.rotvect_quat2psi(quat))

def psi2quat(psi):
    """@brief Convert Cartesian rotation vector into quaternion.
    
    @param psi Cartesian rotation vector
    @returns Quaternion.
    """    
    
    fi = np.linalg.norm(psi)
    if fi<1e-8: nv=np.zeros(3)
    else: nv = psi/fi 
    
    quat=np.zeros(4)
    quat[0] =np.cos(0.5*fi)
    quat[1:]=np.sin(0.5*fi)*nv   
    
    return quat
         

def a1(psi,b):
    """@brief Calculates the A1 operator: delta(T)*b=a1(psi,b)*Delta(psi))
    @param psi cartesian rotation vector.
    @param b vector in R^3.
    @returns A1 matrix (which is size R^3x3)
    @details If magnitude of psi is less than 1e-3 a linearized
    approximation is used to avoid numerical issues.
    """
    return(lib_rotvect.lib_rotvect.rotvect_a1(psi,b))

    
def LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
                 PosIni, PsiIni, PosDefor, PsiDefor, \
                 Force_foll, Force_dead, CAG, iFactor):
    """@brief Assemble separate follower and dead loads.
    @param XBINPUT Xbeam inputs.
    @param XBELEM Xbeam element information.
    @param XBNODE Xbeam nodal information.
    @param XBOPTS Xbeam simulation options.
    @param NumDof numbers of DoF in the problem.
    @param PosInI nodal position at reference configuration.
    @param PsiIni nodal CRV at reference configuration.
    @param PosDefor nodal position at reference configuration.
    @param PsiDefor nodal CRV at reference configuration.
    @param Force_foll Matrix of applied nodal follower forces.
    @param Force_dead Matrix of applied nodal dead forces.
    @param CAG transformation matrix from inertial (G) to body-fixed frame (A).
    @param iFactor factor to account for optional load stepping.
    """
    
    # Initialise function variables
    NumNodes_tot=ct.c_int((XBINPUT.NumNodesElem - 1)*XBINPUT.NumElems + 1)
    
    KglobalFull_foll = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); ksf = ct.c_int()
                            
    FglobalFull_foll = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); fsf = ct.c_int()
                            
    FglobalFull_dead = np.zeros((NumDof.value,NumDof.value),\
                            ct.c_double, 'F'); fsd = ct.c_int()
                            
    Force_foll_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_dead_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    
    # Assemble separate force coefficient matrices
    BeamLib.Cbeam3_Asbly_Fglobal(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                        PosIni, PsiIni, PosDefor, PsiDefor,\
                        Force_foll*iFactor,NumDof,\
                        ksf, KglobalFull_foll, fsf, FglobalFull_foll,\
                        fsd, FglobalFull_dead, CAG)
    
    # Get follower and dead forces on constrained nodes
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                      ct.byref(ct.c_int(6)),\
                      Force_foll.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      ct.byref(NumDof),\
                      Force_foll_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                      ct.byref(ct.c_int(6)),\
                      Force_dead.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      ct.byref(NumDof),\
                      Force_dead_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                      XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    
    # Project follower and dead forces
    Qforces = np.dot(FglobalFull_foll, Force_foll_Dof*iFactor) + np.dot(FglobalFull_dead, Force_dead_Dof*iFactor)  

    #--------------------------------------
    # Repeat in case of rigid-body dynamics
    if XBOPTS.Solution.value == 912:
        FrigidFull_foll = np.zeros((6,NumDof.value+6), ct.c_double, 'F'); frf = ct.c_int()
        FrigidFull_dead = np.zeros((6,NumDof.value+6), ct.c_double, 'F'); frd = ct.c_int()
        
        Force_foll_All = np.zeros(NumDof.value+6, ct.c_double, 'F')
        Force_dead_All = np.zeros(NumDof.value+6, ct.c_double, 'F')
        
        # Assemble separate force coefficient matrices
        BeamLib.Xbeam_Asbly_Frigid(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                                     PosIni, PsiIni, PosDefor, PsiDefor,\
                                     NumDof, frf, FrigidFull_foll,\
                                     frd, FrigidFull_dead, CAG)
    
        # Get follower and dead forces on unconstrained nodes
        BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                                   ct.byref(ct.c_int(6)),\
                                   Force_foll.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                   ct.byref(ct.c_int(NumDof.value+6)),\
                                   Force_foll_All.ctypes.data_as(ct.POINTER(ct.c_double)) )
        BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                                   ct.byref(ct.c_int(6)),\
                                   Force_dead.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                   ct.byref(ct.c_int(NumDof.value+6)),\
                                   Force_dead_All.ctypes.data_as(ct.POINTER(ct.c_double)) )
        
        # Project follower and dead forces
        Qrigid = np.dot(FrigidFull_foll,Force_foll_All*iFactor) + np.dot(FrigidFull_dead,Force_dead_All*iFactor) 
    
    else:
        Qrigid = 0.0

    #------------------
    # Return everything
    return Qforces, KglobalFull_foll, Qrigid

def quat2psi(quat):
    """@brief Convert quaternion into Cartesian rotation vector.
    
    @param quat Quaternion.
    @returns Cartesian rotation vector.
    """
    return(lib_rotvect.lib_rotvect.rotvect_quat2psi(quat))

def a1(psi,b):
    """@brief Calculates the A1 operator: delta(T)*b=a1(psi,b)*Delta(psi))
    @param psi cartesian rotation vector.
    @param b vector in R^3.
    @returns A1 matrix (which is size R^3x3)
    @details If magnitude of psi is less than 1e-3 a linearized
    approximation is used to avoid numerical issues.
    """
    return(lib_rotvect.lib_rotvect.rotvect_a1(psi,b))


if __name__ == '__main__':
    

    # Test quat/CRV conversion
    roll = 15.0 
    pitch= 56.0 
    yaw  = 24.0 
    
    qv = Euler2Quat(roll*np.pi/180.0, pitch*np.pi/180.0, yaw*np.pi/180.0)
    ev = 180.0/np.pi*Quat2Euler(qv)

    print('Original: (roll, pitch, yaw) = (%f,%f,%f)'%(roll, pitch, yaw))
    print('Final: (roll, pitch, yaw) = (%f,%f,%f)'%(ev[0], ev[1], ev[2]))
    

    # Test Cartesian rotation vector.
    
    #create matrices that correspond to Psi of nodes 1 and 2
    Psi1 = np.array([0.1, 0.1, 0.1])
    Psi2 = np.array([0.11, 0.12, 0.13])

    CaB1 = Psi2TransMat(Psi1)
    CaB2 = Psi2TransMat(Psi2)
    
    #Calc local rotation in frame B1
    CB1B2 = np.dot(CaB1.T,CaB2)
    
    #Get skew-symmetric matrix
    phitilda12 = linalg.logm(CB1B2)
    
    #get vector
    phi12 = VectofSkew(phitilda12)
    
    #get C matrix at midpoint
    CaBh = np.dot(CaB1,linalg.expm(0.5 * phitilda12))
    print('Objectively interpolated matrix = \n',CaBh,'\n')
    
    #now try by direct interp of CRV
    Psih_err = 0.5*(Psi2 + Psi1)
    CaBh_err = Psi2TransMat(Psih_err)
    print('Non-objective matrix = \n',CaBh_err,'\n')
    
    #diff
    print('diff = \n',CaBh_err - CaBh,'\n')
    
    #get objectively interped Psi
    Psih = VectofSkew(linalg.logm(CaBh))
    print('Objectively interpolated Psi = \n',Psih,'\n')
    
    #diff
    print('diff = \n',Psih_err - Psih,'\n')
    
    
    # Some test using tangential operator
    
    #Defien simple angular velocity based on CRV rate alone
    # Use a CRV to define rotation in sectional frame
    hingeRot = np.array([45.0*np.pi/180.0, 0.0, 0.0])
    hingeRot2 = np.array([45.0*np.pi/180.0, 0.0, 45.0*np.pi/180.0])
    hingeRotDot = np.array([45.0*np.pi/180.0, 0.0, 0.0])
    hingeRotDot2 = np.array([45.0*np.pi/180.0, 45.0*np.pi/180.0, 0.0])
    # Find lever arm from hinge point
    leverArm = np.array([0.0, 1.0, 0.0]) #Bnew frame
    CBBnew = Psi2TransMat(hingeRot)
    CBBnew2 = Psi2TransMat(hingeRot2)
    # Current implementation of B-frame velocity
    print('planar motion assumption implementation of B-frame velocity: ',np.dot(CBBnew,np.dot(Skew(hingeRotDot),leverArm)),'\n')
    print('General implementation of B-frame velocity: ',np.dot(CBBnew,np.dot(Skew(np.dot(Tangential(hingeRot),hingeRotDot)),leverArm)),'\n')
    
    # Current implementation of B-frame velocity
    print('planar motion assumption implementation of B-frame velocity: ',np.dot(CBBnew2,np.dot(Skew(hingeRotDot2),leverArm)),'\n')
    print('General implementation of B-frame velocity: ',np.dot(CBBnew2,np.dot(Skew(np.dot(Tangential(hingeRot2),hingeRotDot2)),leverArm)),'\n')
    

    Psival = np.pi*30./180.
    nv=np.array([1.,0.,0.])
    Psi = Psival*nv
    C=RotCRV(Psi)
