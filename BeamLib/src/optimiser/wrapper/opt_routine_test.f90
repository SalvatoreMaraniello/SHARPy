!-> Program.- MAIN - 25july2014 Updated: july 2014
!
!-> Author.- Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!   test run for opt_routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program opt_routine_test

 use opt_routine
 use opt_cost
 use xbeam_shared , only: MaxElNod
 use lib_isosec
 use input, only: input_dynforce, input_forcedvel


 implicit none


! ------------------------------------------------------------- Define Variables
    !real(8), pointer, dimension(:) :: pCOST

! Design Variables:
! These are normally inported from the input module through input_setup.
    real(8)  :: BeamLength1, BeamLength2
    real(8)  :: BeamStiffness(6,6)
    real(8)  :: BeamMass(6,6)
    real(8)  :: ExtForce(3)
    real(8)  :: ExtMomnt(3)
    real(8)  :: SectWidth
    real(8)  :: SectHeight

    real(8)  :: ThetaRoot
    real(8)  :: ThetaTip
    real(8)  :: TipMass
    real(8)  :: TipMassY
    real(8)  :: TipMassZ
    real(8)  :: Omega


! Static Loading Input:
     real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.


! Problem Setup:
    integer :: NumElems     ! Number of elements
    integer :: NumNodes     ! Number of nodes in the model.
    real(8), allocatable :: BeamSpanStiffness(:,:,:) ! Local Stiffness Matrix (opt)
    real(8), allocatable :: BeamSpanMass(:,:,:)      ! Local Mass Matrix (opt)
    real(8), allocatable :: PhiNodes(:)         ! Local Twist angle (opt)


 ! Problem Setup Shared:
    integer           :: NumNodesElem   ! Number of nodes on each element.
    character(len=4)  :: ElemType       ! ='STRN','DISP'
    character(len=4)  :: TestCase       ! Define the test case (CANT,ARCH).
    character(len=2)  :: BConds         ! ='CC': Clamped-Clamped

 ! solution type:
    integer :: Solution

 ! Output:
    real(8),      allocatable:: PosIni   (:,:)      ! Initial nodal Coordinates.
    real(8),      allocatable:: PsiIni (:,:,:)      ! Initial element orientation vectors (CRV)
    real(8),      allocatable:: PosDef   (:,:)      ! Current nodal position vector. (sm: local coordinates)
    real(8),      allocatable:: PsiDef (:,:,:)      ! Current element orientation vectors.
    real(8),      allocatable:: InternalForces(:,:) ! Internal force/moments at nodes.

    real(8),      allocatable :: DensityVector (:)  ! Density of each element of the beam. To be used for cost evaluation.
    real(8),      allocatable :: LengthVector  (:)  ! Length of each element of the beam. To be used for cost evaluation.


 ! Dynamic Input set-up
 real(8) :: t0, tfin, dt ! Not to be passed to main routine
 integer :: NumSteps     ! Not necessary to pass to main routine
 real(8), allocatable :: Time(:) ! Passe in input. Array with all time-steps

 ! Structural Dynamic Input/Output
    real(8),      allocatable:: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
    real(8),      allocatable:: ForceTime   (:)   ! Time history of the dynamic nodal forces.
    real(8),      allocatable:: ForcedVel   (:,:) ! Forced velocities at the support.
    real(8),      allocatable:: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.

    real(8),      allocatable:: PosDotDef (:,:)   ! Current nodal position vector.
    real(8),      allocatable:: PsiDotDef (:,:,:) ! Current element orientation vectors.
    real(8),      allocatable:: PosPsiTime(:,:)   ! Position vector/rotation history at beam tip.
    real(8),      allocatable:: VelocTime(:,:)    ! History of velocities.
    real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step

 ! Rigid-body variables
    real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
    real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
    real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.

 ! utilities
   integer :: ii ! counter

 ! --------------------------------------------------------------- Define Input:

! Options and Problem Setup
    NumElems = 10
    Solution = 112

    NumNodesElem = 3
    ElemType='DISP'
    TestCase='lala'
    BConds='CF'

 ! Design
    BeamLength1 = 5.0d0
    BeamLength2 = 0.0d0
    ThetaRoot   = 0.d0
    ThetaTip    = 0.0_8 ! pi/6.d0
    ExtForce=(/ 0.d0, 0.d0, 600.d3 /)
    ExtMomnt=(/ 0.d0, 0.d0,   0.d0 /)

    TipMass =0.0
    TipMassY=0.0
    TipMassZ=0.0
    SectWidth=0.0
    SectHeight=0.0
    Omega=0.0

    call isorect(0.1_8, 0.1_8, 'alumA', BeamMass, BeamStiffness)

    BeamStiffness=0.0_8
    BeamStiffness(1,1)= 4.8d8   ! EA [Nm]
    BeamStiffness(2,2)= 3.231d8 !*0.5 ! GA
    BeamStiffness(3,3)= 3.231d8
    BeamStiffness(4,4)= 1.d6    ! GJ
    BeamStiffness(5,5)= 9.346d6 !*0.5! EI
    BeamStiffness(6,6)= 9.346d6
    !
    BeamMass=0.0_8
    BeamMass(1,1)=100.d0        ! m [kg/m]
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=10.d0         ! J [kgm]
    BeamMass(5,5)=10.d0
    BeamMass(6,6)=10.d0


 ! Build Beam structure
 allocate( BeamSpanStiffness(NumElems,6,6) )
 allocate( BeamSpanMass(NumElems,6,6) )


 do ii = 1,NumElems
    BeamSpanStiffness(ii,:,:)=BeamStiffness
    BeamSpanMass(ii,:,:)     =BeamMass
 end do
 !BeamSpanStiffness(10,:,:)=0.1_8 * BeamStiffness

 ! Determine static output size:
 if ( (NumNodesElem==2) .or. (NumNodesElem==3) ) then
    NumNodes = (NumNodesElem -1)*NumElems + 1
 else
    stop 'NumNodesElem must be equal to 2 or 3'
 end if

 ! Preallocate:
 allocate( PhiNodes(NumNodes) )
 allocate(        PosIni(          NumNodes,3)); PosIni=0.0_8;
 allocate(        PosDef(          NumNodes,3)); PosDef=0.0_8;
 allocate(        PsiIni(Numelems, MaxElNod,3)); PsiIni=0.0_8;
 allocate(        PsiDef(Numelems, MaxElNod,3)); PsiDef=0.0_8;
 allocate(InternalForces(          NumNodes,6)); InternalForces=0.0_8;


 PosIni= 0.d0
 PhiNodes = 0.0_8
 do ii =1,NumNodes
    PosIni(ii,1)=BeamLength1*(dble(ii-1)/dble(NumNodes-1))
    PhiNodes(ii) = ThetaRoot+(ThetaTip-ThetaRoot)*(dble(ii-1)/dble(NumNodes-1))
 end do

 ! Allocate Variables independent of the solution
 allocate( DensityVector (NumElems) )
 allocate(  LengthVector (NumElems) )

 ! Static Forces
 allocate( ForceStatic(NumNodes,6) )
 ForceStatic(NumNodes,1:3)=ExtForce
 ForceStatic(NumNodes,4:6)=ExtMomnt


 ! --------------------------------------------------------------- Dynamic Input

 t0  =  0.0d0
 tfin= 10.0d0
 dt  = 1.0d-2
 NumSteps= ceiling((tfin-t0)/dt)

 Omega=20.d0    ! rad/s
 allocate (Time(NumSteps+1))
 do ii=1,NumSteps+1
     Time(ii)=t0+dt*dble(ii-1)
 end do

 ! Force or velocity input.
 allocate (ForceTime   (NumSteps+1));   ForceTime   = 0.d0
 allocate (ForceDynAmp (NumNodes,6));   ForceDynAmp = 0.d0
 if (Solution == 142) then
    allocate (ForcedVel   (1,6)); ForcedVel   = 0.d0
 else
    allocate (ForcedVel   (NumSteps+1,6)); ForcedVel   = 0.d0
 end if
 allocate (ForcedVelDot(NumSteps+1,6)); ForcedVelDot= 0.d0

 call input_dynforce  (NumNodes,Time,0.0_8*ForceDynAmp,ForceDynAmp,ForceTime)
 call input_forcedvel (NumNodes,Time,ForcedVel,ForcedVelDot)

 allocate (PosDotDef(NumNodes,3));           PosDotDef= 0.d0
 allocate (PsiDotDef(NumElems,MaxElNod,3));  PsiDotDef= 0.d0
 allocate (PosPsiTime(NumSteps+1,6));        PosPsiTime=0.d0
 allocate (VelocTime(NumSteps+1,NumNodes));  VelocTime= 0.d0
 allocate (DynOut((NumSteps+1)*NumNodes,3)); DynOut=0.d0



! --------------------------------------------------- Rigid Body + Dynamic Input
if ((Solution.ge.900).and.(Solution.le.952)) then

    ! Initialize
    allocate (RefVel   (NumSteps+1,6));  RefVel   =ForcedVel;       ! RefVel(1,5)=0.5d0
    allocate (RefVelDot(NumSteps+1,6));  RefVelDot=ForcedVelDot
    allocate (Quat     (4));             Quat     =(/1.d0,0.d0,0.d0,0.d0/)

end if


! ------------------------------------------------------------------ Call Solver
  call opt_main( NumElems, NumNodes,               &
                     & NumNodesElem , ElemType, TestCase, BConds,    & ! Problem Setup Shared
                     & BeamSpanStiffness, BeamSpanMass,  & ! Span Properties
                     & PhiNodes,                         &
                     & ForceStatic,                      &
                     & TipMass, TipMassY, TipMassZ,      &
                     & .false.,.true.,.true.,            & ! Options
                     & .true., .false., 0, 99,           &
                     & 10,Solution, 1.d-5, 1.d-5,  1.d-4,&
                     & PosIni, PsiIni,                   & ! Initial Pos/Rot
                     & PosDef, PsiDef, InternalForces,   & ! Output Static
                     & DensityVector, LengthVector,      & ! Design output up to v1.0
                     & NumSteps, Time,                   & ! input_dynsetup
                     & ForceTime, ForceDynAmp,           & ! input_dynforce
                     & ForcedVel, ForcedVelDot,          & ! input_forcedvel
                     & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! output from sol 202, 212, 302, 312, 322
                     & RefVel, RefVelDot, Quat)

                    ! FollowerForce,FollowerForceRig,PrintInfo,                 & ! Options
                    ! OutInBframe,OutInaframe,ElemProj,MaxIterations,         &
                    ! NumLoadSteps,Solution,DeltaCurved,MinDelta, NewmarkDamp,&

                    !!! NCB1 options
                    !& FollowerForce=.false.,PrintInfo=.false.,                & ! Options
                    !& MaxIterations=99,         &
                    !& NumLoadSteps=10,Solution=112,MinDelta= 1.d-5)


  print *, 'Max. Tip Displ. ', cost_node_disp(PosIni,PosDef,NumNodes)

end program opt_routine_test

