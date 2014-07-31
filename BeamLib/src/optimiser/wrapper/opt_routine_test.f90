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

 implicit none


! ------------------------------------------------------------- Define Variables

    real(8) :: W_COST(2)    ! <--- example of known size array, 2 is NCOSTFUNS into opt_shared
    real(8), pointer, dimension(:) :: pCOST

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

! Problem Setup:
    integer :: NumElems     ! Number of elements
    integer :: NumNodes     ! Number of nodes in the model.

 ! Problem Setup Shared:
    integer           :: NumNodesElem   ! Number of nodes on each element.
    character(len=4)  :: ElemType       ! ='STRN','DISP'
    character(len=4)  :: TestCase       ! Define the test case (CANT,ARCH).
    character(len=2)  :: BConds         ! ='CC': Clamped-Clamped





 ! --------------------------------------------------------------- Define Input:

! Options and Problem Setup
    NumElems = 10
    NumNodes = 3

    NumNodesElem = 3
    ElemType='DISP'
    TestCase='NCB1'
    BConds='CF'

 ! Design
    BeamLength1 = 5.0d0
    BeamLength2 = 0.0d0
    ThetaRoot   = 0.d0
    ThetaTip    = 0.d0
    ExtForce=(/ 0.d0, 0.d0, 600.d3 /)
    ExtMomnt=(/ 0.d0, 0.d0,   0.d0 /)

    TipMass =0.0
    TipMassY=0.0
    TipMassZ=0.0
    SectWidth=0.0
    SectHeight=0.0
    Omega=0.0

    BeamStiffness=0.0_8
    BeamStiffness(1,1)= 4.8d8   ! EA [Nm]
    BeamStiffness(2,2)= 3.231d8 ! GA
    BeamStiffness(3,3)= 3.231d8
    BeamStiffness(4,4)= 1.d6    ! GJ
    BeamStiffness(5,5)= 9.346d6 ! EI
    BeamStiffness(6,6)= 9.346d6

    BeamMass=0.0_8
    BeamMass(1,1)=100.d0        ! m [kg/m]
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=10.d0         ! J [kgm]
    BeamMass(5,5)=10.d0
    BeamMass(6,6)=10.d0


! Solver options.
!  select case (Options%Solution)
!  case (102,112,142,202,212,302,312,322,900,902,910,912,922,952)
!    ElemType= 'DISP'
!  case default
!    STOP 'Error: Wrong solution code (51707)'
!  end select


! -------------------------------------------------------------------- Call Main
    call opt_main( NumElems, NumNodes,  pCOST, W_COST,&
                 & NumNodesElem , ElemType, TestCase, BConds,    & ! Problem Setup Shared
                 & BeamLength1, BeamLength2,         & ! design variables
                 & BeamStiffness, BeamMass,          &
                 & ExtForce, ExtMomnt,               &
                 & SectWidth, SectHeight,            &
                 & ThetaRoot, ThetaTip,              &
                 & TipMass, TipMassY, TipMassZ,      &
                 & Omega,                            &
                 & FollowerForce=.false.,PrintInfo=.false.,                & ! Options
                 & MaxIterations=99,         &
                 & NumLoadSteps=10,Solution=112,MinDelta= 1.d-5)
                 !& FollowerForce,FollowerForceRig,PrintInfo,              & ! Options
                 !& OutInBframe,OutInaframe,ElemProj,MaxIterations,        &
                 !& NumLoadSteps,NumGauss,Solution,DeltaCurved,MinDelta,   &
                 !& NewmarkDamp                                            )
    ! NCB1 options


end program opt_routine_test

