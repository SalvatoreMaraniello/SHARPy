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

 implicit none


! ------------------------------------------------------------- Define Variables

    real(8) :: W_COST(2)    ! <--- example of known size array, 2 is NCOSTFUNS into opt_shared
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

! Problem Setup:
    integer :: NumElems     ! Number of elements
    integer :: NumNodes     ! Number of nodes in the model.

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

 ! --------------------------------------------------------------- Define Input:

! Options and Problem Setup
    NumElems = 10
    Solution = 112

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



 ! Determine output size:
 if ( (NumNodesElem==2) .or. (NumNodesElem==3) ) then
    NumNodes = (NumNodesElem -1)*NumElems + 1
 else
    stop 'NumNodesElem must be equal to 2 or 3'
 end if
 !print *, 'Total Number of nodes', NumNodes

! -------------------------------------------------------------------- Call Main
! The select case is necessary to preallocate the outputs

select case (Solution)

    case (102, 112, 142)

        ! Preallocate:
        allocate(        PosIni(          NumNodes,3)); PosIni=0.0_8;
        allocate(        PosDef(          NumNodes,3)); PosDef=0.0_8;
        allocate(        PsiIni(Numelems, MaxElNod,3)); PsiIni=0.0_8;
        allocate(        PsiDef(Numelems, MaxElNod,3)); PsiDef=0.0_8;
        allocate(InternalForces(          NumNodes,6)); InternalForces=0.0_8;


        call opt_main( NumElems, NumNodes,               &
                     & NumNodesElem , ElemType, TestCase, BConds,    & ! Problem Setup Shared
                     & BeamLength1, BeamLength2,         & ! design variables
                     & BeamStiffness, BeamMass,          &
                     & ExtForce, ExtMomnt,               &
                     & SectWidth, SectHeight,            &
                     & ThetaRoot, ThetaTip,              &
                     & TipMass, TipMassY, TipMassZ,      &
                     & Omega,                            &
                     & .false.,.true.,.false.,           & ! Options
                     & .true., .false., 0, 99,           &
                     & 10,Solution, 1.d-5, 1.d-5,  1.d-4,&
                     & PosIni, PsiIni,                   & ! Initial Pos/Rot
                     & PosDef, PsiDef, InternalForces    ) ! Output Static

                 !!! NCB1 options
                 !& FollowerForce=.false.,PrintInfo=.false.,                & ! Options
                 !& MaxIterations=99,         &
                 !& NumLoadSteps=10,Solution=112,MinDelta= 1.d-5)


    case default
        stop '[opt_routine_test] Check the allocation process is defined for the Solution sought.'


    end select


    ! Correct
    print *, 'Max. Tip Displ. ', cost_node_disp(PosIni,PosDef,NumNodes)


    ! NCB1 options


end program opt_routine_test

