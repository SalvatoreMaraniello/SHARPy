!-> Module.- INPUT Salvatore Maraniello 27/06/2014
!
!-> Author: Salvatore Maraniello, copied from input_rob.f90 (R. Simpson) &
!                             ... input_rafa.f90 (R. Palacios) ...
!                             ... look for sm in the file
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Input data for XBeam test cases (test case x000).
!
!-> Subroutines.-
!
!    -input_setup:    Setup test case:
!    -input_elem :    Compute element information.
!    -input_node :    Compute nodal information.
!    -input_modal:    Compute natural vibration modes.
!    -input_dynsetup: Setup parameters for dynamic solution
!    -input_dynforce: Define time history of applied forces.
!    -input_foredvel: Define time-varying forced velocities.
!    -output_elem:    Write element information in output file.
!
! -> Details & main features.-
!    -input_setup:    Define here geometry, discretisation mass/stiffness ...
!                 ... properties & other solution options (including Gauss
!                 ... integration points)
!    -input_elem :    Define here connectivity and element orientation.
!                 ... Linear elements are used by default. The inverse of the...
!                 ... stiffness/mass matrices are also allocated for each element
!    -input_node :    Define here coordinates, pretwist, BCs and nodal forces.
!
!
!-> Remarks.-
!
!  1) Test cases from Geradin & Cardona's book.
!  2) 1 point for 2-node beam, to prevent shear locking.
!
!-> sm TestCases Summary.-
!
! 'NCB1': from Geradin & Cardona (input_rafa.f90) and used in 'Numerical aspects
!         of nonlinear flexible aircraft flight dynamics modeling', Simpson &
!         Palacios, 2013
! 'BEND': 45-degree bend in Geradin & Cardona, p 133 (input_rafa.f90)
! 'CANT': (from input_henrik.f90)
! 'HALE', 'GOLD', 'TPY0': (input_rob.f90)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module input
 use xbeam_shared
 implicit none

 ! Shared variables.
 real(8),private,save:: BeamLength1, BeamLength2     ! Beam defining lengths.
 real(8),private,save:: BeamStiffness(6,6)           ! Beam element stiffness matrix (assumed constant).
 real(8),private,save:: BeamMass(6,6)                ! Beam element mass matrix (assumed constant).
 real(8),private,save:: ExtForce(3)                  ! Applied forces at the tip.
 real(8),private,save:: ExtMomnt(3)                  ! Applied moments at the tip.
 integer,private,save:: NumNodesElem                 ! Number of nodes on each element.
 real(8),private,save:: SectWidth,SectHeight         ! Height and width of the cross section.
 real(8),private,save:: ThetaRoot=0.d0               ! Pretwist angle at root.
 real(8),private,save:: ThetaTip =0.d0               ! Pretwist angle at tip.
 real(8),private,save:: TipMass  =0.d0               ! Mass at the beam tip.
 real(8),private,save:: TipMassY =0.d0               ! Y offset of the tip mass.
 real(8),private,save:: TipMassZ =0.d0               ! Z offset of the tip mass.

 real(8),private,save:: Omega   =0.d0                ! Frequency of oscillatory motions.

 character(len=4),private,save:: ElemType            ! ='STRN','DISP'
 character(len=4),private,save:: TestCase            ! Define the test case (CANT,ARCH).
 character(len=2),private,save:: BConds              ! ='CC': Clamped-Clamped

 !real(8), parameter :: pi = 4.d0!*datan(1.d0)
                                                     ! ='CF': Clamped-Free'
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_SETUP
!
!-> Description:
!
!    Setup test case:
!    - geometry, discretisation
!    - applied forces
!    - eleements stiffness and mass properties
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_setup (NumElems,OutFile,Options)

! I/O Variables.
  integer,          intent(out):: NumElems       ! Number of elements in the model.
  character(len=25),intent(out):: OutFile        ! Output file.
  type(xbopts),     intent(out):: Options        ! Solution options.

! Local variables.
! (sm: rho is used to define the beam prop. in the BEND case)
  real(8) :: E=0.d0,G=0.d0,rho=0.d0
  real(8) :: sigma

! INPUT START
  TestCase='TPY0'!'PTW2'! 'BEND'!'NCB1'
  Options%Solution=112

! Default values.
  ExtForce(1:3)=0.d0
  ExtMomnt(1:3)=0.d0
  BeamStiffness=0.d0
  BeamMass     =0.d0

  select case (trim(TestCase))

! Results in Pai's book (p. 392). X axis goes along the vertical beam.
  case ('FPAI')
    NumElems   = 40
    ThetaRoot  = 0.d0
    ThetaTip   = 0.d0
    BeamLength1= 479.d-3

    SectWidth= 50.8d-3
    SectHeight=0.45d-3
    E  = 1.27d11
    G  = E/(2.d0*(1.d0+0.36d0))
    rho= 4.43d3

    ExtForce=(/-9.81d0,0.d0,0.d0/)*rho*SectWidth*SectHeight*BeamLength1/dble(NumElems)
    Options%FollowerForce=.false.
    Options%NumLoadSteps = 10
    BConds ='CF'
    Options%MinDelta = 1.d-5

    BeamStiffness(1,1)=SectWidth*SectHeight*E
    BeamStiffness(2,2)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(3,3)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(5,5)=SectWidth*(SectHeight**3.d0)*E/12.d0
    BeamStiffness(6,6)=SectHeight*(SectWidth**3.d0)*E/12.d0
    BeamStiffness(4,4)=0.5d0*(BeamStiffness(5,5)+BeamStiffness(6,6))

    BeamMass(1,1)= rho*SectWidth*SectHeight
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(5,5)=rho*SectWidth*(SectHeight**3.d0)/12.d0
    BeamMass(6,6)=rho*SectHeight*(SectWidth**3.d0)/12.d0
    BeamMass(4,4)=0.5d0*(BeamMass(5,5)+BeamMass(6,6))


     case ('PTW1', 'PTW2') ! Check effect of pre-twist
      ! Inputs are as per NCB1 but:
      ! a. GA and EI of cross-section 2 are halved to make the bean asimetric
        NumElems    = 25
        NumNodesElem= 3
        ThetaRoot   = 0.d0
        ThetaTip    = Pi/6.d0 ! ps: Pi defined in xbeam_shared
        BeamLength1 = 5.0d0

        BConds  ='CF'
        ExtForce=(/ 0.d0, 0.d0, 600.d3 /)
        ExtMomnt=(/ 0.d0, 0.d0,   0.d0 /)
        Options%FollowerForce = .false.

        Options%NumLoadSteps  = 10!NumElems
        Options%MinDelta      = 1.d-5
        Options%MaxIterations = 200

        BeamStiffness(1,1)= 4.8d8   ! EA [Nm]
        BeamStiffness(2,2)= 0.5 * 3.231d8 ! GA
        BeamStiffness(3,3)= 3.231d8
        BeamStiffness(4,4)= 1.d6    ! GJ
        BeamStiffness(5,5)= 0.5 * 9.346d6 ! EI
        BeamStiffness(6,6)= 9.346d6

        BeamMass(1,1)=100.d0        ! m [kg/m]
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(4,4)=10.d0         ! J [kgm]
        BeamMass(5,5)=10.d0
        BeamMass(6,6)=10.d0

     case ('NCB1')
      ! Analytical NatFreqs (rad/s): 43.0(B1), 99.3(T1), 269.4(B2),
      ! 298.0(T2), 496.7(T3), 699.9(A1), 759.5(B3)
        NumElems    = 10
        NumNodesElem= 3
        ThetaRoot   = 0.d0
        ThetaTip    = 0.d0
        BeamLength1 = 5.0d0

        BConds  ='CF'
        ExtForce=(/ 0.d0, 0.d0, 600.d3 /)
        ExtMomnt=(/ 0.d0, 0.d0,   0.d0 /)
        Options%FollowerForce = .false.

        Options%NumLoadSteps  = 10
        Options%MinDelta      = 1.d-5
        Options%MaxIterations = 99

        BeamStiffness(1,1)= 4.8d8   ! EA [Nm]
        BeamStiffness(2,2)= 3.231d8 ! GA
        BeamStiffness(3,3)= 3.231d8
        BeamStiffness(4,4)= 1.d6    ! GJ
        BeamStiffness(5,5)= 9.346d6 ! EI
        BeamStiffness(6,6)= 9.346d6

        BeamMass(1,1)=100.d0        ! m [kg/m]
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(4,4)=10.d0         ! J [kgm]
        BeamMass(5,5)=10.d0
        BeamMass(6,6)=10.d0

    case ('CANT')
        NumElems    =1
        NumNodesElem=2

        ThetaRoot  =  0.d0
        ThetaTip   =  0.d0   ! Pi/6.d0
        BeamLength1= 10.d0

        BConds  ='CF'
        ExtForce=(/1.d0,1.d0,1.d0/)*1.d0
        ExtMomnt=(/1.d0,1.d0,1.d0/)*1.d0

        Options%FollowerForce    = .true.
        Options%FollowerForceRig = .true.
        Options%OutInaframe      = .false.
        Options%NumLoadSteps  = 10
        Options%MinDelta      = 1.d-5
        Options%MaxIterations = 99

        BeamStiffness(1,1)= 10.d3
        BeamStiffness(2,2)= 10.d3
        BeamStiffness(3,3)= 10.d3
        BeamStiffness(4,4)= 500.d0
        BeamStiffness(5,5)= 500.d0
        BeamStiffness(6,6)= 500.d0

        BeamStiffness=1.d0*BeamStiffness

        BeamMass(1,1)=1.d0
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(4,4)=10.d0
        BeamMass(5,5)=10.d0
        BeamMass(6,6)=10.d0

    case ('BEND')
        NumElems= 6
        NumNodesElem=2
        BeamLength1=100.d0

        SectWidth= 1.d0
        SectHeight=1.d0
        E  = 1.0d7
        G  = 1.0d5! 0.5d7 ! sm: as per Cardona this should be zero
        rho= 1.

        BConds  ='CF'
        ExtForce=(/0.d0,0.d0,600.d0/)
        ExtMomnt=(/0.d0,0.d0,0.d0/)
        Options%FollowerForce=.true.

        Options%MinDelta     = 1.d-4     ! 1.d-4
        Options%NumLoadSteps =  20       ! sm 2->10
        Options%MaxIterations = 999

        BeamStiffness(1,1)=SectWidth*SectHeight*E
        BeamStiffness(2,2)=(5.d0/6.d0)*SectWidth*SectHeight*G
        BeamStiffness(3,3)=(5.d0/6.d0)*SectWidth*SectHeight*G
        BeamStiffness(5,5)=SectWidth*(SectHeight**3.d0)*E/12.d0
        BeamStiffness(6,6)=SectHeight*(SectWidth**3.d0)*E/12.d0
        BeamStiffness(4,4)=0.5d0*(BeamStiffness(5,5)+BeamStiffness(6,6))

        BeamMass(1,1)= rho*SectWidth*SectHeight
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(5,5)=rho*SectWidth*(SectHeight**3.d0)/12.d0
        BeamMass(6,6)=rho*SectHeight*(SectWidth**3.d0)/12.d0
        BeamMass(4,4)=0.5d0*(BeamMass(5,5)+BeamMass(6,6))

    case ('HALE')
      ! half wing of HALE model aircraft presented in Table 4 in Murua et al (2011), AIAA-2011-1915
      ! Rob: modified for SharPy PyBeam NonLinearStatic test case 0
        NumElems    = 8
        NumNodesElem= 2
        ThetaRoot   = 0.d0
        ThetaTip    = 0.d0
        BeamLength1 = 16.0d0

        sigma=1.0d0     ! Parameter to change torsion and bending stiffness

        TipMass =0.0d0;
        TipMassY=0.0d0;
        TipMassZ=0.0d0;

        BConds  ='CF'
        ExtForce=(/0.d0,0.d0,800.d0/)*1.d0
        ExtMomnt=(/0.d0,0.d0,800.d0/)*1.d0

        Options%FollowerForce    = .true.
        Options%FollowerForceRig = .true.
        Options%OutInaframe      = .false.
        Options%NumLoadSteps  = 10
        Options%MinDelta      = 1.d-4
        Options%MaxIterations = 99

        BeamMass(1,1)=0.75d0
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(4,4)=1.d-1
        BeamMass(5,5)=1.d-3
        BeamMass(6,6)=1.d-3

        BeamStiffness(1,1)=1.0d9
        BeamStiffness(2,2)=BeamStiffness(1,1)
        BeamStiffness(3,3)=BeamStiffness(1,1)
        BeamStiffness(4,4)=1.0d4
        BeamStiffness(5,5)=2.0d4
        BeamStiffness(6,6)=4.0d6

        BeamStiffness=BeamStiffness*sigma

    case ('GOLD')
      ! half Goland wing taken from Table 3 in Murua et al (2011), AIAA-2011-1915
        NumElems    = 6
        NumNodesElem= 3
        ThetaRoot   = 0.d0
        ThetaTip    = 0.d0
        BeamLength1 = 6.096d0

        TipMass =0.0d0;
        TipMassY=0.0d0;
        TipMassZ=0.0d0;

        BConds  ='CF'
        ExtForce=(/0.d0,0.d0,0.d0/)*1.d5
        ExtMomnt=(/0.d0,0.d0,0.d0/)*1.d0

        Options%FollowerForce    = .false.
        Options%FollowerForceRig = .true.
        Options%OutInaframe      = .false.
        Options%NumLoadSteps  = 10
        Options%MinDelta      = 1.d-4
        Options%MaxIterations = 99

        BeamMass(1,1)=35.709121d0
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(4,4)=8.6405832d0
        BeamMass(5,5)=1.d-3
        BeamMass(6,6)=1.d-3

    !    BeamMass(1,6)= 6.53048405d0
        BeamMass(3,4)=-BeamMass(1,6)
        BeamMass(4:6,1:3)=transpose(BeamMass(1:3,4:6))

        BeamStiffness(1,1)=1.0d9
        BeamStiffness(2,2)=BeamStiffness(1,1)
        BeamStiffness(3,3)=BeamStiffness(1,1)
        BeamStiffness(4,4)=0.987581d6
        BeamStiffness(5,5)=9.77221d6
        BeamStiffness(6,6)=9.77221d8

    case ('TPY0')
      ! half wing of HALE model aircraft presented in Table 4 in Murua et al (2011), AIAA-2011-1915
      ! Rob: modified for SharPy PyBeam NonLinearStatic test case 0
        NumElems    = 8
        NumNodesElem= 3
        ThetaRoot   = 0.d0
        ThetaTip    = 0.d0
        BeamLength1 = 16.0d0

        sigma=1.0d0     ! Parameter to change torsion and bending stiffness

        TipMass =0.0d0;
        TipMassY=0.0d0;
        TipMassZ=0.0d0;

        BConds  ='CF'
        ExtForce=(/0.d0,0.d0,80.d0/)*1.d0
        ExtMomnt=(/0.d0,0.d0,0.d0/)*1.d0

        Options%FollowerForce    = .false.
        Options%FollowerForceRig = .true.
        Options%OutInaframe      = .false.
        Options%NumLoadSteps  = 10
        Options%MinDelta      = 1.d-4
        Options%MaxIterations = 99

        BeamMass(1,1)=0.75d0
        BeamMass(2,2)=BeamMass(1,1)
        BeamMass(3,3)=BeamMass(1,1)
        BeamMass(4,4)=1.d-1
        BeamMass(5,5)=1.d-3
        BeamMass(6,6)=1.d-3

        BeamStiffness(1,1)=1.0d9
        BeamStiffness(2,2)=BeamStiffness(1,1)
        BeamStiffness(3,3)=BeamStiffness(1,1)
        BeamStiffness(4,4)=1.0d4
        BeamStiffness(5,5)=2.0d4
        BeamStiffness(6,6)=4.0d6

        BeamStiffness=BeamStiffness*sigma

    end select

! Set name for output file.
  OutFile(1:8)=trim(TestCase)//'_SOL'
  write (OutFile(9:11),'(I3)') Options%Solution

! Solver options.
  select case (Options%Solution)
  case (102,112,142,202,212,302,312,322,900,902,910,912,922,952)
    ElemType= 'DISP'
  case default
    STOP 'Error: Wrong solution code (51707)'
  end select

! Define number of Gauss points (2-noded displacement-based element needs reduced integration).
  select case (NumNodesElem)
    case (2)
      Options%NumGauss=1
    case (3)
      Options%NumGauss=2
  end select

! Minimum angle for two unit vectors to be parallel. real(8):: Vector  (3)            ! Element orientation vector. It goes along the local Y
  Options%DeltaCurved=1.d-5

! sm check:  stop execution if NumNodesElem is not valid
  if ((NumNodesElem /=2) .and. (NumNodesElem /=3)) then
    STOP 'Error: Input data for Number of Nodes for Each Element not defined.'
  end if ! end sm

  print "(a4,i3,a)", TestCase, Options%Solution, ': input_setup done!' ! sm
  return
 end subroutine input_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_ELEM
!
!-> Description:
!
!    Define element properties.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine input_elem (NumElems,NumNodes,Elem)
 use lib_lu
 use lib_rot

! I/O Variables.
  integer,intent(in)  :: NumElems        ! Number of elements in the model.
  integer,intent(out) :: NumNodes        ! Total Number of nodes in the model.
  type(xbelem),intent(out) :: Elem(:)    ! Element information

! Local variables.
  integer:: i                            ! Counter.
  integer,save :: fl=2,fw=2              ! Multiplier.
  real(8):: BeamInvStiffness(6,6)        ! Inverse of the stiffness matrix.
  real(8):: LocPos(3)                    ! Local position vector of the lumped mass.

! Connectivies.
  select case (trim(TestCase))
    case default
      select case (NumNodesElem)
      case (2)
        do i=1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=i
            Elem(i)%Conn(2)=i+1
            Elem(i)%NumNodes=2 ! sm: this field is redundant (computed into fem_glob2loc_extract
        end do
        NumNodes=NumElems+1
      case (3)
        do i=1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=2*(i-1)+1
            Elem(i)%Conn(2)=2*(i-1)+3
            Elem(i)%Conn(3)=2*(i-1)+2
            Elem(i)%NumNodes=3 ! sm: this field is redundant (computed into fem_glob2loc_extract
        end do
        NumNodes=2*NumElems+1
    end select
  end select

! Store element stiffness/mass (constant)
  BeamInvStiffness=0.d0
  call lu_invers (BeamStiffness, BeamInvStiffness)
  do i=1,NumElems
    Elem(i)%Stiff   = BeamStiffness
    Elem(i)%Mass    = BeamMass
    Elem(i)%InvStiff= BeamInvStiffness
  end do

! Define lumped masses at element nodes.
  do i=1,NumElems
    Elem(i)%RBMass = 0.d0
  end do

  select case (trim(TestCase))

  case default
        LocPos(1)=0.d0
        LocPos(2)=TipMassY
        LocPos(3)=TipMassZ

!       Elem(NumElems/2)%RBMass(2,1:3,1:3)= TipMass*Unit
!       Elem(NumElems/2)%RBMass(2,1:3,4:6)=-TipMass*rot_skew(LocPos)
!       Elem(NumElems/2)%RBMass(2,4:6,1:3)= TipMass*rot_skew(LocPos)
!       Elem(NumElems/2)%RBMass(2,4:6,4:6)=-TipMass*matmul(rot_skew(LocPos),rot_skew(LocPos))
!
!       Elem(NumElems)%RBMass(1,:,:)= Elem(NumElems/2)%RBMass(2,:,:)
  end select

! Element orientation.
! sm: this block defines the local Y axis.
  do i=1,NumElems

    select case (trim(TestCase))

    case default  ! straight beam ('HALE','GOLD','TPY0','CANT','NCB1')
        Elem(i)%Vector(2)= 1.d0

    case ('PTW2') ! rotate the Elem(i)%Vector
        Elem(i)%Vector= 0.d0
        Elem(i)%Vector(2)= dcos( ThetaRoot+(ThetaTip - ThetaRoot)*(dble(i-1)/dble(NumElems)) ) !-cos(ThetaRoot + (ThetaTip - ThetaRoot) * y/L)
        Elem(i)%Vector(3)= dsin( ThetaRoot+(ThetaTip - ThetaRoot)*(dble(i-1)/dble(NumElems)) ) ! sin( * y/L)
        !sm 7 july 2014: lines below commented for twist-orientation vector investigation, case 3
        ThetaTip=0.d0 !ThetaTip/100.d0
        ThetaRoot=0.d0 !ThetaRoot/100.d0
    case ('BEND') ! the beam develops in the Oxy plane
    ! sm: the case is set up as per input_rafa,f90 but I'd set Elem%Vector=(/0,0,1/)
        Elem(i)%Vector= 0.d0
        Elem(i)%Vector(1)=-dcos((Pi/4.d0)*(dble(i-1)/dble(NumElems))) !-cos(Pi/4 * y/L)
        Elem(i)%Vector(2)= dsin((Pi/4.d0)*(dble(i-1)/dble(NumElems))) ! sin(Pi/4 * y/L)
        !Elem(i)%Vector(3)=1.d0
    end select

  end do

! Define element types.
  select case (ElemType)
  case ('DISP')
    Elem(1:NumElems)%MemNo=0
  case ('STRN') ! from input_rafa.f90
    Elem(1:NumElems)%MemNo=1
  end select

  return
 end subroutine input_elem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_NODE
!
!-> Description:
!
!    Define nodal properties.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_node (NumNodes,Elem,BoundConds,Coords,Forces,PhiNodes)

! I/O Variables.
  integer,      intent(in) :: NumNodes            ! Number of nodes in the model.
  type(xbelem), intent(in) :: Elem      (:)       ! Element information
  integer,      intent(out):: BoundConds(:)       ! =0 on free nodes; =1 on clamped nodes; =-1 on free nodes.
  real(8),      intent(out):: Coords    (:,:)     ! Initial nodal Coordinates.
  real(8),      intent(out):: Forces    (:,:)     ! Applied nodal forces.
  real(8),      intent(out):: PhiNodes  (:)       ! Initial twist at grid points.

! Local variables.
  integer      :: i           ! Counters.
  integer,save :: fl=2,fw=2   ! Multiplier.
  real(8)      :: Theta       ! Parameter on curved beams.
  integer, allocatable :: ffvec(:) ! sm: store the indices of nodes with applied forces

! Initial position vector of grid points.
  Coords= 0.d0
    select case (trim(TestCase))

    case default !('HALE','GOLD','TPY0','NCB1')
        do i=1,NumNodes
            Coords(i,1)=BeamLength1*(dble(i-1)/dble(NumNodes-1))
        end do

    case ('BEND')
        do i=1,NumNodes
            Theta=(Pi/4.d0)*(dble(i-1)/dble(NumNodes-1))
            Coords(i,1)=BeamLength1*(1.d0-dcos(Theta))
            Coords(i,2)=BeamLength1*(     dsin(Theta))
        end do

    end select

! Initial pretwist angle.
select case (trim(TestCase))
  case default
    do i=1,NumNodes
      PhiNodes(i)=ThetaRoot+(ThetaTip-ThetaRoot)*(dble(i-1)/dble(NumNodes-1))
    end do
  !case ('PTW2')
  !  PhiNodes(i)=0 ! this leads to error after the execution is terminated...
end select

! Static point forces.
  !ffvec=0
  Forces=0.d0
  select case (trim(TestCase))

  case default !('BEND','CANT','NCB1') ! Forces at the tip only
      allocate(ffvec(1))
      ffvec = (/ NumNodes /)
      Forces(ffvec,1:3)=reshape(ExtForce,(/1,3/))
      Forces(ffvec,4:6)=reshape(ExtMomnt,(/1,3/))
      deallocate(ffvec)

  case ('HALE','GOLD','TPY0')
          Forces(1,1:3)=-ExtForce
          Forces(1,4:6)=-ExtMomnt
          Forces(NumNodes,1:3)=ExtForce
          Forces(NumNodes,4:6)=ExtMomnt

  end select

! Boundary conditions
  BoundConds=0
  select case (trim(TestCase))

  case default
    ! Node 1
    if (BConds(1:1).eq.'C') BoundConds(1)       = 1
    if (BConds(1:1).eq.'F') BoundConds(1)       =-1
    if (BConds(1:1).eq.'S') BoundConds(1)       = 2
    if (BConds(1:1).eq.'T') BoundConds(1)       = 3
    ! Node NumNodes
    if (BConds(2:2).eq.'C') BoundConds(NumNodes)= 1
    if (BConds(2:2).eq.'F') BoundConds(NumNodes)=-1
    if (BConds(2:2).eq.'S') BoundConds(NumNodes)= 2
    if (BConds(2:2).eq.'T') BoundConds(NumNodes)= 3
  end select

  return
 end subroutine input_node


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_DYNSETUP
!
!-> Description:
!
!    Setup test case.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_dynsetup (NumSteps,t0,dt,Options)

! I/O Variables.
  integer,intent(out):: NumSteps             ! Number of time steps.
  real(8),intent(out):: t0                   ! Initial time.
  real(8),intent(out):: dt                   ! Initial, final and delta time.
  type(xbopts),intent(inout):: Options      ! Solution options.

! Local variables.
  real(8):: tfin                             ! Final time.

  select case (trim(TestCase))

  case ('GOLD')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-3
    tfin= 1.0d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.01d0
    Omega=20.d0    ! rad/s

  case ('HALE')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-2
    tfin= 10.0d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.05d0
    Omega=20.d0    ! rad/s

  case ('TPY0')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-2
    tfin= 1.0d0
    NumSteps= ceiling((tfin-t0)/dt) + 2
    Options%NewmarkDamp=0.01d0
    Omega=20.d0    ! rad/s

  case default
    STOP 'Error: Input data for dynamic analysis not defined.'
  end select

  return
 end subroutine input_dynsetup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_DYNFORCE
!
!-> Description:
!
!    Define time-varying forcing terms.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_dynforce (NumNodes,Time,ForceStatic,ForceDynAmp,ForceTime)

! I/O Variables.
  integer,intent(in) :: NumNodes            ! Number of nodes in the model.
  real(8),intent(in) :: Time(:)             ! Time history.
  real(8),intent(in) :: ForceStatic(:,:)    ! Static force.
  real(8),intent(out):: ForceDynAmp(:,:)    ! Nodal force amplitudes.
  real(8),intent(out):: ForceTime  (:)      ! Time history of the applied forces.

! Local variables.
  integer:: i             ! Counter.
  integer:: NumSteps
  real(8):: Time0         ! Time for discontinuity in function.

! Initialize
  NumSteps=size(Time)-1

  select case (trim(TestCase))

  case default
    !ForceDynAmp = 2.d0*ForceStatic
    ForceTime(:)   = 1.d0
    ForceDynAmp(NumNodes,2) = 160.d0

! Ramped harmonic load.
    if (.false.) then
      Time0= Time(NumSteps)/2.d0
      ForceTime=0.d0
      do i=1,NumSteps
        ForceTime(i+1)=sin(Omega*Time(i+1))
        if (Time(i+1) < Time0) ForceTime(i+1)=ForceTime(i+1)*Time(i+1)/Time0
      end do
    end if

! Initial Load
    if (.false.) then
      ForceTime(:) = 0.d0
      ForceTime(1)= 1.d0
      ForceDynAmp(NumNodes,3) = 100.d3
    end if

! 1-Cos load.
    if (.false.) then
     ForceDynAmp(NumNodes,3) = 1.d03
     do i=1,NumSteps+1
         if ((Time(i).ge.0.d0).and.(Time(i).le.1.d-2)) then
             ForceTime(i)=(1.d0-cos(2*Pi*dble(Time(i)-0.d0)/dble(1.d-2)))/2.d0
         end if
     end do
    end if

! Ramped load.
    if (.false.) then
      ForceTime=1.d0
      do i=1,NumSteps+1
        if ((Time(i).ge.0.d0).and.(Time(i).le.2.5d0)) then
            ForceTime(i)=Time(i)/2.5d0
        end if
      end do
    end if

! sinusiodal load.
    if (.false.) then
     do i=1,NumSteps+1
        ! ForceTime(i) = sin(2.*1.*Time(i))
        ForceTime(i)=sin(Pi*dble(Time(i)-0.d0)/dble(0.5))
     end do
    end if
  end select
print *, ForceDynAmp
  return
 end subroutine input_dynforce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_FORCEDVEL
!
!-> Description:
!
!    Define time-varying forcing velocity.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine input_forcedvel (NumNodes,Time,ForcedVel,ForcedVelDot)

! I/O Variables.
  integer,intent(in) :: NumNodes           ! Number of nodes in the model.
  real(8),intent(in) :: Time(:)            ! Time history.
  real(8),intent(out):: ForcedVel(:,:)     ! Forced root velocities.
  real(8),intent(out):: ForcedVelDot(:,:)  ! Time derivatives of the forced root velocities.

! Local variables.
  integer:: i             ! Counter.
  integer:: NumSteps
  real(8):: VelAmp(6)     ! Amplitude of the velocities.

! Initialize.
  NumSteps=size(Time)-1
  ForcedVel   =0.d0
  ForcedVelDot=0.d0
  VelAmp=0.d0

  return
 end subroutine input_forcedvel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine OUTPUT_ELEMS
!
!-> Description:
!
!    Write information in elements.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_elems (iuOut,Elem,Coords,Psi)
  use lib_fem

! I/O Variables.
  integer,      intent(in)   :: iuOut             ! Output file.
  type(xbelem), intent(in)   :: Elem   (:)        ! Element information.
  real(8),      intent(in)   :: Coords (:,:)      ! Coordinates of the grid points.
  real(8),      intent(in)   :: Psi    (:,:,:)    ! CRV of the nodes in the elements.

! Local variables.
  integer:: i                    ! Counter.
  integer:: iElem                ! Counter on the finite elements.
  integer:: NumE                 ! Number of elements in the model.
  integer:: NumNE                ! Number of nodes in an element.
  real(8):: PosElem (MaxElNod,3) ! Coordinates/CRV of nodes in the element.

  NumE=size(Elem)

! Loop in the elements in the model.
  do iElem=1,NumE
    ! associates the coordinates of the nn-th node (global) to the ii-th node
    ! (local) of the element.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,PosElem,NumNE)

    do i=1,NumNE
      write (iuOut,'(2I4,1P6E15.6)') iElem,i,PosElem(i,:),Psi(iElem,i,:)
    end do
  end do

end subroutine output_elems

end module input
