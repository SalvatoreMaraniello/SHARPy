!-> Module.- INPUT Salvatore Maraniello 10/07/2014
!
!-> Author: Salvatore Maraniello,
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Input data for optxbeam. The input for the forward mode need to be input into
!  input_xxx.f90.
!
!
!
!
!-> Subroutines.-
!
!
! -> Details & main features.-
!
!
!-> Remarks.-
!
!
!-> TestCases Summary.-
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module opt_input

    use input
    use opt_perturb
    use opt_cost_utl

    implicit none

    ! Remark: all flags have implicitely the save attribute as an initial value
    !         is assigned to them

    ! Flags to Perturb Shared Variables
    logical :: FLAG_DESIGN_SHARED(8+6+6*6+6*6+2) = .false. ! collects all the flags below into a vector                                                        ! see opt_set_FLAT_DESIGN_SHARED
    ! FLAGS for single inputs.
    logical,private :: FLAG_BeamLength1=.false., FLAG_BeamLength2=.false.! Beam defining lengths.
    logical,private :: FLAG_ThetaRoot=.false.                    ! Pretwist angle at root.
    logical,private :: FLAG_ThetaTip =.false.                    ! Pretwist angle at tip.
    logical,private :: FLAG_TipMass =.false.                     ! Mass at the beam tip.
    logical,private :: FLAG_TipMassY =.false.                    ! Y offset of the tip mass.
    logical,private :: FLAG_TipMassZ=.false.                     ! Z offset of the tip mass.
    logical,private :: FLAG_Omega=.false.                        ! Frequency of oscillatory motions.
    logical,private :: FLAG_ExtForce(3)=.false.                  ! Applied forces at the tip.
    logical,private :: FLAG_ExtMomnt(3)=.false.                  ! Applied moments at the tip.
    ! For constant cross section beam only
    logical,private :: FLAG_BeamStiffness(6,6)=.false.           ! Beam element stiffness matrix (assumed constant).
    logical,private :: FLAG_BeamMass(6,6)=.false.                ! Beam element mass matrix (assumed constant).
    logical,private :: FLAG_SectWidth=.false.,FLAG_SectHeight=.false.    ! Height and width of the cross section.

   ! Deltas for Perturbations
   real(8) :: d_BeamLength1, d_BeamLength2
   real(8) :: d_BeamStiffness(6,6)
   real(8) :: d_BeamMass(6,6)
   real(8) :: d_ExtForce(3)
   real(8) :: d_ExtMomnt(3)
   real(8) :: d_SectWidth,SectHeight
   real(8) :: d_ThetaRoot
   real(8) :: d_ThetaTip
   real(8) :: d_TipMass
   real(8) :: d_TipMassY
   real(8) :: d_TipMassZ
   real(8) :: d_Omega

   real(8), allocatable :: XSH(:,:) ! SHARED DESIGN VARIABLES


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine opt_setup(gradmode,solmode,fdmode,NOPTMAX)

        character(len=3), intent(out) :: gradmode  ! gradient method
        character(len=3), intent(out) :: solmode   ! solution mode
        character(len=3), intent(out) :: fdmode    ! finite differences method

        integer         , intent(out) :: NOPTMAX   ! Max Number of iterations for the optimisation

        ! ----------------------------------------------------------------------
        gradmode='FDF' ! gradient method: DIR: direct
                       !                  ADJ: adjointfinite differences
                       !                  FDF: Finite Differences
        solmode ='OPT' ! solution mode:   FWD: forward
                       !                  OPT: optimisation
                       !                  SNS: sensitivity analysis
        fdmode  ='FWD' ! FD method:       FWD: forward differences
                       !                  BKW: backward differences
                       !                  CNT: central differences
        NOPTMAX= 3     ! Maximum number of iterations for the optimiser

        ! ----------------------------------------------------------------------
        ! Design Parameters: shared variables
        !
        ! Remark: the value of these flags is set to zero at the end of the
        !   subroutine and changed during the code execution (e.g. to update the
        !   the design or perform FDs based sensitivity analysis.
        !   FLAG_DESIGN_SHARED is used to keep memory of the design variables.
        FLAG_BeamLength1 = .true.
        !FLAG_TipMassY =.true.
        FLAG_ExtForce(3)=.true.
        !FLAG_BeamMass(1,1)=.true.
        !FLAG_BeamMass(2,4)=.true.
        !FLAG_ExtMomnt(2)=.true.
        !FLAG_BeamStiffness(6,6)=.true.

        ! ----------------------------------------------------------------------
        ! Design Parameters: element dependent
        !
        ! -> TO BE IMPLEMENTED


        ! ----------------------------------------------------------------------
        ! PreProcessing:

        ! assign all the flags to FLAG_DESIGN_SHARED
        call opt_pack_FLAG_DESIGN_SHARED
        !!!!call opt_print_FLAG_DESIGN_SHARED
        call opt_set_shared_FLAGS_to_false

        ! set up XSH vector
        allocate( XSH(size(FLAG_DESIGN_SHARED), 0:NOPTMAX ) ); XSH=0_8
        call opt_pack_DESIGN_SHARED(0)

        ! TESTING:
        !call print_shared_input
        !call opt_print_XSH(0)
        !call opt_unpack_DESIGN_SHARED(1) ! this deletes everything
        !call print_shared_input
        !call opt_unpack_DESIGN_SHARED(0)   ! this should reassign values
        !call opt_print_XSH(0)
        !call print_shared_input
        !stop

    end subroutine opt_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_input_cost
!
!-> Description:
!    Use this routine to set up the input for the cost/constraint functions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine opt_input_cost(FLAG_COST,FLAG_CONSTR,W_COST,W_CONSTR)

        logical         , intent(out) :: FLAG_COST  (2) ! Flag array for cost funcitons
        logical         , intent(out) :: FLAG_CONSTR(2) ! Flag array for cost funcitons

        real(8)         , intent(out) :: W_COST  (size(FLAG_COST))    ! arrays with weights/scaling factors...
        real(8)         , intent(out) :: W_CONSTR(size(FLAG_COST))   ! ...for cost and constraint functions

        ! input for displacement cost function
        real(8)            ::  WEIGHT          ! wieght to apply to function
        character(len= 4)  ::   ADDTO          ! determines whether to add the function
                                               ! to cost, constraint or none
                                               ! values: 'cost', 'cstr', 'none'
        character(LEN=10)  ::   FUNID          ! link to function
                                               ! 'node_disp': nodal displacement


        ! ------------------------------------------------------- initialise ---
        W_COST=0.0;W_CONSTR=0.0;FLAG_COST=.false.;FLAG_CONSTR=.false.


        ! ------------------------------------------------------------ FLAG_DISP
        ! node displacement
        ! (cost_node_disp in opt_cost module)
        FUNID='node_disp';     ADDTO='cost';     WEIGHT=1.0_8;
        !NODE_DISP=0; ! <--- the passing interface for the option has not been implemented yet!!
        call cost_utl_allocate_flags_and_weights(FUNID,ADDTO,WEIGHT, &
                            & FLAG_COST,FLAG_CONSTR,W_COST,W_CONSTR)

        ! structural mass
        ! (cost_total_mass in opt_cost module)
        FUNID='mass_tot';     ADDTO='cstr';     WEIGHT=1.0_8;
        call cost_utl_allocate_flags_and_weights(FUNID,ADDTO,WEIGHT, &
                            & FLAG_COST,FLAG_CONSTR,W_COST,W_CONSTR)
        ! ----------------------------------------------------------------------


    end subroutine opt_input_cost


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_update_shared_input
!
!-> Description:
!    The routine perturbs one of the shared variables in the module.
!
!-> Remarks:
!    In the current input setup, the values of the shared variables is set in
!    input_setup. This function can be called in a main after input_setup to
!    perturb the design.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine opt_update_shared_input( NOPT )

    integer :: NOPT                        ! value in XSH to allocate



    call opt_unpack_DESIGN_SHARED(NOPT)
    print *, 'unpacked: ',XSH(:,NOPT)


 end subroutine opt_update_shared_input



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_set_FLAG_DESIGN_SHARED
!
!-> Description:
!    This routine sets the FLAG_DESIGN_SHARED logical array.
!    Each element of FLAG_DESIGN_SHARED contains a flag associated to one of the
!    shared design variables contained in the input module.
!    If the flag value is .true., the design variable associated will be used in
!    the optimisation process.
!
!-> Remark:
!    During the optimisation execution, the FLAG_BeamLength1, FLAG_Beam... etc
!    variables are set to .false./.zero. if the sensitivity analysis is performed
!    usinf FDs.
!    FLAG_DESIGN_SHARED keeps, therefore, memory of the design parameters chosen
!    for the analysis.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine opt_pack_FLAG_DESIGN_SHARED

    ! ------------------------------------------------- General Design Variables
    ! scalars
    FLAG_DESIGN_SHARED(1) = FLAG_BeamLength1
    FLAG_DESIGN_SHARED(2) = FLAG_BeamLength2
    FLAG_DESIGN_SHARED(3) = FLAG_ThetaRoot
    FLAG_DESIGN_SHARED(4) = FLAG_ThetaTip
    FLAG_DESIGN_SHARED(5) = FLAG_TipMass
    FLAG_DESIGN_SHARED(6) = FLAG_TipMassY
    FLAG_DESIGN_SHARED(7) = FLAG_TipMassZ
    FLAG_DESIGN_SHARED(8) = FLAG_Omega
    ! vectors
    FLAG_DESIGN_SHARED( 9:11) = FLAG_ExtForce
    FLAG_DESIGN_SHARED(12:14) = FLAG_ExtMomnt

    ! ------------------------------------------- Constant Beam Design Variables
    ! scalars
    FLAG_DESIGN_SHARED(15)=FLAG_SectWidth
    FLAG_DESIGN_SHARED(16)=FLAG_SectHeight
    ! matrices
    FLAG_DESIGN_SHARED(17:17+35)=reshape(FLAG_BeamStiffness,(/36/))
    FLAG_DESIGN_SHARED(53:53+35)=reshape(     FLAG_BeamMass,(/36/))

end subroutine opt_pack_FLAG_DESIGN_SHARED



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_assigns_shared_FLAGS
!
!-> Description:
!    The routine assignes a logical value .true./.false. to all the
!    FLAG_BeamLength1, FLAG_Beam... etc variables.
!    The input value is a logical array which follows the same connectivity rules
!    used to assign FLAG_DESIGN_SHARED.
!
!-> Remark:
!    this routine allows to perturb each single variable alone during the FD
!    sensitivity analysis or to update all the variables when a new design point
!    has to be explored.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine opt_assigns_shared_FLAGS( FLAGS_vector )

    logical, intent(in) :: FLAGS_vector(size( FLAG_DESIGN_SHARED ))
    ! ------------------------------------------------- General Design Variables
    ! scalars
    FLAG_BeamLength1   = FLAGS_vector(1)
    FLAG_BeamLength2   = FLAGS_vector(2)
    FLAG_ThetaRoot     = FLAGS_vector(3)
    FLAG_ThetaTip      = FLAGS_vector(4)
    FLAG_TipMass       = FLAGS_vector(5)
    FLAG_TipMassY      = FLAGS_vector(6)
    FLAG_TipMassZ      = FLAGS_vector(7)
    FLAG_Omega         = FLAGS_vector(8)
    ! vectors
    FLAG_ExtForce      = FLAGS_vector( 9:11)
    FLAG_ExtMomnt      = FLAGS_vector(12:14)

    ! ------------------------------------------- Constant Beam Design Variables
    ! scalars
    FLAG_SectWidth  = FLAGS_vector(15)
    FLAG_SectHeight = FLAGS_vector(16)
    ! matrices
    FLAG_BeamStiffness = reshape( FLAGS_vector(17:17+35), (/6,6/))
    FLAG_BeamMass      = reshape( FLAGS_vector(53:53+35), (/6,6/))

end subroutine opt_assigns_shared_FLAGS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_set_shared_FLAGS_to_false
!
!-> Description:
!    The routine assignes a .false. value to all the FLAG_BeamLength1, FLAG_Beam...
!    etc variables.
!
!-> Remark:
!    this routine allows to ensure that during the FD sensitivity analysis only
!    one design variable will be perturbed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine opt_set_shared_FLAGS_to_false

    logical :: FALSE_vector(size( FLAG_DESIGN_SHARED ))=.false.
    call opt_assigns_shared_FLAGS( FALSE_vector )

end subroutine opt_set_shared_FLAGS_to_false


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_pack_DESIGN_SHARED
!
!-> Description:
!    This routine packs all the shared design variables into the
!    multidimensional array XSH.
!    As XSH is also used to keep track of all the designs explored, the iteration
!    of the optimisation process is also required (NOPT).
!
!-> Remark:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine opt_pack_DESIGN_SHARED(NOPT)

   integer, optional :: NOPT            ! Number of iteration for the optimisation

   ! Shared variables from module input
   real(8) :: BeamLength1, BeamLength2     ! Beam defining lengths.
   real(8) :: BeamStiffness(6,6)           ! Beam element stiffness matrix (assumed constant).
   real(8) :: BeamMass(6,6)                ! Beam element mass matrix (assumed constant).
   real(8) :: ExtForce(3)                  ! Applied forces at the tip.
   real(8) :: ExtMomnt(3)                  ! Applied moments at the tip.
   real(8) :: SectWidth,SectHeight         ! Height and width of the cross section.
   real(8) :: ThetaRoot                    ! Pretwist angle at root.
   real(8) :: ThetaTip                     ! Pretwist angle at tip.
   real(8) :: TipMass                      ! Mass at the beam tip.
   real(8) :: TipMassY                     ! Y offset of the tip mass.
   real(8) :: TipMassZ                     ! Z offset of the tip mass.
   real(8) :: Omega                        ! Frequency of oscillatory motions.


    ! -------------------------------------------- Read values from input module
    call read_shared_input( BeamLength1, BeamLength2,  &
                   & BeamStiffness, BeamMass,          &
                   & ExtForce, ExtMomnt,               &
                   & SectWidth, SectHeight,            &
                   & ThetaRoot, ThetaTip,              &
                   & TipMass, TipMassY, TipMassZ,      &
                   & Omega                             )

    ! ------------------------------------------------- General Design Variables
    ! scalars
    XSH(1,NOPT) = BeamLength1
    XSH(2,NOPT) = BeamLength2
    XSH(3,NOPT) = ThetaRoot
    XSH(4,NOPT) = ThetaTip
    XSH(5,NOPT) = TipMass
    XSH(6,NOPT) = TipMassY
    XSH(7,NOPT) = TipMassZ
    XSH(8,NOPT) = Omega
    ! vectors
    XSH( 9:11,NOPT) = ExtForce
    XSH(12:14,NOPT) = ExtMomnt

    ! ------------------------------------------- Constant Beam Design Variables
    ! scalars
    XSH(15,NOPT)=SectWidth
    XSH(16,NOPT)=SectHeight
    ! matrices
    XSH(17:17+35,NOPT)=reshape(BeamStiffness,(/36/))
    XSH(53:53+35,NOPT)=reshape(     BeamMass,(/36/))

end subroutine opt_pack_DESIGN_SHARED


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_unpack_DESIGN_SHARED
!
!-> Description:
!    This routine unpacks all the shared design variables from the
!    multidimensional array XSH.
!    As XSH is also used to keep track of all the designs explored, the iteration
!    of the optimisation process is also required.
!
!-> Remark:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine opt_unpack_DESIGN_SHARED(NOPT)

   integer, optional :: NOPT               ! Number of iteration for the optimisation

   ! Shared variables from module input
   real(8) :: BeamLength1, BeamLength2     ! Beam defining lengths.
   real(8) :: BeamStiffness(6,6)           ! Beam element stiffness matrix (assumed constant).
   real(8) :: BeamMass(6,6)                ! Beam element mass matrix (assumed constant).
   real(8) :: ExtForce(3)                  ! Applied forces at the tip.
   real(8) :: ExtMomnt(3)                  ! Applied moments at the tip.
   real(8) :: SectWidth,SectHeight         ! Height and width of the cross section.
   real(8) :: ThetaRoot                    ! Pretwist angle at root.
   real(8) :: ThetaTip                     ! Pretwist angle at tip.
   real(8) :: TipMass                      ! Mass at the beam tip.
   real(8) :: TipMassY                     ! Y offset of the tip mass.
   real(8) :: TipMassZ                     ! Z offset of the tip mass.
   real(8) :: Omega                        ! Frequency of oscillatory motions.

    ! ------------------------------------------------- General Design Variables
    ! scalars
    BeamLength1 = XSH(1,NOPT)
    BeamLength2 = XSH(2,NOPT)
    ThetaRoot   = XSH(3,NOPT)
    ThetaTip    = XSH(4,NOPT)
    TipMass     = XSH(5,NOPT)
    TipMassY    = XSH(6,NOPT)
    TipMassZ    = XSH(7,NOPT)
    Omega       = XSH(8,NOPT)
    ! vectors
    ExtForce    = XSH( 9:11,NOPT)
    ExtMomnt    = XSH(12:14,NOPT)

    ! ------------------------------------------- Constant Beam Design Variables
    ! scalars
    SectWidth  = XSH(15,NOPT)
    SectHeight = XSH(16,NOPT)
    ! matrices
    BeamStiffness = reshape( XSH(17:17+35,NOPT), (/6,6/))
    BeamMass      = reshape( XSH(53:53+35,NOPT), (/6,6/))


    ! ------------------------------------Update shared variable in input module
    call update_shared_input( BeamLength1, BeamLength2,&
                   & BeamStiffness, BeamMass,          &
                   & ExtForce, ExtMomnt,               &
                   & SectWidth, SectHeight,            &
                   & ThetaRoot, ThetaTip,              &
                   & TipMass, TipMassY, TipMassZ,      &
                   & Omega                             )

end subroutine opt_unpack_DESIGN_SHARED




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_print_FLAG_DESIGN_SHARED
!
!-> Description:
!    prints FLAG_DESIGN_SHARED
!
!-> Remark:
!   only for testing purposed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine opt_print_FLAG_DESIGN_SHARED

     print *, 'FLAG_DESIGN_SHARED(1:8) '  , FLAG_DESIGN_SHARED(1:8)
     print *, 'FLAG_DESIGN_SHARED(9:14) ' , FLAG_DESIGN_SHARED(9:14)
     print *, 'FLAG_DESIGN_SHARED(15:16) ', FLAG_DESIGN_SHARED(15:16)
     print *, 'FLAG_DESIGN_SHARED(17:52) ', FLAG_DESIGN_SHARED(17:52)
     print *, 'FLAG_DESIGN_SHARED(53:) '  , FLAG_DESIGN_SHARED(53:)

 end subroutine opt_print_FLAG_DESIGN_SHARED




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine opt_print_XSH
!
!-> Description:
!    prints XSH
!
!-> Remark:
!   only for testing purposed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine opt_print_XSH(NOPT)

     integer, optional :: NOPT               ! Number of iteration for the optimisation

     print *, 'NOPT = ', NOPT
     print *, 'XSH(1:8) '  , XSH(1:8,NOPT)
     print *, 'XSH(9:14) ' , XSH(9:14,NOPT)
     print *, 'XSH(15:16) ', XSH(15:16,NOPT)
     print *, 'XSH(17:52) ', XSH(17:52,NOPT)
     print *, 'XSH(53:) '  , XSH(53:,NOPT)

 end subroutine opt_print_XSH


end module opt_input






