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


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine opt_setup(gradmode,solmode)

        character(len=3), intent(out) :: gradmode ! gradient method
        character(len=3), intent(out) :: solmode  ! solution mode
        ! ----------------------------------------------------------------------



        gradmode='FDF' ! gradient method: DIR: direct
                       !                  ADJ: adjointfinite differences
                       !                  FDF: Finite Differences
        solmode ='SNS' ! solution mode:   FWD: forward
                       !                  OPT: optimisation
                       !                  SNS: sensitivity analysis

        ! ----------------------------------------------------------------------
        ! Design Parameters: shared variables
        FLAG_BeamLength1 = .true.
        FLAG_TipMassY =.true.
        FLAG_ExtForce(3)=.true.
        FLAG_BeamMass(1,1)=.true.
        FLAG_BeamMass(2,4)=.true.
        FLAG_ExtMomnt(2)=.true.
        FLAG_BeamStiffness(6,6)=.true.


        ! assign all the flags to FLAG_DESIGN_SHARED
        call opt_pack_FLAG_DESIGN_SHARED
        !print *, 'FLAG_DESIGN_SHARED(1:8) '  , FLAG_DESIGN_SHARED(1:8)
        !print *, 'FLAG_DESIGN_SHARED(9:14) ' , FLAG_DESIGN_SHARED(9:14)
        !print *, 'FLAG_DESIGN_SHARED(15:16) ', FLAG_DESIGN_SHARED(15:16)
        !print *, 'FLAG_DESIGN_SHARED(17:52) ', FLAG_DESIGN_SHARED(17:52)
        !print *, 'FLAG_DESIGN_SHARED(53:) '  , FLAG_DESIGN_SHARED(53:)
        call opt_set_shared_FLAGS_to_false

    end subroutine opt_setup


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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine opt_update_shared_input( FLAGS_vector )

   logical, intent(in) :: FLAGS_vector(size( FLAG_DESIGN_SHARED ))

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

   ! imports input module shared variables
   call read_shared_input( BeamLength1, BeamLength2,   &
                   & BeamStiffness, BeamMass,          &
                   & ExtForce, ExtMomnt,               &
                   & SectWidth, SectHeight,            &
                   & ThetaRoot, ThetaTip,              &
                   & TipMass, TipMassY, TipMassZ,      &
                   & Omega                             )


    call opt_assigns_shared_FLAGS( FLAGS_vector )

    ! Flags to Perturb Shared Variables
    if ( FLAG_BeamLength1 .eqv. .true. ) then
        d_BeamLength1 = fperturb_delta_scalar(BeamLength1)
        BeamLength1 = BeamLength1 + d_BeamLength1
    end if

    ! Flags to Perturb Shared Variables
    if ( FLAG_BeamLength2 .eqv. .true. ) then
        d_BeamLength2 = fperturb_delta_scalar(BeamLength2)
        BeamLength2 = BeamLength2 + d_BeamLength2
    end if

!    if ( FLAG_BeamLength2=.true. ) then
!    if ( FLAG_SectWidth=.true. ) then
!    if ( FLAG_SectHeight=.true. ) then
!    if ( FLAG_ThetaRoot=.true. ) then
!    if ( FLAG_ThetaTip =.true. ) then
!    if ( FLAG_TipMass =.true. ) then
!    if ( FLAG_TipMassY =.true. ) then
!    if ( FLAG_TipMassZ=.true. ) then
!    if ( FLAG_Omega=.true. ) then


     where ( FLAG_BeamMass .eqv. .true. )
         d_BeamMass = fperturb_delta_2d_array(BeamMass)
         BeamMass = BeamMass+d_BeamMass
     end where



!    if ( FLAG_BeamStiffness(6,6)=.true. ) then
!    if ( FLAG_BeamMass(6,6)=.true. ) then
!    if ( FLAG_ExtForce(3)=.true. ) then
!    if ( FLAG_ExtMomnt(3)=.true. ) then


    call update_shared_input( BeamLength1, BeamLength2, &
                   & BeamStiffness, BeamMass,          &
                   & ExtForce, ExtMomnt,               &
                   & SectWidth, SectHeight,            &
                   & ThetaRoot, ThetaTip,              &
                   & TipMass, TipMassY, TipMassZ,      &
                   & Omega                             )

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


end module opt_input
