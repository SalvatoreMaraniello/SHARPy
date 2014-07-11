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
    implicit none

    ! Shared Variables - won;t leave the module if private



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
        solmode ='FWD' ! solution mode:   FWD: forward
                       !                  OPT: optimisation
                       !                  SNS: sensitivity analysis

    end subroutine opt_setup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the Design Parameters
!
!    subroutine opt_design()
!        logical :: OPTL_do=.false. ! optimise for Beam Length
!
!    end subroutine







end module opt_input
