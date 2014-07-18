!-> Module.- opt_fd
!
!-> Author: Salvatore Maraniello 17 july 2014
!
!-> Language: FORTRAN90, Free format.
!
!-> Description: Module to handle the FD gradient evaluation
!
!-> Subroutines.-
!   - fd_main:
!
! -> Remark: this module is placed at main level to allow the call of:
!       a. input module: this is necessary to perturbe each design variable
!          individually;
!       b. fwd_main: to allow the fwd code execution;

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module opt_fd

 use input       ! to perturb the shared design input
 use opt_input   ! to access the FLAGS for design variables and their values XSH
 use fwd_main    ! to run the fwd problem
 use opt_perturb

 implicit none

 ! variables visible from opt_input
 ! real    ::  XSH(size(FLAG_DESIGN_SHARED), 0:NOPT_MAX )
 ! logical ::  FLAG_DESIGN_SHARED(8+6+6*6+6*6+2)

 contains

    subroutine fd_main(NOPT,fdmode )


    character(len=3), intent(in) :: fdmode   ! finite differences method
    integer                      :: NOPT     ! number of iteration for the optimisation - required to understand which column of XSH to modify
    real(8) :: DXSH( size(XSH,1) )           ! deltas for shared design variables


    ! ---------------------- Compute Deltas for Design current design (NOPT) ---
    where (FLAG_DESIGN_SHARED .eqv. .true.)
        DXSH = fperturb_delta_1d_array( XSH(:,NOPT) )
    elsewhere
        DXSH = 0.0_8
    end where


    ! ------------------------------ Testing: result without where statement ---
    !call opt_print_XSH(NOPT)
    !call opt_print_FLAG_DESIGN_SHARED
    !print *, 'Delta only on target design variables: ', DXSH
    !DXSH = fperturb_delta_1d_array( XSH(:,NOPT) )
    !print *, 'Delta without over all array: ', DXSH
    !stop


    ! ------------------------------------ Perturb each variable at the time ---




    end subroutine fd_main

 end module opt_fd


