!-> Module. - fwd_wrapper - 25/jul/2014
!
!-> Author.- Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN95, Free format.
!
!-> Description.-
!   - This module creates a high level wrapper for the beam code solver.
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fwd_wrapper

use fwd_main
use input


 implicit none


contains

! subroutine launch_fwd_solver_test()
! -----------------------------------
    ! Test for wrapper to launch the code from fortran.
    ! Most of the input are stored inside the input module.
    ! this is run at the beginning of the execution, to allocate the default
    ! variables in the input module.
    ! These variables will then be overwritten if they are input of this
    ! subroutine.
    ! this process has been proved to work in the opt_main module.
    subroutine launch_fwd_solver_test()






 end subroutine launch_fwd_solver_test





end module fwd_wrapper
