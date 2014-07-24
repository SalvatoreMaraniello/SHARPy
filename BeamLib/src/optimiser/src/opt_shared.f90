!-> Module.- opt_shared Salvatore Maraniello 14/07/2014
!
!-> Author: Salvatore Maraniello
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Shared Variables for Optimisaiton
!
!-> Subroutines.-
!
! - generate_report: prints parameters values
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_shared

 implicit none

 real(8), parameter :: perturb_ratio = 1e-8_8                  ! 15 digit precision
 real(8), parameter :: perturb_min   = &
                       & 1e+8_8*epsilon(perturb_ratio) ! about 1e-8: below perturb_min, the value is
                                                       ! considered eqaul to zero and the perturbation is added


 integer, parameter :: NCOSTFUNS=2 ! total number of cost/constraints functions implemented

end module opt_shared


