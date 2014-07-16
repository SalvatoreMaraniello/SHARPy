!-> Program.- MAIN - 15july2014 Updated: july 2014
!
!-> Author.- Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!   Runs generic test of optimisation modules.
!
! -> Developer Notes:
!    Test functions/modules are in the same files as the modules to test under
!    the module name xxx_test
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program opt_test

 use opt_shared
 use opt_perturb_test

 implicit none

 character(len=11) :: testmodule ! module to test

 testmodule='opt_perturb'
 print *, 'START TESTING: ', testmodule

 select case (trim(testmodule))

    case ('opt_perturb')
        call perturb_scalar_test
        call perturb_array_test

    case ('opt_shared')
        call generate_report


    case default
        print *, 'TEST CASES NOT IMPLEMENTED FOR: ', testmodule

 end select

 print *, 'END TESTING: ', testmodule

end program opt_test

