!-> Program.- MAIN - 17july2014 Updated: july 2014
!
!-> Author.- Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!   Runs generic test for the optimisation libraries.
!
! -> Developer Notes:
!    Test functions/modules are in the same files as the modules to test under
!    the module name xxx_test
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program lib_test

 use lib_array_test

 implicit none


 character(len=11) :: testmodule ! module to test


 testmodule='lib_array'
 print *, 'START TESTING: ', testmodule


 select case (trim(testmodule))

    case ('lib_array')
        call array_cond_alloc_test

    case default
        print *, 'TEST CASES NOT IMPLEMENTED FOR: ', testmodule

 end select

 print *, 'END TESTING: ', testmodule

end program lib_test











