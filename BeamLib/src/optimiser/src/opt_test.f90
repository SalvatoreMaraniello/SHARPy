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
 use opt_cost_test
 use opt_cost_utl_test
 use opt_driver
 use opt_driver_test

 implicit none

 character(len=15) :: testmodule ! module to test
 !-----------------------------------------------------------------------------


 testmodule='opt_cost'       ! <--- change here!


 !-----------------------------------------------------------------------------
  print *, 'START TESTING: ', testmodule

 select case (trim(testmodule))

    case ('opt_perturb')
        call perturb_scalar_test
        call perturb_array_test

    case ('opt_shared')
        print *, 'Nothing to be done here...'

    case ('opt_cost')
        call cost_node_disp_test
        call cost_mass_total_test

    case ('opt_cost_utl')
        call cost_utl_build_connectivity_test

    case('opt_driver')
        call optimisers_test
        call lin_solver_test

    case default
        print *, 'TEST CASES NOT IMPLEMENTED FOR: ', testmodule

 end select

 print *, 'END TESTING: ', testmodule

end program opt_test

