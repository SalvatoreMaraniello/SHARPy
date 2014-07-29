!-> Program.- MAIN - 25july2014 Updated: july 2014
!
!-> Author.- Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!   test run for opt_routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program opt_routine_test

 use opt_routine

 implicit none


 integer :: NumElems     ! Number of elements
 integer :: NumNodes     ! Number of nodes in the model.

 real(8), pointer, dimension(:) :: pCOST
 real(8) :: W_COST(2)    ! <--- example of known size array, 2 is NCOSTFUNS into opt_shared


 ! Define Input:
 NumElems = 10
 NumNodes = 3

 call opt_main(NumElems, NumNodes, pCOST, W_COST)

 print '(F10.5)', W_COST

end program opt_routine_test

