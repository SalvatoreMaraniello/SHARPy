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

contains

    subroutine generate_report

        print *, 'Perturbation ratio: ',           perturb_ratio
        print *, 'Zero:',          perturb_min

    end subroutine generate_report

end module opt_shared

    ! Perturb BeamLength
    ! 1. if the beam is pretwisted, the root and tip pretwist of the perturbed
    !    beam are the same as those of the inperturbed bone;
    !subroutine perturb_span(ID)
    !    integer, optional :: ID=1 ! Beam to perturb; defalut BeamLength1
    !    print *, 'lala'
    !end subroutine BeamLength

    !subroutine perturb_node(ii)
    !    integer :: ii ! node number
    !    ! a. is a master?
    !    ! b. if yes, apply a displacement
    !    ! c. apply it to all the slaves
    !
    !end subroutine perturb_node
