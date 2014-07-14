!-> Module.- opt_FD_perturb Salvatore Maraniello 14/07/2014
!
!-> Author: Salvatore Maraniello
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Perturbs a design variable
!
!-> Subroutines.-
!
! - fd_perturb_span: increases beam length
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_fd_perturb


 implicit none


 real(8) :: fd_delta=0.00010000 ! public, can be accessed from outside


contains

    ! Perturb BeamLength
    ! 1. if the beam is pretwisted, the root and tip pretwist of the perturbed
    !    beam are the same as those of the inperturbed bone;
    !subroutine perturb_span(ID)
    !    integer, optional :: ID=1 ! Beam to perturb; defalut BeamLength1
    !    print *, 'lala'
    !end subroutine BeamLength

    subroutine perturb_node()
        ! a. is a master?
        ! b. if yes, apply a displacement
        ! c. apply it to all the slaves

    end subroutine perturb_node


end module
