!->Module lib_shape. Salvatore Maraniello. 11/Aug/2014
!
!->Description.-
!
!  This module includes tools for varying the beam properties along the span.
!
!
!->Subroutines:
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lib_shape


implicit none

contains



! subroutine polyshape
! ------------------------------------------------------------------------------
    ! Varies properties according to a polynomial law:
    !   p(eta) = sum( P_i eta**i ), i=0:N
    ! the function is constrained to consider only polynomial up to the 5th
    ! order
    ! --------------------------------------------------------------------------
    subroutine polyshape(Pvec,N,eta,p)

        integer, intent(in)  :: N         ! Number of terms in the polynomial expantion
        real(8), intent(in)  :: eta       ! non-dimensional coordinate for evaluation
        real(8), intent(in)  :: Pvec(0:N) ! Coefficients of polynomial expantion
        real(8), intent(out) :: p         ! Output value

        integer :: ii ! counter


        if (N>5) then
            stop 'Polynomial evaluation may be inaccurate. See polyshape [lib_shape]'
        end if

        p = 0.0_8
        do ii=0,N
            p = p + Pvec(ii) * (eta**ii)
        end do

    end subroutine polyshape

end module lib_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lib_shape_test

use lib_shape

implicit none

contains

    subroutine polyshape_test

        real(8)  :: eta(1:4) ! non-dimensional coordinate for evaluation
        real(8)  :: Pvec(0:5) ! Coefficients of polynomial expantion
        real(8)  :: p         ! Output value

        integer :: ii ! counter

        Pvec=0.0_8

        ! allocate eta
        do ii=1,4
            eta(ii)=real(ii,8)
        end do

        print *, '0 order:'
        Pvec(0)=10.0_8
        do ii=1,4
            call polyshape(Pvec,5,eta(ii),p)
            print '(F20.5)', p
        end do

        print *, 'I order:'
        Pvec=0.0_8; Pvec(1)=10.0_8
        do ii=1,4
            call polyshape(Pvec,5,eta(ii),p)
            print '(F20.5)', p
        end do

        print *, 'II order:'
        Pvec=0.0_8; Pvec(2)=10.0_8
        do ii=1,4
            call polyshape(Pvec,5,eta(ii),p)
            print '(F20.5)', p
        end do

        print *, 'III order:'
        Pvec=0.0_8; Pvec(3)=10.0_8
        do ii=1,4
            call polyshape(Pvec,5,eta(ii),p)
            print '(F20.5)', p
        end do

        print *, 'IV order:'
        Pvec=0.0_8; Pvec(4)=10.0_8
        do ii=1,4
            call polyshape(Pvec,5,eta(ii),p)
            print '(F20.5)', p
        end do

        print *, 'V order:'
        Pvec=0.0_8; Pvec(5)=10.0_8
        do ii=1,4
            call polyshape(Pvec,5,eta(ii),p)
            print '(F20.5)', p
        end do

    end subroutine polyshape_test

end module lib_shape_test
