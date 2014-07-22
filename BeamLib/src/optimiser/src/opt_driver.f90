!-> Module.- opt_driver Salvatore Maraniello 21/07/2014
!
!-> Author: Salvatore Maraniello
!
!-> Language: FORTRAN90, Free format.
!
!-> Description: collects drivers for the optimisation
!
!-> Subroutines.-
!
!   - Cost Functions:
!     a.
!
!
!-> Remarks:
!   a.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_driver


 implicit none

 contains


!-------------------------------------------------------------------------------
!
!------------------------
    !
    subroutine simple_driver(XSH,COST,CONSTR,DCDXSH,DCONDXSH,NOPT,CONN_XSH,CONN_CONSTR)

        real(8), intent(inout) :: XSH(:,0:) ! SHARED DESIGN VARIABLES

        real(8), intent(in) :: COST(0:)                 ! cost function
        real(8), intent(in) :: CONSTR(:,0:)             ! constrain vector
        real(8), intent(in) :: DCDXSH  (:,0:)  ! gradient of cost in respect to shared design
        real(8), intent(in) :: DCONDXSH(:,:,0:)! gradient of constrain in respect to design
        integer, intent(in) :: NOPT

        integer :: CONN_CONSTR(:), CONN_XSH(:), nn, ii   ! connectivity matrix for contrains and design variables array

        real(8)  :: DXSH( size(CONN_XSH) ) ! local dx
        DXSH=0.0




        DXSH=steepest_descent(COST(NOPT),DCDXSH(:,NOPT))
        print *, 'step computed', DXSH



        XSH(:,NOPT+1) = XSH(:,NOPT)
        do nn=1,size(CONN_XSH)
            ii = CONN_XSH(nn)
            print *, 'XSH(', ii,')=',XSH(ii,NOPT)
            XSH( ii, NOPT+1 ) = XSH( ii, NOPT ) + DXSH( nn )
            print *, 'XSH(', ii,')=',XSH(ii,NOPT+1)
        end do



    end subroutine simple_driver


!-------------------------------------------------------------------------------
! steepest descent:
!------------------------
    ! simplified version of steepest_descent algorythm.
    ! The descent direction is chosen to be opposite to the gradient, while the
    ! step is set such that, assuming the function being linear in respect to
    ! all the design parameter, each change in each design parameter gives a
    ! factor reduction of the cost.

    function steepest_descent(cost,gvec,factor)

        real(8), optional   :: factor ! factor for determining step length
        real(8), intent(in) :: cost   ! cost at current design
        real(8), intent(in) :: gvec(:)! gradient vector

        real(8) :: steepest_descent(size(gvec))  ! step for next iteration

        real(8) :: cf
        integer :: ii

        steepest_descent=0.0_8; cf = -1.e-2_8 * cost
        if (present(factor)) then
            cf = -factor*cost
        end if

        where (abs(gvec)<1e-6_8)
            steepest_descent=0.0_8
        elsewhere
            steepest_descent=cf/gvec
        end where

        print *, 'gradeint in input', gvec
        print *, 'descent step', steepest_descent

    end function steepest_descent





end module opt_driver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module opt_driver_test

 use opt_driver
 implicit none

 contains


    subroutine simple_driver_test

        real(8) :: factor ! factor for determining step length
        real(8) :: cost   ! cost at current design
        real(8) :: gvec(5)! gradient vector
        real(8) :: dx(5)  ! step for next iteration


        cost=1.0
        gvec = (/ 1.0, 1.0, 1.0, 1., 0.0001 /)

        dx=steepest_descent(cost,gvec)

        print*, dx


    end subroutine simple_driver_test

end module opt_driver_test







