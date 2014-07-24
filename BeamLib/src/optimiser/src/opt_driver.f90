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
!-> References:
! a. Nocedal, J., Wright, S. J., Numerical Optimization, 2006
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_driver

 use interface_lapack
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

        real(8) :: DXSH(size(CONN_XSH)), xsh_old(size(CONN_XSH)), xsh_now(size(CONN_XSH)) ! local dx
        DXSH=0.0_8; xsh_old=0.0_8; xsh_now=0.0_8

        !if (NOPT>=0) then
            DXSH=steepest_descent(COST(NOPT),DCDXSH(:,NOPT))
        !else
        !    ! extract values of design variables
        !    do nn=1,size(CONN_XSH)
        !        ii = CONN_XSH(nn)
        !        xsh_now(nn)=XSH(ii,NOPT  );
        !        xsh_old(nn)=XSH(ii,NOPT-1);
        !    end do
        !    ! quasi-newton method is wrong! Hessian computation is not correct.
        !    DXSH=quasi_newton(DCDXSH(:,NOPT),DCDXSH(:,NOPT-1),xsh_now,xsh_old)
        !end if

        XSH(:,NOPT+1) = XSH(:,NOPT)
        do nn=1,size(CONN_XSH)
            ii = CONN_XSH(nn)
            !print *, 'XSH(', ii,')=',XSH(ii,NOPT)
            XSH( ii, NOPT+1 ) = XSH( ii, NOPT ) + DXSH( nn )
            !print *, 'XSH(', ii,')=',XSH(ii,NOPT+1)
        end do



    end subroutine simple_driver


!-------------------------------------------------------------------------------
! steepest descent:
!------------------
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

        where (abs(gvec)<1e-10_8)
            steepest_descent=0.0_8
        elsewhere
            steepest_descent=cf/gvec
        end where

    end function steepest_descent



!-------------------------------------------------------------------------------
! quasi-newton: ROUTINE INCOMPLETE
!------------------
    ! simple quasi-newton method
    ! The Newton direction (ref.a, 6.1) is found using an approximation of the
    ! Hessian matrix. This last is recomputed by FDs of the gradient at each
    ! iteration (instead of being updated as in other methods -  e.g. BFGS)
    !
    ! IMPORTANT: the routine is WRONG:
    ! the approximate Hessian can only be used as initial guess for a BFGS or
    ! SR1 method. Even in this case,

    function quasi_newton(gvec,gvec_old,xvec,xvec_old,factor)

        real(8), optional   :: factor ! factor to correct the step length. This
                                      ! is not really rewuired as a natural length
                                      ! of 1 is associated with the method (ref.a)

        real(8), intent(in) :: gvec(:)      ! current gradient
        real(8), intent(in)    :: gvec_old(:)  ! old  gradient
        real(8), intent(in)    :: xvec(:),xvec_old(:) ! design old and current

        integer :: Nvec
        real(8) :: quasi_newton(size(gvec))  ! step for next iteration
        real(8) :: H(size(gvec), size(gvec)) ! Hessian matrix

        real(8) :: cf, dx(size(gvec)), dg(size(gvec))
        integer :: ii, jj

        integer :: LUpivot(size(gvec))  ! pivoting array for LU factorisation
        integer :: Info                 ! Info: for LU factorisation
        real(8) :: Hcp(size(gvec), size(gvec)) ! copy of hessian
        real(8) :: gcp(size(gvec))

        Nvec=size(gvec)
        quasi_newton=0.0_8; cf = 1.0_8;
        gcp=gvec ! copy of gvec necessary cause the lu functions overwrite gvec

        if (present(factor)) then
            cf = factor
        end if

        ! approximate hessian:
        H=0.0
        dg = gvec-gvec_old
        dx = xvec-xvec_old

        print *, dg
        print *, dx

        ! the approx below is wrong
        do 10 jj=1,Nvec
            do 20 ii=1,Nvec
                H(ii,jj)=dg(jj)/dx(ii)
            20 end do
        10 end do
        H = 0.5_8 * (H + transpose(H))

        print *, 'Hessian:'
        print '(F12.6)', H

        ! solve for newton direction
        ! see interface lapack module. The function called are:
        !! http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/dgetrf.html
        !CALL DGETRF( M,    N, A, LDA, IPIV, INFO )
        ! http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/dgetrs.html
        !CALL DGETRS( 'N', Nvec, 1, H, Nvec, IPIV, B, LDB, INFO )

        Hcp = H
        call lapack_lufact(Nvec,H,LUpivot,Info)

        print *, 'Hessian Factorised:'
        print '(F12.6)', H

        call lapack_luback(Nvec,H,LUpivot,gcp,Info) ! solution in gvec
        quasi_newton = gcp

        print *, 'Solution'
        print '(F12.6)', quasi_newton
        print *, 'check solution:'
        print '(F12.6)', matmul(Hcp,gvec)

    end function quasi_newton



!-------------------------------------------------------------------------------
! solve_linear_sys:
!------------------
    ! routine to solve linear system by LU factorisation using LApaCK library
    ! The routine holds on the interface_lapack library
    ! The linear system solved is:
    !   Ax = b
    !
    subroutine solve_linear_sys(A,b,x)

        real(8)  :: A(:,:), b(:)
        real(8), intent(out) :: x(:)
        integer :: N, Info, LUpivot(1:size(b))

        N=size(b)
        LUpivot=0
        x=0.0_8
        Info=0

        call lapack_lufact(N,A,LUpivot,Info)
        call lapack_luback(N,A,LUpivot,b,Info) ! solution in gvec
        x=b

    end subroutine solve_linear_sys



end module opt_driver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module opt_driver_test

 use opt_driver
 implicit none

 contains


    subroutine optimisers_test

        real(8) :: factor ! factor for determining step length
        real(8) :: cost   ! cost at current design
        real(8) :: gvec(3), gvec_old(3) ! gradient vector
        real(8) :: xvec(3), xvec_old(3) ! design vectors
        real(8) :: dx(3)  ! step for next iteration

        cost=1.0
        gvec     = (/ 7.0_8, 12.0_8, 15.0_8 /)
        gvec_old = (/ 0.0_8, 1.0_8, 2.0_8 /)

        xvec     = (/ 5.0_8, 3.0_8, 1.0_8 /)
        xvec_old = (/ 0.0_8, 0.0_8, 0.0_8 /)

        dx=steepest_descent(cost,gvec)

        print '(A20)', 'gradient in input'
        print '(F12.6)', gvec
        print *, ' '
        print '(A20)', 'stepeest descent'
        print '(F12.6)', dx

        dx=quasi_newton(gvec,gvec_old,xvec,xvec_old)
        print *, ' '
        print '(A20)', 'quasi newton'
        print '(F16.8)', dx

    end subroutine optimisers_test





    subroutine lin_solver_test
        use interface_lapack
        real(8)  :: A(3,3), Acp(3,3), Ainv(3,3), b(3), bcp(3)
        real(8) ::          x(3)
        integer :: Info

        A = transpose( reshape( (/  5.0_8 , 10.0_8 ,  2.0_8, &
                                 &  4.0_8 ,  3.0_8 ,  1.0_8, &
                                 &  4.0_8 ,  1.0_8 ,  2.0_8  /) , (/3,3/)) )
        b = (/  10.0_8 , 20.0_8 , 4.0_8  /)
        Acp=A
        bcp=b

        print *, 'Expected Solution: A x =b -> x = [9.4839 -0.3871 -16.7742]'
        call solve_linear_sys(A,b,x)
        print *, 'x ='
        print '(F16.8)', x
        print *, 'check: A*x='
        print  '(F16.8)', matmul(Acp,x)


        print *, 'Expected Solution: AT x=b -> x = [2 0 0]'

        Acp=transpose(Acp)
        A=Acp
        print '(F16.8)', A
        call solve_linear_sys(A,bcp,x)
        print *, 'x ='
        print '(F16.8)', x
        print *, 'check: Acp*x='
        print  '(F16.8)', matmul(Acp,x)

        ! -------------------------------------------- expected LU factorisation
        !L =
        ! 1.0000         0         0
        ! 0.8000    1.0000         0
        ! 0.8000    0.7143    1.0000
        !U =
        !    5.0000   10.0000    2.0000
        !         0   -7.0000    0.4000
        !         0         0   -0.8857
        !P =
        !     1     0     0
        !     0     0     1
        !     0     1     0
        ! ----------------------------------------------------------------------

    end subroutine lin_solver_test



end module opt_driver_test







