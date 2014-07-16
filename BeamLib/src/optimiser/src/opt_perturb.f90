!-> Module.- opt_perturb
!            opt_perturb_test
!
!-> Author: Salvatore Maraniello 15 july 2014
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  opt_perturb: Perturb design Variables for Optimisation
!  opt_perturb_test: testing functions
!
!-> Subroutines.-
!   - perturb_*: the subroutines ??
!

! -> Function.-
!    - fperturb_delta_*: these are implemented in function format to allow them to be
!                  used within a 'where' statement
!       - fperturb_delta_scalar: returns the delta to perturb a scalar
!       - fperturb_delta_1d_array: returns the delta to perturb a 1d array
!       - fperturb_delta_2d_array: returns the delta to perturb a 2d array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_perturb

 use opt_shared
 implicit none

contains


    ! ---------------------------------------------------------- fperturb_delta_scalar
    ! determines the delta to apply to a scalar x.
    function fperturb_delta_scalar(x)

        real(8), intent(in) ::  x               ! scalar to perturb
        real(8)             ::  fperturb_delta_scalar ! delta given to x

        if (abs(x)>perturb_min) then
            fperturb_delta_scalar=perturb_ratio*abs(x)
        else
            fperturb_delta_scalar=perturb_ratio
        end if

    end function fperturb_delta_scalar



    ! -------------------------------------------------------- fperturb_delta_1d_array
    ! determines the delta to apply to a 1d array xv.
    function fperturb_delta_1d_array(xv)

        real(8), intent(in)  ::                       xv(:) ! array to perturb
        real(8)              :: fperturb_delta_1d_array(size(xv)) ! delta given to xv

        where (abs(xv)>perturb_min)
            fperturb_delta_1d_array=perturb_ratio*abs(xv)
        elsewhere
            fperturb_delta_1d_array=perturb_ratio
        end where

    end function fperturb_delta_1d_array


    ! -------------------------------------------------------- fperturb_delta_2d_array
    ! determines the delta to apply to a 2d array xv.
    !
    ! Remark:
    ! old piece of code commented out - discarded for performance
    function fperturb_delta_2d_array(X)

        real(8), intent(in) ::                           X(:,:) ! array to perturb
        real(8)             :: fperturb_delta_2d_array( &
                                   & size(X(:,1)),size(X(1,:))) ! delta given to X

        ! performance: 0.055 sec for 3000x2000 array
        where (abs(X)>perturb_min)
            fperturb_delta_2d_array=perturb_ratio*abs(X)
        elsewhere
            fperturb_delta_2d_array=perturb_ratio
        end where

        ! old code: discarded for low performance
        ! performance: 0.12 sec for a 3000x2000 array
        ! internal variables
        !integer(8)             ::  sh(2)   ! shape of input array
        !integer                ::   N      ! Number of elements in the array
        !real(8),allocatable    ::  xv(:)   ! reshaped 1d array corresponding to x
        !real(8),allocatable    :: dxv(:)   ! reshaped 1d array corresponding to dx
        !sh=shape(X); N=product(sh)
        !allocate(  xv(N), dxv(N) )
        !xv=reshape(X,(/N/)); dxv=0.0_8
        !dxv= fperturb_delta_1d_array(xv)
        !DX=reshape(dxv,sh)

    end function fperturb_delta_2d_array



    !subroutine perturb_




end module opt_perturb



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module opt_perturb_test

 !use opt_shared
 use opt_perturb
 use lib_perf
 implicit none

contains

    ! --------------------------------------------------------------------------
    subroutine perturb_scalar_test

        real(8) :: x, dx

        print *, '----------------------- fperturb_delta_scalar test'

        x=+1e-10_8
        dx = fperturb_delta_scalar(x)
        x = x+dx
        print *, 'Perturb 1e-10 [x, dx]: ', x, dx

        x=+1e-4_8
        dx = fperturb_delta_scalar(x)
        x = x+dx
        print *, 'Perturb 1e-4 [x, dx]: ', x, dx

        x=-1e+6_8
        dx = fperturb_delta_scalar(x)
        x = x+dx
        print *, 'Perturb 1e+6 [x, dx]: ', x, dx

    end subroutine perturb_scalar_test


    ! --------------------------------------------------------------------------
    subroutine perturb_array_test

        real(8) :: A(3),     DA(3)
        real(8) :: B(3,2),   DB(3,2)
        real(8) :: C(3000,2000),   DC(3000,2000)

        print *, '----------------------- fperturb_delta_1d_array test'
        A=4e5; A(1)=-4.0; A(2)=1e-14
        print *, '1D array:'
        print *, 'initial:', A
        DA = fperturb_delta_1d_array(A)
        print *, 'perturbed', A+DA
        print *, 'delta', DA

        print*, ' '
        print *, '----------------------- fperturb_delta_2d_array test'
        B=3.0; B(:,2)=(/ 4e5, -4.0, 1e-14 /); B(2,1)=-1e7
        print *, '2D array:'
        print *, 'initial:', B
        call tic
        DB = fperturb_delta_2d_array(B)
        call toc
        print *, 'perturbed', B+DB
        print *, 'delta', DB

        print*, ' '
        C=3.0; C(:,1)=4e5; C(:,2)=-1e7; C(:,3)=-1e-14
        print *, '2D array of size:', size(C)
        print *, 'initial - C(2,1:4):', C(2,1:4)
        call tic
        DC = fperturb_delta_2d_array(C)
        C = C+DC
        call toc
        print *, 'perturbed  - C(2,1:4):', C(2,1:4)
        print *, 'delta', DC(2,1:4)

    end subroutine perturb_array_test


end module opt_perturb_test




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
