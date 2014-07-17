!->Module lib_array. Salvatore Maraniello. 16/Jul/2014
!
!->Description.-
!
!  This module includes few tools for arrays manipulation
!
!->Subroutines:
! - array1_cond_alloc: conditional allocation for arrays of dimension 1
! - array2_cond_alloc:                           "                    2
! - array3_cond_alloc:                           "                    3                          
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lib_array

    implicit none

contains


! array1_cond_alloc(X,Nrows,IC)
!-------------------------------------------------------------------------------
    ! If the 1D array X is not allocated, the routine will allocate it.
    ! Nrows will specify the length of the array.
    ! Whether the array will be allocated or not, if setzero is equal to .true.
    ! the array will be set to zero
    ! --------------------------------------------------------------------------
    subroutine array1_cond_alloc(X,Nrows,setzero)

        real(8),    intent(inout), allocatable :: X(:)  ! 1D array
        integer,    intent(in)                 :: Nrows ! number of rows (or length)
        logical,    intent(in)   , optional    :: setzero ! real number to assign as IC to the array X if this is allocated

        if (allocated(X) .eqv. .false.) then
            allocate( X(Nrows) )
        end if

        if ( present(setzero) .and. (setzero .eqv. .true.) ) then
            X=0.0_8
        end if

    end subroutine array1_cond_alloc


! array2_cond_alloc(X,Nrows,IC)
!-------------------------------------------------------------------------------
    ! If the 2D array X is not allocated, the routine will allocate it.
    ! Nrows will specify the length of the array.
    ! Whether the array will be allocated or not, if setzero is equal to .true.
    ! the array will be set to zero
    ! --------------------------------------------------------------------------
    subroutine array2_cond_alloc(X,Nrows,Ncols,setzero)

        real(8),    intent(inout), allocatable :: X(:,:)  ! 1D array
        integer,    intent(in)                 :: Nrows, Ncols ! number of rows/columns
        logical,    intent(in)   , optional    :: setzero ! real number to assign as IC to the array X if this is allocated

        if (allocated(X) .eqv. .false.) then
            allocate( X(Nrows, Ncols) )
        end if

        if ( present(setzero) .and. (setzero .eqv. .true.) ) then
            X=0.0_8
        end if

    end subroutine array2_cond_alloc


! array3_cond_alloc(X,Nrows,IC)
!-------------------------------------------------------------------------------
    ! If the 3D array X is not allocated, the routine will allocate it.
    ! Nrows will specify the length of the array.
    ! Whether the array will be allocated or not, if setzero is equal to .true.
    ! the array will be set to zero
    ! --------------------------------------------------------------------------
    subroutine array3_cond_alloc(X,Nrows,Ncols,N3,setzero)

        real(8),    intent(inout), allocatable :: X(:,:,:)  ! 1D array
        integer,    intent(in)                 :: Nrows, Ncols, N3 ! number of rows/columns/length in third dimension
        logical,    intent(in)   , optional    :: setzero ! real number to assign as IC to the array X if this is allocated

        if (allocated(X) .eqv. .false.) then
            allocate( X(Nrows, Ncols, N3) )
        end if

        if ( present(setzero) .and. (setzero .eqv. .true.) ) then
            X=0.0_8
        end if

    end subroutine array3_cond_alloc


end module lib_array



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module lib_array_test

 use lib_array

 implicit none

 contains

! array_cond_alloc_test
!-------------------------------------------------------------------------------
    ! tests the array*_cond_alloc subroutines
    subroutine array_cond_alloc_test

        real(8), allocatable :: X(:)
        real(8), allocatable :: Y(:,:)
        real(8), allocatable :: Z(:,:,:)

        print *, '-------------------------------------- Test array1_cond_alloc'
        print *, 'Allocate X without initialising to zero'
        call array1_cond_alloc(X,3)
        print *, 'done! X= ', X
        print *, 'Try to reallocate with  initialise to zero:'
        call array1_cond_alloc(X,3,.true.)
        print *, 'X should be zero now: X=', X

        print *, '-------------------------------------- Test array2_cond_alloc'
        print *, 'Allocate Y without initialising to zero'
        call array2_cond_alloc(Y,3,2)
        print *, 'done! Y= ', Y
        print *, 'Try to reallocate without  initialise to zero:'
        call array2_cond_alloc(Y,3,2)
        print *, 'Y should be the same: Y=', Y

        print *, '-------------------------------------- Test array3_cond_alloc'
        print *, 'Allocate Z without initialising to zero'
        call array3_cond_alloc(Z,1,3,2)
        print *, 'done! Z= ', Z
        print *, 'Try to reallocate with  initialise to zero:'
        call array3_cond_alloc(Z,1,3,2,.true.)
        print *, 'Z should be zero now: Z=', Z

    end subroutine array_cond_alloc_test

end module lib_array_test
