!-> Module.- opt_cost Salvatore Maraniello 18/07/2014
!
!-> Author: Salvatore Maraniello
!
!-> Language: FORTRAN90, Free format.
!
!-> Description: the module contains a collection of:
!   a. Cost functions for optimisation problem;
!
!-> Subroutines.-
!   - Cost Functions:
!     - cost_node_disp: returns the absolute value of a nodal displacement
!
!
!
!-> Ramarks.-
! - The module required xbeam_shared to access xbelem type
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_cost

 use opt_cost_utl
 use xbeam_shared

 implicit none

 contains


!-------------------------------------------------------------------------------
! function cost_node_disp
!------------------------
    ! the function returns the absolute value of the displacements of the node
    ! Nnode.
    ! If Nnode is not passed, this wil be taken as the higher node in the
    ! model (presumely the tip?)
    function cost_node_disp(PosIni,PosDef,Nnode)

        real(8)                :: cost_node_disp  ! absolute value of node displacement
        real(8),    intent(in) :: PosIni (:,:)    ! Initial nodal Coordinates.
        real(8),    intent(in) :: PosDef (:,:)    ! Current nodal position vector. (sm: local coordinates)
        integer, optional      :: Nnode           ! Number of the node at which to evaluate the displacements.

        real(8)                :: dispvec(3)    ! displacement vector
        integer                :: sh(2)         ! shape of Pos* arrays

        if ( present(Nnode)) then
            dispvec =  PosDef(Nnode,:) - PosIni(Nnode,:)
        else
            !Nnode = size( PosIni(:,1) ) !<-- generates error
            sh = shape(PosIni)
            dispvec =  PosDef( sh(1) ,:) - PosIni( sh(1) ,:)
        end if
        cost_node_disp=sqrt( dot_product(dispvec,dispvec) )

    end function cost_node_disp



!-------------------------------------------------------------------------------
! function cost_total_mass
!-------------------------
    ! the function returns the total mass of the structure.
    ! Geometry data are taken from the undeformed structure configuration
    function cost_total_mass(MassVector,LengthVector)
        real(8),      intent(in) :: MassVector(:)   ! Mass of each element
        real(8),      intent(in) :: LengthVector(:) ! Length of each element
        real(8)                  :: cost_total_mass ! structure mass
        integer                  :: ii

        cost_total_mass=0.0_8
        do ii=1,size(MassVector)
            !cost_total_mass = cost_total_mass + Elem(ii)%Mass(1,1)*Elem(ii)%Length
            cost_total_mass = cost_total_mass + MassVector(ii)*LengthVector(ii)
        end do

    end function cost_total_mass




!-------------------------------------------------------------------------------
! function cost_global:
! ---------------------
    ! This function evaluates the global cost for those cases in which the cost
    ! function is a weighted combination of several cost functions.
    ! the function takes in input all the input of all the cost functions
    ! implemented in the model.

    function cost_global( FLAG_COST, W_COST,        &  ! to determine which functions and weight to evaluate
                        & PosIni,PosDef,            &  ! input for cost_node_disp
                        & MassVector,LengthVector,  &  ! input for cost_total_mass
                        & Nnode                     &  ! optional input for cost_node_disp
                        )                              ! add here other input

        ! see cost_utl_allocate_flags_and_weights
        logical         , intent(in) :: FLAG_COST  (:)    ! determines which functions to evaluate
        real(8)         , intent(in) :: W_COST     (:)    ! arrays with weights/scaling factors...
        ! from cost_node_disp
        real(8),    intent(in) :: PosIni (:,:)    ! Initial nodal Coordinates.
        real(8),    intent(in) :: PosDef (:,:)    ! Current nodal position vector
        integer, optional      :: Nnode           ! Number of the node at which to evaluate the displacements.
        ! from cost_total_mass
        real(8),      intent(in) :: MassVector(:)   ! Mass of each element
        real(8),      intent(in) :: LengthVector(:) ! Length of each element
        ! output
        real(8) :: cost_global
        integer                :: nn              ! counter for FLAG_COST & W_COST

        ! Initialise:
        cost_global=0.0_8

        ! ------------------------------------- Evaluate functions one by one ---

        do nn=1,size(FLAG_COST)
            if (FLAG_COST(nn) .eqv. .true.) then
                cost_global = cost_global +                                    &
                            & W_COST(nn)*eval_function( nn,                     &
                                                      & PosIni,PosDef,          &
                                                      & MassVector,LengthVector,&
                                                      & Nnode                   &
                                                                               )
            end if
       end do

    end function cost_global


!-------------------------------------------------------------------------------
! function constraint_vector:
! --------------------------
    ! This function evaluates the constraint vector for a certain status
    ! the function takes in input all the input of all the cost functions
    ! implemented in the model.
    ! As per cost_global, FLAG and WEIGHT arrays are required.
    ! The connectivity array (cost_utl_build_constraints_connectivity) is also
    ! needed.
    function cost_constraints( W_CONSTR, CONN_CONSTR,   &  ! to determine which functions and weight to use
                            & PosIni,PosDef,            &  ! input for cost_node_disp
                            & MassVector,LengthVector,  &  ! input for cost_total_mass
                            & Nnode                     &  ! optional input for cost_node_disp
                            )                              ! add here other input

        ! see cost_utl_allocate_flags_and_weights
        !logical         , intent(in) :: FLAG_CONSTR  (:)    ! determines which functions to evaluate
        real(8)         , intent(in) :: W_CONSTR     (:)    ! arrays with weights/scaling factors...
        integer         , intent(in) :: CONN_CONSTR  (:)    ! connectivity array for constraints
        ! from cost_node_disp
        real(8),    intent(in) :: PosIni (:,:)    ! Initial nodal Coordinates.
        real(8),    intent(in) :: PosDef (:,:)    ! Current nodal position vector
        integer, optional      :: Nnode           ! Number of the node at which to evaluate the displacements.
        ! from cost_total_mass
        real(8),      intent(in) :: MassVector(:)   ! Mass of each element
        real(8),      intent(in) :: LengthVector(:) ! Length of each element
        ! output
        real(8) :: cost_constraints(size(CONN_CONSTR))

        integer                :: ii        ! counter for cost_constraints
        integer                :: nn        ! counter for FLAG_CONSTR & W_CONSTR


        ! Initialise:
        cost_constraints=0.0_8

        do ii=1,size(cost_constraints)
            nn=CONN_CONSTR(ii)
            cost_constraints(ii) = W_CONSTR(nn)*eval_function( nn,             &
                                                            & PosIni,PosDef,   &
                                                            & MassVector,LengthVector,&
                                                            & Nnode            &
                                                            &                  )
       end do

    end function cost_constraints



!-------------------------------------------------------------------------------
! function eval_function:
! --------------------------
    ! This function evaluates the cost/constraint functions. The function id nn
    ! (see cost_utl_allocate_flags_and_weights) is required for the evaluation.
    function eval_function  ( NN,                   &  ! to determine which function
                            & PosIni,PosDef,        &  ! input for cost_node_disp
                            & MassVector,LengthVector,&! input for cost_total_mass
                            & Nnode                 &  ! optional input for cost_node_disp
                            )                          ! add here other input

        integer                :: NN              ! function ID according to cost_utl_allocate_flags_and_weights

        ! from cost_node_disp
        real(8),    intent(in) :: PosIni (:,:)    ! Initial nodal Coordinates.
        real(8),    intent(in) :: PosDef (:,:)    ! Current nodal position vector
        integer, optional      :: Nnode           ! Number of the node at which to evaluate the displacements.
        ! from cost_total_mass
        real(8),      intent(in) :: MassVector(:)   ! Mass of each element
        real(8),      intent(in) :: LengthVector(:) ! Length of each element
        ! output
        real(8) :: eval_function


        select case (NN)

            case (1)
            ! cost_node_disp
                if ( present(Nnode)) then
                    eval_function = cost_node_disp( PosIni,PosDef,Nnode )
                else
                    eval_function = cost_node_disp( PosIni,PosDef )
                end if

            case (2)
            !cost_mass_total
                eval_function = cost_total_mass(MassVector,LengthVector)

            case default
                stop 'ID not valid! The ID must be according to cost_utl_allocate_flags_and_weights!'

        end select

    end function


end module opt_cost


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module opt_cost_test

 use opt_cost
 implicit none

 character(len=*), parameter :: nfmt='(A,F12.6)'
 character(len=*), parameter :: mfmt='(F12.6,F12.6,F12.6)'


 contains

! Test for cost_node_disp
!-------------------------------------------------------------------------------
    subroutine cost_node_disp_test

        real(8) :: PosIni (4,3)    ! Initial nodal Coordinates.
        real(8) :: PosDef (4,3)    ! Current nodal position vector. (sm: local coordinates)

        PosDef = transpose( reshape( (/    5.5 ,  100. ,  1.7 ,  &
                                      &    4.0 ,  2.90 ,  5.0 ,  &
                                      &    4.4 ,  1.00 ,  -4.2 ,  &
                                      &    4.3 ,  0.00 ,  2.1    /) , (/3,4/)) )

        PosIni = transpose( reshape( (/    3.5 ,   0.9 ,  1.0 ,  &
                                      &    4.0 ,  -0.1 ,  1.0 ,  &
                                      &    4.4 ,  -1.0 ,  -2.2 ,  &
                                      &    4.2 ,   0.1 ,  2.0    /) , (/3,4/)) )

        print *, '---------------------------------------------- cost_node_disp'

        print *, 'Initial:'
        print mfmt, PosIni
        print *, 'Final:'
        print mfmt, PosDef

        print nfmt, 'Node 1:', cost_node_disp(PosIni,PosDef,1)
        print nfmt, 'Node 2:', cost_node_disp(PosIni,PosDef,2)
        print nfmt, 'Node 3:', cost_node_disp(PosIni,PosDef,3)
        print nfmt, 'Node 4:', cost_node_disp(PosIni,PosDef)

        ! other way to define the format is to place it here (after the print!)
        ! print 100, 'Node 1:', cost_node_disp(PosIni,PosDef,1)
        ! 100 FORMAT(A,F10.6)
    end subroutine cost_node_disp_test



! Test for cost_mass-total
!-------------------------------------------------------------------------------
    subroutine cost_mass_total_test

        type(xbelem) :: Elem      (5) ! Element information

        real(8) :: M=100.0_8, L=2.0_8, Mtotal
        integer :: ii
        real(8) :: ratio

        print *, '--------------------------------------------- cost_mass_total'

        print '(A12,A12,A12)', 'Elem.', 'Mass.', 'Length'
        do ii=1,size(Elem)
            ratio = real(ii,8)/real(size(Elem),8)
            Elem(ii)%Mass(1,1)=ratio*M
            Elem(ii)%Length   =(ratio**2)*L
            print '(I12,F12.6,F12.6)', ii, Elem(ii)%Mass(1,1), Elem(ii)%Length
        end do

        Mtotal = cost_total_mass(Elem(:)%Mass(1,1),Elem(:)%Length)
        print '(A,F16.8)', 'Total Mass:', Mtotal
        print '(A,F16.8)', 'Expected:', ( 1.0_8 + (1.0_8 + 2.0_8**3 + 3.0_8**3 + 4.0_8**3)/(5.0_8**3) )*M*L

    end subroutine cost_mass_total_test

end module opt_cost_test







