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
!
!   - Cost Functions:
!     - cost_node_disp: returns the absolute value of a nodal displacement
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_cost

implicit none

contains


!-------------------------------------------------------------------------------
    function cost_node_disp(PosIni,PosDef,Nnode)
    ! the function returns the absolute value of the displacements of the node
    ! Nnode.
    ! If Nnode is not passed, this wil be taken as the higher node in the
    ! model (presumely the tip?)

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

end module opt_cost



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module opt_cost_test

 use opt_cost
 implicit none

 contains


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

        print *, PosIni
        print *, PosDef

        print *, 'Node 1:', cost_node_disp(PosIni,PosDef,1)
        print *, 'Node 2:', cost_node_disp(PosIni,PosDef,2)
        print *, 'Node 3:', cost_node_disp(PosIni,PosDef,3)
        print *, 'Node 4:', cost_node_disp(PosIni,PosDef)


    end subroutine cost_node_disp_test

end module opt_cost_test







