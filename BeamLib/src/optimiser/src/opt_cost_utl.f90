!-> Module.- opt_cost_utl Salvatore Maraniello 18/07/2014
!
!-> Author: Salvatore Maraniello
!
!-> Language: FORTRAN90, Free format.
!
!-> Description: the module contains a collection of methods to organise the
!   input related to cost and constraint function in input from the opt_input module
!
!-> Subroutines.-
!
!   - cost_utl_allocate_flags_and_weights: this is called during the execution
!     of opt_input_cost and allocates the FLAGS nad WEIGHT vectors required to
!     define cost and contraints
!
!-> Functions.-
!
!   - cost_utl_get_fun_position: associates each function to an element of the
!     FLAGS and WEIGHT arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_cost_utl

implicit none

contains


 !------------------------------------------------------------------------------
    ! Reads input from opt_input and allocates FLAG vectors and weights.
    ! the FLAG vectors are used to identify whether a function is associated to
    ! the cost, constraint or is simply not used.
    ! Each element of the FLAG and WEIGHT vector is associated to a function via
    ! the cost_utl_get_fun_position function.
    ! The WEIGHTS vector can be used to scale or to combine different function
    subroutine cost_utl_allocate_flags_and_weights(FUNID,ADDTO,WEIGHT, &   !  input
                                  & FLAG_COST,FLAG_CONSTR,W_COST,W_CONSTR) !  output

        ! --- input ---
        real(8),           intent(in)  ::  WEIGHT   ! wieght to apply to function
        character(len= 4), intent(in)  ::   ADDTO   ! determines whether to add the function
                                                    ! to cost, constraint or none
                                                    ! values: 'cost', 'cstr', 'none'
        character(LEN=10), intent(in)  ::   FUNID   ! link to function
                                                    ! 'node_disp': nodal displacement

        ! --- output ---
        logical         , intent(out) :: FLAG_COST  (1)             ! Flag array for cost functions
        logical         , intent(out) :: FLAG_CONSTR(1)             ! Flag array for cost functions
        real(8)         , intent(out) :: W_COST  (size(FLAG_COST))    ! arrays with weights/scaling factors...
        real(8)         , intent(out) :: W_CONSTR(size(FLAG_COST))    ! ...for cost and constraint functions

        integer            ::           nn          ! counter

        nn = cost_utl_get_fun_position(FUNID)

        select case (ADDTO)
            case ('none')
                FLAG_COST(nn)=.false.;    FLAG_CONSTR(nn)=.false.
                   W_COST(nn)=0.0_8;         W_CONSTR(nn)=0.0_8;
            case ('cost')
                FLAG_COST(nn)=.true.;     FLAG_CONSTR(nn)=.false.
                   W_COST(nn)=WEIGHT;        W_CONSTR(nn)=0.0_8;
            case ('cstr')
                FLAG_COST(nn)=.false.;    FLAG_CONSTR(nn)=.true.
                   W_COST(nn)=0.0_8;         W_CONSTR(nn)=WEIGHT;
        end select

    end subroutine cost_utl_allocate_flags_and_weights


 ! ----------------------------------------------------------------------------
    function cost_utl_get_fun_position(FUNID)
    ! Provides a mapping between cost/constraint function and a position in the
    ! arrays FLAg_COST, FLAG_CONSTR, W_COST, W_CONSTR
        integer                        ::   cost_utl_get_fun_position
        character(LEN=10), intent(in)  ::   FUNID   ! link to function

        select case(FUNID)

            case ('node_disp')
                cost_utl_get_fun_position=1

            case default
                stop 'FUNID not found! Check IDs for cost and constraint functions in the opt_input module'

        end select

    end function cost_utl_get_fun_position


end module opt_cost_utl









