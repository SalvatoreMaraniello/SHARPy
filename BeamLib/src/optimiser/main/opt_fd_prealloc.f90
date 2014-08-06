!-> Module.- opt_fd
!
!-> Author: Salvatore Maraniello 17 july 2014
!
!-> Language: FORTRAN90, Free format.
!
!-> Important:
!   - This module is identical to opt_fd, but does not define the variables at
!     the interface with python as allocatable.
!   - The module fwd_main_prealloc is used.
!
!-> Description: Module to handle the FD gradient evaluation
!
!-> Subroutines.-
!   - fd_main:
!
! -> Remark: this module is placed at main level to allow the call of:
!       a. input module: this is necessary to perturbe each design variable
!          individually;
!       b. fwd_main: to allow the fwd code execution;

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module opt_fd_prealloc

 use input       ! to perturb the shared design input
 use opt_input   ! to access the FLAGS for design variables and their values XSH
 use fwd_main_prealloc    ! to run the fwd problem
 use opt_perturb
 use opt_cost

 implicit none

 ! variables visible from opt_input
 ! real    ::  XSH(size(FLAG_DESIGN_SHARED), 0:NOPT_MAX )
 ! logical ::  FLAG_DESIGN_SHARED(8+6+6*6+6*6+2)

 contains

 ! subroutine fd_main
 !-------------------
    ! The subroutine starts from the current design cost & constraints and,
    ! perturbing each design variable at the time, computes the FD gradient.
    ! Because the forward problem needs to be executed as many time as the number
    ! of design parameters for each cost/constraints, all the variables required
    ! for the forward problem need to be available.
    subroutine fd_main_prealloc( NumElems,OutFile,Options,           &
                      & Elem,                                        &
                      & NumNodes,                                    &
                      & BoundConds,PosIni,ForceStatic,PhiNodes,      &
                      & OutGrids,                                    &
                      & PsiIni,                                      &
                      & Node, NumDof,                                &
                      & PosDef, PsiDef, InternalForces,              &
                      & NumSteps, t0, dt,                            &
                      & ForceDynAmp,ForceTime,                       &
                      & ForcedVel,ForcedVelDot,                      &
                      & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, &
                      & RefVel, RefVelDot, Quat,                      &
                      & NOPT, fdmode, COST, CONSTR,  &  ! cost and constraint at current design point!
                      & FLAG_COST, W_COST, W_CONSTR, &  ! flags and weights
                      & CONN_CONSTR, CONN_XSH,       &  ! connectivity matrices
                      & DCDXSH, DCONDXSH,            &  ! gradients   ! storage for gradients
                      & Nnode                        )  ! optional argument for cost_node_disp

   ! The following variables are allocatable in opt_fd
     real(8) :: InternalForces(:,:)  ! Internal force/moments at nodes.
     real(8) :: PosIni   (:,:)    ! Initial nodal Coordinates.
     real(8) :: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
     real(8) :: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
     real(8) :: PsiDef (:,:,:)    ! Current element orientation vectors.

     real(8):: t0,dt                               ! Initial time and time step.
     integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
     integer:: NumSteps                            ! Number of time steps.
     integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
     type(xbopts)             :: Options            ! Solution options (structure defined in xbeam_shared).
     type(xbelem), allocatable:: Elem(:)            ! Element information.
     type(xbnode), allocatable:: Node(:)            ! Nodal information.
     integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node

     real(8),      allocatable:: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
     real(8),      allocatable:: ForceTime   (:)   ! Time history of the dynamic nodal forces.
     real(8),      allocatable:: ForcedVel   (:,:) ! Forced velocities at the support.
     real(8),      allocatable:: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.
     real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.
     character(len=25)        :: OutFile           ! Output file.

     real(8),      allocatable:: PosDotDef (:,:)   ! Current nodal position vector.
     real(8),      allocatable:: PsiDotDef (:,:,:) ! Current element orientation vectors.

     real(8),      allocatable:: PosPsiTime(:,:)   ! Position vector/rotation history at beam tip.
     real(8),      allocatable:: ForcesTime(:,:)   ! History of the force/moment vector at the beam root element.
     real(8),      allocatable:: VelocTime(:,:)    ! History of velocities.

     real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
     logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.

     ! Rigid-body variables
     real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
     real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
     real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.
     real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step

        character(len=3), intent(in) :: fdmode          ! finite differences method
        integer                      :: NOPT            ! number of iteration for the optimisation - required to understand which column of XSH to modify
        real(8), intent(in)          :: COST(0:)         ! cost function value at current design point
        real(8), intent(in)          :: CONSTR(:,0:)     ! constraint vector at current design point
        real(8), intent(in)          :: W_CONSTR(:), W_COST(:)       ! arrays with weights/scaling factors...
        integer, intent(in)          :: CONN_CONSTR(:), CONN_XSH(:)  ! connectivity array for constraints
        logical, intent(in)          :: FLAG_COST(:)    ! determines which functions to evaluate
        integer, optional            :: Nnode           ! Number of the node at which to evaluate the displacements.

        real(8), intent(inout) :: DCDXSH  (:,0:)  ! gradient of cost in respect to shared design
        real(8), intent(inout) :: DCONDXSH(:,:,0:)! gradient of constraint in respect to design

        real(8)              :: DXSH( size(XSH(:,NOPT)) ), XSH_COPY( size(XSH(:,NOPT)) ) ! deltas for shared design variables and copy of current design
        integer              :: nn, ii              ! counters for target design variables
        real(8)              :: cost_val                        ! value of cost function at perturbed points
        real(8)              :: constr_val(size(CONN_CONSTR))   ! value of constraint array at perturbed points



        ! ------------------ Compute Deltas for Design current design (NOPT) ---
        where (FLAG_DESIGN_SHARED .eqv. .true.)
            DXSH = fperturb_delta_1d_array( XSH(:,NOPT) )
        elsewhere
            DXSH = 0.0_8
        end where
        ! -------------------------- Testing: result without where statement ---
        !call opt_print_XSH(NOPT)
        !call opt_print_FLAG_DESIGN_SHARED
        !print *, 'Delta only on target design variables: ', DXSH
        !DXSH = fperturb_delta_1d_array( XSH(:,NOPT) )
        !print *, 'Delta without over all array: ', DXSH
        !stop


        ! -------------------------------- Perturb each variable at the time ---
        !call print_shared_input
        XSH_COPY=XSH(:,NOPT)

        do nn=1,size(CONN_XSH)
            ii = CONN_XSH(nn)

            ! compute xsh perturbed:
            XSH(  :, NOPT ) = XSH_COPY
            XSH( ii, NOPT ) = XSH_COPY( ii ) + DXSH( ii )
            print  '(A8,I4,A8,I4,A3,$)', 'NOPT:', NOPT,  'Var. No.', ii, '='
            print '(E14.8,A8,E14.8,A8)', XSH_COPY(ii), '(init.)', XSH( ii, NOPT ), '(pert.)'
            call opt_unpack_DESIGN_SHARED(NOPT)

            ! execute forward
            call fwd_problem_prealloc(NumElems,OutFile,Options,    &    ! from input_setup
                 &                        Elem,           &   ! from opt_main_xxx
                 &                    NumNodes,           &   ! from input_elem
                 &  BoundConds,PosIni,ForceStatic,PhiNodes,    &   ! from input_node
                 &                    OutGrids,    &   ! from pt_main_xxx
                 &                      PsiIni,    &   ! from xbeam_undef_geom
                 &                Node, NumDof,    &   ! from xbeam_undef_dofs
                 & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 &            NumSteps, t0, dt,    &   ! input_dyn_setup
                 &       ForceDynAmp,ForceTime,    &   ! input_dynforce
                 &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! to be allocated in fwd_dynamic_presolve and out of dynamic analysis
                 & RefVel, RefVelDot, Quat)

            ! -------------------------------- allocate cost and constraints ---
            cost_val= cost_global( FLAG_COST, W_COST, PosIni, PosDef, Elem(:)%Mass(1,1),Elem(:)%Length )
            constr_val = cost_constraints( W_CONSTR, CONN_CONSTR, PosIni, PosDef, Elem(:)%Mass(1,1),Elem(:)%Length )

            ! -------------------------------------------- Compute Gradients ---
            DCDXSH(nn,NOPT)     = (   cost_val - COST(NOPT)     ) /  DXSH(ii)
            DCONDXSH(:,nn,NOPT) = ( constr_val - CONSTR(:,NOPT) ) /  DXSH(ii)


            !print '(A10,F16.8,A10,F16.8,A10,F16.8,A10)', 'cost: ', COST(NOPT), '(old)', cost_val, '(new)' ,  (cost_val - COST(NOPT)), '(delta)'
            print '(A10,E14.8,A10,E14.8,A10,E14.8,A10)', 'cost: ', COST(NOPT), '(old)', cost_val, '(new)' ,  (cost_val - COST(NOPT)), '(delta)'



            print '(A,E12.6)', 'Cost Gradient:', DCDXSH(nn,NOPT)
            print '(A,$)', 'Constrain Gradient:'
            print '(E12.6)', DCONDXSH(:,nn,NOPT)

        end do

        ! go back to original design
        XSH(:,NOPT)=XSH_COPY

    end subroutine fd_main_prealloc


 end module opt_fd_prealloc


