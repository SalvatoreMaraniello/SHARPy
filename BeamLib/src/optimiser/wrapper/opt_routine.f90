!-> Program.- MAIN - 06Jan2011 Updated: july 2014
!
!-> Author.- Henrik Hesse  (h.hesse09@imperial.ac.uk)
!            Rafa Palacios (rpalacio@imperial.ac.uk)
!            Rob Simpson   (rjs10@imperial.ac.uk) copied from main_andrea.f90
!            Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!   Same as opt_main but in a routine form.
!   The routine is suitable to be called by a main program or a wrapper.
!   The routine can run the direct solution, sensitivity analysis and
!   unconstrained optimisation.
!   Though an optimisation option is present, the main driver for any optimisation
!   should be located 'outside'.
!
! -> Developer Notes:
! a. in function call, in & out are intended as:
!     - in : goes inside the function
!     - out: is returned from the function
! b. Memory usage:
!     - variables for the forward problem are allocated in the fwd_main module.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module opt_routine

 use xbeam_shared   ! Forward Solution Modules
 use xbeam_undef
 use cbeam3_solv
 use xbeam_solv
 use xbeam_perturb
 use input
 use lib_out
 use opt_input      ! Optimisation Modules
 use opt_fd
 use fwd_main
 use lib_perf
 use opt_shared
 use opt_perturb
 use opt_cost
 use opt_cost_utl
 use opt_driver

 implicit none


contains


subroutine opt_main(NumElems, NumNodes, pCOST, W_COST)






 real(8):: t0,dt                               ! Initial time and time step.
 integer:: i,j                                 ! Counter.
 integer, intent(inout) :: NumElems               ! Number of elements
 integer, intent(inout) :: NumNodes               ! Number of nodes in the model.
 integer:: NumSteps                            ! Number of time steps.
 integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
 type(xbopts)             :: Options            ! Solution options (structure defined in xbeam_shared).
 type(xbelem), allocatable:: Elem(:)            ! Element information.
 type(xbnode), allocatable:: Node(:)            ! Nodal information.
 integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node
 real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
 real(8),      allocatable:: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
 real(8),      allocatable:: ForceTime   (:)   ! Time history of the dynamic nodal forces.
 real(8),      allocatable:: ForcedVel   (:,:) ! Forced velocities at the support.
 real(8),      allocatable:: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.
 real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.
 real(8),      allocatable:: InternalForces(:,:)  ! Internal force/moments at nodes.
 logical,      allocatable:: OutGrids(:)       ! Grid nodes where output is written.
 character(len=25)        :: OutFile           ! Output file.

 real(8),      allocatable:: PosIni   (:,:)    ! Initial nodal Coordinates.
 real(8),      allocatable:: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
 real(8),      allocatable:: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
 real(8),      allocatable:: PsiDef (:,:,:)    ! Current element orientation vectors.
 real(8),      allocatable:: PosDotDef (:,:)   ! Current nodal position vector.
 real(8),      allocatable:: PsiDotDef (:,:,:) ! Current element orientation vectors.

 real(8),      allocatable:: PosPsiTime(:,:)   ! Position vector/rotation history at beam tip.
 real(8),      allocatable:: ForcesTime(:,:)   ! History of the force/moment vector at the beam root element.
 real(8),      allocatable:: VelocTime(:,:)    ! History of velocities.
 real(8),      allocatable:: Time(:)           ! Discrete time vector in the dynamic simulation.

 ! Rigid-body variables
 real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
 real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
 real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.
 real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step

 ! Optimisation
 character(len=3) :: gradmode ! gradient method
 character(len=3) :: solmode  ! solution mode
 character(len=3) :: fdmode   ! finite differences method
 integer          :: NOPTMAX  ! Max Number of iterations for the optimisation
 integer          :: NOPT     ! number of iteration for the optimisation

 logical          :: FLAG_COST  (NCOSTFUNS) ! Flag array for cost funcitons
 logical          :: FLAG_CONSTR(NCOSTFUNS) ! Flag array for cost funcitons
real(8)          :: W_COST  (NCOSTFUNS)   ! arrays with weights/scaling factors... ! no need to define size here
 !!!!!!! real(8)          :: W_COST  (:)                                           ! this works just as fine
 real(8)          :: W_CONSTR(NCOSTFUNS)   ! ...for cost and constraint functions
 integer, allocatable :: CONN_CONSTR(:), CONN_XSH(:)   ! connectivity matrix for contrains and design variables array
 real(8), target, allocatable :: COST(:)                 ! cost function
 real(8), allocatable :: CONSTR(:,:)             ! constraint vector

 real(8), allocatable :: DCDXSH  (:,:)  ! gradient of cost in respect to shared design
 real(8), allocatable :: DCONDXSH(:,:,:) ! gradient of constraints in respect to design
 ! gradients ordering:
 !     DCDXSH(   ii,NOPT)           DCONDXSH(nn,ii,NOPT)
 ! where:
 !   nn-th constraint // ii-th design variable // NOPT: optimisation iteration

 ! -----------------------------------------------------------------------------
 ! Python Interface
 ! -----------------------------------------------------------------------------
 ! Create Pointers to arrays
 real(8), pointer, dimension(:), intent(out) :: pCOST


 !------------------------------------------------------------------------------


 call tic
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! set shared design variables for fwd problem
 call input_setup (NumElems,OutFile,Options)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Optimiser Input
 call opt_setup(gradmode,solmode,fdmode,NOPTMAX)
 if ( (solmode == 'OPT') .or. (solmode == 'SNS') ) then
     call opt_input_cost(FLAG_COST,FLAG_CONSTR,W_COST,W_CONSTR)
     ! note: cost_utl_build_connectivity is in the opt_cost_utl module.
     call cost_utl_build_connectivity(FLAG_CONSTR,CONN_CONSTR)
     call cost_utl_build_connectivity(FLAG_DESIGN_SHARED,CONN_XSH)
     !!! for testing...
     !print *, 'FLAG_COST:', FLAG_COST
     !print *, 'FLAG_CONSTR:', FLAG_CONSTR
     !print *, 'Constr. Connectivity: ', CONN_CONSTR
     !print *, 'Constr. XSH: ', CONN_XSH

     if (solmode == 'SNS') then        ! this is only needed for the allocation
        NOPTMAX=0
     end if
     allocate(CONSTR( size(CONN_CONSTR),0:NOPTMAX ) ); CONSTR=0.0_8
     allocate(  COST(                   0:NOPTMAX ) ); COST=0.0_8
     allocate(   DCDXSH(                    size(CONN_XSH),0:NOPTMAX ) ); DCDXSH=0.0_8
     allocate( DCONDXSH( size(CONN_CONSTR), size(CONN_XSH),0:NOPTMAX ) ); DCONDXSH=0.0_8
 end if


!NOPT=0 ! NOPT=0 is assumed inside opt_setup to allocate XSH
!do while (NOPT<=NOPTMAX)
do NOPT=0,NOPTMAX

    call fwd_problem( NumElems,OutFile,Options,    &   ! from input_setup
                 &                        Elem,    &   ! from opt_main_xxx
                 &                    NumNodes,    &   ! from input_elem
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

     ! Store results in text file.
     open (unit=11,file=OutFile(1:17)//'_def.txt',status='replace')
     call output_elems (11,Elem,PosDef,PsiDef)
     close (11)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Loop control
     select case (solmode)

        case default ! ('FWD')
            print *, 'Forward Problem Completed!'
            exit

        case ('OPT','SNS')
            ! -------------- Evaluate cost and Constraints at current design ---
            COST(NOPT) = cost_global( FLAG_COST, W_COST,  &
                                    & PosIni, PosDef,     &
                                    & Elem                   )

            CONSTR(:,NOPT) = cost_constraints( W_CONSTR, CONN_CONSTR, &
                                            & PosIni,PosDef,         &
                                            & Elem                   )
            !PRINT *, 'COST',    COST
            !PRINT *, 'CONSTR:', CONSTR

            ! ----------------------------------------- Sensitivity analysis ---
            print *, 'Sensitivity Analysis started. NOPT=', NOPT
            select case (gradmode)
                case ('FDF')
                    print *, 'Gradients will be computed via Finite Differences'
                    call fd_main( NumElems,OutFile,Options,                    &
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
                                & RefVel, RefVelDot, Quat,                     &
                                & NOPT, fdmode, COST, CONSTR, &  ! cost and constraint at current design point!
                                & FLAG_COST, W_COST, W_CONSTR,&  ! flags and weights
                                & CONN_CONSTR, CONN_XSH,      &  ! connectivity matrices
                                & DCDXSH, DCONDXSH            )  ! gradients
                case default
                    stop 'Only FD method available!'

            end select


            ! ---------------------------------------------------- Save Etc. ---
            !
            ! -> To be Developed



            ! ------------------------------- Terminate sensitivity Analysis ---
            if (solmode == 'SNS') then
                print *, 'Sensitivity Analysis Completed!'
                print *, 'COST',    COST
                print *, 'CONSTR:', CONSTR
                print *, 'cost gradient', DCDXSH
                print *, 'constraints gradient', DCONDXSH
                exit
            end if

     end select


    ! ------------------------------------------------ Set next Design Point ---
    if (NOPT < NOPTMAX) then

        call simple_driver(XSH,COST,CONSTR,DCDXSH,DCONDXSH,NOPT,CONN_XSH,CONN_CONSTR)
        PRINT *, 'UPDATE DESIGN for optimisation loop No.', NOPT


        ! this line shows that the update works correctly and the input change
        ! at each optimisation step
        !XSH(:,NOPT+1)=XSH(:,NOPT)+0.01_8 * XSH(:,NOPT)
        call opt_unpack_DESIGN_SHARED(NOPT+1)

    end if

end do


 print *, 'Optimisation Terminated: '
 print '(A15,$)', 'cost: '
 print '(F10.6,$)', COST
 print *, ' '

 print '(A15,$)', 'constraints: '
 print '(F10.6,$)', CONSTR
 print *, ' '

 print '(A15,$)', 'cost grad.: '
 print '(F10.6,$)', DCDXSH
 print *, ' '

 print '(A15,$)', 'constr. grad.: '
 print '(F10.6,$)', DCONDXSH
 print *, ' '

 call toc

 ! assign pointers to allow python interface
 pCOST => COST
 print *, pCOST

 end subroutine opt_main

end module opt_routine