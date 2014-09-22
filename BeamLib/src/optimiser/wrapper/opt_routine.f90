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
! c. - Design Variables from the original code are updated in the input module
!    using the update_shared_input routine.
!    - Design Variables/Input added for the optimisation are passed as dummy
!    arguments to the forward solver/presolver
!       - BeamSpanStiffness: element by element stiffness matrix
!       - BeamSpanMass: element by element mass matrix
!       - PhiNodes: node by node twist array. This has the same meaning as
!         PhiNodes.
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
 use opt_fd_prealloc
 use fwd_main_prealloc
 use lib_perf
 use opt_shared
 use opt_perturb
 use opt_cost
 use opt_cost_utl
 use opt_driver
 !use test, only: pack_xbopts


 implicit none


contains



!subroutine arg_test(NumElem, pCOST, W_COST)
!   integer, intent(inout) :: NumElems          ! Number of elements





subroutine opt_main( NumElems, NumNodes,                                       & ! Problem SetUp
                   & IN_NumNodesElem , IN_ElemType, IN_TestCase, IN_BConds,    & ! Problem Setup Shared
                   & IN_BeamLength1, IN_BeamLength2,                           & ! Design Variables
                   & BeamSpanStiffness, BeamSpanMass,                          & ! Properties along the span
                   & PhiNodes,                                                 & ! Twist angle along the span
                   & ForceStatic,                                              & ! Static Nodal forces
                   & IN_TipMass, IN_TipMassY, IN_TipMassZ,                     &
                   & IN_Omega,                                                 &
                   & FollowerForce,FollowerForceRig,PrintInfo,                 & ! Options
                   & OutInBframe,OutInaframe,ElemProj,MaxIterations,           &
                   & NumLoadSteps,Solution,DeltaCurved,MinDelta, NewmarkDamp,  &
                   & PosIni, PsiIni,                   & ! Initial Pos/Rot
                   & PosDef, PsiDef, InternalForces,   & ! Output Static
                   & DensityVector, LengthVector,      & ! Design output
                   & NumSteps, Time,                   & ! input_dynsetup
                   & ForceTime, ForceDynAmp,           & ! input_dynforce
                   & ForcedVel, ForcedVelDot,          & ! input_forcedvel
                   & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! output from sol 202, 212, 302, 312, 322
                   & RefVel, RefVelDot, Quat)

! ---------------------------------------------------------------- Removed Input
! these variables have been removed from the input interface as they are managed
! in the python wrapper. They are set to zero and passed to the fwd solver to
! allow consistency with old code.
   real(8)  :: IN_ThetaRoot=0.0_8
   real(8)  :: IN_ThetaTip=0.0_8
   real(8)  :: IN_BeamStiffness(6,6)=0.0_8
   real(8)  :: IN_BeamMass(6,6)=0.0_8
   real(8)  :: IN_SectWidth=0.0_8,IN_SectHeight=0.0_8
   real(8)  :: IN_ExtForce(3)
   real(8)  :: IN_ExtMomnt(3)

! ---------------------------------------------------------------------- Options
        type(xbopts)                      :: Options ! not dummy
        logical     ,intent(in), optional :: FollowerForce
        logical     ,intent(in), optional :: FollowerForceRig
        logical     ,intent(in), optional :: PrintInfo
        logical     ,intent(in), optional :: OutInBframe
        logical     ,intent(in), optional :: OutInaframe
        integer     ,intent(in), optional :: ElemProj
        integer     ,intent(in), optional :: MaxIterations
        integer     ,intent(in), optional :: NumLoadSteps
        !integer     ,intent(in), optional :: NumGauss !NumNodesElem used instead
        integer     ,intent(in), optional :: Solution
        real(8)     ,intent(in), optional :: DeltaCurved
        real(8)     ,intent(in), optional :: MinDelta
        real(8)     ,intent(in), optional :: NewmarkDamp

! ------------------------------------------------------ Design Variables Shared
! These are normally inported from the input module through input_setup.
   real(8), intent(inout)  :: IN_BeamLength1, IN_BeamLength2
   real(8), intent(inout)  :: IN_TipMass
   real(8), intent(inout)  :: IN_TipMassY
   real(8), intent(inout)  :: IN_TipMassZ
   real(8), intent(inout)  :: IN_Omega

! ---------------------------------------------------------------- Problem Setup
! remark: NumNodes is recomputed inside input_elem (input module)
! It is passed as an input to check that the preallocated variables in the python
! wrapper have the same size as those allocated inside the fortran code.
   integer, intent(inout) :: NumElems          ! Number of elements
   integer, intent(inout) :: NumNodes          ! Number of nodes in the model.
   real(8)                :: BeamSpanStiffness(NumElems,6,6) ! Element by Element Stiffness matrix
   real(8)                :: BeamSpanMass(NumElems,6,6)      ! Element by Element Mass matrix
   real(8)                :: PhiNodes(NumNodes)              ! Initial Twist angle (opt). This variable overwrites PhiNodes in the input setup

! --------------------------------------------------------- Problem Setup Shared
! These variables appear as shared in the input module.
    integer,          intent(inout)  :: IN_NumNodesElem   ! Number of nodes on each element.
    character(len=4), intent(inout)  :: IN_ElemType       ! ='STRN','DISP'
    character(len=4), intent(inout)  :: IN_TestCase       ! Define the test case (CANT,ARCH).
    character(len=2), intent(inout)  :: IN_BConds         ! ='CC': Clamped-Clamped


! ----------------------------------------------------------------------- Output
! General
 real(8), intent(inout) :: DensityVector (NumElems)       ! Linear Density of each element of the beam. To be used for cost evaluation.
 real(8), intent(inout) :: LengthVector (NumElems)        ! Length of each element of the beam. To be used for cost evaluation.

! Static
 real(8), intent(inout) :: ForceStatic(NumNodes,6)         ! Internal force/moments at nodes.
 real(8), intent(inout) :: InternalForces(NumNodes,6)      ! Internal force/moments at nodes.
 real(8), intent(inout) :: PosIni (NumNodes,3)             ! Initial nodal Coordinates.
 real(8), intent(inout) :: PsiIni (NumElems,MaxElNod,3)    ! Initial element orientation vectors (CRV)
 real(8), intent(inout) :: PosDef (NumNodes,3)             ! Current nodal position vector. (sm: local coordinates)
 real(8), intent(inout) :: PsiDef (NumElems,MaxElNod,3)    ! Current element orientation vectors.

! Dynamic Input/Output
 integer, intent(in)    :: NumSteps          ! Number of time steps
 real(8), intent(inout) :: Time(NumSteps+1)           ! Discrete time vector in the dynamic simulation.
 real(8), intent(inout) :: ForceDynAmp (NumNodes,3) ! Amplitude of the applied dynamic nodal forces.
 real(8), intent(inout) :: ForceTime   (NumSteps+1)   ! Time history of the dynamic nodal forces.
 real(8), intent(inout) :: ForcedVel   (NumSteps+1,6) ! Forced velocities at the support.
 real(8), intent(inout) :: ForcedVelDot(NumSteps+1,6) ! Derivatives of the forced velocities at the support.

 real(8), intent(inout) :: PosDotDef (NumNodes,3)   ! Current nodal position vector.
 real(8), intent(inout) :: PsiDotDef (NumElems,MaxElNod,3) ! Current element orientation vectors.
 real(8), intent(inout) :: PosPsiTime(NumSteps+1,6)   ! Position vector/rotation history at beam tip.
 real(8), intent(inout) :: VelocTime(NumSteps+1,NumNodes)    ! History of velocities.
 real(8), intent(inout) :: DynOut   ((NumSteps+1)*NumNodes,6)    ! Position of all nodes wrt to global frame a for each time step

 ! Rigid-body variables
 real(8), intent(inout) :: RefVel   (NumSteps+1,6) ! Velocities of reference frame at the support (rigid body).
 real(8), intent(inout) :: RefVelDot(NumSteps+1,6) ! Derivatives of the velocities of reference frame a.
 real(8), intent(inout) :: Quat     (4)            ! Quaternions to describe propagation of reference frame a.


 integer:: i,j                                 ! Counter.
 integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).

 type(xbelem), allocatable:: Elem(:)            ! Element information.
 type(xbnode), allocatable:: Node(:)            ! Nodal information.
 integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node




 logical,      allocatable:: OutGrids(:)       ! Grid nodes where output is written.
 character(len=25)        :: OutFile           ! Output file.


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
 ! real(8), pointer, dimension(:), intent(out) :: pCOST

 call tic

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! pack Options
        if (present( FollowerForce )) then
            Options%FollowerForce = FollowerForce
        end if
        if (present( FollowerForceRig )) then
        Options%FollowerForceRig = FollowerForceRig
        end if
        if (present(  PrintInfo )) then
        Options%PrintInfo = PrintInfo
        end if
        if (present( OutInBframe )) then
        Options%OutInBframe = OutInBframe
        end if
        if (present( OutInaframe )) then
        Options%OutInaframe = OutInaframe
        end if
        if (present( ElemProj )) then
        Options%ElemProj = ElemProj
        end if
        if (present( MaxIterations )) then
        Options%MaxIterations = MaxIterations
        end if
        if (present( NumLoadSteps )) then
        Options%NumLoadSteps = NumLoadSteps
        end if
        !if (present( NumGauss )) then   ! NumNodesElem instead
        !Options%NumGauss = NumGauss
        !end if
        if (present( Solution )) then
        Options%Solution = Solution
        end if
        if (present( DeltaCurved )) then
        Options%DeltaCurved = DeltaCurved
        end if
        if (present( MinDelta )) then
        Options%MinDelta = MinDelta
        end if
    if (present( NewmarkDamp )) then
        Options%NewmarkDamp = NewmarkDamp
    end if

! Solver options.
! this ckeck is inconsistent with the input setup - as both ElemType and
! Options%Solution are input - but has been kept for consistency with the
! original code
   select case (Options%Solution)
       case (102,112,142,202,212,302,312,322,900,902,910,912,922,952)
           IN_ElemType= 'DISP'
       case default
           STOP 'Error: Wrong solution code (51707)'
   end select

  ! Define number of Gauss points (2-noded displacement-based element needs reduced integration).
  select case (IN_NumNodesElem)
    case (2)
      Options%NumGauss=1
    case (3)
      Options%NumGauss=2
  end select

! Set name for output file.
  OutFile(1:14)='./res/'//trim(IN_TestCase)//'_SOL'
  write (OutFile(15:17),'(I3)') Options%Solution

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! set shared design variables for fwd problem
 ! In the original code the shared design variables (see input mode) are set
 ! via the input_setup routine. Here, the input_setup method has been replaced
 ! by the input_wrap method: the shared design variables are updated in the
 ! input module and used as required (without their value being reassigned) by
 ! input_setup_wrap
 !call input_setup (NumElems,OutFile,Options)

 call update_shared_setting(IN_NumNodesElem , IN_ElemType, IN_TestCase, IN_BConds)

 call update_shared_input( IN_BeamLength1, IN_BeamLength2,       &
                       & IN_BeamStiffness, IN_BeamMass,          & ! not used in fwd_code
                       & IN_ExtForce, IN_ExtMomnt,               & ! not used in fwd code
                       & IN_SectWidth, IN_SectHeight,            & ! not used in fwd code
                       & IN_ThetaRoot, IN_ThetaTip,              & ! not used in fwd code
                       & IN_TipMass, IN_TipMassY, IN_TipMassZ,   &
                       & IN_Omega                                ) ! not used in fwd code

 !call input_setup_wrap( NumElems,  OutFile, Options )

 print *, 'part 1 done'
 print *, 'Options%Solution = ', Options%Solution


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Optimiser Input
 call opt_setup(gradmode,solmode,fdmode,NOPTMAX)
 if ( (solmode == 'OPT') .or. (solmode == 'SNS') ) then
     call opt_input_cost(FLAG_COST,FLAG_CONSTR,W_COST,W_CONSTR)
     ! note: cost_utl_build_connectivity is in the opt_cost_utl module.
     call cost_utl_build_connectivity(FLAG_CONSTR,CONN_CONSTR)
     call cost_utl_build_connectivity(FLAG_DESIGN_SHARED,CONN_XSH)

     if (solmode == 'SNS') then        ! this is only needed for the allocation
        NOPTMAX=0
     end if
     allocate(CONSTR( size(CONN_CONSTR),0:NOPTMAX ) ); CONSTR=0.0_8
     allocate(  COST(                   0:NOPTMAX ) ); COST=0.0_8
     allocate(   DCDXSH(                    size(CONN_XSH),0:NOPTMAX ) ); DCDXSH=0.0_8
     allocate( DCONDXSH( size(CONN_CONSTR), size(CONN_XSH),0:NOPTMAX ) ); DCONDXSH=0.0_8
 end if


do NOPT=0,NOPTMAX

    call fwd_problem_prealloc( NumElems,OutFile,Options,    &   ! from input_setup
                 &                        Elem,     &  ! from opt_main_xxx
                 &                    NumNodes,     &  ! from input_elem
                 &  BeamSpanMass, BeamSpanStiffness,&  ! Input added for the optimisation
                 &  BoundConds,PosIni,ForceStatic,PhiNodes, &   ! from input_node
                 &                    OutGrids,    &   ! from pt_main_xxx
                 &                      PsiIni,    &   ! from xbeam_undef_geom
                 &                Node, NumDof,    &   ! from xbeam_undef_dofs
                 & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 &                        Time,    &   ! input_dyn_setup
                 &       ForceTime,ForceDynAmp,    &   ! input_dynforce
                 &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! ! output from sol 202, 212, 302, 312, 322
                 & RefVel, RefVelDot, Quat)

     ! sm 21 aug 2014
     ! all saving moved into python
     !
     ! Store results in text file.
     !open (unit=11,file=OutFile(1:17)//'_def.txt',status='replace')
     !call output_elems (11,Elem,PosDef,PsiDef)
     !close (11)

    ! Compute Design Output
    DensityVector   = Elem(:)%Mass(1,1)
    LengthVector    = Elem(:)%Length

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Loop control
     select case (solmode)

        case default !('FWD')
            print *, 'Forward Problem Completed!'
            exit

        case ('OPT','SNS')
            ! -------------- Evaluate cost and Constraints at current design ---
            COST(NOPT) = cost_global( FLAG_COST, W_COST,  &
                                    & PosIni, PosDef,     &
                                    & Elem(:)%Mass(1,1),Elem(:)%Length         )

            CONSTR(:,NOPT) = cost_constraints( W_CONSTR, CONN_CONSTR, &
                                            & PosIni,PosDef,         &
                                            & Elem(:)%Mass(1,1),Elem(:)%Length )

            ! ----------------------------------------- Sensitivity analysis ---
            print *, 'Sensitivity Analysis started. NOPT=', NOPT
            select case (gradmode)
                case ('FDF')
                    print *, 'Gradients will be computed via Finite Differences'
                    ! interface not updated
                    call fd_main_prealloc( NumElems,OutFile,Options,           &
                                & Elem,                                        &
                                & NumNodes,                                    &
                                & BeamSpanMass, BeamSpanStiffness, &   ! Input added for the optimisation
                                & BoundConds,PosIni,ForceStatic,PhiNodes,      &
                                & OutGrids,                                    &
                                & PsiIni,                                      &
                                & Node, NumDof,                                &
                                & PosDef, PsiDef, InternalForces,              &
                                & Time,                                        &
                                & ForceTime,ForceDynAmp,                       &
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

call toc

 ! IMPORTANT: this gives error if the routine is called more then once from the
 ! python wrapper
 deallocate(XSH)

 end subroutine opt_main

end module opt_routine
