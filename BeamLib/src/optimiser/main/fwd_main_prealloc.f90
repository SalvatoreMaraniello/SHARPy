!-> Module. - Forward - 10Jan2011 Updated: 10/jul/2014
!
!-> Author.- Henrik Hesse  (h.hesse09@imperial.ac.uk)
!            Rafa Palacios (rpalacio@imperial.ac.uk)
!            Rob Simpson   (rjs10@imperial.ac.uk) copied from main_andrea.f90
!            Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Important:
!   - this module is equivalemnt to fwd_main, but all variables used in the
!     fortran-python interface are not allocatable
!
!-> Description.-
!   - This module contains the xbeam forward solver.
!   - It has been recasted into a module structure from the Main xbeam program
!   (proj_root/BeamLib/src/fortran/main/main_xxx.f90) to allow the solver to
!   be called from the outside.
!   - The original Main has been divided into subroutines to allow a more
!   flexible access from the 'outside'
!   - All printing statements have been kept in here as well.
!
! -> Subroutines:
!   - fwd_static_input: reads input for static problem
!   - pre_static_solver: prepares and allocates variables for static solver;
!   - fwd_static_solver: points to the static solver;
!   - pre_dynamic_solver: reads and allocates variables for dynamic solver;
!   - fwd_dynamic_solver: contains the subroutine version of the xbeam Main program;
!   - pre_coupled_solver: reads and allocates variables for dynamic solver;
!   - fwd_coupled_solver: points to the coupled solver;
!   - The input to this function are passed in the same 'order of  appearence'
!     as in the opt_main_xxx file;
!
! -> Bugs/Reminders:
!  - The pre_solver for dynamic and couple dproblems need to be developed. These
!     require a call to the static pre solver
!  - all allocation statements substituted with conditional allocation (lib_array)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fwd_main_prealloc

 use xbeam_shared
 use xbeam_undef
 use cbeam3_solv
 use xbeam_solv
 use xbeam_perturb
 use lib_out
 use input

 use lib_array ! optimisation

 implicit none

contains

! This subroutine executes all the main program with exeption of the input setup
! routine.
 subroutine fwd_problem_prealloc(NumElems,OutFile,Options,  &   ! from input_setup
                 &                        Elem,    &   ! from opt_main_xxx
                 &                    NumNodes,    &   ! from input_elem
                 & BeamSpanMass, BeamSpanStiffness,&   ! Input added for the optimisation
                 & BoundConds,PosIni,ForceStatic,PhiNodes, &   ! from input_node
                 &                    OutGrids,         &   ! from pt_main_xxx
                 &                      PsiIni,    &   ! from xbeam_undef_geom
                 &                Node, NumDof,    &   ! from xbeam_undef_dofs
                 & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 & Time,                             & ! input_dynsetup
                 & ForceTime, ForceDynAmp, ForceDynamic,& ! input_dynforce
                 & ForcedVel, ForcedVelDot,          & ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! output from sol 202, 212, 302, 312, 322
                 & RefVel, RefVelDot, Quat,  &       ! to be allocated in fwd_pre_coupled_solver
                 & SUCCESS  ) ! optional, for error ecceptions




 ! Input added for the optimisaiton
     real(8)   :: BeamSpanStiffness(NumElems,6,6) ! Element by Element Stiffness matrix
     real(8)   :: BeamSpanMass(NumElems,6,6)      ! Element by Element Mass matrix
     real(8)   :: PhiNodes (NumNodes)             ! Initial twist at grid points.

 ! The following variables appear as allocatable in fwd_main
     real(8) :: ForceStatic (NumNodes,6) ! Applied static nodal forces.
     real(8) :: InternalForces(:,:)  ! Internal force/moments at nodes.
     real(8) :: PosIni   (:,:)    ! Initial nodal Coordinates.
     real(8) :: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
     real(8) :: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
     real(8) :: PsiDef (:,:,:)    ! Current element orientation vectors.

 ! Dynamic Input/Output
    real(8), intent(inout) :: Time(:)           ! Discrete time vector in the dynamic simulation.
    real(8), intent(inout) :: ForceDynamic (:,:,:) ! applied dynamic force
    real(8), intent(inout) :: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
    real(8), intent(inout) :: ForceTime   (:)   ! Time history of the dynamic nodal forces.
    real(8), intent(inout) :: ForcedVel   (:,:) ! Forced velocities at the support.
    real(8), intent(inout) :: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.

    real(8), intent(inout) :: PosDotDef  (:,:)   ! Current nodal position vector.
    real(8), intent(inout) :: PsiDotDef  (:,:,:) ! Current element orientation vectors.
    real(8), intent(inout) :: PosPsiTime (:,:)   ! Position vector/rotation history at beam tip.
    real(8), intent(inout) :: VelocTime  (:,:)   ! History of velocities.
    real(8), intent(inout) :: DynOut     (:,:)   ! Position of all nodes wrt to global frame a for each time step

    logical, intent(inout), optional :: SUCCESS  ! Variable to allow python wrapper to handle ecceptions.
                                                 ! If the solution does not converge, the variable is set to .false.


    character(len=25)        :: OutFile           ! Output file.
    integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
    integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).

    type(xbopts)             :: Options            ! Solution options (structure defined in xbeam_shared).
    type(xbelem), allocatable:: Elem(:)            ! Element information.
    type(xbnode), allocatable:: Node(:)            ! Nodal information.
    integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node

     logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.

     ! Rigid-body variables
     real(8), intent(inout) :: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
     real(8), intent(inout) :: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
     real(8), intent(inout) :: Quat     (:)      ! Quaternions to describe propagation of reference frame a.


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Set Up Input for Static Problem
     ! (stage required also for dynamic and coupled solutions)
     call fwd_static_input( NumElems, OutFile, Options, &   ! from input_setup
                          & Elem,                       &   ! from opt_main_xxx
                          & NumNodes,                   &   ! from input_ele
                          & BeamSpanMass, BeamSpanStiffness, &   ! Input added for the optimisation
                          & BoundConds, PosIni,         &   ! from input_node
                          & OutGrids                    )   ! from pt_main_xxx

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Reads Forward Problem Input and allocates the required variables for the
     ! forward static problem solution
     call fwd_presolver(NumElems,OutFile,Options,     &   ! from input_setup
                    &                        Elem,    &   ! from opt_main_xxx
                    &                    NumNodes,    &   ! from input_elem
                    &  BoundConds,PosIni,PhiNodes,    &   ! from input_node
                    &                      PsiIni,    &   ! from xbeam_undef_geom
                    &                Node, NumDof,    &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces  )  ! allocated in fwd_presolve_static and output of static analysis

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Forward Solution.
     call fwd_solver(NumElems,OutFile,Options,    &   ! from input_setup
                     &                        Elem,    &   ! from opt_main_xxx
                     &                    NumNodes,    &   ! from input_elem
                     &  PosIni,ForceStatic,PhiNodes,    &   ! from input_node
                     &                  OutGrids,      &   ! from pt_main_xxx
                     &                      PsiIni,    &   ! from xbeam_undef_geom
                     &                Node, NumDof,    &   ! from xbeam_undef_dofs
                     & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 & Time,                             & ! input_dynsetup
                 & ForceTime, ForceDynAmp, ForceDynamic,   & ! input_dynforce
                 & ForcedVel, ForcedVelDot,          & ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! output from sol 202, 212, 302, 312, 322
                 & RefVel, RefVelDot, Quat,   &      ! to be allocated in fwd_pre_coupled_solver
                 & SUCCESS ) ! optional, for error ecceptions

 end subroutine fwd_problem_prealloc



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fwd_solver
!
! Select the right pre-process sequence of operations according to the
! sought solution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fwd_solver(NumElems,OutFile,Options,    &   ! from input_setup
                 &                        Elem,    &   ! from opt_main_xxx
                 &                    NumNodes,    &   ! from input_elem
                 & PosIni,ForceStatic,PhiNodes,    &   ! from input_node
                 &                  OutGrids,      &   ! from pt_main_xxx
                 &                      PsiIni,    &   ! from xbeam_undef_geom
                 &                Node, NumDof,    &   ! from xbeam_undef_dofs
                 & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 & Time,                           &   ! input_dynsetup
                 & ForceTime, ForceDynAmp, ForceDynamic, & ! input_dynforce
                 & ForcedVel, ForcedVelDot,        &   ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! output from sol 202, 212, 302, 312, 322
                 & RefVel, RefVelDot, Quat,        &   ! to be allocated in fwd_pre_coupled_solver
                 & SUCCESS  ) ! optional, for error ecceptions

    ! The following variables appear as allocatable in fwd_main
    real(8) :: ForceStatic (:,:) ! Applied static nodal forces.
    real(8) :: InternalForces(:,:)  ! Internal force/moments at nodes.
    real(8) :: PosIni   (:,:)    ! Initial nodal Coordinates.
    real(8) :: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
    real(8) :: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
    real(8) :: PsiDef (:,:,:)    ! Current element orientation vectors.
    real(8) :: PhiNodes (:)      ! Initial twist at grid points.

    ! Dynamic Input/Output
    real(8), intent(inout) :: Time(:)           ! Discrete time vector in the dynamic simulation.
    real(8), intent(inout) :: ForceDynamic (:,:,:) ! sm: applied dynamic force
    real(8), intent(inout) :: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
    real(8), intent(inout) :: ForceTime   (:)   ! Time history of the dynamic nodal forces.
    real(8), intent(inout) :: ForcedVel   (:,:) ! Forced velocities at the support.
    real(8), intent(inout) :: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.

    real(8), intent(inout) :: PosDotDef (:,:)   ! Current nodal position vector.
    real(8), intent(inout) :: PsiDotDef (:,:,:) ! Current element orientation vectors.
    real(8), intent(inout) :: PosPsiTime(:,:)   ! Position vector/rotation history at beam tip.
    real(8), intent(inout) :: VelocTime(:,:)    ! History of velocities.
    real(8), intent(inout) :: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step

    logical, intent(inout), optional :: SUCCESS  ! Variable to allow python wrapper to handle ecceptions.
                                                 ! If the solution does not converge, the variable is set to .false.

    integer:: i,j                                 ! Counter.
    integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
    integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
    integer:: NumSteps                            ! Steps for time simulation

    type(xbopts)             :: Options           ! Solution options (structure defined in xbeam_shared).
    type(xbelem), allocatable:: Elem(:)           ! Element information.
    type(xbnode), allocatable:: Node(:)           ! Nodal information.

    character(len=25)        :: OutFile           ! Output file.

    logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.

    ! Rigid-body variables
    real(8), intent(inout) :: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
    real(8), intent(inout) :: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
    real(8), intent(inout) :: Quat     (:)      ! Quaternions to describe propagation of reference frame a.




    ! ---------------------------------------------------------- Static Solution
    select case (Options%Solution)

        case (102,302)
            ! sm: NumDof is 6*NumberIndepNodes
            ! ForceStatic is a matrix (row: global numbering; columns: forces and Moments)
            ! PsiIni: CRV at the nodes [Psi0(NumElems,MaxElNod,3)]
            call cbeam3_solv_linstatic ( NumDof,Elem,Node,ForceStatic,PosIni,PsiIni, &
                                       & PosDef,PsiDef,Options)

        case (112,142,312,322)
            call cbeam3_solv_nlnstatic ( NumDof,Elem,Node,ForceStatic,PosIni,PsiIni, &
                                       & PosDef,PsiDef,Options)

        case default
            print *, 'No static solution will be run'

     end select


    !-------------------------------------------------------- Vibration analysis
    select case (Options%Solution)

        case (142)
        ! Tangent linear vibration analysis (around the current deformed beam)
            ! sm: allocation of ForcedVel moved inside the select case
            !allocate (ForcedVel(1,6)); ForcedVel   = 0.d0
            call cbeam3_solv_modal ( 12,NumDof,Elem,Node,ForcedVel,PosIni,PsiIni, &
                                   & PosDef,PsiDef,Options)
            !deallocate(ForcedVel)

    end select


    ! --------------------------------------------------------- Dynamic solution
    if ((Options%Solution.ge.200).and.(Options%Solution.le.952)) then

        NumSteps = size(Time)-1

        ! ------------------------------------- Structural dynamic analysis only
        select case (Options%Solution)

            case (202,322)
            ! CBEAM3: Tangent linear dynamic (around the current deformed beam).
                call cbeam3_solv_lindyn ( 12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,         &
                                        & ForceTime,ForcedVel,ForcedVelDot,PosIni,PsiIni,                &
                                        & PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut, &
                                        & OutGrids,Options)

            case (302)
            ! CBEAM3: Linear static + linear dynamic.
                ! PosDef=PosIni !!! commented in orginal code
                ! PsiDef=PsiIni !!! commented in orginal code
                call cbeam3_solv_lindyn ( 12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp*0.d0,         &
                                        & ForceTime*0.d0,ForcedVel*0.d0,ForcedVelDot*0.d0,PosIni,PsiIni, &
                                        & PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut, &
                                        & OutGrids,Options)

            case (212,312)
            ! CBEAM3: Nonlinear dynamic (around the current deformed beam).
                !call cbeam3_solv_nlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp,              &
                !                        & ForceTime,ForcedVel,ForcedVelDot,PosIni,PsiIni,                &
                !                        & PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut, &
                !                        & OutGrids,Options)
                call cbeam3_solv_nlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic,ForceDynamic,             &
                                        & ForcedVel,ForcedVelDot,PosIni,PsiIni,                          &
                                        & PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut, &
                                        & OutGrids,Options)

        end select


        ! --------------------------- Coupled Structural and Rigid-Body Dynamics
        if ((Options%Solution.ge.900).and.(Options%Solution.le.952)) then

            ! Initialize (sm: commented out - moved to python wrapper)
            !RefVel   =ForcedVel       ! RefVel(1,5)=0.5d0
            !RefVelDot=ForcedVelDot
            !Quat     =(/1.d0,0.d0,0.d0,0.d0/)

            !print *, PosDotDef
            !print *, PsiDotDef
            !stop

            select case (Options%Solution)

                case (900)
                ! linear rigid body dynamics only
                    call xbeam_solv_rigidlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp, &
                                               & ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,         &
                                               & PosDef,PsiDef,PosDotDef,PsiDotDef,Options)

                case (910)
                ! nonlinear rigid body dynamics only
                    call xbeam_solv_rigidnlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp, &
                                                & ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,         &
                                                & PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

                case (902)
                ! coupled linear rigid body dynamics and linear structural dynamics
                    call xbeam_solv_coupledlindyn ( 12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp, &
                                                  & ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,         &
                                                  & PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

                case (912)
                ! coupled nonlinear rigid body dynamics and nonlinear structural dynamics
                    call xbeam_solv_couplednlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynamic, &
                                                  & RefVel,RefVelDot,Quat,PosIni,PsiIni,         &
                                                  & PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

                case (932)
                ! sm
                ! coupled nonlinear rigid body dynamics and nonlinear structural dynamics with capability
                ! of static load. The static solution is not run first (e.g. case of shperical joint, where
                ! a static solution does not apply)
                    call xbeam_solv_couplednlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic,ForceDynamic, &
                                                  & RefVel,RefVelDot,Quat,PosIni,PsiIni,         &
                                                  & PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options,SUCCESS)



                case (922)
                ! static and then nonlinear rigid body dynamics and nonlinear structural dynamics
                    call cbeam3_solv_nlnstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni,PosDef,PsiDef,Options)
                    PosIni = PosDef
                    PsiIni = PsiDef
                    call xbeam_solv_couplednlndyn ( 12,NumDof,Time,Elem,Node,ForceStatic,ForceDynamic, &
                                                  & RefVel,RefVelDot,Quat,PosIni,PsiIni,    &
                                                  & PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

                case (952)
                ! coupled linear-elastic dynamics
                    call xbeam_perturb_solv ( 12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,    &
                                            & ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,            &
                                            & PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

            end select

            !!! Store rigid body velocities and accelerations of global reference frame
            !open (unit=11,file=OutFile(1:17)//'_rigid.txt',status='replace')
            !    do i=1,NumSteps+1;  write (11,'(1X,1P7E15.6)') Time(i),RefVel(i,:); end do
            !close (11)
            !
            !open (unit=11,file=OutFile(1:17)//'_vreldot.txt',status='replace')
            !    do i=1,NumSteps+1;  write (11,'(1X,1P7E15.6)') Time(i),RefVelDot(i,:); end do
            !close (11)

        end if ! Coupled analysis

    end if     ! Dynamic analysis

end subroutine fwd_solver



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fwd_presolver
!
! Select the right pre-process sequence of operations according to the
! sought solution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fwd_presolver(NumElems,OutFile,Options,    &   ! from input_setup
                    &                        Elem,    &   ! from opt_main_xxx
                    &                    NumNodes,    &   ! from input_elem
                    &  BoundConds,PosIni,PhiNodes,    &   ! from input_node
                    &                      PsiIni,    &   ! from xbeam_undef_geom
                    &                Node, NumDof,    &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces  )  ! allocated in fwd_presolve_static and output of static analysis

    ! The following variables appear as allocatable in fwd_main
    real(8) :: InternalForces(:,:)  ! Internal force/moments at nodes.
    real(8) :: PosIni   (:,:)    ! Initial nodal Coordinates.
    real(8) :: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
    real(8) :: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
    real(8) :: PsiDef (:,:,:)    ! Current element orientation vectors.
    real(8) :: PhiNodes (:)      ! Initial twist at grid points.

    integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
    integer:: NumSteps                            ! Number of time steps.
    integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
    type(xbopts)             :: Options            ! Solution options (structure defined in xbeam_shared).
    type(xbelem), allocatable:: Elem(:)            ! Element information.
    type(xbnode), allocatable:: Node(:)            ! Nodal information.
    integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node

    character(len=25)        :: OutFile           ! Output file.

    ! ------------------------------------ Compute initial (undeformed) geometry
    !call array3_cond_alloc(PsiIni,NumElems,MaxElNod,3,.true.)
    !allocate(PsiIni(NumElems,MaxElNod,3)); PsiIni=0.d0
    call xbeam_undef_geom ( Elem,PosIni,PhiNodes,PsiIni,Options)
    !    xbeam_undef_geom (inout,    in,      in,   out,     in)
    ! - PosIni (in): Coords in unput_xxx.f90
    ! - PhiNodes (in): pretwist
    ! - PsiIni (out): CRV at the node

    ! Identify nodal degrees of freedom.
    if ( allocated(Node) .eqv. .false. ) then
        allocate (Node(NumNodes))
    end if

    ! sm 25 Oct 2014
    ! Options%Solution added in input as, for spherical joint BCs, the number of
    ! Structural dofs depends on the solution (rigid-body+structural or static) as
    ! well.
    call xbeam_undef_dofs (Elem,BoundConds,  Node, NumDof, Options%Solution)
    ! xbeam_undef_dofs    (  in,        in, inout,    out,      in)
    ! Node: Nodal information [type(xbnode)]
    !NumDof Number of independent degrees of freedom.

    ! Avoid overwriting {Pos/Psi}Def to allow previous solution from optimisation
    ! to be passed in input:
    !if ( max( maxval(abs(PosDef)),maxval(abs(PsiDef)) ) < 1.e-16 ) then
    print *, '[frw_main_prealloc, fwd_problem_prealloc, fwd_presolver]'
    if ( max( maxval(abs(PosDef)),maxval(abs(PsiDef)) ) < epsilon(PosDef(1,1)) ) then
        print *, 'Initial Guess Overwritten with Initial Positions/Rotations!!!'
        print *, '(Max. Elem. < epsilon = ', epsilon(PosDef(1,1))
        PosDef= PosIni
        PsiDef= PsiIni
    else
        print *, 'Initial guess for PosDef and PsiDef found!'
    end if

end subroutine fwd_presolver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fwd_static_input
!
! Calls input function for static forward problem. The routine is a
! pre_requisite for dynamic and coupled solutions as well.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fwd_static_input(NumElems,OutFile,Options, &   ! from input_setup
                           & Elem,                    &   ! from opt_main_xxx
                           & NumNodes,                &   ! from input_elem
                           & BeamSpanMass, BeamSpanStiffness,  & ! Input added for the optimisation
                           & BoundConds, PosIni,      & ! from input_node
                           & OutGrids)                  ! from pt_main_xxx

    ! Input added for the optimisation
     real(8)  :: BeamSpanStiffness(NumElems,6,6) ! Element by Element Stiffness matrix
     real(8)  :: BeamSpanMass(NumElems,6,6)      ! Element by Element Mass matrix
     ! Original Input
     real(8)  :: PosIni(:,:)                     ! Initial nodal Coordinates.
     integer  :: NumElems,NumNodes               ! Number of elements/nodes in the model.
     type(xbopts)             :: Options         ! Solution options (structure defined in xbeam_shared).
     type(xbelem), allocatable:: Elem(:)         ! Element information.
     integer,      allocatable:: BoundConds(:)   ! =0: no BC; =1: clamped nodes; =-1: free node
     logical,      allocatable:: OutGrids(:)     ! Grid nodes where output is written.
     character(len=25)        :: OutFile         ! Output file.

    ! sm: moved in opt_main Read input data.
    ! call input_setup (NumElems,OutFile,Options)

    ! if required for optimisation/sensitivity analysis
    if (allocated(Elem) .eqv. .false.) then
        allocate (Elem(NumElems)) ! sm: initial orientation stored here
    end if
    call input_elem_span (NumElems,NumNodes,Elem, BeamSpanMass, BeamSpanStiffness)
    ! conditional allocation
    !!call array2_cond_alloc(     PosIni, NumNodes, 3, .true.)

    if ( allocated(BoundConds) .eqv. .false. ) then
        allocate(BoundConds (NumNodes));   BoundConds = 0
    end if

    call input_node (NumNodes,Elem,BoundConds)
    ! sm: in input_xxx.f90
    ! input_node (NumNodes,Elem,BoundConds,Coords)
    ! input_node (      in,  in,       out,   out)
    if ( allocated(OutGrids) .eqv. .false. ) then
        allocate(OutGrids(NumNodes))
    end if

    OutGrids          = .false.
    OutGrids(NumNodes)= .true.

 end subroutine fwd_static_input



end module fwd_main_prealloc
