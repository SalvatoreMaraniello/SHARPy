!-> Module. - Forward - 10Jan2011 Updated: 10/jul/2014
!
!-> Author.- Henrik Hesse  (h.hesse09@imperial.ac.uk)
!            Rafa Palacios (rpalacio@imperial.ac.uk)
!            Rob Simpson   (rjs10@imperial.ac.uk) copied from main_andrea.f90
!            Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
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
!  - all allocation statemnets substituted with conditional allocation (lib_array)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fwd_main

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
 subroutine fwd_problem(NumElems,OutFile,Options,    &   ! from input_setup
                 &                        Elem,    &   ! from opt_main_xxx
                 &                    NumNodes,    &   ! from input_elem
                 &  BoundConds,PosIni,ForceStatic,PhiNodes,    &   ! from input_node
                 &                  OutGrids,         &   ! from pt_main_xxx
                 &                      PsiIni,    &   ! from xbeam_undef_geom
                 &                Node, NumDof,    &   ! from xbeam_undef_dofs
                 & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 &            NumSteps, t0, dt,    &   ! input_dyn_setup
                 &       ForceDynAmp,ForceTime,    &   ! input_dynforce
                 &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! to be allocated in fwd_dynamic_presolve and out of dynamic analysis
                 & RefVel, RefVelDot, Quat)

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
     real(8),      allocatable:: InternalForces(:,:)  ! Internal force/moments at nodes.
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

     real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
     logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.

     ! Rigid-body variables
     real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
     real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
     real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.
     real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Set Up Input for Static Problem
     ! (stage required also for dynamic and coupled solutions)
     call fwd_static_input( NumElems, OutFile, Options,                &   ! from input_setup
                          &                       Elem,                &   ! from opt_main_xxx
                          &                   NumNodes,                &   ! from input_elem
                          & BoundConds, PosIni, ForceStatic, PhiNodes, &   ! from input_node
                          & OutGrids                                   )   ! from pt_main_xxx

    print *, 'Check on Elem'
    print *, allocated(Elem)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Reads Forward Problem Input and allocates the required variables for the
     ! forward static problem solution
     call fwd_presolver (NumElems,OutFile,Options,    &   ! from input_setup
                    &                        Elem,    &   ! from opt_main_xxx
                    &                    NumNodes,    &   ! from input_elem
                    &  BoundConds,PosIni,PhiNodes,    &   ! from input_node
                    &                      PsiIni,    &   ! from xbeam_undef_geom
                    &                Node, NumDof,    &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                    &            NumSteps, t0, dt,    &   ! input_dyn_setup
                    &       ForceDynAmp,ForceTime,    &   ! input_dynforce
                    &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                    & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! to be allocated in fwd_dynamic_presolve and out of dynamic analysis
                    & RefVel, RefVelDot, Quat)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Forward Solution.
     call fwd_solver(NumElems,OutFile,Options,    &   ! from input_setup
                     &                        Elem,    &   ! from opt_main_xxx
                     &                    NumNodes,    &   ! from input_elem
                     &  BoundConds,PosIni,ForceStatic,PhiNodes,    &   ! from input_node
                     &                  OutGrids,      &   ! from pt_main_xxx
                     &                      PsiIni,    &   ! from xbeam_undef_geom
                     &                Node, NumDof,    &   ! from xbeam_undef_dofs
                     & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                     &            NumSteps, t0, dt,    &   ! input_dyn_setup
                     &       ForceDynAmp,ForceTime,    &   ! input_dynforce
                     &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                     & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! to be allocated in fwd_dynamic_presolve and out of dynamic analysis
                     & RefVel, RefVelDot, Quat)         ! to be allocated in fwd_pre_coupled_solver

 end subroutine fwd_problem






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
                 &  BoundConds,PosIni,ForceStatic,PhiNodes,    &   ! from input_node
                 &                  OutGrids,         &   ! from pt_main_xxx
                 &                      PsiIni,    &   ! from xbeam_undef_geom
                 &                Node, NumDof,    &   ! from xbeam_undef_dofs
                 & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                 &            NumSteps, t0, dt,    &   ! input_dyn_setup
                 &       ForceDynAmp,ForceTime,    &   ! input_dynforce
                 &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                 & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! to be allocated in fwd_dynamic_presolve and out of dynamic analysis
                 & RefVel, RefVelDot, Quat)         ! to be allocated in fwd_pre_coupled_solver


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
 real(8),      allocatable:: InternalForces(:,:)  ! Internal force/moments at nodes.
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

 real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
 logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.

 ! Rigid-body variables
 real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
 real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
 real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.
 real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step


select case (Options%Solution)

    case (102,112,142)
        call fwd_static_solver(NumElems,OutFile,Options,  &   ! from input_setup
                     &                      Elem,         &   ! from opt_main_xxx
                     &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,ForceStatic,PhiNodes, & ! from input_node
                     &                  OutGrids,         &   ! from pt_main_xxx
                     &                    PsiIni,         &   ! from xbeam_undef_geom
                     &              Node, NumDof,         &   ! from xbeam_undef_dofs
                     & PosDef, PsiDef, InternalForces     )   ! OUTPUTS!!!

    case (202,212,302,312,322)
        !call fwd_dynamic_presolver
        print *, 'to be implemented...'

    case (900,910,902,912,922,952)
        !call fwd_coupled_presolver
        print *, 'to be implemented...'

    case default
        stop 'Chec the Options%Solution!!!'

    end select

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
                    & PosDef, PsiDef, InternalForces, &   ! allocated in fwd_presolve_static and output of static analysis
                    &            NumSteps, t0, dt,    &   ! input_dyn_setup
                    &       ForceDynAmp,ForceTime,    &   ! input_dynforce
                    &      ForcedVel,ForcedVelDot,    &   ! input_forcedvel
                    & PosDotDef, PsiDotDef, PosPsiTime, VelocTime, DynOut, & ! to be allocated in fwd_dynamic_presolve and out of dynamic analysis
                    & RefVel, RefVelDot, Quat)         ! to be allocated in fwd_pre_coupled_solver

 real(8):: t0,dt                               ! Initial time and time step.

 integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
 integer:: NumSteps                            ! Number of time steps.
 integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
 type(xbopts)             :: Options            ! Solution options (structure defined in xbeam_shared).
 type(xbelem), allocatable:: Elem(:)            ! Element information.
 type(xbnode), allocatable:: Node(:)            ! Nodal information.
 integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node

 real(8),      allocatable:: ForceDynAmp (:,:) ! Amplituif (allocated(Elem) .eqv. .false.) thende of the applied dynamic nodal forces.
 real(8),      allocatable:: ForceTime   (:)   ! Time history of the dynamic nodal forces.
 real(8),      allocatable:: ForcedVel   (:,:) ! Forced velocities at the support.
 real(8),      allocatable:: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.
 real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.
 real(8),      allocatable:: InternalForces(:,:)  ! Internal force/moments at nodes.
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
 ! real(8),      allocatable:: Time(:)           ! Discrete time vector in the dynamic simulation.

 ! Rigid-body variables
 real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
 real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
 real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.
 real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step


select case (Options%Solution)

     ! Reads Forward Problem Input and allocates the required variables for the
     ! forward static problem solution (102,112) and vibrational analysis (142)
     case (102,112,142)
         call fwd_static_presolver(NumElems,OutFile,Options, &   ! from input_setup
              &                      Elem,              &   ! from opt_main_xxx
              &                  NumNodes,              &   ! from input_elem
              & BoundConds, PosIni, PhiNodes,           &   ! from input_node
              &                    PsiIni,              &   ! from xbeam_undef_geom
              &              Node, NumDof,              &   ! from xbeam_undef_dofs
              & PosDef, PsiDef, InternalForces          )

    case (202,212,302,312,322)
        !call fwd_dynamic_presolver
        print *, 'to be implemented...'

    case (900,910,902,912,922,952)
        !call fwd_coupled_presolver
        print *, 'to be implemented...'

    case default
        stop 'Chec the Options%Solution!!!'

    end select

end subroutine fwd_presolver




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fwd_static_input
!
! Calls input function for static forward problem. The routine is a
! pre_requisite for dynamic and coupled solutions as well.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fwd_static_input(NumElems,OutFile,Options, &      ! from input_setup
                    &                      Elem,         &   ! from opt_main_xxx
                    &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,ForceStatic,PhiNodes, & ! from input_node
                    &                  OutGrids)             ! from pt_main_xxx

 ! Interface
 integer :: NumElems,NumNodes                  ! Number of elements/nodes in the model.
 type(xbopts)             :: Options           ! Solution options (structure defined in xbeam_shared).
 type(xbelem), allocatable:: Elem(:)           ! Element information.
 integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node
 real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
 logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.
 character(len=25)        :: OutFile           ! Output file.
 real(8),      allocatable:: PosIni   (:,:)    ! Initial nodal Coordinates.
 real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! sm: moved in opt_main Read input data.
 !!!call input_setup (NumElems,OutFile,Options)

 ! if required for optimisation/sensitivity analysis
 if (allocated(Elem) .eqv. .false.) then
    allocate (Elem(NumElems)) ! sm: initial orientation stored here
 end if
 call input_elem (NumElems,NumNodes,Elem)

 ! conditional allocation
 call array2_cond_alloc(     PosIni, NumNodes, 3, .true.)
 call array2_cond_alloc(ForceStatic, NumNodes, 6, .true.)
 call array1_cond_alloc(   PhiNodes, NumNodes   , .true.)
 if ( allocated(BoundConds) .eqv. .false. ) then
    allocate(BoundConds (NumNodes));   BoundConds = 0
 end if
 ! original code
 !allocate(PosIni     (NumNodes,3)); PosIni     = 0.d0 ! sm: coord in input_xxx.f90
 !allocate(ForceStatic(NumNodes,6)); ForceStatic= 0.d0
 !allocate(PhiNodes   (NumNodes));   PhiNodes   = 0.d0
 !allocate(BoundConds (NumNodes));   BoundConds = 0

 call input_node (NumNodes,Elem,BoundConds,PosIni,ForceStatic,PhiNodes)
 ! sm: in input_xxx.f90
 !    input_node (NumNodes,Elem,BoundConds,Coords,     Forces,PhiNodes)
 !    input_node (      in,  in,       out,   out,        out,     out)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Open main output file and select grid points where output will be written.
 ! Unit 12 is opened here but subprocesses can access it - not preferred solution
 open (unit=12,file=OutFile(1:11)//'.mrb',status='replace')

 if ( allocated(OutGrids) .eqv. .false. ) then
     allocate(OutGrids(NumNodes))
 end if

 OutGrids          = .false.
 OutGrids(NumNodes)= .true.
 call out_title (12,'GLOBAL CONSTANTS IN THE MODEL:')
 write (12,'(14X,A,I12)')    'Number of Beam DOFs:    ', 6
 call out_title (12,'OUTPUT OPTIONS:')
 write (12,'(14X,A,I12)')    'Number of Output Nodes: ', 1
 write (12,'(14X,A,I12)')    'Print Displacements:    ', 1
 write (12,'(14X,A,I12)')    'Print Velocities:       ', 1
 close(12)

end subroutine fwd_static_input



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fwd_static_pre_solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine fwd_static_presolver(NumElems,OutFile,Options, & ! from input_setup
                    &                      Elem,         &   ! from opt_main_xxx
                    &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,PhiNodes,        & ! from input_node
                    &                    PsiIni,         &   ! from xbeam_undef_geom
                    &              Node, NumDof,         &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces     )

 ! Interface
 integer :: NumElems,NumNodes                  ! Number of elements/nodes in the model.
 integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
 type(xbopts)             :: Options           ! Solution options (structure defined in xbeam_shared).
 type(xbelem), allocatable:: Elem(:)           ! Element information.
 type(xbnode), allocatable:: Node(:)           ! Nodal information.
 integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node

 real(8),      allocatable:: InternalForces(:,:) ! Internal force/moments at nodes.

 character(len=25)        :: OutFile           ! Output file.
 real(8),      allocatable:: PosIni   (:,:)    ! Initial nodal Coordinates.
 real(8),      allocatable:: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
 real(8),      allocatable:: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
 real(8),      allocatable:: PsiDef (:,:,:)    ! Current element orientation vectors.

 real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Compute initial (undeformed) geometry.
 call array3_cond_alloc(PsiIni,NumElems,MaxElNod,3,.true.)
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
 call xbeam_undef_dofs (Elem,BoundConds,  Node,NumDof)
 ! xbeam_undef_dofs    (  in,        in, inout,   out)
 ! Node: Nodal information [type(xbnode)]
 !NumDof Number of independent degrees of freedom.

 call array2_cond_alloc( PosDef,             NumNodes, 3, .true.);  PosDef= PosIni;
 call array3_cond_alloc( PsiDef,   NumElems ,MaxElNod, 3, .true.);  PsiDef= PsiIni;
 call array2_cond_alloc( InternalForces,     NumNodes, 6, .true.)
 !allocate (PosDef(NumNodes,3));          PosDef= PosIni
 !allocate (PsiDef(NumElems,MaxElNod,3)); PsiDef= PsiIni
 !allocate (InternalForces(NumNodes,6));  InternalForces= 0.d0

 ! Store undeformed geometry in external text file.
 open (unit=11,file=OutFile(1:11)//'_und.txt',status='replace')
     call output_elems (11,Elem,PosIni,PsiIni)
 close (11)

end subroutine fwd_static_presolver



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! fwd_static_solver
!
! The routine points to:
! - cbeam3_solv_linstatic
! -
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine fwd_static_solver(NumElems,OutFile,Options,  &   ! from input_setup
                    &                      Elem,         &   ! from opt_main_xxx
                    &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,ForceStatic,PhiNodes, & ! from input_node
                    &                  OutGrids,         &   ! from pt_main_xxx
                    &                    PsiIni,         &   ! from xbeam_undef_geom
                    &              Node, NumDof,         &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces     )   ! OUTPUTS!!!

 ! Interface
 integer :: NumElems,NumNodes                  ! Number of elements/nodes in the model.
 integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
 type(xbopts)             :: Options           ! Solution options (structure defined in xbeam_shared).
 type(xbelem), allocatable:: Elem(:)           ! Element information.
 type(xbnode), allocatable:: Node(:)           ! Nodal information.
 integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node
 real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
 real(8),      allocatable:: InternalForces(:,:) ! Internal force/moments at nodes.
 logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.
 character(len=25)        :: OutFile           ! Output file.
 real(8),      allocatable:: PosIni (:,:)      ! Initial nodal Coordinates.
 real(8),      allocatable:: PsiIni (:,:,:)    ! Initial element orientation vectors (CRV)
 real(8),      allocatable:: PosDef (:,:)      ! Current nodal position vector. (sm: local coordinates)
 real(8),      allocatable:: PsiDef (:,:,:)    ! Current element orientation vectors.
 real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.

 ! Internal Variables
 real(8),      allocatable:: ForcedVel(:,:)  ! Forced velocities at the support.

     select case (Options%Solution)
     case (102,302)
       ! sm: NumDof is 6*NumberIndepNodes
       ! ForceStatic is a matrix (row: global numbering; columns: forces and Moments)
       ! PsiIni: CRV at the nodes [Psi0(NumElems,MaxElNod,3)]
       call cbeam3_solv_linstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni, &
    &                              PosDef,PsiDef,Options)


     case (112,142,312,322)
       call cbeam3_solv_nlnstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni, &
    &                              PosDef,PsiDef,Options)

     case default
       print *, 'No static solution'
     end select


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Vibration analysis.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CBEAM3: Tangent linear vibration analysis (around the current deformed beam).
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      select case (Options%Solution)
      case (142)
         ! sm: allocation of ForcedVel moved insuide the select case
         allocate (ForcedVel(1,6)); ForcedVel   = 0.d0
         call cbeam3_solv_modal (12,NumDof,Elem,Node,ForcedVel,PosIni,PsiIni,     &
    &                            PosDef,PsiDef,Options)
         deallocate(ForcedVel)
      end select

end subroutine fwd_static_solver




 subroutine fwd_dynamic_solver(NumElems,OutFile,Options,         &   ! from input_setup
                    &                      Elem,         &   ! from opt_main_xxx
                    &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,ForceStatic,PhiNodes, & ! from input_node
                    &                  OutGrids,         &   ! from pt_main_xxx
                    &                    PsiIni,         &   ! from xbeam_undef_geom
                    &              Node, NumDof,         &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces     )   ! OUTPUTS!!!

     real(8):: t0,dt                               ! Initial time and time step.
     integer:: i,j                                 ! Counter.
     integer :: NumElems,NumNodes                  ! Number of elements/nodes in the model.
     integer:: NumSteps                            ! Number of time steps.
     integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
     type(xbopts)             :: Options           ! Solution options (structure defined in xbeam_shared).
     type(xbelem), allocatable:: Elem(:)           ! Element information.
     type(xbnode), allocatable:: Node(:)           ! Nodal information.
     integer,      allocatable:: BoundConds(:)     ! =0: no BC; =1: clamped nodes; =-1: free node
     real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
     real(8),      allocatable:: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
     real(8),      allocatable:: ForceTime   (:)   ! Time history of the dynamic nodal forces.
     real(8),      allocatable:: ForcedVel   (:,:) ! Forced velocities at the support.
     real(8),      allocatable:: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.
     real(8),      allocatable::  PhiNodes (:)      ! Initial twist at grid points.
     real(8),      allocatable:: InternalForces(:,:) ! Internal force/moments at nodes.
     logical,      allocatable:: OutGrids(:)        ! Grid nodes where output is written.
     character(len=25)        :: OutFile           ! Output file.

     real(8)                  :: PosIni   (:,:)    ! Initial nodal Coordinates.
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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Dynamic solution.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((Options%Solution.ge.200).and.(Options%Solution.le.952)) then

        ! Input data for transient dynamic solution.
        call input_dynsetup (NumSteps, t0, dt,Options)
        !    input_dynsetup (     out,out,out,  inout)
        ! In Option the damping for the Newmark is assigned.
        allocate (Time(NumSteps+1))
        do i=1,NumSteps+1
          Time(i)=t0+dt*dble(i-1)
        end do

        ! Force or velocity input.
        allocate (ForceTime   (NumSteps+1));   ForceTime   = 0.d0
        allocate (ForceDynAmp (NumNodes,6));   ForceDynAmp = 0.d0
        allocate (ForcedVel   (NumSteps+1,6)); ForcedVel   = 0.d0
        allocate (ForcedVelDot(NumSteps+1,6)); ForcedVelDot= 0.d0

        call input_dynforce  (NumNodes,Time,ForceStatic,ForceDynAmp,ForceTime)
        !    input_dynforce  (      in,  in,         in,        out,      out)
        call input_forcedvel (NumNodes,Time,ForcedVel,ForcedVelDot)
        !    input_forcedvel (      in,  in,      out,         out)

        open (unit=11,file=OutFile(1:11)//'_force.txt',status='replace')
          do i=1,NumSteps
            write (11,'(1X,1P14E13.5)') Time(i), ForceTime(i), ForcedVel(i,:), ForcedVelDot(i,:)
          end do
        close (11)

        allocate (PosDotDef(NumNodes,3));           PosDotDef= 0.d0
        allocate (PsiDotDef(NumElems,MaxElNod,3));  PsiDotDef= 0.d0
        allocate (PosPsiTime(NumSteps+1,6));        PosPsiTime=0.d0
        allocate (VelocTime(NumSteps+1,NumNodes));  VelocTime= 0.d0
        allocate (DynOut((NumSteps+1)*NumNodes,3)); DynOut=0.d0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Structural dynamic analysis only
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        select case (Options%Solution)
        ! CBEAM3: Tangent linear dynamic (around the current deformed beam).
        case (202,322)
          call cbeam3_solv_lindyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,         &
    &                              ForceTime,ForcedVel,ForcedVelDot,PosIni,PsiIni,                &
    &                              PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut, &
    &                              OutGrids,Options)

        ! CBEAM3: Linear static + linear dynamic.
        case (302)
    !      PosDef=PosIni
    !      PsiDef=PsiIni
          call cbeam3_solv_lindyn (12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp*0.d0,             &
    &                              ForceTime*0.d0,ForcedVel*0.d0,ForcedVelDot*0.d0,PosIni,PsiIni,     &
    &                              PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut,     &
    &                              OutGrids,Options)

        ! CBEAM3: Nonlinear dynamic (around the current deformed beam).
        case (212,312)
          call cbeam3_solv_nlndyn (12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp,             &
    &                              ForceTime,ForcedVel,ForcedVelDot,PosIni,PsiIni,                    &
    &                              PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut,     &
    &                              OutGrids,Options)

        end select

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Coupled analysis of structural and rigid-body dynamics.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if ((Options%Solution.ge.900).and.(Options%Solution.le.952)) then

          ! Initialize
          allocate (RefVel   (NumSteps+1,6));      RefVel   =ForcedVel;       ! RefVel(1,5)=0.5d0
          allocate (RefVelDot(NumSteps+1,6));      RefVelDot=ForcedVelDot
          allocate (Quat     (4));                 Quat     =(/1.d0,0.d0,0.d0,0.d0/)

          select case (Options%Solution)
          ! linear rigid body dynamics only
          case (900)
            call xbeam_solv_rigidlndyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,      &
    &                                   ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,              &
    &                                   PosDef,PsiDef,PosDotDef,PsiDotDef,Options)

          ! nonlinear rigid body dynamics only
          case (910)
            call xbeam_solv_rigidnlndyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,     &
    &                                    ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,             &
    &                                    PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

          ! coupled linear rigid body dynamics and linear structural dynamics
          case (902)
            call xbeam_solv_coupledlindyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,   &
    &                                      ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,           &
    &                                      PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

          ! coupled nonlinear rigid body dynamics and nonlinear structural dynamics
          case (912)
            call xbeam_solv_couplednlndyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,   &
    &                                      ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,           &
    &                                      PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

          ! static and then nonlinear rigid body dynamics and nonlinear structural dynamics
          case (922)
            call cbeam3_solv_nlnstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni,PosDef,PsiDef,Options)

            PosIni = PosDef
            PsiIni = PsiDef

            call xbeam_solv_couplednlndyn (12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp,   &
    &                                      ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,      &
    &                                      PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

          ! coupled linear-elastic dynamics
          case (952)
            call xbeam_perturb_solv (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,    &
    &                                ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,            &
    &                                PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

          end select

          ! Store rigid body velocities and accelerations of global reference frame
          open (unit=11,file=OutFile(1:11)//'_rigid.txt',status='replace')
            do i=1,NumSteps+1;  write (11,'(1X,1P7E15.6)') Time(i),RefVel(i,:);  end do
          close (11)

          open (unit=11,file=OutFile(1:11)//'_vreldot.txt',status='replace')
            do i=1,NumSteps+1;  write (11,'(1X,1P7E15.6)') Time(i),RefVelDot(i,:);  end do
          close (11)

        end if      ! Coupled analysis


        ! Store results for general dynamic analysis.
        ! Dynamic response of specified node with respect to global frame a.
        open (unit=11,file=OutFile(1:11)//'_dyn.txt',status='replace')
          do i=1,NumSteps-1;  write (11,'(1X,1P7E15.6)') Time(i),PosPsiTime(i,:);  end do
        close (11)

        open (unit=11,file=OutFile(1:11)//'_vel.txt',status='replace')
          do i=1,NumSteps-1;
            if (Time(i).ge.0.d0) then
              do j=1,NumNodes
                write (11,'(1X,2I8,1PE15.6)') i,j,VelocTime(i,j)
              end do
            end if
          end do
        close (11)

        ! Position vector of every node wrt global frame a at each time step.
        open (unit=11,file=OutFile(1:11)//'_shape.txt',status='replace')
          do i=1,NumSteps-1;
            do j=1,NumNodes
               write (11,'(1X,1P7E15.6)') Time(i), DynOut((i-1)*NumNodes+j,:);
            end do
          end do
        close (11)

        ! Screen output position and CRV of specified node at last time step.
        write (*,'(1P6E12.4)') PosDef(NumNodes,:),PsiDef(NumElems,2,:)

      end if     ! Dynamic analysis

 end subroutine fwd_dynamic_solver
end module fwd_main
