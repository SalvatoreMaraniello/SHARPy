!-> Program.- MAIN - 06Jan2011 Updated: 22/21/2012
!
!-> Author.- Henrik Hesse  (h.hesse09@imperial.ac.uk)
!            Rafa Palacios (rpalacio@imperial.ac.uk)
!            Rob Simpson   (rjs10@imperial.ac.uk) copied from main_andrea.f90
!            Salvatore Maraniello (salvatore.maraniello10@imperial.ac.uk)
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!
!   Main programme for the core routines of the multibeam assembly.
!
!
! -> Developer Notes:
! a. in function call, in & out are intended as:
!     - in : goes inside the function
!     - out: is returned from the function
! b. Memory usage:
!     - variables for the forward problem are allocated in the fed_main module.
!     - this allowed a more clear division of the main into modules
!
! -> Bugs/Iusses:
! a. fwd_main to be moved in src once all input/output modifications are completed
! b. unit(12) is opened here but is modifiend into the solvers. The solvers can,
!    however, access the file unless the unit is closed (somewhere in the code).
!    Changes to file will not be permanent until the unit is closed.
!    To be improved.
! c. NumElemes seems to be overwritten by input_elem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program opt_main

 use xbeam_shared                              ! Forward Solution Modules
 use xbeam_undef
 use cbeam3_solv
 use xbeam_solv
 use xbeam_perturb
 use input
 use lib_out
 use opt_input                                 ! Optimisation Modules
 use fwd_main
 use lib_perf

 implicit none

 real(8):: t0,dt                               ! Initial time and time step.
 integer:: i,j                                 ! Counter.
 integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
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

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call fwd_static_input(NumElems,OutFile,Options, &      ! from input_setup
                    &                      Elem,         &   ! from opt_main_xxx
                    &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,ForceStatic,PhiNodes, & ! from input_node
                    &                  OutGrids)   ! from pt_main_xxx


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Optimiser Input
 call opt_setup(gradmode,solmode)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Define a Pointer to Function:
 ! This is required to access the appropriate pre_solve and solve routines!



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Reads Forward Problem Input and allocates the required variables for the
 ! forward static problem solution
 call fwd_static_presolver(NumElems,OutFile,Options,     &   ! from input_setup
                    &                      Elem,         &   ! from opt_main_xxx
                    &                  NumNodes,         &   ! from input_elem
                    & BoundConds,PosIni,       PhiNodes, &   ! from input_node
                    &                    PsiIni,         &   ! from xbeam_undef_geom
                    &              Node, NumDof,         &   ! from xbeam_undef_dofs
                    & PosDef, PsiDef, InternalForces     )


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Static solution.
 ! Allocate memory for problem variables.
 ! sm note: at this point the arrays have size 1
 ! The allocation can be done in the solver function and the output array will
 ! have whichever size given during the allocation.
 ! Even if the allocation doesn't happen in here (but in the fwd main solver)
 ! variables can be deallocated here!
 !
 ! Moved into fwd_solver or fwd_pre_static_solver
 !allocate (PosDef(NumNodes,3));          PosDef= PosIni
 !allocate (PsiDef(NumElems,MaxElNod,3)); PsiDef= PsiIni
 !allocate (InternalForces(NumNodes,6));  InternalForces= 0.d0

 call tic()
 call fwd_static_solver(NumElems,OutFile,Options,        &   ! from input_setup
             &                      Elem,         &   ! from opt_main_xxx
             &                  NumNodes,         &   ! from input_elem
             & BoundConds,PosIni,ForceStatic,PhiNodes, & ! from input_node
             &                  OutGrids,         &   ! from pt_main_xxx
             &                    PsiIni,         &   ! from xbeam_undef_geom
             &              Node, NumDof,         &   ! from xbeam_undef_dofs
             & PosDef, PsiDef, InternalForces     )   ! OUTPUTS!!!
 call toc()

 ! Store results in text file.
 open (unit=11,file=OutFile(1:11)//'_def.txt',status='replace')
 call output_elems (11,Elem,PosDef,PsiDef)
 close (11)



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Set Next Step
 select case (solmode)

    case ('FWD')


    case ('OPT')


    case ('SNS')
        print *, 'Sensitivity Analysis'

 end select


end program opt_main

