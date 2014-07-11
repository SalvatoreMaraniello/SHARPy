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
 type(xbopts)            :: Options            ! Solution options (structure defined in xbeam_shared).
 type(xbelem),allocatable:: Elem(:)            ! Element information.
 type(xbnode),allocatable:: Node(:)            ! Nodal information.
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
 ! Read input data.
 call input_setup (NumElems,OutFile,Options)
 call opt_setup(gradmode,solmode)

 allocate (Elem(NumElems)) ! sm: initial orientation stored here
 call input_elem (NumElems,NumNodes,Elem)

 allocate(PosIni     (NumNodes,3)); PosIni     = 0.d0 ! sm: coord in input_xxx.f90
 allocate(ForceStatic(NumNodes,6)); ForceStatic= 0.d0
 allocate(PhiNodes   (NumNodes));   PhiNodes   = 0.d0
 allocate(BoundConds (NumNodes));   BoundConds = 0
 call input_node (NumNodes,Elem,BoundConds,PosIni,ForceStatic,PhiNodes)
 ! sm: in input_xxx.f90
 !    input_node (NumNodes,Elem,BoundConds,Coords,     Forces,PhiNodes)
 !    input_node (      in,  in,       out,   out,        out,     out)


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Open main output file and select grid points where output will be written.
 ! Unit 12 is opened here but subprocesses can access it - not preferred solution
 open (unit=12,file=OutFile(1:11)//'.mrb',status='replace')
 allocate(OutGrids(NumNodes))
 OutGrids          = .false.
 OutGrids(NumNodes)= .true.
 call out_title (12,'GLOBAL CONSTANTS IN THE MODEL:')
 write (12,'(14X,A,I12)')    'Number of Beam DOFs:    ', 6
 call out_title (12,'OUTPUT OPTIONS:')
 write (12,'(14X,A,I12)')    'Number of Output Nodes: ', 1
 write (12,'(14X,A,I12)')    'Print Displacements:    ', 1
 write (12,'(14X,A,I12)')    'Print Velocities:       ', 1

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Compute initial (undeformed) geometry.
 allocate(PsiIni(NumElems,MaxElNod,3)); PsiIni=0.d0
 call xbeam_undef_geom ( Elem,PosIni,PhiNodes,PsiIni,Options)
 !    xbeam_undef_geom (inout,    in,      in,   out,     in)
 ! - PosIni (in): Coords in unput_xxx.f90
 ! - PhiNodes (in): pretwist
 ! - PsiIni (out): CRV at the node

 ! Store undeformed geometry in external text file.
 open (unit=11,file=OutFile(1:11)//'_und.txt',status='replace')
     call output_elems (11,Elem,PosIni,PsiIni)
 close (11)

 ! Identify nodal degrees of freedom.
 allocate (Node(NumNodes))
 call xbeam_undef_dofs (Elem,BoundConds,  Node,NumDof)
 ! xbeam_undef_dofs    (  in,        in, inout,   out)
 ! Node: Nodal information [type(xbnode)]
 !NumDof Number of independent degrees of freedom.



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Static solution.
 ! Allocate memory for problem variables.
 !allocate (PosDef(NumNodes,3));          PosDef= PosIni
 !allocate (PsiDef(NumElems,MaxElNod,3)); PsiDef= PsiIni
 !allocate (InternalForces(NumNodes,6));  InternalForces= 0.d0

 call tic()
 call fwd_solver(NumElems,OutFile,Options,        &   ! from input_setup
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

