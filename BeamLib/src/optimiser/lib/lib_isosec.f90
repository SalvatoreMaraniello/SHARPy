!->Module lib_isosec. Salvatore Maraniello. 11/Aug/2014
!
!->Description.-
!
!  This module includes tools for computing mass and stiffness matrices of
!  isoentropic beams cross sections
!
!->Subroutines:
! - getprop: returns properties of material (only alluminium implemented)
! - isorect: returns mass and stiffness matrices of a beam element with full
!   rectangular cross-section.
! - isohollowrect:returns mass and stiffness matrices of a beam element with
!   hollow rectangular cross-section.
! - isoellip: full isentropic elliptical cross-section
! - isohollowellip: isentropic hollow elliptical cross-section
! - isocirc: full isentropic circular cross-section
! - isohollowcirc: isentropic hollow circular cross-section
!
!->Reference:
! - http://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lib_isosec

implicit none

    real(8),parameter:: pi=3.14159265358979

contains

! subroutine getprop
! ------------------------------------------------------------------------------
    ! Return material properties
    !
    ! Remark:
    !   - material names ofther have a letter to allow to add variations of the
    !     same material.
    ! --------------------------------------------------------------------------

    subroutine getprop(material, E, G, rho)

        real(8),          intent(out) :: E        !   E: Young's module [Pa]
        real(8),          intent(out) :: G        !   G: Shear Modulus [Pa]
        real(8),          intent(out) :: rho      ! rho: density [kg/m3]
        real(8)                       :: mu       !  mu: Poisson's ratio
        character(len=5), intent(in)  :: material

        select case(trim(material))

            case ('qiqi1') ! dummy material to reproduce Qiqi stiffer beam

                !!!assuming a full square cross-section of size 0.2 x 0.2 m^2
                !E = 75000.0_8
                !rho = 25.0_8
                !!!assuming a full square cross-section of size 0.1 x 0.1 m^2
                E = 1200000.0_8
                rho = 100.0_8
                mu = 0.36_8
                G = E/( 2.0_8 * (1.0_8+mu) ) ! not tested







            case ('titnA') ! Titanium (as per Pai/Palacios, Cook, Muroa)
                E = 127.d9
                mu = 0.36_8
                G = E/( 2.0_8 * (1.0_8+mu) )
                rho = 4430.0_8

            case ('alumA') ! typical aluminium alloy
                E  = 70.d9
                mu = 0.33_8
                G  = E/( 2.0_8 * (1.0_8+mu) )
                rho = 2700.0_8

            case default ! as per BEND test case - allum in
                E  = 1.0d7
                G  = 1.0d5
                rho= 1.0_8


        end select

    end subroutine getprop



! Subroutine isorect
! ------------------------------------------------------------------------------
    ! Compute Stiffness and Mass matrix properties of full rectanfular beam
    ! cross-section made of isentropic material
    ! Unitary Shear Factor is assumed to get a high shear stiffness
    ! Zero warping is assumed
    ! --------------------------------------------------------------------------
    subroutine isorect(l2, l3, material, M, K)

        real(8),          intent(in)  :: l2, l3     ! length along 2 and 3 axis
        character(len=5), intent(in)  :: material

        real(8), intent(out) :: K(6,6)  ! Beam element stiffness matrix
        real(8), intent(out) :: M(6,6)  ! Beam element mass matrix

        real(8)    :: E, G, rho    ! see getprop
        real(8)    :: J, I2, I3, A ! Polar, intertia moments about 2 and 3 axis and section area

        real(8), parameter :: fact =(1.0_8 / 12.0_8)
        real(8), parameter :: shear_fact = 1.d0!(5.d0/6.d0)

        ! Get material properties
        call getprop(material, E, G, rho)

        ! Geometrical Properties
        A   = l2*l3

        I2  = fact  * l2 * l3**3
        I3  = fact  * l3 * l2**3
        J   = (I2 + I3)           ! fact * A * (l2**2 + l3**2)

        ! Stiffness
        K = 0.0_8
        K(1,1) = E*A
        K(2,2) = shear_fact*G*A; K(3,3) = shear_fact*G*A;
        K(4,4) = G*J
        K(5,5) = E*I2;  K(6,6) = E*I3;

        ! Mass
        M = 0.0_8
        M(1,1) = rho*A;   M(2,2)=M(1,1);   M(3,3)=M(1,1);
        M(4,4) = rho*J
        M(5,5) = rho*I2;  M(6,6) = rho*I3;

    end subroutine isorect



! Subroutine isohollowrect
! ------------------------------------------------------------------------------
    ! Compute Stiffness and Mass matrix properties of a hollow rectanfular beam
    ! cross-section made of isentropic material
    ! Unitary Shear Factor is assumed to get a high shear stiffness
    ! Zero warping is assumed
    ! The tickness is defined such that the internal length of the rectangle is:
    ! l_int = l_out - 2*thickness
    ! --------------------------------------------------------------------------
    subroutine isohollowrect(l2, l3, t2, t3, material, M, K)

        real(8),          intent(in)  :: l2, l3     ! length along 2 and 3 axis
        real(8),          intent(in)  :: t2, t3     ! thickness along 2 and 3 axis
        character(len=5), intent(in)  :: material

        real(8), intent(out) :: K(6,6)  ! Beam element stiffness matrix
        real(8), intent(out) :: M(6,6)  ! Beam element mass matrix

        real(8)    :: Mtmp(6,6), Ktmp(6,6)

        ! Get Outer Rectangle
        call isorect(l2, l3, material, M, K)

        ! Get Inner Rectangle
        call isorect(l2-2.0_8*t2, l3-2.0_8*t3, material, Mtmp, Ktmp)

        ! Compute Final Mass and Stiffness
        M = M - Mtmp
        K = K - Ktmp

    end subroutine isohollowrect



! Subroutine isoellip
! ------------------------------------------------------------------------------
    ! Compute Stiffness and Mass matrix properties of full ellipsic beam
    ! cross-section made of isentropic material
    ! Unitary Shear Factor is assumed to get a high shear stiffness
    ! Zero warping is assumed
    ! --------------------------------------------------------------------------
    subroutine isoellip(l2, l3, material, M, K)

        real(8),          intent(in)  :: l2, l3   ! semilength along 2 and 3 axis
        character(len=5), intent(in)  :: material

        real(8), intent(out) :: K(6,6)  ! Beam element stiffness matrix
        real(8), intent(out) :: M(6,6)  ! Beam element mass matrix

        real(8)    :: E, G, rho    ! see getprop
        real(8)    :: J, I2, I3, A ! Polar, intertia moments about 2 and 3 axis and section area

        real(8), parameter :: fact =(pi / 4.0_8)
        real(8), parameter :: shear_fact = 1.d0!(5.d0/6.d0)

        ! Get material properties
        call getprop(material, E, G, rho)

        ! Geometrical Properties
        A   = pi*l2*l3

        I2  = fact  * l2 * l3**3
        I3  = fact  * l3 * l2**3
        J   = (I2 + I3)

        ! Stiffness
        K = 0.0_8
        K(1,1) = E*A
        K(2,2) = shear_fact*G*A; K(3,3) = shear_fact*G*A;
        K(4,4) = G*J
        K(5,5) = E*I2;  K(6,6) = E*I3;

        ! Mass
        M = 0.0_8
        M(1,1) = rho*A;   M(2,2)=M(1,1);   M(3,3)=M(1,1);
        M(4,4) = rho*J
        M(5,5) = rho*I2;  M(6,6) = rho*I3;

    end subroutine isoellip



! Subroutine isohollowellip
! ------------------------------------------------------------------------------
    ! Compute Stiffness and Mass matrix properties of hollow ellipsic beam
    ! cross-section made of isentropic material
    ! Unitary Shear Factor is assumed to get a high shear stiffness
    ! Zero warping is assumed
    ! The tickness is defined such that the internal semilength of the ellipse is:
    ! l_int = l_out - thickness
    ! --------------------------------------------------------------------------
    subroutine isohollowellip(l2, l3, t2, t3, material, M, K)

        real(8),          intent(in)  :: l2, l3     ! semilength along 2 and 3 axis
        real(8),          intent(in)  :: t2, t3     ! thickness along 2 and 3 axis
        character(len=5), intent(in)  :: material

        real(8), intent(out) :: K(6,6)  ! Beam element stiffness matrix
        real(8), intent(out) :: M(6,6)  ! Beam element mass matrix

        real(8)    :: Mtmp(6,6), Ktmp(6,6)

        ! Get Outer Rectangle
        call isoellip(l2, l3, material, M, K)

        ! Get Inner Rectangle
        call isoellip(l2-t2, l3-t3, material, Mtmp, Ktmp)

        ! Compute Final Mass and Stiffness
        M = M - Mtmp
        K = K - Ktmp

    end subroutine isohollowellip



! Subroutine isocirc
! ------------------------------------------------------------------------------
    ! Compute Stiffness and Mass matrix properties of full circular beam
    ! cross-section made of isentropic material
    ! Unitary Shear Factor is assumed to get a high shear stiffness
    ! Zero warping is assumed
    ! --------------------------------------------------------------------------
    subroutine isocirc(r, material, M, K)

        real(8),          intent(in)  :: r   ! cross-section radious
        character(len=5), intent(in)  :: material

        real(8), intent(out) :: K(6,6)  ! Beam element stiffness matrix
        real(8), intent(out) :: M(6,6)  ! Beam element mass matrix

        real(8)    :: E, G, rho    ! see getprop
        real(8)    :: J, I2, I3, A ! Polar, intertia moments about 2 and 3 axis and section area

        real(8), parameter :: fact =(pi / 4.0_8)
        real(8), parameter :: shear_fact = 1.d0!(5.d0/6.d0)

        ! Get material properties
        call getprop(material, E, G, rho)

        call isoellip(r, r, material, M, K)

    end subroutine isocirc


! Subroutine isohollowcirc
! ------------------------------------------------------------------------------
    ! Compute Stiffness and Mass matrix properties of hollow circular beam
    ! cross-section made of isentropic material
    ! Unitary Shear Factor is assumed to get a high shear stiffness
    ! Zero warping is assumed
    ! The tickness is defined such that the internal radious is:
    ! r_int = r_out - thickness
    ! --------------------------------------------------------------------------
    subroutine isohollowcirc(r, t, material, M, K)

        real(8),          intent(in)  :: r     ! radious
        real(8),          intent(in)  :: t     ! thickness
        character(len=5), intent(in)  :: material

        real(8), intent(out) :: K(6,6)  ! Beam element stiffness matrix
        real(8), intent(out) :: M(6,6)  ! Beam element mass matrix

        real(8)    :: Mtmp(6,6), Ktmp(6,6)

        ! Get Outer Rectangle
        call isocirc(r, material, M, K)

        ! Get Inner Rectangle
        call isocirc(r-t, material, Mtmp, Ktmp)

        ! Compute Final Mass and Stiffness
        M = M - Mtmp
        K = K - Ktmp

    end subroutine isohollowcirc


end module lib_isosec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module lib_isosec_test

 use lib_isosec

 implicit none

    integer :: ii

 contains


! ------------------------------------------------------------------------------
    subroutine getprop_test

        real(8)             :: E, G, rho
        character(len=5)    :: material

        material = 'alumA'

        call getprop(material, E, G, rho)

        print *, 'material: ', material
        print *, 'E ', 'G ', 'rho'
        print '(F20.4)', E, G, rho

    end subroutine getprop_test


    subroutine isorect_test

        real(8)           :: l2, l3     ! length along 2 and 3 axis
        character(len=5)  :: material

        real(8) :: K(6,6)  ! Beam element stiffness matrix
        real(8) :: M(6,6)  ! Beam element mass matrix

        print *, '------------------------------------------------- Test isrect'

        material='alumA'
        l2=0.2_8; l3=0.5_8;

        call isorect(l2, l3, material, M, K)

        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

    end subroutine isorect_test


! ------------------------------------------------------------------------------
    subroutine isohollowrect_test

        real(8)           :: l2, l3     ! length along 2 and 3 axis
        real(8)           :: t2, t3     ! thickness along 2 and 3 axis
        character(len=5)  :: material

        real(8) :: K(6,6)  ! Beam element stiffness matrix
        real(8) :: M(6,6)  ! Beam element mass matrix


        print *, '------------------------------------------ Test isohollowrect'

        material='alumA'

        l2=0.2_8;      l3=0.5_8;
        t2=0.5_8 * l2; t3=0.5_8 * l3
        call isohollowrect(l2, l3, t2, t3, material, M, K)

        print *, 'Expected full...'
        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

        l2=0.2_8; l3=0.5_8;
        t2=0.0_8; t3=0.0_8
        call isohollowrect(l2, l3, t2, t3, material, M, K)

        print *, 'Expected zero...'
        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

    end subroutine isohollowrect_test


! ------------------------------------------------------------------------------
    subroutine isoellip_test

        real(8)           :: l2, l3     ! length along 2 and 3 axis
        character(len=5)  :: material

        real(8) :: K(6,6)  ! Beam element stiffness matrix
        real(8) :: M(6,6)  ! Beam element mass matrix

        print *, '------------------------------------------------- Test isellip'

        material='alumA'
        l2=0.1_8; l3=0.25_8;

        call isoellip(l2, l3, material, M, K)

        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

    end subroutine isoellip_test


! ------------------------------------------------------------------------------
    subroutine isohollowellip_test

        real(8)           :: l2, l3     ! length along 2 and 3 axis
        real(8)           :: t2, t3     ! thickness along 2 and 3 axis
        character(len=5)  :: material

        real(8) :: K(6,6)  ! Beam element stiffness matrix
        real(8) :: M(6,6)  ! Beam element mass matrix


        print *, '------------------------------------------ Test isohollowellip'

        material='alumA'

        l2=0.1_8;      l3=0.25_8;
        t2= l2; t3= l3
        call isohollowellip(l2, l3, t2, t3, material, M, K)

        print *, 'Expected full...'
        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

        l2=0.2_8; l3=0.5_8;
        t2=0.0_8; t3=0.0_8
        call isohollowellip(l2, l3, t2, t3, material, M, K)

        print *, 'Expected zero...'
        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

    end subroutine isohollowellip_test


! ------------------------------------------------------------------------------
    subroutine isocirc_test

        real(8)           :: r
        character(len=5)  :: material

        real(8) :: K(6,6)  ! Beam element stiffness matrix
        real(8) :: M(6,6)  ! Beam element mass matrix

        print *, '------------------------------------------------- Test isocirc'

        material='alumA'
        r=0.2;

        call isocirc(r, material, M, K)

        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

    end subroutine isocirc_test


! ------------------------------------------------------------------------------
    subroutine isohollowcirc_test

        real(8)           :: r     ! radious
        real(8)           :: t     ! thickness
        character(len=5)  :: material

        real(8) :: K(6,6)  ! Beam element stiffness matrix
        real(8) :: M(6,6)  ! Beam element mass matrix


        print *, '------------------------------------------ Test isohollowcirc'

        material='bend'!'alumA'

        r=0.2_8
        t= r;
        call isohollowcirc(r, t, material, M, K)

        print *, 'Expected full...'
        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

        r=0.2_8
        t= 0.0_8;
        call isohollowcirc(r, t , material, M, K)

        print *, 'Expected zero...'
        print *, 'Mass Diagonal'
        call print_diagonal(M)
        print *, 'Stiffness Diagonal'
        call print_diagonal(K)

    end subroutine isohollowcirc_test



! ------------------------------------------------------------------------------
    subroutine print_diagonal(A)

        real(8) :: A(:,:)
        integer :: N

        N=size(A(:,1))

        do ii=1,N
             print '(F16.4)', A(ii,ii)
        end do

    end subroutine print_diagonal


end module lib_isosec_test

