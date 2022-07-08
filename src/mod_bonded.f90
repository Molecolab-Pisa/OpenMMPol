module mod_bonded
    !! Module to handle the bonded part of the FF, it closely follows the 
    !! AMOEBA functional form.
    
    use mod_memory, only: ip, rp
    
    implicit none
    private

    ! Bond section
    integer(ip) :: nbond
    integer(ip), allocatable :: bonda(:), bondb(:)
    real(rp) :: bond_cubic, bond_quartic
    real(rp), allocatable :: kbond(:), l0bond(:)

    public :: bond_init, bond_potential, bonda, bondb, kbond, &
              l0bond, bond_cubic, bond_quartic

    ! Urey-Bradley
    integer(ip) :: nurey
    integer(ip), allocatable ::ureya(:), ureyb(:)
    real(rp) :: urey_cubic, urey_quartic
    real(rp), allocatable :: kurey(:), l0urey(:)
    public :: urey_init, urey_potential, ureya, ureyb, kurey, &
              l0urey, urey_cubic, urey_quartic

    contains

    subroutine bond_init(n) 
        !! Initialize bond stretching arrays

        use mod_memory, only: mallocate

        implicit none

        integer(ip) :: n
        !! Number of bond stretching functions in the potential
        !! energy of the system

        call mallocate('bond_init [bonda]', n, bonda)
        call mallocate('bond_init [bondb]', n, bondb)
        call mallocate('bond_init [kbond]', n, kbond)
        call mallocate('bond_init [l0bond]', n, l0bond)
        nbond = n
        bond_cubic = 0.0_rp
        bond_quartic = 0.0_rp

    end subroutine bond_init

    subroutine bond_potential(V)
        !! Compute the bond-stretching potential

        use mod_constants, only : eps_rp
        use mod_mmpol, only: cmm

        implicit none

        real(rp), intent(inout) :: V
        !! Bond potential, result will be added to V

        integer :: i
        logical :: use_cubic, use_quartic
        real(rp) :: dr(3), l, dl, dl2

        use_cubic = (abs(bond_cubic) > eps_rp)
        use_quartic = (abs(bond_quartic) > eps_rp)

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            do i=1, nbond
                dr = cmm(:,bonda(i)) - cmm(:,bondb(i))
                l = sqrt(dot_product(dr, dr))
                dl = l - l0bond(i)
                V = V + kbond(i) * dl * dl
            end do
        else
            do i=1, nbond
                dr = cmm(:,bonda(i)) - cmm(:,bondb(i))
                l = sqrt(dot_product(dr, dr))
                dl = l - l0bond(i)
                dl2 = dl * dl

                V = V + kbond(i)*dl2 * (1.0_rp + bond_cubic*dl + bond_quartic*dl2)
            end do
        end if
        
    end subroutine bond_potential
    
    
    subroutine urey_init(n) 
        !! Initialize Urey-Bradley potential arrays

        use mod_memory, only: mallocate

        implicit none

        integer(ip) :: n
        !! Number of Urey-Bradley functions in the potential
        !! energy of the system

        call mallocate('urey_init [ureya]', n, ureya)
        call mallocate('urey_init [ureyb]', n, ureyb)
        call mallocate('urey_init [kurey]', n, kurey)
        call mallocate('urey_init [l0urey]', n, l0urey)
        nurey = n
        urey_cubic = 0.0_rp
        urey_quartic = 0.0_rp

    end subroutine urey_init

    subroutine urey_potential(V)
        !! Compute the Urey-Bradley potential

        use mod_constants, only : eps_rp
        use mod_mmpol, only: cmm

        implicit none

        real(rp), intent(inout) :: V
        !! Urey-Bradley potential, result will be added to V

        integer :: i
        logical :: use_cubic, use_quartic
        real(rp) :: dr(3), l, dl, dl2

        use_cubic = (abs(urey_cubic) > eps_rp)
        use_quartic = (abs(urey_quartic) > eps_rp)

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            do i=1, nurey
                dr = cmm(:,ureya(i)) - cmm(:,ureyb(i))
                l = sqrt(dot_product(dr, dr))
                dl = l - l0urey(i)
                V = V + kurey(i) * dl * dl
            end do
        else
            do i=1, nurey
                dr = cmm(:,ureya(i)) - cmm(:,ureyb(i))
                l = sqrt(dot_product(dr, dr))
                dl = l - l0urey(i)
                dl2 = dl * dl

                V = V + kurey(i)*dl2 * (1.0_rp + urey_cubic*dl + urey_quartic*dl2)
            end do
        end if
        
    end subroutine urey_potential

end module mod_bonded
