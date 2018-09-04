module inp_data

    implicit none

    ! double precision definition
    integer, parameter :: dp = selected_real_kind(15, 307)
    !integer, parameter :: dp = selected_real_kind(33, 4931)

    ! constants
    real(dp), parameter :: pi = 3.1415926535897932384626_dp
    real(dp), parameter :: hbarc = 197.326971_dp            ! MeV fm

    ! physical quantities
    real(dp), parameter :: nmass = 939.0_dp
    real(dp), parameter :: hmass = 493.7_dp
    real(dp), parameter :: rhoc  = 0.16_dp
    real(dp) :: mu

    ! nuclei variables
    real(dp) :: Rn, a0
    integer :: Z, A, pqn, aqn
     
    ! chiral potential variables
    real(dp), allocatable, dimension(:) :: r, rrho, pf, VReal, VImag

    ! scattering amplitude variables
    real(dp), allocatable, dimension(:) :: fsqrts, tpr, tpi, tnr, tni

contains

    !----------------------------------------------------------------------
    !  Subroutine to read nuclei data from a file
    !----------------------------------------------------------------------
    subroutine read_nuclei_data ()

        open(11, status = 'old', action = 'read', file = 'inp/input.inp')

        read(11,*) Z, A, Rn, a0, pqn, aqn

        close(11)

    end subroutine read_nuclei_data

    !----------------------------------------------------------------------
    !  Subroutine to read Chiral potential from a file 
    !----------------------------------------------------------------------
    subroutine read_chiral_potential ()

        ! internal variables
        integer, parameter :: vpot_fdim = 50
        integer :: i

        ! allocate the variables
        allocate(r(vpot_fdim))
        allocate(rrho(vpot_fdim))
        allocate(pf(vpot_fdim))
        allocate(VReal(vpot_fdim))
        allocate(VImag(vpot_fdim))

        open(12, status = 'old', action = 'read', file = 'inp/vpot.pika')

        ! it is spected to have heading (otherwise comment next line)
        read(12,*) 

        do i = 1, vpot_fdim
            read(12,*) r(i), rrho(i), pf(i), VReal(i), VImag(i)
        end do

        close(12)

    end subroutine read_chiral_potential

    !----------------------------------------------------------------------
    !  Subroutine to read scattering amplitude from a file 
    !----------------------------------------------------------------------
    subroutine read_scattering_amplitude ()

        ! internal variables
        integer, parameter :: scat_fdim = 201
        integer :: i

        ! allocate the variables
        allocate(fsqrts(scat_fdim))
        allocate(tpr(scat_fdim))
        allocate(tpi(scat_fdim))
        allocate(tnr(scat_fdim))
        allocate(tni(scat_fdim))

        open(13, status = 'old', action = 'read', file = 'inp/amplKN10.dat')

        ! it is spected to have heading (otherwise comment next line)
        read(13,*) 
        
        do i = 1, scat_fdim
            read(13,*) fsqrts(i), tpr(i), tpi(i), tnr(i), tni(i)
        end do

        close(13)

    end subroutine read_scattering_amplitude

end module inp_data
