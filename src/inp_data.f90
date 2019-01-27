!***********************************************************************
! All the variable and parameters defined in this module are global
! quantities for the other modules, subroutines and programs
!***********************************************************************
module inp_data

    use ansi_colors

    implicit none

    ! double precision definition
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! debugin file parameter
    integer, parameter :: FDEB = DEBUG_FILE

    ! constants
    real(dp), parameter :: pi = 3.1415926535897932384626_dp
    real(dp), parameter :: hbarc = 197.326971_dp            ! MeV fm

    ! physical quantities
    real(dp), parameter :: rhoc  = 0.16_dp
    real(dp), parameter :: nmass = 939.0_dp
    real(dp), parameter :: Kmass = 493.7_dp
    real(dp), parameter :: Emass = 547.862_dp
    real(dp) :: hmass, mu

    ! user input parameters
    integer :: part_type, part_charge, pot_type, eta_amp_type
    
    ! nuclei variables
    real(dp) :: Rn, a0
    integer  :: Z, A, pqn, aqn
    integer, parameter :: nnuclei = 8
    integer, dimension(nnuclei) :: x1, x2, x5, x6
    real(dp), dimension(nnuclei) :: x3, x4
     
    ! chiral potential variables
    real(dp), allocatable, dimension(:) :: r, rrho, pf, VReal, VImag

    ! sqrt(s) solution
    real(dp) :: sqrts_out

    ! kaons scattering amplitude variables
    real(dp), allocatable, dimension(:) :: rhos
    real(dp), allocatable, dimension(:) :: fsqrts, tpr, tpi, tnr, tni   !TODO: delete
    real(dp), allocatable, dimension(:,:) :: fsqrts2d, tpr2d, tpi2d, tnr2d, tni2d

    ! eta free scattering amplitudes
    real(dp), allocatable, dimension(:) :: Resqrts, ReFetaN, Imsqrts, ImFetaN

contains

    !----------------------------------------------------------------------
    !  Subroutine to read nuclei data from a file
    !----------------------------------------------------------------------
    subroutine read_nuclei_data ()
        
        ! internal variables
        integer :: i
        character(len=*), parameter :: fpath = 'inp/nuclei/'

        open(11, file = fpath//'input.inp', status = 'old', action = 'read')

        do i = 1, nnuclei
            read(11,*) x1(i), x2(i), x3(i), x4(i), x5(i), x6(i)
        end do

        close(11)

#if (DEBUG >= 3)
        print *, color('Nuclei data file loaded successfully', c_pink)
#endif
    end subroutine read_nuclei_data

    !----------------------------------------------------------------------
    !  Subroutine to read Chiral potential from a file 
    !----------------------------------------------------------------------
    subroutine read_chiral_potential ()

        ! internal variables
        integer, parameter :: vpot_fdim = 50
        integer :: i
        character(len=*), parameter :: fpath = 'inp/potentials/'

        ! allocate the variables
        allocate(r(vpot_fdim))
        allocate(rrho(vpot_fdim))
        allocate(pf(vpot_fdim))
        allocate(VReal(vpot_fdim))
        allocate(VImag(vpot_fdim))

        open(12, file = fpath//'vpot.pika', status = 'old', action = 'read')

        ! it is spected to have heading (otherwise comment next line)
        read(12,*) 

        do i = 1, vpot_fdim
            read(12,*) r(i), rrho(i), pf(i), VReal(i), VImag(i)
        end do

        close(12)

#if (DEBUG >= 3)
        print *, color('Chiral potential file loaded successfully', c_pink)
#endif
    end subroutine read_chiral_potential

    !----------------------------------------------------------------------
    !  Subroutine to read scattering amplitude from a file 
    !----------------------------------------------------------------------
    subroutine read_kaon_scatt_ampl ()

        ! internal variables
        integer, parameter :: scatf_rows = 201, scatf_cols = 5
        integer, parameter :: funit = 333
        integer :: i, j
        character(len=*), parameter :: fpath = 'inp/amplitudes/'

        ! save in a vector de different densities for which the
        ! scattering amplitude is computed
        allocate(rhos(scatf_cols))
        rhos = rhoc * (/ 0._dp, 0.25_dp, 0.5_dp, 0.75_dp, 1._dp /)

        ! allocate scattering amplitude file variables
        allocate(fsqrts(scatf_rows))
        allocate(tpr(scatf_rows))
        allocate(tpi(scatf_rows))
        allocate(tnr(scatf_rows))
        allocate(tni(scatf_rows))
        
        ! allocate in 2d the scattering amplitude variables
        allocate(fsqrts2d(scatf_rows, scatf_cols))
        allocate(tpr2d(scatf_rows, scatf_cols))
        allocate(tpi2d(scatf_rows, scatf_cols))
        allocate(tnr2d(scatf_rows, scatf_cols))
        allocate(tni2d(scatf_rows, scatf_cols))

        open(unit = 13, file = fpath//'amplKN10.inp', status = 'old', action = 'read')
        
        read(13,*)
        do i = 1, scatf_rows
            read(13,*) fsqrts(i), tpr(i), tpi(i), tnr(i), tni(i)
        end do

        close(13)
        
        
        open(unit = funit    , file = fpath//'amplKN00.inp',         status = 'old', action = 'read')
        open(unit = funit + 1, file = fpath//'amplKN025_smooth.inp', status = 'old', action = 'read')
        open(unit = funit + 2, file = fpath//'amplKN05_smooth.inp',  status = 'old', action = 'read')
        open(unit = funit + 3, file = fpath//'amplKN075_smooth.inp', status = 'old', action = 'read')
        open(unit = funit + 4, file = fpath//'amplKN10.inp',         status = 'old', action = 'read')


        do j = 1, scatf_cols
            ! it is spected to have heading (otherwise comment next line)
            read(j + funit - 1,*) 
            do i = 1, scatf_rows
                read(j + funit - 1,*) fsqrts2d(i, j), tpr2d(i, j), tpi2d(i, j), tnr2d(i, j), tni2d(i, j)
            end do
        end do

        close(funit)
        close(funit+1)
        close(funit+2)
        close(funit+3)
        close(funit+4)

#if (DEBUG >= 3)
        print *, color('Kaon scattering amplitudes files loaded successfully', c_pink)
#endif
    end subroutine read_kaon_scatt_ampl

    subroutine read_eta_scatt_ampl ()

        ! internal variables
        integer :: i, ios
        integer, parameter :: funit  = 736
        integer :: Resize
        integer :: Imsize
        character(len=*), parameter :: fpath = 'inp/eta_amplitudes/'


        if (eta_amp_type == 1) then

            Resize = 102
            Imsize = 125
            
            open(unit=funit, file=fpath//'CS_model/free_Re_CS.dat', status='old', action='read', iostat=ios)
            open(unit=funit+1, file=fpath//'CS_model/free_Im_CS.dat', status='old', action='read', err=92, iostat=ios)

        else if (eta_amp_type == 2) then

            Resize = 71
            Imsize = 66

            open(unit=funit, file=fpath//'CS_model/inme_Re_CS.dat', status='old', action='read', err=92, iostat=ios)
            open(unit=funit+1, file=fpath//'CS_model/inme_Im_CS.dat', status='old', action='read', err=92, iostat=ios)

        else if (eta_amp_type == 3) then

            Resize = 99
            Imsize = 110

            open(unit=funit, file=fpath//'GW_model/free_Re_GW.dat', status='old', action='read', err=92, iostat=ios)
            open(unit=funit+1, file=fpath//'GW_model/free_Im_GW.dat', status='old', action='read', err=92, iostat=ios)
        else if (eta_amp_type == 4) then

            Resize = 57
            Imsize = 65

            open(unit=funit, file=fpath//'IOV_model/free_Re_IOV.dat', status='old', action='read', err=92, iostat=ios)
            open(unit=funit+1, file=fpath//'IOV_model/free_Im_IOV.dat', status='old', action='read', err=92, iostat=ios)

        else if (eta_amp_type == 5) then

            Resize = 74
            Imsize = 90

            open(unit=funit, file=fpath//'KSW_model/free_Re_KSW.dat', status='old', action='read', err=92, iostat=ios)
            open(unit=funit+1, file=fpath//'KSW_model/free_Im_KSW.dat', status='old', action='read', err=92, iostat=ios)

        else if (eta_amp_type == 6) then

            Resize = 73
            Imsize = 88

            open(unit=funit, file=fpath//'M2_model/free_Re_M2.dat', status='old', action='read', err=92, iostat=ios)
            open(unit=funit+1, file=fpath//'M2_model/free_Im_M2.dat', status='old', action='read', err=92, iostat=ios)

        else

            print *, color( 'Error: Problems opening eta scattering amplitudes files', c_red)
            stop

        end if

92      continue

        if (ios /= 0) then
            print *, 'ERROR: Problem opening eta scattering amplitude files'
            print *, 'ois =', ios
            STOP
        else
#if (DEBUG >= 1)
            print *, color('Eta scattering amplitude files open successfully', c_pink)
#endif
        end if



        allocate(Resqrts(Resize))
        allocate(ReFetaN(Resize))
        allocate(Imsqrts(Imsize))
        allocate(ImFetaN(Imsize))

        read(funit,*)
        do i = 1, Resize
           read(funit,*) Resqrts(i), ReFetaN(i)
        end do

        read(funit+1,*)
        do i = 1, Imsize
           read(funit+1,*) Imsqrts(i), ImFetaN(i)
        end do


        close(funit)
        close(funit+1)

    end subroutine read_eta_scatt_ampl

end module inp_data
