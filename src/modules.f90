module modules

  use read_files

  implicit none

  ! declare double precision parameter
  integer, parameter :: dp = selected_real_kind(15, 307)

  ! physical constants
  real(dp), parameter :: pi = 3.1415926535897932384626D0    ! pi
  real(dp), parameter :: hbarc = 197.326971D0               ! MeV fm
  real(dp), parameter :: bohrr = 5.29177298D4               ! Bohr radius [fm]
  real(dp), parameter :: emass = 5.109989027D-1             ! electron mass [MeV]
  real(dp), parameter :: nmass = 939.D0                     ! average nucleon mass [MeV]

  ! physical quantities as parameters
  real(dp), parameter :: hmass = 493.7D0    ! Mass of the hadron (kaon)
  
contains

  !--------------------------------------------------------------------
  ! Function which integrates Fermi-type distribution
  !--------------------------------------------------------------------
  double precision function DenSuma (dx)

    implicit none

    ! external
    real(dp), intent(in) :: dx
    real(dp) :: Rn, a0
    integer :: Z, A, pqn, aqn

    ! internal
    real(dp) :: suma, x
    integer :: N    

    suma = 0.D0
    N    = 10              ! times the nuclear radius
    x    = 0.D0

    if (A <= 16) then
        ! light nuclei MHO density normalization (analytic)
        suma = (Rn*dsqrt(pi))**3*(1.D0+1.5D0*a0)
        
        !do while (x < N*Rn)
        !  suma = suma + x*x*dx*(1.D0+a0*(x/Rn)**2)*dexp(-x*x/(Rn*Rn))
        !  x = x + dx
        !end do
    else if (A > 16) then
        ! heavy nuclei Wood-Saxon density normalization (numerical integration)
        do while (x < N*Rn)
            suma = suma + x*x*dx / (1.D0 + dexp((x-Rn)/a0))
            x = x + dx
        end do
        suma = 4.D0*pi*suma
    end if
    
    DenSuma = suma
    PRINT *, "Module DenSuma = ", DenSuma 
    PRINT *, "Here = ", Z, A, Rn, a0, pqn, aqn
  end function DenSuma

end module modules
