module read_files

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  ! data of nucleus
  real(dp) :: Rn, a0
  integer :: Z, A, pqn, aqn
  

contains
  !--------------------------------------------------------------------
  !  Subroutine to read file with nucleus information
  !--------------------------------------------------------------------
  subroutine read_file(Z, A, Rn, a0, pqn, aqn)

    implicit none

    ! external

    open(1, file="input.inp")
    read(1,*) Z, A, Rn, a0, pqn, aqn
    close(1)

  end subroutine read_file

end module read_files
