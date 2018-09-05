! **********************************************************************
!                  - KLEIN-GORDON EQUATION SOLVER -
! **********************************************************************
! Master's Thesis
! Author: Ignacio López de Arbina
!
! Descrition: The program solves the Kleing-Gordon equation for an
!             hadronic atom, particularly for a Kaon, using a Coulomb
!             and an optical potential.
!
! Algorithm:  Numerov. Implemented in subroutine "Numerov".
!		  	  The subroutine inetegrates the (radial) Klein-Gordon
!		  	  eq. from the origin forward to the classical turning 
!		  	  point and from the 'infinite' backwards to the classical
!		      turning point, for a given eigenvalue.
!		  	  The program studies the 'matching' of both wave functions
!		      and change the eigenvalue acordingly to perform a new 
!		      iteration until reach convergence.
! **********************************************************************

program nuc_klein_gordon

    use ansi_colors
    use inp_data
    use potentials

    implicit none

    ! external variables
    complex(dp) :: BE_in

    ! internal variables
    real(dp) :: t1, t2
    integer :: part_charge, pot_type
    character(*), parameter :: head_color = c_lblue
    character(*), parameter :: text_color = c_lcyan
    character(*), parameter :: res_color  = c_lcyan
    character(len=15) :: str

    print *, color(repeat('=',60), head_color)
    print *, repeat(' ',14), color('* Klein-Gordon Equation Solver *', head_color), repeat(' ',14)
    print *, color(repeat('=',60), head_color)

    ! asks the user the charge of the particle y/n
72  continue
    print *, color("Select charged or non-charged particle (Ctrl+C to stop):", text_color)
    print *, color(" 1. Charged (K⁻)", text_color)
    print *, color(" 2. Non-charged (K⁰)", text_color)
    read(*,*,err = 72) part_charge

    ! if the input is not correct ask another time
    if (part_charge /= 1.and.part_charge /= 2) then
        print *, color("You shall select 1 or 2!", c_red)
        call sleep (2)
        go to 72
    end if

    ! user choose the potential to use    
73	continue
    print *, color("Which potential would you like to use? (Ctrl+C to stop)", text_color)
    print *, color(" 1. Phenomenological", text_color)
    print *, color(" 2. Phenomenological density dependent",text_color)
    print *, color(" 3. Chiral theory based", text_color)
    print *, color(" 4. Scattering amplitude based (Ramos-Oset)", text_color)
    print *, color(" 5. Scattering amplitude based (Mares et.al.)", text_color)
    print *, color(" 6. No nuclear potential", text_color)
    read(*,*,err = 73) pot_type 
    
    ! if the input is not correct ask another time
    if (pot_type /= 1.and.pot_type /= 2&
        &.and.pot_type /= 3.and.pot_type /= 4&
        &.and.pot_type /= 5.and. pot_type /= 6) then
        print *, "You shall choose between 1, 2, 3, 4, 5 or 6!"
        call sleep(2)
        go to 73
    end if

    call cpu_time (t1)
    call Numerov (part_charge, pot_type, BE_in)
    call cpu_time (t2)

    print *,
    print *, repeat(' ', 1), color(repeat('*', 54), res_color)
    print '(2X,A6,5X,A6,10X,A12,16X,A16)', &
            color('Z', res_color),       &
            color('A', res_color),       &
            color('B (keV)', res_color), &
            color('Width (keV)', res_color)
    print *, repeat(' ', 1), '-', repeat(' ', 5), '-',&
             repeat(' ', 6), repeat('-', 16),          &
             repeat(' ', 7), repeat('-', 18)
    PRINT '(1X,I2,5X,I2,3X,F18.14,7X,F18.14)', Z, A, real(BE_in), 2.0_dp*aimag(BE_in)
    print *,

    print *, repeat(' ', 1), color("CPU computing time =", text_color), t2 - t1, color("s", text_color)
    print *, repeat(' ', 1), color(repeat('*', 54), res_color)

    PRINT *,
    PRINT *, color("===== END OF PROGRAM: Klein-Gordon Equation Solver", head_color)
    PRINT *,

end program nuc_klein_gordon
