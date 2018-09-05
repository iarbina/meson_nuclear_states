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

    use inp_data
    use potentials

    implicit none

    ! external variables
    complex(dp) :: BE_in

    ! internal variables
    real(dp) :: t1, t2
    integer :: part_charge, pot_type

    print *, "----------------------------------------------------------"
    print *, "            * Klein-Gordon Equation Solver *              "
    print *, "----------------------------------------------------------"

    ! asks the user the charge of the particle y/n
72  continue
    print *, "Select charged or non-charged particle (Ctrl+C to stop):"
    print *, " 1. Charged (K⁻)"
    print *, " 2. Non-charged (K⁰)"
    read(*,*,err = 72) part_charge

    ! if the input is not correct ask another time
    if (part_charge /= 1.and.part_charge /= 2) then
        print *, "You shall select 1 or 2!"
        call sleep (2)
        go to 72
    end if

    ! user choose the potential to use    
73	continue
    print *, "Which potential would you like to use? (Ctrl+C to stop)"
    print *, " 1. Phenomenological"
    print *, " 2. Phenomenological density dependent"
    print *, " 3. Chiral theory based"
    print *, " 4. Scattering amplitude based (A.R.)"
    print *, " 5. Scattering amplitude based (Mares et.al.)"
    print *, " 6. No nuclear potential"
    read(*,*,err = 73) pot_type 
    
    ! if the input is not correct ask another time
    if (pot_type /= 1.and.pot_type /= 2&
        &.and.pot_type /= 3.and.pot_type /= 4&
        &.and.pot_type /= 5.and. pot_type /= 6) then
        print *, "You shall choose between 1, 2, 3, 4, 5 or 6!"
        call sleep(2)
        go to 73
    end if

    print '(6X,A1,3X,A2,5X,A11,5X,A11)', "Z", "A", "Ec-En (keV)", "Width (keV)"
    call cpu_time (t1)
    call Numerov (part_charge, pot_type, BE_in)
    call cpu_time (t2)
   
    PRINT *, Z, A, real(BE_in), 2.0_dp*aimag(BE_in)
 
    print *, "CPU computing time = ", t2 - t1, " s"

    PRINT *,
    PRINT *, "===== END OF PROGRAM: Klein-Gordon Equation Solver"
    PRINT *,

end program nuc_klein_gordon
