! **********************************************************************
!                  - KLEIN-GORDON EQUATION SOLVER -
! **********************************************************************
! Master's Thesis
! Author: Ignacio López de Arbina
! Insti.: Univerisity of Barcelona
!
! Descrition: The program solves the Kleing-Gordon equation for an
!             mesic nuclei, particularly for Kaons and etas, using 
!             a Coulomb and a wide range of optical potentials.
!
! Algorithm:  Numerov. Implemented in subroutine "Numerov".
!		  	  The subroutine inetegrates the (radial) Klein-Gordon
!		  	  eq. from the origin forward until the nucleus radius 
!		  	  and from the 'infinite' backwards until the nuclear
!		      radius for a given eigenvalue.
!		  	  Then the subroutine studies the 'matching' of both 
!             wave functions (left and right). If the eigenvalue
!             is not solution of KG eq. there would be no matching
!             between wave functions. A new eigenvalue is determined
!             using a stepest descent method in the complex plane.
!		      the process is iterated until the left and right hand
!             side wave functions match smoothly.
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
    integer :: i, answer
    character(*), parameter :: head_color = c_lblue
    character(*), parameter :: text_color = c_lcyan
    character(*), parameter :: res_color  = c_lcyan
    logical :: all_nuc

    
    ! read the data IMPORTANT
    call read_nuclei_data ()            ! required for any calculation
    call read_chiral_potential ()       ! optional: only chiral potential
    call read_kaon_scatt_ampl ()   ! optional: only scatt potential
    call read_eta_scatt_ampl ()

    
    
    
    ! Print the heading
    print *, color(repeat('*',60), head_color)
    print *, color('|', head_color),    &
             repeat(' ',9),            &
             color('* KN and ηN nuclear quasi-bound states *', head_color), &
             repeat(' ',9), &
             color('|', head_color)
    print *, color(repeat('*',60), head_color)


    ! print nuclei data base storaged
    print *, color('|', head_color),  color(repeat(' ',58), head_color), color('|', head_color)
    write(*,332) color( "Nuclei data base", head_color )
    write(*,333) color( repeat('-', 32), head_color )
    write(*,334) color("Z", head_color),   &
                 color("A", head_color),   &
                 color("Rn",head_color),   &
                 color("a", head_color),   &
                 color("n", head_color),   &
                 color("l", head_color)
    write(*,333) color( repeat('-', 32), head_color )

    do i = 1, nnuclei
        print '(15X, I2, 3X, I3, 3X, F5.3, 3X, F5.3, 3X, I1, 3X, I1)', &
                x1(i), x2(i), x3(i), x4(i), x5(i), x6(i)
    end do
    write(*,333) color( repeat('-', 32), head_color )


    write(*,336) color('|', head_color), color('|', head_color)
    write(*,337) color('|', head_color), color('|', head_color)
    write(*,338) color('|', head_color), color('|', head_color)
    write(*,339) color('|', head_color), color('|', head_color)
    write(*,340) color('|', head_color), color('|', head_color)
    write(*,341) color('|', head_color), color('|', head_color)
    print *, color('|', head_color),  color(repeat(' ',58), head_color), color('|', head_color)
    print *, color(repeat('*',60), head_color) 


332 format(' ', 22X, A21)
333 format(' ', 14X, A37)
334 format(' ', 15X, A6, 5X, A6, 5X, A7, 6X, A6, 5X, A6, 3X, A6)
!335 format('1X, I2, 1X, I3, 1X, F5.3, 6X, F5.3, 6X, I1, I1')
336 format(' ', A6, 14X, 'Z: nuclear charge', 27X, A6)
337 format(' ', A6, 14X, 'A: number of nucleons', 23X, A6)
338 format(' ', A6, 14X, 'Rn: nuclear radius', 26X, A6)
339 format(' ', A6, 14X, 'a: WS surface thickness', 21X, A6)
340 format(' ', A6, 14X, 'n: principal quantum number', 17X, A6)
341 format(' ', A6, 14X, 'l: angular quantum number', 19X, A6)



    ! asks the user the type of the particle (K⁻, K⁰ or η)
    print *,
72  continue
    print *, color("Select type of particle (Ctrl+C to stop):", text_color)

    print *, color("    1. Anti-kaon (K⁻)", text_color)
    print *, color("    2. Nuetral kaon (K⁰)", text_color)
    print *, color("    3. Eta (η)", text_color)
    read(*,*, err = 72) part_type

    ! if the input is not correct ask another time
    if (part_type /= 1 .and. part_type /= 2 .and. part_type /= 3) then
        
        print *, color("You shall type 1, 2 or 3!", c_red)
        call sleep (2)
        go to 72
    
    end if

    ! define the charge of the particle (1 = charge, 2 = no charge)
    ! and the mass
    if (part_type == 1) then
        part_charge = 1
        hmass = Kmass
    else if (part_type == 2) then
        part_charge = 2
        hmass = Kmass
    else if (part_type == 3) then
        part_charge = 2
        hmass = Emass
    else
        print *, color('ERROR:', c_red), 'Particle type not defined'
        STOP
    end if

    ! user choose the potential to use    
73	continue
    print *, color("Which potential would you like to use? (Ctrl+C to stop)", text_color)
    print *, color("    1. Phenomenological", text_color)
    print *, color("    2. Phenomenological density dependent",text_color)
    print *, color("    3. Chiral theory based", text_color)
    print *, color("    4. Scattering amplitude based for K⁻ and K⁰", text_color)
    print *, color("    5. Scattering amplitude based for η", text_color)
    print *, color("    6. No nuclear potential", text_color)
    read(*,*,err = 73) pot_type 
    
    ! if the input is not correct ask another time
    if (pot_type /= 1.and.pot_type /= 2&
        &.and.pot_type /= 3.and.pot_type /= 4&
        &.and.pot_type /= 5.and. pot_type /= 6) then
        
        print *, "ERROR: You shall choose between 1, 2, 3, 4, 5 or 6!"
        call sleep(2)
        go to 73
    
    end if


    ! last question
911 continue
    if (pot_type == 5 .and. part_type == 3) then
        
        print *, color('    Select a scattering amplitude', text_color)
        print *, color('        1. Free scattering amplitude', text_color)
        print *, color('        2. In-medium scattering amplitude', text_color)
        read(*,*, err = 911) eta_amp_type

        if (eta_amp_type /= 1 .and. eta_amp_type /= 2) then
            
            print *, 'ERROR: You shall choose between 1 and 2!'
            call sleep(2)
            go to 911
        
        end if

    end if


    ! ask if the calculation is made for all nucleus or for one
70  continue
    print *, color( "Perform the calculation for all nuclei or just for one?", text_color )
    print *, color( "   1. one (recomended)", text_color )
    print *, color( "   2. all (long calculation)", text_color )
    read(*,*, err = 70) answer
    print *, answer

    if (answer /= 1 .and. answer /= 2) then
        print *, "You have to choose between type 1 or 2"
        call sleep (2)
        go to 70
    else if (answer == 2) then
        i = 1
        go to 74
    end if


    ! ask which nuclei is wanted
    print *,
71  continue
    i = 0
    print *, color("Select one of the nuclei in the data base ordered numerically", text_color)
    print *, color("in increasing values, i.e., 1, 2, 3... (Ctrl+C to stop):", text_color)
    read(*,*, err = 71) i

74  continue
    if (i .eq. 0 .or. i > nnuclei) then
        
        print *, 'Sorry, the nucleus selected is not in the data base'
        call sleep(2)
        go to 71
    
    else

        Z = x1(i); A = x2(i); Rn = x3(i); a0 = x4(i); pqn = x5(i); aqn = x6(i)
        print *, Z, A, Rn, a0, pqn, aqn
    
    end if


    ! call subroutine Numerov to perfom the calculations
    call cpu_time (t1)
    call Numerov (BE_in)
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

    ! Save in file result
    open(unit = 983, file = 'out/solutions.out', &
         status = 'old', position = 'append', action = 'write')
    WRITE(983,*) Z, A, real(BE_in), 2.0_dp*aimag(BE_in), part_charge, pot_type, sqrts_out
    close(unit = 983)


    ! if a calculation for all nuclei was selected repeat the calculation for the next one
    if (answer == 2 .and. i < nnuclei) then
        i = i + 1
        go to 74
    end if


    PRINT *,
    PRINT *, color("===== END OF PROGRAM: Klein-Gordon Equation Solver", head_color)
    PRINT *,

end program nuc_klein_gordon
