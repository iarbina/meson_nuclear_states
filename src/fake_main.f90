program fake_main

    use inp_data
    use modules

    implicit none

    call read_nuclei_data ()
    call read_chiral_potential ()
    call read_scattering_amplitude ()

    print *, "DenSuma =", DenSuma(1.0E-1_dp)

end program fake_main
