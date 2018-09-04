# meson_nuclear_states

Description:
    This program solves Klein-Gordon equation describing 
    the bound or quasi-bound states of meson particles
    whitin the nucleus. For that, it integrates different
    types of nuclear potentials (optical potentials-like).

Conent:
    /src/ directory contains the source files:
        - inp_data.f90 (module): contains global parameters
          and physical constants. Some of the parameters are
          read from files.

        - potentials.f90 (module): contains the different
          nuclear potentials as well as the nuclear density.

        - Numerov.f90 (subroutine): contians the numerical
          method (Numerov) to solve KG equation.

        - spls3.f (subroutine): cubic spline interpolation 
          subroutine. It's used by some of the potentials.

        - nuc-modular.f90 (main): main program. Just read
          some user control inputs and calls Numerov.f90.

    /viz/ directory contains data visualization programs.py
        - wavef.py: plots the squared radial wave function
        
        - wavef2.py: plots the Real and Imaginary parts
          of the radial wave fucntion.

        Note: useful to check the mathcing of both left
        and right integraton of the wave fucntion.
