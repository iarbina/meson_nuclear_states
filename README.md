# Mesic in-medium nuclear states
 
## Description
    This program solves Klein-Gordon equation describing 
    the bound or quasi-bound states of meson particles
    whitin the nucleus. For that, it integrates different
    types of nuclear potentials (optical potentials-like).

## Conent
### /src/ directory contains the source files:
        - inp_data.f90 (module): contains global variables
          parameters and physical constants. Some of the
          variables are read from files in different 
          subroutines below the 'contains' section.

        - potentials.f90 (module): contains the different 
          functions and subrotuines decribing the nuclear
          potentials as well as the nuclear density.

        - Numerov.f90 (subroutine): contians the numerical
          method (Numerov) to solve KG equation. It calls
          to the different subroutines (in the same script).

        - spls3.f (subroutine): cubic spline interpolation 
          subroutine. It's used by some of the potentials
          and scattering amplitudes.

        - nuc-modular.f90 (main): main program. Just read
          some user control inputs and calls Numerov.f90.

### /obj/ directory contains the .o and .mod files to run the program

### /inp/ directory contains the input data read by the program:
        - /nuclei/ contains nuclei parameters (Z, A, Rn, a, n, l)

        - /potentials/ contains data of the potential.

        - /amplitudes/ contains data of the scattering amplitudes.

### /out/ is an auxiliary directory where to save outputs from the program

### /bin/ directory contains the binaray (executable) file of the program

### /viz/ directory contains data visualization programs (.py) of the output data:
        - wavef.py: plots the squared radial wave function.
        
        - wavef2.py: plots the Real and Imaginary parts
          of the radial wave fucntion.

            -- Note: useful to check the mathcing of both left
               and right integraton of the wave fucntion.
### /python/ directory contains Python programs for data analysis and visualization:
        - /data_smoothing/ contains a .py script for data 
          smoothing of the scattering amplitudes.

        - /data_visualization/ contains .py scripts to visualize
          scattering amplitudes data.

### Makefile
        Script that contains the steps for the correct compilation
        of all the modules, subroutines and programs, in GNU 
        makefile language.

### README.md
        Current script.
