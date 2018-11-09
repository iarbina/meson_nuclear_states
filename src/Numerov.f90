subroutine Numerov (BE_out)

    use ansi_colors
    use inp_data
    use potentials

    implicit none

    ! variables to return
    complex(dp), intent(out) :: BE_out

    ! global variables
    complex(dp) :: BE
    real(dp) :: dens

    ! internal variables
    real(dp), parameter :: dx = 1.0E-1_dp
    real(dp), parameter :: delta = 1.0E-3_dp
    real(dp) :: xturn, rho0
    real(dp) :: Er0, Ei0, Er, Ei 
    real(dp) :: G0, G1, Gx, Gy, MG


    ! - Physical parameters required for the calulation:
    ! reduced mass
    mu = hmass*nmass*A/(hmass+nmass*A)

    ! normalization density constant
    rho0 = A/DenSuma(dx)

    ! define match/turn point (in this case the nuclear radius)
    xturn = Rn

    ! make the turning point match with a point in the grid
    xturn = dx*int(xturn/dx)+dx

    ! Steepest descent method to maximize function G:
    ! initial values of the energy
    if (pot_type == 1) then
        Er0 = abs( real(Vopt(0.0_dp, Density(0.0_dp, dx, rho0))))
        Ei0 = abs(aimag(Vopt(0.0_dp, Density(0.0_dp, dx, rho0))))
    else if (pot_type == 2) then
        Er0 = abs( real(DDVopt(0.0_dp, Density(0.0_dp, dx, rho0))))
        Ei0 = abs(aimag(DDVopt(0.0_dp, Density(0.0_dp, dx, rho0))))
    else if (pot_type == 3) then
        Er0 = abs( real(ChVopt(0.0_dp, Density(0.0_dp, dx, rho0))))
        Ei0 = abs(aimag(ChVopt(0.0_dp, Density(0.0_dp, dx, rho0))))
    else if (pot_type == 4) then
        !Er0 = abs( real(SelfEnergy(0._dp, Density(0._dp, dx, rho0), BE, part_charge, pot_type)) )*0.5_dp/hmass
        !Ei0 = abs(aimag(SelfEnergy(0._dp, Density(0._dp, dx, rho0), BE, part_charge, pot_type)) )*0.5_dp/hmass
        !Er0 = 40._dp
        !Ei0 = 40.578666687011719_dp
        Er0 = 30._dp
        Ei0 = 30._dp
    else if (pot_type == 5) then
        Er0 = 14._dp
        Ei0 = 1.5_dp
    end if
    

    ! initial values of the function to maximize
    G1 = G(Er0,Ei0)
    G0 = G1-1.D0
    
    ! steepest descent method in action
    do while (G1 > G0)
        G0 = G1
        Gx = (G(Er0+delta,Ei0)-G0)/delta
        Gy = (G(Er0,Ei0+delta)-G0)/delta
        MG = sqrt(Gx*Gx+Gy*Gy)
        Er  = Er0+delta*Gx/MG
        Ei  = Ei0+delta*Gy/MG
        G1  = G(Er,Ei)
        Er0 = Er
        Ei0 = Ei
    end do

    BE_out = cmplx(Er, Ei)

    ! print the solution
#if (DEBUG >= 1)
    WRITE(*,*) Z, A, Er, 2.D0*Ei, sqrts_out
#endif

contains

    !----------------------------------------------------------------------
    !  Integration of the equation:
    !  - Left/right integration with logarithmic derivatives matching at the
    !    turning point
    !  - Returns the inverse of the difference of the logarithmic derivatives
    !    to be maximized by the steepest descent methods
    !----------------------------------------------------------------------
    real(dp) function G(Ereal, Eimag)

        ! external variables
        real(dp), intent(in) :: Ereal, Eimag 

        ! internal variables
        complex(dp) :: y0, y1, ynew, ynext
        complex(dp) :: ymatch, yprev, scaling
        complex(dp) :: dul, dur
        real(dp) :: x, xend
        integer :: i, imatch

        ! internal variables to store wave function
        integer, parameter :: Nend = 5
        real(dp), dimension(int(xturn/dx)+1) :: x_l, Wf_l, ReWf_l, ImWf_l
        real(dp), dimension(int(Nend*xturn/dx)+1) :: x_r, Wf_r, ReWf_r, ImWf_r
        complex(dp), dimension(int(Nend*xturn/dx)+1) :: CWf_r


        open(99, file = "viz/wf2.dat", status = 'replace', action = 'write')
        open(98, file = "viz/lwf.dat", action = 'write')
        open(97, file = "viz/rwf.dat", action = 'write')

        BE  = cmplx(Ereal, Eimag)

        ! Forward integration: from origin to turn/match point
        ! first two points must be given -> initial conditions

        ! fist point
        i = 1
        x = 0.0_dp
        y1 = 0.0_dp

        Wf_l(i) = abs(y1)*abs(y1)
        ReWf_l(i) = real(y1)
        ImWf_l(i) = aimag(y1)
        x_l(i) = x
        write(99,*) i, x_l(i), Wf_l(i), xturn
        write(98,*) i, x_l(i), ReWf_l(i), ImWf_l(i), xturn
        
        ! second point
        i = i + 1
        x = x + dx
        y0 = dx**(aqn+1)
        
        Wf_l(i) = abs(y0)*abs(y0)
        ReWf_l(i) = real(y0)
        ImWf_l(i) = aimag(y0)
        x_l(i) = x
        write(99,*) i, x_l(i), Wf_l(i), xturn
        write(98,*) i, x_l(i), ReWf_l(i), ImWf_l(i), xturn

        ! define the matching/turn point as a index in the grid
        imatch = int(xturn/dx) + 1

        ! integrate by Numerov's algorith KGE in forward direction
        do while (i < imatch)
            x = x + dx
            ynew = y(x, y0, y1, dx)
            y1 = y0
            y0 = ynew
            i = i + 1
            Wf_l(i) = abs(ynew)**2
            ReWf_l(i) = real(ynew)
            ImWf_l(i) = aimag(ynew)
            x_l(i) = x
            write(99,*) i, x_l(i), Wf_l(i), xturn
            write(98,*) i, x_l(i), ReWf_l(i), ImWf_l(i), xturn
        end do

        ! store wave function values at match and previous-to-match points
        ymatch = ynew
        yprev  = y1
        
! Backwards integration: from 'infinity' to turning/match point
! first two points must be given -> initial conditions
        
        ! define infinity (first point)
        xend = Nend*xturn
        x = xend
        
        i = int(xend/dx)+1
        y1 = exp(-(x/Rn)**2)
        Wf_r(i) = abs(y1)**2
        ReWf_r(i) = real(y1)
        ImWf_r(i) = aimag(y1)
        x_r(i) = x

        ! second point
        i = i - 1
        x = x - dx
        y0 = exp(-(x/Rn)**2)
        Wf_r(i) = abs(y0)**2
        ReWf_r(i) = real(y1)
        ImWf_r(i) = aimag(y1)
        x_r(i) = x


        ! integrate by Numerv's algorith KGE in backward direction
        do while (i > imatch)
            x = x - dx
            ynew = y(x, y0, y1, -dx)
            y1 = y0
            y0 = ynew
            i = i - 1
            Wf_r(i) = abs(ynew)**2
            CWf_r(i) = ynew
            x_r(i) = x
        end do
        
        scaling = ymatch/ynew
        ynext = y1*scaling

        ! store wave function data from backward integrations
        do i = imatch, int(xend/dx)+1
            CWf_r(i) = CWf_r(i)*scaling
            ReWf_r(i) = real(CWf_r(i))
            ImWf_r(i) = aimag(CWf_r(i))
            write(99,*) i, x_r(i), Wf_r(i)*abs(scaling)**2, xturn
            write(97,*) i, x_r(i), ReWf_r(i), ImWf_r(i), xturn
        end do
        

        ! close file units
        close(97)
        close(98)
        close(99)

        ! define logarithmic derivative
        dur = (y1 - ynew)/dx
        dul = (ymatch - yprev)/dx

        ! define function to maximize as inverse of the logarithmic
        ! derivatives difference
        G = dx/abs(dur/ynew - dul/ymatch) 

#if (DEBUG >= 3)
    write(FDEB,*) 'Solution: ', BE, G
#endif

#if (DEBUG >= 1)
        PRINT *, BE, G, sqrts_out
#endif

        return
    end function G

    !----------------------------------------------------------------------
    !  Numerov's wave function
    !----------------------------------------------------------------------
    complex(dp) function y(x, y0, y1, dxx)

        ! external variables
        complex(dp), intent(in) :: y0, y1
        real(dp),    intent(in) :: x, dxx

        y = ((12.0_dp-10.0_dp*f(x))*y0-f(x-dxx)*y1)/f(x-dxx)
        return
    end function y
 
    !----------------------------------------------------------------------
    !  Numerov's function f
    !----------------------------------------------------------------------
    complex(dp) function f(x)

        ! external variables
        real(dp), intent(in) :: x

        ! internal variables
        complex(dp) :: g, SE
        real(dp) :: dens

        dens = Density(x, dx, rho0)
        
        ! choose the corresponding Self-Energy for each potential type
        if (pot_type == 1) then
            SE = 2.0_dp*mu*Vopt(x, dens)
        else if (pot_type == 2) then
            SE = 2.0_dp*mu*DDVopt(x, dens)
        else if (pot_type == 3) then
            SE = 2.0_dp*mu*ChVopt(x, dens)
        else if (pot_type == 4) then
            SE = SelfEnergy(x, dens, BE)
        else if (pot_type == 5) then
            SE = SelfEnergy(x, dens, BE)
        else if (pot_type == 6) then
            SE = 0.0_dp
        else
            print *, color("ERROR: Potential type not defined", c_red)
            stop
        end if

        if (part_charge == 1) then
            g = (mu-BE-Vc(x))**2-mu**2-SE-aqn*(aqn+1.0_dp)/(x/hbarc)**2
        else if (part_charge == 2) then
            g = (mu-BE)**2-mu**2-SE-aqn*(aqn+1.0_dp)/(x/hbarc)**2
        else
            print *, color("ERROR: Particle charge not defined", c_red)
            stop
        end if

        f = 1.0_dp+g*(dx/hbarc)**2/12.0_dp
        return
    end function f


end subroutine Numerov
