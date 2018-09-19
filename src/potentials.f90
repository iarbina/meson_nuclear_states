module potentials

    use ansi_colors
    use inp_data

    implicit none

contains

    !----------------------------------------------------------------------
    !  Coulomb Potential
    !----------------------------------------------------------------------
    real(dp) function Vc(x)

        ! external variables
        real(dp), intent(in) :: x

        ! internal variables
        real(dp), parameter :: finest = 137.0359991_dp
        
        if (x.lt.Rn) then
            Vc = -0.5_dp*Z*hbarc*(3.0_dp-(x/Rn)**2)/(Rn*finest)
        else
            Vc = -Z*hbarc/(x*finest)
        end if

        return
    end function Vc

    !----------------------------------------------------------------------
    !  Phenomenological Optical Potential
    !----------------------------------------------------------------------
    complex(dp) function Vopt(x, dens)

        ! external variables
        real(dp), intent(in) :: x, dens

        ! internal variable
        complex(dp), parameter :: b0 = (0.52_dp, 0.8_dp)

        Vopt = -2.0_dp*pi*(1.0_dp+mu/nmass)*b0*dens*hbarc**2/mu

        return
    end function Vopt
    
    !----------------------------------------------------------------------
    !  Density Dependent Phenomenological Optical Potential
    !----------------------------------------------------------------------
    complex(dp) function DDVopt(x, dens)

        ! external variables
        real(dp), intent(in) :: x, dens

        ! internal variables
        complex(dp), parameter :: b0exp = (-0.15_dp,0.62_dp)
        complex(dp), parameter :: B0 = (1.62_dp,-0.028_dp)
        real(dp),    parameter :: alpha = 0.27_dp
 
        DDVopt = -2.0_dp*pi*(1.0_dp+mu/nmass)*dens*hbarc**2/mu&
                &*(b0exp+B0*(dens/rhoc)**(alpha))

        return
    end function DDVopt
    
    !----------------------------------------------------------------------
    !  Theoretical Optical Potential Based on Chiral Theory 
    !----------------------------------------------------------------------
    complex(dp) function ChVopt(x, dens)

        ! tell the compuler to use the external subroutine for interporalate
        external :: SPLS3

        ! external variables
        real(dp), intent(in) :: x, dens 

        ! internal variables
        integer :: i, vdim

        ! reverse file parameter to interpolate
        real(dp), allocatable, dimension(:) :: rrho_inv, VReal_inv, VImag_inv
 
        ! the interpolation subroutine requires one-dim arrays for input/output
        real(dp) :: inter_dens, inter_VReal, inter_VImag
        
        ! and other arrays with same dimension of file
        real(dp), allocatable, dimension(:) :: Q, AU

        ! get the number of lines in the file -> dimension of the variables
        vdim = size(r)

        ! allocate the other parameter for the interpolation
        allocate(rrho_inv(vdim))
        allocate(VReal_inv(vdim))
        allocate(VImag_inv(vdim))
        allocate(Q(vdim))
        allocate(AU(vdim))

        ! Since the file is given for decreasing values of the density
        ! and the interpolation has to be done for increasing values
        ! we re-store the variables of the file for incresing values of
        ! the density

        rrho_inv (1:vdim) = (/ ( rrho (vdim-i+1)*mu/hmass, i=1, vdim ) /)
        VReal_inv(1:vdim) = (/ ( VReal(vdim-i+1)*mu/hmass, i=1, vdim ) /)
        VImag_inv(1:vdim) = (/ ( VImag(vdim-i+1)*mu/hmass, i=1, vdim ) /)

        ! the potential written in the file is not well defined at the origin
        ! we impose to be equal to the potential at the maximum density value
        if (dens >= rrho_inv(vdim)*rhoc) then
            ChVopt = cmplx(VReal_inv(vdim), VImag_inv(vdim))
        else if (dens < rrho_inv(vdim)*rhoc .and. dens >= rrho_inv(1)*rhoc) then
            inter_dens = dens/rhoc

            ! call interpolation subroutine for the real and imag parts
            call SPLS3 (rrho_inv,VReal_inv,vdim,inter_dens,inter_VReal,1,Q,AU,1,0)
            call SPLS3 (rrho_inv,VImag_inv,vdim,inter_dens,inter_VImag,1,Q,AU,1,0)

            ChVopt = cmplx(inter_VReal, inter_VImag)
        else
            ChVopt = cmplx(0.0_dp, 0.0_dp)
        end if
        
        ! allocate the other parameter for the interpolation
        deallocate(rrho_inv)
        deallocate(VReal_inv)
        deallocate(VImag_inv)
        deallocate(Q)
        deallocate(AU)

        return
    end function ChVopt
    
    !----------------------------------------------------------------------
    !  Scattering amplitude based Self-Energy (SE)
    !----------------------------------------------------------------------
    complex(dp) function SelfEnergy(x, dens, Bh, part_charge, pot_type)

        ! external variables
        complex(dp), intent(in) :: Bh
        real(dp),    intent(in) :: x, dens
        integer,     intent(in) :: part_charge, pot_type

        ! internal variables
        complex(dp) :: wh, tp, tn, ThN, FhN, Vh
        real(dp), parameter :: Bn = 8.5_dp
        real(dp) :: xih, xin, Eth
        real(dp) :: sqrts_a, sqrts_b, sqrts_out

        ! interpolation variables
        real(dp) :: Itpr, Itpi, Itnr, Itni
        real(dp), allocatable, dimension(:) :: Q, AU
        integer :: sdim

        ! get the dimension of the file containing the scattering amplitude
        sdim = size(fsqrts)

        ! allocate interpolation variables
        allocate(Q(sdim))
        allocate(AU(sdim))

        ! declaration of needed physical variables (constats) 
        xih = hmass/(nmass + hmass)
        xin = nmass/(nmass + hmass)
        Eth = hmass + nmass
        if (part_charge == 1) then
            wh  = hmass - Bh - Vc(x)
        else if (part_charge == 2) then
            wh  = hmass - Bh
        else
            PRINT *, color("ERROR at SelfEnergy: particle charge note defined", c_red)
            STOP
        end if

        sqrts_a = fsqrts(1)
        sqrts_b = fsqrts(sdim)
        
        !call regula_falsi (sqrts_a, sqrts_b, sqrts_out)
        !call bisection (sqrts_a, sqrts_b, sqrts_out)
        sqrts_out = random_guess()

        if (pot_type == 4) then
            SelfEnergy = ThN*dens*hbarc**3
        else if (pot_type == 5) then
            SelfEnergy = -4._dp*pi*sqrts_out*FhN*dens*hbarc**2/nmass
        else
            PRINT *, color("ERROR at SelfEnergy: potential type not defined", c_red)
            STOP
        end if

        ! deallocate variables
        deallocate(Q)
        deallocate(AU)

        !PRINT *, "Calculated sqrts value =", sqrts_out

        return
    contains

        real(dp) function random_guess() result (sqrts1)

            ! external variables

            ! internal variables
            real(dp) :: sqrts2
            real(dp), parameter :: epsi = 1.E-2_dp
            real(dp), parameter :: Dsqrts = 1.E-2_dp
            integer :: counter

            sqrts1 = 1.0E2_dp
            sqrts2 = 0._dp

            do while (abs(sqrts1 - sqrts2) > 0._dp)
                
                counter = counter + 1
                sqrts1 = sqrts1 + Dsqrts

                if (sqrts1 > fsqrts(sdim)) then
                    Itpr = tpr(sdim)
                    Itpi = tpi(sdim)
                    Itnr = tnr(sdim)
                    Itni = tni(sdim)
                else if (sqrts1 < fsqrts(1)) then
                    Itpr = tpr(1)
                    Itpi = tpi(1)
                    Itnr = tnr(1)
                    Itni = tni(1)
                else
                    ! interpolation of the scattering amplitude for a given sqrt(s) (cm energy)
                    call SPLS3 (fsqrts,tpr,sdim,sqrts1,Itpr,1,Q,AU,1,0)
                    call SPLS3 (fsqrts,tpi,sdim,sqrts1,Itpi,1,Q,AU,1,0)
                    call SPLS3 (fsqrts,tnr,sdim,sqrts1,Itnr,1,Q,AU,1,0)
                    call SPLS3 (fsqrts,tni,sdim,sqrts1,Itni,1,Q,AU,1,0)
                end if

                if (part_charge == 1) then
                    tp = cmplx(Itpr,Itpi)               ! proton  scattering amplitude
                    tn = cmplx(Itnr,Itni)               ! neutron scattering amplitude
                    ThN = 0.5_dp*(tp + tn)              ! total scattering amplitude
                else if (part_charge == 2) then
                    tn = cmplx(Itpr,Itpi)               ! proton  scattering amplitude
                    tp = cmplx(Itnr,Itni)               ! neutron scattering amplitude
                    ThN = 0.5_dp*(tp + tn)              ! total scattering amplitude
                end if

                Vh = 0.5_dp*(1.0_dp + nmass/wh)*ThN*dens*hbarc**3/sqrts1
            
                if (part_charge == 1) then
                    sqrts2 = Eth - Bn - xin*(real(Bh)+Vc(x)) - 15.1_dp*(dens/rhoc)**(2./3) + xih*real(Vh)
                else if (part_charge == 2) then
                    sqrts2 = Eth - Bn - xin*real(Bh) - 15.1_dp*(dens/rhoc)**(2./3) + xih*real(Vh)
                end if
                
                PRINT *, "Loop", counter, sqrts1, sqrts2, abs(sqrts1 - sqrts2)
                if (sqrts2 < 0._dp) STOP
            end do

            PRINT *, "OUT OF THE LOOP!"
            STOP

        end function random_guess

        subroutine regula_falsi (xa_in, xb_in, xt)

            ! external variables
            real(dp), intent(in)  :: xa_in, xb_in
            real(dp), intent(out) :: xt

            ! internal variables
            real(dp) :: ffunc, f_xa, f_xb
            real(dp) :: xa, xb
            real(dp), parameter :: epsi = 1.0E-6_dp

            ! define variables
            xa = xa_in
            xb = xb_in

            ! method condition: f(xa)*f(xd) < 0
            !if (func(xa) * func(xb) > 0) then
            !    print *, "Method condition: f(xa)*f(xd) < 0 not fulfiled!"
            !    STOP
            !end if

            ! 'Regula falsi' method to get sqrts
            f_xa = func(xa)
            f_xb = func(xb)
            xt = xa - f_xa * (xb - xa) / (f_xb - f_xa)
            do while (abs(xa - xb) > epsi*xt)
                ffunc = f_xa * func(xt)
                if (ffunc < 0) then
                    xb = xt
                else if (ffunc > 0) then
                    xa = xt
                else
                    print *, color("ERROR at 'Regula Falsi' subroutine: function not greater nor smaller than 0", c_red)
                end if
                !print *, func(xt), xt
                !if (func(xt) == 0._dp) STOP
                !if (xt == xa .or. xt == xb) then
                !    PRINT *, "'Regula falsi' didn't coverged! Function value =", ffunc, "when it should be 0."
                !    go to 738
                !end if
                f_xa = func(xa)
                f_xb = func(xb)
                xt = xa - f_xa * (xb - xa) / (f_xb - f_xa)
            end do
738         continue

            return
        end subroutine regula_falsi

        subroutine bisection (xa_in, xb_in, xt)

            ! external variables
            real(dp), intent(in)  :: xa_in, xb_in
            real(dp), intent(out) :: xt

            ! internal variables
            real(dp) :: ffunc
            real(dp) :: xa, xb
            real(dp), parameter :: epsi = 1.0E-6_dp

            ! define variables
            xa = xa_in
            xb = xb_in

            ! method condition: f(xa)*f(xd) < 0
            !if (func(xa) * func(xb) > 0) then
            !    print *, "Method condition: f(xa)*f(xd) < 0 not fulfiled!"
            !    STOP
            !end if

            ! Bisection method to get sqrts
            xt = 0.5_dp*(xa + xb)
 !           PRINT *, "***", xa, func(xa), xb, func(xb)
            do while (abs(xa - xb) > epsi*xt)
                ffunc = func(xa) * func(xt)
                if (ffunc < 0) then
                    xb = xt
                else if (ffunc > 0) then
                    xa = xt
                else
                    print *, color("ERROR at Bisection subroutine:function not greater nor smaller than 0", c_red)
                    STOP
                end if
                !if (xt == xa .or. xt == xb) then
                !    PRINT *, "Bisection didn't coverged! Function value =", ffunc, "when it should be 0."
                !    go to 739
                !end if
                xt = 0.5_dp*(xa + xb)
!                PRINT *, "***", xa, func(xa), xb, func(xb)
            end do
739         continue

            return
        end subroutine bisection

        real(dp) function func(sqrts1)

            ! external subroutine
            external :: SPLS3

            ! external variables
            real(dp), intent(in) :: sqrts1

            ! internal variables
            real(dp) :: sqrts2
            
            ! interpolation of the scattering amplitude for a given sqrt(s) (cm energy)
            call SPLS3 (fsqrts,tpr,sdim,sqrts1,Itpr,1,Q,AU,1,0)
            call SPLS3 (fsqrts,tpi,sdim,sqrts1,Itpi,1,Q,AU,1,0)
            call SPLS3 (fsqrts,tnr,sdim,sqrts1,Itnr,1,Q,AU,1,0)
            call SPLS3 (fsqrts,tni,sdim,sqrts1,Itni,1,Q,AU,1,0)

            if (part_charge == 1) then
                tp = cmplx(Itpr,Itpi)               ! proton  scattering amplitude
                tn = cmplx(Itnr,Itni)               ! neutron scattering amplitude
                ThN = 0.5_dp*(tp + tn)              ! total scattering amplitude
            else if (part_charge == 2) then
                tn = cmplx(Itpr,Itpi)               ! proton  scattering amplitude
                tp = cmplx(Itnr,Itni)               ! neutron scattering amplitude
                ThN = 0.5_dp*(tp + tn)              ! total scattering amplitude
            end if

            Vh = 0.5_dp*(1.0_dp + nmass/wh)*ThN*dens*hbarc**3/sqrts1
            
            if (part_charge == 1) then
                sqrts2 = Eth - Bn - xin*(real(Bh)+Vc(x)) - 15.1_dp*(dens/rhoc)**(2./3) + xih*real(Vh)
            else if (part_charge == 2) then
                sqrts2 = Eth - Bn - xin*real(Bh) - 15.1_dp*(dens/rhoc)**(2./3) + xih*real(Vh)
            end if
            func = sqrts1 - sqrts2
            !PRINT *, "SQRTS", sqrts1, sqrts2
            return
        end function func

    end function SelfEnergy

    !----------------------------------------------------------------------
    ! Nuclear density
    !----------------------------------------------------------------------
    real(dp) function Density(x, dx, rho0)

        ! external variables
        real(dp), intent(in) :: x, dx, rho0
        
        if (A <= 16) then
            ! Density for light nuclei
            Density = rho0*(1.0_dp+a0*(x/Rn)**2)*exp(-x*x/(Rn*Rn))
        else if (A > 16) then
            ! Density for heavy nuclei
            Density = rho0/(1.0_dp+exp((x-Rn)/a0))
        end if

    end function Density

    !----------------------------------------------------------------------
    ! Function which integrates Woods-Saxon distribution
    !----------------------------------------------------------------------
    real(dp) function DenSuma(dx)
    
        ! external variables
        real(dp), intent(in) :: dx

        ! internal variables
        real(dp) :: x, suma
        integer :: N

        suma = 0.0_dp
        N = 10              ! times the nuclear radius
        x = 0.0_dp

        if (A <= 16) then
            ! Light nuclei MHO density normalization (analytic)
            suma = (Rn*sqrt(pi))**3*(1.0_dp+1.5_dp*a0)
        else if (A > 16) then
            ! Heavy nuclei Wood-Saxon density normalization (numerical integration)
            do while (x < N*Rn)
                suma = suma + x*x*dx / (1.0_dp + exp((x-Rn)/a0))
                x = x + dx
            end do
            suma = 4.0_dp*pi*suma
        end if

        DenSuma = suma
    end function DenSuma

end module potentials
