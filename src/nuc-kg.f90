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

    use modules
!
    implicit none
!
    integer, parameter :: scat_fdim = 201
!
    complex*16 :: BE, Vc, Vopt, Eigen, SelfEnergy
   	real*8  :: den0, den, x, dx, xturn, const, dxl, dxr
   	real*8  :: t1, t2, mu, rho0, beta
   	real*8, dimension(50) :: r, rrho, pf, VReal, VImag
    real*8, dimension(scat_fdim) :: fsqrts, tpr, tpi, tnr, tni
    integer :: i, inc, pot_type, inp_type, le, ri
!
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
    common /intevar/ r, rrho, pf, VReal, VImag
    common /scatvar/ fsqrts, tpr, tpi, tnr, tni
!
! 	Read the Chiral potential from the file
!
    open(1, file='vpot.pika', status='old', action='read')
    read(1,*)
    do i = 1, 50
       read(1,*) r(i), rrho(i), pf(i), VReal(i), VImag(i)
    end do
    close(1)
!
!   Read scattering amplitudes from the file
!
    open(4, file='amplKN10.dat', status='old', action='read')
    read(4,*)
    do i = 1, scat_fdim
       read(4,*) fsqrts(i), tpr(i), tpi(i), tnr(i), tni(i)
    end do
    close(4)

!
    print *, "----------------------------------------------------------"
    print *, "            * Klein-Gordon Equation Solver *              "
    print *, "----------------------------------------------------------"
!
73	continue
    print *, "Which potential would you like to use? (Ctrl+C to stop)"
    print *, "1. Phenomenological"
    print *, "2. Phenomenological density dependent"
    print *, "3. Chiral theory based"
    print *, "4. Real square well"
    print *, "5. Scattering amplitude based"
    read(*,*,err=73) pot_type 
    print *,  pot_type
    
    if (pot_type /= 1.and.pot_type /= 2&
        &.and.pot_type /= 3.and.pot_type /= 4.and.pot_type /= 5) then
        print *, "You shall choose between 1, 2, 3, 4 or 5!"
        call sleep(2)
        go to 73
    end if
!	
!		print *, "Enter forward integration step #: 1E-#"
!   	read(*,*) le
!   	print *, "Enter backward integration step #: 1E-#"
!   	read(*,*) ri
!   	dxl = 10.D0**(-le)
!   	dxr = 10.D0**(-ri)
!   	print *, dxl, dxr

    call read_file(Z, A, Rn, a0, pqn, aqn)

!    open(3,file="input.inp")
!    open(2,file="output3.txt",action="write")
!    print '(6X,A1,3X,A2,5X,A11,5X,A11)', "Z", "A", "Ec-En (keV)", "Width (keV)"
!    do
!        read(3,*,end=9) Z, A, Rn, a0, pqn, aqn
        call cpu_time(t1)
        call Numerov()
        call cpu_time(t2)
!    end do
!9  	continue
!    close(3)
!    close(2)
!
    print '(A20,F20.17,A2)', "CPU computing time = ", t2 - t1, " s"
end program nuc_klein_gordon
! **********************************************************************
! Subroutine Numerov's Method
!
subroutine Numerov()
    
    use modules

    implicit none
    complex*16 :: BE, Vc, Vopt, DDVopt, ChVopt, SelfEnergy
   	real*8  :: den0, den, x, dx, xturn, const, ECoul, Er0, Ei0,MCE, dxl, dxr
   	real*8  :: t1, t2,  mu, rho0, Er, Ei, CEx,CEy, beta
   	real*8  :: Eigen, RombergSuma, Density, CEnergy, CE0, CE1,delta, Vsquare
    integer :: inc, counter, mm, pot_type
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta  

    mu = hmass*nmass*A/(hmass+nmass*A)        ! Reduced mass
    
    dx = 1.D-1
    xturn = 1.D0*Rn
! 	What interest us from xturn is its the point within our mesh
! 	where the functions are to be evaluated
    xturn = dx*int(xturn/dx)+dx
!
    rho0 = A/DenSuma(dx)
    PRINT *, "Firs I am here = ", rho0, DenSuma(dx)
!
    go to 7
!  	First compute the eigenvalue for only Coulomb potential
    inc = 0
!
    beta = -500.D0
    delta = 1.D-3
!
!  	Initial conditions
    !Er0 = real(Vopt(0.D0))/2.D0
    Er0 = 540.D0
    CE1 = Eigen(Er0)
    CE0 = CE1-1.D0
!
!  	Maximization CEenergy() for the derivative matching
    counter = 0
    do while (CE1 > CE0)
        CE0 = CE1
        CEx = (Eigen(Er0+delta)-CE0)/delta
        MCE = dsqrt(CEx*CEx)
        Er  = Er0+delta*CEx/MCE
        CE1 = Eigen(Er)
        Er0 = Er
    !	counter = counter+1
        !print *, counter, Er, CE1
    end do
!
    print *, "ECoul = ", Er, xturn
!
!  	Now compute the energy eigenvalue with the optical potential (complex)
7   inc   = 1
    beta  = 1.D0
    delta = 1.D-3
!
!  	Initial conditions for the steepest descent maximization method

    if (pot_type == 1) then
        Er0 = abs(real(Vopt(0.D0)))
        Ei0 = abs(aimag(Vopt(0.D0)))
    else if (pot_type == 2) then
        !Er0 = 0.D0
        Er0 = abs(real(DDVopt(0.D0)))/2.D0
        Ei0 = abs(aimag(DDVopt(0.D0)))
    else if (pot_type == 3) then
        Er0 = abs(real(ChVopt(0.D0,Density(0.D0))))
        Ei0 = abs(aimag(ChVopt(0.D0,Density(0.D0))))
    else if (pot_type == 5) then
        Er0 = 100.D0
        Ei0 = 100.D0
    end if
!
    CE1 = CEnergy(Er0,Ei0)
    CE0 = CE1-1.D0
!  	Maximization CEenergy() for the derivative matching
    do while (CE1 > CE0)
        CE0 = CE1
        CEx = (CEnergy(Er0+delta,Ei0)-CE0)/delta
        CEy = (CEnergy(Er0,Ei0+delta)-CE0)/delta
        MCE = dsqrt(CEx*CEx+CEy*CEy)
        Er  = Er0+delta*CEx/MCE
        Ei  = Ei0+delta*CEy/MCE
        CE1 = CEnergy(Er,Ei)
        Er0 = Er
        Ei0 = Ei
        !print *, counter, Er*1.D3, ECoul*1.D3, (Er-ECoul)*1.D3, CE1-CE0
    !	counter = counter+1
    end do
!   
    write(*,*) Z, A, Er, 2.D0*Ei
    !write(*,'(4X,I3,2X,I3,5X,F7.3,8X,F7.3)') Z, A, (Er-ECoul)*1.D3, 2.D0*Ei*1.D3
    return
end subroutine Numerov
! ······································································
! Real Binding energy - Eigenvalue
!
real*8 function Eigen(Ereal)
    implicit none
    complex*16 :: BE, f, yl, yr, y0, y1, ynew, yrnext
	real*8  :: Rn, a0, dx, xturn, x, Emax, Emin, rho0, Density, dxl, dxr, Ereal
   	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, Dderi, tol, aux, aux2
   	real*8  :: ymatch, yright, const, scaling, xmatch, ylprev, beta, Vsquare
   	real*8, dimension(int(xturn/dx)+2) :: lwf
   	real*8, dimension(int(xturn/dx)+2) :: lx
   	real*8, dimension(int(3*xturn/dx)+1) :: rwf
   	real*8, dimension(int(3*xturn/dx)+1) :: rx
    !real*8, dimension(int(xturn/dx)-1+int(3*xturn/dx)) :: wf
    integer :: Z, A, pqn, aqn, nodes, inodes, counter, inc, pot_type, i, k, j
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!  	
    BE = Ereal
    const = 5.29177298D4*5.109989027D-1
!
    open(99,file="wf.dat")
!
    dxl = dx
    i = 1
    x = 0.D0
    y1 = 0.D0
    lwf(i) = y1
    lx(i) = x
    write(99,*) i, lx(i), lwf(i), xturn
!
    i = i + 1
    x = x + dxl
    y0 = dexp(-(x/Rn)**2)
    !print *, y1, y0
    !y0 = dxl**(aqn+1)
    lwf(i) = y0
    lx(i) = x
    write(99,*) i, lx(i), lwf(i), xturn
!
    i = i + 1
    do while (x <= xturn)
        x = x + dxl
        ynew = yl(x,y0,y1,dxl)
        y1 = y0
        y0 = ynew
        i = i + 1
        lwf(i) = real(ynew)**2+aimag(ynew)**2
        lx(i) = x
        !print *,  i, lx(i), lwf(i), xturn
        write(99,*) i, lx(i), lwf(i), xturn
    end do
    xmatch = x
    ymatch = real(ynew)
    ylprev = real(y1)
!
!	Integration from 'infinity' to turning point (no nodes)
    x = 4.D0*xturn                  ! our infinity
    dxr = dx
!
    i = int(3*xturn/dx)+1
    y1 = 1.D-3
    rwf(i) = y1
    rx(i) = x
    !print *, i, x, rwf(i), xturn
!
    i = i - 1
    x = x - dxr
!	Coulomb potential at infinity
    y0 = dexp(-(x/Rn)**2)
    rwf(i) = y0
    rx(i) = x
    !print *, i, x, rwf(i), xturn
!
    do while (int(x*10) > int(xturn*10))
        x = x - dxr
        ynew = yl(x,y0,y1,-dxr)
        y1 = y0
        y0 = ynew
        i = i - 1
        rwf(i) = real(ynew)**2+aimag(ynew)**2
        rx(i) = x
        !print *, i, rx(i), rwf(i), xturn
    end do
!	Matching of the wf and rescaling
    scaling = ymatch/real(ynew)
    yright = real(y1)*scaling
    
    do i=1, int(3*xturn/dx)+1
        if (i == int(3*xturn/dx)+1.or.i == int(3*xturn/dx)) then
            write(99,*) i, rx(i), rwf(i), xturn
        else
            write(99,*) i, rx(i), rwf(i)*scaling**2, xturn
        end if
    end do
    close(99)
    
    Eigen = dxr/abs((ylprev+yright-2.D0*ymatch))
    !print *, BE, Eigen
    return
end function Eigen
! ······································································
! Complex energy function
!
real*8 function CEnergy(Ereal,Eimag)
    implicit none
    complex*16 :: BE, f, yl, yr, y0, y1, ynew, dur, dul
    complex*16 :: ymatch, yright, ylprev, yrnext, scaling
   	real*8  :: Rn, a0, dx, xturn, x, Ereal, Eimag, const, y0real, y1real, dxl, dxr, k
   	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0, xmatch, Density, beta, eta
   	real*8  :: Reyright, Imyright, Rescal, Imscal
    integer, parameter :: Nend = 5
   	real*8, dimension(int(xturn/dx)+1) :: lwf
    complex*16, dimension(int(xturn/dx)+1) :: clwf
   	real*8, dimension(int(xturn/dx)+1) :: lx
   	real*8, dimension(int(xturn/dx)+1) :: Relwf
   	real*8, dimension(int(xturn/dx)+1) :: Imlwf
   	real*8, dimension(int((Nend)*xturn/dx)+1) :: Rerwf
   	real*8, dimension(int((Nend)*xturn/dx)+1) :: Imrwf
   	real*8, dimension(int((Nend)*xturn/dx)+1) :: rwf
   	real*8, dimension(int((Nend)*xturn/dx)+1) :: rx
    complex*16, dimension(int(Nend*xturn/dx)+1) :: crwf
    integer :: Z, A, pqn, aqn, inc, counter, pot_type, i, imatch
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!
    counter = 0
    open(99,file="wf2.dat")
    open(98,file="lwf.dat")
    open(97,file="rwf.dat")
!
    BE = dcmplx(Ereal,Eimag)
! 	Integrate from 0 to turning point
   
    PRINT *, "I am here = ", BE
 
    i = 1
    x = 0.D0
    y1 = 0.D0
    lwf(i) = abs(y1)*abs(y1)
    Relwf(i) = real(y1)
    Imlwf(i) = aimag(y1)
    lx(i) = x
    write(99,*) i, lx(i), lwf(i), xturn
    write(98,*) i, lx(i), Relwf(i), Imlwf(i), xturn
    !print *, i, lx(i), lwf(i), xturn
!
    i = i + 1
    x = x + dx
    !y0 = dexp(-(x/Rn)**2)
    y0 = dx**(aqn+1)
    lwf(i) = abs(y0)*abs(y0)
    Relwf(i) = real(y0)
    Imlwf(i) = aimag(y0)
    lx(i) = x
    write(99,*) i, lx(i), lwf(i), xturn
    write(98,*) i, lx(i), Relwf(i), Imlwf(i), xturn
    !print *, i, lx(i), lwf(i), xturn
!
    imatch = int(xturn/dx)+1
!
    do while (i < imatch)
        x = x + dx
        ynew = yl(x,y0,y1,dx)
        y1 = y0
        y0 = ynew
        i = i + 1
        lwf(i) = abs(ynew)**2
        clwf(i) = ynew
        Relwf(i) = real(ynew)
        Imlwf(i) = aimag(ynew)
        lx(i) = x
        !print *,  i, lx(i), lwf(i), xturn
        write(99,*) i, lx(i), lwf(i), xturn
        write(98,*) i, lx(i), Relwf(i), Imlwf(i), xturn
    end do
!
    ymatch = ynew
    ylprev = y1
!
! 	Integrate from 'infinity' to turning point
    x = Nend*xturn
!
    i = int(Nend*xturn/dx)+1
    y1 = dexp(-(x/Rn)**2)
    rwf(i) = abs(y1)**2
    Rerwf(i) = real(y1)
    Imrwf(i) = aimag(y1)
    rx(i) = x
    !print *, i, rx(i), rwf(i), xturn
    
    i = i - 1
    x = x - dx
!	Coulomb potential at infinity
    y0 = dexp(-(x/Rn)**2)
    rwf(i) = abs(y0)**2
    Rerwf(i) = real(y1)
    Imrwf(i) = aimag(y1)
    rx(i) = x
    !print *, i, rx(i), rwf(i), xturn
    
    do while (i > imatch)
        x = x - dxr
        ynew = yl(x,y0,y1,-dx)
        y1 = y0
        y0 = ynew
        i = i - 1
        !if (i == 0) then 
        !	go to 337
        !end if
        rwf(i) = abs(ynew)**2
        crwf(i) = ynew
        rx(i) = x
        !print *, i, rx(i), rwf(i), xturn
        !print *, x, xturn
    end do
337	continue
!
! 	Define logarithmic dereivatives of the wf for the matching
!
    dur = (y1 - ynew)/dx
    dul = (ymatch - ylprev)/dx
!
    scaling = ymatch/ynew
    yright = y1*scaling
!
    !Rescal = real(ymatch)/real(ynew)
    !Imscal = aimag(ymatch)/aimag(ynew)
    !Reyright = real(y1)*Rescal
    !Imyright = aimag(y1)*Imscal
!
    do i=imatch, int((Nend)*xturn/dx)+1
        crwf(i) = crwf(i)*scaling
        Rerwf(i) = real(crwf(i))
        Imrwf(i) = aimag(crwf(i))
        write(99,*) i, rx(i), rwf(i)*abs(scaling)**2, xturn
        write(97,*) i, rx(i), Rerwf(i), Imrwf(i), xturn
    end do
!
    close(99)
    close(98)
    close(97)
!
! 	CEnergy is defined as the inverse of the derivatives around the 
! 	turning point
    CEnergy = dx/abs(dur/ynew-dul/ymatch)
!
    print *, BE, CEnergy
    return
end function CEnergy
! ······································································
! Wave Function integration from left
!
complex*16 function yl(x,y0,y1,dxx)
    implicit none
    complex*16 :: BE, f, y0, y1
	real*8  :: Rn, a0, dx, dxx, xturn, x, dxl, dxr, beta
	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0
    integer :: Z, A, pqn, aqn, inc, pot_type
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
    yl = ((12.D0-10.D0*f(x))*y0-f(x-dxx)*y1)/f(x+dxx)
    return
end function yl
! ······································································
! Function containing the energies
!
complex*16 function f(x)
    implicit none
    complex*16 :: BE, Vc, Vopt, g, ChVopt, VVopt, DDVopt, SelfEnergy, SE
   	real*8  :: Rn, a0, dx, xturn, x, Density, dxl, dxr, beta
   	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0, Vsquare
    integer :: Z, A, pqn, aqn, inc, pot_type
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!
    if (pot_type == 1) then
        VVopt = Vopt(x)
    else if (pot_type == 2) then
        VVopt = DDVopt(x)
    else if (pot_type == 3) then
        VVopt = ChVopt(x,Density(x))
    else if (pot_type == 4) then
        VVopt = Vsquare(x)
    else if (pot_type == 5) then
        SE = SelfEnergy(x, Density(x), real(BE))
    end if
!
    if (pot_type == 1.or.pot_type == 2.or.pot_type == 3.or.pot_type == 4) then
       g = (mu-BE-Vc(x))**2-mu**2-2.D0*mu*VVopt-aqn*(aqn+1.D0)/(x/hbarc)**2
       !g = (mu-BE)**2-mu**2-2.D0*mu*VVopt-aqn*(aqn+1.D0)/(x/hbarc)**2
    else if (pot_type == 5) then
       g = (mu-BE-Vc(x))**2-mu**2-SE-aqn*(aqn+1.D0)/(x/hbarc)**2
       !g = (mu-BE)**2-mu**2+SE-aqn*(aqn+1.D0)/(x/hbarc)**2
    end if
    f = 1.D0+g*(dx/hbarc)**2/12.D0
    return
end function f
! ······································································
!
! Coulomb potential
!
complex*16 function Vc(x)
    implicit none
    complex*16 :: BE
	real*8  :: Rn, a0, dx, xturn, x, beta
	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0
    integer :: Z, A, pqn, aqn, inc, pot_type
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
!    
!	finesc = 1.D0/137.0359991 ! Fine structure constant
!	  
!   Spherical distribution of the charge in the nucleus
!   Different potentials for the inner and outer contributions
!
    if (x.lt.Rn) then
        Vc = -0.5D0*Z*hbarc*(3.D0-(x/Rn)**2)/(Rn*137.0359991D0)
    else
        Vc = -Z*hbarc/(x*137.0359991D0)
    end if
!	  
    return
end function Vc
!-----------------------------------------------------------------------
!
! Scattering amplitude based potential
!
complex*16 function SelfEnergy(x, rho, Bh)
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)
!    
    ! external variables
    integer, parameter :: ndim = 201
    real(kind=dp), intent(in) :: x, rho, Bh
    complex(kind=dp), external :: Vc
!
    ! commons
    real(kind=dp), dimension(ndim) :: fsqrts, tpr, tpi, tnr, tni
    real(kind=dp) :: pi, hbarc, hmass, nmass, mu, rho0
    real(kind=dp) :: BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!
    ! internal variables
    complex(kind=dp) :: Vh, wh, thN, tp, tn
    real(kind=dp), parameter :: rhoc = 0.16D0
    real(kind=dp) :: epsi, aux
    real(kind=dp) :: Ecm, Ecm1, Ecm2, Eth, xih, xin, Bn
    real(kind=dp) :: Itpr, Itpi, Itnr, Itni
    real(kind=dp), dimension(ndim) :: Q, AU
    integer :: i
!
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /scatvar/ fsqrts, tpr, tpi, tnr, tni
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!
    ! declaration of needed physical variables (constats) 
    xih = hmass/(nmass + hmass)
    xin = nmass/(nmass + hmass)
    Eth = hmass + nmass
    Bn  = 8.5D0
!
    ! bisection method to get sqrts
    Ecm1 = fsqrts(1)
    Ecm2 = fsqrts(201)
    Ecm  = 0.5D0*(Ecm1 + Ecm2)
    epsi = 1.D-8
    do while (abs(func(Ecm)) > epsi)
       if (func(Ecm1) * func(Ecm) < 0) then
          Ecm2 = Ecm
       else if (func(Ecm1) * func(Ecm) > 0) then
          Ecm1 = Ecm
       end if
       Ecm  = 0.5D0*(Ecm1 + Ecm2)
    end do
!
    SelfEnergy = thN*rho*hbarc**3
!
    ! This would be in Mâres et.al. notation, and dimensions should be checked
    !SelfEnergy = 4.D0*Ecm*thN*rho/nmass
!
    return
contains
    
    real(dp) function func(sqrts)
      implicit none
!
      external :: SPLS3
!      
      ! internal variables
      real(kind=dp) :: sqrts, aux
!
      ! interpolation of the scattering amplitude for a given sqrt(s) (cm energy)
      call SPLS3 (fsqrts,tpr,ndim,sqrts,Itpr,1,Q,AU,1,0)
      call SPLS3 (fsqrts,tpi,ndim,sqrts,Itpi,1,Q,AU,1,0)
      call SPLS3 (fsqrts,tnr,ndim,sqrts,Itnr,1,Q,AU,1,0)
      call SPLS3 (fsqrts,tni,ndim,sqrts,Itni,1,Q,AU,1,0)
 
      ! computing sqrts self-consistently
      wh = hmass - BE - Vc(x)
      tp = cmplx(Itpr,Itpi)
      tn = cmplx(Itnr,Itni)
      thN = 0.5D0*(tp + tn)
      Vh = -2.D0*pi*(1.D0 + wh/nmass)*thN*rho/wh
      aux = Eth - Bn - xin*(Bh+Vc(x)) - 15.1D0*(rho/rhoc)**(2./3) + xih*real(Vh)
      func = sqrts - aux
    end function func
 
end function SelfEnergy
! ······································································
!
! Square well potential
!
real*8 function Vsquare(x)
    implicit none
    complex*16 :: BE
	real*8  :: Rn, a0, dx, xturn, x, Density, dxl, dxr, beta
	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0, V0
    integer :: Z, A, pqn, aqn, inc, pot_type
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
    V0 = beta
    if (x <= Rn) then
        Vsquare = - V0
    else if (x > Rn) then
        Vsquare = 0.D0
    end if
    return
end function Vsquare
! ······································································
!
! Phenomenological Optical potential
!
complex*16 function Vopt(x)
    implicit none
    complex*16 :: BE
    complex*16, parameter :: b0 = (0.52D0,0.8D0)
	real*8  :: Rn, a0, dx, xturn, x, Density, dxl, dxr, beta
	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0
    integer :: inc, pot_type
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!   	
    if (inc == 0) then
        Vopt = 0.D0
    else if (inc == 1) then
        Vopt = -2.D0*pi*(1.D0+mu/nmass)*b0*Density(x)*hbarc**2/mu
    end if
    Vopt = dcmplx(real(Vopt),beta*aimag(Vopt))
    return
end function Vopt
! ······································································
!
! Phenomenological density dependent Optical potential
!
complex*16 function DDVopt(x)
    implicit none
    complex*16 :: BE
    complex*16, parameter :: b0exp = (-0.15D0,0.62D0), B0 = (1.62D0,-0.028D0)
	real*8 :: Density, x, dx, xturn, pi, hbarc, hmass, nmass, mu, rho0, dxl, dxr, beta
	real*8, parameter :: rhoc = 0.16D0, alpha = 0.270D0
    integer :: inc, pot_type
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!	
    if (inc == 0) then
        DDVopt = 0.D0
    else if (inc == 1) then
        DDVopt = -2.D0*pi*(1.D0+mu/nmass)*Density(x)*hbarc**2/mu&
                &*(b0exp+B0*(Density(x)/rhoc)**(alpha))
    end if
    DDVopt = dcmplx(real(DDVopt),beta*aimag(DDVopt))
    return
end function DDVopt
! ······································································
!
! Chiral (kaon) Optical potential
!
complex*16 function ChVopt(x,dens)
    implicit none   
   
    complex*16 :: b0, BE
	real*8  :: Rn, a0, dx, xturn, x, Density, dxl, dxr, beta
	real*8  :: pi, hbarc, finesc, hmass, nmass, mu, rho0
    integer :: Z, A, pqn, aqn, inc, pot_type, i
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
    common /ctrlvar/ BE, dx, xturn, inc, pot_type, dxl, dxr, beta
!	
    real*8 :: dens, svre, svim
	real*8, parameter     :: rhoc = 0.16D0
	real*8, dimension(1)  :: xx, den, vre, vim
	real*8, dimension(50) :: r, rrho, pf, VReal, VImag, Q, AU, &
                           & rrhoinv, VRealinv, VImaginv
    common /intevar/ r, rrho, pf, VReal, VImag
!

    do i=1,50
        rrhoinv(i)  =  rrho(50-i+1)*mu/hmass
        VRealinv(i) = VReal(50-i+1)*mu/hmass
        VImaginv(i) = VImag(50-i+1)*mu/hmass
    end do
!    
!   
    if (inc == 0) then
        ChVopt = 0.D0
    else if (inc == 1) then
        if (dens >= rrho(1)*rhoc) then
            ChVopt = dcmplx(VReal(1),VImag(1))
        elseif (dens < rrho(1)*rhoc .and. dens >= rrho(50)*rhoc) then
            xx(1) = dens
!	    
            !call SPLS3(r,rrho*rhoc,50,xx,den,1,Q,AU,1,0)
            call SPLS3(rrhoinv,VRealinv,50,xx/rhoc,vre,1,Q,AU,1,0)
            call SPLS3(rrhoinv,VImaginv,50,xx/rhoc,vim,1,Q,AU,1,0)
!    
            svre = vre(1)
            svim = vim(1)
!    
            ChVopt = dcmplx(svre,svim)
        else
            ChVopt = dcmplx(0.D0,0.D0)
        end if
    end if
    ChVopt = dcmplx(real(ChVopt),beta*aimag(ChVopt))
    return
end function ChVopt
! ······································································
!
! Nuclear density
!
double precision function Density(x)
    implicit none
    complex*16 :: BE
	real*8 :: Rn, a0, dx, xturn, B, RombergSuma, DenSuma, x
	real*8 :: pi, hbarc, hmass, nmass, mu, rho0
    integer A, Z, pqn, aqn, inc, pot_type
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
    common /physvar/ pi, hbarc, hmass, nmass, mu, rho0
!
    if (A <= 16) then
! 		Density for light nuclei
        Density = rho0*(1.D0+a0*(x/Rn)**2)*dexp(-x*x/(Rn*Rn))
    else if (A > 16) then
! 		Density for heavy nuclei
        Density = rho0/(1.D0+dexp((x-Rn)/a0))
    end if
    return
end function Density
! ······································································
!
! Romberg integration formula
!
double precision function RombergSuma()
    implicit none
    complex*16 :: BE
	real*8  :: Rn, a0, dx, xturn, FermiDist
	real*8  :: pi, hbarc, hmass, nmass, mu, rho0
	real*8  :: suma, xi, xf, ddx, dom, step, func, auxn, aux2n
    integer :: A, Z, N, nstep, i, j, pqn, aqn, inc, pot_type
!
    N = 10           ! Times the nuclear radius
    xi = 0.D0        ! Initial integration point
    xf = N*Rn        ! Final integration point
    dom = xf - xi    ! Integration domain
    nstep = 1000     ! Number of integration steps
    step = dom/nstep ! Integration step
    do j = 1, 2
        suma = 0.D0
        do i = 0, j*nstep
            if (i == 0) then
                suma = suma + 0.5D0*step*FermiDist(xi+i*step)
            else if (i == nstep) then
                suma = suma + 0.5D0*step*FermiDist(xi+i*step)
            else
                suma = suma + step*FermiDist(xi+i*step)
            end if
        end do
        if (j == 1) then
            auxn = suma
        else if (j == 2) then
            aux2n = suma
        end if
    end do
    RombergSuma = 4.D0*pi*(4.D0*aux2n-auxn)/3.D0
    return
end function RombergSuma
! ······································································
!
! Fermi distribution
!
double precision function FermiDist(x)
    implicit none
    complex*16 :: BE
    real*8  :: Rn, a0, dx, xturn, x
	real*8  :: pi, hbarc, hmass, nmass, mu, rho0
    integer :: A, Z, pqn, aqn, inc, pot_type
    common /fileinp/ Z, A, Rn, a0, pqn, aqn
!	
    FermiDist = x*x/(1.D0 + dexp((x-Rn)/a0))
    return
end function FermiDist
! ·······································································
!
! SPLINE interpolation subroutine (Courtesy of Angels Ramos, coded in F77) 
!
      !SUBROUTINE SPLS3(X,Y,N,XI,FI,M,Q,AU,IGO,ISPL)                     
!                                                                       
!     ******************************************************************
!                                                                       
!     CUBIC SPLINE, STARTING WITH ZERO SECOND DERIVATIVES AT THE        
!     BOUNDARIES OF THE APPROXIlopez1991!IMATION INTERVAL.                         
!                                                                       
!     IGO = 0      BUILD UP SPLINE ONLY.                                
!     IGO = 1      BUILD UP SPLINE AND INTERPOLATE FOR M POINTS.        
!     IGO = 2      BUILD UP SPLINE AND COMPUTE DERIVATIVES AT M POINTS. 
!                                                                       
!     REAL*8 VERSION.        J.GALONSKA, 15.12.1971           
!     ******************************************************************
      !IMPLICIT REAL*8 (A-H,O-Z)                                       
      !DIMENSION X(1),Y(1),XI(1),Q(1),AU(1),FI(1)
!                                                                       
      !ZERO=0.D0                                                         
      !THREE=3.D0                                                        
      !SIX=6.D0                                                          
      !FACT=0.1666666666667D0                                            
      !IF (ISPL.NE.0)  GO TO 30                                          
!                                                                       
!                                                                       
      !AU(1) = ZERO                                                      
      !AU(N) = ZERO                                                      
      !Q(1) = ZERO                                                       
      !HK = X(2) - X(1)                                                  
      !YSAVE = (Y(2)-Y(1)) / HK                                          
      !AUX = ZERO                                                        
      !NN = N - 1                                                        
      !DO 10  K = 2,NN                                                   
      !HX = X(K+1) - X(K-1)                                              
      !DIVQ = (HK*Q(K-1)+HX+HX)                                          
      !HK = X(K+1) - X(K)                                                
      !YK = (Y(K+1)-Y(K)) / HK                                           
      !Q(K) = - HK / DIVQ                                                
      !AU(K) = (SIX*(YK-YSAVE)-AUX) / DIVQ                               
      !YSAVE = YK                                                        
      !AUX = AU(K) * HK                                                  
   !10 CONTINUE                                                          
!                                                                       
      !NN2 = NN + 2                                                      
      !DO 20  KK = 2,NN                                                  
      !K = NN2 - KK                                                      
   !20 AU(K) = Q(K) * AU(K+1) + AU(K)                                    
!                                                                       
      !IF (IGO.EQ.0)  RETURN                                             
!                                                                       
!     ******************************************************************
!                                                                       
!     INTERPOLATION OR COMPUTATION OF DERIVATIVES.                      
!                                                                       
!     IGO = 1      INTERPOLATE FOR M POINTS.                            
!     IGO = 2      COMPUTE DERIVATIVES AT M POINTS.                     
!                                                                       
!     ******************************************************************
!                                                                       
   !30 DO 100  J = 1,M                                                   
      !IF (X(1).GT.XI(J))  THEN                                         
      !M1=1
      !M2=2
      !GO TO 90
      !ELSE
      !ENDIF
      !IF (XI(J).GT.X(N))  THEN                                       
      !M1=N-1
      !M2=N
      !GO TO 90
      !ELSE
      !ENDIF
      !M1 = 1                                                            
      !M2 = N                                                            
   !50 M3 = (M2+M1)/2                                                    
      !IF (XI(J).GE.X(M3))  GO TO 70                                     
      !M2 = M3                                                           
      !GO TO 80                                                          
   !70 M1 = M3                                                           
   !80 IF (M1+1-M2.NE.0)  GO TO 50                                       
   !90 DIJ = X(M2) - XI(J)                                               
      !DIM1J = X(M1) - XI(J)                                             
      !HI = X(M2) - X(M1)                                                
      !HI2 = HI * HI                                                     
      !IF (IGO.GE.2)  GO TO 95                                           
      !DIJ3 = DIJ * DIJ * DIJ                                            
      !DIM1J3 = DIM1J * DIM1J * DIM1J                                    
      !FI(J) = FACT * (AU(M1)*DIJ3-AU(M2)*DIM1J3+(SIX*Y(M1)-HI2*AU(M1))  &
     !&        *DIJ-(SIX*Y(M2)-HI2*AU(M2))*DIM1J) / HI                   
      !GO TO 100                                                         
   !95 FI(J) = FACT * (THREE*(AU(M2)*DIM1J*DIM1J-AU(M1)*DIJ*DIJ)&         
     !&       -SIX*(Y(M1)-Y(M2))+HI2*(AU(M1)-AU(M2))) / HI               
  !100 CONTINUE                                                          
      !RETURN                                                            
!                                                                       
! 110 M1 = 1                                                            
!     M2 = 2                                                            
!     GO TO 90                                                          
!                                                                       
! 120 M1 = N - 1                                                        
!     M2 = N                                                            
!     GO TO 90                                                          
!                                                                       
      !END   
