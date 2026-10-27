program project1
    implicit none
    
    integer :: j
    
    integer, parameter :: dp = selected_real_kind(15,307) !!use this precision
    integer, parameter :: jmax = 5000   !! dimension of the grid
    
    ! constant
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: kpc = 3.084e21_dp

    real(dp), parameter :: Ggrav = 6.67430e-8_dp
    real(dp), parameter :: kBoltz = 1.380649e-16_dp
    real(dp), parameter :: mp = 1.6726219e-24_dp
    
    real(dp), parameter :: mu = 0.60_dp
    real(dp), parameter :: Tgas = 8.9e7_dp
    
    real(dp), parameter :: Tmg = 8.9e7_dp
    
    ! Part 1: Grid and dark matter (NFW) profile
    real(dp) :: r(jmax)
    real(dp) :: rHalf(jmax-1), shellVol(jmax-1)
    real(dp) :: rhoNFW(jmax-1)
    real(dp) :: massDM(jmax), theoreticalMassDM(jmax) !! variables, rhoNFW is DM density NFW
    
    real(dp) :: rmin, rmax, rs, rhoNFW0, Msun, dr, xs
    ! Part 2: Isothermal gas hydrostatic equilibrium
    integer :: jVir
    
    real(dp) :: rhoGas(jmax-1), rhoGasAnalytic(jmax-1)
    real(dp) :: gasMass(jmax), rVir, fb

    real(dp) :: rhoGas0, hseFactor, xs0, f0, f

    ! Part 3: BCG contribution to the gravitational potential
    real(dp) :: massBCG(jmax), massTotal(jmax)
    real(dp) :: MBCG, aBCG
    real(dp) :: targetFb 
    real(dp) :: gasMassTarget

    ! Part 4: Non-isothermal temperature profile
    real(dp) :: temperature(jmax-1)
    real(dp) :: r500 , x500
    
    ! Part 5: Non-isothermal temperature profile with turb
    real(dp) :: rhoGasTurb(jmax-1), Pturb(jmax-1)
    real(dp) :: velocity_turb

    rmin    = 0.0_dp
    rmax    = 3000.0_dp * kpc
    rs      = 435.7_dp * kpc
    rhoNFW0 = 7.35e-26_dp  
    Msun    = 1.9885e33_dp   !  (cgs)

    rhoGas0 = 4.0e-26_dp
    targetFb = 0.16_dp
    rVir = 2800.0_dp * kpc
    
    MBCG =  1.0e12_dp * Msun
    aBCG = (12.0_dp / (1.0_dp + sqrt(2.0_dp) ) )* kpc 
    
    r500 = 1400.0_dp * kpc

    velocity_turb = 5.0e7_dp 

    ! step 1: building the grids (normal and half stepped)
    dr = (rmax - rmin) / real(jmax - 1, dp)

    do j = 1, jmax
        r(j) = rmin + real(j - 1, dp) * dr
    end do
    
    ! Step 2: shell centers, volumes, and NFW density
    do j = 1, jmax - 1
        rHalf(j) = 0.5_dp * (r(j) + r(j+1))

        shellVol(j) = (4.0_dp * pi / 3.0_dp) * (r(j+1)**3 - r(j)**3)

        xs = rHalf(j) / rs
        rhoNFW(j) = rhoNFW0 / (xs * (1.0_dp + xs)**2)
    end do

    ! Step 3 :Grid output with left boundary and shell center, in kpc
    open(unit=10, file='grid.dat', status='replace', action='write')
    do j = 1, jmax -1
        write(10,'(2(1pe16.8))') r(j) / kpc, rHalf(j) / kpc
    end do
    close(10)

    ! Step 4: enclosed DM mass at each grid boundary
    massDM(1) = 0.0_dp
    theoreticalMassDM(1) = 0.0_dp

    do j = 1, jmax - 1
        massDM(j+1) = massDM(j) + rhoNFW(j) * shellVol(j)  !! numerical DM

        xs = r(j+1) / rs
        theoreticalMassDM(j+1) = 4.0_dp * pi * rhoNFW0 * rs**3 * (log(1.0_dp + xs) - xs / (1.0_dp + xs))  !! theoretical DM
    end do

    ! Step 5: Output Columns with radius [kpc], numerical mass [Msun], theoretical mass [Msun]
    open(unit=20, file='mass.dat', status='replace', action='write')
        do j = 1, jmax
        write(20, '(3(1pe16.8))') r(j) / kpc, massDM(j) / Msun, theoreticalMassDM(j) / Msun
    end do
    close(20)

    ! Step 6 :Isothermal HSE solution with numerical and analytic profiles
    hseFactor = 4.0_dp * pi * Ggrav * rhoNFW0 * rs**2 * mu * mp / (kBoltz * Tgas)

    ! Set the first shell-center density using the central density rhoGas0
    rhoGas(1) = rhoGas0

    xs0 = rHalf(1) / rs
    f0 = 1.0_dp - log(1.0_dp + xs0) / xs0
    rhoGasAnalytic(1) = rhoGas0 * exp(-hseFactor * f0)

    ! Integrate HSE outward. Gas densities are at shell centers;
    ! enclosed DM masses are at grid boundaries.
    do j = 1, jmax - 2
        rhoGas(j+1) = rhoGas(j) * exp( -dr * Ggrav * massDM(j+1) * mu * mp / (kBoltz * Tgas * r(j+1)**2) )

        xs = rHalf(j+1) / rs
        f = 1.0_dp - log(1.0_dp + xs) / xs
        rhoGasAnalytic(j+1) = rhoGas0 * exp(-hseFactor * f)
    end do

    ! Integrate gas mass over shells
    gasMass(1) = 0.0_dp
    do j = 1, jmax - 1
        gasMass(j+1) = gasMass(j) + rhoGas(j) * shellVol(j)
    end do

    !Step 7: Baryon fraction at the nearest grid point to rVir
    jVir = nint(rVir / dr) + 1
    fb = gasMass(jVir) / (massDM(jVir) + gasMass(jVir))

    write(*,*) "r at baryon-fraction check [kpc] =", r(jVir) / kpc
    write(*,*) "Baryon fraction at rVir =", fb

    open(unit=30, file='gas.dat', status='replace', action='write')
    do j = 1, jmax - 1
        write(30, '(4(1pe24.16))') rHalf(j) / kpc, rhoGas(j), rhoGasAnalytic(j) ,rhoNFW(j)
    end do
    close(30)

    ! Step 8 :Isothermal HSE solution with BCG
    do j = 1, jmax
        massBCG(j) = MBCG * r(j)**2 / (r(j) + aBCG)**2
        massTotal(j) = massDM(j) + massBCG(j)
    end do

    do j = 1, jmax - 2
        rhoGas(j+1) = rhoGas(j) * exp( -dr * Ggrav * massTotal(j+1) * mu * mp / (kBoltz * Tgas * r(j+1)**2) )
    end do

    ! Integrate gas mass over shells
    gasMass(1) = 0.0_dp
    do j = 1, jmax - 1
        gasMass(j+1) = gasMass(j) + rhoGas(j) * shellVol(j)
    end do
    
    !Scale the gas profile to reach the target baryon fraction
    gasMassTarget = targetFb / (1.0_dp - targetFb) * massDM(jVir) - massBCG(jVir)
    
    rhoGas0 = rhoGas0 * (gasMassTarget / gasMass(jVir))
    rhoGas = rhoGas * (gasMassTarget / gasMass(jVir))
    gasMass = gasMass * (gasMassTarget / gasMass(jVir))
    
    fb = (gasMass(jVir) + massBCG(jVir)) / (massDM(jVir) + massBCG(jVir) + gasMass(jVir))
    
    write(*,*) "Baryon fraction with BCG =", fb

    open(unit=30, file='gas_BCG.dat', status='replace', action='write')
    do j = 1, jmax - 1
        write(30, '(3(1pe24.16))') rHalf(j) / kpc, rhoGas(j) , rhoNFW(j)
    end do
    close(30)

    ! step 9 Non-isothermal temperature profile 
    ! reset initial value
    rhoGas0 = 4.0e-26_dp
    rhoGas(1) = rhoGas0

    do j = 1, jmax - 1
    x500 = rHalf(j) / r500
    temperature(j) = Tmg * 1.35_dp * ((x500 / 0.045_dp)**1.9_dp + 0.45_dp) / ((x500 / 0.045_dp)**1.9_dp + 1.0_dp) / (1.0_dp + (x500 / 0.6_dp)**2)**0.45_dp
    end do
    do j = 1 , jmax - 2
        rhoGas(j+1) = rhoGas(j) * temperature(j) / temperature(j+1) * exp(-dr * Ggrav * massTotal(j+1) * mu * mp / (kBoltz * temperature(j) * r(j+1)**2))
    end do
    
    ! Integrate gas mass over shells
    gasMass(1) = 0.0_dp
    do j = 1, jmax - 1
        gasMass(j+1) = gasMass(j) + rhoGas(j) * shellVol(j)
    end do

    !Scale the gas profile to reach the target baryon fraction
    gasMassTarget = targetFb / (1.0_dp - targetFb) * massDM(jVir) - massBCG(jVir)
    
    rhoGas0 = rhoGas0 * (gasMassTarget / gasMass(jVir))
    rhoGas = rhoGas * (gasMassTarget / gasMass(jVir))
    gasMass = gasMass * (gasMassTarget / gasMass(jVir))
    
    fb = (gasMass(jVir) + massBCG(jVir)) / (massDM(jVir) + massBCG(jVir) + gasMass(jVir))
    
    write(*,*) "Baryon fraction with BCG =", fb

    open(unit=30, file='gas_BCG_non-iso.dat', status='replace', action='write')
    do j = 1, jmax - 1
        write(30, '(3(1pe24.16))') rHalf(j) / kpc, rhoGas(j) , rhoNFW(j)
    end do
    close(30)

    ! step 10 Non-isothermal density profile with turburlence
    ! reset initial value
    rhoGas0 = 4.0e-26_dp
    rhoGas(1) = rhoGas0
    rhoGasTurb(1) = rhoGas(1)
    
    do j = 1, jmax - 1
    x500 = rHalf(j) / r500
    temperature(j) = Tmg * 1.35_dp * ((x500 / 0.045_dp)**1.9_dp + 0.45_dp) / ((x500 / 0.045_dp)**1.9_dp + 1.0_dp) / (1.0_dp + (x500 / 0.6_dp)**2)**0.45_dp
    end do
    do j = 1 , jmax - 2
        rhoGas(j+1) = rhoGas(j) * temperature(j) / temperature(j+1) * exp(-dr * Ggrav * massTotal(j+1) * mu * mp / (kBoltz * temperature(j) * r(j+1)**2))
    end do
    ! add turb term
    do j = 1, jmax - 1
        Pturb(j) = 1.0_dp / 3.0_dp * rhoGas(j) * velocity_turb **2
    end do
    do j = 1 , jmax - 2
        rhoGasTurb(j+1) = mu * mp / (kBoltz * temperature(j+1)) * (rhoGasTurb(j) * kBoltz * temperature(j) / (mu * mp) - dr * rhoGasTurb(j) * Ggrav * massTotal(j+1) / r(j+1)**2 - (Pturb(j+1) - Pturb(j)))
    end do


    ! Integrate gas mass over shells
    gasMass(1) = 0.0_dp
    do j = 1, jmax - 1
        gasMass(j+1) = gasMass(j) + rhoGasTurb(j) * shellVol(j)
    end do

    ! Scale the gas profile to reach the target baryon fraction
    gasMassTarget = targetFb / (1.0_dp - targetFb) * massDM(jVir) - massBCG(jVir)
    
    rhoGas0 = rhoGas0 * (gasMassTarget / gasMass(jVir))
    rhoGas  = rhoGas * (gasMassTarget / gasMass(jVir))
    Pturb  = Pturb * (gasMassTarget / gasMass(jVir))
    rhoGasTurb = rhoGasTurb * (gasMassTarget / gasMass(jVir))
    gasMass = gasMass * (gasMassTarget / gasMass(jVir))
    
    fb = (gasMass(jVir) + massBCG(jVir)) / (massDM(jVir) + massBCG(jVir) + gasMass(jVir))
    
    write(*,*) "Baryon fraction with BCG =", fb

    open(unit=30, file='gas_BCG_non-iso-turb.dat', status='replace', action='write')
    do j = 1, jmax - 1
        write(30, '(3(1pe24.16))') rHalf(j) / kpc, rhoGasTurb(j) , rhoNFW(j)
    end do
    close(30)

end program project1