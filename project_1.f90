program project1
    implicit none

    integer, parameter :: dp = selected_real_kind(15,307) !!use this precision
    integer, parameter :: jmax = 5000   !! dimension of the grid
    
    real(dp), parameter :: pi  = 3.141592653589793_dp

    real(dp) :: r(jmax)
    real(dp) :: rHalf(jmax-1), shellVol(jmax-1)
    real(dp) :: rhoNFW(jmax-1)
    real(dp) :: massDM(jmax), theoreticalMassDM(jmax) !! variables, rhoNFW is DM density NFW
 
    real(dp) :: kpc, rmin, rmax, rs, rhoNFW0, Msun, dr, x
    integer :: j

    kpc     = 3.084e21_dp
    rmin    = 0.0_dp
    rmax    = 3000.0_dp * kpc
    rs      = 435.7_dp * kpc
    rhoNFW0 = 7.35e-26_dp  
    Msun    = 1.9885e33_dp   !! constants (cgs)

    !! step 1: building the grids (normal and half stepped)
    dr = (rmax - rmin) / real(jmax - 1, dp)

    do j = 1, jmax
        r(j) = rmin + real(j - 1, dp) * dr
    end do
    
    ! Step 2: shell centers, volumes, and NFW density
    do j = 1, jmax - 1
        rHalf(j) = 0.5_dp * (r(j) + r(j+1))

        shellVol(j) = (4.0_dp * pi / 3.0_dp) * (r(j+1)**3 - r(j)**3)

        x = rHalf(j) / rs
        rhoNFW(j) = rhoNFW0 / (x * (1.0_dp + x)**2)
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

        x = r(j+1) / rs
        theoreticalMassDM(j+1) = 4.0_dp * pi * rhoNFW0 * rs**3 * (log(1.0_dp + x) - x / (1.0_dp + x))  !! theoretical DM
    end do
    
    ! Step 5: Output Columns with radius [kpc], numerical mass [Msun], theoretical mass [Msun]
    open(unit=20, file='mass.dat', status='replace', action='write')
        do j = 1, jmax
        write(20, '(3(1pe16.8))') r(j) / kpc, &
                                  massDM(j) / Msun, &
                                  theoreticalMassDM(j) / Msun
    end do
    close(20)
end program project1