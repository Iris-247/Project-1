program project2
    implicit none

    ! Numerical precision
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: jmax = 5000
    integer, parameter :: nstep = 20000
    integer, parameter :: jcell = jmax - 1

    ! Physical constants
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: kpc = 3.084e21_dp
    real(dp), parameter :: year = 365.25_dp * 24.0_dp * 3600.0_dp
    real(dp), parameter :: ZFe_sun = 1.8e-3_dp

    integer :: j, n
    real(dp) :: dt, tmax, time

    ! Radial grid:
    real(dp) :: r(jmax), rHalf(jcell)
    real(dp) :: rhoGas(jcell)
    real(dp) :: velocity_turb, length_turb, D
    real(dp) :: rGasKpc


    real(dp) :: ZFe_old(jcell), ZFe_new(jcell)
    real(dp) :: F(jcell-1), rhoFace(jcell-1)

    integer :: ios

    ! Set the turbulent using diffusion coefficient D = v_t * l_t / 3
    velocity_turb = 2.0e7_dp
    length_turb = 10.0_dp * kpc
    D = velocity_turb * length_turb / 3.0_dp

    ! Set the time step and total simulated time
    dt = 0.05e6_dp * year
    tmax = dt * real(nstep, dp)

    ! Read the radial grid from grid.dat.
    open(unit=10, file='grid.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot open grid.dat'

    do j = 1, jmax - 1
        read(10, *, iostat=ios) r(j), rHalf(j)
        if (ios /= 0) error stop 'Error reading grid.dat'
    end do

    close(10)

    ! Reconstruct the outermost grid boundary using the grid spacing
    r(jmax) = r(jmax - 1) + r(2) - r(1)

    ! Convert from kpc to cm
    r = r * kpc
    rHalf = rHalf * kpc

    ! Read the background gas density from gas.dat.
    open(unit=10, file='gas.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot open gas.dat'

    do j = 1, jcell
        read(10, *, iostat=ios) rGasKpc, rhoGas(j)
        if (ios /= 0) error stop 'Error reading gas.dat'
    end do

    close(10)
    ! Set the initial profile 
    do j = 1, jcell
        ZFe_old(j) = 0.3_dp * (2.2_dp + ((rHalf(j) / kpc) / 80.0_dp)**3) / (1.0_dp + ((rHalf(j) / kpc) / 80.0_dp)**3) * ZFe_sun
    end do

    open(unit=20, file='ZFe_all_times.dat', status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) error stop 'Cannot open ZFe_all_times.dat'

    ! step number, time [Myr], radius [kpc], iron abundance same below
    do j = 1, jcell
        write(20, '(I8,1X,3(1PE16.8,1X))') 0, 0.0_dp, rHalf(j) / kpc, ZFe_old(j)
    end do

    ! Advance the iron abundance through time
    do n = 1, nstep

        ! Calculate the gas density and diffusive flux at each cell interface
        do j = 1, jmax - 2
            rhoFace(j) = 0.5_dp * (rhoGas(j) + rhoGas(j+1))

            F(j) = r(j+1)**2 * D * rhoFace(j) * (ZFe_old(j+1) - ZFe_old(j)) / (rHalf(j+1) - rHalf(j))
        end do

        ! Update the iron abundance in the interior cells using the
        do j = 2, jcell - 1
            ZFe_new(j) = ZFe_old(j) + dt * (F(j) - F(j-1)) / (1.4_dp * rhoGas(j) * (r(j+1)**3 - r(j)**3) / 3.0_dp)
        end do

        ! Apply the inner and outer boundary conditions Inner boundary approximate zero gradient at the center
        ZFe_new(1) = ZFe_new(2)
        ZFe_new(jcell) = ZFe_old(jcell)
        time = real(n, dp) * dt

        do j = 1, jcell
            write(20, '(I8,1X,3(1PE16.8,1X))') n, time / (1.0e6_dp * year), rHalf(j) / kpc, ZFe_new(j)
        end do

        ZFe_old = ZFe_new
    end do

    close(20)

end program project2