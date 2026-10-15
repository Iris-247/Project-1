program project1
    implicit none
    integer, parameter :: jmax = 5000   !! dimension of the grid
    real*8 :: r(jmax), hr(jmax), Vol(jmax)
    real*8 :: rhonfw(jmax), Mdm(jmax), theoricalMdm(jmax)  !! variables, rhonwf is DM density NFW, Mdm is dark matter mass
    real*8 :: x(jmax), y(jmax)
    real*8 :: kpc, pi, rs, rhonfw0, Msun   !! constants (cgs)
    real*8 :: rmin, rmax
    integer :: j

    kpc     = 3.084d21
    pi      = 3.141592653589793d0
    rmin    = 0.d0
    rmax    = 3000.d0 * kpc
    rs      = 435.7d0 * kpc
    rhonfw0 = 7.35d-26
    Msun    = 1.9885d33

    !! step 1: building the grids (normal and half stepped)
    do j = 1, jmax
        r(j) = rmin + (j-1)*rmax/(jmax-1)
    end do
    do j = 1, jmax - 1
        hr(j) = r(j) + 0.5d0*(r(j+1) - r(j))
    end do
    hr(jmax) = r(jmax-1) + 0.5d0*(r(jmax-2) - r(jmax -1))
    open(unit=10, file='grid.dat', status='replace', action='write')
    do j = 1, jmax
        write(10,*) real(r(j)/kpc), real(hr(j)/kpc)
    end do
    close(10)

    !! step 2.1 calculating the volume
    do j = 1, jmax-1
        Vol(j) = (4.d0*pi/3.d0) * (r(j+1)**3 - r(j)**3)  !! centering at half steps hr(j+1)
    end do
    Vol(jmax) = (4.d0*pi/3.d0) * r(jmax)**3

    !! step 2.2 calculate the NFW DM density profile (we can directly use halfgrid)
    do j = 1, jmax
        x(j) = hr(j) / rs    !! just for the formula to be clearer
        rhonfw(j) = rhonfw0 / (x(j)*(1.d0 + x(j))**2)
    end do

    !! step 2.3
    open(20, file='mass.dat')
    Mdm(1) = 0.d0
    do j = 2, jmax
        x(j) = r(j)/rs
        y(j) = r(j)/(r(j)+rs)
        Mdm(j) = Mdm(j-1) + rhonfw(j-1)*Vol(j-1)  !! numerical DM
        theoricalMdm(j) = 4.d0*pi*rhonfw0*rs**3 * (log(1.d0+x(j)) - y(j))  !! theoretical DM
        write(20,1000) r(j)/kpc, Mdm(j)/Msun, theoricalMdm(j)/Msun
    end do
    close(20)

    1000 format(3(1pe12.4))
end program project1