!
! © 2024. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

module traveltime_iso

    use libflit
    use parameters

    implicit none

    real, allocatable, dimension(:, :, :) :: t0, pdxt0, pdyt0, pdzt0, tt, lambda
    logical, allocatable, dimension(:, :, :) :: unknown, recrflag

    private
    public :: forward_iso_fast_march
    public :: forward_iso_fast_sweep
    public :: adjoint_iso

contains

    function smaller_traveltime(node1, node2) result(smaller)

        real, dimension(:), intent(in) :: node1, node2
        logical :: smaller

        smaller = node1(1) < node2(1)

    end function smaller_traveltime

    !
    !> Solve quadratic equation associated with the factorized eikonal equation
    !
    subroutine solve_quadratic(i, j, k, dx, dy, dz, vx, vy, vz, f)

        integer, intent(in) :: i, j, k
        real, intent(in) :: dx, dy, dz, vx, vy, vz, f

        double precision :: timeleft1, timeright1
        double precision :: timefront1, timeback1
        double precision :: timetop1, timebottom1
        double precision :: t0ijk
        double precision :: a, b, c
        double precision :: ax, bx, cx
        double precision :: ay, by, cy
        double precision :: az, bz, cz
        double precision :: sx, sy, sz, s
        double precision :: taux, tauy, tauz
        integer :: dirx, diry, dirz
        double precision :: tau, alpha, beta

        t0ijk = t0(i, j, k)
        timeleft1 = t0(i - 1, j, k)*tt(i - 1, j, k)
        timeright1 = t0(i + 1, j, k)*tt(i + 1, j, k)
        timefront1 = t0(i, j - 1, k)*tt(i, j - 1, k)
        timeback1 = t0(i, j + 1, k)*tt(i, j + 1, k)
        timetop1 = t0(i, j, k - 1)*tt(i, j, k - 1)
        timebottom1 = t0(i, j, k + 1)*tt(i, j, k + 1)

        taux = float_huge
        tauy = float_huge
        tauz = float_huge

        dirx = 0
        if (.not. unknown(i - 1, j, k)) then
            taux = timeleft1
            dirx = -1
        end if
        if (.not. unknown(i + 1, j, k) .and. taux > timeright1) then
            taux = timeright1
            dirx = 1
        end if

        diry = 0
        if (.not. unknown(i, j - 1, k)) then
            tauy = timefront1
            diry = -1
        end if
        if (.not. unknown(i, j + 1, k) .and. tauy > timeback1) then
            tauy = timeback1
            diry = 1
        end if

        dirz = 0
        if (.not. unknown(i, j, k - 1)) then
            tauz = timetop1
            dirz = -1
        end if
        if (.not. unknown(i, j, k + 1) .and. tauz > timebottom1) then
            tauz = timebottom1
            dirz = 1
        end if

        a = 0.0
        b = 0.0
        c = 0.0
        ax = 0.0
        bx = 0.0
        cx = 0.0
        if (dirx /= 0) then
            tau = t0(i + 2*dirx, j, k)*tt(i + 2*dirx, j, k)
            alpha = 0.0
            beta = 0.0
            if (taux > tau .and. .not. unknown(i + 2*dirx, j, k)) then
                alpha = pdxt0(i, j, k) - 1.5*t0ijk*dirx/dx
                beta = -t0ijk*dirx*(tt(i + 2*dirx, j, k) - 4*tt(i + dirx, j, k))/(2*dx)
            else
                alpha = pdxt0(i, j, k) - t0ijk*dirx/dx
                beta = t0ijk*dirx*tt(i + dirx, j, k)/dx
            end if
            ax = vx*alpha**2
            bx = vx*2*alpha*beta
            cx = vx*beta**2
            a = a + ax
            b = b + bx
            c = c + cx
        end if

        ay = 0.0
        by = 0.0
        cy = 0.0
        if (diry /= 0) then
            tau = t0(i, j + 2*diry, k)*tt(i, j + 2*diry, k)
            alpha = 0.0
            beta = 0.0
            if (tauy > tau .and. .not. unknown(i, j + 2*diry, k)) then
                alpha = pdyt0(i, j, k) - 1.5*t0ijk*diry/dy
                beta = -t0ijk*diry*(tt(i, j + 2*diry, k) - 4*tt(i, j + diry, k))/(2*dy)
            else
                alpha = pdyt0(i, j, k) - t0ijk*diry/dy
                beta = t0ijk*diry*tt(i, j + diry, k)/dy
            end if
            ay = vy*alpha**2
            by = vy*2*alpha*beta
            cy = vy*beta**2
            a = a + ay
            b = b + by
            c = c + cy
        end if

        az = 0.0
        bz = 0.0
        cz = 0.0
        if (dirz /= 0) then
            tau = t0(i, j, k + 2*dirz)*tt(i, j, k + 2*dirz)
            alpha = 0.0
            beta = 0.0
            if (tauz > tau .and. .not. unknown(i, j, k + 2*dirz)) then
                alpha = pdzt0(i, j, k) - 1.5*t0ijk*dirz/dz
                beta = -t0ijk*dirz*(tt(i, j, k + 2*dirz) - 4*tt(i, j, k + dirz))/(2*dz)
            else
                alpha = pdzt0(i, j, k) - t0ijk*dirz/dz
                beta = t0ijk*dirz*tt(i, j, k + dirz)/dz
            end if
            az = vz*alpha**2
            bz = vz*2*alpha*beta
            cz = vz*beta**2
            a = a + az
            b = b + bz
            c = c + cz
        end if

        c = c - f
        cx = cx - f
        cy = cy - f
        cz = cz - f

        sx = (-bx + sqrt(bx**2 - 4*ax*cx))/(2*ax)
        sy = (-by + sqrt(by**2 - 4*ay*cy))/(2*ay)
        sz = (-bz + sqrt(bz**2 - 4*az*cz))/(2*az)
        s = (-b + sqrt(b**2 - 4*a*c))/(2*a)
        if (isnan(sx) .or. sx <= 0) then
            sx = float_huge
        end if
        if (isnan(sy) .or. sy <= 0) then
            sy = float_huge
        end if
        if (isnan(sz) .or. sz <= 0) then
            sz = float_huge
        end if
        if (isnan(s) .or. s <= 0) then
            s = float_huge
        end if

        tt(i, j, k) = min(sx, sy, sz, s)
        unknown(i, j, k) = .false.

    end subroutine solve_quadratic

    !
    !> Fast marching eikonal solver
    !
    subroutine forward_iso_fast_march(v, d, o, geom, t, trec)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        integer :: nx, ny, nz, l, isf
        integer :: i, j, k, hi, hj, hk, ii, jj, kk
        integer :: isx, isy, isz, irx, iry, irz
        type(heap_float) :: wavefront
        real :: node(1:4), v2, dx, dy, dz, ox, oy, oz
        real, allocatable, dimension(:, :, :) :: vp
        real :: dsx, dsy, dsz, vsource

        dx = d(1)
        dy = d(2)
        dz = d(3)
        ox = o(1)
        oy = o(2)
        oz = o(3)

        vp = permute(v, 321)
        nx = size(vp, 1)
        ny = size(vp, 2)
        nz = size(vp, 3)

        ! Source location
        do i = 1, geom%ns
            if (geom%srcr(i)%amp == 1) then
                isf = i
                exit
            end if
        end do
        isx = nint((geom%srcr(isf)%x - ox)/dx) + 1
        isy = nint((geom%srcr(isf)%y - oy)/dy) + 1
        isz = nint((geom%srcr(isf)%z - oz)/dz) + 1

        ! Background time field
        vsource = vp(isx, isy, isz)
        call alloc_array(t0, [1, nx, 1, ny, 1, nz], pad=2)
        call alloc_array(pdxt0, [1, nx, 1, ny, 1, nz], pad=2)
        call alloc_array(pdyt0, [1, nx, 1, ny, 1, nz], pad=2)
        call alloc_array(pdzt0, [1, nx, 1, ny, 1, nz], pad=2)
        do k = 1 - 2, nz + 2
            do j = 1 - 2, ny + 2
                do i = 1 - 2, nx + 2
                    dsx = (i - isx)*dx
                    dsy = (j - isy)*dy
                    dsz = (k - isz)*dz
                    t0(i, j, k) = sqrt(dsx**2 + dsy**2 + dsz**2)/vsource
                    pdxt0(i, j, k) = dsx/vsource**2/t0(i, j, k)
                    pdyt0(i, j, k) = dsy/vsource**2/t0(i, j, k)
                    pdzt0(i, j, k) = dsz/vsource**2/t0(i, j, k)
                end do
            end do
        end do
        t0(isx, isy, isz) = 0.0
        pdxt0(isx, isy, isz) = 0.0
        pdyt0(isx, isy, isz) = 0.0
        pdzt0(isx, isy, isz) = 0.0

        ! Multiplicative time field
        call alloc_array(unknown, [1, nx, 1, ny, 1, nz], pad=2)
        call alloc_array(tt, [1, nx, 1, ny, 1, nz], pad=2)
        tt = sqrt(float_huge)
        tt(isx, isy, isz) = 1.0
        unknown = .true.
        unknown(isx, isy, isz) = .false.

        ! Initialize wavefront heap structure
        call wavefront%init(nx*ny*nz, 4, smaller_traveltime)

        i = isx
        j = isy
        k = isz

        l = 1
        ! The source point is known, therefore l < nx*ny*nz suffices
        do while (l < nx*ny*nz)

            ! Compute traveltiem at neighbour points
            do hk = -1, 1
                do hj = -1, 1
                    do hi = -1, 1
                        ii = clip(i + hi, 1, nx)
                        jj = clip(j + hj, 1, ny)
                        kk = clip(k + hk, 1, nz)
                        if (unknown(ii, jj, kk) .and. abs(hi) + abs(hj) + abs(hk) == 1) then
                            v2 = vp(ii, jj, kk)**2
                            call solve_quadratic(ii, jj, kk, dx, dy, dz, v2, v2, v2, 1.0)
                            call wavefront%insert([t0(ii, jj, kk)*tt(ii, jj, kk), real(ii), real(jj), real(kk)])
                        end if
                    end do
                end do
            end do

            ! Remove the smallest traveltime
            call wavefront%pop(node)
            i = int(node(2))
            j = int(node(3))
            k = int(node(4))

            l = l + 1

        end do

        call wavefront%free

        t = t0(1:nx, 1:ny, 1:nz)*tt(1:nx, 1:ny, 1:nz) + geom%srcr(isf)%t0

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                iry = nint((geom%recr(i)%y - oy)/dy) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                trec(i, 1) = t(irx, iry, irz) + geom%recr(i)%t0
            end if
        end do

        ! Transpose to output
        t = permute(t, 321)

    end subroutine forward_iso_fast_march

    !
    !> Solve the fast sweep scheme for the factorized eikonal equation
    !

#ifdef parallel_sweeping

    !
    !> Parallel fast sweeping implementing
    !> Detrixhe et al., 2013, JCP,
    !> A parallel fast sweeping method for the Eikonal equation
    !> doi: 10.1016/j.jcp.2012.11.042
    !> with adaptions and modifications for 3D scenario
    !
    subroutine fast_sweep_forward(nx, ny, nz, nxd, nyd, nzd, dx, dy, dz, vp, t0, px0, py0, pz0, tau)

        integer, intent(in) :: nx, ny, nz, nxd, nyd, nzd
        real, intent(in) :: dx, dy, dz
        real, dimension(:, :, :), intent(in) :: vp, t0, px0, py0, pz0
        real, dimension(:, :, :), intent(inout) :: tau

        integer :: nxbeg, nxend, nybeg, nyend, nzbeg, nzend
        integer :: i, j, k
        double precision :: signx, signy, signz
        double precision :: root
        double precision :: dxdt, dydt, dzdt
        double precision :: taux, tauy, tauz, t0x, t0y, t0z
        double precision :: tauc, t0c
        double precision :: taucx, taucy, taucz
        double precision :: px0c, py0c, pz0c
        double precision :: v2
        double precision :: u(3), d(3), array_tauc(3), array_tau(3), array_p0(3), array_t0(3), array_sign(3)
        double precision :: da, db
        double precision :: signa, signb
        double precision :: dadt, dbdt
        double precision :: taua, taub, t0a, t0b
        double precision :: tauca, taucb
        double precision :: pa0c, pb0c
        double precision :: travel
        integer :: i1, i2, j1, j2, k1, k2
        double precision :: dada, dbdb, dxdx, dydy, dzdz
        integer :: level, kbeg, kend, jbeg, jend

        if (nxd > 0) then
            nxbeg = 1
            nxend = nx
        else
            nxbeg = nx
            nxend = 1
        end if

        if (nyd > 0) then
            nybeg = 1
            nyend = ny
        else
            nybeg = ny
            nyend = 1
        end if

        if (nzd > 0) then
            nzbeg = 1
            nzend = nz
        else
            nzbeg = nz
            nzend = 1
        end if

        dxdx = dx**2
        dydy = dy**2
        dzdz = dz**2

        ! In 3D, the maximum level is nz + (ny - 1) + (nx - 1)
        do level = 1, nx + ny + nz - 2

            kbeg = ifthen(level <= nx + ny, nzbeg, nzbeg + (level - (nx + ny - 1))*nzd)
            kend = ifthen(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            jbeg = ifthen(level <= nx + ny, nybeg, nybeg + (level - (nx + ny - 1))*nyd)
            jend = ifthen(level <= ny, nybeg + (level - 1)*nyd, nyend)

            !$omp parallel do private(i, j, k, &
                !$omp signx, signy, signz, &
                !$omp root, &
                !$omp dxdt, dydt, dzdt, &
                !$omp taux, tauy, tauz, t0x, t0y, t0z, &
                !$omp tauc, t0c, &
                !$omp taucx, taucy, taucz, &
                !$omp px0c, py0c, pz0c, &
                !$omp v2, &
                !$omp u, d, array_tauc, array_tau, array_p0, array_t0, array_sign, &
                !$omp da, db, &
                !$omp signa, signb, &
                !$omp dadt, dbdt, &
                !$omp taua, taub, t0a, t0b, &
                !$omp tauca, taucb, &
                !$omp pa0c, pb0c, &
                !$omp travel, &
                !$omp i1, i2, j1, j2, k1, k2, &
                !$omp dada, dbdb) collapse(2) schedule(dynamic)
            do k = kbeg, kend, nzd
                do j = jbeg, jend, nyd

                    ! The condition is
                    ! abs(i - nxbeg) + abs(j - nybeg) + abs(k - nzbeg) + 1 = level
                    i = ((level - 1) - abs(k - nzbeg) - abs(j - nybeg))/nxd + nxbeg

                    ! Different from 2D, in 3D, the above condition can satisfy with out-of-range i values.
                    ! This can potentially cause load unbalance, but is the best I can do at this moment.
                    if (i < 1 .or. i > nx) then
                        cycle
                    end if

                    v2 = vp(i, j, k)**2
                    tauc = tau(i, j, k)
                    t0c = t0(i, j, k)
                    px0c = px0(i, j, k)
                    py0c = py0(i, j, k)
                    pz0c = pz0(i, j, k)

                    if (i == 1) then

                        i1 = i
                        i2 = i + 1

                        taux = tau(i2, j, k)
                        t0x = t0(i2, j, k)
                        signx = -1.0

                    else if (i == nx) then

                        i1 = i - 1
                        i2 = i

                        taux = tau(i1, j, k)
                        t0x = t0(i1, j, k)
                        signx = 1.0

                    else

                        i1 = i - 1
                        i2 = i + 1

                        if (tau(i1, j, k)*t0(i1, j, k) <= tau(i2, j, k)*t0(i2, j, k)) then
                            taux = tau(i1, j, k)
                            t0x = t0(i1, j, k)
                            signx = 1.0
                        else
                            taux = tau(i2, j, k)
                            t0x = t0(i2, j, k)
                            signx = -1.0
                        end if

                    end if

                    if (j == 1) then

                        j1 = j
                        j2 = j + 1

                        tauy = tau(i, j2, k)
                        t0y = t0(i, j2, k)
                        signy = -1.0

                    else if (j == ny) then

                        j1 = j - 1
                        j2 = j

                        tauy = tau(i, j1, k)
                        t0y = t0(i, j1, k)
                        signy = 1.0

                    else

                        j1 = j - 1
                        j2 = j + 1

                        if (tau(i, j1, k)*t0(i, j1, k) <= tau(i, j2, k)*t0(i, j2, k)) then
                            tauy = tau(i, j1, k)
                            t0y = t0(i, j1, k)
                            signy = 1.0
                        else
                            tauy = tau(i, j2, k)
                            t0y = t0(i, j2, k)
                            signy = -1.0
                        end if

                    end if

                    if (k == 1) then

                        k1 = k
                        k2 = k + 1

                        tauz = tau(i, j, k2)
                        t0z = t0(i, j, k2)
                        signz = -1.0

                    else if (k == nz) then

                        k1 = k - 1
                        k2 = k

                        tauz = tau(i, j, k1)
                        t0z = t0(i, j, k1)
                        signz = 1.0

                    else

                        k1 = k - 1
                        k2 = k + 1

                        if (tau(i, j, k1)*t0(i, j, k1) <= tau(i, j, k2)*t0(i, j, k2)) then
                            tauz = tau(i, j, k1)
                            t0z = t0(i, j, k1)
                            signz = 1.0
                        else
                            tauz = tau(i, j, k2)
                            t0z = t0(i, j, k2)
                            signz = -1.0
                        end if

                    end if

                    if (taux == float_huge .and. tauy == float_huge .and. tauz == float_huge) then
                        cycle
                    end if

                    taucx = (t0c*taux + dx/vp(i, j, k))/(t0c + abs(px0c)*dx)
                    taucy = (t0c*tauy + dy/vp(i, j, k))/(t0c + abs(py0c)*dy)
                    taucz = (t0c*tauz + dz/vp(i, j, k))/(t0c + abs(pz0c)*dz)

                    u(1) = min(tau(i1, j, k)*t0(i1, j, k), tau(i2, j, k)*t0(i2, j, k))
                    u(2) = min(tau(i, j1, k)*t0(i, j1, k), tau(i, j2, k)*t0(i, j2, k))
                    u(3) = min(tau(i, j, k1)*t0(i, j, k1), tau(i, j, k2)*t0(i, j, k2))
                    d = [dx, dy, dz]
                    array_tauc = [taucx, taucy, taucz]
                    array_p0 = [px0c, py0c, pz0c]
                    array_tau = [taux, tauy, tauz]
                    array_t0 = [t0x, t0y, t0z]
                    array_sign = [signx, signy, signz]

                    if (u(1) == float_huge .and. u(2) == float_huge .and. u(3) == float_huge) then
                        cycle
                    end if

                    if (u(1) > u(3)) then
                        call swap(u(1), u(3))
                        call swap(d(1), d(3))
                        call swap(array_tau(1), array_tau(3))
                        call swap(array_t0(1), array_t0(3))
                        call swap(array_sign(1), array_sign(3))
                        call swap(array_p0(1), array_p0(3))
                        call swap(array_tauc(1), array_tauc(3))
                    end if
                    if (u(1) > u(2)) then
                        call swap(u(1), u(2))
                        call swap(d(1), d(2))
                        call swap(array_tau(1), array_tau(2))
                        call swap(array_t0(1), array_t0(2))
                        call swap(array_sign(1), array_sign(2))
                        call swap(array_p0(1), array_p0(2))
                        call swap(array_tauc(1), array_tauc(2))
                    end if
                    if (u(2) > u(3)) then
                        call swap(u(2), u(3))
                        call swap(d(2), d(3))
                        call swap(array_tau(2), array_tau(3))
                        call swap(array_t0(2), array_t0(3))
                        call swap(array_sign(2), array_sign(3))
                        call swap(array_p0(2), array_p0(3))
                        call swap(array_tauc(2), array_tauc(3))
                    end if

                    travel = u(1) + d(1)/vp(i, j, k)
                    tauc = travel/t0(i, j, k)
                    if (travel > u(2)) then
                        ! Consider the 2D plane case
                        ! Here, new notations (a, b) are used

                        da = d(1)
                        db = d(2)
                        dada = da**2
                        dbdb = db**2
                        taua = array_tau(1)
                        taub = array_tau(2)
                        t0a = array_t0(1)
                        t0b = array_t0(2)
                        signa = array_sign(1)
                        signb = array_sign(2)
                        pa0c = array_p0(1)
                        pb0c = array_p0(2)
                        tauca = array_tauc(1)
                        taucb = array_tauc(2)

                        ! Sovling the quadratic equation analbtically
                        root = (da*dbdb*pa0c*signa*t0c*taua*v2 + dbdb*signa**2*t0c**2*taua*v2 + &
                            dada*(signb*t0c*(db*pb0c + signb*t0c)*taub*v2 + &
                            dbdb*sqrt((v2*(dada*dbdb*(pa0c**2 + pb0c**2) + 2*da*dbdb*pa0c*signa*t0c + 2*dada*db*pb0c*signb*t0c + &
                            dbdb*signa**2*t0c**2 + dada*signb**2*t0c**2 - &
                            (signa*t0c*(db*pb0c + signb*t0c)*taua - signb*t0c*(da*pa0c + signa*t0c)*taub)**2*v2))/(dada*dbdb))))/ &
                            ((dada*dbdb*(pa0c**2 + pb0c**2) + 2*da*db*(db*pa0c*signa + da*pb0c*signb)*t0c + &
                            (dbdb*signa**2 + dada*signb**2)*t0c**2)*v2)

                        dadt = (root*t0c - taua*t0a)/(signa*da)
                        dbdt = (root*t0c - taub*t0b)/(signb*db)

                        if (dadt*signa > 0 .and. dbdt*signb > 0) then
                            tauc = root
                        else
                            if (tauca*t0a < taucb*t0b) then
                                tauc = tauca
                            else
                                tauc = taucb
                            end if
                        end if

                        travel = tauc*t0(i, j, k)
                        if (travel > u(3)) then
                            ! Consider the 3D cuboid case
                            ! In this case, the order of u no long matters, and for convenience
                            ! here I just use the x-y-z notation to avoid more new symbols

                            ! Sovling the quadratic equation analytically
                            root = (dx*dydy*dzdz*px0c*signx*t0c*taux*v2 + dydy*dzdz*signx**2*t0c**2*taux*v2 + &
                                dxdx*t0c*(dzdz*signy*(dy*py0c + signy*t0c)*tauy + dydy*signz*(dz*pz0c + signz*t0c)*tauz)*v2 + &
                                dxdx*dydy*dzdz*sqrt(v2*(((dzdz*t0c*(dydy*signx*(dx*px0c + signx*t0c)*taux + &
                                dxdx*signy*(dy*py0c + signy*t0c)*tauy) + dxdx*dydy*signz*t0c*(dz*pz0c + signz*t0c)*tauz)**2*v2)/ &
                                (dxdx*dydy*dzdz)**2 - (px0c**2 + py0c**2 + pz0c**2 + (2*px0c*signx*t0c)/dx + &
                                t0c*((2*py0c*signy)/dy + (signx**2*t0c)/dxdx + (signy**2*t0c)/dydy + (signz*(2*dz*pz0c + signz*t0c))/dzdz))* &
                                (-1 + t0c**2*((signx**2*taux**2)/dxdx + (signy**2*tauy**2)/dydy + (signz**2*tauz**2)/dzdz)*v2))))/ &
                                ((dxdx*dydy*dzdz*(px0c**2 + py0c**2 + pz0c**2) + &
                                2*dx*dy*dz*(dy*dz*px0c*signx + dx*dz*py0c*signy + dx*dy*pz0c*signz)*t0c + &
                                (dxdx*dzdz*signy**2 + dydy*(dzdz*signx**2 + dxdx*signz**2))*t0c**2)*v2)

                            dxdt = (root*t0c - taux*t0x)/(signx*dx)
                            dydt = (root*t0c - tauy*t0y)/(signy*dy)
                            dzdt = (root*t0c - tauz*t0z)/(signz*dz)

                            if (dxdt*signx > 0 .and. dydt*signy > 0 .and. dzdt*signz > 0) then
                                tauc = root
                            else
                                if (taucx*t0x < taucz*t0z .and. taucx*t0x < taucy*t0y) then
                                    tauc = taucx
                                else if (taucy*t0y < taucz*t0z .and. taucy*t0y < taucx*t0x) then
                                    tauc = taucy
                                else
                                    tauc = taucz
                                end if
                            end if

                        end if
                    end if

                    tau(i, j, k) = min(tau(i, j, k), tauc)

                end do
            end do
            !$omp end parallel do

        end do

    end subroutine fast_sweep_forward

    !
    !> Fast sweeping to solve the adjoint state equation
    !
    subroutine fast_sweep_adjoint(n1, n2, n3, n1d, n2d, n3d, d1, d2, d3, lambda)

        integer, intent(in) :: n1, n2, n3, n1d, n2d, n3d
        real, intent(in) :: d1, d2, d3
        real, dimension(:, :, :), intent(inout) :: lambda

        integer :: i, j, k
        double precision :: app, amp, apm, amm
        double precision :: bpp, bmp, bpm, bmm
        double precision :: cpp, cmp, cpm, cmm
        double precision :: ap, am, bp, bm, cp, cm
        double precision :: lhs, rhs, t
        integer :: n1beg, n1end, n2beg, n2end, n3beg, n3end
        integer :: i1, i2, j1, j2, k1, k2
        integer :: level, kbeg, kend, jbeg, jend

        if (n1d > 0) then
            n1beg = 1
            n1end = n1
        else
            n1beg = n1
            n1end = 1
        end if

        if (n2d > 0) then
            n2beg = 1
            n2end = n2
        else
            n2beg = n2
            n2end = 1
        end if

        if (n3d > 0) then
            n3beg = 1
            n3end = n3
        else
            n3beg = n3
            n3end = 1
        end if

        ! In 3D, the maximum level is n3 + (n2 - 1) + (n1 - 1)
        do level = 1, n1 + n2 + n3 - 2

            kbeg = ifthen(level <= n1 + n2, n3beg, n3beg + (level - (n1 + n2 - 1))*n3d)
            kend = ifthen(level <= n3, n3beg + (level - 1)*n3d, n3end)

            jbeg = ifthen(level <= n1 + n2, n2beg, n2beg + (level - (n1 + n2 - 1))*n2d)
            jend = ifthen(level <= n2, n2beg + (level - 1)*n2d, n2end)

            !$omp parallel do private(i, j, k, &
                !$omp app, amp, apm, amm, &
                !$omp bpp, bmp, bpm, bmm, &
                !$omp cpp, cmp, cpm, cmm, &
                !$omp ap, am, bp, bm, cp, cm, &
                !$omp lhs, rhs, t) collapse(2) schedule(dynamic)
            do k = kbeg, kend, n3d
                do j = jbeg, jend, n2d

                    ! The condition is
                    ! abs(i - n1beg) + abs(j - n2beg) + abs(k - n3beg) + 1 = level
                    i = ((level - 1) - abs(k - n3beg) - abs(j - n2beg))/n1d + n1beg

                    ! Different from 2D, in 3D, the above condition can satisfy with out-of-range i values.
                    ! This can potentially cause load unbalance, but is the best I can do at this moment.
                    if (i < 1 .or. i > n1) then
                        cycle
                    end if

                    if (.not. recrflag(i, j, k)) then

                        if (i == 1) then
                            i1 = i
                        else
                            i1 = i - 1
                        end if
                        if (i == n1) then
                            i2 = i
                        else
                            i2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = j
                        else
                            j1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = j
                        else
                            j2 = j + 1
                        end if

                        if (k == 1) then
                            k1 = k
                        else
                            k1 = k - 1
                        end if
                        if (k == n3) then
                            k2 = k
                        else
                            k2 = k + 1
                        end if

                        ! Solve equation (A-9) in Taillander et al (2009)
                        ap = (tt(i2, j, k) - tt(i, j, k))/d1
                        am = (tt(i, j, k) - tt(i1, j, k))/d1

                        bp = (tt(i, j2, k) - tt(i, j, k))/d2
                        bm = (tt(i, j, k) - tt(i, j1, k))/d2

                        cp = (tt(i, j, k2) - tt(i, j, k))/d3
                        cm = (tt(i, j, k) - tt(i, j, k1))/d3

                        app = (ap + abs(ap))/2.0
                        apm = (ap - abs(ap))/2.0

                        amp = (am + abs(am))/2.0
                        amm = (am - abs(am))/2.0

                        bpp = (bp + abs(bp))/2.0
                        bpm = (bp - abs(bp))/2.0

                        bmp = (bm + abs(bm))/2.0
                        bmm = (bm - abs(bm))/2.0

                        cpp = (cp + abs(cp))/2.0
                        cpm = (cp - abs(cp))/2.0

                        cmp = (cm + abs(cm))/2.0
                        cmm = (cm - abs(cm))/2.0

                        lhs = (apm - amp)/d1 + (bpm - bmp)/d2 + (cpm - cmp)/d3
                        rhs = (amm*lambda(i1, j, k) - app*lambda(i2, j, k))/d1 &
                            + (bmm*lambda(i, j1, k) - bpp*lambda(i, j2, k))/d2 &
                            + (cmm*lambda(i, j, k1) - cpp*lambda(i, j, k2))/d3

                        if (lhs == 0) then
                            t = 0
                        else
                            t = rhs/lhs
                        end if

                        lambda(i, j, k) = min(lambda(i, j, k), t)

                    end if

                end do
            end do
        end do

    end subroutine fast_sweep_adjoint

#else

    !
    !> Serial fast sweeping
    !
    subroutine fast_sweep_forward(nx, ny, nz, nxd, nyd, nzd, dx, dy, dz, vp, t0, px0, py0, pz0, tau)

        integer, intent(in) :: nx, ny, nz, nxd, nyd, nzd
        real, intent(in) :: dx, dy, dz
        real, dimension(:, :, :), intent(in) :: vp, t0, px0, py0, pz0
        real, dimension(:, :, :), intent(inout) :: tau

        integer :: nxbeg, nxend, nybeg, nyend, nzbeg, nzend
        integer :: i, j, k
        double precision :: signx, signy, signz
        double precision :: root
        double precision :: dxdt, dydt, dzdt
        double precision :: taux, tauy, tauz, t0x, t0y, t0z
        double precision :: tauc, t0c
        double precision :: taucx, taucy, taucz
        double precision :: px0c, py0c, pz0c
        double precision :: v2
        double precision :: u(3), d(3), array_tauc(3), array_tau(3), array_p0(3), array_t0(3), array_sign(3)
        double precision :: da, db !, dc
        double precision :: signa, signb !, signc
        double precision :: dadt, dbdt !, dcdt
        double precision :: taua, taub, t0a, t0b
        double precision :: tauca, taucb !, taucd
        double precision :: pa0c, pb0c !, pc0c
        double precision :: travel
        integer :: i1, i2, j1, j2, k1, k2
        double precision :: dada, dbdb, dxdx, dydy, dzdz

        if (nxd > 0) then
            nxbeg = 1
            nxend = nx
        else
            nxbeg = nx
            nxend = 1
        end if

        if (nyd > 0) then
            nybeg = 1
            nyend = ny
        else
            nybeg = ny
            nyend = 1
        end if

        if (nzd > 0) then
            nzbeg = 1
            nzend = nz
        else
            nzbeg = nz
            nzend = 1
        end if

        dxdx = dx**2
        dydy = dy**2
        dzdz = dz**2

        do k = nzbeg, nzend, nzd
            do j = nybeg, nyend, nyd
                do i = nxbeg, nxend, nxd

                    v2 = vp(i, j, k)**2
                    tauc = tau(i, j, k)
                    t0c = t0(i, j, k)
                    px0c = px0(i, j, k)
                    py0c = py0(i, j, k)
                    pz0c = pz0(i, j, k)

                    if (i == 1) then

                        i1 = i
                        i2 = i + 1

                        taux = tau(i2, j, k)
                        t0x = t0(i2, j, k)
                        signx = -1.0

                    else if (i == nx) then

                        i1 = i - 1
                        i2 = i

                        taux = tau(i1, j, k)
                        t0x = t0(i1, j, k)
                        signx = 1.0

                    else

                        i1 = i - 1
                        i2 = i + 1

                        if (tau(i1, j, k)*t0(i1, j, k) <= tau(i2, j, k)*t0(i2, j, k)) then
                            taux = tau(i1, j, k)
                            t0x = t0(i1, j, k)
                            signx = 1.0
                        else
                            taux = tau(i2, j, k)
                            t0x = t0(i2, j, k)
                            signx = -1.0
                        end if

                    end if

                    if (j == 1) then

                        j1 = j
                        j2 = j + 1

                        tauy = tau(i, j2, k)
                        t0y = t0(i, j2, k)
                        signy = -1.0

                    else if (j == ny) then

                        j1 = j - 1
                        j2 = j

                        tauy = tau(i, j1, k)
                        t0y = t0(i, j1, k)
                        signy = 1.0

                    else

                        j1 = j - 1
                        j2 = j + 1

                        if (tau(i, j1, k)*t0(i, j1, k) <= tau(i, j2, k)*t0(i, j2, k)) then
                            tauy = tau(i, j1, k)
                            t0y = t0(i, j1, k)
                            signy = 1.0
                        else
                            tauy = tau(i, j2, k)
                            t0y = t0(i, j2, k)
                            signy = -1.0
                        end if

                    end if

                    if (k == 1) then

                        k1 = k
                        k2 = k + 1

                        tauz = tau(i, j, k2)
                        t0z = t0(i, j, k2)
                        signz = -1.0

                    else if (k == nz) then

                        k1 = k - 1
                        k2 = k

                        tauz = tau(i, j, k1)
                        t0z = t0(i, j, k1)
                        signz = 1.0

                    else

                        k1 = k - 1
                        k2 = k + 1

                        if (tau(i, j, k1)*t0(i, j, k1) <= tau(i, j, k2)*t0(i, j, k2)) then
                            tauz = tau(i, j, k1)
                            t0z = t0(i, j, k1)
                            signz = 1.0
                        else
                            tauz = tau(i, j, k2)
                            t0z = t0(i, j, k2)
                            signz = -1.0
                        end if

                    end if

                    if (taux == float_huge .and. tauy == float_huge .and. tauz == float_huge) then
                        cycle
                    end if

                    taucx = (t0c*taux + dx/vp(i, j, k))/(t0c + abs(px0c)*dx)
                    taucy = (t0c*tauy + dy/vp(i, j, k))/(t0c + abs(py0c)*dy)
                    taucz = (t0c*tauz + dz/vp(i, j, k))/(t0c + abs(pz0c)*dz)

                    u(1) = min(tau(i1, j, k)*t0(i1, j, k), tau(i2, j, k)*t0(i2, j, k))
                    u(2) = min(tau(i, j1, k)*t0(i, j1, k), tau(i, j2, k)*t0(i, j2, k))
                    u(3) = min(tau(i, j, k1)*t0(i, j, k1), tau(i, j, k2)*t0(i, j, k2))
                    d = [dx, dy, dz]
                    array_tauc = [taucx, taucy, taucz]
                    array_p0 = [px0c, py0c, pz0c]
                    array_tau = [taux, tauy, tauz]
                    array_t0 = [t0x, t0y, t0z]
                    array_sign = [signx, signy, signz]

                    if (u(1) == float_huge .and. u(2) == float_huge .and. u(3) == float_huge) then
                        cycle
                    end if

                    if (u(1) > u(3)) then
                        call swap(u(1), u(3))
                        call swap(d(1), d(3))
                        call swap(array_tau(1), array_tau(3))
                        call swap(array_t0(1), array_t0(3))
                        call swap(array_sign(1), array_sign(3))
                        call swap(array_p0(1), array_p0(3))
                        call swap(array_tauc(1), array_tauc(3))
                    end if
                    if (u(1) > u(2)) then
                        call swap(u(1), u(2))
                        call swap(d(1), d(2))
                        call swap(array_tau(1), array_tau(2))
                        call swap(array_t0(1), array_t0(2))
                        call swap(array_sign(1), array_sign(2))
                        call swap(array_p0(1), array_p0(2))
                        call swap(array_tauc(1), array_tauc(2))
                    end if
                    if (u(2) > u(3)) then
                        call swap(u(2), u(3))
                        call swap(d(2), d(3))
                        call swap(array_tau(2), array_tau(3))
                        call swap(array_t0(2), array_t0(3))
                        call swap(array_sign(2), array_sign(3))
                        call swap(array_p0(2), array_p0(3))
                        call swap(array_tauc(2), array_tauc(3))
                    end if

                    travel = u(1) + d(1)/vp(i, j, k)
                    tauc = travel/t0(i, j, k)
                    if (travel > u(2)) then
                        ! Consider the 2D plane case
                        ! Here, new notations (a, b) are used

                        da = d(1)
                        db = d(2)
                        dada = da**2
                        dbdb = db**2
                        taua = array_tau(1)
                        taub = array_tau(2)
                        t0a = array_t0(1)
                        t0b = array_t0(2)
                        signa = array_sign(1)
                        signb = array_sign(2)
                        pa0c = array_p0(1)
                        pb0c = array_p0(2)
                        tauca = array_tauc(1)
                        taucb = array_tauc(2)

                        ! Sovling the quadratic equation analbtically
                        root = (da*dbdb*pa0c*signa*t0c*taua*v2 + dbdb*signa**2*t0c**2*taua*v2 + &
                            dada*(signb*t0c*(db*pb0c + signb*t0c)*taub*v2 + &
                            dbdb*sqrt((v2*(dada*dbdb*(pa0c**2 + pb0c**2) + 2*da*dbdb*pa0c*signa*t0c + 2*dada*db*pb0c*signb*t0c + &
                            dbdb*signa**2*t0c**2 + dada*signb**2*t0c**2 - &
                            (signa*t0c*(db*pb0c + signb*t0c)*taua - signb*t0c*(da*pa0c + signa*t0c)*taub)**2*v2))/(dada*dbdb))))/ &
                            ((dada*dbdb*(pa0c**2 + pb0c**2) + 2*da*db*(db*pa0c*signa + da*pb0c*signb)*t0c + &
                            (dbdb*signa**2 + dada*signb**2)*t0c**2)*v2)

                        dadt = (root*t0c - taua*t0a)/(signa*da)
                        dbdt = (root*t0c - taub*t0b)/(signb*db)

                        if (dadt*signa > 0 .and. dbdt*signb > 0) then
                            tauc = root
                        else
                            if (tauca*t0a < taucb*t0b) then
                                tauc = tauca
                            else
                                tauc = taucb
                            end if
                        end if

                        travel = tauc*t0(i, j, k)
                        if (travel > u(3)) then
                            ! Consider the 3D cuboid case
                            ! In this case, the order of u no long matters, and for convenience
                            ! here I just use the x-y-z notation to avoid more new symbols

                            ! Sovling the quadratic equation analytically
                            root = (dx*dydy*dzdz*px0c*signx*t0c*taux*v2 + dydy*dzdz*signx**2*t0c**2*taux*v2 + &
                                dxdx*t0c*(dzdz*signy*(dy*py0c + signy*t0c)*tauy + dydy*signz*(dz*pz0c + signz*t0c)*tauz)*v2 + &
                                dxdx*dydy*dzdz*sqrt(v2*(((dzdz*t0c*(dydy*signx*(dx*px0c + signx*t0c)*taux + &
                                dxdx*signy*(dy*py0c + signy*t0c)*tauy) + dxdx*dydy*signz*t0c*(dz*pz0c + signz*t0c)*tauz)**2*v2)/ &
                                (dxdx*dydy*dzdz)**2 - (px0c**2 + py0c**2 + pz0c**2 + (2*px0c*signx*t0c)/dx + &
                                t0c*((2*py0c*signy)/dy + (signx**2*t0c)/dxdx + (signy**2*t0c)/dydy + (signz*(2*dz*pz0c + signz*t0c))/dzdz))* &
                                (-1 + t0c**2*((signx**2*taux**2)/dxdx + (signy**2*tauy**2)/dydy + (signz**2*tauz**2)/dzdz)*v2))))/ &
                                ((dxdx*dydy*dzdz*(px0c**2 + py0c**2 + pz0c**2) + &
                                2*dx*dy*dz*(dy*dz*px0c*signx + dx*dz*py0c*signy + dx*dy*pz0c*signz)*t0c + &
                                (dxdx*dzdz*signy**2 + dydy*(dzdz*signx**2 + dxdx*signz**2))*t0c**2)*v2)

                            dxdt = (root*t0c - taux*t0x)/(signx*dx)
                            dydt = (root*t0c - tauy*t0y)/(signy*dy)
                            dzdt = (root*t0c - tauz*t0z)/(signz*dz)

                            if (dxdt*signx > 0 .and. dydt*signy > 0 .and. dzdt*signz > 0) then
                                tauc = root
                            else
                                if (taucx*t0x < taucz*t0z .and. taucx*t0x < taucy*t0y) then
                                    tauc = taucx
                                else if (taucy*t0y < taucz*t0z .and. taucy*t0y < taucx*t0x) then
                                    tauc = taucy
                                else
                                    tauc = taucz
                                end if
                            end if

                        end if
                    end if

                    tau(i, j, k) = min(tau(i, j, k), tauc)

                end do
            end do
        end do

    end subroutine fast_sweep_forward

    !
    !> Fast sweeping to solve the adjoint state equation
    !
    subroutine fast_sweep_adjoint(n1, n2, n3, n1d, n2d, n3d, d1, d2, d3, lambda)

        integer, intent(in) :: n1, n2, n3, n1d, n2d, n3d
        real, intent(in) :: d1, d2, d3
        real, dimension(:, :, :), intent(inout) :: lambda

        integer :: i, j, k
        double precision :: app, amp, apm, amm
        double precision :: bpp, bmp, bpm, bmm
        double precision :: cpp, cmp, cpm, cmm
        double precision :: ap, am, bp, bm, cp, cm
        double precision :: lhs, rhs, t
        integer :: n1beg, n1end, n2beg, n2end, n3beg, n3end
        integer :: i1, i2, j1, j2, k1, k2

        if (n1d > 0) then
            n1beg = 1
            n1end = n1
        else
            n1beg = n1
            n1end = 1
        end if

        if (n2d > 0) then
            n2beg = 1
            n2end = n2
        else
            n2beg = n2
            n2end = 1
        end if

        if (n3d > 0) then
            n3beg = 1
            n3end = n3
        else
            n3beg = n3
            n3end = 1
        end if

        ! Sweep over finite-difference grids
        do k = n3beg, n3end, n3d
            do j = n2beg, n2end, n2d
                do i = n1beg, n1end, n1d

                    if (.not. recrflag(i, j, k)) then

                        if (i == 1) then
                            i1 = i
                        else
                            i1 = i - 1
                        end if
                        if (i == n1) then
                            i2 = i
                        else
                            i2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = j
                        else
                            j1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = j
                        else
                            j2 = j + 1
                        end if

                        if (k == 1) then
                            k1 = k
                        else
                            k1 = k - 1
                        end if
                        if (k == n3) then
                            k2 = k
                        else
                            k2 = k + 1
                        end if

                        ! Solve equation (A-9) in Taillander et al (2009)
                        ap = (tt(i2, j, k) - tt(i, j, k))/d1
                        am = (tt(i, j, k) - tt(i1, j, k))/d1

                        bp = (tt(i, j2, k) - tt(i, j, k))/d2
                        bm = (tt(i, j, k) - tt(i, j1, k))/d2

                        cp = (tt(i, j, k2) - tt(i, j, k))/d3
                        cm = (tt(i, j, k) - tt(i, j, k1))/d3

                        app = (ap + abs(ap))/2.0
                        apm = (ap - abs(ap))/2.0

                        amp = (am + abs(am))/2.0
                        amm = (am - abs(am))/2.0

                        bpp = (bp + abs(bp))/2.0
                        bpm = (bp - abs(bp))/2.0

                        bmp = (bm + abs(bm))/2.0
                        bmm = (bm - abs(bm))/2.0

                        cpp = (cp + abs(cp))/2.0
                        cpm = (cp - abs(cp))/2.0

                        cmp = (cm + abs(cm))/2.0
                        cmm = (cm - abs(cm))/2.0

                        lhs = (apm - amp)/d1 + (bpm - bmp)/d2 + (cpm - cmp)/d3
                        rhs = (amm*lambda(i1, j, k) - app*lambda(i2, j, k))/d1 &
                            + (bmm*lambda(i, j1, k) - bpp*lambda(i, j2, k))/d2 &
                            + (cmm*lambda(i, j, k1) - cpp*lambda(i, j, k2))/d3

                        if (lhs == 0) then
                            t = 0
                        else
                            t = rhs/lhs
                        end if

                        lambda(i, j, k) = min(lambda(i, j, k), t)

                    end if

                end do
            end do
        end do

    end subroutine fast_sweep_adjoint

#endif

    !
    !> Fast sweeping factorized eikonal solver
    !
    subroutine forward_iso_fast_sweep(v, d, o, geom, t, trec)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        real, allocatable, dimension(:, :, :) :: tt, tt_prev, vp
        real :: dx, dy, dz, ox, oy, oz
        integer :: nx, ny, nz, niter, i, j, k, l, isx, isy, isz, irx, iry, irz
        real, allocatable, dimension(:) :: t1, t2, t3, t4
        real :: ttdiff, vsource, dsx, dsy, dsz
        integer :: imin

        dx = d(1)
        dy = d(2)
        dz = d(3)
        ox = o(1)
        oy = o(2)
        oz = o(3)

        vp = permute(v, 321)
        nx = size(vp, 1)
        ny = size(vp, 2)
        nz = size(vp, 3)

        t0 = zeros(nx, ny, nz)
        pdxt0 = zeros(nx, ny, nz)
        pdyt0 = zeros(nx, ny, nz)
        pdzt0 = zeros(nx, ny, nz)

        ! Source location
        t1 = zeros(geom%ns) + sqrt(float_huge)
        t2 = zeros_like(t1)
        t3 = zeros_like(t1)
        t4 = zeros_like(t1)

        ! Background time field
        !$omp parallel do private(i, j, k, l, isx, isy, isz, vsource, dsx, dsy, dsz, imin, t1, t2, t3, t4)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    do l = 1, geom%ns

                        if (geom%srcr(l)%amp == 1) then

                            isx = nint((geom%srcr(l)%x - ox)/dx) + 1
                            isy = nint((geom%srcr(l)%y - oy)/dy) + 1
                            isz = nint((geom%srcr(l)%z - oz)/dz) + 1
                            vsource = vp(isx, isy, isz)
                            dsx = (i - isx)*dx
                            dsy = (j - isy)*dy
                            dsz = (k - isz)*dz

                            t1(l) = sqrt(dsx**2 + dsy**2 + dsz**2)/vsource
                            if (t1(l) == 0) then
                                t2(l) = 0
                                t3(l) = 0
                                t4(l) = 0
                            else
                                t2(l) = dsx/vsource**2/t1(l)
                                t3(l) = dsy/vsource**2/t1(l)
                                t4(l) = dsz/vsource**2/t1(l)
                            end if

                            t1(l) = t1(l) + geom%srcr(l)%t0

                        end if

                    end do

                    imin = as_scalar(minloc(t1))
                    t0(i, j, k) = t1(imin)
                    pdxt0(i, j, k) = t2(imin)
                    pdyt0(i, j, k) = t3(imin)
                    pdzt0(i, j, k) = t4(imin)

                end do
            end do
        end do
        !$omp end parallel do

        tt_prev = zeros(nx, ny, nz)
        tt = zeros(nx, ny, nz) + sqrt(float_huge)
        !$omp parallel do private(i, isx, isy, isz)
        do i = 1, geom%ns
            if (geom%srcr(i)%amp == 1) then
                isx = nint((geom%srcr(i)%x - ox)/dx) + 1
                isy = nint((geom%srcr(i)%y - oy)/dy) + 1
                isz = nint((geom%srcr(i)%z - oz)/dz) + 1
                tt(isx, isy, isz) = 1.0
            end if
        end do
        !$omp end parallel do

        ! Fast sweeping via Gauss-Seidel iterations
        niter = 0
        ttdiff = sqrt(float_huge)
        do while (ttdiff >= sweep_stop_threshold .and. niter < sweep_niter_max)

            tt_prev = tt
            call fast_sweep_forward(nx, ny, nz, 1, 1, 1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, 1, 1, -1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, 1, -1, 1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, 1, -1, -1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, 1, 1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, 1, -1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, -1, 1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, -1, -1, dx, dy, dz, vp, t0, pdxt0, pdyt0, pdzt0, tt)

            ttdiff = mean(abs(tt - tt_prev))
            niter = niter + 1

        end do

        t = t0*tt

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        !$omp parallel do private(i, irx, iry, irz)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                iry = nint((geom%recr(i)%y - oy)/dy) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                trec(i, 1) = t(irx, iry, irz) + geom%recr(i)%t0
            end if
        end do
        !$omp end parallel do

        t = permute(t, 321)

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(niter)// &
            ', relative diff = '//num2str(ttdiff, '(es)'))

    end subroutine forward_iso_fast_sweep

    !
    !> Compute adjoint field for FATT
    !
    subroutine adjoint_iso(v, d, o, geom, t, tresidual, tadj)

        real, dimension(:, :, :), intent(in) :: v, t
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :, :), allocatable, intent(out) :: tadj

        real, allocatable, dimension(:, :, :) :: vp, lambda_prev
        integer :: nx, ny, nz, iter, i, irx, iry, irz
        real :: lambda_diff, dx, dy, dz, ox, oy, oz

        dx = d(1)
        dy = d(2)
        dz = d(3)
        ox = o(1)
        oy = o(2)
        oz = o(2)

        tt = permute(t, 321)
        vp = permute(v, 321)
        nx = size(tt, 1)
        ny = size(tt, 2)
        nz = size(tt, 3)

        recrflag = falses(nx, ny, nz)
        lambda = zeros(nx, ny, nz) + sqrt(float_huge)
        lambda_prev = zeros(nx, ny, nz)
        lambda_diff = sqrt(float_huge)

        ! Get the misfit values at the receivers
        !$omp parallel do private(i, irx, iry, irz)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                iry = nint((geom%recr(i)%y - oy)/dy) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                lambda(irx, iry, irz) = tresidual(i, 1)
                recrflag(irx, iry, irz) = .true.
            end if
        end do
        !$omp end parallel do

        ! Fast sweep to solve the adjoint-state equation
        iter = 0
        do while (lambda_diff >= sweep_stop_threshold .and. iter < sweep_niter_max)

            ! Save previous step
            lambda_prev = lambda

            ! Fast sweeps
            call fast_sweep_adjoint(nx, ny, nz, 1, 1, 1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, 1, 1, -1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, 1, -1, 1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, 1, -1, -1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, 1, 1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, 1, -1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, -1, 1, dx, dy, dz, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, -1, -1, dx, dy, dz, lambda)

            ! Check threshold
            lambda_diff = mean(abs(lambda - lambda_prev))

            iter = iter + 1

        end do

        tadj = permute(lambda, 321)

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(iter)// &
            ', relative diff = '//num2str(lambda_diff, '(es)'))

    end subroutine adjoint_iso

end module traveltime_iso
