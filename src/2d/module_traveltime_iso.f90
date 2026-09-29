!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

module traveltime_iso

    ! To ensure accuracy, computations in the internal subroutines/functions of this module
    ! are double-precision, but the input/output variables are single-precision.

    use libflit
    use parameters
    use utility
    use omp_lib

    implicit none

    real, parameter :: huge_value = sqrt(float_huge)

    double precision, allocatable, dimension(:, :) :: t0, pdxt0, pdzt0, tt, lambda
    logical, allocatable, dimension(:, :) :: unknown, recrflag

    private
    public :: forward_iso
    public :: adjoint_iso

contains

#ifdef legacy_solver

    !
    !> Parallel fast sweeping implementing
    !> Detrixhe et al., 2013, JCP,
    !> A parallel fast sweeping method for the Eikonal equation
    !> doi: 10.1016/j.jcp.2012.11.042
    !> with modifications
    !
    subroutine fast_sweep_forward(nx, nz, nxd, nzd, dx, dz, vp, t0, px0, pz0, tau)

        integer, intent(in) :: nx, nz, nxd, nzd
        double precision, intent(in) :: dx, dz
        double precision, dimension(:, :), intent(in) :: vp, t0, px0, pz0
        double precision, dimension(:, :), intent(inout) :: tau

        integer :: nxbeg, nxend, nzbeg, nzend
        integer :: i, j
        double precision :: signx, signz
        double precision :: root
        double precision :: taux, tauz, t0x, t0z
        double precision :: tauc, t0c
        double precision :: taucx, taucz
        double precision :: px0c, pz0c
        double precision :: v2
        double precision :: u(2), d(2), array_tauc(2), array_tau(2), array_p0(2), array_t0(2), array_sign(2)
        double precision :: da, db
        double precision :: signa, signb
        double precision :: dadt, dbdt
        double precision :: taua, taub, t0a, t0b
        double precision :: tauca, taucb
        double precision :: pa0c, pb0c
        double precision :: travel
        integer :: i1, i2, j1, j2
        double precision :: dada, dbdb
        integer :: level, jbeg, jend

        if (nxd > 0) then
            nxbeg = 1
            nxend = nx
        else
            nxbeg = nx
            nxend = 1
        end if

        if (nzd > 0) then
            nzbeg = 1
            nzend = nz
        else
            nzbeg = nz
            nzend = 1
        end if

        do level = 1, nx + nz - 1

            jbeg = ifelse(level <= nx, nzbeg, nzbeg + (level - nx)*nzd)
            jend = ifelse(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            !$omp parallel do private(i, j, v2, tauc, t0c, px0c, pz0c, i1, i2, j1, j2, &
                !$omp taux, t0x, signx, tauz, t0z, signz, taucx, taucz, u, d, &
                !$omp array_tauc, array_tau, array_p0, array_t0, array_sign, &
                !$omp da, db, signa, signb, dadt, dbdt, taua, taub, t0a, t0b, tauca, taucb, &
                !$omp  pa0c, pb0c, travel, dada, dbdb, root) schedule(dynamic)
            do j = jbeg, jend, nzd

                ! The condition is
                ! abs(i - nxbeg) + abs(j - nzbeg) + 1 = level
                i = ((level - 1) - abs(j - nzbeg))/nxd + nxbeg

                v2 = vp(i, j)**2
                tauc = tau(i, j)
                t0c = t0(i, j)
                px0c = px0(i, j)
                pz0c = pz0(i, j)

                if (i == 1) then

                    i1 = i
                    i2 = i + 1

                    taux = tau(i2, j)
                    t0x = t0(i2, j)
                    signx = -1.0

                else if (i == nx) then

                    i1 = i - 1
                    i2 = i

                    taux = tau(i1, j)
                    t0x = t0(i1, j)
                    signx = 1.0

                else

                    i1 = i - 1
                    i2 = i + 1

                    if (tau(i1, j)*t0(i1, j) <= tau(i2, j)*t0(i2, j)) then
                        taux = tau(i1, j)
                        t0x = t0(i1, j)
                        signx = 1.0
                    else
                        taux = tau(i2, j)
                        t0x = t0(i2, j)
                        signx = -1.0
                    end if

                end if

                if (j == 1) then

                    j1 = j
                    j2 = j + 1

                    tauz = tau(i, j2)
                    t0z = t0(i, j2)
                    signz = -1.0

                else if (j == nz) then

                    j1 = j - 1
                    j2 = j

                    tauz = tau(i, j1)
                    t0z = t0(i, j1)
                    signz = 1.0

                else

                    j1 = j - 1
                    j2 = j + 1

                    if (tau(i, j1)*t0(i, j1) <= tau(i, j2)*t0(i, j2)) then
                        tauz = tau(i, j1)
                        t0z = t0(i, j1)
                        signz = 1.0
                    else
                        tauz = tau(i, j2)
                        t0z = t0(i, j2)
                        signz = -1.0
                    end if

                end if

                if (taux == huge_value .and. tauz == huge_value) then
                    cycle
                end if

                taucx = (t0c*taux + dx/vp(i, j))/(t0c + abs(px0c)*dx)
                taucz = (t0c*tauz + dz/vp(i, j))/(t0c + abs(pz0c)*dz)

                u(1) = min(tau(i1, j)*t0(i1, j), tau(i2, j)*t0(i2, j))
                u(2) = min(tau(i, j1)*t0(i, j1), tau(i, j2)*t0(i, j2))
                d = [dx, dz]
                array_tauc = [taucx, taucz]
                array_p0 = [px0c, pz0c]
                array_tau = [taux, tauz]
                array_t0 = [t0x, t0z]
                array_sign = [signx, signz]

                if (u(1) == huge_value .and. u(2) == huge_value) then
                    cycle
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

                travel = u(1) + d(1)/vp(i, j)
                tauc = travel/t0(i, j)
                if (travel > u(2)) then

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

                    ! Solve the quadratic equation analytically
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
                end if

                tau(i, j) = min(tau(i, j), tauc)

            end do
            !$omp end parallel do

        end do

    end

#else

    !
    ! Fast sweeping local solver
    !
    function local_solver_2d(t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz, dx, dz, signx, signz, v) result(txz)

        double precision :: t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz
        double precision :: dx, dz, v
        integer :: signx, signz
        double precision :: txz

        double precision :: a1, a2, b1, b2
        double precision :: a, b, c
        double precision :: taucx, taucz
        double precision :: tau1, tau2
        logical :: causality1, causality2

        a1 = pdxt0c + t0c/dx*signx
        a2 = pdzt0c + t0c/dz*signz
        b1 = -taux*t0c/dx*signx
        b2 = -tauz*t0c/dz*signz

        a = a1**2 + a2**2
        b = 2*(a1*b1 + a2*b2)
        c = b1**2 + b2**2 - 1.0/v**2

        if (b**2 - 4*a*c > 0) then

            tau1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
            tau2 = (-b - sqrt(b**2 - 4*a*c))/(2*a)

            causality1 = tau1*t0c >= max(taux*t0x, tauz*t0z)
            causality2 = tau2*t0c >= max(taux*t0x, tauz*t0z)

            if (causality1 .and. causality2) then
                txz = min(tau1, tau2)
            else if (causality1 .and. .not. causality2) then
                txz = tau1
            else if (.not. causality1 .and. causality2) then
                txz = tau2
            else
                taucx = max((t0c*taux + dx/v)/(t0c + pdxt0c*dx*signx), taux*t0x/t0c)
                taucz = max((t0c*tauz + dz/v)/(t0c + pdzt0c*dz*signz), tauz*t0z/t0c)
                txz = min(taucx, taucz)
            end if

        else

            taucx = max((t0c*taux + dx/v)/(t0c + pdxt0c*dx*signx), taux*t0x/t0c)
            taucz = max((t0c*tauz + dz/v)/(t0c + pdzt0c*dz*signz), tauz*t0z/t0c)
            txz = min(taucx, taucz)

        end if

    end function local_solver_2d

    !
    !> Parallel fast sweeping implementing
    !> Detrixhe et al., 2013, JCP,
    !> A parallel fast sweeping method for the eikonal equation
    !> doi: 10.1016/j.jcp.2012.11.042
    !> with modifications
    !
    subroutine fast_sweep_forward(nx, nz, nxd, nzd, dx, dz, vp, t0, px0, pz0, tau)

        integer, intent(in) :: nx, nz, nxd, nzd
        double precision, intent(in) :: dx, dz
        double precision, dimension(:, :), intent(in) :: vp, t0, px0, pz0
        double precision, dimension(:, :), intent(inout) :: tau

        integer :: nxbeg, nxend, nzbeg, nzend
        integer :: i, j
        integer :: level, jbeg, jend
        double precision :: txz1, txz2, txz3, txz4
        logical :: valid1, valid2, valid3, valid4

        if (nxd > 0) then
            nxbeg = 1
            nxend = nx
        else
            nxbeg = nx
            nxend = 1
        end if

        if (nzd > 0) then
            nzbeg = 1
            nzend = nz
        else
            nzbeg = nz
            nzend = 1
        end if

        do level = 1, nx + nz - 1

            jbeg = ifelse(level <= nx, nzbeg, nzbeg + (level - nx)*nzd)
            jend = ifelse(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            !$omp parallel do private(i, j, txz1, txz2, txz3, txz4, valid1, valid2, valid3, valid4) schedule(dynamic)
            do j = jbeg, jend, nzd

                ! The condition is
                ! abs(i - nxbeg) + abs(j - nzbeg) + 1 = level
                i = ((level - 1) - abs(j - nzbeg))/nxd + nxbeg

                txz1 = huge_value
                txz2 = huge_value
                txz3 = huge_value
                txz4 = huge_value

                valid1 = .true.
                valid2 = .true.
                valid3 = .true.
                valid4 = .true.

                if (i == 1) then
                    valid1 = .false.
                    valid4 = .false.
                else if (i == nx) then
                    valid2 = .false.
                    valid3 = .false.
                end if
                if (j == 1) then
                    valid1 = .false.
                    valid2 = .false.
                else if (j == nz) then
                    valid3 = .false.
                    valid4 = .false.
                end if

                if (valid1) then
                    txz1 = local_solver_2d(t0(i, j), t0(i - 1, j), t0(i, j - 1), &
                        px0(i, j), pz0(i, j), tau(i - 1, j), tau(i, j - 1), dx, dz, +1, +1, vp(i, j))
                end if

                if (valid2) then
                    txz2 = local_solver_2d(t0(i, j), t0(i + 1, j), t0(i, j - 1), &
                        px0(i, j), pz0(i, j), tau(i + 1, j), tau(i, j - 1), dx, dz, -1, +1, vp(i, j))
                end if

                if (valid3) then
                    txz3 = local_solver_2d(t0(i, j), t0(i + 1, j), t0(i, j + 1), &
                        px0(i, j), pz0(i, j), tau(i + 1, j), tau(i, j + 1), dx, dz, -1, -1, vp(i, j))
                end if

                if (valid4) then
                    txz4 = local_solver_2d(t0(i, j), t0(i - 1, j), t0(i, j + 1), &
                        px0(i, j), pz0(i, j), tau(i - 1, j), tau(i, j + 1), dx, dz, +1, -1, vp(i, j))
                end if

                tau(i, j) = min(tau(i, j), txz1, txz2, txz3, txz4)

            end do
            !$omp end parallel do

        end do

    end subroutine fast_sweep_forward

#endif

    !
    !> Fast sweeping for adjoint state equation
    !
    subroutine fast_sweep_adjoint(n1, n2, n1d, n2d, d1, d2, lambda)

        integer, intent(in) :: n1, n2, n1d, n2d
        double precision, intent(in) :: d1, d2
        double precision, dimension(:, :), intent(inout) :: lambda

        integer :: i, j
        double precision :: app, amp, apm, amm
        double precision :: bpp, bmp, bpm, bmm
        double precision :: ap, am, bp, bm
        double precision :: lhs, rhs, t
        integer :: n1beg, n1end, n2beg, n2end
        integer :: i1, i2, j1, j2
        integer :: level, jbeg, jend

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

        do level = 1, n1 + n2 - 1

            jbeg = ifelse(level <= n1, n2beg, n2beg + (level - n1)*n2d)
            jend = ifelse(level <= n2, n2beg + (level - 1)*n2d, n2end)

            !$omp parallel do private(i, j, i1, i2, j1, j2, &
                !$omp app, amp, apm, amm, &
                !$omp bpp, bmp, bpm, bmm, &
                !$omp ap, am, bp, bm, &
                !$omp lhs, rhs, t)
            do j = jbeg, jend, n2d

                ! The condition is
                ! abs(i - n1beg) + abs(j - n2beg) + 1 = level
                i = ((level - 1) - abs(j - n2beg))/n1d + n1beg

                if (.not. recrflag(i, j)) then

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

                    ! Solve equation (A-9) in Taillandier et al. (2009)
                    ap = (tt(i2, j) - tt(i, j))/d1
                    am = (tt(i, j) - tt(i1, j))/d1

                    bp = (tt(i, j2) - tt(i, j))/d2
                    bm = (tt(i, j) - tt(i, j1))/d2

                    app = (ap + abs(ap))/2.0
                    apm = (ap - abs(ap))/2.0

                    amp = (am + abs(am))/2.0
                    amm = (am - abs(am))/2.0

                    bpp = (bp + abs(bp))/2.0
                    bpm = (bp - abs(bp))/2.0

                    bmp = (bm + abs(bm))/2.0
                    bmm = (bm - abs(bm))/2.0

                    lhs = (apm - amp)/d1 + (bpm - bmp)/d2
                    rhs = (amm*lambda(i1, j) - app*lambda(i2, j))/d1 &
                        + (bmm*lambda(i, j1) - bpp*lambda(i, j2))/d2

                    if (lhs == 0) then
                        t = 0
                    else
                        t = rhs/lhs
                    end if

                    lambda(i, j) = min(lambda(i, j), t)

                end if

            end do
        end do

    end subroutine fast_sweep_adjoint

    !
    !> Background time field of a source of one or many points: at every
    !> grid point the analytic time from the source point that arrives
    !> first, t0 = min_l (|x - s_l|/v_l + t0_l), and its derivatives, into
    !> the module arrays t0, pdxt0, pdzt0.
    !
    !> The minimum used to re-interpolate every source point's velocity at
    !> every grid point and fill three ns-long arrays per grid point:
    !> nx*nz*ns interpolations, with the source points taken one after
    !> another. A reflector ensemble (forward_iso_reflection) has one for
    !> every reflector node. Now each source point's velocity is
    !> interpolated once, and an ensemble larger than background_direct is
    !> binned in space: a bin whose lower bound -- distance to its box over
    !> its fastest velocity, plus its earliest start time -- is above the
    !> best time found so far cannot hold the minimum, and is skipped. The
    !> pruning is exact: the same source point wins, by minloc's rule for
    !> ties (the lowest index), and its time and derivatives are the same
    !> expressions, so the field is identical to evaluating every point.
    !> The loop runs over grid rows in parallel; along a row, each grid
    !> point starts from its neighbour's winner, which is usually its own.
    !
    subroutine background_time(geom, source_inside, nx, nz, ox, oz, dx, dz, vp)

        type(source_receiver_geometry), intent(in) :: geom
        logical, dimension(:), intent(in) :: source_inside
        integer, intent(in) :: nx, nz
        double precision, intent(in) :: ox, oz, dx, dz
        double precision, dimension(:, :), intent(in) :: vp

        ! Up to this many source points, all are evaluated at every grid point
        integer, parameter :: background_direct = 16
        double precision, allocatable, dimension(:) :: sx, sz, sv, st, gx, gz
        double precision, allocatable, dimension(:) :: bxl, bxh, bzl, bzh, bvmax, btmin
        integer, allocatable, dimension(:) :: member, bin_of, bin_count, bin_start, order, bfirst, blast
        integer :: ns, nb, nbin, nbx, nbz, npos, i, j, l, m, b, h, win, prev, first
        double precision :: ex, ez, cell, best, tval, lb, lbmin, tol

        t0 = zeros(nx, nz) + huge_value
        pdxt0 = zeros(nx, nz)
        pdzt0 = zeros(nx, nz)

        ! The source points that take part, in index order, with what the
        ! minimum needs of each
        ns = count(source_inside)
        if (ns == 0) then
            return
        end if
        member = pack([(l, l = 1, geom%ns)], source_inside)
        allocate (sx(1:ns), sz(1:ns), sv(1:ns), st(1:ns))
        !$omp parallel do private(m, l)
        do m = 1, ns
            l = member(m)
            sx(m) = geom%srcr(l)%x
            sz(m) = geom%srcr(l)%z
            ! Interpolate to get the velocity corresponding to the source point
            ! inside a rectange or triangle
            sv(m) = get_point_value_inside(dble([geom%srcr(l)%x, geom%srcr(l)%z]), [nx, nz], [ox, oz], [dx, dz], vp)
            st(m) = geom%srcr(l)%t0
            !                ! In comparison, the nearest grid point approach is
            !                isx = nint((geom%srcr(l)%x - ox)/dx) + 1
            !                isz = nint((geom%srcr(l)%z - oz)/dz) + 1
            !                vsource = vp(isx, isz)
            !                dsx = (i - isx)*dx
            !                dsz = (j - isz)*dz
        end do
        !$omp end parallel do

        allocate (gx(1:nx), gz(1:nz))
        do i = 1, nx
            gx(i) = ox + (i - 1)*dx
        end do
        do j = 1, nz
            gz(j) = oz + (j - 1)*dz
        end do

        if (ns <= background_direct) then

            !$omp parallel do private(i, j, m, best, win, tval)
            do j = 1, nz
                do i = 1, nx
                    best = huge_value
                    win = 0
                    do m = 1, ns
                        tval = source_time(m, i, j)
                        if (tval < best) then
                            best = tval
                            win = m
                        end if
                    end do
                    call set_background(win, i, j)
                end do
            end do
            !$omp end parallel do
            return

        end if

        ! Bins: a uniform grid over the box of the source points, about
        ! 2*sqrt(ns) cells -- per grid point, one bound per non-empty bin
        ! plus the points of the few bins the bound does not rule out
        ex = maxval(sx) - minval(sx)
        ez = maxval(sz) - minval(sz)
        npos = count([ex, ez] > 0)
        nbx = 1
        nbz = 1
        if (npos > 0) then
            cell = (product(pack([ex, ez], [ex, ez] > 0))/(2.0d0*sqrt(dble(ns))))**(1.0d0/npos)
            if (ex > 0) nbx = max(1, min(ns, ceiling(ex/cell)))
            if (ez > 0) nbz = max(1, min(ns, ceiling(ez/cell)))
        end if
        nbin = nbx*nbz
        allocate (bin_of(1:ns))
        do m = 1, ns
            bin_of(m) = bin_index(sx(m), minval(sx), ex, nbx) &
                + nbx*(bin_index(sz(m), minval(sz), ez, nbz) - 1)
        end do

        ! A counting sort by bin, stable, so each bin lists its points in
        ! index order
        allocate (bin_count(1:nbin), bin_start(1:nbin + 1))
        bin_count = 0
        do m = 1, ns
            bin_count(bin_of(m)) = bin_count(bin_of(m)) + 1
        end do
        bin_start(1) = 1
        do b = 1, nbin
            bin_start(b + 1) = bin_start(b) + bin_count(b)
        end do
        allocate (order(1:ns))
        bin_count = 0
        do m = 1, ns
            b = bin_of(m)
            order(bin_start(b) + bin_count(b)) = m
            bin_count(b) = bin_count(b) + 1
        end do

        ! The non-empty bins: their members, tight box, fastest velocity
        ! and earliest start time
        nb = count(bin_count > 0)
        allocate (bfirst(1:nb), blast(1:nb), bxl(1:nb), bxh(1:nb), bzl(1:nb), bzh(1:nb), &
            bvmax(1:nb), btmin(1:nb))
        h = 0
        do b = 1, nbin
            if (bin_count(b) == 0) then
                cycle
            end if
            h = h + 1
            bfirst(h) = bin_start(b)
            blast(h) = bin_start(b + 1) - 1
            associate (mm => order(bin_start(b):bin_start(b + 1) - 1))
                bxl(h) = minval(sx(mm))
                bxh(h) = maxval(sx(mm))
                bzl(h) = minval(sz(mm))
                bzh(h) = maxval(sz(mm))
                bvmax(h) = maxval(sv(mm))
                btmin(h) = minval(st(mm))
            end associate
        end do

        !$omp parallel do schedule(dynamic) &
            !$omp private(i, j, m, b, h, best, win, prev, first, tval, lb, lbmin, tol)
        do j = 1, nz
            prev = 0
            do i = 1, nx

                best = huge_value
                win = 0
                if (prev > 0) then
                    ! the neighbour's winner, which is usually this
                    ! point's as well: a tight bound from the start
                    best = source_time(prev, i, j)
                    win = prev
                else
                    ! otherwise the bin with the lowest bound goes first
                    lbmin = huge_value
                    first = 1
                    do b = 1, nb
                        lb = bin_bound(b, i, j)
                        if (lb < lbmin) then
                            lbmin = lb
                            first = b
                        end if
                    end do
                    do h = bfirst(first), blast(first)
                        m = order(h)
                        tval = source_time(m, i, j)
                        if (tval < best .or. (tval == best .and. m < win)) then
                            best = tval
                            win = m
                        end if
                    end do
                end if

                ! Every bin that could hold a time at or below the best so
                ! far; the margin keeps rounding in the bound from ever
                ! ruling out the true minimum
                tol = 1.0d-12*abs(best)
                do b = 1, nb
                    if (bin_bound(b, i, j) > best + tol) then
                        cycle
                    end if
                    do h = bfirst(b), blast(b)
                        m = order(h)
                        tval = source_time(m, i, j)
                        if (tval < best .or. (tval == best .and. m < win)) then
                            best = tval
                            win = m
                            tol = 1.0d-12*abs(best)
                        end if
                    end do
                end do

                call set_background(win, i, j)
                prev = win

            end do
        end do
        !$omp end parallel do

    contains

        ! The loop indices come in as arguments: inside a parallel loop a
        ! contained procedure would see the original i, j through host
        ! association, not the thread's private copies

        !> The time from source point m to grid point (i, j), as the
        !> original loop computed it
        double precision function source_time(m, i, j)
            integer, intent(in) :: m, i, j
            double precision :: dsx, dsz
            dsx = gx(i) - sx(m)
            dsz = gz(j) - sz(m)
            source_time = sqrt(dsx**2 + dsz**2)/sv(m)
            source_time = source_time + st(m)
        end function source_time

        !> A lower bound of the time from any source point of bin b
        double precision function bin_bound(b, i, j)
            integer, intent(in) :: b, i, j
            double precision :: cx, cz
            cx = max(bxl(b) - gx(i), 0.0d0, gx(i) - bxh(b))
            cz = max(bzl(b) - gz(j), 0.0d0, gz(j) - bzh(b))
            bin_bound = sqrt(cx**2 + cz**2)/bvmax(b) + btmin(b)
        end function bin_bound

        !> The winner's time and derivatives at (i, j), as the original
        !> loop computed them
        subroutine set_background(m, i, j)
            integer, intent(in) :: m, i, j
            double precision :: dsx, dsz, t1
            if (m == 0) then
                return
            end if
            dsx = gx(i) - sx(m)
            dsz = gz(j) - sz(m)
            t1 = sqrt(dsx**2 + dsz**2)/sv(m)
            if (t1 == 0) then
                pdxt0(i, j) = 0
                pdzt0(i, j) = 0
            else
                pdxt0(i, j) = dsx/sv(m)**2/t1
                pdzt0(i, j) = dsz/sv(m)**2/t1
            end if
            t0(i, j) = t1 + st(m)
        end subroutine set_background

    end subroutine background_time

    !
    !> Which of n bins a coordinate falls in, over [low, low + extent]
    !
    pure integer function bin_index(x, low, extent, n)
        double precision, intent(in) :: x, low, extent
        integer, intent(in) :: n
        if (extent > 0 .and. n > 1) then
            bin_index = min(n, int((x - low)/extent*n) + 1)
        else
            bin_index = 1
        end if
    end function bin_index

    !
    !> A point's grid cell as one integer: the floor and ceiling indices of
    !> each coordinate, exactly what within_the_same_grid compares
    !
    function cell_key(p, n, o, d) result(key)
        double precision, dimension(:), intent(in) :: p, o, d
        integer, dimension(:), intent(in) :: n
        integer(kind=8) :: key
        integer(kind=8) :: c1, c2
        double precision :: p1, p2
        p1 = p(1) - o(1)
        p2 = p(2) - o(2)
        c1 = 2*int(floor(p1/d(1)), 8) + (ceiling(p1/d(1)) - floor(p1/d(1))) + 4
        c2 = 2*int(floor(p2/d(2)), 8) + (ceiling(p2/d(2)) - floor(p2/d(2))) + 4
        key = c1 + (2*int(n(1), 8) + 8)*c2
    end function cell_key

    !
    !> Fast sweeping factorized eikonal solver
    !
    subroutine forward_iso(v, d, o, geom, t, trec)

        real, dimension(:, :), intent(in) :: v
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        double precision, allocatable, dimension(:, :) :: tt_prev, vp
        double precision :: dx, dz, ox, oz
        integer :: nx, nz, niter, i, l, isx, isz, irx, irz, itx, itz
        double precision :: ttdiff, vsource, dsx, dsz
        integer :: sw
        integer(kind=8) :: key
        integer(kind=8), allocatable, dimension(:) :: cellkey
        logical, allocatable, dimension(:) :: source_inside, receiver_inside
        double precision :: time1, time2

        dx = d(1)
        dz = d(2)
        ox = o(1)
        oz = o(2)

        vp = transpose(v)
        nx = size(vp, 1)
        nz = size(vp, 2)

        ! Check if sources are at zero-velocity points
        source_inside = falses(geom%ns)
        !$omp parallel do private(l)
        do l = 1, geom%ns
            if (geom%srcr(l)%amp == 1) then
                source_inside(l) = point_in_domain(dble([geom%srcr(l)%x, geom%srcr(l)%z]), [nx, nz], [ox, oz], [dx, dz], vp)
            end if
        end do
        !$omp end parallel do

        ! Check if receivers are at zero-velocity points
        receiver_inside = falses(geom%nr)
        !$omp parallel do private(i)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                receiver_inside(i) = point_in_domain(dble([geom%recr(i)%x, geom%recr(i)%z]), [nx, nz], [ox, oz], [dx, dz], vp)
            end if
        end do
        !$omp end parallel do

        ! Background time field: the analytic time from the source point
        ! that arrives first, and its derivatives
        call background_time(geom, source_inside, nx, nz, ox, oz, dx, dz, vp)

        ! Initialize the multiplicative time field
        tt_prev = zeros(nx, nz)
        tt = zeros(nx, nz) + huge_value
        ! Increase the width of the multiplicative field around the source
        ! in the initialization seems to significantly improve the resulting accuracy,
        ! especially when the source point is not on an integer grid point
        sw = 3
        !$omp parallel do private(i, isx, isz, itx, itz)
        do i = 1, geom%ns
            if (source_inside(i)) then
                isx = max(floor((geom%srcr(i)%x - ox)/dx) + 1 - sw, 1)
                isz = max(floor((geom%srcr(i)%z - oz)/dz) + 1 - sw, 1)
                itx = min(ceiling((geom%srcr(i)%x - ox)/dx) + 1 + sw, nx)
                itz = min(ceiling((geom%srcr(i)%z - oz)/dz) + 1 + sw, nz)
                tt(isx:itx, isz:itz) = 1.0
            end if
        end do
        !$omp end parallel do

        !        ! In comparison, the nearest grid point approach is
        !        !$omp parallel do private(i, isx, isz, itx, itz)
        !        do i = 1, geom%ns
        !            if (source_inside(i)) then
        !                isx = nint((geom%srcr(i)%x - ox)/dx) + 1
        !                isz = nint((geom%srcr(i)%z - oz)/dz) + 1
        !                tt(isx, isz) = 1.0
        !            end if
        !        end do
        !        !$omp end parallel do

        ! Fast sweeping via Gauss-Seidel iterations
        niter = 0
        ttdiff = huge_value
        where (vp == 0)
            vp = float_tiny
        end where

        ! For timing the modeling, commented out for clarity
        time1 = omp_get_wtime()

        do while (ttdiff >= sweep_stop_threshold .and. niter < sweep_niter_max)

            tt_prev = tt
            call fast_sweep_forward(nx, nz, 1, 1, dx, dz, vp, t0, pdxt0, pdzt0, tt)
            call fast_sweep_forward(nx, nz, 1, -1, dx, dz, vp, t0, pdxt0, pdzt0, tt)
            call fast_sweep_forward(nx, nz, -1, 1, dx, dz, vp, t0, pdxt0, pdzt0, tt)
            call fast_sweep_forward(nx, nz, -1, -1, dx, dz, vp, t0, pdxt0, pdzt0, tt)

            ttdiff = mean(abs(tt - tt_prev))
            niter = niter + 1

        end do

        ! For timing the modeling, commented out for clarity
        time2 = omp_get_wtime()
        call warn(' Wall-clock time = '//num2str(time2 - time1, '(es)'))

        tt = t0*tt
        where (vp == float_tiny)
            tt = 0
            vp = 0
        end where

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        ! A receiver in the same grid cell as a source point takes the
        ! analytic time from it -- from the last such source point, as the
        ! pairwise test within_the_same_grid always did. That test ran for
        ! every receiver against every source point; each source point's
        ! cell is now one integer (cell_key), computed once.
        allocate (cellkey(1:geom%ns))
        cellkey = -1
        !$omp parallel do private(l)
        do l = 1, geom%ns
            if (source_inside(l)) then
                cellkey(l) = cell_key(dble([geom%srcr(l)%x, geom%srcr(l)%z]), [nx, nz], [ox, oz], [dx, dz])
            end if
        end do
        !$omp end parallel do
        !$omp parallel do private(i, irx, irz, l, key, vsource, dsx, dsz)
        do i = 1, geom%nr
            if (receiver_inside(i)) then

                ! Linear interpolation; same as for the source points,
                ! here the implementation automatically handles the case of
                ! a receiver falling on integer grid points.
                trec(i, 1) = get_point_value_inside(dble([geom%recr(i)%x, geom%recr(i)%z]), [nx, nz], [ox, oz], [dx, dz], tt)
                key = cell_key(dble([geom%recr(i)%x, geom%recr(i)%z]), [nx, nz], [ox, oz], [dx, dz])
                do l = geom%ns, 1, -1
                    if (cellkey(l) == key) then
                        vsource = get_point_value_inside(dble([geom%srcr(l)%x, geom%srcr(l)%z]), [nx, nz], [ox, oz], [dx, dz], vp)
                        dsx = geom%recr(i)%x - geom%srcr(l)%x
                        dsz = geom%recr(i)%z - geom%srcr(l)%z
                        trec(i, 1) = sqrt(dsx**2 + dsz**2)/vsource + geom%srcr(l)%t0
                        exit
                    end if
                end do

                ! Add virtual receiver time for TLOC
                trec(i, 1) = trec(i, 1) + geom%recr(i)%t0

                !                ! In comparison, the nearest grid point approach is
                !                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                !                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                !                trec(i, 1) = t(irx, irz) + geom%recr(i)%t0

            end if
        end do
        !$omp end parallel do

        t = transpose(tt)

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(niter)// &
            ', relative diff = '//num2str(ttdiff, '(es)'))

    end subroutine forward_iso

    !
    !> Compute adjoint field for FATT
    !
    subroutine adjoint_iso(v, d, o, geom, t, tresidual, tadj)

        real, dimension(:, :), intent(in) :: v, t
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :), allocatable, intent(out) :: tadj

        double precision, allocatable, dimension(:, :) :: vp, lambda_prev
        integer :: nx, nz, iter, i, irx, irz !, j
        double precision :: lambda_diff, dx, dz, ox, oz

        dx = d(1)
        dz = d(2)
        ox = o(1)
        oz = o(2)

        tt = transpose(t)
        vp = transpose(v)
        nx = size(vp, 1)
        nz = size(vp, 2)

        recrflag = falses(nx, nz)
        lambda = zeros(nx, nz)
        lambda_prev = zeros(nx, nz)
        lambda_diff = huge_value

        ! Get the misfit values at the receivers; here we cannot use parallel as
        ! in loc + tomo, initial virtual receivers may be at the same position.
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                lambda(irx, irz) = lambda(irx, irz) + tresidual(i, 1)
                recrflag(irx, irz) = .true.
            end if
        end do
        where (.not. recrflag)
            lambda = huge_value
        end where

        !        ! This is to set lambda = 0 on boundaries, but does not make a difference
        !        ! when both source and receivers are inside and could be ignored.
        !        ! But when this is set, the preconditioner input residual must be negative to
        !        ! avoid unwanted boundary artifacts.
        !        !$omp parallel do private(i, j)
        !        do j = 1, nz
        !            do i = 1, nx
        !                if ((i == 1 .or. i == nx .or. j == 1 .or. j == nz) .and. (.not.recrflag(i, j))) then
        !                    lambda(i, j) = 0.0
        !                end if
        !            end do
        !        end do
        !        !$omp end parallel do

        ! Fast sweep to solve the adjoint-state equation
        iter = 0
        do while (lambda_diff >= sweep_stop_threshold .and. iter < sweep_niter_max)

            ! Save previous step
            lambda_prev = lambda

            ! Fast sweeps
            call fast_sweep_adjoint(nx, nz, 1, 1, dx, dz, lambda)
            call fast_sweep_adjoint(nx, nz, 1, -1, dx, dz, lambda)
            call fast_sweep_adjoint(nx, nz, -1, 1, dx, dz, lambda)
            call fast_sweep_adjoint(nx, nz, -1, -1, dx, dz, lambda)

            ! Check threshold
            lambda_diff = mean(abs(lambda - lambda_prev))

            iter = iter + 1

        end do

        tadj = transpose(lambda)

        call warn(date_time_compact()//' Fast sweeping adjoint niter = '//num2str(iter)// &
            ', relative diff = '//num2str(lambda_diff, '(es)'))

    end subroutine adjoint_iso

end module traveltime_iso
