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

    real, allocatable, dimension(:, :) :: t0, pdxt0, pdzt0, tt, lambda
    logical, allocatable, dimension(:, :) :: unknown, recrflag

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
    subroutine solve_quadratic(i, j, dx, dz, vx, vz, f)

        integer, intent(in) :: i, j
        real, intent(in) :: dx, dz, vx, vz, f

        double precision :: timeleft1, timeright1
        double precision :: timetop1, timebottom1
        double precision :: t0ij, taux, tauz, tau
        double precision :: a, b, c
        double precision :: ax, bx, cx
        double precision :: az, bz, cz
        double precision :: sx, sz, s
        integer :: dirx, dirz
        double precision :: alpha, beta

        t0ij = t0(i, j)
        timeleft1 = t0(i - 1, j)*tt(i - 1, j)
        timeright1 = t0(i + 1, j)*tt(i + 1, j)
        timetop1 = t0(i, j - 1)*tt(i, j - 1)
        timebottom1 = t0(i, j + 1)*tt(i, j + 1)

        taux = float_huge
        tauz = float_huge

        dirx = 0
        if (.not. unknown(i - 1, j)) then
            taux = timeleft1
            dirx = -1
        end if
        if (.not. unknown(i + 1, j) .and. taux > timeright1) then
            taux = timeright1
            dirx = 1
        end if

        dirz = 0
        if (.not. unknown(i, j - 1)) then
            tauz = timetop1
            dirz = -1
        end if
        if (.not. unknown(i, j + 1) .and. tauz > timebottom1) then
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
            tau = t0(i + 2*dirx, j)*tt(i + 2*dirx, j)
            alpha = 0.0
            beta = 0.0
            if (taux > tau .and. .not. unknown(i + 2*dirx, j)) then
                alpha = pdxt0(i, j) - 1.5d0*t0ij*dirx/dx
                beta = -t0ij*dirx*(tt(i + 2*dirx, j) - 4*tt(i + dirx, j))/(2*dx)
            else
                alpha = pdxt0(i, j) - t0ij*dirx/dx
                beta = t0ij*dirx*tt(i + dirx, j)/dx
            end if
            ax = vx*alpha**2
            bx = vx*2*alpha*beta
            cx = vx*beta**2
            a = a + ax
            b = b + bx
            c = c + cx
        end if

        az = 0.0
        bz = 0.0
        cz = 0.0
        if (dirz /= 0) then
            tau = t0(i, j + 2*dirz)*tt(i, j + 2*dirz)
            alpha = 0.0
            beta = 0.0
            if (tauz > tau .and. .not. unknown(i, j + 2*dirz)) then
                alpha = pdzt0(i, j) - 1.5d0*t0ij*dirz/dz
                beta = -t0ij*dirz*(tt(i, j + 2*dirz) - 4*tt(i, j + dirz))/(2*dz)
            else
                alpha = pdzt0(i, j) - t0ij*dirz/dz
                beta = t0ij*dirz*tt(i, j + dirz)/dz
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
        cz = cz - f

        sx = (-bx + sqrt(bx**2 - 4*ax*cx))/(2*ax)
        sz = (-bz + sqrt(bz**2 - 4*az*cz))/(2*az)
        s = (-b + sqrt(b**2 - 4*a*c))/(2*a)
        if (isnan(sx) .or. sx <= 0) then
            sx = float_huge
        end if
        if (isnan(sz) .or. sz <= 0) then
            sz = float_huge
        end if
        if (isnan(s) .or. s <= 0) then
            s = float_huge
        end if

        tt(i, j) = min(sx, sz, s)
        unknown(i, j) = .false.

    end subroutine solve_quadratic

    !
    !> Fast marching eikonal solver
    !
    subroutine forward_iso_fast_march(v, d, o, geom, t, trec)

        real, dimension(:, :), intent(in) :: v
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        integer :: nx, nz, l, isf
        integer :: i, j, hi, hj, ii, jj
        integer :: isx, isz, irx, irz
        type(heap_float) :: wavefront
        real :: node(1:3), v2, dx, dz, ox, oz
        real, allocatable, dimension(:, :) :: vp
        real :: dsx, dsz, vsource

        dx = d(1)
        dz = d(2)
        ox = o(1)
        oz = o(2)

        vp = transpose(v)
        nx = size(vp, 1)
        nz = size(vp, 2)

        ! Source location
        do i = 1, geom%ns
            if (geom%srcr(i)%amp == 1) then
                isf = i
                exit
            end if
        end do
        isx = nint((geom%srcr(isf)%x - ox)/dx) + 1
        isz = nint((geom%srcr(isf)%z - oz)/dz) + 1

        ! Background time field
        vsource = vp(isx, isz)
        call alloc_array(t0, [1, nx, 1, nz], pad=2)
        call alloc_array(pdxt0, [1, nx, 1, nz], pad=2)
        call alloc_array(pdzt0, [1, nx, 1, nz], pad=2)
        do j = 1 - 2, nz + 2
            do i = 1 - 2, nx + 2
                dsx = (i - isx)*dx
                dsz = (j - isz)*dz
                t0(i, j) = sqrt(dsx**2 + dsz**2)/vsource
                pdxt0(i, j) = dsx/vsource**2/t0(i, j)
                pdzt0(i, j) = dsz/vsource**2/t0(i, j)
            end do
        end do

        ! Multiplicative time field
        call alloc_array(unknown, [1, nx, 1, nz], pad=2)
        call alloc_array(tt, [1, nx, 1, nz], pad=2)
        tt = sqrt(float_huge)
        unknown = .true.
        tt(isx, isz) = 1.0
        unknown(isx, isz) = .false.
        i = isx
        j = isz

        ! Initialize heap for the wavefront
        call wavefront%init(nx*nz, 3, smaller_traveltime)

        l = 1
        ! The source point is known, therefore l < nx*nz suffices
        do while (l < nx*nz)

            ! Compute traveltime at neighbour points
            do hj = -1, 1
                do hi = -1, 1
                    ii = clip(i + hi, 1, nx)
                    jj = clip(j + hj, 1, nz)
                    if (unknown(ii, jj) .and. abs(hi) /= abs(hj)) then
                        v2 = vp(ii, jj)**2
                        call solve_quadratic(ii, jj, dx, dz, v2, v2, 1.0)
                        call wavefront%insert([t0(ii, jj)*tt(ii, jj), real(ii), real(jj)])
                    end if
                end do
            end do

            ! Remove the smallest traveltime
            call wavefront%pop(node)
            i = int(node(2))
            j = int(node(3))

            l = l + 1

        end do

        call wavefront%free

        ! Select and resample to the original domain, i.e, [ox, oz] --> [0, 0] origin
        t = t0(1:nx, 1:nz)*tt(1:nx, 1:nz) + geom%srcr(isf)%t0

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                ! No interpolation here; just nearest grid point; add virtual receiver time for TLOC
                trec(i, 1) = t(irx, irz) + geom%recr(i)%t0
            end if
        end do

        ! Transpose to output
        t = transpose(t)

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
    !> with modifications
    !
    subroutine fast_sweep_forward(nx, nz, nxd, nzd, dx, dz, vp, t0, px0, pz0, tau)

        integer, intent(in) :: nx, nz, nxd, nzd
        real, intent(in) :: dx, dz
        real, dimension(:, :), intent(in) :: vp, t0, px0, pz0
        real, dimension(:, :), intent(inout) :: tau

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

            jbeg = ifthen(level <= nx, nzbeg, nzbeg + (level - nx)*nzd)
            jend = ifthen(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            !$omp parallel do private(i, j, v2, tauc, t0c, px0c, pz0c, i1, i2, j1, j2, &
                !$omp taux, t0x, signx, tauz, t0z, signz, taucx, taucz, u, d, &
                !$omp array_tauc, array_tau, array_p0, array_t0, array_sign, &
                !$omp da, db, signa, signb, dadt, dbdt, taua, taub, t0a, t0b, tauca, taucb, &
                !$omp  pa0c, pb0c, travel, dada, dbdb, root)
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

                if (taux == float_huge .and. tauz == float_huge) then
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

                if (u(1) == float_huge .and. u(2) == float_huge) then
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
                end if

                tau(i, j) = min(tau(i, j), tauc)

            end do
            !$omp end parallel do

        end do

    end subroutine fast_sweep_forward

    !
    !> Fast sweeping for adjoint state equation
    !
    subroutine fast_sweep_adjoint(n1, n2, n1d, n2d, d1, d2, lambda)

        integer, intent(in) :: n1, n2, n1d, n2d
        real, intent(in) :: d1, d2
        real, dimension(:, :), intent(inout) :: lambda

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

            jbeg = ifthen(level <= n1, n2beg, n2beg + (level - n1)*n2d)
            jend = ifthen(level <= n2, n2beg + (level - 1)*n2d, n2end)

            !$omp parallel do private(i, j, &
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

#else

    !
    !> Serial fast sweeping
    !
    subroutine fast_sweep_forward(nx, nz, nxd, nzd, dx, dz, vp, t0, px0, pz0, tau)

        integer, intent(in) :: nx, nz, nxd, nzd
        real, intent(in) :: dx, dz
        real, dimension(:, :), intent(in) :: vp, t0, px0, pz0
        real, dimension(:, :), intent(inout) :: tau

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

        do j = nzbeg, nzend, nzd
            do i = nxbeg, nxend, nxd

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

                if (taux == float_huge .and. tauz == float_huge) then
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

                if (u(1) == float_huge .and. u(2) == float_huge) then
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
                end if

                tau(i, j) = min(tau(i, j), tauc)

            end do
        end do

    end subroutine fast_sweep_forward

    !
    !> Fast sweeping for adjoint state equation
    !
    subroutine fast_sweep_adjoint(n1, n2, n1d, n2d, d1, d2, lambda)

        integer, intent(in) :: n1, n2, n1d, n2d
        real, intent(in) :: d1, d2
        real, dimension(:, :), intent(inout) :: lambda

        integer :: i, j
        double precision :: app, amp, apm, amm
        double precision :: bpp, bmp, bpm, bmm
        double precision :: ap, am, bp, bm
        double precision :: lhs, rhs, t
        integer :: n1beg, n1end, n2beg, n2end
        integer :: i1, i2, j1, j2

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

        ! Sweep over finite-difference grids
        do j = n2beg, n2end, n2d
            do i = n1beg, n1end, n1d

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

#endif

    !
    !> Fast sweeping factorized eikonal solver
    !
    subroutine forward_iso_fast_sweep(v, d, o, geom, t, trec)

        real, dimension(:, :), intent(in) :: v
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        real, allocatable, dimension(:, :) :: tt_prev, vp
        real :: dx, dz, ox, oz
        integer :: nx, nz, niter, i, j, l, isx, isz, irx, irz
        real, allocatable, dimension(:) :: t1, t2, t3
        real :: ttdiff, vsource, dsx, dsz
        integer :: imin

        dx = d(1)
        dz = d(2)
        ox = o(1)
        oz = o(2)

        vp = transpose(v)
        nx = size(vp, 1)
        nz = size(vp, 2)

        t0 = zeros(nx, nz)
        pdxt0 = zeros(nx, nz)
        pdzt0 = zeros(nx, nz)

        ! Source location
        t1 = zeros(geom%ns) + sqrt(float_huge)
        t2 = zeros_like(t1)
        t3 = zeros_like(t1)

        ! Background time field
        !$omp parallel do private(i, j, l, isx, isz, vsource, dsx, dsz, t1, t2, t3, imin)
        do j = 1, nz
            do i = 1, nx

                do l = 1, geom%ns

                    if (geom%srcr(l)%amp == 1) then

                        isx = nint((geom%srcr(l)%x - ox)/dx) + 1
                        isz = nint((geom%srcr(l)%z - oz)/dz) + 1
                        vsource = vp(isx, isz)
                        dsx = (i - isx)*dx
                        dsz = (j - isz)*dz

                        t1(l) = sqrt(dsx**2 + dsz**2)/vsource
                        if (t1(l) == 0) then
                            t2(l) = 0
                            t3(l) = 0
                        else
                            t2(l) = dsx/vsource**2/t1(l)
                            t3(l) = dsz/vsource**2/t1(l)
                        end if

                        t1(l) = t1(l) + geom%srcr(l)%t0

                    end if

                end do

                imin = as_scalar(minloc(t1))
                t0(i, j) = t1(imin)
                pdxt0(i, j) = t2(imin)
                pdzt0(i, j) = t3(imin)

            end do
        end do
        !$omp end parallel do

        tt_prev = zeros(nx, nz)
        tt = zeros(nx, nz) + sqrt(float_huge)
        !$omp parallel do private(i, isx, isz)
        do i = 1, geom%ns
            if (geom%srcr(i)%amp == 1) then
                isx = nint((geom%srcr(i)%x - ox)/dx) + 1
                isz = nint((geom%srcr(i)%z - oz)/dz) + 1
                tt(isx, isz) = 1.0
            end if
        end do
        !$omp end parallel do

        ! Fast sweeping via Gauss-Seidel iterations
        niter = 0
        ttdiff = sqrt(float_huge)
        do while (ttdiff >= sweep_stop_threshold .and. niter < sweep_niter_max)

            tt_prev = tt
            call fast_sweep_forward(nx, nz, 1, 1, dx, dz, vp, t0, pdxt0, pdzt0, tt)
            call fast_sweep_forward(nx, nz, 1, -1, dx, dz, vp, t0, pdxt0, pdzt0, tt)
            call fast_sweep_forward(nx, nz, -1, 1, dx, dz, vp, t0, pdxt0, pdzt0, tt)
            call fast_sweep_forward(nx, nz, -1, -1, dx, dz, vp, t0, pdxt0, pdzt0, tt)

            ttdiff = mean(abs(tt - tt_prev))
            niter = niter + 1

        end do

        t = t0*tt

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        !$omp parallel do private(i, irx, irz)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                ! No interpolation here; just nearest grid point; add virtual receiver time for TLOC
                trec(i, 1) = t(irx, irz) + geom%recr(i)%t0
            end if
        end do
        !$omp end parallel do

        t = transpose(t)

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(niter)// &
            ', relative diff = '//num2str(ttdiff, '(es)'))

    end subroutine forward_iso_fast_sweep

    !
    !> Compute adjoint field for FATT
    !
    subroutine adjoint_iso(v, d, o, geom, t, tresidual, tadj)

        real, dimension(:, :), intent(in) :: v, t
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :), allocatable, intent(out) :: tadj

        real, allocatable, dimension(:, :) :: vp, lambda_prev
        integer :: nx, nz, iter, i, irx, irz
        real :: lambda_diff, dx, dz, ox, oz

        dx = d(1)
        dz = d(2)
        ox = o(1)
        oz = o(2)

        tt = transpose(t)
        vp = transpose(v)
        nx = size(vp, 1)
        nz = size(vp, 2)

        recrflag = falses(nx, nz)
        lambda = zeros(nx, nz) + sqrt(float_huge)
        lambda_prev = zeros(nx, nz)
        lambda_diff = sqrt(float_huge)

        ! Get the misfit values at the receivers
        !$omp parallel do private(i, irx, irz)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                lambda(irx, irz) = tresidual(i, 1)
                recrflag(irx, irz) = .true.
            end if
        end do
        !$omp end parallel do

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

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(iter)// &
            ', relative diff = '//num2str(lambda_diff, '(es)'))

    end subroutine adjoint_iso

end module traveltime_iso
