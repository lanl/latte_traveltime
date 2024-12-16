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

    real, parameter :: huge_value = sqrt(float_huge)

    real, allocatable, dimension(:, :, :) :: t0, pdxt0, pdyt0, pdzt0, tt, lambda
    logical, allocatable, dimension(:, :, :) :: unknown, recrflag

    private
    public :: forward_iso
    public :: adjoint_iso

contains

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

            kbeg = ifelse(level <= nx + ny, nzbeg, nzbeg + (level - (nx + ny - 1))*nzd)
            kend = ifelse(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            jbeg = ifelse(level <= nx + ny, nybeg, nybeg + (level - (nx + ny - 1))*nyd)
            jend = ifelse(level <= ny, nybeg + (level - 1)*nyd, nyend)

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

                    if (taux == huge_value .and. tauy == huge_value .and. tauz == huge_value) then
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

                    if (u(1) == huge_value .and. u(2) == huge_value .and. u(3) == huge_value) then
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

                            ! Solve the quadratic equation analytically
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

            kbeg = ifelse(level <= n1 + n2, n3beg, n3beg + (level - (n1 + n2 - 1))*n3d)
            kend = ifelse(level <= n3, n3beg + (level - 1)*n3d, n3end)

            jbeg = ifelse(level <= n1 + n2, n2beg, n2beg + (level - (n1 + n2 - 1))*n2d)
            jend = ifelse(level <= n2, n2beg + (level - 1)*n2d, n2end)

            !$omp parallel do private(i, j, k, i1, i2, j1, j2, k1, k2, &
                !$omp app, amp, apm, amm, &
                !$omp bpp, bmp, bpm, bmm, &
                !$omp cpp, cmp, cpm, cmm, &
                !$omp ap, am, bp, bm, cp, cm, &
                !$omp lhs, rhs, t) collapse(2)
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
        double precision :: da, db
        double precision :: signa, signb
        double precision :: dadt, dbdt
        double precision :: taua, taub, t0a, t0b
        double precision :: tauca, taucb
        double precision :: pa0c, pb0c
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

                    if (taux == huge_value .and. tauy == huge_value .and. tauz == huge_value) then
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

                    if (u(1) == huge_value .and. u(2) == huge_value .and. u(3) == huge_value) then
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

                            ! Solve the quadratic equation analytically
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
    subroutine forward_iso(v, d, o, geom, t, trec)

        use omp_lib

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        real, allocatable, dimension(:, :, :) :: tt, tt_prev, vp
        real :: dx, dy, dz, ox, oy, oz
        integer :: nx, ny, nz, niter, i, j, k, l, isx, isy, isz, irx, iry, irz, itx, ity, itz
        real, allocatable, dimension(:) :: t1, t2, t3, t4
        real :: ttdiff, vsource, dsx, dsy, dsz
        integer :: imin, sw
        logical, allocatable, dimension(:) :: source_inside, receiver_inside
        !        real :: time1, time2

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

        ! Check if sources are at zero-velocity points
        source_inside = falses(geom%ns)
        !$omp parallel do private(l)
        do l = 1, geom%ns
            if (geom%srcr(l)%amp == 1) then
                source_inside(l) = point_in_domain([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z], &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
            end if
        end do
        !$omp end parallel do

        ! Check if receivers are at zero-velocity points
        receiver_inside = falses(geom%nr)
        !$omp parallel do private(i)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                receiver_inside(i) = point_in_domain([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z], &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
            end if
        end do
        !$omp end parallel do

        t0 = zeros(nx, ny, nz)
        pdxt0 = zeros(nx, ny, nz)
        pdyt0 = zeros(nx, ny, nz)
        pdzt0 = zeros(nx, ny, nz)

        ! Source location
        t1 = zeros(geom%ns) + huge_value
        t2 = zeros_like(t1)
        t3 = zeros_like(t1)
        t4 = zeros_like(t1)

        ! Background time field
        !$omp parallel do private(i, j, k, l, isx, isy, isz, vsource, dsx, dsy, dsz, imin, t1, t2, t3, t4)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    do l = 1, geom%ns

                        if (source_inside(l)) then

                            ! Linear interpolation to get velocity at the source position
                            vsource = get_point_value_inside([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z], &
                                [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)

                            dsx = ox + (i - 1)*dx - geom%srcr(l)%x
                            dsy = oy + (j - 1)*dy - geom%srcr(l)%y
                            dsz = oz + (k - 1)*dz - geom%srcr(l)%z

                            !                            ! In comparison, the nearest grid point apporach
                            !                            isx = nint((geom%srcr(l)%x - ox)/dx) + 1
                            !                            isy = nint((geom%srcr(l)%y - oy)/dy) + 1
                            !                            isz = nint((geom%srcr(l)%z - oz)/dz) + 1
                            !                            vsource = vp(isx, isy, isz)
                            !                            dsx = (i - isx)*dx
                            !                            dsy = (j - isy)*dy
                            !                            dsz = (k - isz)*dz

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

        ! Initialize the multiplicative time field
        tt_prev = zeros(nx, ny, nz)
        tt = zeros(nx, ny, nz) + huge_value
        ! Increase the width of the multiplicative field around the source
        ! in the initialization seems to significantly improve the resulting accuracy,
        ! especially when the source point is not on an integer grid point
        sw = 3
        !$omp parallel do private(i, isx, isy, isz, itx, ity, itz)
        do i = 1, geom%ns
            if (source_inside(i)) then
                isx = max(floor((geom%srcr(i)%x - ox)/dx) + 1 - sw, 1)
                isy = max(floor((geom%srcr(i)%y - oy)/dy) + 1 - sw, 1)
                isz = max(floor((geom%srcr(i)%z - oz)/dz) + 1 - sw, 1)
                itx = min(ceiling((geom%srcr(i)%x - ox)/dx) + 1 + sw, nx)
                ity = min(ceiling((geom%srcr(i)%y - oy)/dy) + 1 + sw, ny)
                itz = min(ceiling((geom%srcr(i)%z - oz)/dz) + 1 + sw, nz)
                tt(isx:itx, isy:ity, isz:itz) = 1.0
            end if
        end do
        !$omp end parallel do

        !        ! As a comparison, the nearest grid point approach is
        !        !$omp parallel do private(i, isx, isy, isz, itx, ity, itz)
        !        do i = 1, geom%ns
        !            if (source_inside(i)) then
        !                isx = nint((geom%srcr(i)%x - ox)/dx) + 1
        !                isy = nint((geom%srcr(i)%y - oy)/dy) + 1
        !                isz = nint((geom%srcr(i)%z - oz)/dz) + 1
        !                tt(isx, isy, isz) = 1.0
        !            end if
        !        end do
        !        !$omp end parallel do

        ! Fast sweeping via Gauss-Seidel iterations
        niter = 0
        ttdiff = huge_value
        where (vp == 0)
            vp = float_tiny
        end where

        !        ! For timing the modeling, commented out for clarity
        !        call cpu_time(time1)

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

        !        ! For timing the modeling, commented out for clarity
        !        call cpu_time(time2)
        !        if (file_exists('./time3.txt')) then
        !            open(3, file='./time3.txt', position='append')
        !        else
        !            open(3, file='./time3.txt')
        !        end if
        !        write (3, *) time2 - time1
        !        close(3)

        t = t0*tt
        where (vp == float_tiny)
            t = 0
            vp = 0
        end where

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        !$omp parallel do private(i, irx, iry, irz, l, vsource, dsx, dsy, dsz)
        do i = 1, geom%nr
            if (receiver_inside(i)) then

                ! Linear interpolation; same as for the source points,
                ! here the implementation automatically handles the case of
                ! a receiver falling on integer grid points.
                trec(i, 1) = get_point_value_inside([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z], &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], t)
                do l = 1, geom%ns
                    if (source_inside(l)) then
                        if (within_the_same_grid([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z], &
                                [geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z], [ox, oy, oz], [dx, dy, dz])) then
                            vsource = get_point_value_inside([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z], &
                                [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
                            dsx = geom%recr(i)%x - geom%srcr(l)%x
                            dsy = geom%recr(i)%y - geom%srcr(l)%y
                            dsz = geom%recr(i)%z - geom%srcr(l)%z
                            trec(i, 1) = sqrt(dsx**2 + dsy**2 + dsz**2)/vsource + geom%srcr(l)%t0
                        end if
                    end if
                end do

                ! Add virtual receiver time for TLOC
                trec(i, 1) = trec(i, 1) + geom%recr(i)%t0

                !                ! As a comparison, the nearest grid point approach is
                !                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                !                iry = nint((geom%recr(i)%y - oy)/dy) + 1
                !                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                !                trec(i, 1) = t(irx, iry, irz) + geom%recr(i)%t0

            end if
        end do
        !$omp end parallel do

        t = permute(t, 321)

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(niter)// &
            ', relative diff = '//num2str(ttdiff, '(es)'))

    end subroutine forward_iso

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
        integer :: nx, ny, nz, iter, i, irx, iry, irz !, j, k
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
        lambda = zeros(nx, ny, nz)
        lambda_prev = zeros(nx, ny, nz)
        lambda_diff = huge_value

        ! Get the misfit values at the receivers
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                iry = nint((geom%recr(i)%y - oy)/dy) + 1
                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                lambda(irx, iry, irz) = lambda(irx, iry, irz) + tresidual(i, 1)
                recrflag(irx, iry, irz) = .true.
            end if
        end do
        where (.not.recrflag)
            lambda = huge_value
        end where

        !        ! This is to set lambda = 0 on boundaries, but does not make a difference
        !        ! when both source and receivers are inside and could be ignored.
        !        ! But when this is set, the preconditioner input residual must be negative to
        !        ! avoid unwanted boundary artifacts.
        !        !$omp parallel do private(i, j, k)
        !        do k = 1, nz
        !            do j = 1, ny
        !                do i = 1, nx
        !                    if ((i == 1 .or. i == nx .or. j == 1 .or. j == ny &
            !                        .or. k == 1 .or. k == nz) .and. (.not.recrflag(i, j, k))) then
        !                        lambda(i, j, k) = 0.0
        !                    end if
        !                end do
        !            end do
        !        end do
        !        !$omp end parallel do

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

        call warn(date_time_compact()//' Fast sweeping adjoint niter = '//num2str(iter)// &
            ', relative diff = '//num2str(lambda_diff, '(es)'))

    end subroutine adjoint_iso

    !
    ! Check if a point is inside of the model domain
    !
    function point_in_domain(p, n, o, d, v) result(f)

        real, dimension(:) :: p, o, d
        integer, dimension(:) :: n
        real, dimension(:, :, :) :: v
        logical :: f

        real :: p1, p2, p3
        integer :: n1, n2, n3
        real :: d1, d2, d3, o1, o2, o3
        integer, dimension(1:2) :: ii, jj, kk
        integer :: i, j, k

        f = .true.

        n1 = n(1)
        n2 = n(2)
        n3 = n(3)
        d1 = d(1)
        d2 = d(2)
        d3 = d(3)
        o1 = o(1)
        o2 = o(2)
        o3 = o(3)
        p1 = p(1) - o1
        p2 = p(2) - o2
        p3 = p(3) - o3

        if (p1 < 0 .or. p1 > (n1 - 1)*d1 &
                .or. p2 < 0 .or. p2 > (n2 - 1)*d2 &
                .or. p3 < 0 .or. p3 > (n3 - 1)*d3) then
            f = .false.
            return
        end if

        ii = [floor(p1/d1) + 1, ceiling(p1/d1) + 1]
        jj = [floor(p2/d2) + 1, ceiling(p2/d2) + 1]
        kk = [floor(p3/d3) + 1, ceiling(p3/d3) + 1]
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    if (v(ii(i), jj(j), kk(k)) == 0) then
                        f = .false.
                        return
                    end if
                end do
            end do
        end do

    end function point_in_domain

    !
    ! Linear interpolation to get the value of a point inside a grid
    !
    function get_point_value_inside(p, n, o, d, v) result(f)

        real, dimension(:) :: p, o, d
        integer, dimension(:) :: n
        real, dimension(:, :, :) :: v
        real :: f

        real :: p1, p2, p3
        integer :: n1, n2, n3
        real :: d1, d2, d3, o1, o2, o3
        integer :: i1beg, i1end
        integer :: i2beg, i2end
        integer :: i3beg, i3end

        f = .true.

        n1 = n(1)
        n2 = n(2)
        n3 = n(3)
        d1 = d(1)
        d2 = d(2)
        d3 = d(3)
        o1 = o(1)
        o2 = o(2)
        o3 = o(3)
        p1 = p(1) - o1
        p2 = p(2) - o2
        p3 = p(3) - o3

        i1beg = floor(p1/d1) + 1
        i1end = ceiling(p1/d1) + 1
        i2beg = floor(p2/d2) + 1
        i2end = ceiling(p2/d2) + 1
        i3beg = floor(p3/d3) + 1
        i3end = ceiling(p3/d3) + 1

        if (i1beg == i1end .and. i2beg == i2end .and. i3beg == i3end) then

            f = v(i1beg, i2beg, i3beg)

        else if (i1beg == i1end .and. i2beg == i2end .and. i3beg /= i3end) then

            f = point_interp_linear( &
                ([i3beg, i3end] - 1)*d3, &
                v(i1beg, i2beg, i3beg:i3end), &
                p3)

        else if (i1beg == i1end .and. i2beg /= i2end .and. i3beg == i3end) then

            f = point_interp_linear( &
                ([i2beg, i2end] - 1)*d2, &
                v(i1beg, i2beg:i2end, i3beg), &
                p2)

        else if (i1beg /= i1end .and. i2beg == i2end .and. i3beg == i3end) then

            f = point_interp_linear( &
                ([i1beg, i1end] - 1)*d1, &
                v(i1beg:i1end, i2beg, i3beg), &
                p1)

        else if (i1beg == i1end .and. i2beg /= i2end .and. i3beg /= i3end) then

            f = point_interp_linear( &
                ([i2beg, i2end] - 1)*d2, &
                ([i3beg, i3end] - 1)*d3, &
                v(i1beg, i2beg:i2end, i3beg:i3end), &
                p2, p3)

        else if (i1beg /= i1end .and. i2beg == i2end .and. i3beg /= i3end) then

            f = point_interp_linear( &
                ([i1beg, i1end] - 1)*d1, &
                ([i3beg, i3end] - 1)*d3, &
                v(i1beg:i1end, i2beg, i3beg:i3end), &
                p1, p3)

        else if (i1beg /= i1end .and. i2beg /= i2end .and. i3beg == i3end) then

            f = point_interp_linear( &
                ([i1beg, i1end] - 1)*d1, &
                ([i2beg, i2end] - 1)*d2, &
                v(i1beg:i1end, i2beg:i2end, i3beg), &
                p1, p2)

        else

            f = point_interp_linear( &
                ([i1beg, i1end] - 1)*d1, &
                ([i2beg, i2end] - 1)*d2, &
                ([i3beg, i3end] - 1)*d3, &
                v(i1beg:i1end, i2beg:i2end, i3beg:i3end), &
                p1, p2, p3)

        end if

    end function get_point_value_inside

    !
    ! Check if two points are in the same grid
    !
    function within_the_same_grid(pa, pb, o, d) result(f)

        real, dimension(:) :: pa, pb, o, d
        logical :: f

        real :: p1_a, p2_a, p3_a
        real :: p1_b, p2_b, p3_b
        real :: d1, d2, d3, o1, o2, o3
        integer :: i1beg_a, i1end_a
        integer :: i2beg_a, i2end_a
        integer :: i3beg_a, i3end_a
        integer :: i1beg_b, i1end_b
        integer :: i2beg_b, i2end_b
        integer :: i3beg_b, i3end_b

        d1 = d(1)
        d2 = d(2)
        d3 = d(3)
        o1 = o(1)
        o2 = o(2)
        o3 = o(3)
        p1_a = pa(1) - o1
        p2_a = pa(2) - o2
        p3_a = pa(3) - o3
        p1_b = pb(1) - o1
        p2_b = pb(2) - o2
        p3_b = pb(3) - o3

        i1beg_a = floor(p1_a/d1) + 1
        i1end_a = ceiling(p1_a/d1) + 1
        i2beg_a = floor(p2_a/d2) + 1
        i2end_a = ceiling(p2_a/d2) + 1
        i3beg_a = floor(p3_a/d3) + 1
        i3end_a = ceiling(p3_a/d3) + 1

        i1beg_b = floor(p1_b/d1) + 1
        i1end_b = ceiling(p1_b/d1) + 1
        i2beg_b = floor(p2_b/d2) + 1
        i2end_b = ceiling(p2_b/d2) + 1
        i3beg_b = floor(p3_b/d3) + 1
        i3end_b = ceiling(p3_b/d3) + 1

        f = i1beg_a == i1beg_b &
            .and. i1end_a == i1end_b &
            .and. i2beg_a == i2beg_b &
            .and. i2end_a == i2end_b &
            .and. i3beg_a == i3beg_b &
            .and. i3end_a == i3end_b

    end function within_the_same_grid

end module traveltime_iso
