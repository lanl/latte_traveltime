!
! © 2025. Triad National Security, LLC. All rights reserved.
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

    double precision, allocatable, dimension(:, :, :) :: t0, pdxt0, pdyt0, pdzt0, tt, lambda
    logical, allocatable, dimension(:, :, :) :: unknown, recrflag

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
    !> with adaptions and modifications for 3D scenario
    !
    subroutine fast_sweep_forward(nx, ny, nz, nxd, nyd, nzd, dx, dy, dz, vp, t0, px0, py0, pz0, tau)

        integer, intent(in) :: nx, ny, nz, nxd, nyd, nzd
        double precision, intent(in) :: dx, dy, dz
        double precision, dimension(:, :, :), intent(in) :: vp, t0, px0, py0, pz0
        double precision, dimension(:, :, :), intent(inout) :: tau

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

            causality1 = tau1*t0c >= taux*t0x .and. tau1*t0c >= tauz*t0z
            causality2 = tau2*t0c >= taux*t0x .and. tau2*t0c >= tauz*t0z

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

    function local_solver_3d(t0c, t0x, t0y, t0z, pdxt0c, pdyt0c, pdzt0c, taux, tauy, tauz, dx, dy, dz, signx, signy, signz, v) result(txyz)

        double precision :: t0c, t0x, t0y, t0z, pdxt0c, pdyt0c, pdzt0c, taux, tauy, tauz
        double precision :: dx, dy, dz, v
        integer :: signx, signy, signz
        double precision :: txyz

        double precision :: a1, a2, a3, b1, b2, b3
        double precision :: a, b, c
        double precision :: tau1, tau2
        logical :: causality1, causality2

        a1 = pdxt0c + t0c/dx*signx
        a2 = pdyt0c + t0c/dy*signy
        a3 = pdzt0c + t0c/dz*signz
        b1 = -taux*t0c/dx*signx
        b2 = -tauy*t0c/dy*signy
        b3 = -tauz*t0c/dz*signz

        a = a1**2 + a2**2 + a3**2
        b = 2*(a1*b1 + a2*b2 + a3*b3)
        c = b1**2 + b2**2 + b3**2 - 1.0/v**2

        if (b**2 - 4*a*c > 0) then

            tau1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
            tau2 = (-b - sqrt(b**2 - 4*a*c))/(2*a)

            causality1 = tau1*t0c >= max(taux*t0x, tauy*t0y, tauz*t0z)
            causality2 = tau2*t0c >= max(taux*t0x, tauy*t0y, tauz*t0z)

            if (causality1 .and. causality2) then
                txyz = min(tau1, tau2)
            else if (causality1 .and. .not. causality2) then
                txyz = tau1
            else if (.not. causality1 .and. causality2) then
                txyz = tau2
            else
                txyz = min( &
                    local_solver_2d(t0c, t0x, t0y, pdxt0c, pdyt0c, taux, tauy, dx, dy, signx, signy, v), &
                    local_solver_2d(t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz, dx, dz, signx, signz, v), &
                    local_solver_2d(t0c, t0y, t0z, pdyt0c, pdzt0c, tauy, tauz, dy, dz, signy, signz, v))
            end if

        else

            txyz = min( &
                local_solver_2d(t0c, t0x, t0y, pdxt0c, pdyt0c, taux, tauy, dx, dy, signx, signy, v), &
                local_solver_2d(t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz, dx, dz, signx, signz, v), &
                local_solver_2d(t0c, t0y, t0z, pdyt0c, pdzt0c, tauy, tauz, dy, dz, signy, signz, v))

        end if

    end function local_solver_3d

    !
    !> Parallel fast sweeping
    !>
    !> The global scheme is based on
    !>      Detrixhe et al., 2013, JCP,
    !>      A parallel fast sweeping method for the eikonal equation
    !>      doi: 10.1016/j.jcp.2012.11.042
    !>
    !> The local scheme is based on
    !>      Fomel et al., 2009, JCP,
    !>      Fast sweeping method for the factored eikonal equation
    !>      doi: 10.1016/j.jcp.2009.05.029
    !>
    !> I made modifications to these schemes for efficiency
    !
    subroutine fast_sweep_forward(nx, ny, nz, nxd, nyd, nzd, dx, dy, dz, vp, t0, px0, py0, pz0, tau)

        integer, intent(in) :: nx, ny, nz, nxd, nyd, nzd
        double precision, intent(in) :: dx, dy, dz
        double precision, dimension(:, :, :), intent(in) :: vp, t0, px0, py0, pz0
        double precision, dimension(:, :, :), intent(inout) :: tau

        integer :: nxbeg, nxend, nybeg, nyend, nzbeg, nzend
        integer :: i, j, k
        integer :: level, jbeg, jend, kbeg, kend
        double precision :: txyz1, txyz2, txyz3, txyz4, txyz5, txyz6, txyz7, txyz8
        logical :: valid1, valid2, valid3, valid4, valid5, valid6, valid7, valid8

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

        ! In 3D, the maximum level is nz + (ny - 1) + (nx - 1)
        do level = 1, nx + ny + nz - 2

            kbeg = ifelse(level <= nx + ny, nzbeg, nzbeg + (level - (nx + ny - 1))*nzd)
            kend = ifelse(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            jbeg = ifelse(level <= nx + ny, nybeg, nybeg + (level - (nx + ny - 1))*nyd)
            jend = ifelse(level <= ny, nybeg + (level - 1)*nyd, nyend)

            !$omp parallel do private(i, j, k, txyz1, txyz2, txyz3, txyz4, txyz5, txyz6, txyz7, txyz8, &
                !$omp valid1, valid2, valid3, valid4, valid5, valid6, valid7, valid8) schedule(dynamic)
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

                    txyz1 = huge_value
                    txyz2 = huge_value
                    txyz3 = huge_value
                    txyz4 = huge_value
                    txyz5 = huge_value
                    txyz6 = huge_value
                    txyz7 = huge_value
                    txyz8 = huge_value

                    valid1 = .true.
                    valid2 = .true.
                    valid3 = .true.
                    valid4 = .true.
                    valid5 = .true.
                    valid6 = .true.
                    valid7 = .true.
                    valid8 = .true.

                    if (i == 1) then
                        valid1 = .false.
                        valid4 = .false.
                        valid5 = .false.
                        valid8 = .false.
                    else if (i == nx) then
                        valid2 = .false.
                        valid3 = .false.
                        valid6 = .false.
                        valid7 = .false.
                    end if
                    if (j == 1) then
                        valid1 = .false.
                        valid2 = .false.
                        valid5 = .false.
                        valid6 = .false.
                    else if (j == ny) then
                        valid3 = .false.
                        valid4 = .false.
                        valid7 = .false.
                        valid8 = .false.
                    end if
                    if (k == 1) then
                        valid1 = .false.
                        valid2 = .false.
                        valid3 = .false.
                        valid4 = .false.
                    else if (k == nz) then
                        valid5 = .false.
                        valid6 = .false.
                        valid7 = .false.
                        valid8 = .false.
                    end if

                    if (valid1) then
                        txyz1 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j - 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j - 1, k), tau(i, j, k - 1), &
                            dx, dy, dz, +1, +1, +1, vp(i, j, k))
                    end if

                    if (valid2) then
                        txyz2 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j - 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j - 1, k), tau(i, j, k - 1), &
                            dx, dy, dz, -1, +1, +1, vp(i, j, k))
                    end if

                    if (valid3) then
                        txyz3 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j + 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j + 1, k), tau(i, j, k - 1), &
                            dx, dy, dz, -1, -1, +1, vp(i, j, k))
                    end if

                    if (valid4) then
                        txyz4 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j + 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j + 1, k), tau(i, j, k - 1), &
                            dx, dy, dz, +1, -1, +1, vp(i, j, k))
                    end if

                    if (valid5) then
                        txyz5 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j - 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j - 1, k), tau(i, j, k + 1), &
                            dx, dy, dz, +1, +1, -1, vp(i, j, k))
                    end if

                    if (valid6) then
                        txyz6 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j - 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j - 1, k), tau(i, j, k + 1), &
                            dx, dy, dz, -1, +1, -1, vp(i, j, k))
                    end if

                    if (valid7) then
                        txyz7 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j + 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j + 1, k), tau(i, j, k + 1), &
                            dx, dy, dz, -1, -1, -1, vp(i, j, k))
                    end if

                    if (valid8) then
                        txyz8 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j + 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j + 1, k), tau(i, j, k + 1), &
                            dx, dy, dz, +1, -1, -1, vp(i, j, k))
                    end if

                    tau(i, j, k) = min(tau(i, j, k), txyz1, txyz2, txyz3, txyz4, txyz5, txyz6, txyz7, txyz8)

                end do
            end do
            !$omp end parallel do

        end do

    end subroutine fast_sweep_forward

#endif

    !
    !> Fast sweeping to solve the adjoint state equation
    !
    subroutine fast_sweep_adjoint(n1, n2, n3, n1d, n2d, n3d, d1, d2, d3, lambda)

        integer, intent(in) :: n1, n2, n3, n1d, n2d, n3d
        double precision, intent(in) :: d1, d2, d3
        double precision, dimension(:, :, :), intent(inout) :: lambda

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

    !
    !> Fast sweeping factorized eikonal solver
    !
    subroutine forward_iso(v, d, o, geom, t, trec)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        double precision, allocatable, dimension(:, :, :) :: tt, tt_prev, vp
        double precision :: dx, dy, dz, ox, oy, oz
        integer :: nx, ny, nz, niter, i, j, k, l, isx, isy, isz, irx, iry, irz, itx, ity, itz
        double precision, allocatable, dimension(:) :: t1, t2, t3, t4
        double precision :: ttdiff, vsource, dsx, dsy, dsz
        integer :: imin, sw
        logical, allocatable, dimension(:) :: source_inside, receiver_inside
        double precision :: time1, time2

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
                source_inside(l) = point_in_domain(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
            end if
        end do
        !$omp end parallel do

        ! Check if receivers are at zero-velocity points
        receiver_inside = falses(geom%nr)
        !$omp parallel do private(i)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                receiver_inside(i) = point_in_domain(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
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
                            vsource = get_point_value_inside(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
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

        ! For timing the modeling, commented out for clarity
        time1 = omp_get_wtime()

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
        !$omp parallel do private(i, irx, iry, irz, l, vsource, dsx, dsy, dsz)
        do i = 1, geom%nr
            if (receiver_inside(i)) then

                ! Linear interpolation; same as for the source points,
                ! here the implementation automatically handles the case of
                ! a receiver falling on integer grid points.
                trec(i, 1) = get_point_value_inside(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], tt)
                do l = 1, geom%ns
                    if (source_inside(l)) then
                        if (within_the_same_grid(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
                                dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), [ox, oy, oz], [dx, dy, dz])) then
                            vsource = get_point_value_inside(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
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

        t = permute(tt, 321)

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

        double precision, allocatable, dimension(:, :, :) :: vp, lambda_prev
        integer :: nx, ny, nz, iter, i, irx, iry, irz !, j, k
        double precision :: lambda_diff, dx, dy, dz, ox, oy, oz

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
        where (.not. recrflag)
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

end module traveltime_iso
