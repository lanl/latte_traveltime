!
! Test driver for adjoint_receiver = nearest and linear of LATTE (3D)
!
!   kernel3d fd   model nx ny nz d geometry delta ix0 ix1 iy0 iy1 iz0 iz1 out
!       Sensitivity of every receiver traveltime of forward_iso to the
!       velocity at the grid points of the window ix0:ix1, iy0:iy1, iz0:iz1
!       (1-based), by central differences with the relative perturbation
!       delta; out holds (nr, nz, ny, nx) float32, zero outside the window,
!       and out.trec the receiver traveltimes.
!   kernel3d each model nx ny nz d geometry mode out
!       Adjoint field of each receiver alone with a unit residual, (nr, nz, ny, nx).
!   kernel3d adj  model nx ny nz d geometry residuals mode out
!       Adjoint field of all receivers with the given residuals, (nz, ny, nx).
!
! The geometry file has LATTE's shot-file format: id, ns, ns lines x y z t0,
! nr, nr lines x y z weight. The grid spacing d is the same in x, y and z,
! and the origin is (0, 0, 0).
!
program kernel3d

    use libflit
    use parameters
    use traveltime_iso

    implicit none

    character(len=1024) :: task, file_model, file_geom, file_out, file_res, arg, mode
    integer :: gnx, gny, gnz, i, j, k, l, m, nr, w(6), nw
    real :: dd, delta
    real, allocatable, dimension(:, :, :) :: v, vv, t, lam, buf
    real, allocatable, dimension(:, :) :: trec, trec_p, trec_m, res
    real, allocatable, dimension(:, :, :, :) :: kern
    type(source_receiver_geometry) :: geom, g1

    call mpistart

    call get_command_argument(1, task)
    call get_command_argument(2, file_model)
    call get_command_argument(3, arg)
    read (arg, *) gnx
    call get_command_argument(4, arg)
    read (arg, *) gny
    call get_command_argument(5, arg)
    read (arg, *) gnz
    call get_command_argument(6, arg)
    read (arg, *) dd
    call get_command_argument(7, file_geom)

    ! Converge the sweeps fully
    sweep_stop_threshold = 1.0e-12
    sweep_niter_max = 200

    v = load(file_model, gnz, gny, gnx)
    call read_geometry(file_geom, geom)
    nr = geom%nr

    select case (task)

        case ('fd')
            call get_command_argument(8, arg)
            read (arg, *) delta
            do m = 1, 6
                call get_command_argument(8 + m, arg)
                read (arg, *) w(m)
            end do
            call get_command_argument(15, file_out)
            allocate (kern(1:nr, 1:gnz, 1:gny, 1:gnx))
            kern = 0
            nw = (w(2) - w(1) + 1)*(w(4) - w(3) + 1)*(w(6) - w(5) + 1)
            do m = rankid + 1, nw, nrank
                ! z fastest, then y, then x within the window
                k = w(5) + mod(m - 1, w(6) - w(5) + 1)
                j = w(3) + mod((m - 1)/(w(6) - w(5) + 1), w(4) - w(3) + 1)
                i = w(1) + (m - 1)/((w(6) - w(5) + 1)*(w(4) - w(3) + 1))
                vv = v
                vv(k, j, i) = v(k, j, i)*(1 + delta)
                call forward_iso(vv, [dd, dd, dd], [0.0, 0.0, 0.0], geom, t, trec_p)
                vv(k, j, i) = v(k, j, i)*(1 - delta)
                call forward_iso(vv, [dd, dd, dd], [0.0, 0.0, 0.0], geom, t, trec_m)
                kern(:, k, j, i) = real((dble(trec_p(:, 1)) - dble(trec_m(:, 1)))/(2*dble(v(k, j, i))*delta))
            end do
            do l = 1, nr
                buf = kern(l, :, :, :)
                call allreduce_array(buf)
                kern(l, :, :, :) = buf
            end do
            if (rankid == 0) then
                call forward_iso(v, [dd, dd, dd], [0.0, 0.0, 0.0], geom, t, trec)
                call output_array(kern, file_out)
                call output_array(trec, tidy(file_out)//'.trec')
            end if

        case ('each')
            call get_command_argument(8, mode)
            call get_command_argument(9, file_out)
            call forward_iso(v, [dd, dd, dd], [0.0, 0.0, 0.0], geom, t, trec)
            allocate (kern(1:nr, 1:gnz, 1:gny, 1:gnx))
            kern = 0
            res = ones(1, 1)
            do l = rankid + 1, nr, nrank
                g1 = geom
                g1%nr = 1
                g1%recr = geom%recr(l:l)
                call adjoint_iso(v, [dd, dd, dd], [0.0, 0.0, 0.0], g1, t, res, lam, tidy(mode))
                kern(l, :, :, :) = lam
            end do
            do l = 1, nr
                buf = kern(l, :, :, :)
                call allreduce_array(buf)
                kern(l, :, :, :) = buf
            end do
            if (rankid == 0) then
                call output_array(kern, file_out)
            end if

        case ('adj')
            call get_command_argument(8, file_res)
            call get_command_argument(9, mode)
            call get_command_argument(10, file_out)
            res = load(file_res, nr, 1)
            call forward_iso(v, [dd, dd, dd], [0.0, 0.0, 0.0], geom, t, trec)
            call adjoint_iso(v, [dd, dd, dd], [0.0, 0.0, 0.0], geom, t, res, lam, tidy(mode))
            if (rankid == 0) then
                call output_array(lam, file_out)
            end if

    end select

    call mpibarrier
    call mpiend

contains

    subroutine read_geometry(file, g)

        character(len=*), intent(in) :: file
        type(source_receiver_geometry), intent(out) :: g

        integer :: u, n

        open (newunit=u, file=file, action='read')
        read (u, *) g%id
        read (u, *) g%ns
        allocate (g%srcr(1:g%ns))
        do n = 1, g%ns
            read (u, *) g%srcr(n)%x, g%srcr(n)%y, g%srcr(n)%z, g%srcr(n)%t0
            g%srcr(n)%amp = 1
        end do
        read (u, *) g%nr
        allocate (g%recr(1:g%nr))
        do n = 1, g%nr
            read (u, *) g%recr(n)%x, g%recr(n)%y, g%recr(n)%z, g%recr(n)%weight
            g%recr(n)%aoff = 0
            g%recr(n)%t0 = 0
        end do
        close (u)

    end subroutine read_geometry

end program kernel3d
