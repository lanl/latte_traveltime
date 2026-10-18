!
! Test driver for adjoint_receiver = nearest and linear of LATTE (2D)
!
!   kernel2d fd   model nx nz dx dz geometry delta out
!       Sensitivity of every receiver traveltime of forward_iso to the
!       velocity at every grid point, by central differences with the
!       relative perturbation delta; out holds (nr, nz, nx) float32,
!       out.t the traveltime field and out.trec the receiver traveltimes.
!   kernel2d fwd  model nx nz dx dz geometry out
!       Receiver traveltimes of forward_iso, (nr).
!   kernel2d each model nx nz dx dz geometry mode out
!       Adjoint field of each receiver alone with a unit residual, (nr, nz, nx).
!   kernel2d adj  model nx nz dx dz geometry residuals mode out
!       Adjoint field of all receivers with the given residuals, (nz, nx).
!
! The geometry file has LATTE's shot-file format: id, ns, ns lines x y z t0,
! nr, nr lines x y z weight. The grid origin is (0, 0).
!
program kernel2d

    use libflit
    use parameters
    use traveltime_iso

    implicit none

    character(len=1024) :: task, file_model, file_geom, file_out, file_res, arg, mode
    integer :: gnx, gnz, i, j, k, l, nr
    real :: ddx, ddz, delta
    real, allocatable, dimension(:, :) :: v, vv, t, trec, trec_p, trec_m, lam, res
    real, allocatable, dimension(:, :, :) :: kern
    type(source_receiver_geometry) :: geom, g1

    call mpistart

    call get_command_argument(1, task)
    call get_command_argument(2, file_model)
    call get_command_argument(3, arg)
    read (arg, *) gnx
    call get_command_argument(4, arg)
    read (arg, *) gnz
    call get_command_argument(5, arg)
    read (arg, *) ddx
    call get_command_argument(6, arg)
    read (arg, *) ddz
    call get_command_argument(7, file_geom)

    ! Converge the sweeps fully
    sweep_stop_threshold = 1.0e-12
    sweep_niter_max = 200

    v = load(file_model, gnz, gnx)
    call read_geometry(file_geom, geom)
    nr = geom%nr

    select case (task)

        case ('fd')
            call get_command_argument(8, arg)
            read (arg, *) delta
            call get_command_argument(9, file_out)
            allocate (kern(1:nr, 1:gnz, 1:gnx))
            kern = 0
            do k = rankid + 1, gnz*gnx, nrank
                i = mod(k - 1, gnz) + 1
                j = (k - 1)/gnz + 1
                vv = v
                vv(i, j) = v(i, j)*(1 + delta)
                call forward_iso(vv, [ddx, ddz], [0.0, 0.0], geom, t, trec_p)
                vv(i, j) = v(i, j)*(1 - delta)
                call forward_iso(vv, [ddx, ddz], [0.0, 0.0], geom, t, trec_m)
                kern(:, i, j) = real((dble(trec_p(:, 1)) - dble(trec_m(:, 1)))/(2*dble(v(i, j))*delta))
            end do
            call allreduce_array(kern)
            if (rankid == 0) then
                call forward_iso(v, [ddx, ddz], [0.0, 0.0], geom, t, trec)
                call output_array(kern, file_out)
                call output_array(t, tidy(file_out)//'.t')
                call output_array(trec, tidy(file_out)//'.trec')
            end if

        case ('fwd')
            call get_command_argument(8, file_out)
            call forward_iso(v, [ddx, ddz], [0.0, 0.0], geom, t, trec)
            if (rankid == 0) then
                call output_array(trec, file_out)
            end if

        case ('each')
            call get_command_argument(8, mode)
            call get_command_argument(9, file_out)
            call forward_iso(v, [ddx, ddz], [0.0, 0.0], geom, t, trec)
            allocate (kern(1:nr, 1:gnz, 1:gnx))
            kern = 0
            res = ones(1, 1)
            do l = rankid + 1, nr, nrank
                g1 = geom
                g1%nr = 1
                g1%recr = geom%recr(l:l)
                call adjoint_iso(v, [ddx, ddz], [0.0, 0.0], g1, t, res, lam, tidy(mode))
                kern(l, :, :) = lam
            end do
            call allreduce_array(kern)
            if (rankid == 0) then
                call output_array(kern, file_out)
            end if

        case ('adj')
            call get_command_argument(8, file_res)
            call get_command_argument(9, mode)
            call get_command_argument(10, file_out)
            res = load(file_res, nr, 1)
            call forward_iso(v, [ddx, ddz], [0.0, 0.0], geom, t, trec)
            call adjoint_iso(v, [ddx, ddz], [0.0, 0.0], geom, t, res, lam, tidy(mode))
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

        integer :: u, m
        real :: y

        open (newunit=u, file=file, action='read')
        read (u, *) g%id
        read (u, *) g%ns
        allocate (g%srcr(1:g%ns))
        do m = 1, g%ns
            read (u, *) g%srcr(m)%x, y, g%srcr(m)%z, g%srcr(m)%t0
            g%srcr(m)%y = 0
            g%srcr(m)%amp = 1
        end do
        read (u, *) g%nr
        allocate (g%recr(1:g%nr))
        do m = 1, g%nr
            read (u, *) g%recr(m)%x, y, g%recr(m)%z, g%recr(m)%weight
            g%recr(m)%y = 0
            g%recr(m)%aoff = 0
            g%recr(m)%t0 = 0
        end do
        close (u)

    end subroutine read_geometry

end program kernel2d
