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


module gradient

    use libflit
    use parameters
    use vars
    use utility

    implicit none

    real :: refltaper_x, refltaper_y, reflwindow_x, reflwindow_y, reflmin, reflmax
    real, allocatable, dimension(:, :, :, :) :: refl_all

    interface
        module subroutine compute_gradient_shots_fatt
        end subroutine
        module subroutine compute_gradient_shots_trtt
        end subroutine
        module subroutine compute_gradient_shots_tloc
        end subroutine
        module subroutine compute_gradient_shots_tloc_dd
        end subroutine
        module subroutine compute_image_shots
        end subroutine
    end interface

    private

    public :: zero_gradient
    public :: compute_gradient_shots
    public :: process_gradient
    public :: output_gradient
    public :: compute_image_shots

contains

    !
    !> Initialize gradient arrays
    !
    subroutine zero_gradient

        integer :: i

        do i = 1, nmodel
            select case (model_name(i))
                case default
                    model_grad(i)%array = zeros(nz, ny, nx)
                case ('sx', 'sy', 'sz', 'st0')
                    model_grad(i)%array = zeros(nr_virtual, 1, 1)
            end select
        end do

        if (yn_update_reflector) then
            refl_all = zeros(nz, ny, nx, nrefl)
        end if

    end subroutine zero_gradient

    !
    !> Compute gradients shot by shot
    !
    subroutine compute_gradient_shots

        select case (which_program)
            case ('fatt')
                call compute_gradient_shots_fatt
            case ('trtt')
                call compute_gradient_shots_trtt
            case ('tloc')
                if (yn_dd_no_st0) then
                    call compute_gradient_shots_tloc_dd
                else
                    call compute_gradient_shots_tloc
                end if
        end select

    end subroutine compute_gradient_shots

    !
    !> Process gradient
    !
    subroutine process_gradient

        integer :: i, l

        do i = 1, nmodel

            if (.not. any(model_m(i)%name == ['sx', 'sy', 'sz', 'st0'])) then

                if (model_m(i)%name /= 'refl') then

                    if (uniform_processing) then
                        call process_model_single_parameter(model_grad(i)%array, 'grad', param_name=model_m(i)%name)
                    else
                        call process_model_single_parameter(model_grad(i)%array, 'grad_'//tidy(model_name(i)), param_name=model_m(i)%name)
                    end if

                else

                    do l = 1, nrefl
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_taper_length_x', refltaper_x, 10*dx, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_taper_length_y', refltaper_y, 10*dy, iter*1.0)
                        call assert(refltaper_x < nx*dx/3.0 .and. refltaper_y < ny*dy/3.0, ' <compute_image_shots> Error: tapers are too long. ')
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_smooth_window_x', reflwindow_x, 0.2*(nx*dx - 2*refltaper_x), iter*1.0)
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_smooth_window_y', reflwindow_y, 0.2*(ny*dy - 2*refltaper_y), iter*1.0)
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_min', reflmin, oz, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_max', reflmax, oz + (nz - 1)*dz, iter*1.0)
                        call smooth_reflector(refl_all(:, :, :, l))
                        model_grad(i)%array = model_grad(i)%array - refl_all(:, :, :, l)*l
                    end do
                    call process_model_single_parameter(model_grad(i)%array, 'refl')

                end if

            end if

        end do

    end subroutine process_gradient

    !
    !> Output gradient
    !
    subroutine output_gradient

        integer :: i

        if (rankid == 0) then

            do i = 1, nmodel
                call output_array(model_grad(i)%array, &
                    dir_iter_model(iter)//'/grad_'//tidy(model_name(i))//'.bin')
            end do

            if (nmodel > 1) then
                call warn(date_time_compact()//' >>>>>>>>>> Gradients are saved. ')
            else
                call warn(date_time_compact()//' >>>>>>>>>> Gradient is saved. ')
            end if

        end if

    end subroutine output_gradient

    !
    !> Specially process reflector image
    !
    subroutine smooth_reflector(img)

        real, dimension(:, :, :), intent(inout) :: img

        integer :: tplenx, tpleny
        real :: windowx, windowy
        real, allocatable, dimension(:, :) :: x, y, z, xx, yy, zz
        integer :: n1, n2, n3, i, j, k
        integer :: imin
        real, allocatable, dimension(:, :, :) :: v

        n1 = size(img, 1)
        n2 = size(img, 2)
        n3 = size(img, 3)
        tplenx = nint(refltaper_x/dx)
        tpleny = nint(refltaper_y/dy)
        windowx = reflwindow_x/dx/n3
        windowy = reflwindow_y/dy/n2

        v = img

        x = reshape(meshgrid([n2, n3], [1.0, 1.0], [0.0, 0.0], dim=1), [n2, n3])
        y = reshape(meshgrid([n2, n3], [1.0, 1.0], [0.0, 0.0], dim=2), [n2, n3])
        z = zeros(n2, n3)

        do k = 1, n3
            do j = 1, n2
                imin = as_scalar(maxloc(v(:, j, k)))
                do i = 1, n1
                    if (i == imin) then
                        v(imin, j, k) = 1.0
                        z(j, k) = imin
                    else
                        v(i, j, k) = 0.0
                    end if
                end do
            end do
        end do
        z = clip(z, nint((reflmin - oz)/dz) + 1.0, nint((reflmax - oz)/dz) + 1.0)

        xx = x(tpleny + 1:n2 - tpleny, tplenx + 1:n3 - tplenx)
        yy = y(tpleny + 1:n2 - tpleny, tplenx + 1:n3 - tplenx)
        zz = z(tpleny + 1:n2 - tpleny, tplenx + 1:n3 - tplenx)

        ! The window should be a fraction of the total length of xx
        zz = reshape(lowess_filt(flatten(yy), flatten(xx), flatten(zz), flatten(yy), flatten(xx), &
            window=[windowy, windowx], robust=.true.), [n2 - 2*tpleny, n3 - 2*tplenx])
        z = interp(zz, [n2 - 2*tpleny, n3 - 2*tplenx], [1.0, 1.0], [tpleny + 1.0, tplenx + 1.0], &
            [n2, n3], [1.0, 1.0], [1.0, 1.0], ['linear', 'linear'])

        z = clip(z, nint((reflmin - oz)/dz) + 1.0, nint((reflmax - oz)/dz) + 1.0)

        v = 0
        do k = 1, n3
            do j = 1, n2
                v(nint(z(j, k)), j, k) = 1.0
            end do
        end do

        img = v

    end subroutine smooth_reflector

    !
    !> Process gradient of one shot
    !
    subroutine process_model_single_shot(ishot, w, ws, name)

        integer, intent(in) :: ishot
        real, dimension(:, :, :), intent(inout) :: w, ws
        character(len=*), intent(in) :: name

        integer :: i, j, k
        character(len=1024) :: shot_prefix, dir_mask
        real :: shot_w_movingbalx, shot_w_movingbaly, shot_w_movingbalz
        real, allocatable, dimension(:) :: shot_w_taperx, shot_w_tapery, shot_w_taperz
        real :: shot_w_medianfiltx, shot_w_medianfilty, shot_w_medianfiltz
        real :: shot_w_smoothx, shot_w_smoothy, shot_w_smoothz
        real, allocatable, dimension(:) :: wavenums, wamps, fkdips, fkdipamps
        real, allocatable, dimension(:, :, :) :: w_mask
        character(len=32), allocatable, dimension(:) :: process_shot_grad_w
        real, allocatable, dimension(:) :: st
        real :: tp
        real :: conemuter, conemutez, conemutepower, conemutetaper
        character(len=24) :: conemuteorigin
        real :: srcx, srcy, srcz, depth, ds, dist
        integer :: pr, pt, ix, iy
        character(len=1024) :: file_mask

        call readpar_nstring(file_parameter, 'process_shot_'//tidy(name)//'', process_shot_grad_w, [''])

        do i = 1, size(process_shot_grad_w)

            select case (process_shot_grad_w(i))

                case ('smooth')
                    ! Gaussian smooth
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothx', shot_w_smoothx, 3*dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothy', shot_w_smoothy, 3*dy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothz', shot_w_smoothz, 3*dz)
                    ws = gauss_filt(ws, [shot_w_smoothz/dz, shot_w_smoothy/dy, shot_w_smoothx/dx])

                case ('maxbal')
                    if (maxval(ws) /= 0) then
                        ws = ws/maxval(ws)
                    end if

                case ('rmsbal')
                    ! Normalize with shot image energy
                    ws = ws/mean(ws, 2)

                case ('movingbal')
                    ! Moving balance
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbalx', shot_w_movingbalx, 6*dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbaly', shot_w_movingbaly, 6*dy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbalz', shot_w_movingbalz, 6*dz)
                    ws = balance_filt(ws, nint([0.5*shot_w_movingbalz/dz, 0.5*shot_w_movingbaly/dy, 0.5*shot_w_movingbalx/dx]), 0.01)

                case ('medianfilt')
                    ! Median filtering
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfiltx', shot_w_medianfiltx, dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfilty', shot_w_medianfilty, dy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfiltz', shot_w_medianfiltz, dz)
                    ws = median_filt(ws, nint([shot_w_medianfiltz/dz, shot_w_medianfilty/dy, shot_w_medianfiltx/dx]))

                case ('dipfilt')
                    ! Dip filtering
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzx', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzx_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        ws = dip_filt(ws, [1.0, dy/dz, dx/dz], fkdips, fkdipamps, 13)
                    end if

                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzy', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzy_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        ws = dip_filt(ws, [1.0, dy/dz, dx/dz], fkdips, fkdipamps, 12)
                    end if

                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltyx', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltyx_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        ws = dip_filt(ws, [1.0, dy/dz, dx/dz], fkdips, fkdipamps, 23)
                    end if

                case ('removenan')
                    ! Remove NaN
                    ws = return_normal(ws)

                case ('laplacefilt')
                    ws = laplace_filt(ws)

                case ('wavenumfilt')
                    ! Wavenumber-domain filtering in x-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumx', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumx_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j, k)
                        do k = 1, size(ws, 1)
                            do j = 1, size(ws, 2)
                                ws(k, j, :) = fourier_filt(ws(k, j, :), dx, wavenums, wamps)
                            end do
                        end do
                        !$omp end parallel do
                    end if

                    ! Wavenumber-domain filtering in x-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumy', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumy_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j, k)
                        do k = 1, size(ws, 1)
                            do j = 1, size(ws, 3)
                                ws(k, :, j) = fourier_filt(ws(k, :, j), dy, wavenums, wamps)
                            end do
                        end do
                        !$omp end parallel do
                    end if

                    ! Wavenumber-domain filtering in z-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumz', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumz_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j, k)
                        do k = 1, size(ws, 2)
                            do j = 1, size(ws, 3)
                                ws(:, k, j) = fourier_filt(ws(:, k, j), dz, wavenums, wamps)
                            end do
                        end do
                        !$omp end parallel do
                    end if

                case ('mask')
                    ! Masking
                    call readpar_string(file_parameter, 'dir_shot_'//tidy(name)//'_mask', dir_mask, '')
                    if (dir_mask == '') then
                        call readpar_string(file_parameter, 'shot_'//tidy(name)//'_mask', file_mask, '')
                    else
                        file_mask = tidy(dir_mask)//'/'//tidy(shot_prefix)//'_mask.bin'
                    end if
                    call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    call alloc_array(w_mask, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                        source=w_mask(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                    ws = ws*w_mask

                case ('taper')
                    ! Tapering
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taperx', shot_w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_tapery', shot_w_tapery, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taperz', shot_w_taperz, [0.0, 0.0])
                    if (size(shot_w_taperx) == 1) then
                        call alloc_array(shot_w_taperx, [1, 2], source=[shot_w_taperx(1), shot_w_taperx(1)])
                    end if
                    if (size(shot_w_tapery) == 1) then
                        call alloc_array(shot_w_tapery, [1, 2], source=[shot_w_tapery(1), shot_w_tapery(1)])
                    end if
                    if (size(shot_w_taperz) == 1) then
                        call alloc_array(shot_w_taperz, [1, 2], source=[shot_w_taperz(1), shot_w_taperz(1)])
                    end if
                    ws = taper(ws, nint([shot_w_taperz/dz, shot_w_tapery/dy, shot_w_taperx/dx]), &
                        ['blackman', 'blackman', 'blackman', 'blackman', 'blackman', 'blackman'])

                case ('conemute')
                    call readpar_string(file_parameter, 'shot_'//tidy(name)//'_conemuteorigin', conemuteorigin, 'source')
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemuter', conemuter, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutez', conemutez, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutepower', conemutepower, 2.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutetaper', conemutetaper, 10*mean([dx, dy, dz]))
                    if (conemuter >= 0 .and. conemutez >= 0) then

                        conemutepower = max(1.0, conemutepower)

                        select case (conemuteorigin)
                            case ('source')
                                srcx = mean(gmtr(ishot)%srcr(:)%x - shot_xbeg)
                                srcy = mean(gmtr(ishot)%srcr(:)%y - shot_ybeg)
                                srcz = mean(gmtr(ishot)%srcr(:)%z)
                            case ('receiver')
                                srcx = mean(gmtr(ishot)%recr(:)%x - shot_xbeg)
                                srcy = mean(gmtr(ishot)%recr(:)%y - shot_ybeg)
                                srcz = minval(gmtr(ishot)%recr(:)%z)
                        end select

                        do k = 1, shot_nz

                            depth = (k - 1)*dz + shot_zbeg

                            if (depth < srcz) then
                                ! When the depth is shallower than source then set zero
                                ws(k, :, :) = 0.0
                            else
                                ! Otherwise, compute the circular taper at each depth
                                ! A sufficiently fine sampling interval to capture variation
                                ds = 0.25*min(dy, dx)
                                ! The radius of muting and tapering
                                pr = nint(conemuter*(min(depth - srcz, conemutez)/conemutez)**(1.0/conemutepower)/ds)
                                pt = nint(conemutetaper/ds)
                                ! Taper
                                call alloc_array(st, [1, pr + pt])
                                st = 1.0
                                st = taper(st, [0, pt], ['', 'nuttall'])
                                ! Do ther tapering based on distance of each spatial grid point
                                !$omp parallel do private(ix, iy, dist, tp)
                                do iy = 1, shot_ny
                                    do ix = 1, shot_nx
                                        ! Distance to source location (or its vertical projection)
                                        dist = sqrt(((ix - 1)*dx - srcx)**2 + ((iy - 1)*dy - srcy)**2)
                                        ! Weight
                                        if (dist <= (pr - 1)*ds) then
                                            tp = 1.0
                                        else if (dist >= (pr + pt - 1)*ds) then
                                            tp = 0.0
                                        else
                                            tp = 1 - (dist - (pr - 1)*ds)/(pt*ds)
                                        end if
                                        ws(k, iy, ix) = ws(k, iy, ix)*tp
                                    end do
                                end do
                                !$omp end parallel do
                            end if

                        end do

                    end if

            end select

            if (process_shot_grad_w(i) /= '') then
                call warn(' Model value range = '//num2str(minval(ws), '(es)')//', '//num2str(maxval(ws), '(es)'))
                call warn(date_time_compact()//' shot '//num2str(gmtr(ishot)%id)//' '//tidy(name) &
                    //' processing ('//tidy(process_shot_grad_w(i))//') finished. ')
            end if

        end do

        ! Merge image
        w(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend) = w(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend) + ws

        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)//' '//tidy(name)//' is merged. ')

    end subroutine process_model_single_shot

    !
    !> Processing gradients for a single model parameter
    !
    !> Smoothing gradients can remove high-wavenumber noises
    !> generated during gradient calculations. This functiionality is
    !> not as fancy as anisotropic diffusion denoising, but
    !> is much more less computational expensive.
    !
    subroutine process_model_single_parameter(w, name, param_name)

        real, dimension(:, :, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: param_name

        integer :: i, j, k
        real, allocatable, dimension(:) :: w_taperx, w_tapery, w_taperz
        real :: w_movingbalx, w_movingbaly, w_movingbalz
        real :: w_medianfiltx, w_medianfilty, w_medianfiltz
        real :: w_smoothx, w_smoothy, w_smoothz
        type(andf_param) :: param
        real, allocatable, dimension(:, :, :) :: andfaux, andfcoh
        character(len=32), allocatable, dimension(:) :: process_grad_w
        real :: w_signed_power, w_scalar
        real, allocatable, dimension(:, :, :) :: wt, w_mask
        integer :: wrx, wry, wrz
        real :: w_rmsbalx, w_rmsbaly, w_rmsbalz
        integer, allocatable, dimension(:) :: update_iter
        character(len=1024) :: file_mask

        call readpar_nstring(file_parameter, 'process_'//tidy(name), process_grad_w, [''])

        if (present(param_name)) then
            call readpar_nint(file_parameter, tidy(param_name)//'_update_iter', update_iter, [1, niter_max])
            if (size(update_iter) == 1) then
                update_iter = [update_iter(1), niter_max]
            end if
            if (iter < update_iter(1) .or. iter > update_iter(2)) then
                w = 0
                return
            end if
        end if

        ! Process gradient
        do i = 1, size(process_grad_w)

            select case (process_grad_w(i))

                case ('scale')
                    call readpar_xfloat(file_parameter, tidy(name)//'_scale', w_scalar, 1.0, iter*1.0)
                    w = w*w_scalar

                case ('signed_power')
                    call readpar_xfloat(file_parameter, tidy(name)//'_signed_power', w_signed_power, 0.5, iter*1.0)
                    w = sign(1.0, w)*(abs(w))**w_signed_power

                case ('movingbal')
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbalx', w_movingbalx, 3*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbaly', w_movingbaly, 3*dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbalz', w_movingbalz, 3*dz, iter*1.0)
                    w = balance_filt(w, nint([0.5*w_movingbalz/dz, 0.5*w_movingbaly/dy, 0.5*w_movingbalx/dx]), 0.01)

                case ('taper')
                    call readpar_nfloat(file_parameter, tidy(name)//'_taperx', w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_tapery', w_tapery, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_taperz', w_taperz, [0.0, 0.0])
                    if (size(w_taperx) == 1) then
                        call alloc_array(w_taperx, [1, 2], source=[w_taperx(1), w_taperx(1)])
                    end if
                    if (size(w_tapery) == 1) then
                        call alloc_array(w_tapery, [1, 2], source=[w_tapery(1), w_tapery(1)])
                    end if
                    if (size(w_taperz) == 1) then
                        call alloc_array(w_taperz, [1, 2], source=[w_taperz(1), w_taperz(1)])
                    end if
                    w = taper(w, nint([w_taperz/dz, w_tapery/dy, w_taperx/dx]), &
                        ['blackman', 'blackman', 'blackman', 'blackman', 'blackman', 'blackman'])

                case ('smooth')
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothx', w_smoothx, 3*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothy', w_smoothy, 3*dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothz', w_smoothz, 3*dz, iter*1.0)
                    w = gauss_filt(w, [w_smoothz/dz, w_smoothy/dy, w_smoothx/dx])

                case ('andf')
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothx', param%smooth3, 2*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothy', param%smooth2, 2*dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothz', param%smooth1, 8*dz, iter*1.0)
                    param%smooth3 = param%smooth3/dx
                    param%smooth2 = param%smooth2/dy
                    param%smooth1 = param%smooth1/dz
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_powerm', param%powerm, 1.0, iter*1.0)
                    call readpar_xint(file_parameter, tidy(name)//'_andf_t', param%niter, 5, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_sigma', param%sigma, 6*max(dx, dy, dz), iter*1.0)
                    param%sigma = param%sigma/max(dx, dy, dz)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_alpha', param%lambda1, 1.0e-3, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_gamma', param%lambda3, 1.0, iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_aux', file_andfaux, '', iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_coh', file_andfcoh, '', iter*1.0)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rankx', rank3, 1)
                    call readpar_int(file_parameter, tidy(name)//'_andf_ranky', rank2, 1)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rankz', rank1, 1)
                    if (file_andfaux == '' .and. file_andfcoh == '') then
                        w = andf_filt_mpi(w, param)
                    else if (file_andfaux /= '' .and. file_andfcoh == '') then
                        call prepare_model_single_parameter(andfaux, 'andfaux', file_andfaux, update=.false.)
                        w = andf_filt_mpi(w, param, aux=andfaux)
                    else if (file_andfaux == '' .and. file_andfcoh /= '') then
                        call prepare_model_single_parameter(andfcoh, 'andfcoh', file_andfcoh, update=.false.)
                        w = andf_filt_mpi(w, param, acoh=andfcoh)
                    else if (file_andfaux /= '' .and. file_andfcoh /= '') then
                        call prepare_model_single_parameter(andfaux, 'andfaux', file_andfaux, update=.false.)
                        call prepare_model_single_parameter(andfcoh, 'andfcoh', file_andfcoh, update=.false.)
                        w = andf_filt_mpi(w, param, aux=andfaux, acoh=andfcoh)
                    end if

                case ('medianfilt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfiltx', w_medianfiltx, dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfilty', w_medianfilty, dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfiltz', w_medianfiltz, dz, iter*1.0)
                    w = median_filt(w, nint([w_medianfiltz/dz, w_medianfilty/dy, w_medianfiltx/dx]))

                case ('rmsbal')
                    if (mean(w, 2) /= 0) then
                        w = w/mean(w, 2)
                    end if

                case ('rmsbalx')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbalx', w_rmsbalx, 1.0*dx, iter*1.0)
                    wrx = nint(w_rmsbalx/dx)
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, 0, 0, wrx, wrx])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 3)
                        w(:, :, j) = w(:, :, j)/norm2(wt(:, :, j - (wrx - 1)/2:j + (wrx - 1)/2))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rmsbaly')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbaly', w_rmsbaly, 1.0*dy, iter*1.0)
                    wry = nint(w_rmsbaly/dy)
                    if (mod(wry, 2) == 0) then
                        wry = wry - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wry, wry, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 2)
                        w(:, j, :) = w(:, j, :)/norm2(wt(:, j - (wry - 1)/2:j + (wry - 1)/2, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rmsbalxy')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbalx', w_rmsbalx, 1.0*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbaly', w_rmsbaly, 1.0*dy, iter*1.0)
                    wrx = nint(w_rmsbalx/dx)
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wry = nint(w_rmsbaly/dy)
                    if (mod(wry, 2) == 0) then
                        wry = wry - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wry, wry, wrx, wrx])
                    !$omp parallel do private(j, k)
                    do k = 1, size(w, 3)
                        do j = 1, size(w, 2)
                            w(:, j, k) = w(:, j, k)/norm2(wt(:, j - (wry - 1)/2:j + (wry - 1)/2, k - (wrx - 1)/2:k + (wrx - 1)/2))
                        end do
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rmsbalz')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbalz', w_rmsbalz, 1.0*dz, iter*1.0)
                    wrz = nint(w_rmsbalz/dz)
                    if (mod(wrz, 2) == 0) then
                        wrz = wrz - 1
                    end if
                    wt = w
                    call pad_array(wt, [wrz, wrz, 0, 0, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 1)
                        w(j, :, :) = w(j, :, :)/norm2(wt(j - (wrz - 1)/2:j + (wrz - 1)/2, :, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('mask')
                    call readpar_xstring(file_parameter, tidy(name)//'_mask', file_mask, file_mask, iter*1.0)
                    call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    w = w*w_mask

            end select

            if (rankid == 0) then
                if (process_grad_w(i) /= '') then
                    call warn(' Model value range = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
                    call warn(date_time_compact()//' '//tidy(name)//' processing ('//tidy(process_grad_w(i))//') finished. ')
                end if
            end if

        end do

    end subroutine process_model_single_parameter

end module gradient
