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

    real :: refltaper, reflwindow, reflmin, reflmax
    real, allocatable, dimension(:, :, :) :: refl_all

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
                    model_grad(i)%array = zeros(nz, nx)
                case ('sx', 'sz', 'st0')
                    model_grad(i)%array = zeros(nr_virtual, 1)
            end select
        end do

        if (yn_update_reflector) then
            refl_all = zeros(nz, nx, nrefl)
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

            if (.not. any(model_m(i)%name == ['sx', 'sz', 'st0'])) then

                if (model_m(i)%name /= 'refl') then

                    if (uniform_processing) then
                        call process_model_single_parameter(model_grad(i)%array, 'grad', param_name=model_m(i)%name)
                    else
                        call process_model_single_parameter(model_grad(i)%array, 'grad_'//tidy(model_name(i)), param_name=model_m(i)%name)
                    end if

                else

                    do l = 1, nrefl
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_taper_length', refltaper, 10*dx, iter*1.0)
                        call assert(refltaper < nx*dx/3.0, ' <compute_image_shots> Error: taper is too long. ')
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_smooth_window', reflwindow, 0.2*(nx*dx - 2*refltaper), iter*1.0)
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_min', reflmin, oz, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_max', reflmax, oz + (nz - 1)*dz, iter*1.0)
                        call smooth_reflector(refl_all(:, :, l))
                        model_grad(i)%array = model_grad(i)%array - refl_all(:, :, l)*l
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

        real, dimension(:, :), intent(inout) :: img

        integer :: tplen
        real :: window
        real, allocatable, dimension(:) :: x, y, xx, yy
        integer :: n1, n2, i, j
        integer :: imin
        real, allocatable, dimension(:, :) :: v

        n1 = size(img, 1)
        n2 = size(img, 2)
        tplen = nint(refltaper/dx)
        window = reflwindow/dx/n2

        v = img

        x = regspace(0.0, 1.0, n2 - 1.0)
        y = zeros(n2)

        do j = 1, n2
            imin = as_scalar(maxloc(v(:, j)))
            do i = 1, n1
                if (i == imin) then
                    v(imin, j) = 1.0
                    y(j) = imin
                else
                    v(i, j) = 0.0
                end if
            end do
        end do
        y = clip(y, nint((reflmin - oz)/dz) + 1.0, nint((reflmax - oz)/dz) + 1.0)

        xx = x(tplen + 1:n2 - tplen)
        yy = y(tplen + 1:n2 - tplen)

        ! The window should be a fraction of the total length of xx
        yy = lowess_filt(xx, yy, xx, window=window, robust=.true.)
        y = ginterp(xx, yy, x, 'linear')

        y = clip(y, nint((reflmin - oz)/dz) + 1.0, nint((reflmax - oz)/dz) + 1.0)

        v = 0
        do j = 1, n2
            v(nint(y(j)), j) = 1.0
        end do

        img = v

    end subroutine smooth_reflector

    !
    !> Process 2D model for each shot
    !
    subroutine process_model_single_shot(ishot, w, ws, name)

        integer, intent(in) :: ishot
        real, dimension(:, :), intent(inout) :: w, ws
        character(len=*), intent(in) :: name

        integer :: i, j
        character(len=1024) :: shot_prefix, dir_mask, file_mask
        real :: shot_w_movingbalx, shot_w_movingbalz
        real :: shot_w_medianfiltx, shot_w_medianfiltz
        real, allocatable, dimension(:) :: shot_w_taperx, shot_w_taperz
        real :: shot_w_smoothx, shot_w_smoothz
        real :: shot_w_adpmutex, shot_w_adpmutez
        real :: recmin, recmax, srcmin, srcmax
        integer :: lb, ub, px, l
        real, allocatable, dimension(:, :) :: w_mask
        real, allocatable, dimension(:) :: st, tp
        real, allocatable, dimension(:) :: wavenums, wamps, fkdips, fkdipamps
        real :: conemutex, conemutez, conemutepower, conemutetaper
        real :: srcx, srcz, depth
        character(len=32), allocatable, dimension(:) :: process_shot_grad_w

        call readpar_nstring(file_parameter, 'process_shot_'//tidy(name), process_shot_grad_w, [''])

        do i = 1, size(process_shot_grad_w)

            select case (process_shot_grad_w(i))

                case ('smooth')
                    ! Gaussian smooth
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothx', shot_w_smoothx, 3.0*dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothz', shot_w_smoothz, 3.0*dz)
                    ws = gauss_filt(ws, [shot_w_smoothz/dz, shot_w_smoothx/dx])

                case ('maxbal')
                    if (maxval(ws) /= 0) then
                        ws = ws/maxval(ws)
                    end if

                case ('rmsbal')
                    ! Normalize with shot image energy
                    ws = ws/mean(ws, 2)

                case ('movingbal')
                    ! Moving balance
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbalx', shot_w_movingbalx, 3*dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbalz', shot_w_movingbalz, 3*dz)
                    ws = balance_filt(ws, nint([0.5*shot_w_movingbalz/dz, 0.5*shot_w_movingbalx/dx]), 0.01)

                case ('medianfilt')
                    ! Median filtering
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfiltx', shot_w_medianfiltx, dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfiltz', shot_w_medianfiltz, dz)
                    ws = median_filt(ws, nint([shot_w_medianfiltz/dz, shot_w_medianfiltx/dx]))

                case ('dipfilt')
                    ! Dip filtering
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzx', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzx_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        ws = dip_filt(ws, [1.0, dx/dz], fkdips, fkdipamps)
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
                        !$omp parallel do private(j)
                        do j = 1, size(ws, 1)
                            ws(j, :) = fourier_filt(ws(j, :), dx, wavenums, wamps)
                        end do
                        !$omp end parallel do
                    end if

                    ! Wavenumber-domain filtering in z-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumz', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumz_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j)
                        do j = 1, size(ws, 2)
                            ws(:, j) = fourier_filt(ws(:, j), dz, wavenums, wamps)
                        end do
                        !$omp end parallel do
                    end if

                case ('taper')
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taperx', shot_w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taperz', shot_w_taperz, [0.0, 0.0])
                    if (size(shot_w_taperx) == 1) then
                        call alloc_array(shot_w_taperx, [1, 2], source=[shot_w_taperx(1), shot_w_taperx(1)])
                    end if
                    if (size(shot_w_taperz) == 1) then
                        call alloc_array(shot_w_taperz, [1, 2], source=[shot_w_taperz(1), shot_w_taperz(1)])
                    end if
                    ws = taper(ws, nint([shot_w_taperz/dz, shot_w_taperx/dx]), &
                        ['blackman', 'blackman', 'blackman', 'blackman'])

                case ('mask')
                    ! Masking
                    call readpar_string(file_parameter, 'dir_shot_'//tidy(name)//'_mask', dir_mask, '')
                    if (dir_mask == '') then
                        call readpar_string(file_parameter, 'shot_'//tidy(name)//'_mask', file_mask, '')
                    else
                        file_mask = tidy(dir_mask)//'/'//tidy(shot_prefix)//'_mask.bin'
                    end if
                    call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    call alloc_array(w_mask, [1, shot_nz, 1, shot_nx], &
                        source=w_mask(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                    ws = ws*w_mask

                case ('adpmute')
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_adpmutex', shot_w_adpmutex, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_adpmutez', shot_w_adpmutez, -1.0)
                    if (shot_w_adpmutex >= 0) then
                        ! find source-receiver widest possible range
                        recmin = minval(gmtr(ishot)%recr(:)%x - shot_xbeg)
                        recmax = maxval(gmtr(ishot)%recr(:)%x - shot_xbeg)
                        srcmin = minval(gmtr(ishot)%srcr(:)%x - shot_xbeg)
                        srcmax = maxval(gmtr(ishot)%srcr(:)%x - shot_xbeg)
                        ! ... and their integer grid point positions in the computed image
                        lb = nint(min(recmin, srcmin)/dx + 1)
                        ub = nint(max(recmax, srcmax)/dx + 1)
                        ! the taper length
                        px = nint(shot_w_adpmutex/dx)
                        ! create the taper
                        call alloc_array(st, [lb, ub], pad=px)
                        st = 1.0
                        st = taper(st, [px, px], ['blackman', 'blackman'])
                        ! put the taper in the whole x range, which can be longer than the taper
                        call alloc_array(tp, [1, shot_nx])
                        do l = lb - px, ub + px
                            if (l >= 1 .and. l <= shot_nx) then
                                tp(l) = st(l)
                            end if
                        end do
                        ! now the taper has the same length with the image, and do the tapering
                        ! the resulting tapered image is now restricted to the region cropped by
                        ! the largest possible source-receiver offset
                        do l = 1, shot_nz
                            ws(l, :) = ws(l, :)*tp
                        end do
                    end if

                case ('conemute')
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutex', conemutex, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutez', conemutez, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutepower', conemutepower, 2.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutetaper', conemutetaper, 10*dx)
                    if (conemutex >= 0 .and. conemutez >= 0) then

                        conemutepower = max(1.0, conemutepower)
                        srcx = mean(gmtr(ishot)%srcr(:)%x - shot_xbeg)
                        srcz = mean(gmtr(ishot)%srcr(:)%z)
                        call alloc_array(tp, [1, shot_nx])

                        do j = 1, shot_nz

                            depth = (j - 1)*dz + shot_zbeg

                            if (depth < srcz) then
                                ! When the depth is shallower than source then set zero
                                ws(j, :) = 0.0
                            else
                                ! Depth is deeper than the source
                                tp = 0.0
                                ! lower and upper spatial range
                                px = conemutex*(min(depth - srcz, conemutez)/conemutez)**(1.0/conemutepower)
                                lb = nint((srcx - px)/dx + 1)
                                ub = nint((srcx + px)/dx + 1)
                                ! the taper length
                                px = nint(conemutetaper/dx)
                                ! create the taper
                                call alloc_array(st, [lb, ub], pad=px)
                                st = 1.0
                                st = taper(st, [px, px], ['blackman', 'blackman'])
                                ! put the taper in the whole x range, which can be longer than the taper
                                do l = lb - px, ub + px
                                    if (l >= 1 .and. l <= shot_nx) then
                                        tp(l) = st(l)
                                    end if
                                end do
                                ws(j, :) = ws(j, :)*tp
                            end if

                        end do

                    end if

            end select

            if (process_shot_grad_w(i) /= '') then
                call warn(' Model value range = '//num2str(minval(ws), '(es)')//', '//num2str(maxval(ws), '(es)'))
                call warn(date_time_compact()//' shot '//num2str(gmtr(ishot)%id)//' '//tidy(name) &
                    //' processing ('//tidy(process_shot_grad_w(i))//') is done. ')
            end if

        end do

        ! Merge w
        w(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend) = w(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend) + ws

        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)//' '//tidy(name)//' is merged. ')

    end subroutine process_model_single_shot

    !
    !> Processing model for a single parameter
    !
    !> Smoothing gradients can remove high-wavenumber noises
    !> generated during gradient calculations. This functiionality is
    !> not as fancy as anisotropic diffusion denoising, but
    !> is much more less computational expensive.
    !
    subroutine process_model_single_parameter(w, name, param_name)

        real, dimension(:, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: param_name

        integer :: i, j
        real, allocatable, dimension(:) :: w_taperx, w_taperz
        real :: w_movingbalx, w_movingbalz
        real :: w_medianfiltx, w_medianfiltz
        real :: w_smoothx, w_smoothz
        type(andf_param) :: param
        real, allocatable, dimension(:, :) :: andfaux, andfcoh, w_mask
        character(len=32), allocatable, dimension(:) :: process_grad_w
        real :: w_signed_power
        real :: w_scalar
        real :: w_rmsbalx, w_rmsbalz
        integer :: wrx, wrz
        real, allocatable, dimension(:, :) :: wt
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

        ! Process wient
        do i = 1, size(process_grad_w)

            if (rankid == 0) then
                call warn(' Model value range before processing = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
            end if

            select case (process_grad_w(i))

                case ('scale')
                    call readpar_xfloat(file_parameter, tidy(name)//'_scale', w_scalar, 1.0, iter*1.0)
                    w = w*w_scalar

                case ('signed_power')
                    call readpar_xfloat(file_parameter, tidy(name)//'_signed_power', w_signed_power, 0.5, iter*1.0)
                    w = sign(1.0, w)*(abs(w))**w_signed_power

                case ('maxbal')
                    if (maxval(w) /= 0) then
                        w = w/maxval(w)
                    end if

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
                    call pad_array(wt, [0, 0, wrx, wrx])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 2)
                        w(:, j) = w(:, j)/norm2(wt(:, j - (wrx - 1)/2:j + (wrx - 1)/2))
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
                    call pad_array(wt, [wrz, wrz, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 1)
                        w(j, :) = w(j, :)/norm2(wt(j - (wrz - 1)/2:j + (wrz - 1)/2, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('movingbal')
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbalx', w_movingbalx, 6*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbalz', w_movingbalz, 6*dz, iter*1.0)
                    w = balance_filt(w, nint([0.5*w_movingbalz/dz, 0.5*w_movingbalx/dx]), 0.01)

                case ('taper')
                    call readpar_nfloat(file_parameter, tidy(name)//'_taperx', w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_taperz', w_taperz, [0.0, 0.0])
                    if (size(w_taperx) == 1) then
                        call alloc_array(w_taperx, [1, 2], source=[w_taperx(1), w_taperx(1)])
                    end if
                    if (size(w_taperz) == 1) then
                        call alloc_array(w_taperz, [1, 2], source=[w_taperz(1), w_taperz(1)])
                    end if
                    w = taper(w, nint([w_taperz/dz, w_taperx/dx]), ['blackman', 'blackman', 'blackman', 'blackman'])

                case ('smooth')
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothx', w_smoothx, 3*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothz', w_smoothz, 3*dz, iter*1.0)
                    w = gauss_filt(w, [w_smoothz/dz, w_smoothx/dx])

                case ('andf')
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothx', param%smooth2, 2.0*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothz', param%smooth1, 8.0*dz, iter*1.0)
                    param%smooth2 = param%smooth2/dx
                    param%smooth1 = param%smooth1/dz
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_powerm', param%powerm, 1.0, iter*1.0)
                    call readpar_xint(file_parameter, tidy(name)//'_andf_t', param%niter, 5, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_sigma', param%sigma, 6.0*max(dx, dz), iter*1.0)
                    param%sigma = param%sigma/max(dx, dz)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_alpha', param%lambda1, 1.0e-3, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_aux', file_andfaux, '', iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_coh', file_andfcoh, '', iter*1.0)
                    if (file_andfaux == '' .and. file_andfcoh == '') then
                        w = andf_filt(w, param)
                    else if (file_andfaux /= '' .and. file_andfcoh == '') then
                        call prepare_model_single_parameter(andfaux, 'andfaux', file_andfaux, update=.false.)
                        w = andf_filt(w, param, aux=andfaux)
                    else if (file_andfaux == '' .and. file_andfcoh /= '') then
                        call prepare_model_single_parameter(andfcoh, 'andfcoh', file_andfcoh, update=.false.)
                        w = andf_filt(w, param, acoh=andfcoh)
                    else if (file_andfaux /= '' .and. file_andfcoh /= '') then
                        call prepare_model_single_parameter(andfaux, 'andfaux', file_andfaux, update=.false.)
                        call prepare_model_single_parameter(andfcoh, 'andfcoh', file_andfcoh, update=.false.)
                        w = andf_filt(w, param, aux=andfaux, acoh=andfcoh)
                    end if

                case ('medianfilt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfiltx', w_medianfiltx, dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfiltz', w_medianfiltz, dz, iter*1.0)
                    w = median_filt(w, nint([w_medianfiltz/dz, w_medianfiltx/dx]))

                case ('mask')
                    call readpar_xstring(file_parameter, tidy(name)//'_mask', file_mask, file_mask, iter*1.0)
                    call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    w = w*w_mask

            end select

            if (rankid == 0) then
                call warn(' Model value range after processing = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
                if (process_grad_w(i) /= '') then
                    call warn(date_time_compact()//' '//tidy(name)//' processing ('//tidy(process_grad_w(i))//') finished. ')
                end if
            end if

        end do

    end subroutine process_model_single_parameter

end module gradient
