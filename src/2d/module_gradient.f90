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


module gradient

    use libflit
    use parameters
    use vars
    use utility

    implicit none

    real :: refltaper, reflwindow, reflmin, reflmax
    ! Gauss-Newton step of the reflector depths (refine_reflector, and the
    ! joint update with step_max_refl > 0): row l of model_grad(refl)
    ! accumulates the numerator of reflector l by column, and refl_den the
    ! denominator (see reflector_direction)
    real, allocatable, dimension(:, :) :: refl_den

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
        module subroutine compute_reflector_step_shots
        end subroutine
    end interface

    private

    public :: zero_gradient
    public :: compute_gradient_shots
    public :: process_gradient
    public :: output_gradient
    public :: compute_image_shots
    public :: move_reflector
    public :: compute_reflector_step_shots
    public :: zero_reflector_step
    public :: reflector_step

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

        if (reflector_joint()) then
            refl_den = zeros(nz, nx)
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

        integer :: i
        real, allocatable, dimension(:, :) :: dzr

        do i = 1, nmodel

            if (.not. any(model_m(i)%name == ['sx', 'sz', 'st0'])) then

                if (model_m(i)%name /= 'refl') then

                    if (uniform_processing) then
                        call process_model_single_parameter(model_grad(i)%array, 'grad', param_name=model_m(i)%name)
                    else
                        call process_model_single_parameter(model_grad(i)%array, 'grad_'//tidy(model_name(i)), param_name=model_m(i)%name)
                    end if

                else if (reflector_joint()) then

                    ! Joint update of the reflector depths: row l of the gradient
                    ! holds minus the Gauss-Newton depth change of reflector l by
                    ! column, so that the search direction is that change
                    call reflector_direction(model_grad(i)%array, refl_den, dzr)
                    model_grad(i)%array = -dzr

                else

                    ! The reflector model is set by imaging and refinement
                    ! (update_reflector_image), not by a gradient step
                    model_grad(i)%array = 0.0

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

            call warn(date_time_compact()//' >>>>>>>>>> Gradient saved. ')

        end if

    end subroutine output_gradient

    !
    !> Build the model of reflector l from the per-shot images zs(nx, ns): the
    !> fractional 1-based z index of the reflector in each column for each
    !> shot, or 0 where the shot has no crossing. Every isochron lies above
    !> the reflector except at its specular point, so a crossing that is not
    !> specular (the isochron of an end receiver beyond the lit part of the
    !> reflector, or the small negative region around a receiver) lies above
    !> the reflector, and the reflector is the deepest crossing over the
    !> shots. A column therefore takes the median of the shots whose crossing
    !> is within half a grid cell of the deepest one. The columns without a
    !> crossing are interpolated from the others. The curve is smoothed when
    !> reflector_<l>_smooth_window > 0, clipped to the depth range of the
    !> reflector, and stored at the nearest grid point of each column as
    !> l + dfrac, where dfrac in [-0.5, 0.5) is the offset of the reflector
    !> from the grid point in grid cells.
    !
    subroutine build_reflector(zs, l, refl)

        real, dimension(:, :), intent(in) :: zs
        integer, intent(in) :: l
        real, dimension(:, :), intent(inout) :: refl

        real, allocatable, dimension(:) :: zbar, x, xx, yy
        logical, allocatable, dimension(:) :: has
        integer :: i, j, k, n, j1, j2, tplen
        real :: dfrac, window, ztol, zmax

        ! The crossings within half a grid cell of the deepest one
        ztol = 0.5
        n = size(zs, 1)
        zbar = zeros(n)
        has = falses(n)
        do j = 1, n
            if (any(zs(j, :) > 0)) then
                zmax = maxval(zs(j, :))
                zbar(j) = median(pack(zs(j, :), zs(j, :) >= zmax - ztol))
                has(j) = .true.
            end if
        end do
        if (rankid == 0) then
            call warn(date_time_compact()//' Reflector '//num2str(l)//' is imaged in ' &
                //num2str(count(has))//' of '//num2str(n)//' columns. ')
        end if
        if (.not. any(has)) then
            return
        end if

        ! Columns without an image: linear interpolation between imaged
        ! columns, and the nearest imaged column beyond the ends
        j1 = 1
        do while (.not. has(j1))
            j1 = j1 + 1
        end do
        j2 = n
        do while (.not. has(j2))
            j2 = j2 - 1
        end do
        zbar(1:j1 - 1) = zbar(j1)
        zbar(j2 + 1:n) = zbar(j2)
        j = j1
        do while (j < j2)
            if (has(j + 1)) then
                j = j + 1
            else
                i = j + 1
                do while (.not. has(i))
                    i = i + 1
                end do
                do k = j + 1, i - 1
                    zbar(k) = zbar(j) + (zbar(i) - zbar(j))*real(k - j)/real(i - j)
                end do
                j = i
            end if
        end do

        ! Optional LOWESS smoothing, excluding a taper at each end
        if (reflwindow > 0) then
            tplen = min(max(nint(refltaper/dx), 0), (n - 2)/2)
            window = reflwindow/dx/n
            x = regspace(0.0, 1.0, n - 1.0)
            xx = x(tplen + 1:n - tplen)
            yy = zbar(tplen + 1:n - tplen)
            yy = lowess_filt(xx, yy, xx, window=window, robust=.true.)
            zbar = ginterp(xx, yy, x, 'linear')
        end if

        ! Depth range, and the label-plus-offset model
        zbar = clip(zbar, (reflmin - oz)/dz + 1.0, (reflmax - oz)/dz + 1.0)
        zbar = clip(zbar, 1.0, nz*1.0)
        do j = 1, n
            i = nint(zbar(j))
            dfrac = zbar(j) - i
            i = min(max(i, 1), nz)
            refl(i, j) = l + dfrac
        end do

    end subroutine build_reflector

    !
    !> The fractional 1-based z index of reflector l in each column of a
    !> label-plus-offset reflector model, or 0 where the column has no point
    !> of the reflector
    !
    subroutine reflector_depth_index(refl, l, zidx)

        real, dimension(:, :), intent(in) :: refl
        integer, intent(in) :: l
        real, dimension(:), intent(out) :: zidx

        integer :: i, j

        zidx = 0.0
        do j = 1, size(refl, 2)
            do i = 1, size(refl, 1)
                if (nint(refl(i, j)) == l) then
                    zidx(j) = i + (refl(i, j) - l)
                    exit
                end if
            end do
        end do

    end subroutine reflector_depth_index

    !
    !> Zero the accumulators of the Gauss-Newton step of the reflector depths
    !
    subroutine zero_reflector_step(irefl)

        integer, intent(in) :: irefl

        model_grad(irefl)%array = zeros(nz, nx)
        refl_den = zeros(nz, nx)

    end subroutine zero_reflector_step

    !
    !> The Gauss-Newton step of the reflector depths from the accumulators of
    !> compute_reflector_step_shots (see reflector_direction), capped at
    !> step_max: row l of dzr holds the depth change of reflector l in each
    !> column, in length units.
    !
    subroutine reflector_step(irefl, step_max, dzr)

        integer, intent(in) :: irefl
        real, intent(in) :: step_max
        real, allocatable, dimension(:, :), intent(out) :: dzr

        call reflector_direction(model_grad(irefl)%array, refl_den, dzr)
        dzr = clip(dzr, -step_max, step_max)

    end subroutine reflector_step

    !
    !> The Gauss-Newton depth change of the reflectors from the accumulated
    !> numerator num and denominator den: row l of dzr holds the change of
    !> reflector l in each column, in length units, without a cap. num and den
    !> are smoothed along x with a Gaussian of length refl_refine_smooth
    !> before the division, which gives the least-squares depth change that is
    !> locally constant over the smoothing length. The damping is 1% of the
    !> largest smoothed den of each reflector, so weakly lit columns move
    !> little.
    !
    subroutine reflector_direction(num, den, dzr)

        real, dimension(:, :), intent(in) :: num, den
        real, allocatable, dimension(:, :), intent(out) :: dzr

        real, allocatable, dimension(:, :) :: sn, sd
        real :: smooth
        integer :: l

        call readpar_xfloat(file_parameter, 'refl_refine_smooth', smooth, 5*dx, iter*1.0)
        sn = num
        sd = den
        if (smooth > 0) then
            sn = gauss_filt(sn, [0.0, smooth/dx])
            sd = gauss_filt(sd, [0.0, smooth/dx])
        end if
        dzr = zeros(nz, nx)
        do l = 1, nrefl
            where (sd(l, :) > 0)
                dzr(l, :) = -sn(l, :)/(sd(l, :) + 0.01*maxval(sd(l, :)))
            end where
        end do

    end subroutine reflector_direction

    !
    !> Whether the reflector depths take part in the joint update with the
    !> velocities, which step_max_refl > 0 selects in the current iteration
    !
    logical function reflector_joint()
        real :: step
        reflector_joint = .false.
        if (yn_update_reflector) then
            call readpar_xfloat(file_parameter, 'step_max_refl', step, 0.0, iter*1.0)
            reflector_joint = step > 0
        end if
    end function reflector_joint

    !
    !> Move the reflectors of a label-plus-offset model: row l of dzr holds the
    !> depth change of reflector l in each column, in length units. The depths
    !> stay within reflector_<l>_depth_min and _max.
    !
    subroutine move_reflector(refl, dzr)

        real, dimension(:, :), intent(inout) :: refl
        real, dimension(:, :), intent(in) :: dzr

        real, allocatable, dimension(:, :) :: zidx
        real :: zlo, zhi, dfrac
        integer :: i, j, l, n

        n = size(refl, 2)
        zidx = zeros(n, nrefl)
        do l = 1, nrefl
            call reflector_depth_index(refl, l, zidx(:, l))
        end do
        refl = 0.0
        do l = 1, nrefl
            call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_min', zlo, oz, iter*1.0)
            call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_max', zhi, oz + (nz - 1)*dz, iter*1.0)
            zlo = (zlo - oz)/dz + 1.0
            zhi = (zhi - oz)/dz + 1.0
            do j = 1, n
                if (zidx(j, l) > 0) then
                    zidx(j, l) = zidx(j, l) + dzr(l, j)/dz
                    zidx(j, l) = min(max(zidx(j, l), max(zlo, 1.0)), min(zhi, nz*1.0))
                    i = nint(zidx(j, l))
                    dfrac = zidx(j, l) - i
                    i = min(max(i, 1), nz)
                    refl(i, j) = l + dfrac
                end if
            end do
        end do

    end subroutine move_reflector

    !
    !> Process 2D model for each shot
    !
    subroutine process_model_single_shot(ishot, w, ws, name)

        integer, intent(in) :: ishot
        real, dimension(:, :), intent(inout) :: w, ws
        character(len=*), intent(in) :: name

        integer :: i, j
        character(len=1024) :: dir_mask, file_mask
        real :: shot_w_movingbalx, shot_w_movingbalz
        real :: shot_w_medianfiltx, shot_w_medianfiltz
        real, allocatable, dimension(:) :: shot_w_taperx, shot_w_taperz
        real :: shot_w_smoothx, shot_w_smoothz
        integer :: lb, ub, px, l
        real, allocatable, dimension(:, :) :: w_mask
        real, allocatable, dimension(:) :: st, tp
        real, allocatable, dimension(:) :: wavenums, wamps, fkdips, fkdipamps
        real :: conemutex, conemutez, conemutepower, conemutetaper
        real :: srcx, srcz, depth
        character(len=32), allocatable, dimension(:) :: process_shot_grad_w
        integer, dimension(1:4) :: ntaper

        call readpar_nstring(file_parameter, 'process_shot_'//tidy(name), process_shot_grad_w, [''])

        ws = return_normal(ws)

        do i = 1, size(process_shot_grad_w)

            select case (process_shot_grad_w(i))

                case ('smooth')
                    ! Gaussian smooth
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothx', shot_w_smoothx, 3.0*dx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothz', shot_w_smoothz, 3.0*dz)
                    ws = gauss_filt(ws, [shot_w_smoothz/dz, shot_w_smoothx/dx])

                case ('maxbal')
                    if (maxval(abs(ws)) /= 0) then
                        ws = ws/maxval(abs(ws))
                    end if

                case ('rmsbal')
                    ! Normalize with shot image energy
                    if (mean(ws, 2) /= 0) then
                        ws = ws/mean(ws, 2)
                    end if

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
                        file_mask = tidy(dir_mask)//'/shot_'//num2str(gmtr(ishot)%id)//'_mask.bin'
                    end if
                    if (file_mask /= '') then
                        call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    else
                        w_mask = ones(nz, nx)
                    end if 
                    call alloc_array(w_mask, [1, shot_nz, 1, shot_nx], &
                        source=w_mask(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                    ws = ws*w_mask

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
                    //' processing ('//tidy(process_shot_grad_w(i))//') completed. ')
            end if

        end do

        ! Taper the sides of the shot range that lie inside the model, so that
        ! the merged gradient has no sharp edges along the adaptive range;
        ! an axis without adaptive range is never tapered
        if (name /= 'refl') then
            ntaper = [ &
                merge(min(nint(adptaperz/dz), shot_nz/2), 0, yn_adpz .and. shot_nzbeg > 1), &
                merge(min(nint(adptaperz/dz), shot_nz/2), 0, yn_adpz .and. shot_nzend < nz), &
                merge(min(nint(adptaperx/dx), shot_nx/2), 0, yn_adpx .and. shot_nxbeg > 1), &
                merge(min(nint(adptaperx/dx), shot_nx/2), 0, yn_adpx .and. shot_nxend < nx)]
            if (any(ntaper > 0)) then
                ws = taper(ws, ntaper, ['blackman', 'blackman', 'blackman', 'blackman'])
            end if
        end if

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
                    if (maxval(abs(w)) /= 0) then
                        w = w/maxval(abs(w))
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
                    call readpar_xstring(file_parameter, tidy(name)//'_mask', file_mask, '', iter*1.0)
                    if (file_mask /= '') then
                        call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    else
                        w_mask = ones_like(w)
                    end if 
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
