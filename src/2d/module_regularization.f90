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


module regularization

    use libflit
    use parameters
    use vars
    use utility

    implicit none

contains

    !
    !> Solve the Tikhonov regularization in a hard way:
    !> I decompose the overall inversion problem into two subproblems, and then solve
    !> the Tikhonov regularization problem as
    !>             $$||\nabla u||_2^2 = ||u-m^{k+1}||_2^2$$
    !> The first-order optimality condition of the equation is
    !>             $$\nabla^T \nabla u = u - m^{k+1}$$
    !> which can be solved using the Gauss-Seidel method.
    !
    function tikhonov_denoise(m, tikhonov_lambda) result(u)

        real, dimension(:, :) :: m
        real :: tikhonov_lambda
        real, allocatable, dimension(:, :) :: u

        real, allocatable, dimension(:, :) :: pu
        integer :: i, j, iter
        real :: sumu
        integer :: n1, n2, niter
        integer :: i1, i2, j1, j2
        integer :: ii1, ii2, jj1, jj2
        real :: mu

        u = m
        mu = 1.0
        niter = 200
        n1 = size(u, 1)
        n2 = size(u, 2)

        ! TGpV iteration
        do iter = 1, niter

            pu = u

            do j = 1, n2
                do i = 1, n1

                    if (i == 1) then
                        i1 = 0
                        ii1 = i
                    else
                        i1 = 1
                        ii1 = i - 1
                    end if

                    if (i == n1) then
                        i2 = 0
                        ii2 = i
                    else
                        i2 = 1
                        ii2 = i + 1
                    end if

                    if (j == 1) then
                        j1 = 0
                        jj1 = j
                    else
                        j1 = 1
                        jj1 = j - 1
                    end if
                    if (j == n2) then
                        j2 = 0
                        jj2 = j
                    else
                        j2 = 1
                        jj2 = j + 1
                    end if

                    sumu = tikhonov_lambda*( &
                        +i2*u(ii2, j) + i1*u(ii1, j) &
                        + j2*u(i, jj2) + j1*u(i, jj1)) &
                        + mu*m(i, j)

                    u(i, j) = sumu/(mu + (4.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2)) &
                        *tikhonov_lambda)

                end do
            end do

            ! progress
            if (rankid == 0 .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
                call warn(date_time_compact()//'>> Tikhonov iteration '//tidy(num2str(iter, '(i)')) &
                    //' of '//tidy(num2str(niter, '(i)')// &
                    ' relative norm2 diff = '//tidy(num2str(norm2(u - pu), '(es)'))))
            end if

        end do

    end function tikhonov_denoise

    !
    !> Apply regularization to a single model parameter
    !
    subroutine regularize_single_parameter(m, mr, name, default_value, smooth_inverse)

        !> The model to be regularized
        real, dimension(:, :), intent(in) :: m
        !> The regularized model
        real, dimension(:, :), intent(inout) :: mr
        !> Name of the model
        character(len=*), intent(in) :: name
        !> Default value of \(\mu\) in TGpV regularization
        real, intent(in) :: default_value
        !> For smooth regularization, regularize the inverse of the model or not
        logical, intent(in), optional :: smooth_inverse

        integer :: i, n1, n2
        real, allocatable, dimension(:, :) :: mt, maux, mcoh
        character(len=1024) :: file_aux, file_coh
        real :: reg_smoothx, reg_smoothz
        type(andf_param) :: param
        real :: tv_mu
        real :: tv_lambda1
        real :: tv_lambda2
        real :: tv_norm
        integer :: tv_niter
        real :: tikhonov_lambda
        real :: sim_smoothx, sim_smoothz
        real, allocatable, dimension(:, :) :: vp, vs, r
        real :: rmin, rmax

        n1 = size(m, 1)
        n2 = size(m, 2)
        mt = m

        call readpar_int(file_parameter, 'rankx', rank2, floor(nrank**0.5d0))
        call readpar_int(file_parameter, 'rankz', rank1, floor(nrank**0.5d0))

        do i = 1, size(model_regularization_method)

            select case (model_regularization_method(i))

                case ('Tikhonov', 'tikhonov')
                    call readpar_xfloat(file_parameter, 'reg_tikhonov_lambda', tikhonov_lambda, 10.0, iter*1.0)
                    mt = tikhonov_denoise(mt, tikhonov_lambda)
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Tikhonov regularization finished. ')
                    end if

                case ('smooth')
                    call readpar_xfloat(file_parameter, 'reg_smoothx', reg_smoothx, -1.0, iter*1.0)
                    if (reg_smoothx < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smoothx_'//tidy(name), reg_smoothx, 1.0*dx, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smoothz', reg_smoothz, -1.0, iter*1.0)
                    if (reg_smoothz < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smoothz_'//tidy(name), reg_smoothz, 1.0*dz, iter*1.0)
                    end if

                    if (present(smooth_inverse) .and. smooth_inverse) then
                        mt = 1.0/gauss_filt(1.0/mt, [reg_smoothz/dz, reg_smoothx/dx])
                    else
                        mt = gauss_filt(mt, [reg_smoothz/dz, reg_smoothx/dx])
                    end if

                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Smoothing regularization finished. ')
                    end if

                case ('TGpV', 'tgpv')
                    call readpar_xfloat(file_parameter, 'reg_tv_mu_'//tidy(name), tv_mu, default_value, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_lambda1', tv_lambda1, -1.0, iter*1.0)
                    if (tv_lambda1 == -1.0) then
                        call readpar_xfloat(file_parameter, 'reg_tv_lambda1_'//tidy(name), tv_lambda1, 1.0, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_tv_lambda2', tv_lambda2, -1.0, iter*1.0)
                    if (tv_lambda2 == -1.0) then
                        call readpar_xfloat(file_parameter, 'reg_tv_lambda2_'//tidy(name), tv_lambda2, 1.0, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_tv_norm', tv_norm, 0.5, iter*1.0)
                    call readpar_xint(file_parameter, 'reg_tv_niter', tv_niter, 50, iter*1.0)
                    mt = tgpv_filt_mpi(mt, tv_mu, tv_lambda1, tv_lambda2, tv_niter, tv_norm)

                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> TGpV regularization finished. ')
                    end if

                case ('structure')
                    call readpar_xfloat(file_parameter, 'reg_andf_alpha', param%lambda1, 0.001, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_smoothx', param%smooth2, 2.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_smoothz', param%smooth1, 8.0, iter*1.0)
                    call readpar_xint(file_parameter, 'reg_andf_t', param%niter, 10, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_sigma', param%sigma, 10.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_powerm', param%powerm, 4.0, iter*1.0)
                    call readpar_xstring(file_parameter, 'reg_andf_aux', file_aux, '', iter*1.0)
                    call readpar_xstring(file_parameter, 'reg_andf_coh', file_coh, '', iter*1.0)
                    if (file_aux /= '' .and. file_coh == '') then
                        maux = load(file_aux, n1, n2)
                        mt = andf_filt_mpi(mt, param, aux=maux)
                    else if (file_aux == '' .and. file_coh /= '') then
                        mcoh = load(file_coh, n1, n2)
                        mt = andf_filt_mpi(mt, param, acoh=mcoh)
                    else if (file_aux /= '' .and. file_coh /= '') then
                        maux = load(file_aux, n1, n2)
                        mcoh = load(file_coh, n1, n2)
                        mt = andf_filt_mpi(mt, param, aux=maux, acoh=mcoh)
                    else
                        mt = andf_filt_mpi(mt, param)
                    end if
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Structure-oriented regularization finished. ')
                    end if

                case ('vpvs_similarity')
                    if (name == 'vs') then

                        call readpar_xfloat(file_parameter, 'reg_similarity_smoothx', sim_smoothx, 5*dx, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_smoothz', sim_smoothz, 5*dz, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_vpvs_ratio_min', rmin, min_vpvsratio, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_vpvs_ratio_max', rmax, max_vpvsratio, iter*1.0)

                        if (any('vp' == model_name)) then
                            vp = get_meta_array_core(model_m, 'vp')
                        else if (any('vp' == model_name_aux)) then
                            vp = get_meta_array_core(model_aux, 'vp')
                        end if
                        vs = mt

                        mt = ones_like(vs)
                        where (vs == 0)
                            mt = 0
                        end where

                        r = zeros_like(vp)
                        where (mt /= 0)
                            r = vp/vs
                        end where
                        r = clip(r, rmin, rmax)
                        r = median_filt(r)
                        r = gauss_filt(r, [sim_smoothz, sim_smoothx]/[dz, dx])

                        where (mt /= 0)
                            vs = vp/r
                        end where
                        mt = vs

                    end if

            end select

        end do

        call mpibarrier
        mr = mt

        ! Output regularized model
        if (rankid == 0) then
            call output_array(mr, dir_iter_model(iter)//'/reg_'//tidy(name)//'.bin')
        end if

    end subroutine regularize_single_parameter

    !
    !> Apply regularization to model parameters
    !
    subroutine model_regularization

        integer :: i

        do i = 1, nmodel

            if (.not. any(model_m(i)%name == ['sx', 'sy', 'sz', 'st0'])) then

                if (model_m(i)%name == 'vp' .or. model_m(i)%name == 'vs') then
                    call regularize_single_parameter(model_m(i)%array, &
                        model_reg(i)%array, model_name(i), 100.0/(maxval(abs(model_m(i)%array)) + float_tiny), &
                        smooth_inverse=.true.)
                else
                    call regularize_single_parameter(model_m(i)%array, &
                        model_reg(i)%array, model_name(i), 100.0/(maxval(abs(model_m(i)%array)) + float_tiny))
                end if

            end if

        end do

    end subroutine model_regularization

    !
    !> Apply regularization to model parameters
    !
    subroutine source_regularization

        real, allocatable, dimension(:, :) :: sx, sz, st0
        real, allocatable, dimension(:, :) :: sxr, szr
        integer :: i, ic
        type(hdbscan_param) :: p
        integer, allocatable, dimension(:) :: unindex, sindex
        integer :: nc
        real, allocatable, dimension(:) :: pr1, pr2
        character(len=32) :: fit_method
        real :: fit_smooth
        integer :: fit_order
        real, allocatable, dimension(:, :) :: meq, f
        integer :: dist, i1, i2, l, isx, isz, sp
        real, allocatable, dimension(:) :: fx, fz
        logical, allocatable, dimension(:) :: fi
        real :: d, dd
        character(len=1024) :: ml_python, ml_src
        character(len=1024) :: ml_model_infer, ml_model_refine
        integer :: ml_niter_refine
        real, allocatable, dimension(:) :: ml_xyz_weight
        real :: ml_max_dist

        if (any_in(['sx', 'sy', 'sz', 'st0'], model_name)) then

            sx = zeros(nr_virtual, 1)
            if (.not. any(model_name == 'sx')) then
                sx(:, 1) = gmtr(1)%recr(:)%x
            else
                sx = get_meta_array_core(model_m, 'sx')
            end if

            sz = zeros(nr_virtual, 1)
            if (.not. any(model_name == 'sz')) then
                sz(:, 1) = gmtr(1)%recr(:)%z
            else
                sz = get_meta_array_core(model_m, 'sz')
            end if

            st0 = zeros(nr_virtual, 1)
            if (.not. any(model_name == 'st0')) then
                st0(:, 1) = gmtr(1)%recr(:)%t0
            else
                st0 = get_meta_array_core(model_m, 'st0')
            end if

            sxr = sx
            szr = sz
            do i = 1, size(source_regularization_method)

                select case(source_regularization_method(i))

                    case('clustering-fit')

                        ! Clustering
                        p%n = nr_virtual
                        p%data = zeros(nr_virtual, 2)
                        p%data(:, 1) = szr(:, 1) - oz
                        p%data(:, 2) = sxr(:, 1) - ox
                        p%nd = 2
                        call readpar_xint(file_parameter, 'clustering_fit_min_sample', p%min_sample, nint(0.25*nr_virtual), iter*1.0)
                        call readpar_xint(file_parameter, 'clustering_fit_min_cluster_size', p%min_cluster_size, p%min_sample, iter*1.0)
                        call hdbscan(p)

                        ! Find unique labels
                        unindex = sort(unique(p%labels))
                        nc = size(unindex)

                        ! ic = 1 -> label = 0 -> noisy data
                        call readpar_xstring(file_parameter, 'clustering_fit_method', fit_method, 'polynomial', iter*1.0)
                        call readpar_xfloat(file_parameter, 'clustering_fit_smooth', fit_smooth, 0.5, iter*1.0)
                        call readpar_xint(file_parameter, 'clustering_fit_order', fit_order, 1, iter*1.0)

                        if (rankid == 0) then
                            call warn(date_time_compact()//' >>>>>>>>>> Number of clusters = '//num2str(nc))
                        endif

                        do ic = 2, nc

                            ! Select points associated with label_i
                            sindex = pack(regspace(1, 1, p%n), mask=(p%labels==unindex(ic)))

                            ! Fit the curve
                            call fit_curve(p%data(sindex, 1), p%data(sindex, 2), pr1, pr2, method=fit_method, smooth=fit_smooth, order=fit_order)

                            sxr(sindex, 1) = pr2
                            szr(sindex, 1) = pr1
                            ! st0r(sindex, 1) = mean(st0(sindex, 1))

                            call warn(date_time_compact()//' Source regularization for cluster '//num2str(ic)//' is done.')

                        end do

                    case ('ml')

                        meq = zeros(nz, nx)
                        f = zeros(nz, nx)

                        call readpar_string(file_parameter, 'reg_ml_python', ml_python, 'python', required=.true.)
                        call readpar_string(file_parameter, 'reg_ml_src', ml_src, '$HOME/src/latte/ml/main2.py', required=.true.)
                        call readpar_string(file_parameter, 'reg_ml_model_infer', ml_model_infer, '$HOME/src/latte/ml/infer2.model', required=.true.)
                        call readpar_int(file_parameter, 'reg_ml_niter_refine', ml_niter_refine, 3, required=.false.)
                        call readpar_string(file_parameter, 'reg_ml_model_refine', ml_model_refine, '$HOME/src/latte/ml/refine2.model', required=.true.)
                        call readpar_nfloat(file_parameter, 'reg_ml_xyz_weight', ml_xyz_weight, [1.0, 1.0])
                        call assert(size(ml_xyz_weight) == 2, ' <source_regularization> Error: size(ml_xyz_weight) must = 2')
                        call readpar_xfloat(file_parameter, 'reg_ml_max_dist', ml_max_dist, float_huge, iter*1.0)

                        if (rankid == 0) then

                            ! ML fault detection
                            dist = 3
                            do l = 1, nr_virtual
                                isx = nint((sxr(l, 1) - ox)/dx + 1)
                                isz = nint((szr(l, 1) - oz)/dz + 1)
                                do i2 = -2*dist  - 1, 2*dist + 1
                                    do i1 = -2*dist  - 1, 2*dist + 1
                                        if (isz + i1 >= 1 .and. isz + i1 <= nz &
                                                .and. isx + i2 >= 1 .and. isx + i2 <= nx) then
                                            meq(isz + i1, isx + i2) = &
                                                max(meq(isz + i1, isx + i2), exp(-0.3*(i1**2 + i2**2)))
                                        end if
                                    end do
                                end do
                            end do

                            call output_array(meq, dir_iter_model(iter)//'/source_image.bin')
                            call execute_command_line(tidy(ml_python)//" "//tidy(ml_src) &
                                //' --n1='//num2str(nz)//' --n2='//num2str(nx) &
                                //' --model_infer='//tidy(ml_model_infer) &
                                //' --model_refine='//tidy(ml_model_refine) &
                                //' --niter='//num2str(ml_niter_refine) &
                                //' --input='//dir_iter_model(iter)//'/source_image.bin' &
                                //' --output='//dir_iter_model(iter)//'/source_image.bin')
                            f = load(dir_iter_model(iter)//'/source_image.bin.fsem', nz, nx)

                        end if

                        call mpibarrier

                        call bcast_array(f)

                        ! Move source location to nearest fault position
                        fx = meshgrid([nz, nx], [dz, dx], [oz, ox], dim=2)
                        fz = meshgrid([nz, nx], [dz, dx], [oz, ox], dim=1)

                        fi = falses(nx*nz)
                        l = 1
                        do i2 = 1, nx
                            do i1 = 1, nz
                                if (f(i1, i2) >= 0.5) then
                                    fi(l) = .true.
                                end if
                                l = l + 1
                            end do
                        end do

                        fx = pack(fx, fi)
                        fz = pack(fz, fi)
                        do l = 1, nr_virtual
                            d = float_huge
                            do ic = 1, size(fx)
                                dd = ml_xyz_weight(1)*(fx(ic) - sxr(l, 1))**2 &
                                    + ml_xyz_weight(2)*(fz(ic) - szr(l, 1))**2
                                if (dd < d) then
                                    d = dd
                                    sp = ic
                                end if
                            end do
                            ! Only regularize if close enough
                            if (d <= ml_max_dist) then
                                sxr(l, 1) = fx(sp)
                                szr(l, 1) = fz(sp)
                            end if
                        end do

                end select

            end do

            ! Assign back
            call set_meta_array_core(model_reg, 'sx', sxr)
            call set_meta_array_core(model_reg, 'sz', szr)

            ! Output regularized model
            if (rankid == 0) then
                call output_array(get_meta_array_core(model_reg, 'sx'), dir_iter_model(iter)//'/reg_sx.bin')
                call output_array(get_meta_array_core(model_reg, 'sz'), dir_iter_model(iter)//'/reg_sz.bin')
            end if

        end if

        call mpibarrier

    end subroutine source_regularization

end module regularization
