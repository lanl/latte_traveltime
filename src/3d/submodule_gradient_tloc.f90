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


submodule(gradient) gradient_tloc

    use traveltime_iso

contains

    !
    !> Compute gradients shot by shot for TLOC-AD
    !
    module subroutine compute_gradient_shots_tloc

        real, allocatable, dimension(:, :, :) :: ttp, tts
        real, allocatable, dimension(:, :) :: ttp_syn, ttp_obs, ttp_residual
        real, allocatable, dimension(:, :) :: tts_syn, tts_obs, tts_residual
        real, allocatable, dimension(:, :, :) :: lam, energy
        integer :: i, j, irx, iry, irz
        real :: dxt, dyt, dzt
        real, allocatable, dimension(:, :, :) :: vp, vs, sx, sy, sz, st0
        character(len=1024) :: dir_field

        ! temporary directory
        dir_scratch = tidy(dir_working)//'/scratch'
        dir_synthetic = dir_iter_synthetic(iter)
        dir_field = dir_iter_record(iter)
        if (rankid == 0) then
            call make_directory(dir_synthetic)
            call make_directory(dir_scratch)
        end if
        call mpibarrier

        ! Step misfit set to zero
        step_misfit = 0.0d0

        ! Set default values for real source parameters
        sx = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'sx')) then
            sx(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%x
        end if

        sy = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'sy')) then
            sy(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%y
        end if

        sz = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'sz')) then
            sz(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%z
        end if

        st0 = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'st0')) then
            st0(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%t0
        end if

        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
                case ('sx')
                    sx = model_m(i)%array
                case ('sy')
                    sy = model_m(i)%array
                case ('sz')
                    sz = model_m(i)%array
                case ('st0')
                    st0 = model_m(i)%array
            end select
        end do
        do i = 1, nmodel_aux
            select case (model_aux(i)%name)
                case ('vp')
                    vp = model_aux(i)%array
                case ('vs')
                    vs = model_aux(i)%array
                case ('sx')
                    sx = model_aux(i)%array
                case ('sy')
                    sy = model_aux(i)%array
                case ('sz')
                    sz = model_aux(i)%array
                case ('st0')
                    st0 = model_aux(i)%array
            end select
        end do

        if (rankid == 0) then
            do i = 1, nmodel
                call plot_histogram(model_m(i)%array, &
                    label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
            end do
        end if

        call mpibarrier

        ! Forward modeling to compute misfit
        do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

            ! Updated virtual receiver (real source) information
            do i = 1, gmtr(ishot)%nr
                if (any(model_name == 'sx')) then
                    gmtr(ishot)%recr(i)%x = sx(i, 1, 1)
                end if
                if (any(model_name == 'sy')) then
                    gmtr(ishot)%recr(i)%y = sy(i, 1, 1)
                end if
                if (any(model_name == 'sz')) then
                    gmtr(ishot)%recr(i)%z = sz(i, 1, 1)
                end if
                if (any(model_name == 'st0')) then
                    gmtr(ishot)%recr(i)%t0 = st0(i, 1, 1)
                end if
            end do

            ! Set adaptive range
            call set_adaptive_model_range(gmtr(ishot))

            select case (which_medium)

                case ('acoustic-iso')
                    select case (forward_eikonal_method)
                        case ('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                        case ('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                    end select
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')

                case ('elastic-iso')
                    select case (forward_eikonal_method)
                        case ('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                            call forward_iso_fast_march( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, tts_syn)
                        case ('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                            call forward_iso_fast_sweep( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, tts_syn)
                    end select
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    call output_array(tts_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')

            end select

            call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id) &
                //' traveltime computation is done. ')

            if (sum(data_misfit) == 0) then
                ! Copy initial iteration synthetic data to dir_working/iteration_0
                call make_directory(dir_iter_synthetic(0))
                select case (which_medium)
                    case ('acoustic-iso', 'acoustic-tti')
                        call copy_file( &
                            tidy(dir_synthetic)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', &
                            tidy(dir_iter_synthetic(0))//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    case ('elastic-iso', 'elastic-tti')
                        call copy_file( &
                            tidy(dir_synthetic)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', &
                            tidy(dir_iter_synthetic(0))//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                        call copy_file( &
                            tidy(dir_synthetic)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin', &
                            tidy(dir_iter_synthetic(0))//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                end select
                ! Copy initial models to dir_working/iteration_0
                if (rankid == 0) then
                    call make_directory(tidy(dir_working)//'/iteration_0/model')
                    do i = 1, nmodel
                        call output_array(model_m(i)%array, tidy(dir_working)//'/iteration_0/model/' &
                            //tidy(model_name(i))//'.bin')
                    end do
                end if
            end if

            step_misfit(ishot) = 0
            select case (which_medium)
                case ('acoustic-iso', 'acoustic-tti')
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot))
                case ('elastic-iso', 'elastic-tti')
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot))
                    call compute_shot_misfit(ishot, dir_field, 's', tts_syn, tts_obs, tts_residual, step_misfit(ishot))
            end select

            call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                //' misfit = '//num2str(step_misfit(ishot), '(es)'))

            if (yn_misfit_only) then
                cycle
            end if

            ! Adjoint modeling to compute gradient
            select case (which_medium)

                case ('acoustic-iso')

                    do i = 1, nmodel

                        ! Update model parameters, if necessary
                        if (model_name(i) == 'vp') then

                            ! Compute adjoint-state field
                            call adjoint_iso( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                ttp_residual, lam)

                            if (yn_precond) then

                                ! Compute preconditioner
                                call adjoint_iso( &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                    ones_like(ttp_residual), energy)

                                lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                            else

                                lam = -lam/vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

                            end if

                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(1))
                            end if

                            ! Update source parameters
                        else if (model_name(i) == 'sx') then
                            ! grad_sx = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂x, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dxt = (ttp(irz, iry, clip(irx + 1, 1, shot_nx)) - ttp(irz, iry, clip(irx - 1, 1, shot_nx)))/(2*dx)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dxt
                            end do

                        else if (model_name(i) == 'sy') then
                            ! grad_sy = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂y, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dyt = (ttp(irz, clip(iry + 1, 1, shot_ny), irx) - ttp(irz, clip(iry - 1, 1, shot_ny), irx))/(2*dy)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dyt
                            end do

                        else if (model_name(i) == 'sz') then
                            ! grad_sz = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂z, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dzt = (ttp(clip(irz + 1, 1, shot_nz), iry, irx) - ttp(clip(irz - 1, 1, shot_nz), iry, irx))/(2*dz)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dzt
                            end do

                        else if (model_name(i) == 'st0') then
                            ! grad_st0 = (t_syn - (t_obs - t_0))*δ(x - x_vr)

                            do j = 1, gmtr(ishot)%nr
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)
                            end do

                        end if

                    end do

                    call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                        //' gradient computation is done. ')

                case ('elastic-iso')

                    do i = 1, nmodel

                        ! Update model parameters, if necessary
                        if (model_name(i) == 'vp') then

                            ! For vp
                            ! Compute adjoint-state field
                            call adjoint_iso( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                ttp_residual, lam)

                            if (yn_precond) then

                                ! Compute preconditioner
                                call adjoint_iso( &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                    ones_like(ttp_residual), energy)

                                lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                            else

                                lam = -lam/vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

                            end if

                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(i))
                            end if

                        else if (model_name(i) == 'vs') then

                            ! For vs
                            ! Compute adjoint-state field
                            call adjoint_iso( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, &
                                tts_residual, lam)

                            if (yn_precond) then

                                ! Compute preconditioner
                                call adjoint_iso( &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, &
                                    ones_like(tts_residual), energy)

                                lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                            else

                                lam = -lam/vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

                            end if

                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(i))
                            end if

                            ! Update source parameters
                        else if (model_name(i) == 'sx') then
                            ! grad_sx = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂x, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dxt = (ttp(irz, iry, clip(irx + 1, 1, shot_nx)) - ttp(irz, iry, clip(irx - 1, 1, shot_nx)))/(2*dx)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dxt
                                dxt = (tts(irz, iry, clip(irx + 1, 1, shot_nx)) - tts(irz, iry, clip(irx - 1, 1, shot_nx)))/(2*dx)
                                model_grad(i)%array(j, 1, 1) = model_grad(i)%array(j, 1, 1) + tts_residual(j, 1)*dxt
                            end do

                        else if (model_name(i) == 'sy') then
                            ! grad_sy = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂y, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dyt = (ttp(irz, clip(iry + 1, 1, shot_ny), irx) - ttp(irz, clip(iry - 1, 1, shot_ny), irx))/(2*dy)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dyt
                                dyt = (tts(irz, clip(iry + 1, 1, shot_ny), irx) - tts(irz, clip(iry - 1, 1, shot_ny), irx))/(2*dy)
                                model_grad(i)%array(j, 1, 1) = model_grad(i)%array(j, 1, 1) + tts_residual(j, 1)*dyt
                            end do

                        else if (model_name(i) == 'sz') then
                            ! grad_sz = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂z, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dzt = (ttp(clip(irz + 1, 1, shot_nz), iry, irx) - ttp(clip(irz - 1, 1, shot_nz), iry, irx))/(2*dz)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dzt
                                dzt = (tts(clip(irz + 1, 1, shot_nz), iry, irx) - tts(clip(irz - 1, 1, shot_nz), iry, irx))/(2*dz)
                                model_grad(i)%array(j, 1, 1) = model_grad(i)%array(j, 1, 1) + tts_residual(j, 1)*dzt
                            end do

                        else if (model_name(i) == 'st0') then
                            ! grad_st0 = (t_syn - (t_obs - t_0))*δ(x - x_vr)

                            do j = 1, gmtr(ishot)%nr
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1) &
                                    + tts_residual(j, 1)
                            end do

                        end if

                    end do

                    call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                        //' gradient computation is done. ')

            end select

        end do

        call mpibarrier

        ! collect misfit
        call allreduce_array(step_misfit)
        if (sum(data_misfit) == 0) then
            data_misfit(0) = sum(step_misfit)
            shot_misfit(:, 0) = step_misfit
        end if
        data_misfit(iter) = sum(step_misfit)

        if (yn_misfit_only) then
            return
        end if

        call mpibarrier

        do i = 1, nmodel
            call allreduce_array(model_grad(i)%array)
        end do

        call mpibarrier

    end subroutine compute_gradient_shots_tloc

    !
    !> Compute gradients shot by shot for TLOC-DD
    !
    module subroutine compute_gradient_shots_tloc_dd

        real, allocatable, dimension(:, :, :) :: ttp, tts
        real, allocatable, dimension(:, :) :: ttp_syn, ttp_obs, ttp_residual
        real, allocatable, dimension(:, :) :: tts_syn, tts_obs, tts_residual
        real, allocatable, dimension(:, :, :) :: lam, energy
        integer :: i, j, irx, iry, irz
        real :: dxt, dyt, dzt
        real, allocatable, dimension(:, :, :) :: vp, vs, sx, sy, sz, st0
        character(len=1024) :: dir_field
        real, allocatable, dimension(:, :) :: tpsyn_all, tssyn_all, tpobs_all, tsobs_all

        ! temporary directory
        dir_scratch = tidy(dir_working)//'/scratch'
        dir_synthetic = dir_iter_synthetic(iter)
        dir_field = dir_iter_record(iter)
        if (rankid == 0) then
            call make_directory(dir_synthetic)
            call make_directory(dir_scratch)
        end if
        call mpibarrier

        ! Step misfit set to zero
        step_misfit = 0.0d0

        ! Set default values for real source parameters
        sx = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'sx')) then
            sx(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%x
        end if

        sy = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'sy')) then
            sy(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%y
        end if

        sz = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'sz')) then
            sz(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%z
        end if

        st0 = zeros(nr_virtual, 1, 1)
        if (.not. any(model_name == 'st0')) then
            st0(:, 1, 1) = gmtr(shot_in_rank(rankid, 1))%recr(:)%t0
        end if

        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
                case ('sx')
                    sx = model_m(i)%array
                case ('sy')
                    sy = model_m(i)%array
                case ('sz')
                    sz = model_m(i)%array
                case ('st0')
                    st0 = model_m(i)%array
            end select
        end do
        do i = 1, nmodel_aux
            select case (model_aux(i)%name)
                case ('vp')
                    vp = model_aux(i)%array
                case ('vs')
                    vs = model_aux(i)%array
                case ('sx')
                    sx = model_aux(i)%array
                case ('sy')
                    sy = model_aux(i)%array
                case ('sz')
                    sz = model_aux(i)%array
                case ('st0')
                    st0 = model_aux(i)%array
            end select
        end do

        if (rankid == 0) then
            do i = 1, nmodel
                call plot_histogram(model_m(i)%array, &
                    label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
            end do
        end if

        call mpibarrier

        select case (which_medium)
            case ('acoustic-iso')
                tpsyn_all = zeros(nr_virtual, ns)
                tpobs_all = zeros(nr_virtual, ns)
            case ('elastic-iso')
                tpsyn_all = zeros(nr_virtual, ns)
                tssyn_all = zeros(nr_virtual, ns)
                tpobs_all = zeros(nr_virtual, ns)
                tsobs_all = zeros(nr_virtual, ns)
        end select

        ! Forward modeling to compute misfit
        do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

            ! Updated virtual receiver (real source) information
            do i = 1, gmtr(ishot)%nr
                if (any(model_name == 'sx')) then
                    gmtr(ishot)%recr(i)%x = sx(i, 1, 1)
                end if
                if (any(model_name == 'sy')) then
                    gmtr(ishot)%recr(i)%y = sy(i, 1, 1)
                end if
                if (any(model_name == 'sz')) then
                    gmtr(ishot)%recr(i)%z = sz(i, 1, 1)
                end if
                if (any(model_name == 'st0')) then
                    gmtr(ishot)%recr(i)%t0 = st0(i, 1, 1)
                end if
            end do

            ! Set adaptive range
            call set_adaptive_model_range(gmtr(ishot))

            select case (which_medium)

                case ('acoustic-iso')
                    select case (forward_eikonal_method)
                        case ('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                        case ('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                    end select
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')

                    tpsyn_all(:, ishot) = ttp_syn(:, 1)
                    tpobs_all(:, ishot) = load(tidy(dir_field)//'/shot_'//num2str(ishot)//'_traveltime_p.bin', nr_virtual)

                case ('elastic-iso')
                    select case (forward_eikonal_method)
                        case ('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                            call forward_iso_fast_march( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, tts_syn)
                        case ('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                            call forward_iso_fast_sweep( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, tts_syn)
                    end select
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    call output_array(tts_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')

                    tpsyn_all(:, ishot) = ttp_syn(:, 1)
                    tpobs_all(:, ishot) = load(tidy(dir_field)//'/shot_'//num2str(ishot)//'_traveltime_p.bin', nr_virtual)
                    tssyn_all(:, ishot) = tts_syn(:, 1)
                    tsobs_all(:, ishot) = load(tidy(dir_field)//'/shot_'//num2str(ishot)//'_traveltime_s.bin', nr_virtual)

            end select

            call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id) &
                //' traveltime computation is done. ')

            if (sum(data_misfit) == 0) then
                ! Copy initial iteration synthetic data to dir_working/iteration_0
                call make_directory(dir_iter_synthetic(0))
                select case (which_medium)
                    case ('acoustic-iso', 'acoustic-tti')
                        call copy_file( &
                            tidy(dir_synthetic)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', &
                            tidy(dir_iter_synthetic(0))//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    case ('elastic-iso', 'elastic-tti')
                        call copy_file( &
                            tidy(dir_synthetic)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', &
                            tidy(dir_iter_synthetic(0))//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                        call copy_file( &
                            tidy(dir_synthetic)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin', &
                            tidy(dir_iter_synthetic(0))//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                end select
                ! Copy initial models to dir_working/iteration_0
                if (rankid == 0) then
                    call make_directory(tidy(dir_working)//'/iteration_0/model')
                    do i = 1, nmodel
                        call output_array(model_m(i)%array, tidy(dir_working)//'/iteration_0/model/' &
                            //tidy(model_name(i))//'.bin')
                    end do
                end if
            end if

        end do

        call mpibarrier

        ! Get all data
        select case (which_medium)
            case ('acoustic-iso', 'acoustic-tti')
                call allreduce_array(tpsyn_all)
                call allreduce_array(tpobs_all)
                if (rankid == 0) then
                    call output_array(tpsyn_all, tidy(dir_synthetic)//'/tp_all.bin')
                    if (sum(data_misfit) == 0) then
                        call output_array(tpsyn_all, tidy(dir_iter_synthetic(0))//'/tp_all.bin')
                    end if
                    if (iter == 1) then
                        call output_array(tpobs_all, tidy(dir_field)//'/tp_all.bin')
                    end if
                end if
            case ('elastic-iso', 'elastic-tti')
                call allreduce_array(tpsyn_all)
                call allreduce_array(tpobs_all)
                call allreduce_array(tssyn_all)
                call allreduce_array(tsobs_all)
                if (rankid == 0) then
                    call output_array(tpsyn_all, tidy(dir_synthetic)//'/tp_all.bin')
                    call output_array(tssyn_all, tidy(dir_synthetic)//'/ts_all.bin')
                    if (sum(data_misfit) == 0) then
                        call output_array(tpsyn_all, tidy(dir_iter_synthetic(0))//'/tp_all.bin')
                        call output_array(tssyn_all, tidy(dir_iter_synthetic(0))//'/ts_all.bin')
                    end if
                    if (iter == 1) then
                        call output_array(tpobs_all, tidy(dir_field)//'/tp_all.bin')
                        call output_array(tsobs_all, tidy(dir_field)//'/ts_all.bin')
                    end if
                end if
        end select

        call mpibarrier

        ! Resume misfit computation and if necessary, gradient computation
        do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

            ! Set adaptive range
            call set_adaptive_model_range(gmtr(ishot))

            step_misfit(ishot) = 0
            select case (which_medium)
                case ('acoustic-iso', 'acoustic-tti')
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot), tpsyn_all, tpobs_all)
                case ('elastic-iso', 'elastic-tti')
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot), tpsyn_all, tpobs_all)
                    call compute_shot_misfit(ishot, dir_field, 's', tts_syn, tts_obs, tts_residual, step_misfit(ishot), tssyn_all, tsobs_all)
            end select

            call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                //' misfit = '//num2str(step_misfit(ishot), '(es)'))

            if (yn_misfit_only) then
                cycle
            end if

            ! Adjoint modeling to compute gradient
            select case (which_medium)

                case ('acoustic-iso')

                    do i = 1, nmodel

                        ! Update model parameters, if necessary
                        if (model_name(i) == 'vp') then

                            ! Compute adjoint-state field
                            call adjoint_iso( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                ttp_residual, lam)

                            if (yn_precond) then

                                ! Compute preconditioner
                                call adjoint_iso( &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                    ones_like(ttp_residual), energy)

                                lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                            else

                                lam = -lam/vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

                            end if

                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(1))
                            end if

                            ! Update source parameters
                        else if (model_name(i) == 'sx') then
                            ! grad_sx = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂x, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dxt = (ttp(irz, iry, clip(irx + 1, 1, shot_nx)) - ttp(irz, iry, clip(irx - 1, 1, shot_nx)))/(2*dx)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dxt
                            end do

                        else if (model_name(i) == 'sy') then
                            ! grad_sy = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂y, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dyt = (ttp(irz, clip(iry + 1, 1, shot_ny), irx) - ttp(irz, clip(iry - 1, 1, shot_ny), irx))/(2*dy)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dyt
                            end do

                        else if (model_name(i) == 'sz') then
                            ! grad_sz = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂z, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dzt = (ttp(clip(irz + 1, 1, shot_nz), iry, irx) - ttp(clip(irz - 1, 1, shot_nz), iry, irx))/(2*dz)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dzt
                            end do

                        else if (model_name(i) == 'st0') then
                            ! grad_st0 = (t_syn - (t_obs - t_0))*δ(x - x_vr)

                            do j = 1, gmtr(ishot)%nr
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)
                            end do

                        end if

                    end do

                    call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                        //' gradient computation is done. ')

                case ('elastic-iso')

                    do i = 1, nmodel

                        ! Update model parameters, if necessary
                        if (model_name(i) == 'vp') then

                            ! For vp
                            ! Compute adjoint-state field
                            call adjoint_iso( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                ttp_residual, lam)

                            if (yn_precond) then

                                ! Compute preconditioner
                                call adjoint_iso( &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, &
                                    ones_like(ttp_residual), energy)

                                lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                            else

                                lam = -lam/vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

                            end if

                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(i))
                            end if

                        else if (model_name(i) == 'vs') then

                            ! For vs
                            ! Compute adjoint-state field
                            call adjoint_iso( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, &
                                tts_residual, lam)

                            if (yn_precond) then

                                ! Compute preconditioner
                                call adjoint_iso( &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, &
                                    ones_like(tts_residual), energy)

                                lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                            else

                                lam = -lam/vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

                            end if

                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(i))
                            end if

                            ! Update source parameters
                        else if (model_name(i) == 'sx') then
                            ! grad_sx = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂x, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dxt = (ttp(irz, iry, clip(irx + 1, 1, shot_nx)) - ttp(irz, iry, clip(irx - 1, 1, shot_nx)))/(2*dx)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dxt
                                dxt = (tts(irz, iry, clip(irx + 1, 1, shot_nx)) - tts(irz, iry, clip(irx - 1, 1, shot_nx)))/(2*dx)
                                model_grad(i)%array(j, 1, 1) = model_grad(i)%array(j, 1, 1) + tts_residual(j, 1)*dxt
                            end do

                        else if (model_name(i) == 'sy') then
                            ! grad_sy = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂y, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dyt = (ttp(irz, clip(iry + 1, 1, shot_ny), irx) - ttp(irz, clip(iry - 1, 1, shot_ny), irx))/(2*dy)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dyt
                                dyt = (tts(irz, clip(iry + 1, 1, shot_ny), irx) - tts(irz, clip(iry - 1, 1, shot_ny), irx))/(2*dy)
                                model_grad(i)%array(j, 1, 1) = model_grad(i)%array(j, 1, 1) + tts_residual(j, 1)*dyt
                            end do

                        else if (model_name(i) == 'sz') then
                            ! grad_sz = (t_syn - (t_obs - t_0))*δ(x - x_vr)*∂T_syn/∂z, is a scalar value at each virtual receiver

                            do j = 1, gmtr(ishot)%nr
                                irx = nint((sx(j, 1, 1) - shot_xbeg)/dx) + 1
                                iry = nint((sy(j, 1, 1) - shot_ybeg)/dy) + 1
                                irz = nint((sz(j, 1, 1) - shot_zbeg)/dz) + 1
                                dzt = (ttp(clip(irz + 1, 1, shot_nz), iry, irx) - ttp(clip(irz - 1, 1, shot_nz), iry, irx))/(2*dz)
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1)*dzt
                                dzt = (tts(clip(irz + 1, 1, shot_nz), iry, irx) - tts(clip(irz - 1, 1, shot_nz), iry, irx))/(2*dz)
                                model_grad(i)%array(j, 1, 1) = model_grad(i)%array(j, 1, 1) + tts_residual(j, 1)*dzt
                            end do

                        else if (model_name(i) == 'st0') then
                            ! grad_st0 = (t_syn - (t_obs - t_0))*δ(x - x_vr)

                            do j = 1, gmtr(ishot)%nr
                                model_grad(i)%array(j, 1, 1) = ttp_residual(j, 1) + tts_residual(j, 1)
                            end do

                        end if

                    end do

                    call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                        //' gradient computation is done. ')

            end select

        end do

        call mpibarrier

        ! collect misfit
        call allreduce_array(step_misfit)
        if (sum(data_misfit) == 0) then
            data_misfit(0) = sum(step_misfit)
            shot_misfit(:, 0) = step_misfit
        end if
        data_misfit(iter) = sum(step_misfit)

        if (yn_misfit_only) then
            return
        end if

        call mpibarrier

        do i = 1, nmodel
            call allreduce_array(model_grad(i)%array)
        end do

        call mpibarrier

    end subroutine compute_gradient_shots_tloc_dd

end submodule gradient_tloc
