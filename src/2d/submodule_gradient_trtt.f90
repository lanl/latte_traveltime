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


submodule(gradient) gradient_trtt

    use traveltime_iso_reflection

contains

    !
    !> Compute gradients shot by shot for TRTT
    !
    module subroutine compute_gradient_shots_trtt

        real, allocatable, dimension(:, :, :) :: ttp_all, tts_all
        real, allocatable, dimension(:, :) :: ttp_syn, ttp_obs, ttp_residual
        real, allocatable, dimension(:, :) :: tts_syn, tts_obs, tts_residual
        real, allocatable, dimension(:, :, :) :: lamp_all, lams_all, energyp_all, energys_all
        real, allocatable, dimension(:, :) :: lamp, lams
        real, allocatable, dimension(:, :, :) :: img
        integer :: i, l
        real, allocatable, dimension(:, :) :: vp, vs, refl
        character(len=1024) :: dir_field
        real :: eps

        ! temporary directory
        dir_scratch = tidy(dir_working)//'/scratch'
        dir_synthetic = dir_iter_synthetic(iter)
        dir_field = dir_iter_record(iter)
        if (rankid == 0) then
            call make_directory(dir_synthetic)
            call make_directory(dir_scratch)
        end if
        call mpibarrier

        ! step misfit set to zero
        step_misfit = 0.0d0

        ! Copy models
        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
                case ('refl')
                    refl = model_m(i)%array
                    ! Limit reflector image to the maximum number of reflectors
                    where (refl > nrefl*1.0)
                        refl = 0.0
                    end where
            end select
        end do
        do i = 1, nmodel_aux
            select case (model_aux(i)%name)
                case ('vp')
                    vp = model_aux(i)%array
                case ('vs')
                    vs = model_aux(i)%array
                case ('refl')
                    refl = model_aux(i)%array
                    ! Limit reflector image to the maximum number of reflectors
                    where (refl > nrefl*1.0)
                        refl = 0.0
                    end where
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

            call set_adaptive_model_range(gmtr(ishot))

            ! Compute transmission and reflection traveltimes
            select case (which_medium)

                case ('acoustic-iso')
                    call forward_iso_fast_sweep_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                        refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), ttp_all, ttp_syn)
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')

                case ('elastic-iso')
                    select case (incident_wave)
                        case ('p')
                            call forward_iso_fast_sweep_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                ttp_all, tts_all, ttp_syn, tts_syn)
                        case ('s')
                            call forward_iso_fast_sweep_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                tts_all, ttp_all, tts_syn, ttp_syn)
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

                    ! Compute adjoint-state field
                    call adjoint_iso_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_all, &
                        refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        ttp_residual, lamp_all)

                    if (yn_precond) then

                        ! Compute preconditioner
                        call adjoint_iso_reflection( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_all, &
                            refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            ones_like(ttp_residual), energyp_all)

                        lamp = zeros(shot_nz, shot_nx)
                        do i = 1, nrefl + 1
                            eps = precond_eps*maxval(abs(energyp_all(:, :, i)))
                            if (eps > 0) then
                                lamp = lamp - misfit_weight(i)*lamp_all(:, :, i)/(energyp_all(:, :, i) + eps)
                            end if
                        end do

                    else

                        lamp = zeros(shot_nz, shot_nx)
                        do i = 1, nrefl + 1
                            lamp = lamp + misfit_weight(i)*lamp_all(:, :, i)
                        end do
                        lamp = -lamp/vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)**3

                    end if

                    do i = 1, nmodel
                        if (model_m(i)%name == 'vp') then
                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lamp, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lamp, 'grad_'//tidy(model_m(i)%name))
                            end if
                        end if
                    end do

                    if (yn_update_reflector) then
                        ! Compute gradient of reflector image
                        call image_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_obs(:, 2:), img)
                        do i = 1, nmodel
                            if (model_m(i)%name == 'refl') then
                                do l = 1, nrefl
                                    call process_model_single_shot(ishot, refl_all(:, :, l), img(:, :, l), model_m(i)%name)
                                end do
                            end if
                        end do
                    end if

                case ('elastic-iso')

                    select case (incident_wave)
                        case ('p')
                            ! Compute adjoint-state field
                            call adjoint_iso_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_all, tts_all, &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                ttp_residual, tts_residual, lamp_all, lams_all)

                            if (yn_precond) then
                                ! Compute preconditioner
                                call adjoint_iso_reflection_elastic( &
                                    vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                    vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                    [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_all, tts_all, &
                                    refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                    ones_like(ttp_residual), ones_like(tts_residual), energyp_all, energys_all)
                            end if

                        case ('s')
                            ! Compute adjoint-state field
                            call adjoint_iso_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts_all, ttp_all, &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                tts_residual, ttp_residual, lams_all, lamp_all)

                            if (yn_precond) then
                                ! Compute preconditioner
                                call adjoint_iso_reflection_elastic( &
                                    vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                    vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                    [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts_all, ttp_all, &
                                    refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                    ones_like(tts_residual), ones_like(ttp_residual), energys_all, energyp_all)
                            end if

                    end select

                    if (yn_precond) then

                        lamp = zeros(shot_nz, shot_nx)
                        lams = zeros(shot_nz, shot_nx)
                        do i = 1, nrefl + 1
                            eps = maxval(abs(energyp_all(:, :, i)))*precond_eps
                            if (eps > 0) then
                                lamp = lamp - misfit_weight(i)*lamp_all(:, :, i)/(energyp_all(:, :, i) + eps)
                            end if
                            eps = maxval(abs(energys_all(:, :, i)))*precond_eps
                            if (eps > 0) then
                                lams = lams - misfit_weight(i)*lams_all(:, :, i)/(energys_all(:, :, i) + eps)
                            end if
                        end do

                    else

                        lamp = zeros(shot_nz, shot_nx)
                        lams = zeros(shot_nz, shot_nx)
                        do i = 1, nrefl + 1
                            lamp = lamp + misfit_weight(i)*lamp_all(:, :, i)
                            lams = lams + misfit_weight(i)*lams_all(:, :, i)
                        end do
                        lamp = -lamp/vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)**3
                        lams = -lams/vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)**3

                    end if

                    do i = 1, nmodel
                        if (model_m(i)%name == 'vp') then
                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lamp, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lamp, 'grad_'//tidy(model_m(i)%name))
                            end if
                        end if
                        if (model_m(i)%name == 'vs') then
                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lams, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lams, 'grad_'//tidy(model_m(i)%name))
                            end if
                        end if
                    end do

                    if (yn_update_reflector) then
                        ! Compute gradient of reflector image
                        call image_iso_elastic( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_obs(:, 2:), tts_obs(:, 2:), img)
                        do i = 1, nmodel
                            if (model_m(i)%name == 'refl') then
                                do l = 1, nrefl
                                    call process_model_single_shot(ishot, refl_all(:, :, l), img(:, :, l), model_m(i)%name)
                                end do
                            end if
                        end do
                    end if

            end select

            call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                //' gradient computation is done. ')

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
        if (yn_update_reflector) then
            call allreduce_array(refl_all)
        end if

        call mpibarrier

    end subroutine compute_gradient_shots_trtt

end submodule gradient_trtt
