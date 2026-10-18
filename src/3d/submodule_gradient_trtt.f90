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


submodule(gradient) gradient_trtt

    use traveltime_iso_reflection

contains

    !
    !> Compute gradients shot by shot for TRTT
    !
    module subroutine compute_gradient_shots_trtt

        real, allocatable, dimension(:, :, :, :) :: ttp_all, tts_all
        real, allocatable, dimension(:, :) :: ttp_syn, ttp_obs, ttp_residual
        real, allocatable, dimension(:, :) :: tts_syn, tts_obs, tts_residual
        real, allocatable, dimension(:, :, :, :) :: lamp_all, lams_all, energyp_all, energys_all
        real, allocatable, dimension(:, :, :, :) :: lamr, lamr2, energyr, energyr2
        real, allocatable, dimension(:, :, :) :: lamp, lams
        real, allocatable, dimension(:, :) :: wgtp, wgts
        integer :: i, irefl
        logical :: joint
        real, allocatable, dimension(:, :, :) :: vp, vs, refl
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
        irefl = 0
        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
                case ('refl')
                    irefl = i
                    refl = model_m(i)%array
                    ! Ignore the reflector labels above nrefl
                    where (nint(refl) > nrefl)
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
                    ! Ignore the reflector labels above nrefl
                    where (nint(refl) > nrefl)
                        refl = 0.0
                    end where
            end select
        end do

        ! With step_max_refl > 0, the reflector depths take part in the joint
        ! update: accumulate their Gauss-Newton step, as refine_reflector does,
        ! from the reflector-point values of the adjoint fields below
        joint = reflector_joint()

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

            select case (which_medium)

                case ('acoustic-iso')
                    ! Compute transmission and reflection traveltimes
                    call forward_iso_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), ttp_all, ttp_syn)
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')

                case ('elastic-iso')
                    ! Compute transmission and reflection traveltimes
                    select case (incident_wave)
                        case ('p')
                            call forward_iso_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                ttp_all, tts_all, ttp_syn, tts_syn)
                        case ('s')
                            call forward_iso_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                tts_all, ttp_all, tts_syn, ttp_syn)
                    end select
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    call output_array(tts_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')

            end select

            call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id) &
                //' traveltime computation completed. ')

            if (sum(data_misfit) == 0 .and. .not. yn_misfit_only) then
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
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot), &
                        weight_used=wgtp)
                case ('elastic-iso', 'elastic-tti')
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot), &
                        weight_used=wgtp)
                    call compute_shot_misfit(ishot, dir_field, 's', tts_syn, tts_obs, tts_residual, step_misfit(ishot), &
                        weight_used=wgts)
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
                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, &
                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        ttp_residual, lamp_all, lamr)

                    if (yn_precond) then

                        ! Compute preconditioner
                        call adjoint_iso_reflection( &
                            vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, &
                            refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            ones_like(ttp_residual), energyp_all, energyr)

                        lamp = zeros(shot_nz, shot_ny, shot_nx)
                        do i = 1, nrefl + 1
                            eps = precond_eps*maxval(abs(energyp_all(:, :, :, i)))
                            if (eps > 0) then
                                lamp = lamp - misfit_weight(i)*lamp_all(:, :, :, i)/(energyp_all(:, :, :, i) + eps)
                            end if
                        end do

                    else

                        lamp = zeros(shot_nz, shot_ny, shot_nx)
                        do i = 1, nrefl + 1
                            lamp = lamp + misfit_weight(i)*lamp_all(:, :, :, i)
                        end do
                        lamp = -lamp/vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

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

                    if (joint) then
                        if (.not. yn_precond .or. any(wgtp(:, 2:) /= 1)) then
                            call adjoint_iso_reflection_at_reflector( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), wgtp, energyr)
                        end if
                        call accumulate_reflector_step(irefl, ttp_all(:, :, :, 1), &
                            vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), lamr, energyr)
                    end if

                case ('elastic-iso')

                    select case (incident_wave)
                        case ('p')
                            ! Compute adjoint-state field
                            call adjoint_iso_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, tts_all, &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                ttp_residual, tts_residual, lamp_all, lams_all, lamr, lamr2)

                            if (yn_precond) then
                                ! Compute preconditioner
                                call adjoint_iso_reflection_elastic( &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, tts_all, &
                                    refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    ones_like(ttp_residual), ones_like(tts_residual), energyp_all, energys_all, &
                                    energyr, energyr2)
                            end if

                            if (joint) then
                                if (.not. yn_precond .or. any(wgtp(:, 2:) /= 1) .or. any(wgts(:, 2:) /= 1)) then
                                    call adjoint_iso_reflection_at_reflector( &
                                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, &
                                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), wgtp, energyr)
                                    call adjoint_iso_reflection_at_reflector( &
                                        vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts_all, &
                                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), wgts, energyr2)
                                end if
                                ! PP and PS reflections
                                call accumulate_reflector_step(irefl, ttp_all(:, :, :, 1), &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), lamr, energyr)
                                call accumulate_reflector_step(irefl, ttp_all(:, :, :, 1), &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), lamr2, energyr2, &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                            end if

                        case ('s')
                            ! Compute adjoint-state field
                            call adjoint_iso_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts_all, ttp_all, &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                tts_residual, ttp_residual, lams_all, lamp_all, lamr, lamr2)

                            if (yn_precond) then
                                ! Compute preconditioner
                                call adjoint_iso_reflection_elastic( &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts_all, ttp_all, &
                                    refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    ones_like(tts_residual), ones_like(ttp_residual), energys_all, energyp_all, &
                                    energyr, energyr2)
                            end if

                            if (joint) then
                                if (.not. yn_precond .or. any(wgtp(:, 2:) /= 1) .or. any(wgts(:, 2:) /= 1)) then
                                    call adjoint_iso_reflection_at_reflector( &
                                        vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts_all, &
                                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), wgts, energyr)
                                    call adjoint_iso_reflection_at_reflector( &
                                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_all, &
                                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), wgtp, energyr2)
                                end if
                                ! SS and SP reflections
                                call accumulate_reflector_step(irefl, tts_all(:, :, :, 1), &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), lamr, energyr)
                                call accumulate_reflector_step(irefl, tts_all(:, :, :, 1), &
                                    vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                    refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), lamr2, energyr2, &
                                    vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                            end if

                    end select

                    if (yn_precond) then

                        lamp = zeros(shot_nz, shot_ny, shot_nx)
                        lams = zeros(shot_nz, shot_ny, shot_nx)
                        do i = 1, nrefl + 1
                            eps = precond_eps*maxval(abs(energyp_all(:, :, :, i)))
                            if (eps > 0) then
                                lamp = lamp - misfit_weight(i)*lamp_all(:, :, :, i)/(energyp_all(:, :, :, i) + eps)
                            end if
                            eps = precond_eps*maxval(abs(energys_all(:, :, :, i)))
                            if (eps > 0) then
                                lams = lams - misfit_weight(i)*lams_all(:, :, :, i)/(energys_all(:, :, :, i) + eps)
                            end if
                        end do

                    else

                        lamp = zeros(shot_nz, shot_ny, shot_nx)
                        lams = zeros(shot_nz, shot_ny, shot_nx)
                        do i = 1, nrefl + 1
                            lamp = lamp + misfit_weight(i)*lamp_all(:, :, :, i)
                            lams = lams + misfit_weight(i)*lams_all(:, :, :, i)
                        end do
                        lamp = -lamp/vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3
                        lams = -lams/vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)**3

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

            end select

        end do

        call mpibarrier

        ! collect misfit
        call allreduce_array(step_misfit)
        if (sum(data_misfit) == 0 .and. .not. yn_misfit_only) then
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
        if (joint) then
            call allreduce_array(refl_den)
        end if
        call mpibarrier

    end subroutine compute_gradient_shots_trtt

    !
    !> Accumulate the Gauss-Newton step of the reflector depths at the current
    !> models (see accumulate_reflector_step). For each shot, it computes the
    !> transmission and reflection traveltimes, their residuals, and the
    !> receiver-side adjoint fields of the reflected legs, driven by the
    !> residuals and by the data weights of the misfit. Converted reflections
    !> (PS, SP) take part in elastic media. The misfit of the current models
    !> is left in step_misfit. Per shot, it solves 1 + nrefl eikonal equations
    !> (2*nrefl + 1 in elastic media) and 2*nrefl adjoint equations (4*nrefl
    !> in elastic media).
    !
    module subroutine compute_reflector_step_shots

        real, allocatable, dimension(:, :, :, :) :: ttp_all, tts_all, lam, energy
        real, allocatable, dimension(:, :) :: ttp_syn, ttp_obs, ttp_residual
        real, allocatable, dimension(:, :) :: tts_syn, tts_obs, tts_residual
        real, allocatable, dimension(:, :, :) :: vp, vs, refl, ws, wp, rs
        real, allocatable, dimension(:, :) :: wgtp, wgts
        real, dimension(1:3) :: d, o
        integer :: i, irefl
        character(len=1024) :: dir_field

        dir_field = dir_iter_record(iter)
        step_misfit = 0.0
        d = [dx, dy, dz]

        irefl = 0
        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
                case ('refl')
                    irefl = i
                    refl = model_m(i)%array
                    ! Ignore the reflector labels above nrefl
                    where (nint(refl) > nrefl)
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
            end select
        end do
        call assert(irefl > 0, ' <compute_reflector_step_shots> Error: refl is not in model_update. ')

        do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

            call set_adaptive_model_range(gmtr(ishot))
            o = [shot_xbeg, shot_ybeg, shot_zbeg]
            step_misfit(ishot) = 0

            ! The models of the shot's range: wp and ws are the velocities of
            ! the incident wave and of the other wave, rs the reflectors
            rs = refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
            if (which_medium == 'elastic-iso' .and. incident_wave == 's') then
                wp = vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
                ws = vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
            else
                wp = vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
                if (which_medium == 'elastic-iso') then
                    ws = vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
                end if
            end if

            select case (which_medium)

                case ('acoustic-iso')
                    call forward_iso_reflection(wp, d, o, gmtr(ishot), rs, ttp_all, ttp_syn)
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot), &
                        weight_used=wgtp)
                    call adjoint_iso_reflection_at_reflector(wp, d, o, gmtr(ishot), ttp_all, rs, ttp_residual, lam)
                    call adjoint_iso_reflection_at_reflector(wp, d, o, gmtr(ishot), ttp_all, rs, wgtp, energy)
                    call accumulate_reflector_step(irefl, ttp_all(:, :, :, 1), wp, rs, lam, energy)

                case ('elastic-iso')
                    ! Same-wave (PP or SS) reflections in the p or s data of the
                    ! incident wave, converted ones (PS or SP) in the other
                    select case (incident_wave)
                        case ('p')
                            call forward_iso_reflection_elastic(wp, ws, d, o, gmtr(ishot), rs, ttp_all, tts_all, ttp_syn, tts_syn)
                        case ('s')
                            call forward_iso_reflection_elastic(wp, ws, d, o, gmtr(ishot), rs, tts_all, ttp_all, tts_syn, ttp_syn)
                    end select
                    call compute_shot_misfit(ishot, dir_field, 'p', ttp_syn, ttp_obs, ttp_residual, step_misfit(ishot), &
                        weight_used=wgtp)
                    call compute_shot_misfit(ishot, dir_field, 's', tts_syn, tts_obs, tts_residual, step_misfit(ishot), &
                        weight_used=wgts)
                    select case (incident_wave)
                        case ('p')
                            call adjoint_iso_reflection_at_reflector(wp, d, o, gmtr(ishot), ttp_all, rs, ttp_residual, lam)
                            call adjoint_iso_reflection_at_reflector(wp, d, o, gmtr(ishot), ttp_all, rs, &
                                wgtp, energy)
                            call accumulate_reflector_step(irefl, ttp_all(:, :, :, 1), wp, rs, lam, energy)
                            call adjoint_iso_reflection_at_reflector(ws, d, o, gmtr(ishot), tts_all, rs, tts_residual, lam)
                            call adjoint_iso_reflection_at_reflector(ws, d, o, gmtr(ishot), tts_all, rs, &
                                wgts, energy)
                            call accumulate_reflector_step(irefl, ttp_all(:, :, :, 1), wp, rs, lam, energy, ws)
                        case ('s')
                            call adjoint_iso_reflection_at_reflector(wp, d, o, gmtr(ishot), tts_all, rs, tts_residual, lam)
                            call adjoint_iso_reflection_at_reflector(wp, d, o, gmtr(ishot), tts_all, rs, &
                                wgts, energy)
                            call accumulate_reflector_step(irefl, tts_all(:, :, :, 1), wp, rs, lam, energy)
                            call adjoint_iso_reflection_at_reflector(ws, d, o, gmtr(ishot), ttp_all, rs, ttp_residual, lam)
                            call adjoint_iso_reflection_at_reflector(ws, d, o, gmtr(ishot), ttp_all, rs, &
                                wgtp, energy)
                            call accumulate_reflector_step(irefl, tts_all(:, :, :, 1), wp, rs, lam, energy, ws)
                    end select

            end select

        end do

        call mpibarrier
        call allreduce_array(step_misfit)
        call allreduce_array(model_grad(irefl)%array)
        call allreduce_array(refl_den)

    end subroutine compute_reflector_step_shots

    !
    !> Accumulate the Gauss-Newton step of the reflector depths of one shot.
    !> At a reflector point, the residual-driven and the weight-driven
    !> receiver-side adjoint fields lam and energy give the weighted mean
    !> residual lam/energy of the rays that end there. A depth change dz
    !> changes their traveltimes by g*dz, where g is the two-way vertical
    !> slowness: g = 2*dT1/dz for the incident field T1, or the sum of the
    !> vertical slownesses of the two legs for converted waves (v2). The
    !> least-squares depth change of a column over all shots is then
    !> -sum(g*lam)/sum(g*g*energy), accumulated here into row l of the
    !> reflector gradient (numerator) and refl_den (denominator).
    !
    subroutine accumulate_reflector_step(irefl, t1, v, refl, lam, energy, v2)

        integer, intent(in) :: irefl
        real, dimension(:, :, :), intent(in) :: t1, v, refl
        real, dimension(:, :, :, :), intent(in) :: lam
        real, dimension(:, :, :, :), intent(in) :: energy
        !> Velocity of the reflected leg for converted waves; the default is v
        real, dimension(:, :, :), intent(in), optional :: v2

        integer :: i, j, k, l, i1, i2, jg, kg
        real :: g, p1z, p1x2, p2z, w2

        do l = 1, nrefl
            do k = 1, size(refl, 3)
                do j = 1, size(refl, 2)
                    do i = 1, size(refl, 1)
                        if (nint(refl(i, j, k)) == l) then
                            i1 = max(i - 1, 1)
                            i2 = min(i + 1, size(refl, 1))
                            ! Vertical slowness of the incident leg, at least that of a
                            ! 60-degree incidence; the reflected leg keeps the horizontal
                            ! slowness, so g is the two-way vertical slowness
                            p1z = max((t1(i2, j, k) - t1(i1, j, k))/((i2 - i1)*dz), 0.5/max(v(i, j, k), float_tiny))
                            if (present(v2)) then
                                w2 = max(v2(i, j, k), float_tiny)
                                p1x2 = max(1.0/max(v(i, j, k), float_tiny)**2 - p1z**2, 0.0)
                                p2z = max(sqrt(max(1.0/w2**2 - p1x2, 0.0)), 0.5/w2)
                                g = p1z + p2z
                            else
                                g = 2.0*p1z
                            end if
                            jg = shot_nybeg + j - 1
                            kg = shot_nxbeg + k - 1
                            model_grad(irefl)%array(l, jg, kg) = model_grad(irefl)%array(l, jg, kg) + g*lam(i, j, k, l)
                            refl_den(l, jg, kg) = refl_den(l, jg, kg) + g*g*energy(i, j, k, l)
                        end if
                    end do
                end do
            end do
        end do

    end subroutine accumulate_reflector_step

end submodule gradient_trtt
