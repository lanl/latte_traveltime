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


submodule(gradient) gradient_fatt

    use traveltime_iso

contains

    !
    !> Compute gradients shot by shot for FATT
    !
    module subroutine compute_gradient_shots_fatt

        real, allocatable, dimension(:, :) :: ttp, tts
        real, allocatable, dimension(:, :) :: ttp_syn, ttp_obs, ttp_residual
        real, allocatable, dimension(:, :) :: tts_syn, tts_obs, tts_residual
        real, allocatable, dimension(:, :) :: lam, energy
        integer :: i
        real, allocatable, dimension(:, :) :: vp, vs
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

        ! step misfit set to zero
        step_misfit = 0.0d0

        ! Copy models
        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
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
                    select case (forward_eikonal_method)
                        case ('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                        case ('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                    end select
                    call output_array(ttp_syn, tidy(dir_synthetic)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')

                case ('elastic-iso')
                    select case (forward_eikonal_method)
                        case ('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                            call forward_iso_fast_march( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, tts_syn)
                        case ('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttp_syn)
                            call forward_iso_fast_sweep( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, tts_syn)
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
                    call adjoint_iso( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, &
                        ttp_residual, lam)

                    if (yn_precond) then

                        ! Compute preconditioner
                        call adjoint_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, &
                            ones_like(ttp_residual), energy)

                        lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                    else

                        lam = -lam/vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)**3

                    end if

                    if (uniform_processing) then
                        call process_model_single_shot(ishot, model_grad(1)%array, lam, 'grad')
                    else
                        call process_model_single_shot(ishot, model_grad(1)%array, lam, 'grad_'//model_name(1))
                    end if

                    call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                        //' gradient computation is done. ')

                case ('elastic-iso')

                    ! For vp
                    ! Compute adjoint-state field
                    call adjoint_iso( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, &
                        ttp_residual, lam)

                    if (yn_precond) then

                        ! Compute preconditioner
                        call adjoint_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, &
                            ones_like(ttp_residual), energy)

                        lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                    else

                        lam = -lam/vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)**3

                    end if

                    do i = 1, nmodel
                        if (model_name(i) == 'vp') then
                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(i))
                            end if
                        end if
                    end do

                    ! For vs
                    ! Compute adjoint-state field
                    call adjoint_iso( &
                        vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, &
                        tts_residual, lam)

                    if (yn_precond) then

                        ! Compute preconditioner
                        call adjoint_iso( &
                            vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, &
                            ones_like(tts_residual), energy)

                        lam = -lam/(energy + precond_eps*maxval(abs(energy)))

                    else

                        lam = -lam/vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)**3

                    end if

                    do i = 1, nmodel
                        if (model_name(i) == 'vs') then
                            if (uniform_processing) then
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad')
                            else
                                call process_model_single_shot(ishot, model_grad(i)%array, lam, 'grad_'//model_name(i))
                            end if
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

    end subroutine compute_gradient_shots_fatt

end submodule gradient_fatt
