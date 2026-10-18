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


module inversion_gradient

    use parameters
    use vars
    use gradient
    use inversion_regularization

    implicit none

contains

    !
    !> Compute gradient by adjoint-state method
    !
    subroutine compute_gradient

        ! Initialization
        call zero_gradient

        ! Compute gradients
        yn_misfit_only = .false.
        call compute_gradient_shots
        call mpibarrier

        ! Process and regularize gradients
        call process_gradient
        call regularize_gradient
        call output_gradient
        call mpibarrier

        ! Print progress
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Gradient computation completed. ')
        end if

    end subroutine compute_gradient

    !
    !> The data misfit of the current models, without a model update and
    !> without touching the misfit history
    !
    function misfit_of_current_model() result(misfit)

        real :: misfit

        logical :: scratch_saved
        real :: misfit_saved

        scratch_saved = put_synthetic_in_scratch
        misfit_saved = data_misfit(iter)
        put_synthetic_in_scratch = .true.
        yn_misfit_only = .true.
        call compute_gradient_shots
        call mpibarrier
        misfit = sum(step_misfit)
        yn_misfit_only = .false.
        put_synthetic_in_scratch = scratch_saved
        data_misfit(iter) = misfit_saved

    end function misfit_of_current_model

    !
    !> Update the reflectors for the current velocity models. Image the
    !> reflectors in the iterations where yn_refl_image is on. The image
    !> replaces a nonzero reflector model only when it does not raise the data
    !> misfit. Then refine the reflector depths by least squares
    !> (refine_reflector).
    !
    subroutine update_reflector_image

        integer :: i, irefl, l
        logical :: yn_image
        real :: misfit_new, misfit_prev
#ifdef dim2
        real, allocatable, dimension(:, :) :: refl_prev, refl_new
#endif
#ifdef dim3
        real, allocatable, dimension(:, :, :) :: refl_prev, refl_new
#endif

        irefl = 0
        do i = 1, nmodel
            if (model_m(i)%name == 'refl') then
                irefl = i
            end if
        end do
        call assert(irefl > 0, ' <update_reflector_image> Error: refl is not in model_update. ')

        ! The least-squares refinement follows the velocity models, so imaging
        ! in every iteration is optional. Without imaging, the refinement starts
        ! from the reflector model in use, which must then hold every reflector.
        call readpar_xlogical(file_parameter, 'yn_refl_image', yn_image, .true., iter*1.0)
        if (.not. yn_image) then
            do l = 1, nrefl
                call assert(any(nint(model_m(irefl)%array) == l), ' <update_reflector_image> Error: the reflector model ' &
                    //'does not have reflector '//num2str(l)//'; image the reflectors with yn_refl_image. ')
            end do
            call refine_reflector(irefl)
            return
        end if

        refl_prev = model_m(irefl)%array

        ! The imaging and the guard use the reflections up to offset_max_refl_depth
        yn_refl_depth_offsets = .true.
        call compute_image_shots
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Reflector image computation completed. ')
        end if

        ! The image replaces a nonzero reflector model only when it does not
        ! raise the data misfit
        if (maxval(abs(refl_prev)) > 0) then

            refl_new = model_m(irefl)%array
            misfit_new = misfit_of_current_model()
            if (iter >= 2 .and. .not. refl_depth_offsets_differ()) then
                ! The misfit of the models in use, from the end of the last iteration
                misfit_prev = data_misfit(iter - 1)
            else
                model_m(irefl)%array = refl_prev
                misfit_prev = misfit_of_current_model()
            end if

            if (misfit_new > misfit_prev) then
                model_m(irefl)%array = refl_prev
                if (rankid == 0) then
                    call warn(date_time_compact()//' >>>>>>>>>> Reflector image rejected: misfit ' &
                        //num2str(misfit_new, '(es)')//' > '//num2str(misfit_prev, '(es)')//'; the previous reflector model stays. ')
                end if
            else
                model_m(irefl)%array = refl_new
                if (rankid == 0) then
                    call warn(date_time_compact()//' >>>>>>>>>> Reflector image accepted: misfit ' &
                        //num2str(misfit_new, '(es)')//' <= '//num2str(misfit_prev, '(es)'))
                end if
            end if

        end if
        yn_refl_depth_offsets = .false.
        call mpibarrier

        ! Least-squares refinement of the reflector depths at the current velocity models
        call refine_reflector(irefl)

    end subroutine update_reflector_image

    !
    !> Whether the reflector depths use fewer reflections than the velocity
    !> inversion, i.e., offset_max_refl_depth < offset_max_refl
    !
    logical function refl_depth_offsets_differ()
        refl_depth_offsets_differ = offset_max_refl_depth < offset_max_refl
    end function refl_depth_offsets_differ

    !
    !> Refine the reflector depths at the current velocity models by
    !> Gauss-Newton steps (refl_refine_niter). A step moves each column of each
    !> reflector by the depth change that zeros, in the least-squares sense,
    !> the residuals of the reflections there, capped at refl_refine_step_max.
    !> The step is halved, at most twice, until it lowers the data misfit, and
    !> dropped otherwise, which ends the refinement. The imaging condition is
    !> an envelope over the receivers and the shots, which places the
    !> reflectors too deep when the picks are noisy. The refinement removes
    !> this bias and keeps the reflectors at the depths that fit the data best
    !> in the current velocity models, so that the velocity gradient does not
    !> compensate for a misplaced reflector.
    !
    subroutine refine_reflector(irefl)

        integer, intent(in) :: irefl

        integer :: k, ntry, niter_refine
        real :: misfit0, misfit1, scale, step_max
        logical :: accepted
#ifdef dim2
        real, allocatable, dimension(:, :) :: refl0, dzr
#endif
#ifdef dim3
        real, allocatable, dimension(:, :, :) :: refl0, dzr
#endif

        call readpar_xint(file_parameter, 'refl_refine_niter', niter_refine, 2, iter*1.0)
        call readpar_xfloat(file_parameter, 'refl_refine_step_max', step_max, 5*dz, iter*1.0)

        ! The steps and their acceptance use the reflections in the offset
        ! range of the reflector depths
        yn_refl_depth_offsets = .true.

        do k = 1, niter_refine

            call zero_reflector_step(irefl)
            call compute_reflector_step_shots
            misfit0 = sum(step_misfit)
            call reflector_step(irefl, step_max, dzr)
            refl0 = model_m(irefl)%array

            accepted = .false.
            scale = 1.0
            do ntry = 1, 3
                model_m(irefl)%array = refl0
                call move_reflector(model_m(irefl)%array, scale*dzr)
                misfit1 = misfit_of_current_model()
                if (misfit1 < misfit0) then
                    accepted = .true.
                    exit
                end if
                scale = 0.5*scale
            end do

            if (accepted) then
                if (rankid == 0) then
                    call warn(date_time_compact()//' >>>>>>>>>> Reflector refinement accepted: misfit ' &
                        //num2str(misfit1, '(es)')//' < '//num2str(misfit0, '(es)') &
                        //', largest depth change '//num2str(scale*maxval(abs(dzr)), '(es)'))
                end if
            else
                model_m(irefl)%array = refl0
                if (rankid == 0) then
                    call warn(date_time_compact()//' >>>>>>>>>> Reflector refinement rejected: misfit ' &
                        //num2str(misfit1, '(es)')//' >= '//num2str(misfit0, '(es)'))
                end if
                exit
            end if

        end do
        yn_refl_depth_offsets = .false.
        call mpibarrier

    end subroutine refine_reflector

end module inversion_gradient
