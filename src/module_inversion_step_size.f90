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

module inversion_step_size

    use libflit
    use parameters
    use vars
    use utility
    use gradient

    implicit none

#ifdef dim2
#define model_dimension dimension(:, :)
#endif

#ifdef dim3
#define model_dimension dimension(:, :, :)
#endif

    private

    ! Maximum number of step size search
    ! In most cases, one (1) search is adequate
    integer :: nsearch_max = 6
    logical :: step_suitable

    public :: compute_step_size

contains

    !
    !> Check if a step size is suitable
    !
    subroutine check_step_range_single_model(step_model, &
            step_scale, srch, step_max)

        real, model_dimension, intent(in) :: srch
        real, intent(in) :: step_model, step_scale, step_max

        if (step_max == 0) then
            step_suitable = step_suitable .and. .true.
        else
            if (maxval(abs(step_model*srch*step_scale)) < step_max) then
                ! require all steps are within the prefined variation range
                step_suitable = step_suitable .and. .true.
            else
                ! require all steps are within the prefined variation range
                step_suitable = step_suitable .and. .false.
            end if
        end if

    end subroutine check_step_range_single_model

    !
    !> Compute the initial step size based on search direction values
    !
    subroutine initial_step_size(cij, srch, init_step, max_perturbation)

        ! arguments
        real, model_dimension, intent(in) :: cij, srch
        real, intent(inout) :: init_step
        real, intent(in) :: max_perturbation

        ! local variables
        real :: max_srch, max_cij

        max_cij = 0.0
        max_srch = 0.0

        ! Find out the maximum value of the search direciton
        max_srch = maxval(abs(srch))
        max_cij = maxval(abs(cij))

        ! Initial step size is maximum initial step size or 1% or maximum model value
        if (max_srch == 0) then
            init_step = 0.0
            return
        end if

        init_step = max_perturbation/max_srch

    end subroutine initial_step_size

    !
    !> Get min max of model purturbation
    !
    subroutine min_max_step_size(step, srch, step_scaling_factor, mname)

        ! Arguments
        real, model_dimension, intent(in) :: srch
        real, intent(in) :: step, step_scaling_factor
        character(len=*) :: mname

        ! Local variables
        real, allocatable, model_dimension :: deltam
        real :: min_deltam, max_deltam

        ! Allocate memory
        deltam = srch*step*step_scaling_factor

        ! Find out the maximum value of the search direciton
        min_deltam = minval(deltam)
        max_deltam = maxval(deltam)

        if (rankid == 0) then
            call warn(date_time_compact()//' Perturbation '//tidy(mname)//' range: ' &
                //num2str(min_deltam, '(es12.4)')//' ' &
                //num2str(max_deltam, '(es12.4)'))
        end if

    end subroutine min_max_step_size

    !
    !> Output updated model
    !
    subroutine output_updated_model

        integer :: i

        do i = 1, nmodel
            call output_array(model_m(i)%array, &
                dir_iter_model(iter)//'/updated_'//tidy(model_name(i))//'.bin')
        end do

    end subroutine output_updated_model

    !
    !> Check if search step size is within the desired range
    !
    subroutine check_step_range(s)

        real, intent(in) :: s
        integer :: i

        step_suitable = .true.

        do i = 1, nmodel
            call check_step_range_single_model(model_step(i), &
                s, model_srch(i)%array, model_step_max(i))
        end do

    end subroutine check_step_range

    !
    !> Backup current model
    !
    subroutine backup_current_model

        integer :: i

        do i = 1, nmodel
            model_m_backup(i)%array = model_m(i)%array
        end do

    end subroutine backup_current_model

    !
    !> Restore current model
    !
    subroutine restore_current_model

        integer :: i

        do i = 1, nmodel
            model_m(i)%array = model_m_backup(i)%array
        end do

    end subroutine restore_current_model

    !
    !> Update model
    !
    subroutine update_model(step)

        real, intent(in) :: step
        integer :: i

        do i = 1, nmodel

            model_m(i)%array = model_m(i)%array + model_step(i)*model_srch(i)%array*step

            ! Box-clip model
            if (model_m(i)%name == 'vp' .or. model_m(i)%name == 'vs') then
                where (model_m(i)%array /= 0)
                    model_m(i)%array = clip(model_m(i)%array, model_min(i), model_max(i))
                end where
            else
                model_m(i)%array = clip(model_m(i)%array, model_min(i), model_max(i))
            end if

            ! For TRTT, the search direction is the reflector's image
            if (model_m(i)%name == 'refl') then
                model_m(i)%array = model_srch(i)%array
            end if

        end do

        ! Enforce Vp/Vs ratio in an appropriate range
        call clip_vpvsratio

    end subroutine update_model

    !
    !> Print perturbation information
    !
    subroutine print_perturb_info

        integer :: i

        do i = 1, nmodel
            call min_max_step_size(model_step(i), model_srch(i)%array, step_scaling_factor, model_name(i))
        end do

    end subroutine print_perturb_info

    !
    !> Calcualte initial step size for all parametes
    !
    !> Note that some parameters may be zero at the very beginning,
    !> then the initial step size is set to step_max_cij
    !> under such circumstances
    !
    subroutine compute_initial_step_size

        integer :: i
        real :: valstep

        do i = 1, nmodel

            select case (remove_string_after(model_name(i), ['c', 'C']))
                case ('vp', 'vs', 'rho')
                    valstep = 100.0
                case ('epsilon', 'delta', 'gamma', 'eps', 'del', 'gam')
                    valstep = 0.1
                case ('c', 'C')
                    valstep = 1.0e9
                case ('refl')
                    valstep = 1.0
                case ('sx')
                    valstep = 0.1*(nx - 1)*dx
                case ('sy')
                    valstep = 0.1*(ny - 1)*dy
                case ('sz')
                    valstep = 0.1*(nz - 1)*dz
                case ('st0')
                    valstep = 0.1
            end select

            call readpar_xfloat(file_parameter, 'step_max_'//tidy(model_name(i)), model_step_max(i), valstep, iter*1.0)
            ! In most cases, step_max_scale_factor = 1
            model_step_max(i) = model_step_max(i)*step_max_scale_factor
            call initial_step_size(model_m(i)%array, model_srch(i)%array, model_step(i), model_step_max(i))

        end do

    end subroutine compute_initial_step_size

    !
    !> Calculate optimal step size coefficient for some parameter
    !
    subroutine compute_step_coef_single_component(srcindex, component_name, sum1, sum2)

        ! Arguments
        integer, intent(in) :: srcindex
        character(len=*) :: component_name
        real, intent(inout) :: sum1, sum2

        ! Local variables
        character(len=1024) :: file_recorded
        character(len=1024) :: file_synthetic
        character(len=1024) :: file_synthetic_prev
        character(len=256) :: filename
        integer :: nr, i, j
        real, allocatable, dimension(:, :) :: seis_obs, seis_syn, seis_syn_prev
        real, allocatable, dimension(:, :) :: weight
        real, allocatable, dimension(:, :) :: r1, r2, m1, m2
        real, allocatable, dimension(:, :) :: tobs_all, tsyn_all, tsyn_all_prev

        if (yn_dd_no_st0) then

            ! Setup filename prefix and filenames
            filename = '/t'//tidy(component_name)//'_all.bin'

            ! Observed data
            file_recorded = dir_iter_record(iter)//tidy(filename)

            ! Synthetic data last iteration
            file_synthetic_prev = dir_iter_synthetic(iter - 1)//tidy(filename)

            ! Synthetic data current iteration
            file_synthetic = dir_iter_synthetic(iter)//tidy(filename)

            tobs_all = load(file_recorded, nr_virtual, ns)
            tsyn_all_prev = load(file_synthetic_prev, nr_virtual, ns)
            tsyn_all = load(file_synthetic, nr_virtual, ns)

            ! For now only considers transmission data (no reflection)
            nr = gmtr(srcindex)%nr
            r1 = zeros(nr, 1)
            r2 = zeros(nr, 1)
            !$omp parallel do private(i, j)
            do j = 1, nr
                if (maxval(tsyn_all(j, :)) > 0 .and. maxval(tsyn_all_prev(j, :)) > 0 .and. maxval(tobs_all(j, :)) > 0) then
                    do i = 1, ns
                        if (tobs_all(j, srcindex) >= 0 .and. tobs_all(j, i) >= 0) then
                            r1(j, 1) = r1(j, 1) + (tsyn_all_prev(j, srcindex) - tsyn_all_prev(j, i)) - (tsyn_all(j, srcindex) - tsyn_all(j, i))
                            r2(j, 1) = r2(j, 1) + (tsyn_all_prev(j, srcindex) - tsyn_all_prev(j, i)) - (tobs_all(j, srcindex) - tobs_all(j, i))
                        end if
                    end do
                end if
            end do
            !$omp end parallel do

        else

            ! Setup filename prefix and filenames
            filename = '/shot_'//num2str(gmtr(srcindex)%id)//'_'//'traveltime_'//tidy(component_name)//'.bin'

            ! Observed data
            file_recorded = dir_iter_record(iter)//tidy(filename)

            ! Synthetic data last iteration
            file_synthetic_prev = dir_iter_synthetic(iter - 1)//tidy(filename)

            ! Synthetic data current iteration
            file_synthetic = dir_iter_synthetic(iter)//tidy(filename)

            ! Read data
            nr = gmtr(srcindex)%nr
            seis_obs = load(file_recorded, nr, nrefl + 1)
            seis_syn = load(file_synthetic, nr, nrefl + 1)
            seis_syn_prev = load(file_synthetic_prev, nr, nrefl + 1)

            call traveltime_residual(seis_syn, seis_syn_prev, r1, m1)
            call traveltime_residual(seis_obs, seis_syn_prev, r2, m2)

        end if

        weight = zeros(nr, nrefl + 1)
        !$omp parallel do private(i, j)
        do i = 1, nr
            do j = 1, nrefl + 1

                weight(i, j) = gmtr(srcindex)%recr(i)%weight*misfit_weight(j)

                if (abs(r1(i, j)) > misfit_threshold .or. abs(r2(i, j)) > misfit_threshold) then
                    weight(i, j) = 0
                end if

                ! For reflection data
                if (j > 1 .and. (gmtr(srcindex)%recr(i)%aoff < offset_min_refl &
                        .or. gmtr(srcindex)%recr(i)%aoff > offset_max_refl)) then
                    weight(i, j) = 0
                end if

            end do
        end do
        !$omp end parallel do

        sum1 = sum1 + sum(abs(r1)*abs(r2)*weight)
        sum2 = sum2 + sum(abs(r1)**2*weight)

    end subroutine compute_step_coef_single_component

    !
    !> Calculate optimal step size coefficient for all parameters to be updated
    !
    subroutine compute_step_coef(srcindex, sum1, sum2)

        integer, intent(in) :: srcindex
        real, intent(inout) :: sum1, sum2

        integer :: i

        do i = 1, ndata
            call compute_step_coef_single_component(srcindex, data_name(i), sum1, sum2)
        end do

    end subroutine compute_step_coef

    !
    !> Compute step misift associated with an updated model
    !
    subroutine compute_step_misfit(step_scaling_factor, misfit)

        real, intent(in) :: step_scaling_factor
        real, intent(inout) :: misfit

        ! update model
        call backup_current_model
        call update_model(step_scaling_factor)

        ! forward modeling
        yn_misfit_only = .true.
        call compute_gradient_shots
        call mpibarrier

        misfit = sum(step_misfit)

        ! backdate model
        call restore_current_model
        call mpibarrier

    end subroutine compute_step_misfit

    !
    !> Compute the optimal step size for a quasi-linear inversion
    !> The part mostly stemmed from the discussion with Dr. Benxin Chi.
    !
    subroutine compute_step_size_linear

        ! Local variables
        real :: trial_misfit, sum1, sum2, data_misfit_current
        integer :: cnt, i
        character(len=1024) :: dir_from, dir_to

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Computing optimal step size...')
        end if

        ! Calculate the initial step size
        call compute_initial_step_size

        ! default initial step scalar = 0.1
        step_scaling_factor = 0.2

        ! check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact()//' Error: Cannot find a suitable initial step size. Exit. ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        ! Step size coefficients
        sum1 = 0.0
        sum2 = 0.0

        ! Update model
        call backup_current_model
        call update_model(step_scaling_factor)

        ! forward modeling and compute misfit
        yn_misfit_only = .true.

        iter = iter + 1
        put_synthetic_in_scratch = .false.
        call make_iter_dir
        call compute_gradient_shots

        call mpibarrier
        do i = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)
            call compute_step_coef(i, sum1, sum2)
        end do

        iter = iter - 1
        call mpibarrier

        ! Restore current iteration model
        call restore_current_model

        ! Calculate optimal step size
        call mpibarrier
        call allreduce(sum1)
        call allreduce(sum2)
        if (isnan(sum1) .or. isnan(sum2)) then
            call warn(date_time_compact()//' Step size coef sum 1 = '//num2str(sum1, '(es)'))
            call warn(date_time_compact()//' Step size coef sum 2 = '//num2str(sum2, '(es)'))
            call mpibarrier
            call mpiend
        else
            step_scaling_factor = sum1/(sum2 + float_tiny)
        end if
        if (rankid == 0) then
            call warn(date_time_compact()//' Step size coef sum 1 = '//num2str(sum1, '(es)'))
            call warn(date_time_compact()//' Step size coef sum 2 = '//num2str(sum2, '(es)'))
            call warn(date_time_compact()//' Step size coef = '//num2str(step_scaling_factor, '(es)'))
        end if

        if (.not. (step_scaling_factor <= float_huge)) then
            if (rankid == 0) then
                call warn(date_time_compact()//' <compute_step_size_linear> Error: Step size is NaN. Exit. ')
                call mpi_abort(mpi_comm_world, mpi_err_other, mpi_ierr)
            end if
        end if

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size scaling factor is computed. ')
        end if

        ! Check step size
        step_scaling_factor = step_scaling_factor*2.0
        data_misfit_current = data_misfit(iter)
        trial_misfit = data_misfit_current*2.0
        cnt = 1

        ! Enforce update or not
        call readpar_xlogical(file_parameter, 'yn_enforce_update', yn_enforce_update, .false., iter*1.0)

        if (yn_enforce_update) then
            ! To enforce update regardless of misfit increase,
            ! use only one trial. This can be enabled after certain number of
            ! iterations, not necessarily from the begining

            do while (cnt < 2)

                step_scaling_factor = step_scaling_factor/2.0

                call check_step_range(step_scaling_factor)
                if (.not. step_suitable) then
                    cycle
                end if

                call compute_step_misfit(step_scaling_factor, trial_misfit)

                call print_step_info(cnt, step_scaling_factor, trial_misfit)

                cnt = cnt + 1

            end do

        else
            ! Otherwise, use conventional trial-error approach

            do while (trial_misfit > data_misfit_current .and. cnt < nsearch_max)

                step_scaling_factor = step_scaling_factor/2.0

                ! make sure step size and model values in suitable range
                ! if not in range, then half until suitable
                call check_step_range(step_scaling_factor)
                if (.not. step_suitable) then
                    cycle
                end if

                put_synthetic_in_scratch = .true.
                call compute_step_misfit(step_scaling_factor, trial_misfit)

                ! Print step size information
                call print_step_info(cnt, step_scaling_factor, trial_misfit)

                cnt = cnt + 1

            end do

            if (trial_misfit > data_misfit_current) then
                ! When searching the optimal step size, if the trial step size
                ! at the final trial cannot produce a smaller data misfit than
                ! that of the last iteration, then set step size to zero
                ! and the trial misfit to that in the last iteration
                step_scaling_factor = 0.0
                trial_misfit = data_misfit_current
            else
                ! Otherwise, copy the synthetic and synthetic processed data residing in scratch
                ! to the current synthetic and synthetic processed data
                ! so that they are the synthetic data of the current iteration
                if (rankid == 0) then
                    ! This only executed by rank0
                    put_synthetic_in_scratch = .true.
                    dir_from = dir_iter_synthetic(iter)
                    put_synthetic_in_scratch = .false.
                    dir_to = dir_iter_synthetic(iter)
                    call delete_directory(dir_to)
                    call move_directory(dir_from, dir_to)
                end if
            end if

        end if

        ! print perturbation infomation
        call print_perturb_info

        ! Current absolute data misfit
        data_misfit(iter) = trial_misfit

        ! Make all process ready
        call mpibarrier

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation is done.')
        end if

        ! update model and output updated model for current iteration
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if
        put_synthetic_in_scratch = .false.

    end subroutine compute_step_size_linear

    !
    !> Calculate the optimal step size to update the model
    !> with a quadratic + interval bounding + bisection line search method
    !
    subroutine compute_step_size_linesearch

        integer :: cnt
        real :: al, ar, ac, fl, fr, fc, a0
        real :: acr, arl, alc, bcr, brl, blc, num, den, an, fn, a2, f2
        !        real :: step0,
        real :: misfit0, step1, misfit1, step2, misfit2
        character(len=1024) :: dir_from, dir_to

        ! Print info
        if (rankid == 0) then
            call warn(' >>>>>>>>>> Computing step size... ')
        end if

        ! Calculate the initial step size
        call compute_initial_step_size

        ! default initial step scalar = 0.1
        step_scaling_factor = 1.0

        ! check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' Error: Cannot find a suitable initial step size. Exit. ')
                ! call warn(date_time_compact()//' one of more parameters have reached their assigned limits ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        a0 = step_scaling_factor

        ! The directory of synthetic
        put_synthetic_in_scratch = .true.

        ! Prepare the initial values
        al = 0.0
        fl = data_misfit(iter)
        misfit0 = fl
        call print_step_info(0, al, fl)
        call mpibarrier

        ar = a0
        call compute_step_misfit(ar, fr)
        call print_step_info(0, ar, fr)
        step1 = ar
        misfit1 = fr
        call mpibarrier

        ! If this step produces a smaller misfit, then skip trial
        if (fr < fl) then
            ac = ar
            fc = fr
            call mpibarrier
            goto 123
        end if

        ac = 0.5*(al + ar)
        call compute_step_misfit(ac, fc)
        call print_step_info(0, ac, fc)
        step2 = ac
        misfit2 = fc
        call mpibarrier

        ! If this step produces a smaller misfit, then skip trial
        if (fc < fl) then
            call mpibarrier
            goto 123
        end if

        cnt = 3

        call mpibarrier

        ! Quadratic fit + bisection line search
        do while (cnt < nsearch_max)

            if ((fc < fl) .and. (fc < fr)) then

                acr = ac - ar
                bcr = ac**2 - ar**2

                arl = ar - al
                brl = ar**2 - al**2

                alc = al - ac
                blc = al**2 - ac**2

                num = bcr*fl + brl*fc + blc*fr
                den = acr*fl + arl*fc + alc*fr

                if (den == 0.0) then
                    exit
                end if

                an = 0.5*num/den
                call compute_step_misfit(an, fn)
                call print_step_info(cnt - 2, an, fn)
                call mpibarrier

                cnt = cnt + 1

                if (an > ac) then
                    if (fn >= fc) then
                        ar = an
                        fr = fn
                    else
                        al = ac
                        fl = fc
                        ac = an
                        fc = fn
                    end if
                else
                    if (fn >= fc) then
                        al = an
                        fl = fn
                    else
                        ar = ac
                        fr = fc
                        ac = an
                        fc = fn
                    end if
                end if

            else
                ! Bisection safe switchover

                a2 = ac + 1.0e-1*(ar - al)
                call compute_step_misfit(a2, f2)
                call print_step_info(cnt - 2, a2, f2)
                call mpibarrier

                cnt = cnt + 1

                if (fc < f2) then
                    ar = a2
                    fr = f2
                else
                    al = ac
                    fl = fc
                end if

                ac = 0.5*(al + ar)
                call compute_step_misfit(ac, fc)
                call print_step_info(cnt - 2, ac, fc)
                call mpibarrier

                cnt = cnt + 1

            end if

            ! Stop search if range too small
            if (abs(ar - al) <= 0.05) then
                if (fc >= misfit1) then
                    ac = step1
                    fc = misfit1
                end if
                if (fc >= misfit2) then
                    ac = step2
                    fc = misfit2
                end if
                exit
            end if

        end do

        123 continue

        if (fc >= misfit0) then
            ! When searching optimal step size, if the trial step size
            ! at the final trial cannot produce a smaller data misfit than
            ! that of the last iteration, then set step size to zero
            ! and the trial misfit to that in the last iteration
            step_scaling_factor = 0.0
            data_misfit(iter) = misfit0
        else
            ! Otherwise, copy the synthetic and synthetic processed data residing in scratch
            ! to the current synthetic and synthetic processed data
            ! so that they are the synthetic data of the current iteration
            step_scaling_factor = ac
            data_misfit(iter) = fc
            if (rankid == 0) then
                ! This only executed by rank0
                put_synthetic_in_scratch = .true.
                dir_from = dir_iter_synthetic(iter)
                put_synthetic_in_scratch = .false.
                dir_to = dir_iter_synthetic(iter)
                call delete_directory(dir_to)
                call move_directory(dir_from, dir_to)
            end if
        end if

        ! Make all process ready
        call mpibarrier

        ! Print info
        !    call date_time_compact(string_date_time)
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation is done. ')
        end if

        ! update model and output
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if

        put_synthetic_in_scratch = .false.

    end subroutine compute_step_size_linesearch

    !
    !> Calculate the optimal step size to update the model
    !> with a quadratic search method
    !
    subroutine compute_step_size_quadratic

        integer :: cnt
        real :: step0, misfit0, step1, misfit1, step2, misfit2
        character(len=1024) :: dir_from, dir_to

        ! Print info
        if (rankid == 0) then
            call warn(' >>>>>>>>>> Computing step size...')
        end if

        ! Calculate the initial step size
        call compute_initial_step_size

        ! default initial step scalar = 0.1
        step_scaling_factor = 1.0

        ! check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact()//' Error: Cannot find a suitable initial step size. Exit. ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        ! The directory of synthetic
        put_synthetic_in_scratch = .true.

        step0 = 0.0
        misfit0 = data_misfit(iter)
        call print_step_info(0, step0, misfit0)
        call mpibarrier

        step1 = 0.1*step_scaling_factor
        call compute_step_misfit(step1, misfit1)
        call print_step_info(1, step1, misfit1)
        call mpibarrier

        step2 = step_scaling_factor
        call compute_step_misfit(step2, misfit2)
        call print_step_info(2, step2, misfit2)
        call mpibarrier

        step_scaling_factor = &
            0.5*((misfit0 - misfit2)*step1**2 + (-misfit0 + misfit1)*step2**2) &
            /(-(misfit2*step1) + misfit0*(step1 - step2) + misfit1*step2)
        call compute_step_misfit(step_scaling_factor, data_misfit(iter))
        call print_step_info(3, step_scaling_factor, data_misfit(iter))
        call mpibarrier

        if (rankid == 0) then
            ! This only executed by rank0
            put_synthetic_in_scratch = .true.
            dir_from = dir_iter_synthetic(iter)
            put_synthetic_in_scratch = .false.
            dir_to = dir_iter_synthetic(iter)
            call delete_directory(dir_to)
            call move_directory(dir_from, dir_to)
        end if

        ! Make all process ready
        call mpibarrier

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation is done. ')
        end if

        ! update model and output
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if

        put_synthetic_in_scratch = .false.

    end subroutine compute_step_size_quadratic

    !
    !> Compute step size
    !
    subroutine compute_step_size

        select case (step_size_method)
            case ('linear')
                call compute_step_size_linear
            case ('quadratic')
                call compute_step_size_quadratic
            case ('line_search')
                call compute_step_size_linesearch
        end select

        call mpibarrier

    end subroutine compute_step_size

end module inversion_step_size
