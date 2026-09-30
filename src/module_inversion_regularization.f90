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


!
!> Regularization of the inversion by alternating optimization:
!>
!>     m^(l+1)  = argmin_m  { J(m) + lambda_m ||m - mu^(l)||_2^2 },
!>     mu^(l+1) = argmin_mu { T(mu) + lambda_mu ||m^(l+1) - mu||_2^2 },
!>
!> where l is the iteration index, m is the model, J is the data misfit, mu is the
!> regularized model (model_reg), T is the regularization functional of the listed
!> methods, and lambda_m and lambda_mu are regularization weights. Each iteration
!> takes one descent step of the first problem, with the correction
!> lambda_m*(m - mu^(l)) added to the gradient (add_l2reg), and then solves the
!> second problem by denoising the updated model (compute_regularization).
!>
!> The strength of a parameter (regularization_strength) is reg_scale_<name>, or
!> reg_lambda_<name> when const_reg is true, and falls back to reg_scale or
!> reg_lambda. A zero strength switches off the regularization of the parameter.
!>
!> In LATTE, the traveltime computation is fast, while the denoising (TGpV, and
!> especially ML-based source regularization) can be slow, and it is rarely needed
!> in early iterations. Therefore, compute_regularization denoises only when the next
!> iteration, which uses the result, has a positive strength. An iteration without an
!> up-to-date mu (the first iteration of a run, or the first one after a zero strength)
!> adds no correction. add_l2reg must only read model_reg: if it stored m - mu
!> there, an iteration after a skipped denoising would use that difference as mu.
!>
!> To save the cost in early iterations, a schedule must stay exactly zero there,
!> for instance reg_scale_vp = 1~20:0, 21:0.4. Schedules are interpolated linearly,
!> so a ramp such as reg_scale_vp = 1:0, 30:0.5 is positive from iteration 2 on.
!>
!> These cost considerations are the reasons why LATTE differs here from OWL
!> (https://github.com/lanl/owl), which uses the same scheme. In OWL, the
!> wave-equation modeling costs far more than the denoising, so OWL denoises in
!> every iteration once a method is listed. Its mu is therefore always up to date,
!> at the cost of denoising while the strength is zero. OWL also defaults the
!> strength to 0.2, while LATTE requires an explicit strength.
!
module inversion_regularization

    use parameters
    use regularization

#ifdef dim2
#define model_dimension dimension(:, :)
#endif

#ifdef dim3
#define model_dimension dimension(:, :, :)
#endif

    implicit none

    ! Whether compute_regularization refreshed the regularized models (model_reg)
    ! in the last iteration: for the medium parameters vp and vs, and for the
    ! source positions sx, sy and sz
    logical :: model_reg_ready = .false.
    logical :: source_reg_ready = .false.

contains

    !
    !> Regularization strength of a parameter in an iteration, by default the current one:
    !> reg_lambda_<name> when const_reg is true, and reg_scale_<name> otherwise
    !
    function regularization_strength(name, iteration) result(s)

        character(len=*), intent(in) :: name
        integer, intent(in), optional :: iteration
        real :: s

        logical :: const_reg
        real :: s0, it

        if (present(iteration)) then
            it = iteration*1.0
        else
            it = iter*1.0
        end if

        ! A parameter without its own value takes the value for all parameters
        call readpar_xlogical(file_parameter, 'const_reg', const_reg, .false., it)
        if (const_reg) then
            call readpar_xfloat(file_parameter, 'reg_lambda', s0, 0.0, it)
            call readpar_xfloat(file_parameter, 'reg_lambda_'//tidy(name), s, s0, it)
        else
            call readpar_xfloat(file_parameter, 'reg_scale', s0, 0.0, it)
            call readpar_xfloat(file_parameter, 'reg_scale_'//tidy(name), s, s0, it)
        end if

    end function regularization_strength

    !
    !> Add L2-norm regularization term to gradient for a single parameter
    !
    subroutine add_l2reg_single_parameter(reg, model, grad, name)

        real, model_dimension, intent(in) :: model, reg
        real, model_dimension, intent(inout) :: grad
        character(len=*), intent(in) :: name
        character(len=1024) :: file_mask

        real :: reg_coef, reg_scale
        logical :: const_reg
        real, allocatable, model_dimension :: r, grad_mask
        character(len=64) :: process_name
        character(len=32), allocatable, dimension(:) :: process_list

        ! Difference between the model and its regularized version; reg keeps
        ! the regularized version, which only compute_regularization updates
        r = model - reg

        ! Set regularization weights
        call readpar_xlogical(file_parameter, 'const_reg', const_reg, .false., iter*1.0)
        if (const_reg) then

            reg_coef = regularization_strength(name)

        else

            reg_scale = regularization_strength(name)
            if (rankid == 0) then
                call warn(date_time_compact()//' Regularization scale for '//tidy(name)//' = '//num2str(reg_scale, '(es)'))
            end if

            ! Compute regularization weight for each parameter
            if (maxval(abs(r)) == 0) then
                reg_coef = 0.0
            else
                reg_coef = mean(grad, 2)/mean(r, 2)*reg_scale
            end if

        end if

        if (rankid == 0) then
            call warn(date_time_compact()//' Regularization coefficient for '//tidy(name)//' = '//num2str(reg_coef, '(es)'))
        end if

        ! modify gradients
        grad = grad + reg_coef*r

        ! Mask gradient again, with the mask that gradient processing applied to
        ! this parameter; source parameters are neither processed nor masked
        if (.not. any(name == ['sx', 'sy', 'sz', 'st0'])) then
            if (uniform_processing) then
                process_name = 'grad'
            else
                process_name = 'grad_'//tidy(name)
            end if
            call readpar_nstring(file_parameter, 'process_'//tidy(process_name), process_list, [''])
            if (any(process_list == 'mask')) then
                call readpar_xstring(file_parameter, tidy(process_name)//'_mask', file_mask, '', iter*1.0)
                if (file_mask /= '') then
                    call prepare_model_single_parameter(grad_mask, 'mask', file_mask, update=.false.)
                else
                    grad_mask = ones_like(grad)
                end if 
                grad = grad*grad_mask
            end if
        end if

    end subroutine add_l2reg_single_parameter

    !
    !> Add regularization
    !
    subroutine add_l2reg

        integer :: i

        do i = 1, nmodel
            ! Only parameters whose regularized versions are up to date; no source
            ! regularization method regularizes st0
            if ((any(model_name(i) == ['vp', 'vs']) .and. model_reg_ready) &
                    .or. (any(model_name(i) == ['sx', 'sy', 'sz']) .and. source_reg_ready)) then
                call add_l2reg_single_parameter(model_reg(i)%array, model_m(i)%array, &
                    model_grad(i)%array, model_name(i))
            end if
        end do

        call mpibarrier

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> L2 regularization finished ')
        end if

    end subroutine add_l2reg

    !
    !> Regularize gradient
    !
    subroutine regularize_gradient

        if (model_reg_ready .or. source_reg_ready) then
            call add_l2reg
        end if

    end subroutine regularize_gradient

    !
    !> Update regularization variables
    !
    subroutine compute_regularization

        integer :: i

        ! The regularized models computed here enter the gradient of the next
        ! iteration, so the strengths of the next iteration decide whether to compute them

        ! For tomography, regularize the medium parameters when vp or vs
        ! has a positive regularization strength in the next iteration
        model_reg_ready = .false.
        if (yn_regularize_model) then
            do i = 1, nmodel
                if (any(model_name(i) == ['vp', 'vs'])) then
                    if (regularization_strength(model_name(i), iter + 1) > 0) then
                        model_reg_ready = .true.
                        exit
                    end if
                end if
            end do
            if (model_reg_ready) then
                call model_regularization
            end if
        end if

        ! For location, regularize the source positions when sx, sy or sz
        ! has a positive regularization strength in the next iteration
        source_reg_ready = .false.
        if (yn_regularize_source) then
            do i = 1, nmodel
                if (any(model_name(i) == ['sx', 'sy', 'sz'])) then
                    if (regularization_strength(model_name(i), iter + 1) > 0) then
                        source_reg_ready = .true.
                        exit
                    end if
                end if
            end do
            if (source_reg_ready) then
                call source_regularization
            end if
        end if

    end subroutine compute_regularization

end module inversion_regularization
