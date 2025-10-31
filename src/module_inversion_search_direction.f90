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
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module inversion_search_direction

    use libflit
    use parameters
    use utility
    use vars

    implicit none

    ! Set L-BFGS maximum layers
    integer :: nstep_lbfgs = 10
    integer :: ng

contains

    !
    !> Calculate search direction in Conjugate-Gradient framework
    !
    function compute_search_direction_cg_single_parameter(parameter_name) result(srch)

        character(len=*), intent(in) :: parameter_name
        real, allocatable, dimension(:) :: srch

        real, allocatable, dimension(:) :: prev_grad, prev_srch, grad
        real :: s1, s2, beta

        ! Initialize search directions
        prev_grad = zeros(ng)
        prev_srch = zeros(ng)
        grad = zeros(ng)
        srch = zeros(ng)

        if (iter == 1) then
            ! Iteration 1 uses steepes descent

            ! Input current gradient
            call input_array(grad, dir_iter_model(iter)//'/grad_'//tidy(parameter_name)//'.bin')

            ! Calculate current search direction
            srch = -grad

        else
            ! For other iterations update search direction

            ! Input previous gradient
            call input_array(prev_grad, dir_iter_model(iter - 1)//'/grad_'//tidy(parameter_name)//'.bin')

            ! Input previous search direction
            call input_array(prev_srch, dir_iter_model(iter - 1)//'/srch_'//tidy(parameter_name)//'.bin')

            ! Input current gradient
            call input_array(grad, dir_iter_model(iter)//'/grad_'//tidy(parameter_name)//'.bin')

            s1 = 0.0
            s2 = 0.0

            ! Polak-Ribière formula
            s1 = sum(prev_grad*prev_grad)
            s2 = sum(grad*(grad - prev_grad))

            if (s1 == 0) then
                beta = 0.0
            else
                beta = max(0.0, s2/s1)
            end if

            ! Update search direction
            srch = -grad + beta*prev_srch

        end if

    end function compute_search_direction_cg_single_parameter

    !
    !> Caclulate search directions in Limited-memory BFGS framework
    !
    function compute_search_direction_lbfgs_single_parameter(parameter_name) result(srch)

        character(len=*), intent(in) :: parameter_name
        real, allocatable, dimension(:) :: srch

        character(len=1024) :: infile
        real, allocatable, dimension(:) :: grad, q, xi, xii, gi, gii, si, yi
        real, allocatable, dimension(:) :: alpha
        real :: beta
        real :: s1, s2
        integer :: i, l

        ! Allocate memory
        srch = zeros(ng)

        if (iter == 1) then
            ! For iter = 1, search direction the steepest-decent direction

            ! Allocate memory
            grad = zeros(ng)

            ! Input current gradient
            call input_array(grad, dir_iter_model(iter)//'/grad_'//tidy(parameter_name)//'.bin')

            ! Calculate current search direction
            srch = -grad

        else
            ! For iter >= 2, update search direction based on two-loop recursive l-BFGS formulas

            ! Allocate memory
            q = zeros(ng)
            xi = zeros(ng)
            xii = zeros(ng)
            gi = zeros(ng)
            gii = zeros(ng)
            si = zeros(ng)
            yi = zeros(ng)
            alpha = zeros(nstep_lbfgs)

            ! Current gradient
            infile = dir_iter_model(iter)//'/grad_'//tidy(parameter_name)//'.bin'
            call input_array(q, infile)

            ! Loop 1
            l = 1
            do i = iter - 1, max(iter - nstep_lbfgs, 1), -1

                ! model at iter = i+1 is the updated model at iter = i
                call input_array(xii, dir_iter_model(i)//'/updated_'//tidy(parameter_name)//'.bin')

                ! model at iter = i is the updated model at iter = i-1
                if (i - 1 == 0) then
                    call input_array(xi, dir_iter_model(i - 1)//'/'//tidy(parameter_name)//'.bin')
                else
                    call input_array(xi, dir_iter_model(i - 1)//'/updated_'//tidy(parameter_name)//'.bin')
                end if

                ! grad at iter = i+1
                call input_array(gii, dir_iter_model(i + 1)//'/grad_'//tidy(parameter_name)//'.bin')

                ! grad at iter = i
                call input_array(gi, dir_iter_model(i)//'/grad_'//tidy(parameter_name)//'.bin')

                s1 = 0.0
                s2 = 0.0

                yi = gii - gi
                si = xii - xi
                s1 = sum(si*q)
                s2 = sum(yi*si)
                if (s2 == 0.0) then
                    alpha(l) = 0.0
                else
                    alpha(l) = s1/s2
                end if

                q = q - alpha(l)*yi
                l = l + 1

            end do

            ! Initialize search direction

            ! model at iter = iter is the updated model at iter = iter-1
            call input_array(xii, dir_iter_model(iter - 1)//'/updated_'//tidy(parameter_name)//'.bin')

            ! model at iter = iter-1 is the updated model at iter = iter - 2
            if (iter - 2 == 0) then
                call input_array(xi, dir_iter_model(iter - 2)//'/'//tidy(parameter_name)//'.bin')
            else
                call input_array(xi, dir_iter_model(iter - 2)//'/updated_'//tidy(parameter_name)//'.bin')
            end if

            ! grad at iter = iter
            call input_array(gii, dir_iter_model(iter)//'/grad_'//tidy(parameter_name)//'.bin')

            ! grad at iter = iter-1
            call input_array(gi, dir_iter_model(iter - 1)//'/grad_'//tidy(parameter_name)//'.bin')

            s1 = 0.0
            s2 = 0.0

            yi = gii - gi
            si = xii - xi
            s1 = sum(si*yi)
            s2 = sum(yi*yi)
            if (s2 == 0.0) then
                srch = 0.0*q
            else
                srch = (s1/s2)*q
            end if

            ! Loop 2
            l = l - 1
            do i = max(iter - nstep_lbfgs, 1), iter - 1

                ! model at iter = i+1 is the updated model at iter = i
                call input_array(xii, dir_iter_model(i)//'/updated_'//tidy(parameter_name)//'.bin')

                ! model at iter = i is the updated model at iter = i - 1
                if (i - 1 == 0) then
                    call input_array(xi, dir_iter_model(i - 1)//'/'//tidy(parameter_name)//'.bin')
                else
                    call input_array(xi, dir_iter_model(i - 1)//'/updated_'//tidy(parameter_name)//'.bin')
                end if

                ! grad at iter = i+1
                call input_array(gii, dir_iter_model(i + 1)//'/grad_'//tidy(parameter_name)//'.bin')

                ! grad at iter = i
                call input_array(gi, dir_iter_model(i)//'/grad_'//tidy(parameter_name)//'.bin')

                s1 = 0.0
                s2 = 0.0

                yi = gii - gi
                si = xii - xi
                s1 = sum(yi*srch)
                s2 = sum(yi*si)
                if (s2 == 0.0) then
                    beta = 0.0
                else
                    beta = s1/s2
                end if

                srch = srch + si*(alpha(l) - beta)
                l = l - 1

            end do

            srch = -srch

        end if

    end function compute_search_direction_lbfgs_single_parameter

    !
    !> For L-BFGS framework, copy initial models to iteration 0 directories
    !
    subroutine init_lbfgs

        integer :: i

        if (to_lower(search_method) == 'l-bfgs' .and. rankid == 0) then
            call make_directory(dir_iter_model(0))
            do i = 1, nmodel
                call output_array(model_m(i)%array, &
                    dir_iter_model(0)//'/updated_'//tidy(model_name(i))//'.bin')
            end do
        end if

    end subroutine init_lbfgs

    !
    !> Calculate search directions
    !
    subroutine compute_search_direction

        real, allocatable, dimension(:) :: tmp
        integer :: i

        do i = 1, nmodel

            select case(model_name(i))

                case default

                    ng = nx*ny*nz
                    tmp = zeros(ng)

                    select case(model_search_method(i))
                        case ('SD', 'sd', 'steepest-descent')
                            tmp = -flatten(model_grad(i)%array)
                        case ('CG', 'cg', 'conjugate-gradient')
                            tmp = compute_search_direction_cg_single_parameter(model_name(i))
                        case ('L-BFGS', 'l-bfgs', 'l-BFGS')
                            tmp = compute_search_direction_lbfgs_single_parameter(model_name(i))
                    end select

#ifdef dim2
                    model_srch(i)%array = reshape(tmp, [nz, nx])
#endif
#ifdef dim3
                    model_srch(i)%array = reshape(tmp, [nz, ny, nx])
#endif

                case('sx', 'sy', 'sz', 'st0')

                    ng = nr_virtual
                    tmp = zeros(ng)

                    select case(model_search_method(i))
                        case ('SD', 'sd', 'steepest-descent')
                            tmp = -flatten(model_grad(i)%array)
                        case ('CG', 'cg', 'conjugate-gradient')
                            tmp = compute_search_direction_cg_single_parameter(model_name(i))
                        case ('L-BFGS', 'l-bfgs', 'l-BFGS')
                            tmp = compute_search_direction_lbfgs_single_parameter(model_name(i))
                    end select

#ifdef dim2
                    model_srch(i)%array = reshape(tmp, [nr_virtual, 1])
#endif
#ifdef dim3
                    model_srch(i)%array = reshape(tmp, [nr_virtual, 1, 1])
#endif

            end select

            if (rankid == 0) then
                call output_array(model_srch(i)%array, dir_iter_model(iter)//'/srch_'//tidy(model_name(i))//'.bin')
            end if

        end do

        call mpibarrier

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Search direction computation completed. ')
        end if

    end subroutine compute_search_direction

end module inversion_search_direction
