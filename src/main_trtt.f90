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


program main

    use libflit
    use parameters
    use utility
    use inversion_gradient
    use inversion_search_direction
    use inversion_step_size

    implicit none

#ifdef dim2
    type(meta_array2_real) :: m
#endif
#ifdef dim3
    type(meta_array3_real) :: m
#endif

    which_program = 'trtt'

    ! Initialization
    call mpistart

    if (command_argument_count() == 0) then
        if (rankid == 0) then
            call warn('')
            call warn(date_time_compact()//' Error: Parameter file does not exist. Exit. ')
            call warn('')
        end if
        call mpibarrier
        call mpistop
    end if

    iter = 0

    ! Read parameters
    call read_parameters

    if (rankid == 0) then
        call warn('')
        call warn(tile('=', 80))
        call warn(center_substring('Joint Transmission-Reflection Tomography Starts', 80))
        call warn('')
        call print_date_time
        call warn('')
    end if

    ! Set space dimensions
    call set_regular_space

    ! Load geometry
    call load_geometry

    ! Divide shots
    call divide_shots

    ! Load models
    call prepare_model

    ! Initialize regularization models if necessary
    call init_reg

    ! For l-BFGS, prepare iter = 0 models if necessary
    call init_lbfgs

    ! Initialize misfit arrays
    call print_misfit

    call mpibarrier

    ! Update
    iter = resume_from_iter
    do while (iter <= niter_max)

        ! Print info
        if (rankid == 0) then
            call warn('')
            call warn('================================================')
            call warn(' Iteration '//num2str(iter))
            call warn('')
        end if

        ! Make update directories
        call make_iter_dir

        ! Update reflector image
        if (iter == 1 .and. any_in(['refl'], model_name)) then
            m = get_meta_array(model_m, 'refl')
            if (maxval(abs(m%array)) == 0) then
                call compute_reflector_image
            end if
        end if

        ! Compute gradients
        call compute_gradient

        ! Compute search directions
        call compute_search_direction

        ! Compute inversion step size
        call compute_step_size

        ! Update ADMM variables if necessary
        call compute_regularization

        ! Print misfit information
        call print_misfit

        ! If not satisfied then update
        iter = iter + 1

    end do

    call mpibarrier

    ! End
    if (rankid == 0) then
        call warn('')
        call print_date_time
        call warn('')
        call warn(center_substring('TRTT ends', 80))
        call warn(tile('=', 80))
        call warn('')
    end if

    call mpiend

end program main
