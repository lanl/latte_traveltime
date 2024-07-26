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


module vars

    use libflit
    use parameters

    implicit none

    type(meta_array2_real), allocatable, dimension(:) :: model_aux
    type(meta_array2_real), allocatable, dimension(:) :: model_m, model_m_backup
    type(meta_array2_real), allocatable, dimension(:) :: model_grad, model_srch, model_reg

contains

    !
    !> Initialize a model
    !
    subroutine prepare_model_single_parameter(w, name, file_w, const, source, update)

        real, allocatable, dimension(:, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: file_w
        real, intent(in), optional :: const
        real, dimension(:, :), intent(in), optional :: source
        logical, intent(in), optional :: update

        logical :: update_this_model
        character(len=1024) :: file_update
        real, allocatable, dimension(:, :) :: wt
        character(len=32) :: interpolation_method

        select case(name)

            case default
                ! For model parameter update

                if (name == 'refl') then
                    interpolation_method = 'nearest'
                else
                    interpolation_method = 'linear'
                end if

                w = zeros(nz, nx)

                ! if read in or assign const value
                if (present(file_w) .and. file_w /= '') then

                    if (require_model_interp) then
                        ! if it is necessary to resample

                        call alloc_array(wt, [1, nz0, 1, nx0])
                        call input_array(wt, file_w)
                        w = interp(wt, [nz0, nx0], [dz0, dx0], [oz0, ox0], &
                            [nz, nx], [dz, dx], [oz, ox], [interpolation_method, interpolation_method])
                        deallocate (wt)

                    else
                        ! if not resampling required
                        call input_array(w, file_w)

                    end if

                else

                    if (present(const)) then
                        w = const
                    end if
                    if (file_w == '') then
                        if (rankid == 0) then
                            call warn(' Warning: Model '//tidy(name)//' is empty and is set to '//num2str(w(1, 1)))
                        end if
                    end if

                end if

            case ('sx', 'sy', 'sz', 'st0')
                ! For source parameter update

                w = zeros(nr_virtual, 1)

                ! if read in or assign const value
                if (present(file_w) .and. file_w /= '') then

                    call input_array(w, file_w)

                else

                    if (present(const)) then
                        w = const
                    end if
                    if (file_w == '') then
                        if (rankid == 0) then
                            call warn(' Warning: Source parameter '//tidy(name)//' is empty and is set to '//num2str(w(1, 1)))
                        end if
                    end if

                end if

        end select

        ! If the source is given
        if (present(source)) then
            w = source
        end if

        ! If the inversion starts from certain iteration other than 1
        update_this_model = .true.
        if (present(update)) then
            if (.not.update) then
                update_this_model = .false.
            end if
        end if

        if (resume_from_iter > 1 .and. update_this_model) then
            file_update = tidy(dir_working)//'/iteration_'// &
                num2str(resume_from_iter - 1)//'/model/updated_'//tidy(name)//'.bin'
            if (file_exists(file_update)) then
                call input_array(w, file_update)
            else
                if (rankid == 0) then
                    call warn(' Error: Updated model '//tidy(name)//' does not exist. Exit')
                end if
                stop
            end if
        end if

    end subroutine prepare_model_single_parameter

    !
    !> Input medium parameter models for inversion
    !
    subroutine prepare_model

        integer :: i
        character(len=1024) :: fname

        select case (which_program)

            case('fatt', 'trtt', 'tloc')

                ! Model static
                if (nmodel_aux >= 1) then
                    allocate (model_aux(1:nmodel_aux))
                    do i = 1, nmodel_aux
                        ! Models that will not be updated
                        call readpar_string(file_parameter, 'file_'//tidy(model_name_aux(i)), fname, '', .true.)
                        call prepare_model_single_parameter(model_aux(i)%array, model_name_aux(i), fname, update=.false.)
                        model_aux(i)%name = model_name_aux(i)
                    end do
                end if

                ! Model to update
                allocate (model_m(1:nmodel))
                allocate (model_grad(1:nmodel))
                allocate (model_srch(1:nmodel))
                allocate (model_m_backup(1:nmodel))

                do i = 1, nmodel
                    ! Models to be updated
                    call readpar_string(file_parameter, 'file_'//tidy(model_name(i)), fname, '')
                    call prepare_model_single_parameter(model_m(i)%array, model_name(i), fname, update=.true.)

                    model_m(i)%name = model_name(i)
                    model_grad(i)%name = model_name(i)
                    model_grad(i)%array = zeros_like(model_m(i)%array)
                    model_srch(i)%name = model_name(i)
                    model_srch(i)%array = zeros_like(model_m(i)%array)
                    model_m_backup(i)%name = model_name(i)
                    model_m_backup(i)%array = zeros_like(model_m(i)%array)

                end do

            case ('eikonal')

                allocate(model_m(1:nmodel))
                do i = 1, nmodel
                    call readpar_string(file_parameter, 'file_'//tidy(model_name(i)), fname, '')
                    call prepare_model_single_parameter(model_m(i)%array, model_name(i), fname)
                    model_m(i)%name = model_name(i)
                end do

        end select

    end subroutine prepare_model

    !
    !> Clip Vp or Vs to ensure that Vp/Vs ratio is within
    !> a reasonable range ([1.2, 3] by default).
    !
    subroutine clip_vpvsratio

        integer :: i
        real, allocatable, dimension(:, :) :: mvp, mvs, r
        logical :: vp_from_aux, vs_from_aux

        if (which_medium == 'elastic-iso' .or. which_medium == 'elastic-tti') then

            vp_from_aux = .false.
            vs_from_aux = .false.

            ! Vp or Vs might come from static models
            do i = 1, nmodel_aux
                if (model_aux(i)%name == 'vp') then
                    mvp = model_aux(i)%array
                    vp_from_aux = .true.
                end if
                if (model_aux(i)%name == 'vs') then
                    mvs = model_aux(i)%array
                    vs_from_aux = .true.
                end if
            end do

            ! Vp or/and Vs might come from update models
            do i = 1, nmodel
                if (model_m(i)%name == 'vp') then
                    mvp = model_m(i)%array
                    vp_from_aux = .false.
                end if
                if (model_m(i)%name == 'vs') then
                    mvs = model_m(i)%array
                    vs_from_aux = .false.
                end if
            end do

            ! If Vp is static model but Vs is update model, or both Vp and Vs are update model
            ! Then use Vp/Vs ratio to constrain Vs
            if ((vp_from_aux .and. .not. vs_from_aux) .or. (.not. vp_from_aux .and. .not. vs_from_aux)) then

                do i = 1, nmodel
                    if (model_m(i)%name == 'vs') then
                        where (mvs /= 0)
                            mvs = clip(mvs, model_min(i), model_max(i))
                        end where
                    end if
                end do

                r = clip(mvp/mvs, min_vpvsratio, max_vpvsratio)
                where (mvs == 0)
                    r = 0
                end where

                if (vpvsratio_smoothx /= 0 .or. vpvsratio_smoothz /= 0) then
                    r = gauss_filt(r, [vpvsratio_smoothz/dz, vpvsratio_smoothx/dx])
                end if
                r = clip(r, min_vpvsratio, max_vpvsratio)

                do i = 1, nmodel
                    if (model_m(i)%name == 'vs') then
                        model_m(i)%array = mvp/r
                        where (mvs == 0)
                            model_m(i)%array = 0
                        end where
                    end if
                end do

            end if

            ! If Vs is static but Vp is update,
            ! Then use Vp/Vs ratio to contrain Vp
            if (vs_from_aux .and. .not. vp_from_aux) then

                do i = 1, nmodel
                    if (model_m(i)%name == 'vp') then
                        where (mvp /= 0)
                            mvp = clip(mvp, model_min(i), model_max(i))
                        end where
                    end if
                end do

                r = clip(mvp/mvs, min_vpvsratio, max_vpvsratio)
                where (mvp == 0)
                    r = 0
                end where

                if (vpvsratio_smoothx /= 0 .or. vpvsratio_smoothz /= 0) then
                    r = gauss_filt(r, [vpvsratio_smoothz/dz, vpvsratio_smoothx/dx])
                end if
                r = clip(r, min_vpvsratio, max_vpvsratio)

                do i = 1, nmodel
                    if (model_m(i)%name == 'vp') then
                        model_m(i)%array = mvs*r
                        where (mvp == 0)
                            model_m(i)%array = 0
                        end where
                    end if
                end do

            end if

        end if

    end subroutine clip_vpvsratio

    !
    !> Zero regularization auxiliary variables
    !
    subroutine init_reg

        integer :: i

        if (yn_regularize_model .or. yn_regularize_source) then

            allocate (model_reg(1:nmodel))

            do i = 1, nmodel
                model_reg(i)%name = model_name(i)
                model_reg(i)%array = model_m(i)%array
            end do

        end if

    end subroutine init_reg

end module vars

