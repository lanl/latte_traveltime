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


submodule(gradient) reflector_imaging

    use traveltime_iso_reflection

contains

    !
    !> Compute reflector image shot by shot
    !
    module subroutine compute_image_shots

        real, allocatable, dimension(:, :) :: ttp_obs, tts_obs
        real, allocatable, dimension(:, :, :) :: vp, vs
        real, allocatable, dimension(:, :, :, :) :: img
        integer :: i, l
        character(len=1024) :: dir_field
        real, allocatable, dimension(:, :) :: weight

        ! temporary directory
        dir_scratch = tidy(dir_working)//'/scratch'
        dir_synthetic = dir_iter_synthetic(iter)
        dir_field = dir_iter_record(iter)
        if (rankid == 0) then
            call make_directory(dir_synthetic)
            call make_directory(dir_scratch)
        end if
        call mpibarrier

        ! Copy models
        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
            end select
        end do
        refl_all = zeros(nz, ny, nx, nrefl)

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

            weight = ones(gmtr(ishot)%nr, nrefl + 1)
            do i = 1, gmtr(ishot)%nr
                if (gmtr(ishot)%recr(i)%aoff < offset_min_refl .or. gmtr(ishot)%recr(i)%aoff > offset_max_refl) then
                    weight(i, :) = 0.0
                end if
            end do

            ! Read observed data
            ttp_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', gmtr(ishot)%nr, nrefl + 1)
            ttp_obs = ttp_obs*weight

            select case (which_medium)
                case ('acoustic-iso')
                    ! Here I only select reflection traveltimes for simplicity -- 2:
                    call image_iso( &
                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_obs(:, 2:), img)
                case ('elastic-iso')
                    ! Read observed data
                    ttp_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', gmtr(ishot)%nr, nrefl + 1)
                    tts_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin', gmtr(ishot)%nr, nrefl + 1)
                    ttp_obs = ttp_obs*weight
                    tts_obs = tts_obs*weight
                    ! Here I only select reflection traveltimes for simplicity -- 2:
                    call image_iso_elastic( &
                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp_obs(:, 2:), tts_obs(:, 2:), img)
            end select

            do i = 1, nmodel
                if (model_m(i)%name == 'refl') then
                    do l = 1, nrefl
                        call process_model_single_shot(ishot, refl_all(:, :, :, l), img(:, :, :, l), 'refl')
                    end do
                end if
            end do

            call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                //' reflector imaging is done. ')

        end do

        call mpibarrier

        do i = 1, nmodel
            if (model_m(i)%name == 'refl') then
                ! It is very challenging to do the reflector fitting simultaneously, therefore
                ! here I do reflectors separately
                do l = 1, nrefl
                    call allreduce_array(refl_all(:, :, :, l))
                    call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_taper_length_x', refltaper_x, 10*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_taper_length_y', refltaper_y, 10*dy, iter*1.0)
                    call assert(refltaper_x < nx*dx/3.0 .and. refltaper_y < ny*dy/3.0, ' <compute_image_shots> Error: tapers are too long. ')
                    call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_smooth_window_x', reflwindow_x, 0.2*(nx*dx - 2*refltaper_x), iter*1.0)
                    call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_smooth_window_y', reflwindow_y, 0.2*(ny*dy - 2*refltaper_y), iter*1.0)
                    call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_min', reflmin, oz, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_max', reflmax, oz + (nz - 1)*dz, iter*1.0)
                    call smooth_reflector(refl_all(:, :, :, l))
                    model_m(i)%array = model_m(i)%array + refl_all(:, :, :, l)*l
                end do
                call process_model_single_parameter(model_m(i)%array, 'refl')
            end if
        end do

    end subroutine compute_image_shots

end submodule reflector_imaging
