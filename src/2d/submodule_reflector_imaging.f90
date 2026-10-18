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


submodule(gradient) reflector_imaging

    use traveltime_iso_reflection

contains

    !
    !> Image the reflectors shot by shot in the current velocity models, and
    !> build the reflector model from the per-shot images
    !
    module subroutine compute_image_shots

        real, allocatable, dimension(:, :) :: ttp_obs, tts_obs
        real, allocatable, dimension(:, :) :: vp, vs, refl_prev
        real, allocatable, dimension(:, :) :: zidx
        real, allocatable, dimension(:, :) :: weight
        ! Per-shot reflector images: the fractional z index of each reflector
        ! in each column, by shot, or 0 where the shot has no crossing
        real, allocatable, dimension(:, :, :) :: refl_depth
        integer :: i, j, l, irefl
        real :: omin, omax
        character(len=1024) :: dir_field

        dir_field = dir_iter_record(iter)

        ! Velocity models, whether they are updated or fixed
        irefl = 0
        do i = 1, nmodel
            select case (model_m(i)%name)
                case ('vp')
                    vp = model_m(i)%array
                case ('vs')
                    vs = model_m(i)%array
                case ('refl')
                    irefl = i
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
        call assert(irefl > 0, ' <compute_image_shots> Error: refl is not in model_update. ')
        call assert(allocated(vp), ' <compute_image_shots> Error: vp is neither in model_update nor in model_aux. ')
        if (which_medium == 'elastic-iso') then
            call assert(allocated(vs), ' <compute_image_shots> Error: vs is neither in model_update nor in model_aux. ')
        end if

        refl_depth = zeros(nx, ns, nrefl)
        call reflection_offset_range(omin, omax)

        do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

            call set_adaptive_model_range(gmtr(ishot))

            ! Reflection data of the receivers in use, in the offset range of
            ! the reflector depths (omin, omax); the others are zeroed
            weight = ones(gmtr(ishot)%nr, nrefl + 1)
            do i = 1, gmtr(ishot)%nr
                if (gmtr(ishot)%recr(i)%weight == 0 .or. gmtr(ishot)%recr(i)%aoff < omin &
                        .or. gmtr(ishot)%recr(i)%aoff > omax) then
                    weight(i, :) = 0.0
                end if
            end do

            select case (which_medium)
                case ('acoustic-iso')
                    ttp_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', gmtr(ishot)%nr, nrefl + 1)
                    ttp_obs = ttp_obs*weight
                    call image_iso( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_obs(:, 2:), zidx)
                case ('elastic-iso')
                    ttp_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin', gmtr(ishot)%nr, nrefl + 1)
                    tts_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin', gmtr(ishot)%nr, nrefl + 1)
                    ttp_obs = ttp_obs*weight
                    tts_obs = tts_obs*weight
                    ! The incident leg is P for incident_wave = p and S for incident_wave = s;
                    ! the reflected leg of the p data is P, that of the s data S
                    select case (incident_wave)
                        case ('p')
                            call image_iso_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp_obs(:, 2:), tts_obs(:, 2:), zidx)
                        case ('s')
                            call image_iso_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts_obs(:, 2:), ttp_obs(:, 2:), zidx)
                    end select
            end select

            ! From the shot's range to the model's columns and rows
            do l = 1, nrefl
                do j = 1, shot_nx
                    if (zidx(j, l) > 0) then
                        refl_depth(shot_nxbeg + j - 1, ishot, l) = zidx(j, l) + shot_nzbeg - 1
                    end if
                end do
            end do

            call warn(date_time_compact()//' >> Shot '//num2str(gmtr(ishot)%id) &
                //' reflector imaging completed. ')

        end do

        call mpibarrier
        call allreduce_array(refl_depth)

        ! Build the reflectors one by one from the per-shot crossings. A
        ! reflector that no shot images keeps its previous depths.
        refl_prev = model_m(irefl)%array
        model_m(irefl)%array = 0.0
        do l = 1, nrefl
            call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_taper_length', refltaper, 10*dx, iter*1.0)
            call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_smooth_window', reflwindow, 0.0, iter*1.0)
            call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_min', reflmin, oz, iter*1.0)
            call readpar_xfloat(file_parameter, 'reflector_'//num2str(l)//'_depth_max', reflmax, oz + (nz - 1)*dz, iter*1.0)
            call build_reflector(refl_depth(:, :, l), l, model_m(irefl)%array)
            if (.not. any(nint(model_m(irefl)%array) == l)) then
                call assert(any(nint(refl_prev) == l), ' <compute_image_shots> Error: no shot images reflector ' &
                    //num2str(l)//', and the reflector model does not have it. ')
                where (nint(refl_prev) == l)
                    model_m(irefl)%array = refl_prev
                end where
                if (rankid == 0) then
                    call warn(date_time_compact()//' Reflector '//num2str(l)//' keeps its previous depths. ')
                end if
            end if
        end do

    end subroutine compute_image_shots

end submodule reflector_imaging
