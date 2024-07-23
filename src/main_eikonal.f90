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
    use vars
    use traveltime_iso
    use traveltime_iso_reflection

    implicit none

#ifdef dim2
    real, allocatable, dimension(:, :) :: ttp, tts, vp, vs, refl
    real, allocatable, dimension(:, :, :) :: ttp_all, tts_all
#endif
#ifdef dim3
    real, allocatable, dimension(:, :, :) :: ttp, tts, vp, vs, refl
    real, allocatable, dimension(:, :, :, :) :: ttp_all, tts_all
#endif
    real, allocatable, dimension(:, :) :: ttprecr, ttsrecr
    integer :: i

    which_program = 'eikonal'

    call mpistart

    if (command_argument_count() == 0) then
        if (rankid == 0) then
            call warn('')
            call warn(date_time_compact()//' Error: Parameter file is missing. Exit. ')
            call warn('')
        end if
        call mpibarrier
        call mpistop
    end if

    if (rankid == 0) then
        call warn('')
        call warn(tile('=', 80))
        call warn(center_substring('Eikonal Traveltime Computation Starts', 80))
        call warn('')
        call print_date_time
        call warn('')
    end if

    ! Read parameters
    call read_parameters

    ! Set dimensions
    call set_regular_space

    ! Load geometry
    call load_geometry

    ! Divide shots to MPI ranks
    call divide_shots

    ! Read models
    call prepare_model

    ! Get model
    do i = 1, nmodel
        select case (model_m(i)%name)
            case ('vp')
                vp = model_m(i)%array
            case ('vs')
                vs = model_m(i)%array
            case ('refl')
                refl = model_m(i)%array
        end select
    end do

    ! Make directory
    call make_directory(dir_synthetic)

    call mpibarrier

    ! Show the medium parameter statistics
    if (rankid == 0) then
        do i = 1, nmodel
            call plot_histogram(model_m(i)%array, &
                label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
        end do
    end if
    call mpibarrier

    if (sum(snaps) > 0 .and. rankid == 0) then
        call make_directory(dir_snapshot)
    end if

    call mpibarrier

    ! Traveltime computation source by source
    do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

        ! Determine ranges of model selected for this source, if necessary
        call set_adaptive_model_range(gmtr(ishot))

        ! Solving eikional equation for traveltime
        select case (which_medium)

            case ('acoustic-iso')

#ifdef dim2
                if (allocated(refl)) then

                    call forward_iso_fast_sweep_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                        refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), ttp_all, ttprecr)

                else

                    select case (forward_eikonal_method)
                        case('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                        case('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                    end select

                end if

#endif
#ifdef dim3
                if (allocated(refl)) then

                    call forward_iso_fast_sweep_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), ttp_all, ttprecr)

                else

                    select case (forward_eikonal_method)
                        case('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                        case('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                    end select

                end if
#endif
                call output_array(ttprecr, tidy(dir_synthetic)// &
                    '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                if (sum(snaps) > 0) then
                    if (allocated(refl)) then
                        call output_array(ttp_all, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    else
                        call output_array(ttp, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    end if
                end if

            case ('elastic-iso')

#ifdef dim2
                if (allocated(refl)) then

                    select case (incident_wave)
                        case ('p')
                            call forward_iso_fast_sweep_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                ttp_all, tts_all, ttprecr, ttsrecr)
                        case ('s')
                            call forward_iso_fast_sweep_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                tts_all, ttp_all, ttsrecr, ttprecr)
                    end select

                else

                    select case (forward_eikonal_method)
                        case('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                            call forward_iso_fast_march( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, ttsrecr)
                        case('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                            call forward_iso_fast_sweep( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, ttsrecr)
                    end select

                end if

#endif
#ifdef dim3
                if (allocated(refl)) then

                    select case (incident_wave)
                        case ('p')
                            call forward_iso_fast_sweep_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                ttp_all, tts_all, ttprecr, ttsrecr)
                        case ('s')
                            call forward_iso_fast_sweep_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                tts_all, ttp_all, ttsrecr, ttprecr)
                    end select

                else

                    select case (forward_eikonal_method)
                        case('fast_march')
                            call forward_iso_fast_march( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                            call forward_iso_fast_march( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, ttsrecr)
                        case('fast_sweep')
                            call forward_iso_fast_sweep( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                            call forward_iso_fast_sweep( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, ttsrecr)
                    end select

                end if
#endif
                call output_array(ttprecr, tidy(dir_synthetic)// &
                    '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                call output_array(ttsrecr, tidy(dir_synthetic)// &
                    '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                if (sum(snaps) > 0) then
                    if (allocated(refl)) then
                        call output_array(ttp_all, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                        call output_array(tts_all, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                    else
                        call output_array(ttp, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                        call output_array(tts, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                    end if
                end if

        end select

        call warn(date_time_compact()//' Source '//num2str(gmtr(ishot)%id)//' traveltime computation is done. ')

    end do

    call mpibarrier

    if (rankid == 0) then
        call warn('')
        call print_date_time
        call warn('')
        call warn(center_substring('EIKONAL ends', 80))
        call warn(tile('=', 80))
        call warn('')
    end if

    call mpiend

end program main
