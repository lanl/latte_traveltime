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


module parameters

    use libflit

    implicit none

    character(len=1024) :: file_parameter

    character(len=32) :: which_program

    character(len=1024) :: dir_record
    character(len=1024) :: dir_working
    character(len=1024) :: dir_scratch
    character(len=1024) :: dir_synthetic
    character(len=1024) :: dir_snapshot
    character(len=1024) :: dir_adjoint

    character(len=1204) :: file_geometry

    integer :: ns = 1

    ! meditum type
    character(len=24) :: which_medium

    ! snapshot times
    real, allocatable, dimension(:) :: snaps

    ! space taper when adaptive range required
    integer :: space_taperx = 10
    integer :: space_tapery = 10
    integer :: space_taperz = 10

    ! shot index and division related
    integer, allocatable, dimension(:, :) :: shot_in_rank
    character(len=1024) :: shot_prefix
    integer :: ishot
    integer, allocatable, dimension(:) :: shot_index, rec_index
    integer, allocatable, dimension(:) :: src_exclude, sid_exclude, rec_exclude

    real :: sxmin = -float_huge
    real :: sxmax = +float_huge
    real :: symin = -float_huge
    real :: symax = +float_huge
    real :: szmin = -float_huge
    real :: szmax = +float_huge
    real :: rxmin = -float_huge
    real :: rxmax = +float_huge
    real :: rymin = -float_huge
    real :: rymax = +float_huge
    real :: rzmin = -float_huge
    real :: rzmax = +float_huge

    ! only compute misfit
    logical :: yn_misfit_only
    real :: misfit0, misfit
    real :: val_misfit

    ! search direction computation method in inversion
    character(len=32), allocatable, dimension(:) :: model_search_method
    character(len=32) :: search_method

    ! regularization for inversion
    character(len=32), allocatable, dimension(:) :: model_regularization_method
    character(len=32), allocatable, dimension(:) :: source_regularization_method
    logical :: yn_regularize_model = .false.
    logical :: yn_regularize_source = .false.

    ! inversion parameters
    integer :: iter
    integer :: niter_max
    character(len=32) :: misfit_type

    ! vp/vs ratio smoothing
    real :: vpvsratio_smoothx = 0
    real :: vpvsratio_smoothy = 0
    real :: vpvsratio_smoothz = 0

    ! method of computing step size: linear or line search
    character(len=32) :: step_size_method

    ! step
    real :: step_scaling_factor
    real :: step_max_scale_factor = 1.0

    ! data misfit arrays
    real, allocatable, dimension(:) :: step_misfit
    real, allocatable, dimension(:) :: data_misfit
    real, allocatable, dimension(:, :) :: shot_misfit
    real, allocatable, dimension(:) :: misfit_weight

    ! maximum and minimum offset
    real :: offset_min = 0.0d0
    real :: offset_max = float_huge
    real :: offset_min_refl = 0.0d0
    real :: offset_max_refl = float_huge

    ! source id
    integer :: sid_min = 1
    integer :: sid_max = int4_huge

    ! maximum and minimum shot number
    integer :: shot_min = 1
    integer :: shot_max = 1
    integer :: shot_every = 1

    ! maximum and minimum receiver number
    integer :: rec_min = 1
    integer :: rec_max = int4_huge
    integer :: rec_every = 1

    ! necessity of adaptivity
    logical :: yn_adpx
    logical :: yn_adpy
    logical :: yn_adpz

    real :: adpextrax = 0.0
    real :: adpextray = 0.0
    real :: adpextraz = 0.0

    integer, allocatable, dimension(:) :: sid_select, src_select

    character(len=24) :: aniso_param = 'thomsen'

    real :: min_vpvsratio, max_vpvsratio

    logical :: verbose
    integer :: sweep_niter_max = 10
    real :: sweep_stop_threshold = 1.0e-4
    character(len=32) :: forward_eikonal_method
    real :: misfit_threshold
    real :: energybal_power = 1.0

    character(len=1024) :: dir_gradmask
    character(len=1024) :: file_andfaux, file_andfcoh

    integer :: resume_from_iter

    character(len=1024) :: file_datamisfit, file_shotmisfit
    logical :: yn_continue_inv = .false.
    logical :: yn_enforce_update = .false.
    logical :: put_synthetic_in_scratch = .false.

    character(len=32), allocatable, dimension(:) :: record_processing
    character(len=32), allocatable, dimension(:) :: encoded_record_processing
    character(len=32), allocatable, dimension(:) :: synthetic_processing
    character(len=32), allocatable, dimension(:) :: gradient_processing
    character(len=32), allocatable, dimension(:) :: adjoint_source_processing

    character(len=32), allocatable, dimension(:) :: process_shot_grad, process_grad

    character(len=24), allocatable, dimension(:) :: model_name, model_name_aux
    integer :: nmodel, nmodel_aux

    real, allocatable, dimension(:) :: model_step, model_step_max
    real, allocatable, dimension(:) :: model_min, model_max

    character(len=24), allocatable, dimension(:) :: data_name
    integer :: ndata

    logical :: require_model_interp = .false.
    logical :: uniform_processing = .true.

    logical :: yn_flat_stop = .false.

    character(len=32) :: data_prefix = 'traveltime'

    integer :: nrefl = 0
    real :: reflector_imaging_threshold = 1.0e-9
    logical :: yn_update_reflector
    logical :: yn_precond = .true.

    character(len=12) :: incident_wave = 'p'
    real :: precond_eps = 1.0e-4

    ! Dimension
    integer :: nx0, ny0, nz0
    real :: dx0, dy0, dz0
    real :: ox0, oy0, oz0

    integer :: nx, ny, nz
    real :: dx, dy, dz
    real :: ox, oy, oz

    real :: xmin = -float_large
    real :: ymin = -float_large
    real :: zmin = -float_large
    real :: xmax = +float_large
    real :: ymax = +float_large
    real :: zmax = +float_large

    integer :: nxbeg1, nxbeg2, nxend1, nxend2, nxbeg_e, nxend_e
    integer :: nybeg1, nybeg2, nyend1, nyend2, nybeg_e, nyend_e
    integer :: nzbeg1, nzbeg2, nzend1, nzend2, nzbeg_e, nzend_e

    real :: shot_xbeg, shot_ybeg, shot_zbeg
    real :: shot_xend, shot_yend, shot_zend
    integer :: shot_nxbeg, shot_nybeg, shot_nzbeg
    integer :: shot_nxend, shot_nyend, shot_nzend
    integer :: shot_nx, shot_ny, shot_nz

    integer :: regular_nx, regular_ny, regular_nz
    real :: regular_dx, regular_dy, regular_dz
    real :: regular_ox, regular_oy, regular_oz

    type source
        real :: x, y, z
        real :: t0 = 0.0
        real :: amp = 1.0
    end type source

    type receiver
        real :: x, y, z
        real :: aoff
        real :: weight = 1.0
        real :: t0 = 0.0
    end type receiver

    type source_receiver_geometry
        integer :: id
        integer :: ns
        type(source), allocatable, dimension(:) :: srcr
        integer :: nr
        type(receiver), allocatable, dimension(:) :: recr
    end type source_receiver_geometry

    type(source_receiver_geometry), allocatable, dimension(:) :: gmtr

    ! For TLOC
    type(source_receiver_geometry), allocatable, dimension(:) :: gmtr_real
    logical :: yn_exchange_sr
    integer :: nr_virtual
    logical :: yn_dd_no_st0

contains

    !
    ! Read global parameters for forward modeling or inversion
    !
    subroutine read_parameters

        character(len=1024) :: tmpparam
        integer :: i
        real :: valmin, valmax

        ! Get name of the parameter file
        call get_command_argument(1, file_parameter)

        ! Check existence of the parameter file
        if (.not. file_exists(file_parameter)) then
            call warn(date_time_compact()//' <read_parameters> Error: Parameter file = ' &
                //tidy(file_parameter)//' does not exist. Exit. ')
            call mpibarrier
            call mpistop
        end if

        ! Copy parameter file to working directory and use it as the working parameter file
        select case (which_program)
            case ('eikonal')
                call readpar_string(file_parameter, 'dir_synthetic', dir_synthetic, './data_synthetic')
                call make_directory(dir_synthetic)
                if (rankid == 0) then
                    tmpparam = tidy(dir_synthetic)//'/parameters.'//tidy(which_program)//'.'//date_time_string()
                    call copy_file(file_parameter, tmpparam)
                    file_parameter = tmpparam
                    call warn(' Parameter file: '//tidy(file_parameter))
                end if
            case ('fatt', 'trtt', 'tloc')
                call readpar_string(file_parameter, 'dir_working', dir_working, './test')
                call make_directory(dir_working)
                if (rankid == 0) then
                    tmpparam = tidy(dir_working)//'/parameters.'//tidy(which_program)//'.'//date_time_string()
                    call copy_file(file_parameter, tmpparam)
                    file_parameter = tmpparam
                    call warn(' Parameter file: '//tidy(file_parameter))
                end if
        end select
        call mpibarrier

        call bcast_array(file_parameter)

        ! which kind of medium
        call readpar_string(file_parameter, 'which_medium', which_medium, 'acoustic-iso')
        call readpar_float(file_parameter, 'min_vpvsratio', min_vpvsratio, 1.1)
        call readpar_float(file_parameter, 'max_vpvsratio', max_vpvsratio, 9.0)
        call readpar_float(file_parameter, 'vpvsratio_smoothx', vpvsratio_smoothx, 0.0)
        call readpar_float(file_parameter, 'vpvsratio_smoothy', vpvsratio_smoothy, 0.0)
        call readpar_float(file_parameter, 'vpvsratio_smoothz', vpvsratio_smoothz, 0.0)
        call readpar_string(file_parameter, 'incident_wave', incident_wave, 'p')

        call readpar_string(file_parameter, 'data_prefix', data_prefix, 'traveltime')
        call readpar_float(file_parameter, 'energybal_power', energybal_power, 1.0)
        call readpar_float(file_parameter, 'misfit_threshold', misfit_threshold, float_huge)

        call readpar_int(file_parameter, 'ns', ns, 1)
        call readpar_string(file_parameter, 'file_geometry', file_geometry, '', required=.true.)
        call readpar_nint(file_parameter, 'src_index', shot_index, [1, 1, ns])
        shot_min = max(shot_index(1), 1)
        shot_every = max(shot_index(2), 1)
        shot_max = min(shot_index(3), ns)
        call readpar_nint(file_parameter, 'rec_index', rec_index, [1, 1, int4_huge])
        rec_min = max(rec_index(1), 1)
        rec_every = max(rec_index(2), 1)
        rec_max = min(rec_index(3), int4_huge)
        call readpar_nint(file_parameter, 'rec_exclude', rec_exclude, [0])
        call readpar_float(file_parameter, 'offset_min', offset_min, 0.0)
        call readpar_float(file_parameter, 'offset_max', offset_max, float_huge)
        call readpar_nint(file_parameter, 'sid_select', sid_select, [0])
        call readpar_nint(file_parameter, 'src_select', src_select, [0])
        call readpar_nint(file_parameter, 'sid_exclude', sid_exclude, [0])
        call readpar_nint(file_parameter, 'src_exclude', src_exclude, [0])
        call readpar_int(file_parameter, 'sid_min', sid_min, 0)
        call readpar_int(file_parameter, 'sid_max', sid_max, int4_huge)
        call readpar_float(file_parameter, 'sxmin', sxmin, -float_huge)
        call readpar_float(file_parameter, 'sxmax', sxmax, +float_huge)
        call readpar_float(file_parameter, 'symin', symin, -float_huge)
        call readpar_float(file_parameter, 'symax', symax, +float_huge)
        call readpar_float(file_parameter, 'szmin', szmin, -float_huge)
        call readpar_float(file_parameter, 'szmax', szmax, +float_huge)
        call readpar_float(file_parameter, 'rxmin', rxmin, -float_huge)
        call readpar_float(file_parameter, 'rxmax', rxmax, +float_huge)
        call readpar_float(file_parameter, 'rymin', rymin, -float_huge)
        call readpar_float(file_parameter, 'rymax', rymax, +float_huge)
        call readpar_float(file_parameter, 'rzmin', rzmin, -float_huge)
        call readpar_float(file_parameter, 'rzmax', rzmax, +float_huge)
        call readpar_logical(file_parameter, 'yn_exchange_sr', yn_exchange_sr, ifthen(which_program == 'tloc', .true., .false.))

        call readpar_string(file_parameter, 'dir_working', dir_working, './test')
        call readpar_string(file_parameter, 'dir_record', dir_record, './data')

        call readpar_nstring(file_parameter, 'model_regularization_method', model_regularization_method, [''])
        call readpar_nstring(file_parameter, 'source_regularization_method', source_regularization_method, [''])
        if (model_regularization_method(1) /= '') then
            yn_regularize_model = .true.
        end if
        if (source_regularization_method(1) /= '') then
            yn_regularize_source = .true.
        end if
        call readpar_int(file_parameter, 'niter_max', niter_max, 100)

        call readpar_string(file_parameter, 'search_method', search_method, 'CG')

        call readpar_string(file_parameter, 'misfit_type', misfit_type, 'ad')
        call assert(any_in([misfit_type], ['absolute-difference', 'ad', 'double-difference', 'dd']), &
            ' <read_parameter> Error: misfit_type must be one of absolute-difference, ad, double-difference, dd')

        call readpar_logical(file_parameter, 'yn_adpx', yn_adpx, .false.)
        call readpar_logical(file_parameter, 'yn_adpy', yn_adpy, .false.)
        call readpar_logical(file_parameter, 'yn_adpz', yn_adpz, .false.)
        call readpar_float(file_parameter, 'adp_extrax', adpextrax, 0.0)
        call readpar_float(file_parameter, 'adp_extray', adpextray, 0.0)
        call readpar_float(file_parameter, 'adp_extraz', adpextraz, 0.0)
        call readpar_int(file_parameter, 'adp_taperx', space_taperx, 0)
        call readpar_int(file_parameter, 'adp_tapery', space_tapery, 0)
        call readpar_int(file_parameter, 'adp_taperz', space_taperz, 0)

        call readpar_string(file_parameter, 'step_size_method', step_size_method, 'linear')
        call readpar_logical(file_parameter, 'yn_precond', yn_precond, .true.)
        call readpar_float(file_parameter, 'precond_eps', precond_eps, 1.0e-4)

        call readpar_nfloat(file_parameter, 'snaps', snaps, [-1.0])
        call readpar_string(file_parameter, 'dir_snapshot', dir_snapshot, './snapshot')
        call readpar_string(file_parameter, 'dir_synthetic', dir_synthetic, './data_synthetic')
        call assert(dir_snapshot /= dir_synthetic, ' <read_parameters> Error: dir_snapshot and dir_synthetic cannot be the same. ')

        call readpar_logical(file_parameter, 'verbose', verbose, .false.)
        call readpar_float(file_parameter, 'sweep_stop_threshold', sweep_stop_threshold, 1.0e-4)
        call readpar_int(file_parameter, 'sweep_niter_max', sweep_niter_max, 10)
        call readpar_string(file_parameter, 'forward_eikonal_method', forward_eikonal_method, 'fast_sweep')

        call readpar_string(file_parameter, 'file_data_misfit', file_datamisfit, tidy(dir_working)//'/data_misfit.txt')
        call readpar_string(file_parameter, 'file_shot_misfit', file_shotmisfit, tidy(dir_working)//'/shot_misfit.bin')

        call readpar_logical(file_parameter, 'yn_continue', yn_continue_inv, .false.)
        if (yn_continue_inv) then
            resume_from_iter = max(1, count_nonempty_lines(file_datamisfit))
        else
            call readpar_int(file_parameter, 'resume_from_iter', resume_from_iter, 1)
        end if

        call readpar_nstring(file_parameter, 'process_record', record_processing, [''])
        call readpar_nstring(file_parameter, 'process_synthetic', synthetic_processing, [''])
        call readpar_logical(file_parameter, 'yn_flat_stop', yn_flat_stop, .false.)

        if (which_program == 'eikonal') then

            call readpar_nstring(file_parameter, 'model_name', model_name, [''], required=.true.)
            nmodel = size(model_name)

        else

            call readpar_nstring(file_parameter, 'model_update', model_name, ['vp'])
            call assert(model_name(1) /= '', ' <read_parameter> Error: model_update cannot be empty or start with null')
            nmodel = size(model_name)
            call readpar_nstring(file_parameter, 'model_aux', model_name_aux, [''])
            if (model_name_aux(1) == '') then
                nmodel_aux = 0
            else
                nmodel_aux = size(model_name_aux)
            end if
            model_min = zeros(nmodel)
            model_max = zeros(nmodel)
            model_step_max = zeros(nmodel)
            model_step = zeros(nmodel)
            allocate(model_search_method(1:nmodel))
            do i = 1, nmodel
                select case (remove_string_after(model_name(i), ['c', 'C']))
                    case ('vp', 'vs')
                        valmin = 0.0
                        valmax = 1.0e5
                    case ('epsilon', 'delta', 'gamma', 'eps', 'del', 'gam')
                        valmin = 0.0
                        valmax = 0.5
                    case ('c', 'C')
                        valmin = 0.0
                        valmax = 1.0e9
                    case ('refl')
                        valmin = 0.0
                        valmax = 1.0e9
                end select
                call readpar_float(file_parameter, 'min_'//tidy(model_name(i)), model_min(i), valmin)
                call readpar_float(file_parameter, 'max_'//tidy(model_name(i)), model_max(i), valmax)
                ! The program allows using different search method for different parameters
                call readpar_string(file_parameter, 'search_method_'//tidy(model_name(i)), model_search_method(i), search_method)
                ! For reflector, the only available option now is SD
                if (model_name(i) == 'refl') then
                    model_search_method(i) = 'SD'
                end if
            end do

            if (which_program == 'tloc' .and. .not.any_in(['st0'], model_name) .and. (misfit_type == 'dd' .or. misfit_type == 'double-difference')) then
                yn_dd_no_st0 = .true.
            else
                yn_dd_no_st0 = .false.
            end if

        end if

        select case (which_medium)
            case ('acoustic-iso', 'acoustic-tti')
                call readpar_nstring(file_parameter, 'data_name', data_name, ['p'])
                call assert(size(data_name) == 1 .and. data_name(1) == 'p', &
                    ' <read_parameter> Error: For acoustic, data_name must = p' )
            case('elastic-iso', 'elastic-tti')
                call readpar_nstring(file_parameter, 'data_name', data_name, ['p', 's'])
                call assert(size(data_name) == 2 .and. data_name(1) == 'p' .and. data_name(2) == 's', &
                    ' <read_parameter> Error: For acoustic, data_name must = p, s' )
        end select
        ndata = size(data_name)

        if (which_program == 'trtt') then
            call readpar_float(file_parameter, 'reflector_imaging_threshold', reflector_imaging_threshold, 1.0e-9)
            call readpar_int(file_parameter, 'nrefl', nrefl, 1)
            call readpar_nfloat(file_parameter, 'misfit_weight', misfit_weight, ones(nrefl + 1))
            call assert(size(misfit_weight) == nrefl + 1, ' <read_parameter> Error: size(misfit_weight) must be nrefl + 1. Exit. ')
            if (any_in(['refl'], model_name)) then
                yn_update_reflector = .true.
            end if
            call readpar_float(file_parameter, 'offset_min_refl', offset_min_refl, 0.0)
            call readpar_float(file_parameter, 'offset_max_refl', offset_max_refl, float_huge)
        else
            call readpar_nfloat(file_parameter, 'misfit_weight', misfit_weight, ones(1))
        end if

    end subroutine read_parameters

    subroutine read_additional_parameters_tloc

        integer :: i
        real :: valmin, valmax

        do i = 1, nmodel
            if (any(model_name(i) == ['sx', 'sy', 'sz', 'st0'])) then
                select case (model_name(i))
                    case ('sx')
                        valmin = ox
                        valmax = ox + (nx - 1)*dx
                    case ('sy')
                        valmin = oy
                        valmax = oy + (ny - 1)*dy
                    case ('sz')
                        valmin = oz
                        valmax = oz + (nz - 1)*dz
                    case ('st0')
                        valmin = 0.0
                        valmax = 1.0e6
                end select
                call readpar_float(file_parameter, 'min_'//tidy(model_name(i)), model_min(i), valmin)
                call readpar_float(file_parameter, 'max_'//tidy(model_name(i)), model_max(i), valmax)
                ! The program allows using different search method for different parameters
                call readpar_string(file_parameter, 'search_method_'//tidy(model_name(i)), model_search_method(i), search_method)
            end if
        end do

    end subroutine read_additional_parameters_tloc

    !
    !> Load geometry
    !
    subroutine load_geometry

        integer :: i, j, l
        character(len=1024) :: shot_file_geometry, dir_geometry
        logical, allocatable, dimension(:) :: qs
        integer, allocatable, dimension(:) :: usid
        type(source_receiver_geometry), allocatable, dimension(:) :: g

        if (rankid == 0) then
            call warn(date_time_compact()//' Loading geometry... ')
        end if

        if (allocated(gmtr)) then
            deallocate (gmtr)
        end if
        allocate (gmtr(1:ns))
        if (.not. file_exists(file_geometry)) then
            call warn(date_time_compact()//' @'//get_hostname() //' Error: Geometry file does not exist. ')
            stop
        else
            dir_geometry = get_file_directory(file_geometry)
        end if

        open (3, file=tidy(file_geometry), action='read', status='old')
        do i = 1, ns

            read(3, *) shot_file_geometry

            open (4, file=tidy(dir_geometry)//tidy(shot_file_geometry), action='read', status='old')

            ! Read geometry ID
            read(4, *) gmtr(i)%id

            ! Read source parameters
            read(4, *) gmtr(i)%ns
            allocate(gmtr(i)%srcr(1:gmtr(i)%ns))
            do j = 1, gmtr(i)%ns
                read (4, *) gmtr(i)%srcr(j)%x, gmtr(i)%srcr(j)%y, gmtr(i)%srcr(j)%z, gmtr(i)%srcr(j)%t0
            end do

            ! Read receiver parameters
            read(4, *) gmtr(i)%nr
            allocate(gmtr(i)%recr(1:gmtr(i)%nr))
            do j = 1, gmtr(i)%nr
                read (4, *) gmtr(i)%recr(j)%x, gmtr(i)%recr(j)%y, gmtr(i)%recr(j)%z, gmtr(i)%recr(j)%weight
            end do

            close(4)

            if (ny == 1) then
                gmtr(i)%srcr(:)%y = 0.5*(ymin + ymax)
                gmtr(i)%recr(:)%y = 0.5*(ymin + ymax)
            end if

            if ((mod(i, max(nint(ns/10.0), 1)) == 0 .or. i == ns .or. i == 1) .and. &
                    rankid == 0) then
                call warn(date_time_compact()//' Loading geometry '//num2str(i)//' of '//num2str(ns))
            end if

        end do
        close(3)

        ! Check validity of geometry
        do i = 1, ns

            ! Check source
            l = 0
            do j =  1, gmtr(i)%ns
                if ( &
                        gmtr(i)%srcr(j)%x >= xmin .and. gmtr(i)%srcr(j)%x <= xmax .and. &
                        gmtr(i)%srcr(j)%y >= ymin .and. gmtr(i)%srcr(j)%y <= ymax .and. &
                        gmtr(i)%srcr(j)%z >= zmin .and. gmtr(i)%srcr(j)%z <= zmax .and. &
                        gmtr(i)%srcr(j)%x >= sxmin .and. gmtr(i)%srcr(j)%x <= sxmax .and. &
                        gmtr(i)%srcr(j)%y >= symin .and. gmtr(i)%srcr(j)%y <= symax .and. &
                        gmtr(i)%srcr(j)%z >= szmin .and. gmtr(i)%srcr(j)%z <= szmax .and. &
                        gmtr(i)%id >= sid_min .and. gmtr(i)%id <= sid_max) then
                    l = l + 1
                else
                    gmtr(i)%srcr(j)%amp = 0
                end if
            end do

            if (l == 0) then
                gmtr(i)%ns = 0
            end if

            ! Check receiver
            l = 0
            do j = 1, gmtr(i)%nr
                ! Compute absolute offset
                gmtr(i)%recr(j)%aoff = sqrt( &
                    (gmtr(i)%recr(j)%x - mean(gmtr(i)%srcr(:)%x))**2 &
                    + (gmtr(i)%recr(j)%y - mean(gmtr(i)%srcr(:)%y))**2 &
                    + (gmtr(i)%recr(j)%z - mean(gmtr(i)%srcr(:)%z))**2)

                ! Check if receiver in the range
                if ( &
                        gmtr(i)%recr(j)%x >= xmin .and. gmtr(i)%recr(j)%x <= xmax .and. &
                        gmtr(i)%recr(j)%y >= ymin .and. gmtr(i)%recr(j)%y <= ymax .and. &
                        gmtr(i)%recr(j)%z >= zmin .and. gmtr(i)%recr(j)%z <= zmax .and. &
                        gmtr(i)%recr(j)%x >= rxmin .and. gmtr(i)%recr(j)%x <= rxmax .and. &
                        gmtr(i)%recr(j)%y >= rymin .and. gmtr(i)%recr(j)%y <= rymax .and. &
                        gmtr(i)%recr(j)%z >= rzmin .and. gmtr(i)%recr(j)%z <= rzmax .and. &
                        gmtr(i)%recr(j)%aoff >= offset_min .and. gmtr(i)%recr(j)%aoff <= offset_max .and. &
                        j >= rec_min .and. j <= rec_max .and. mod(j - rec_min, rec_every) == 0 .and. &
                        gmtr(i)%recr(j)%weight /= 0) then
                    l = l + 1
                else
                    gmtr(i)%recr(j)%weight = 0.0
                end if

            end do
            close (4)

            ! check if any receiver in the model
            if (l == 0) then
                gmtr(i)%nr = 0
            end if

        end do

        ! Remove common-shot gathers that are not qualified
        qs = falses(ns)
        ! Select sources based on offset, gain and representational id restrictions
        do i = shot_min, shot_max, shot_every
            if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0) then
                qs(i) = .true.
            end if
        end do
        if (allocated(src_select)) then
            if (any(src_select /= 0)) then
                ! Select sources based on their representional id
                qs = .false.
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. any(src_select == i)) then
                        qs(i) = .true.
                    end if
                end do
            end if
        end if
        if (allocated(sid_select)) then
            if (any(sid_select /= 0)) then
                ! Select source based on their field id
                qs = .false.
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. any(sid_select == gmtr(i)%id)) then
                        qs(i) = .true.
                    end if
                end do
            end if
        end if
        if (allocated(src_exclude)) then
            if (any(src_exclude /= 0)) then
                ! Exclude sources based on their representional id
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. qs(i) .and. any(src_exclude == i)) then
                        qs(i) = .false.
                    end if
                end do
            end if
        end if
        if (allocated(sid_exclude)) then
            if (any(sid_exclude /= 0)) then
                ! Exclude source based on their field id
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. qs(i) .and. any(sid_exclude == gmtr(i)%id)) then
                        qs(i) = .false.
                    end if
                end do
            end if
        end if

        ! If no qualified shots, then stop
        if (.not. any(qs)) then
            call warn(date_time_compact()//' @'//get_hostname()//' Error: No shot or receiver in the model. ')
            stop
        end if

        ! Select all qualified shots
        allocate(g(1:ns))
        l = 1
        do i = 1, ns
            if (qs(i)) then
                g(l) = gmtr(i)
                l = l + 1
            end if
        end do

        gmtr = g(1:l - 1)

        ! The following pack(gmtr) works with ifort, but not new ifx...
        ! To avoid potential segmentation fault, use the above hard way to pack
        !        gmtr = pack(gmtr, mask=qs)
        ns = size(gmtr)

        ! If any duplicate source id, then stop
        if (ns >= 2) then
            call alloc_array(sid_select, [1, ns])
            do i = 1, ns
                sid_select(i) = gmtr(i)%id
            end do
            usid = unique(sid_select)
            if (size(usid) < ns) then
                call warn(date_time_compact()//' <load_geometry> Erorr: Duplicate source ID found. ')
                stop
            end if
        end if

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' Geometry is loaded with '//num2str(ns)//' sources. ')
        end if
        call mpibarrier

    end subroutine load_geometry

    !
    !> Get and set model dimensions
    !
    subroutine set_regular_space

        ! Original model dimensions
        call readpar_int(file_parameter, 'nx', nx0, 0, required=.true.)
        call readpar_int(file_parameter, 'ny', ny0, 1)
        call readpar_int(file_parameter, 'nz', nz0, 0, required=.true.)
        call readpar_float(file_parameter, 'dx', dx0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'dy', dy0, 1.0)
        call readpar_float(file_parameter, 'dz', dz0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'ox', ox0, 0.0)
        call readpar_float(file_parameter, 'oy', oy0, 0.0)
        call readpar_float(file_parameter, 'oz', oz0, 0.0)

        ! Target model dimensions
        call readpar_int(file_parameter, 'nnx', nx, nx0)
        call readpar_int(file_parameter, 'nny', ny, ny0)
        call readpar_int(file_parameter, 'nnz', nz, nz0)
        call readpar_float(file_parameter, 'ddx', dx, dx0)
        call readpar_float(file_parameter, 'ddy', dy, dy0)
        call readpar_float(file_parameter, 'ddz', dz, dz0)
        call readpar_float(file_parameter, 'oox', ox, ox0)
        call readpar_float(file_parameter, 'ooy', oy, oy0)
        call readpar_float(file_parameter, 'ooz', oz, oz0)

        ! Target model ranges
        call readpar_float(file_parameter, 'xmin', xmin, ox)
        call readpar_float(file_parameter, 'xmax', xmax, ox + (nx - 1)*dx)
        call assert(xmax >= xmin, ' <set_regular_space> Error: xmax must >= xmin')
        call readpar_float(file_parameter, 'ymin', ymin, oy)
        call readpar_float(file_parameter, 'ymax', ymax, oy + (ny - 1)*dy)
        call assert(ymax >= ymin, ' <set_regular_space> Error: ymax must >= ymin')
        call readpar_float(file_parameter, 'zmin', zmin, oz)
        call readpar_float(file_parameter, 'zmax', zmax, oz + (nz - 1)*dz)
        call assert(zmax >= zmin, ' <set_regular_space> Error: zmax must >= zmin')

        ! Target model origins
        ox = xmin
        oy = ymin
        oz = zmin
        nx = ceiling((xmax - xmin)/dx) + 1
        ny = ceiling((ymax - ymin)/dy) + 1
        nz = ceiling((zmax - zmin)/dz) + 1

        ! Check if requires model interpolation
        if (nx /= nx0 .or. ny /= ny0 .or. nz /= nz0 .or. &
                dx /= dx0 .or. dy /= dy0 .or. dz /= dz0 .or. &
                ox /= ox0 .or. oy /= oy0 .or. oz /= oz0) then
            require_model_interp = .true.
        end if

        ! Check non-negativeness
        call assert(nx >= 1 .and. ny >= 1 .and. nz >= 1, ' <set_regular_space> Error: n# must >= 1')
        call assert(dx > 0 .and. dy > 0 .and. dz > 0, ' <set_regular_space> Error: d# must > 0')

        ! Print original and target model dimensions
        if (rankid == 0) then
            if (ny == 1) then
                call warn(' ')
                call warn(' Original model dimensions: ')
                call warn('      nx, nz = '//num2str(nx0)//', '//num2str(nz0))
                call warn('      dx, dz = '//num2str(dx0, '(es)')//', '//num2str(dz0, '(es)'))
                call warn('      ox, oz = '//num2str(ox0, '(es)')//', '//num2str(oz0, '(es)'))
                call warn('      ex, ez = ' &
                    //num2str(ox0 + (nx0 - 1)*dx0, '(es)')//', ' &
                    //num2str(oz0 + (nz0 - 1)*dz0, '(es)'))
                call warn(' ')
                call warn(' Target model dimensions: ')
                call warn('      nx, nz = '//num2str(nx)//', '//num2str(nz))
                call warn('      dx, dz = '//num2str(dx, '(es)')//', '//num2str(dz, '(es)'))
                call warn('      ox, oz = '//num2str(ox, '(es)')//', '//num2str(oz, '(es)'))
                call warn('      ex, ez = ' &
                    //num2str(ox + (nx - 1)*dx, '(es)')//', ' &
                    //num2str(oz + (nz - 1)*dz, '(es)'))
                call warn(' ')
            else
                call warn(' ')
                call warn(' Original model dimensions: ')
                call warn('      nx, ny, nz = '//num2str(nx0)//', '//num2str(ny0)//', '//num2str(nz0))
                call warn('      dx, dy, dz = '//num2str(dx0, '(es)')//', '//num2str(dy0, '(es)')//', '//num2str(dz0, '(es)'))
                call warn('      ox, oy, oz = '//num2str(ox0, '(es)')//', '//num2str(oy0, '(es)')//', '//num2str(oz0, '(es)'))
                call warn('      ex, ey, ez = ' &
                    //num2str(ox0 + (nx0 - 1)*dx0, '(es)')//', ' &
                    //num2str(oy0 + (ny0 - 1)*dy0, '(es)')//', ' &
                    //num2str(oz0 + (nz0 - 1)*dz0, '(es)'))
                call warn(' ')
                call warn(' Target model dimensions: ')
                call warn('      nx, ny, nz = '//num2str(nx)//', '//num2str(ny)//', '//num2str(nz))
                call warn('      dx, dy, dz = '//num2str(dx, '(es)')//', '//num2str(dy, '(es)')//', '//num2str(dz, '(es)'))
                call warn('      ox, oy, oz = '//num2str(ox, '(es)')//', '//num2str(oy, '(es)')//', '//num2str(oz, '(es)'))
                call warn('      ex, ey, ez = ' &
                    //num2str(ox + (nx - 1)*dx, '(es)')//', ' &
                    //num2str(oy + (ny - 1)*dy, '(es)')//', ' &
                    //num2str(oz + (nz - 1)*dz, '(es)'))
                call warn(' ')
            end if
        end if

    end subroutine set_regular_space

    !
    !> Set adaptive range
    !
    subroutine set_adaptive_model_range(geom)

        type(source_receiver_geometry), intent(inout) :: geom

        ! range selection
        if (yn_adpx) then
            shot_xbeg = max(xmin, min( &
                minval(geom%recr(:)%x, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%x)) - adpextrax)
            shot_xend = min(xmax, max( &
                maxval(geom%recr(:)%x, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%x)) + adpextrax)
            shot_xbeg = floor((shot_xbeg - ox)/dx)*dx + ox
            shot_xend = ceiling((shot_xend - ox)/dx)*dx + ox
        else
            shot_xbeg = ox
            shot_xend = ox + (nx - 1)*dx
        end if
        where (geom%srcr(:)%x < shot_xbeg) geom%srcr(:)%amp = 0.0
        where (geom%srcr(:)%x > shot_xend) geom%srcr(:)%amp = 0.0
        where (geom%recr(:)%x < shot_xbeg) geom%recr(:)%weight = 0.0
        where (geom%recr(:)%x > shot_xend) geom%recr(:)%weight = 0.0

        if (yn_adpy) then
            shot_ybeg = max(ymin, min( &
                minval(geom%recr(:)%y, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%y)) - adpextray)
            shot_yend = min(ymax, max( &
                maxval(geom%recr(:)%y, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%y)) + adpextray)
            shot_ybeg = floor((shot_ybeg - oy)/dy)*dy + oy
            shot_yend = ceiling((shot_yend - oy)/dy)*dy + oy
        else
            shot_ybeg = oy
            shot_yend = oy + (ny - 1)*dy
        end if
        where (geom%srcr(:)%y < shot_ybeg) geom%srcr(:)%amp = 0.0
        where (geom%srcr(:)%y > shot_yend) geom%srcr(:)%amp = 0.0
        where (geom%recr(:)%y < shot_ybeg) geom%recr(:)%weight = 0.0
        where (geom%recr(:)%y > shot_yend) geom%recr(:)%weight = 0.0

        if (yn_adpz) then
            shot_zbeg = max(zmin, min( &
                minval(geom%recr(:)%z, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%z)) - adpextraz)
            shot_zend = min(zmax, max( &
                maxval(geom%recr(:)%z, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%z)) + adpextraz)
            shot_zbeg = floor((shot_zbeg - oz)/dz)*dz + oz
            shot_zend = ceiling((shot_zend - oz)/dz)*dz + oz
        else
            shot_zbeg = oz
            shot_zend = oz + (nz - 1)*dz
        end if
        where (geom%srcr(:)%z < shot_zbeg) geom%srcr(:)%amp = 0.0
        where (geom%srcr(:)%z > shot_zend) geom%srcr(:)%amp = 0.0
        where (geom%recr(:)%z < shot_zbeg) geom%recr(:)%weight = 0.0
        where (geom%recr(:)%z > shot_zend) geom%recr(:)%weight = 0.0

        ! Range for each shot
        shot_nxbeg = clip(int((shot_xbeg - ox)/dx) + 1, 1, nx)
        shot_nxend = clip(int((shot_xend - ox)/dx) + 1, 1, nx)
        shot_nybeg = clip(int((shot_ybeg - oy)/dy) + 1, 1, ny)
        shot_nyend = clip(int((shot_yend - oy)/dy) + 1, 1, ny)
        shot_nzbeg = clip(int((shot_zbeg - oz)/dz) + 1, 1, nz)
        shot_nzend = clip(int((shot_zend - oz)/dz) + 1, 1, nz)

        shot_nx = shot_nxend - shot_nxbeg + 1
        shot_ny = shot_nyend - shot_nybeg + 1
        shot_nz = shot_nzend - shot_nzbeg + 1

        if (verbose) then
            call warn(date_time_compact()//' Source '//num2str(gmtr(ishot)%id)//' dimensions:')
            call warn(date_time_compact()//'      xmin, xmax = ' &
                //num2str(shot_xbeg, '(es)')//', '//num2str(shot_xend, '(es)'))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      ymin, ymax = ' &
                    //num2str(shot_ybeg, '(es)')//', '//num2str(shot_yend, '(es)'))
            end if
            call warn(date_time_compact()//'      zmin, zmax = ' &
                //num2str(shot_zbeg, '(es)')//', '//num2str(shot_zend, '(es)'))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      nx, ny, nz = ' &
                    //num2str(shot_nx)//', '//num2str(shot_ny)//', '//num2str(shot_nz))
            else
                call warn(date_time_compact()//'      nx, nz = ' &
                    //num2str(shot_nx)//', '//num2str(shot_nz))
            end if
            call warn(date_time_compact()//'      nxmin, nxmax = ' &
                //num2str(shot_nxbeg)//', '//num2str(shot_nxend))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      nymin, nymax = ' &
                    //num2str(shot_nybeg)//', '//num2str(shot_nyend))
            end if
            call warn(date_time_compact()//'      nzmin, nzmax = ' &
                //num2str(shot_nzbeg)//', '//num2str(shot_nzend))
        end if

    end subroutine set_adaptive_model_range

    function source_receiver_coincide(s, r, threshold) result(y)

        type(source) :: s
        type(receiver) :: r
        real :: threshold
        logical :: y

        y = norm2([s%x, s%y, s%z] - [r%x, r%y, r%z]) <= threshold

    end function source_receiver_coincide

end module parameters
