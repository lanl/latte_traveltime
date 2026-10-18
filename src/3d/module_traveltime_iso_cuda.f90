!
! © 2026. Triad National Security, LLC. All rights reserved.
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

module traveltime_iso_cuda

    ! Fast sweeping on a GPU for forward_iso and adjoint_iso in 3D, with the
    ! CUDA code in traveltime_iso_cuda.cu. The GPU runs the same iterations as
    ! the CPU and gives the same fields.

    use iso_c_binding
    use libflit
    use parameters

    implicit none

    interface

        function latte_cuda_init(name, len) result(id) bind(c, name='latte_cuda_init')
            import :: c_char, c_int
            character(kind=c_char), dimension(*) :: name
            integer(kind=c_int), value :: len
            integer(kind=c_int) :: id
        end function latte_cuda_init

        subroutine latte_cuda_message(text, len) bind(c, name='latte_cuda_message')
            import :: c_char, c_int
            character(kind=c_char), dimension(*) :: text
            integer(kind=c_int), value :: len
        end subroutine latte_cuda_message

        function latte_cuda_forward_3d(nx, ny, nz, dx, dy, dz, vp, t0, px0, py0, pz0, tau, &
                huge_value, threshold, niter_max, niter, diff) result(status) bind(c, name='latte_cuda_forward_3d')
            import :: c_int, c_double
            integer(kind=c_int), value :: nx, ny, nz, niter_max
            real(kind=c_double), value :: dx, dy, dz, huge_value, threshold
            real(kind=c_double), dimension(*) :: vp, t0, px0, py0, pz0, tau
            integer(kind=c_int) :: niter
            real(kind=c_double) :: diff
            integer(kind=c_int) :: status
        end function latte_cuda_forward_3d

        function latte_cuda_adjoint_3d(nx, ny, nz, dx, dy, dz, tt, recrflag, src, use_src, lambda, &
                huge_value, threshold, niter_max, niter, diff) result(status) bind(c, name='latte_cuda_adjoint_3d')
            import :: c_int, c_double, c_signed_char
            integer(kind=c_int), value :: nx, ny, nz, niter_max, use_src
            real(kind=c_double), value :: dx, dy, dz, huge_value, threshold
            real(kind=c_double), dimension(*) :: tt, src, lambda
            integer(kind=c_signed_char), dimension(*) :: recrflag
            integer(kind=c_int) :: niter
            real(kind=c_double) :: diff
            integer(kind=c_int) :: status
        end function latte_cuda_adjoint_3d

    end interface

    ! Whether this process has looked for a GPU, and whether it found one
    logical :: gpu_checked = .false.
    logical :: gpu_found = .false.

    private
    public :: use_gpu
    public :: fast_sweep_forward_gpu
    public :: fast_sweep_adjoint_gpu

contains

    !
    !> Whether the fast sweeping of the forward and adjoint solvers runs on a
    !> GPU, i.e., whether the process has one. The first call selects the GPU.
    !
    function use_gpu() result(yn)

        logical :: yn

        character(len=256) :: name
        integer :: id

        if (.not. gpu_checked) then
            gpu_checked = .true.
            id = latte_cuda_init(name, len(name))
            gpu_found = id >= 0
            if (gpu_found) then
                call warn(date_time_compact()//' Rank '//num2str(rankid)// &
                    ' runs the fast sweeping on GPU '//num2str(id)//', '//tidy(name)//'. ')
            else
                call warn(date_time_compact()//' Rank '//num2str(rankid)// &
                    ' finds no GPU and runs the fast sweeping on the CPU. ')
            end if
        end if

        yn = gpu_found

    end function use_gpu

    !
    !> The fast sweeping iterations of forward_iso on the GPU. On success, tau
    !> is the final multiplicative time field, niter the number of iterations,
    !> and diff the mean change of tau in the last iteration. On a failure,
    !> e.g., too little GPU memory, tau, niter and diff stay unchanged.
    !
    subroutine fast_sweep_forward_gpu(nx, ny, nz, dx, dy, dz, vp, t0, px0, py0, pz0, tau, huge_value, niter, diff)

        integer, intent(in) :: nx, ny, nz
        double precision, intent(in) :: dx, dy, dz, huge_value
        double precision, dimension(nx, ny, nz), intent(in) :: vp, t0, px0, py0, pz0
        double precision, dimension(nx, ny, nz), intent(inout) :: tau
        integer, intent(inout) :: niter
        double precision, intent(inout) :: diff

        integer(kind=c_int) :: n
        real(kind=c_double) :: r

        if (latte_cuda_forward_3d(nx, ny, nz, dx, dy, dz, vp, t0, px0, py0, pz0, tau, &
                huge_value, dble(sweep_stop_threshold), sweep_niter_max, n, r) == 0) then
            niter = n
            diff = r
        else
            call warn_failure
        end if

    end subroutine fast_sweep_forward_gpu

    !
    !> The fast sweeping iterations of adjoint_iso on the GPU, with the same
    !> conventions as fast_sweep_forward_gpu; src holds the point sources of
    !> the adjoint-state equation, if any
    !
    subroutine fast_sweep_adjoint_gpu(nx, ny, nz, dx, dy, dz, tt, recrflag, lambda, huge_value, niter, diff, src)

        integer, intent(in) :: nx, ny, nz
        double precision, intent(in) :: dx, dy, dz, huge_value
        double precision, dimension(nx, ny, nz), intent(in) :: tt
        logical, dimension(nx, ny, nz), intent(in) :: recrflag
        double precision, dimension(nx, ny, nz), intent(inout) :: lambda
        integer, intent(inout) :: niter
        double precision, intent(inout) :: diff
        double precision, dimension(nx, ny, nz), intent(in), optional :: src

        integer(kind=c_signed_char), allocatable, dimension(:, :, :) :: flag
        integer(kind=c_int) :: n, status
        real(kind=c_double) :: r
        real(kind=c_double), dimension(1) :: nosrc

        flag = merge(1_c_signed_char, 0_c_signed_char, recrflag)

        if (present(src)) then
            status = latte_cuda_adjoint_3d(nx, ny, nz, dx, dy, dz, tt, flag, src, 1, lambda, &
                huge_value, dble(sweep_stop_threshold), sweep_niter_max, n, r)
        else
            nosrc = 0
            status = latte_cuda_adjoint_3d(nx, ny, nz, dx, dy, dz, tt, flag, nosrc, 0, lambda, &
                huge_value, dble(sweep_stop_threshold), sweep_niter_max, n, r)
        end if
        if (status == 0) then
            niter = n
            diff = r
        else
            call warn_failure
        end if

    end subroutine fast_sweep_adjoint_gpu

    !
    !> Report why the GPU failed
    !
    subroutine warn_failure

        character(len=512) :: text

        call latte_cuda_message(text, len(text))
        call warn(date_time_compact()//' Rank '//num2str(rankid)//' GPU fast sweeping failed: ' &
            //tidy(text)//'. Fast sweeping runs on the CPU instead. ')

    end subroutine warn_failure

end module traveltime_iso_cuda
