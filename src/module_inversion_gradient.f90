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


module inversion_gradient

    use parameters
    use vars
    use gradient
    use inversion_regularization

    implicit none

contains

    !
    !> Compute gradient by adjoint-state method
    !
    subroutine compute_gradient

        ! Initialization
        call zero_gradient

        ! Compute gradients
        yn_misfit_only = .false.
        call compute_gradient_shots
        call mpibarrier

        ! Process and regularize gradients
        call process_gradient
        call regularize_gradient
        call output_gradient
        call mpibarrier

        ! Print progress
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Gradient computation completed. ')
        end if

    end subroutine compute_gradient

    !
    !> Compute reflector image by consistent time imaging condition
    !
    subroutine compute_reflector_image

        call compute_image_shots

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Reflector image computation completed. ')
        end if

    end subroutine compute_reflector_image

end module inversion_gradient
