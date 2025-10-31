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


module traveltime_iso_reflection

    use libflit
    use parameters
    use traveltime_iso

    implicit none

    public :: forward_iso_reflection
    public :: adjoint_iso_reflection
    public :: image_iso

    public :: forward_iso_reflection_elastic
    public :: adjoint_iso_reflection_elastic
    public :: image_iso_elastic

contains

    !===============================================================
    ! Acoustic domain
    !

    !
    !> Fast sweeping factorized eikonal solver in 2D isotropic acoustic media
    !
    subroutine forward_iso_reflection(v, d, o, geom, refl, tall, trec)

        real, dimension(:, :), intent(in) :: v
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: refl
        real, allocatable, dimension(:, :, :), intent(out) :: tall
        real, allocatable, dimension(:, :), intent(out) :: trec

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: tt0, tt, ttrec
        integer :: nrefl, i, j, l, h, mx, mz !imin,

        nrefl = int(maxval(refl))
        trec = zeros(geom%nr, nrefl + 1)
        mz = size(v, 1)
        mx = size(v, 2)
        !        t = zeros(mz, mx) + float_huge
        tall = zeros(mz, mx, nrefl + 1)

        ! Compute transmission traveltime
        sg = geom
        call forward_iso(v, d, o, sg, tt0, ttrec)
        trec(:, 1) = ttrec(:, 1)
        tall(:, :, 1) = tt0(:, :)
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission traveltime computation completed. ')

        ! Compute reflection traveltime
        do l = 1, nrefl

            sg = geom
            sg%ns = count(nint(refl) == l)
            deallocate(sg%srcr)
            allocate(sg%srcr(1:sg%ns))
            h = 1
            do j = 1, mx
                do i = 1, mz
                    if (nint(refl(i, j)) == l) then
                        sg%srcr(h)%z = (i - 1)*dz
                        sg%srcr(h)%y = 0
                        sg%srcr(h)%x = (j - 1)*dx
                        sg%srcr(h)%t0 = tt0(i, j)
                        h = h + 1
                    end if
                end do
            end do

            call forward_iso(v, d, o, sg, tt, ttrec)
            trec(:, l + 1) = ttrec(:, 1)
            tall(:, :, l + 1) = tt(:, :)

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' traveltime computation completed. ')

        end do

    end subroutine forward_iso_reflection

    !
    !> @brife Compute adjoint field for TRTT in 2D isotropic acoustic media
    !
    subroutine adjoint_iso_reflection(v, d, o, geom, tall, refl, tresidual, tadj)

        real, dimension(:, :), intent(in) :: v, refl
        real, dimension(:, :, :), intent(in) :: tall
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :, :), allocatable, intent(out) :: tadj

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: ta, tr
        integer :: nrefl, i, j, l, h, mx, mz
        integer :: mf = 10

        nrefl = int(maxval(refl))
        mz = size(v, 1)
        mx = size(v, 2)
        tadj = zeros(mz, mx, nrefl + 1)

        ! Compute receiver-side transmission adjoint field
        sg = geom
        call adjoint_iso(v, d, o, sg, tall(:, :, 1), tresidual(:, 1:1), ta)
        tadj(:, :, 1) = ta
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission adjoint field computation completed. ')

        ! Compute reflection adjoint field
        do l = 1, nrefl

            ! Receiver-side reflection adjoint field
            sg = geom
            call adjoint_iso(v, d, o, sg, tall(:, :, l + 1), tresidual(:, l + 1:l + 1), ta)
            ta = median_filt(ta)
            tadj(:, :, l + 1) = ta

            ! Source-side reflection adjoint field
            ! Allocate a much longer array, because along z, a refletor can quasi-vertical,
            ! resulting in more pixels in total than mx
            sg%nr = mf*mx
            deallocate(sg%recr)
            allocate(sg%recr(1:sg%nr))
            tr = zeros(mf*mx, 1)
            h = 1
            do j = 1, mx
                do i = 1, mz
                    if (nint(refl(i, j)) == l) then
                        ! The locations of effective receivers are reflectors
                        sg%recr(h)%z = (i - 1)*dz
                        sg%recr(h)%y = 0
                        sg%recr(h)%x = (j - 1)*dx
                        ! The residual is the value of receiver-side reflection adjoint field
                        tr(h, 1) = ta(i, j)
                        h = h + 1
                    end if
                end do
            end do
            sg%nr = h - 1
            sg%recr = sg%recr(1:sg%nr)
            tr = tr(1:sg%nr, :)
            call adjoint_iso(v, d, o, sg, tall(:, :, 1), tr, ta)
            tadj(:, :, l + 1) = tadj(:, :, l + 1) + ta

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' adjoint field computation completed. ')

        end do

    end subroutine adjoint_iso_reflection

    !
    !> Fast sweeping factorized eikonal solver in 2D isotropic acoustic media
    !
    subroutine image_iso(v, d, o, geom, tobs, refl)

        real, dimension(:, :), intent(in) :: v
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tobs
        real, allocatable, dimension(:, :, :), intent(out) :: refl

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: ttf, ttr, img
        real, allocatable, dimension(:, :) :: ttrec
        integer :: i, l
        integer :: n1, n2
        real :: tmax

        ! Compute transmission traveltime
        sg = geom
        call forward_iso(v, d, o, sg, ttf, ttrec)

        n1 = size(v, 1)
        n2 = size(v, 2)
        refl = zeros(n1, n2, nrefl)

        ! Compute reflection traveltime
        do l = 1, nrefl

            sg = geom
            sg%ns = sg%nr
            deallocate(sg%srcr)
            allocate(sg%srcr(1:sg%ns))
            ! Eikonal solver won't solve negative values, so here to backpropagate
            ! we have to first use the max of tobs to subtract tobs ...
            tmax = maxval(tobs(:, l))
            do i = 1, sg%ns
                sg%srcr(i)%x = sg%recr(i)%x
                sg%srcr(i)%y = 0
                sg%srcr(i)%z = sg%recr(i)%z
                sg%srcr(i)%t0 = tmax - tobs(i, l)
            end do
            call forward_iso(v, d, o, sg, ttr, ttrec)
            ! ... and then get it back
            ttr = tmax - ttr
            ! Get image by finding thresholded match
            img = rescale(abs(ttf - ttr), [0.0, 1.0])**4.0
            img = binarize(img, reflector_imaging_threshold, [1.0, 0.0])

            refl(:, :, l) = img

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' imaging completed. ')

        end do

    end subroutine image_iso

    !===============================================================
    ! Elastic domain
    !

    !
    !> Fast sweeping factorized eikonal solver in 2D isotropic elastic media
    !
    subroutine forward_iso_reflection_elastic(vp, vs, d, o, geom, refl, tpall, tsall, tprec, tsrec)

        real, dimension(:, :), intent(in) :: vp, vs
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: refl
        real, allocatable, dimension(:, :, :), intent(out) :: tpall, tsall
        real, allocatable, dimension(:, :), intent(out) :: tprec, tsrec

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: ttp0, ttp, ttprec, tts, ttsrec
        integer :: nrefl, i, j, l, h, mx, mz

        nrefl = int(maxval(refl))
        tprec = zeros(geom%nr, nrefl + 1)
        tsrec = zeros(geom%nr, nrefl + 1)
        mz = size(vp, 1)
        mx = size(vp, 2)
        tpall = zeros(mz, mx, nrefl + 1)
        tsall = zeros(mz, mx, nrefl + 1)

        ! Compute transmission traveltime
        sg = geom
        call forward_iso(vp, d, o, sg, ttp0, ttprec)
        tprec(:, 1) = ttprec(:, 1)
        tpall(:, :, 1) = ttp0(:, :)
        !        ! Based on the definition of TRTT, there should be no direct S (incident-P case) or direct P (incident-S case)
        !        call forward_iso(vs, d, o, sg, tts0, ttsrec)
        !        tsrec(:, 1) = ttsrec(:, 1)
        !        tsall(:, :, 1) = tts0(:, :)
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission traveltime computation completed. ')

        ! Compute reflection traveltime
        do l = 1, nrefl

            sg = geom
            sg%ns = count(nint(refl) == l)
            deallocate(sg%srcr)
            allocate(sg%srcr(1:sg%ns))
            h = 1
            do j = 1, mx
                do i = 1, mz
                    if (nint(refl(i, j)) == l) then
                        sg%srcr(h)%z = (i - 1)*dz
                        sg%srcr(h)%y = 0
                        sg%srcr(h)%x = (j - 1)*dx
                        ! The reflector FAT is with p, not s
                        sg%srcr(h)%t0 = ttp0(i, j)
                        h = h + 1
                    end if
                end do
            end do

            ! p, reflector time = tp, velocity = vp
            call forward_iso(vp, d, o, sg, ttp, ttprec)
            tprec(:, l + 1) = ttprec(:, 1)
            tpall(:, :, l + 1) = ttp(:, :)
            ! s, reflector time = tp, velocity = vs
            call forward_iso(vs, d, o, sg, tts, ttsrec)
            tsrec(:, l + 1) = ttsrec(:, 1)
            tsall(:, :, l + 1) = tts(:, :)

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' traveltime computation completed. ')

        end do

    end subroutine forward_iso_reflection_elastic

    !
    !> @brife Compute adjoint field for TRTT in 2D isotropic elastic media
    !
    subroutine adjoint_iso_reflection_elastic(vp, vs, d, o, geom, tpall, tsall, refl, tpresidual, tsresidual, tpadj, tsadj)

        real, dimension(:, :), intent(in) :: vp, vs, refl
        real, dimension(:, :, :), intent(in) :: tpall, tsall
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tpresidual, tsresidual
        real, dimension(:, :, :), allocatable, intent(out) :: tpadj, tsadj

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: tpa, tsa, tr
        integer :: nrefl, i, j, l, h, mx, mz
        integer :: mf = 10

        nrefl = int(maxval(refl))
        mz = size(vp, 1)
        mx = size(vs, 2)
        tpadj = zeros(mz, mx, nrefl + 1)
        tsadj = zeros(mz, mx, nrefl + 1)

        ! Compute receiver-side transmission adjoint field
        sg = geom
        call adjoint_iso(vp, d, o, sg, tpall(:, :, 1), tpresidual(:, 1:1), tpa)
        tpadj(:, :, 1) = tpa
        !        ! Based on the definition of TRTT, there should be no direct S (incident-P case) or direct P (incident-S case)
        !        call adjoint_iso(vs, d, o, sg, tsall(:, :, 1), tsresidual(:, 1:1), tsa)
        !        tsadj(:, :, 1) = tsa
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission adjoint field computation completed. ')

        ! Compute reflection adjoint field
        do l = 1, nrefl

            ! PP reflection
            ! Receiver-side reflection adjoint field -- for vp
            sg = geom
            call adjoint_iso(vp, d, o, sg, tpall(:, :, l + 1), tpresidual(:, l + 1:l + 1), tpa)
            tpa = median_filt(tpa)
            tpadj(:, :, l + 1) = tpa

            ! Source-side reflection adjoint field -- for vp
            ! Allocate a much longer array, because along z, a refletor can quasi-vertical,
            ! resulting in more pixels in total than mx
            sg%nr = mf*mx
            deallocate(sg%recr)
            allocate(sg%recr(1:sg%nr))
            tr = zeros(mf*mx, 1)
            h = 1
            do j = 1, mx
                do i = 1, mz
                    if (nint(refl(i, j)) == l) then
                        ! The locations of effective receivers are reflectors
                        sg%recr(h)%z = (i - 1)*dz
                        sg%recr(h)%y = 0
                        sg%recr(h)%x = (j - 1)*dx
                        ! The residual is the value of receiver-side reflection adjoint field
                        tr(h, 1) = tpa(i, j)
                        h = h + 1
                    end if
                end do
            end do
            sg%nr = h - 1
            sg%recr = sg%recr(1:sg%nr)
            tr = tr(1:sg%nr, :)
            ! The base traveltime field here is the transmission field
            call adjoint_iso(vp, d, o, sg, tpall(:, :, 1), tr, tpa)
            tpadj(:, :, l + 1) = tpadj(:, :, l + 1) + tpa

            ! PS reflection
            ! Receiver-side reflection adjoint field -- for vs
            sg = geom
            call adjoint_iso(vs, d, o, sg, tsall(:, :, l + 1), tsresidual(:, l + 1:l + 1), tsa)
            tsa = median_filt(tsa)
            tsadj(:, :, l + 1) = tsa

            !            ! Source-side reflection adjoint field -- for vp
            !            ! Allocate a much longer array, because along z, a refletor can quasi-vertical,
            !            ! resulting in more pixels in total than mx
            !            sg%nr = mf*mx
            !            deallocate(sg%recr)
            !            allocate(sg%recr(1:sg%nr))
            !            tr = zeros(mf*mx, 1)
            !            h = 1
            !            do j = 1, mx
            !                do i = 1, mz
            !                    if (nint(refl(i, j)) == l) then
            !                        ! The locations of effective receivers are reflectors
            !                        sg%recr(h)%z = (i - 1)*dz
            !                        sg%recr(h)%y = 0
            !                        sg%recr(h)%x = (j - 1)*dx
            !                        ! The residual is the value of receiver-side reflection adjoint field
            !                        tr(h, 1) = tsa(i, j)
            !                        h = h + 1
            !                    end if
            !                end do
            !            end do
            !            sg%nr = h - 1
            !            sg%recr = sg%recr(1:sg%nr)
            !            tr = tr(1:sg%nr, :)
            !            ! The base traveltime field here is the transmission field
            !            call adjoint_iso(vp, d, o, sg, tpall(:, :, 1), tr, tpa)
            !            tpadj(:, :, l + 1) = tpadj(:, :, l + 1) + tpa

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' adjoint field computation completed. ')

        end do

    end subroutine adjoint_iso_reflection_elastic

    !
    !> Fast sweeping factorized eikonal solver in 2D isotropic elastic media
    !
    subroutine image_iso_elastic(vp, vs, d, o, geom, tpobs, tsobs, refl)

        real, dimension(:, :), intent(in) :: vp, vs
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tpobs, tsobs
        real, allocatable, dimension(:, :, :), intent(out) :: refl

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: ttf, ttr, imgp, imgs
        real, allocatable, dimension(:, :) :: ttrec
        integer :: i, l
        integer :: n1, n2
        real :: tmax

        ! Compute transmission traveltime -- for p
        sg = geom
        call forward_iso(vp, d, o, sg, ttf, ttrec)

        n1 = size(vp, 1)
        n2 = size(vp, 2)
        refl = zeros(n1, n2, nrefl)

        ! Compute reflection traveltime -- for p and s
        do l = 1, nrefl

            sg = geom
            sg%ns = sg%nr
            deallocate(sg%srcr)
            allocate(sg%srcr(1:sg%ns))
            ! Eikonal solver won't solve negative values, so here to backpropagate
            ! we have to first use the max of tobs to subtract tobs ...
            tmax = maxval(tpobs(:, l))
            do i = 1, sg%ns
                sg%srcr(i)%x = sg%recr(i)%x
                sg%srcr(i)%y = 0
                sg%srcr(i)%z = sg%recr(i)%z
                sg%srcr(i)%t0 = tmax - tpobs(i, l)
            end do
            call forward_iso(vp, d, o, sg, ttr, ttrec)
            ! ... and then get it back
            ttr = tmax - ttr
            ! Get image by finding thresholded match
            imgp = rescale(abs(ttf - ttr), [0.0, 1.0])**4.0
            imgp = binarize(imgp, reflector_imaging_threshold, [1.0, 0.0])

            ! Eikonal solver won't solve negative values, so here to backpropagate
            ! we have to first use the max of tobs to subtract tobs ...
            tmax = maxval(tsobs(:, l))
            do i = 1, sg%ns
                sg%srcr(i)%x = sg%recr(i)%x
                sg%srcr(i)%y = 0
                sg%srcr(i)%z = sg%recr(i)%z
                sg%srcr(i)%t0 = tmax - tsobs(i, l)
            end do
            call forward_iso(vs, d, o, sg, ttr, ttrec)
            ! ... and then get it back
            ttr = tmax - ttr
            ! Get image by finding thresholded match
            imgs = rescale(abs(ttf - ttr), [0.0, 1.0])**4.0
            imgs = binarize(imgs, reflector_imaging_threshold, [1.0, 0.0])

            refl(:, :, l) = imgp + imgs

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' imaging completed. ')

        end do

    end subroutine image_iso_elastic

end module traveltime_iso_reflection
