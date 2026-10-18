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


module traveltime_iso_reflection

    use libflit
    use parameters
    use traveltime_iso

    implicit none

    public :: forward_iso_reflection
    public :: adjoint_iso_reflection
    public :: image_iso
    public :: adjoint_iso_reflection_at_reflector

    public :: forward_iso_reflection_elastic
    public :: adjoint_iso_reflection_elastic
    public :: image_iso_elastic

    private :: reflector_point, backpropagate_time, zero_crossings

contains

    !===============================================================
    ! Acoustic domain
    !

    !
    !> Fast sweeping factorized eikonal solver in 3D isotropic acoustic media
    !
    subroutine forward_iso_reflection(v, d, o, geom, refl, tall, trec)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :, :), intent(in) :: refl
        real, allocatable, dimension(:, :, :, :), intent(out) :: tall
        real, allocatable, dimension(:, :), intent(out) :: trec

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :, :) :: tt0, tt
        real, allocatable, dimension(:, :) :: ttrec
        integer :: nrefl, i, j, k, l, h, mx, my, mz
        !        imin,
        real :: zf, tf

        nrefl = nint(maxval(refl))
        trec = zeros(geom%nr, nrefl + 1)
        mz = size(v, 1)
        my = size(v, 2)
        mx = size(v, 3)
        !        t = zeros(mz, my, mx) + float_huge
        tall = zeros(mz, my, mx, nrefl + 1)

        ! Compute transmission traveltime
        sg = geom
        call forward_iso(v, d, o, sg, tt0, ttrec)
        trec(:, 1) = ttrec(:, 1)
        tall(:, :, :, 1) = tt0(:, :, :)
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission traveltime computation completed. ')

        ! Compute reflection traveltime
        do l = 1, nrefl

            sg = geom
            sg%ns = count(nint(refl) == l)
            deallocate(sg%srcr)
            allocate(sg%srcr(1:sg%ns))
            h = 1
            do k = 1, mx
                do j = 1, my
                    do i = 1, mz
                        if (nint(refl(i, j, k)) == l) then
                            ! The reflector lies refl - l grid cells below the grid point,
                            ! or above it when negative
                            call reflector_point(refl(i, j, k) - l, i, mz, tt0(:, j, k), zf, tf)
                            sg%srcr(h)%z = o(3) + zf*dz
                            sg%srcr(h)%y = o(2) + (j - 1)*dy
                            sg%srcr(h)%x = o(1) + (k - 1)*dx
                            sg%srcr(h)%t0 = tf
                            h = h + 1
                            cycle
                        end if
                    end do
                end do
            end do

            call forward_iso(v, d, o, sg, tt, ttrec)
            trec(:, l + 1) = ttrec(:, 1)
            tall(:, :, :, l + 1) = tt(:, :, :)

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' traveltime computation completed. ')

        end do

    end subroutine forward_iso_reflection

    !
    !> @brife Compute adjoint field for TRTT in 3D isotropic acoustic media
    !
    subroutine adjoint_iso_reflection(v, d, o, geom, tall, refl, tresidual, tadj, trefl)

        real, dimension(:, :, :), intent(in) :: v, refl
        real, dimension(:, :, :, :), intent(in) :: tall
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :, :, :), allocatable, intent(out) :: tadj
        !> The receiver-side adjoint field of the reflected leg of each reflector
        !> at the points of that reflector, zero elsewhere, as in
        !> adjoint_iso_reflection_at_reflector
        real, dimension(:, :, :, :), allocatable, intent(out), optional :: trefl

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :, :) :: ta
        real, allocatable, dimension(:, :) :: tr
        integer :: nrefl, i, j, k, l, h, mx, my, mz
        integer :: mf = 10

        nrefl = nint(maxval(refl))
        mz = size(v, 1)
        my = size(v, 2)
        mx = size(v, 3)
        tadj = zeros(mz, my, mx, nrefl + 1)
        if (present(trefl)) then
            trefl = zeros(mz, my, mx, max(nrefl, 1))
        end if

        ! Compute receiver-side transmission adjoint field
        sg = geom
        call adjoint_iso(v, d, o, sg, tall(:, :, :, 1), tresidual(:, 1:1), ta)
        tadj(:, :, :, 1) = ta
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission adjoint field computation completed. ')

        ! Compute reflection adjoint field
        do l = 1, nrefl

            ! Receiver-side reflection adjoint field
            sg = geom
            call adjoint_iso(v, d, o, sg, tall(:, :, :, l + 1), tresidual(:, l + 1:l + 1), ta)
            ta = median_filt(ta)
            tadj(:, :, :, l + 1) = ta

            ! Source-side reflection adjoint field
            ! Allocate a much longer array, because along z, a refletor can quasi-vertical,
            ! resulting in more pixels in total than mx
            sg%nr = mf*mx*my
            deallocate(sg%recr)
            allocate(sg%recr(1:sg%nr))
            tr = zeros(mf*mx*my, 1)
            h = 1
            do k = 1, mx
                do j = 1, my
                    do i = 1, mz
                        if (nint(refl(i, j, k)) == l) then
                            ! The locations of effective receivers are reflectors
                            sg%recr(h)%z = o(3) + (i - 1)*dz
                            sg%recr(h)%y = o(2) + (j - 1)*dy
                            sg%recr(h)%x = o(1) + (k - 1)*dx
                            ! The residual is the value of receiver-side reflection adjoint field
                            tr(h, 1) = ta(i, j, k)
                            if (present(trefl)) then
                                trefl(i, j, k, l) = ta(i, j, k)
                            end if
                            h = h + 1
                        end if
                    end do
                end do
            end do
            sg%nr = h - 1
            sg%recr = sg%recr(1:sg%nr)
            tr = tr(1:sg%nr, :)
            call adjoint_iso(v, d, o, sg, tall(:, :, :, 1), tr, ta)
            tadj(:, :, :, l + 1) = tadj(:, :, :, l + 1) + ta

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' adjoint field computation completed. ')

        end do

    end subroutine adjoint_iso_reflection

    !
    !> The receiver-side adjoint fields of the reflection legs at the points of
    !> each reflector (zero elsewhere), driven by tresidual: the reflection
    !> residuals, or the data weights of the misfit. tall(:, :, :, l + 1) is the
    !> traveltime field of the reflected leg of reflector l and v its velocity.
    !> It solves only these legs, one adjoint equation per reflector, and
    !> applies the same median filter as adjoint_iso_reflection.
    !
    subroutine adjoint_iso_reflection_at_reflector(v, d, o, geom, tall, refl, tresidual, trefl)

        real, dimension(:, :, :), intent(in) :: v, refl
        real, dimension(:, :, :, :), intent(in) :: tall
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :, :, :), allocatable, intent(out) :: trefl

        real, allocatable, dimension(:, :, :) :: ta
        integer :: nrefl, l

        nrefl = nint(maxval(refl))
        trefl = zeros(size(v, 1), size(v, 2), size(v, 3), max(nrefl, 1))
        do l = 1, nrefl
            call adjoint_iso(v, d, o, geom, tall(:, :, :, l + 1), tresidual(:, l + 1:l + 1), ta)
            ta = median_filt(ta)
            where (nint(refl) == l)
                trefl(:, :, :, l) = ta
            end where
        end do

    end subroutine adjoint_iso_reflection_at_reflector

    !
    !> Image the reflectors of one shot in 3D isotropic acoustic media.
    !>
    !> With the source-side first-arrival time Ts, the receiver-side time Tr_r
    !> of receiver r and its observed reflection traveltime tobs_r, the excess
    !> time
    !>     E(x) = Ts(x) - max_r [tobs_r - Tr_r(x)]
    !>          = min_r [Ts(x) + Tr_r(x) - tobs_r]
    !> is negative above a reflector and positive below it, so the reflector is
    !> the zero crossing of E along z. The receiver-side term comes from one
    !> eikonal solve (backpropagate_time). Where the receiver that attains the
    !> minimum lies at the edge of the receiver array, the zero set of E is
    !> that receiver's isochron and not a reflector. Such a crossing lies
    !> above the reflector, and the lower envelope over the shots removes it
    !> (build_reflector).
    !>
    !> zidx(j, k, l) is the fractional 1-based z index of the first zero
    !> crossing of reflector l in the column (j, k), j along y and k along x,
    !> or 0 where there is none.
    !
    subroutine image_iso(v, d, o, geom, tobs, zidx)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tobs
        real, allocatable, dimension(:, :, :), intent(out) :: zidx

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :, :) :: ttf, ttr
        real, allocatable, dimension(:, :) :: ttrec
        logical, allocatable, dimension(:) :: live
        integer :: l, n2, n3

        ! Source-side first-arrival time
        sg = geom
        call forward_iso(v, d, o, sg, ttf, ttrec)

        n2 = size(v, 2)
        n3 = size(v, 3)
        zidx = zeros(n2, n3, nrefl)

        do l = 1, nrefl

            ! Receivers with reflection data
            live = tobs(:, l) > 0
            if (count(live) < 2) then
                cycle
            end if

            ! Receiver-side time, and the first zero crossing of E in each column
            call backpropagate_time(v, d, o, geom, tobs(:, l), live, ttr)
            call zero_crossings(ttf - ttr, zidx(:, :, l))

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' imaging completed. ')

        end do

    end subroutine image_iso

    !===============================================================
    ! Elastic domain
    !

    !
    !> Fast sweeping factorized eikonal solver in 3D isotropic elastic media
    !
    subroutine forward_iso_reflection_elastic(vp, vs, d, o, geom, refl, tpall, tsall, tprec, tsrec)

        real, dimension(:, :, :), intent(in) :: vp, vs
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :, :), intent(in) :: refl
        real, allocatable, dimension(:, :, :, :), intent(out) :: tpall, tsall
        real, allocatable, dimension(:, :), intent(out) :: tprec, tsrec

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :, :) :: ttp0, ttp, tts
        real, allocatable, dimension(:, :) :: ttprec, ttsrec
        integer :: nrefl, i, j, k, l, h, mx, my, mz
        real :: zf, tf

        nrefl = nint(maxval(refl))
        tprec = zeros(geom%nr, nrefl + 1)
        tsrec = zeros(geom%nr, nrefl + 1)
        mz = size(vp, 1)
        my = size(vp, 2)
        mx = size(vp, 3)
        tpall = zeros(mz, my, mx, nrefl + 1)
        tsall = zeros(mz, my, mx, nrefl + 1)

        ! Compute transmission traveltime
        sg = geom
        call forward_iso(vp, d, o, sg, ttp0, ttprec)
        tprec(:, 1) = ttprec(:, 1)
        tpall(:, :, :, 1) = ttp0(:, :, :)
        !        ! Based on the definition of TRTT, there should be no direct S (incident-P case) or direct P (incident-S case)
        !        call forward_iso(vs, d, o, sg, tts0, ttsrec)
        !        tsrec(:, 1) = ttsrec(:, 1)
        !        tsall(:, :, :, 1) = tts0(:, :, :)
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission traveltime computation completed. ')

        ! Compute reflection traveltime
        do l = 1, nrefl

            sg = geom
            sg%ns = count(nint(refl) == l)
            deallocate(sg%srcr)
            allocate(sg%srcr(1:sg%ns))
            h = 1
            do k = 1, mx
                do j = 1, my
                    do i = 1, mz
                        if (nint(refl(i, j, k)) == l) then
                            ! The reflector FAT is with p, not s; the reflector lies refl - l
                            ! grid cells below the grid point, or above it when negative
                            call reflector_point(refl(i, j, k) - l, i, mz, ttp0(:, j, k), zf, tf)
                            sg%srcr(h)%z = o(3) + zf*dz
                            sg%srcr(h)%y = o(2) + (j - 1)*dy
                            sg%srcr(h)%x = o(1) + (k - 1)*dx
                            sg%srcr(h)%t0 = tf
                            h = h + 1
                            cycle
                        end if
                    end do
                end do
            end do

            ! p, reflector time = tp, velocity = vp
            call forward_iso(vp, d, o, sg, ttp, ttprec)
            tprec(:, l + 1) = ttprec(:, 1)
            tpall(:, :, :, l + 1) = ttp(:, :, :)
            ! s, reflector time = tp, velocity = vs
            call forward_iso(vs, d, o, sg, tts, ttsrec)
            tsrec(:, l + 1) = ttsrec(:, 1)
            tsall(:, :, :, l + 1) = tts(:, :, :)

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' traveltime computation completed. ')

        end do

    end subroutine forward_iso_reflection_elastic

    !
    !> @brife Compute adjoint field for TRTT in 3D isotropic elastic media
    !
    subroutine adjoint_iso_reflection_elastic(vp, vs, d, o, geom, tpall, tsall, refl, tpresidual, tsresidual, tpadj, tsadj, &
            tprefl, tsrefl)

        real, dimension(:, :, :), intent(in) :: vp, vs, refl
        real, dimension(:, :, :, :), intent(in) :: tpall, tsall
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tpresidual, tsresidual
        real, dimension(:, :, :, :), allocatable, intent(out) :: tpadj, tsadj
        !> The receiver-side adjoint fields of the reflected legs at the points
        !> of each reflector, zero elsewhere: tprefl of the same-wave
        !> reflection (velocity vp) and tsrefl of the converted one (velocity vs)
        real, dimension(:, :, :, :), allocatable, intent(out), optional :: tprefl, tsrefl

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :, :) :: tpa, tsa
        real, allocatable, dimension(:, :) :: tr
        integer :: nrefl, i, j, k, l, h, mx, my, mz
        integer :: mf = 10

        nrefl = nint(maxval(refl))
        mz = size(vp, 1)
        my = size(vp, 2)
        mx = size(vp, 3)
        tpadj = zeros(mz, my, mx, nrefl + 1)
        tsadj = zeros(mz, my, mx, nrefl + 1)
        if (present(tprefl)) then
            tprefl = zeros(mz, my, mx, max(nrefl, 1))
        end if
        if (present(tsrefl)) then
            tsrefl = zeros(mz, my, mx, max(nrefl, 1))
        end if

        ! Compute receiver-side transmission adjoint field
        sg = geom
        call adjoint_iso(vp, d, o, sg, tpall(:, :, :, 1), tpresidual(:, 1:1), tpa)
        tpadj(:, :, :, 1) = tpa
        !        ! Based on the definition of TRTT, there should be no direct S (incident-P case) or direct P (incident-S case)
        !        call adjoint_iso(vs, d, o, sg, tsall(:, :, :, 1), tsresidual(:, 1:1), tsa)
        !        tsadj(:, :, :, 1) = tsa
        call warn(date_time_compact()//' Shot '//num2str(geom%id) &
            //' transmission adjoint field computation completed. ')

        ! Compute reflection adjoint field
        do l = 1, nrefl

            ! PP reflection
            ! Receiver-side reflection adjoint field -- for vp
            sg = geom
            call adjoint_iso(vp, d, o, sg, tpall(:, :, :, l + 1), tpresidual(:, l + 1:l + 1), tpa)
            tpa = median_filt(tpa)
            tpadj(:, :, :, l + 1) = tpa

            ! Source-side reflection adjoint field -- for vp
            ! Allocate a much longer array, because along z, a refletor can quasi-vertical,
            ! resulting in more pixels in total than mx
            sg%nr = mf*mx*my
            deallocate(sg%recr)
            allocate(sg%recr(1:sg%nr))
            tr = zeros(mf*mx*my, 1)
            h = 1
            do k = 1, mx
                do j = 1, my
                    do i = 1, mz
                        if (nint(refl(i, j, k)) == l) then
                            ! The locations of effective receivers are reflectors
                            sg%recr(h)%z = o(3) + (i - 1)*dz
                            sg%recr(h)%y = o(2) + (j - 1)*dy
                            sg%recr(h)%x = o(1) + (k - 1)*dx
                            ! The residual is the value of receiver-side reflection adjoint field
                            tr(h, 1) = tpa(i, j, k)
                            if (present(tprefl)) then
                                tprefl(i, j, k, l) = tpa(i, j, k)
                            end if
                            h = h + 1
                        end if
                    end do
                end do
            end do
            sg%nr = h - 1
            sg%recr = sg%recr(1:sg%nr)
            tr = tr(1:sg%nr, :)
            call adjoint_iso(vp, d, o, sg, tpall(:, :, :, 1), tr, tpa)
            tpadj(:, :, :, l + 1) = tpadj(:, :, :, l + 1) + tpa

            ! PS reflection
            ! Receiver-side reflection adjoint field -- for vs
            sg = geom
            call adjoint_iso(vs, d, o, sg, tsall(:, :, :, l + 1), tsresidual(:, l + 1:l + 1), tsa)
            tsa = median_filt(tsa)
            tsadj(:, :, :, l + 1) = tsa
            if (present(tsrefl)) then
                where (nint(refl) == l)
                    tsrefl(:, :, :, l) = tsa
                end where
            end if

            !            ! Source-side reflection adjoint field -- for vp
            !            ! Allocate a much longer array, because along z, a refletor can quasi-vertical,
            !            ! resulting in more pixels in total than mx
            !            sg%nr = mf*mx*my
            !            deallocate(sg%recr)
            !            allocate(sg%recr(1:sg%nr))
            !            tr = zeros(mf*mx*my, 1)
            !            h = 1
            !            do k = 1, mx
            !                do j = 1, my
            !                    do i = 1, mz
            !                        if (nint(refl(i, j, k)) == l) then
            !                            ! The locations of effective receivers are reflectors
            !                            sg%recr(h)%z = (i - 1)*dz
            !                            sg%recr(h)%y = (j - 1)*dy
            !                            sg%recr(h)%x = (k - 1)*dx
            !                            ! The residual is the value of receiver-side reflection adjoint field
            !                            tr(h, 1) = tsa(i, j, k)
            !                            h = h + 1
            !                        end if
            !                    end do
            !                end do
            !            end do
            !            sg%nr = h - 1
            !            sg%recr = sg%recr(1:sg%nr)
            !            tr = tr(1:sg%nr, :)
            !            call adjoint_iso(vp, d, o, sg, tpall(:, :, :, 1), tr, tpa)
            !            tpadj(:, :, :, l + 1) = tpadj(:, :, :, l + 1) + tpa

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' adjoint field computation completed. ')

        end do

    end subroutine adjoint_iso_reflection_elastic

    !
    !> Image the reflectors of one shot in 3D isotropic elastic media. vp is
    !> the velocity of the incident wave and of the reflected leg of the
    !> same-wave reflections tpobs, and vs that of the reflected leg of the
    !> converted reflections tsobs: PP and PS for an incident P wave, or, with
    !> vp and vs swapped by the caller, SS and SP for an incident S wave. Each
    !> data set gives a zero crossing (see image_iso), and a column takes their
    !> mean where both exist, otherwise the one that exists.
    !
    subroutine image_iso_elastic(vp, vs, d, o, geom, tpobs, tsobs, zidx)

        real, dimension(:, :, :), intent(in) :: vp, vs
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tpobs, tsobs
        real, allocatable, dimension(:, :, :), intent(out) :: zidx

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :, :) :: ttf, ttr
        real, allocatable, dimension(:, :) :: ttrec, zp, zs
        logical, allocatable, dimension(:) :: live
        integer :: l, n2, n3

        ! Source-side first-arrival time of the incident wave
        sg = geom
        call forward_iso(vp, d, o, sg, ttf, ttrec)

        n2 = size(vp, 2)
        n3 = size(vp, 3)
        zidx = zeros(n2, n3, nrefl)
        zp = zeros(n2, n3)
        zs = zeros(n2, n3)

        do l = 1, nrefl

            ! Same-wave reflections
            live = tpobs(:, l) > 0
            zp = 0.0
            if (count(live) >= 2) then
                call backpropagate_time(vp, d, o, geom, tpobs(:, l), live, ttr)
                call zero_crossings(ttf - ttr, zp)
            end if

            ! Converted reflections
            live = tsobs(:, l) > 0
            zs = 0.0
            if (count(live) >= 2) then
                call backpropagate_time(vs, d, o, geom, tsobs(:, l), live, ttr)
                call zero_crossings(ttf - ttr, zs)
            end if

            where (zp > 0 .and. zs > 0)
                zidx(:, :, l) = 0.5*(zp + zs)
            elsewhere (zp > 0)
                zidx(:, :, l) = zp
            elsewhere
                zidx(:, :, l) = zs
            end where

            call warn(date_time_compact()//' Shot '//num2str(geom%id)//' reflector ' &
                //num2str(l)//' imaging completed. ')

        end do

    end subroutine image_iso_elastic

    !===============================================================
    ! Reflector helpers
    !

    !
    !> Position and start time of a reflector point. The point is stored at
    !> grid row i with the fractional offset dfrac in [-0.5, 0.5) in grid cells;
    !> zf is its z index (0-based, fractional) and tf the source-side time
    !> interpolated at that depth from the column tcol. With dfrac = 0, the
    !> point is the grid point itself and tf = tcol(i).
    !
    subroutine reflector_point(dfrac, i, mz, tcol, zf, tf)

        real, intent(in) :: dfrac
        integer, intent(in) :: i, mz
        real, dimension(:), intent(in) :: tcol
        real, intent(out) :: zf, tf

        integer :: k
        real :: w

        zf = i - 1 + dfrac
        zf = min(max(zf, 0.0), mz - 1.0)
        k = floor(zf) + 1
        if (k >= mz) then
            tf = tcol(mz)
        else
            w = zf - (k - 1)
            tf = (1.0 - w)*tcol(k) + w*tcol(k + 1)
        end if

    end subroutine reflector_point

    !
    !> Receiver-side time ttr(x) = max_r [tobs_r - Tr_r(x)] over the receivers
    !> in mask, from one eikonal solve with these receivers as source points
    !> that start at tmax - tobs_r, where tmax is the largest tobs_r of these
    !> receivers
    !
    subroutine backpropagate_time(v, d, o, geom, tobs, mask, ttr)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:), intent(in) :: tobs
        logical, dimension(:), intent(in) :: mask
        real, allocatable, dimension(:, :, :), intent(out) :: ttr

        type(source_receiver_geometry) :: sg
        real, allocatable, dimension(:, :) :: ttrec
        integer :: i, h
        real :: tmax

        sg = geom
        sg%ns = count(mask)
        deallocate(sg%srcr)
        allocate(sg%srcr(1:sg%ns))
        ! The eikonal solver takes nonnegative start times
        tmax = maxval(tobs, mask=mask)
        h = 1
        do i = 1, geom%nr
            if (mask(i)) then
                sg%srcr(h)%x = geom%recr(i)%x
                sg%srcr(h)%y = geom%recr(i)%y
                sg%srcr(h)%z = geom%recr(i)%z
                sg%srcr(h)%t0 = tmax - tobs(i)
                h = h + 1
            end if
        end do
        call forward_iso(v, d, o, sg, ttr, ttrec)
        ttr = tmax - ttr

    end subroutine backpropagate_time

    !
    !> The first zero crossing of e from negative to positive along z in each
    !> column, as a fractional 1-based index by linear interpolation, or 0
    !> without one
    !
    subroutine zero_crossings(e, z)

        real, dimension(:, :, :), intent(in) :: e
        real, dimension(:, :), intent(out) :: z

        integer :: i, j, k

        z = 0.0
        !$omp parallel do private(i, j, k)
        do k = 1, size(e, 3)
            do j = 1, size(e, 2)
                do i = 1, size(e, 1) - 1
                    if (e(i, j, k) < 0 .and. e(i + 1, j, k) >= 0) then
                        z(j, k) = i + e(i, j, k)/(e(i, j, k) - e(i + 1, j, k))
                        exit
                    end if
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine zero_crossings

end module traveltime_iso_reflection
