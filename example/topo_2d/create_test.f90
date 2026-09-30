
program test

    use libflit
    use librgm

    type(rgm2_curved) :: p
    real, allocatable, dimension(:) :: tp
    type(fractal_noise_1d) :: xx
    real, allocatable, dimension(:, :) :: vp, vs, m
    integer :: n1, n2
    integer :: pml = 15
    integer, allocatable, dimension(:) :: t

    call make_directory('./model')
    call make_directory('./geometry')

    n1 = 101
    n2 = 301

    ! Topography
    xx%n1 = n2
    xx%seed = 1122
    xx%octaves = 7
    xx%periods1 = 7

    tp = xx%generate()
    call pad_array(tp, [pml + 1, pml + 1])

    tp = gauss_filt(tp, 1.0)
    tp = rescale(tp, [0.0, 300.0])

    ! Mask
    m = ones(n1, n2)
    t = zeros(n2)
    do j = 1, n2
        do i = 1, n1
            if (i - 1.0 < (maxval(tp) - tp(j))/10.0) then
                m(i, j) = 0.0
            end if
        end do
        do i = 1, n1
            if (m(i, j) == 1) then
                t(j) = i
                exit
            end if
        end do
    end do
    call output_array(m, './model/mask.bin')

    ! Models
    p%n1 = n1
    p%n2 = n2
    p%refl_shape = 'perlin'
    p%refl_shape_top = 'perlin'
    p%refl_smooth = 3
    p%refl_smooth_top = 0
    p%seed = 2341
    p%lwv = 0.5
    p%lwh = 0.5
    p%refl_height = [0, 80]
    p%refl_height_top = [0, 2]
    p%nl = 15
    p%disp = [20.0, 40.0]
    p%yn_elastic = .true.
    p%nf = 6
    call p%generate

    p%vp = rescale(p%vp, [2000.0, 4000.0])
    p%vs = rescale(p%vs, [1200.0, 2400.0])

    ! Increase the velocities by 2000 m/s per kilometer below the topography, so that the
    ! first arrivals reach deeper; vs changes by the same factor as vp
    do j = 1, n2
        do i = t(j), n1
            p%vs(i, j) = p%vs(i, j)*(1.0 + 2.0*(i - t(j))*10.0/p%vp(i, j))
            p%vp(i, j) = p%vp(i, j) + 2.0*(i - t(j))*10.0
        end do
    end do

    vp = 1.0/gauss_filt(1.0/p%vp, [6.0, 20.0])
    vs = 1.0/gauss_filt(1.0/p%vs, [6.0, 20.0])

    where (m == 0)
        p%vp = 0.0
        p%vs = 0.0
        vp = 0
        vs = 0
    end where

    call output_array(p%vp, './model/vp.bin')
    call output_array(p%vs, './model/vs.bin')

    call output_array(return_normal(p%vp/p%vs), './model/vpvs_ratio.bin')

    call output_array(vp, './model/vp_init.bin')
    call output_array(vs, './model/vs_init.bin')

    ! Geometry: 30 sources on the topographic surface, 100 m apart; every source records
    ! the receivers at all surface grid points
    open (3, file='./geometry/geometry.txt')
    do i = 1, 30
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'
    end do

    close (3)

    do ishot = 1, 30

        open (3, file='./geometry/shot_'//num2str(ishot)//'_geometry.txt')
        write (3, *) ishot
        write (3, *)
        write (3, *) 1

        j = nint((50.0 + (ishot - 1)*100.0)/10.0 + 1)

        write (3, '(3es, es)') 50.0 + (ishot - 1)*100.0, 0.0, (t(j) - 1.0)*10.0, 0.0
        write (3, *)

        write (3, *) p%n2
        do i = 1, p%n2
            write (3, '(3es, es)') (i - 1)*10.0, 0.0, (t(i) - 1.0)*10.0, 1.0
        end do

        close (3)
    end do

end program test
