
program test

    use libflit
    use geological_model_2d

    real, allocatable, dimension(:, :) :: vp, r
    integer :: n1, n2, ns, i, j
    real, allocatable, dimension(:) :: sx, sz, st0
    type(rgm2) :: p
    integer, allocatable, dimension(:) :: si
    logical, allocatable, dimension(:) :: x
    integer :: nrx, nrz

    n1 = 301
    n2 = 501

    call make_directory('./model')
    call make_directory('./geometry')

    p%n1 = n1
    p%n2 = n2
    p%nf = 17
    p%yn_regular_fault = .true.
    p%nl = 40
    p%unconf = 1
    p%unconf_z = 0.5
    p%disp = [10, 20]
    p%lwv = 0.1
    p%dip = [rand(range=[100.0, 120.0], seed=123), rand(range=[60.0, 80.0], seed=567)]
    p%disp = [10.0, -10.0]
    p%secondary_refl_amp = 0.1
    p%refl_smooth = 10.0
    p%yn_facies = .true.
    p%yn_rgt = .true.
    p%fwidth = 3
    p%seed = 111
    call p%generate

    vp = rescale(p%facies, [1000.0, 3000.0])
    call output_array(vp, './model/vp.bin')

    vp = 1.0/gauss_filt(1.0/vp, [4.0, 6.0])
    call output_array(vp, './model/vp_init.bin')

    sx = meshgrid([n1, n2], [10.0, 10.0], [0.0, 0.0], dim=2)
    sz = meshgrid([n1, n2], [10.0, 10.0], [0.0, 0.0], dim=1)
    x = falses(n1*n2)
    l = 1
    do j = 1, n2
        do i = 1, n1
            if (p%fault(i, j) /= 0 .and. j >= 100 .and. j <= n2 - 100) then
                x(l) = .true.
            end if
            l = l + 1
        end do
    end do

    print *, count(x)

    sx = pack(sx, x)
    sz = pack(sz, x)

    ns = 1200
    si = irandom(ns, range=[1, size(sx)], seed=1234)

    sx = sx(si)
    sz = sz(si)

    open(3, file='./model/sxz.txt')
    do i = 1, size(sx)
        write(3, *) sx(i)/1000.0, sz(i)/1000.0
    end do
    close(3)

    r = zeros(n1, n2)
    do i = 1, ns
        r(nint(sz(i)/10.0 + 1), nint(sx(i)/10.0 + 1)) = 1.0
    end do
    call output_array(r, './model/source_location.bin')

    ! geometry
    st0 = random(ns, range=[0.0, 10.0], seed=789)

    call output_array(sx, './model/sx.bin')
    call output_array(sz, './model/sz.bin')
    call output_array(st0, './model/st0.bin')

    nrx = 50
    nrz = 30

    open (33, file='./model/rxz.txt')
    open (3, file='./geometry/geometry.txt')
    do i = 1, ns
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) sx(i), 0.0, sz(i), st0(i)
        write (4, *)
        write (4, *) nrx + nrz
        do j = 1, nrx
            write (4, *) (j - 1)*100.0 + 50.0, 0.0, 0.0, 1.0
            if(i == 1)  write (33, *) ((j - 1)*100.0 + 50.0)/1000.0, 20.0/1000.0
        end do
        do j = 1, nrz
            write (4, *) 2000.0, 0.0, (j - 1)*60.0 + 100.0, 1.0
           if(i == 1)  write (33, *) 2.0, ((j - 1)*60.0 + 100.0)/1000.0
        end do
        close (4)

    end do
    close (3)
    close(33)

    open (3, file='./geometry/geometry_no_t0.txt')
    do i = 1, size(sx)
        write (3, *) 'shot_'//num2str(i)//'_geometry_no_t0.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry_no_t0.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) sx(i), 0.0, sz(i), 0.0
        write (4, *)
        write (4, *) nrx + nrz
        do j = 1, nrx
            write (4, *) (j - 1)*100.0 + 50.0, 0.0, 0.0, 1.0
        end do
        do j = 1, nrz
            write (4, *) 2000.0, 0.0, (j - 1)*60.0 + 100.0, 1.0
        end do
        close (4)

    end do
    close (3)

    sx = zeros(ns) + mean(sx)
    sz = zeros(ns) + mean(sz)
    st0 = zeros(ns) + mean(st0)
    call output_array(sx, './model/sx_init.bin')
    call output_array(sz, './model/sz_init.bin')
    call output_array(st0, './model/st0_init.bin')

end program test
