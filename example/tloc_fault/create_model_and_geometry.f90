
program test

    use libflit
    use librgm

    real, allocatable, dimension(:, :, :) :: w
    real, allocatable, dimension(:, :, :) :: vp, vs, r
    integer :: n1, n2, n3, ns, i, j, k
    real, allocatable, dimension(:) :: sx, sy, sz, st0
    type(rgm3) :: p
    integer, allocatable, dimension(:) :: si
    logical, allocatable, dimension(:) :: x
    integer :: nrx, nry

    n1 = 101
    n2 = 201
    n3 = 301

    call make_directory('./model')
    call make_directory('./geometry')

    vp = zeros(n1, n2, n3) + 2000
    vs = zeros(n1, n2, n3) + 2000/sqrt(3.0)
    r = zeros(n1, n2, n3)

    p%n1 = n1
    p%n2 = n2
    p%n3 = n3
    p%nf = 4
    p%nl = 20
    p%lwv = 0.5
    p%disp = [10, 20]
    p%lwv = 0.5
    p%ng = 4
    p%refl_shape = 'gaussian'
    p%refl_sigma2 = [200.0, 400.0]
    p%refl_sigma3 = [200.0, 400.0]
    p%secondary_refl_height_ratio = 0.0
    p%refl_smooth = 10.0
    p%yn_rgt = .true.
    p%fwidth = 2
    p%seed = 23342412
    call p%generate

    w = p%rgt
    w = gauss_filt(w, [4.0, 4.0, 4.0])
    vp = rescale(w, [1000.0, 3000.0])
    vs = rescale(w, [500.0, 2000.0])
    call output_array(vp, './model/vp_init.bin')
    call output_array(vs, './model/vs_init.bin')

    r = load('checkerboard.bin', n1, n2, n3)
    r = gauss_filt(r, [5.0, 4.0, 4.0])
    r = rescale(r, [-300.0, 300.0])
    vp = vp + r
    vs = vs + r/sqrt(3.0)
    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')
    
    call output_array(clip(p%fault, 0.0, 1.0), './model/fault.bin')
    call output_array(clip(p%fault_strike/180.0, 0.0, 1.0), './model/strike.bin')

    ! Create sources
    sx = meshgrid([n1, n2, n3], [10.0, 10.0, 10.0], [0.0, 0.0, 0.0], dim=3)
    sy = meshgrid([n1, n2, n3], [10.0, 10.0, 10.0], [0.0, 0.0, 0.0], dim=2)
    sz = meshgrid([n1, n2, n3], [10.0, 10.0, 10.0], [0.0, 0.0, 0.0], dim=1)
    x = falses(n1*n2*n3)
    l = 1
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                if (p%fault(i, j, k) /= 0) then
                    x(l) = .true.
                end if
                l = l + 1
            end do
        end do
    end do

    print *, count(x)

    sx = pack(sx, x)
    sy = pack(sy, x)
    sz = pack(sz, x)

    ns = 1200
    si = irandom(ns, range=[1, size(sx)], seed=1234)

    sx = sx(si)
    sy = sy(si)
    sz = sz(si)
    
    sx = clip(sx, 0.0, (n3 - 1.0)*10)
    sy = clip(sy, 0.0, (n2 - 1.0)*10)
    sz = clip(sz, 0.0, (n1 - 1.0)*10)
	
    r = zeros(n1, n2, n3)
    do i = 1, ns
        r(nint(sz(i)/10.0 + 1), nint(sy(i)/10.0 + 1), nint(sx(i)/10.0 + 1)) = 1.0
    end do
    call output_array(r, './model/source_location.bin')

    ! geometry
    st0 = random(ns, range=[0.0, 10.0], seed=789)

    call output_array(sx, './model/sx.bin')
    call output_array(sy, './model/sy.bin')
    call output_array(sz, './model/sz.bin')
    call output_array(st0, './model/st0.bin')

    nrx = 15
    nry = 10

    open (3, file='./geometry/geometry.txt')
    do i = 1, ns
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) sx(i), sy(i), sz(i), st0(i)
        write (4, *)
        write (4, *) nrx*nry
        do k = 1, nrx
            do j = 1, nry
                write (4, *) (k - 1)*200.0 + 100.0, (j - 1)*200.0 + 100.0, 0.0, 1.0
            end do
        end do
        close (4)

    end do
    close (3)

    open (3, file='./geometry/geometry_no_t0.txt')
    do i = 1, size(sx)
        write (3, *) 'shot_'//num2str(i)//'_geometry_no_t0.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry_no_t0.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) sx(i), sy(i), sz(i), 0.0
        write (4, *)
        write (4, *) nrx*nry
        do k = 1, nrx
            do j = 1, nry
                write (4, *) (k - 1)*200.0 + 100.0, (j - 1)*200.0 + 100.0, 0.0, 1.0
            end do
        end do
        close(4)

    end do
    close (3)

    sx = zeros(ns) + mean(sx)
    sy = zeros(ns) + mean(sy)
    sz = zeros(ns) + mean(sz)
    st0 = zeros(ns) + mean(st0)
    call output_array(sx, './model/sx_init.bin')
    call output_array(sy, './model/sy_init.bin')
    call output_array(sz, './model/sz_init.bin')
    call output_array(st0, './model/st0_init.bin')

    !    sz = zeros(ns) + 980.0
    !    call output_array(sz, './model/sz_init_deep.bin')

end program test
