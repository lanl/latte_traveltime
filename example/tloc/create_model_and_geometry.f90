
program test

    use libflit

    real, allocatable, dimension(:, :) :: vp, vs, r
    integer :: n1, n2, ns, i
    real, allocatable, dimension(:) :: sx, sz, st0

    n1 = 301
    n2 = 501
    ns = 300

    call make_directory('./model')
    call make_directory('./geometry')

    vp = zeros(n1, n2) + 2000
    vs = zeros(n1, n2) + 2000/sqrt(3.0)
    r = zeros(n1, n2)

    call output_array(vp, './model/vp_init.bin')
    call output_array(vs, './model/vs_init.bin')

    r = zeros(n1, n2)
    do j = 1, n2
        do i = 1, n1
            r(i, j) = sin(0.05*i)*cos(0.055*j)
        end do
    end do

    r = r(1:n1, 1:n2)
    r = gauss_filt(r, [7.0, 7.0])
    r = rescale(r, [-200.0, 200.0])
    vp = vp + r
    vs = vs + r/3.0

    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')

    sx = random(ns, range=[1000.0, 4000.0], seed=123)
    sz = random(ns, range=[100.0, 2900.0], seed=456)
    st0 = random(ns, range=[0.0, 100.0], seed=789)

    call output_array(sx, './model/sx.bin')
    call output_array(sz, './model/sz.bin')
    call output_array(st0, './model/st0.bin')

    open (3, file='./geometry/geometry.txt')
    do i = 1, size(sx)
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) sx(i), 0.0, sz(i), st0(i)
        write (4, *)
        write (4, *) 50
        do j = 1, 50
            write (4, *) (j - 1)*100.0 + 50.0, 0.0, 10.0, 1.0
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
        write (4, *) sx(i), 0.0, sz(i), 0.0
        write (4, *)
        write (4, *) 50
        do j = 1, 50
            write (4, *) (j - 1)*100.0 + 50.0, 0.0, 10.0, 1.0
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

    sz = zeros(ns) + 2980.0
    call output_array(sz, './model/sz_init_deep.bin')

end program test
