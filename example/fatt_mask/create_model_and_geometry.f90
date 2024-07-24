
program test

    use libflit
    use geological_model_2d

    real, allocatable, dimension(:, :) :: vp, r
    integer :: n1, n2, ns, nr, i, j
    type(rgm2) :: p

    n1 = 201
    n2 = 801

    call make_directory('./model')
    call make_directory('./geometry')

    p%n1 = n1
    p%n2 = n2
    p%nf = 17
    p%yn_regular_fault = .false.
    p%nl = 40
    p%unconf = 1
    p%unconf_z = 0.3
    p%disp = [20, 30]
    p%lwv = 0.4
    p%dip = [rand(range=[100.0, 120.0], seed=123), rand(range=[60.0, 80.0], seed=567)]
    p%secondary_refl_amp = 0.1
    p%refl_smooth = 10.0
    p%yn_facies = .true.
    p%yn_rgt = .true.
    p%fwidth = 3
    p%seed = 101010
    call p%generate

    vp = rescale(p%facies, [1000.0, 3500.0])
    r = vp
    vp(1:50, :) = 1500
    call output_array(vp, './model/vp.bin')

    vp = zeros_like(r)
    do i = 1, n1
        vp(i, :) = mean(r(i, :))
    end do
    vp = 1.0/gauss_filt(1.0/vp, [3.0, 1.0])
    vp(1:50, :) = 1500
    call output_array(vp, './model/vp_init.bin')

    r = ones_like(vp)
    r(1:50, :) = 0.0
    call output_array(r, './model/mask.bin')

    r = zeros_like(vp)
    r(190, :) = 1.0
    call output_array(r, './model/refl.bin')

    ns = 160
    nr = 160

    open (33, file='./model/rxz.txt')
    open (3, file='./geometry/geometry.txt')
    do i = 1, ns
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) (i - 1)*50.0 + 10.0, 0.0, 0.0, 0.0
        write (4, *)
        write (4, *) nr
        do j = 1, nr
            write (4, *) (j - 1)*50.0 + 40.0, 0.0, 500.0, 1.0
            if(i == 1)  then
                write (33, *) 500.0, (j - 1)*50.0 + 40.0
            end if
        end do
        close (4)

    end do
    close (3)
    close(33)

end program test
