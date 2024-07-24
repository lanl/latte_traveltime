
program test

    use libflit
    use geological_model_2d

    real, allocatable, dimension(:, :) :: v, r, f
    integer :: n1, n2
    type(rgm2) :: p

    n1 = 51
    n2 = 401

    call make_directory('./model')
    call make_directory('./geometry')

    v = zeros(n1, n2)
    do j = 1, n2
        do i = 1, n1
            v(i, j) = 1000 + (i - 1)*15
        end do
    end do
    call output_array(v, './model/vp_init.bin')

    p%n1 = n1
    p%n2 = n2
    p%nf = 3
    p%nl = 10
    p%lwv = 0.5
    p%disp = [10, 20]
    p%secondary_refl_amp = 0.1
    p%yn_facies = .true.
    p%fwidth = 5
    p%seed = 1202
    call p%generate
    r = p%facies
    r = rescale(r, [-500.0, 500.0])
    v = v + r

    f = clip(p%fault, 0.0, 1.0)
    f = f*600
    where (f /= 0)
        v = f
    end where

    call output_array(v, './model/vp.bin')

    !    print *, norm2(v - load('./test_d/iteration_50/model/updated_vp.bin', n1, n2))/norm2(v)
    !    print *, norm2(v - load('./test_dd/iteration_50/model/updated_vp.bin', n1, n2))/norm2(v)

    open (3, file='./geometry/geometry.txt')
    do i = 1, 40
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) (i - 1)*100 + 50, 0, 0, 0
        write (4, *)
        write (4, *) n2
        do j = 1, n2
            write (4, *) (j - 1)*10, 0, 0, 1
        end do
        close (4)

    end do
    close (3)

end program test
