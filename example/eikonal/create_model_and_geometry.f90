
program test

    use libflit

    real, allocatable, dimension(:, :) :: v
    integer :: n1, n2

    n1 = 201
    n2 = 401

    call make_directory('./geometry')
    call make_directory('./model')

    v = zeros(n1, n2)
    do j = 1, n2
        do i = 1, n1
            v(i, j) = 1000 + (i - 1)*10
        end do
    end do
    call output_array(v, './model/vp_init.bin')
    call output_array(v/sqrt(3.0), './model/vs_init.bin')

    do j = 1, n2
        do i = 1, n1
            v(i, j) = v(i, j) &
                + 200*exp(-((i - 50)**2/3.0 + (j - 100)**2)/(2*20.0**2)) &
                - 400*exp(-((i - 50)**2/10.0 + (j - 200)**2)/(2*20.0**2)) &
                + 500*exp(-((i - 75)**2/5.0 + (j - 300)**2)/(2*10.0**2))
        end do
    end do

    call output_array(v, './model/vp.bin')
    call output_array(v/sqrt(3.0), './model/vs.bin')

    v = v*0
    do j = 1, n2
        v(80 + 10*sin(0.015*j), j) = 1.0
        v(160 - 20*exp(-(j - 100.0)**2/200**2), j) = 2.0
    end do
    call output_array(v, './model/refl.bin')

    open (3, file='./geometry/geometry.txt')
    do i = 1, 20
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) (i - 1)*200 + 100, 0, 0, 0
        write (4, *)
        write (4, *) n2
        do j = 1, n2
            write (4, *) (j - 1)*10, 0, 0, 1
        end do
        close (4)

    end do
    close (3)

end program test
