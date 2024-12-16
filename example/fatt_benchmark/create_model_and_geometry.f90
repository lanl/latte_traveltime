
program test

    use libflit
    use librgm

    real, allocatable, dimension(:, :) :: v
    integer :: n1, n2

    n1 = 201
    n2 = 201

    call make_directory('./model')
    call make_directory('./geometry')

    v = zeros(n1, n2) + 2000
    call output_array(v, './model/vp_homo.bin')

    do j = 1, n2
        do i = 1, n1
            if ((i - n1/2.0)**2 + (j - n2/2.0)**2 <= 50.0**2) then
                v(i, j) = 1000.0
            end if
        end do
    end do
    call output_array(v, './model/vp_low.bin')

        do j = 1, n2
        do i = 1, n1
            if ((i - n1/2.0)**2 + (j - n2/2.0)**2 <= 50.0**2) then
                v(i, j) = 3000.0
            end if
        end do
    end do
    call output_array(v, './model/vp_high.bin')

    ! geometry
    open (3, file='./geometry/geometry.txt')
    do i = 1, 1
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open (4, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write (4, *) i
        write (4, *)
        write (4, *) 1
        write (4, *) 250.0, 0.0, 1000.0, 0.0
        write (4, *)
        write (4, *) 5
        do j = 1, 5
            write (4, *) 1750.0, 0, (j - 1)*250.0 + 500.0, 1
        end do
        close (4)

    end do
    close (3)

end program test
