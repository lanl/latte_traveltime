
program test

    use libflit

    implicit none

    integer :: ishot, i

    call make_directory('./geometry')

    open (3, file='./geometry/geometry.txt')
    do i = 1, 49
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'
    end do

    close (3)

    do ishot = 1, 49

        open (3, file='./geometry/shot_'//num2str(ishot)//'_geometry.txt')
        write (3, *) ishot
        write (3, *)
        write (3, *) 1

        write (3, '(3es, es)') 0.0 + (ishot - 1)*2.0, 0.0, 0.0, 0.0
        write (3, *)

        write (3, *) 48
        do i = 1, 48
            write (3, '(3es, es)') 1.0 + (i - 1)*2.0, 0.0, 0.0, 1.0
        end do

        close (3)
    end do

end program test
