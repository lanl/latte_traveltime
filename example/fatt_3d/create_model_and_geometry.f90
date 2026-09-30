
program test

    use libflit
    use librgm

    implicit none

    type(rgm3_curved) :: p
    real, allocatable, dimension(:, :, :) :: vp, vp_init
    integer :: n1, n2, n3, nsx, nsy, nrx, nry, i, j, k, m, l
    real :: vtop, vbot

    ! Model size along z, y and x; the grid spacing is 10 m
    n1 = 61
    n2 = 201
    n3 = 201

    call make_directory('./model')
    call make_directory('./geometry')

    ! True model: a layered and faulted random model, whose velocity increases with depth
    p%n1 = n1
    p%n2 = n2
    p%n3 = n3
    p%nf = 4
    p%nl = 10
    p%lwv = 0.5
    p%lwh = 0.2
    p%disp = [5.0, 10.0]
    p%fwidth = 2
    p%vmin = 1000.0
    p%vmax = 3000.0
    p%delta_v = 700.0
    p%delta_strike = [20, 40]
    p%refl_height = [0, 50]
    p%refl_smooth = 10.0
    p%seed = 1122
    call p%generate
    vp = p%vp
    call output_array(vp, './model/vp.bin')

    ! Initial model: a linear increase with depth between the mean velocities at the top and bottom
    vtop = mean(vp(1, :, :))
    vbot = mean(vp(n1, :, :))
    vp_init = zeros(n1, n2, n3)
    do i = 1, n1
        vp_init(i, :, :) = vtop + (vbot - vtop)*(i - 1.0)/(n1 - 1.0)
    end do
    call output_array(vp_init, './model/vp_init.bin')

    ! Refraction geometry on the surface: 4 x 4 sources, 500 m apart, and 41 x 41 receivers,
    ! 50 m apart; every source records all receivers
    nsx = 4
    nsy = 4
    nrx = 41
    nry = 41

    open (3, file='./geometry/geometry.txt')
    l = 0
    do j = 1, nsy
        do i = 1, nsx

            l = l + 1
            write (3, *) 'shot_'//num2str(l)//'_geometry.txt'

            open (4, file='./geometry/shot_'//num2str(l)//'_geometry.txt')
            write (4, *) l
            write (4, *)
            write (4, *) 1
            write (4, *) (i - 1)*500.0 + 250.0, (j - 1)*500.0 + 250.0, 0.0, 0.0
            write (4, *)
            write (4, *) nrx*nry
            do m = 1, nry
                do k = 1, nrx
                    write (4, *) (k - 1)*50.0, (m - 1)*50.0, 0.0, 1.0
                end do
            end do
            close (4)

        end do
    end do
    close (3)

end program test
