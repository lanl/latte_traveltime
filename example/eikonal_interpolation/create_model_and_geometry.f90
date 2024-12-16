
program test

    use libflit
    use librgm

    implicit none

    block

        integer :: l, i, j, k
        integer :: n1, n2
        real, dimension(:), allocatable :: rx, rz, r, t, topo
        real, allocatable, dimension(:, :) :: rr, vp, tref
        real :: sx, sz
        real :: v0
        type(fractal_noise_1d) :: p
        integer :: iz, nrx, nrz
        real :: dx, dz
        real :: g(1:2)
        real :: s0

        call make_directory('./geometry')
        call make_directory('./model')

        ! ============================================================
        ! Test with random source receiver locations

        n1 = 101
        n2 = 101
        v0 = 1500

        rx = random(400, range=[0.0, 1000.0], seed=121)
        rz = random(400, range=[0.0, 1000.0], seed=333)

        sz = 67.77
        sx = 453.91

        !    sz = 65
        !    sx = 455

        open (3, file='./geometry/geometry.txt')
        write (3, *) 'shot_1_geometry.txt'
        close (3)

        open (4, file='./geometry/shot_1_geometry.txt')
        write (4, *) 1
        write (4, *)
        write (4, *) 1
        write (4, *) sx, 0.0, sz, 0.0
        write (4, *)

        r = linspace(0.0, 1000.0, 600)
        rr = zeros(size(r) + size(rx), 2)
        rr(1:size(r), 1) = r
        rr(1:size(r), 2) = sz
        rr(size(r) + 1:, 1) = rx
        rr(size(r) + 1:, 2) = rz

        write (4, *) size(rr, 1)

        do i = 1, size(rr, 1)
            write (4, *) rr(i, 1), 0.0, rr(i, 2), 1.0
        end do
        close (4)

        ! Reference solution on receivers
        t = zeros(size(rr, 1))
        do i = 1, size(rr, 1)
            t(i) = norm2(rr(i, :) - [sx, sz])/v0
        end do
        call output_array(t, './tref_homo.bin')

        ! Output receiver locations for plotting
        open (3, file='./recr_homo.txt')
        do i = 1, size(rr, 1)
            write (3, *) rr(i, 1), rr(i, 2), 1.0
        end do
        close (3)
        open (3, file='./srcr_homo.txt')
        write (3, *) sx, sz, 1.0
        close (3)

        ! Reference snapshot
        rr = zeros(101, 101)
        do j = 1, n2
            do i = 1, n1
                rr(i, j) = norm2([(j - 1)*10.0, (i - 1)*10.0] - [sx, sz])/v0
            end do
        end do
        call output_array(rr, './sref_homo.bin')

        ! Create velocity model
        rr = ones(101, 101)*v0
        call output_array(rr, './model/v_homo.bin')

        ! ============================================================
        ! Another test for topography

        n1 = 31
        n2 = 51
        dx = 10.0
        dz = 10.0
        v0 = 2500

        vp = ones(n1, n2)*v0
        p%n1 = n2
        p%octaves = 5
        p%seed = 12123
        topo = p%generate()
        topo = rescale(topo, [0.0, 7.0])

        do i = 1, n2
            vp(1:nint(topo(i)) + 1, i) = 0.0
        end do
        call output_array(vp, './model/v_topo.bin')

        sx = 257.2
        sz = 253.1

        open (3, file='./geometry/geometry_topo.txt')
        write (3, *) 'shot_1_geometry_topo.txt'
        close (3)

        open (33, file='./topo.txt')

        open (4, file='./geometry/shot_1_geometry_topo.txt')
        write (4, *) 1
        write (4, *)
        write (4, *) 1
        write (4, *) sx, 0.0, sz, 0.0
        write (4, *)

        nrx = 5
        nrz = 5

        rr = zeros((n2 - 1)*nrz*nrx, 2)

        write (4, *) size(rr, 1)

        l = 1
        do i = 1, n2 - 1

            iz = index_first_nonzero(vp(:, i)) - 1

            do j = 1, nrx
                do k = 1, nrz
                    rr(l, 1) = (i - 1 - 0.5)*dx + dx/(nrx - 1.0)*(j - 1)
                    rr(l, 2) = (iz - 1 + 0.5)*dz + dz/(nrz - 1.0)*(k - 1)
                    l = l + 1
                end do
            end do

            write (33, *) (i - 1)*dx, iz*dz

        end do

        iz = index_first_nonzero(vp(:, n2)) - 1
        write (33, *) (n2 - 1)*dx, iz*dz
        close (33)

        do i = 1, size(rr, 1)
            write (4, *) rr(i, 1), 0.0, rr(i, 2), 1.0
        end do

        close (4)

        print *, size(rr, 1)

        ! Reference solution on receivers
        t = zeros(size(rr, 1))
        do i = 1, size(rr, 1)
            t(i) = norm2(rr(i, :) - [sx, sz])/v0
        end do

        open (3, file='./tref_topo.txt')
        do i = 1, size(rr, 1)
            write (3, *) rr(i, :), t(i)
        end do
        close (3)

        ! Reference snapshot
        rr = zeros(n1, n2)
        do j = 1, n2
            do i = 1, n1
                rr(i, j) = norm2([(j - 1)*10.0, (i - 1)*10.0] - [sx, sz])/v0
            end do
        end do
        call output_array(rr, './sref_topo.bin')

        ! ============================================================
        ! Test with random source receiver locations

        n1 = 101
        n2 = 101
        g = [0.2, 0.7]
        s0 = 1.0/3000.0

        sz = 505.0
        sx = 505.0

        open (3, file='./geometry/geometry_gradient.txt')
        write (3, *) 'shot_1_geometry_gradient.txt'
        close (3)

        open (4, file='./geometry/shot_1_geometry_gradient.txt')
        write (4, *) 1
        write (4, *)
        write (4, *) 1
        write (4, *) sx, 0.0, sz, 0.0
        write (4, *)

        rx = random(500, range=[0.0, 1000.0], seed=444)
        rz = random(500, range=[0.0, 1000.0], seed=555)

        write (4, *) size(rx)

        do i = 1, size(rx, 1)
            write (4, *) rx(i), 0.0, rz(i), 1.0
        end do
        close (4)

        rr = zeros(n1, n2)
        tref = zeros(n1, n2)
        do j = 1, n2
            do i = 1, n1
                rr(i, j) = 1.0/(1.0/s0 + dot_product(g, [(i - 1)*10.0 - sz, (j - 1)*10.0 - sx]))
                tref(i, j) = 1.0/norm2(g)*acosh(1.0 + 0.5*rr(i, j)*s0*norm2(g)**2*norm2([(i - 1)*10.0 - sz, (j - 1)*10.0 - sx])**2)
            end do
        end do
        call output_array(1.0/rr, './model/v_gradient.bin')

        ! Reference solution on receivers
        t = zeros(size(rx))
        do i = 1, size(rx)
            v0 = 1.0/(1.0/s0 + dot_product(g, [rz(i) - sz, rx(i) - sx]))
            t(i) = 1.0/norm2(g)*acosh(1.0 + 0.5*v0*s0*norm2(g)**2*norm2([rz(i) - sz, rx(i) - sx])**2)
        end do
        call output_array(t, './tref_gradient.bin')

        ! Output receiver locations for plotting
        open (3, file='./recr_gradient.txt')
        do i = 1, size(rx)
            write (3, *) rx(i), rz(i), 1.0
        end do
        close (3)
        open (3, file='./srcr_gradient.txt')
        write (3, *) sx, sz, 1.0
        close (3)

        ! Reference snapshot
        call output_array(tref, './sref_gradient.bin')

    end block

    block

        integer :: i, j, k
        integer :: n1, n2, n3
        real, dimension(:), allocatable :: rx, ry, rz, t
        real, allocatable, dimension(:, :, :) :: rr, tref
        real :: sx, sy, sz
        real :: v0
        real :: g(1:3)
        real :: s0

        n1 = 101
        n2 = 121
        n3 = 141
        g = [0.2, 0.5, 0.7]
        s0 = 1.0/3000.0

        sz = 505.0
        sy = 505.0
        sx = 505.0

        open (3, file='./geometry/geometry3_gradient.txt')
        write (3, *) 'shot_1_geometry3_gradient.txt'
        close (3)

        open (4, file='./geometry/shot_1_geometry3_gradient.txt')
        write (4, *) 1
        write (4, *)
        write (4, *) 1
        write (4, *) sx, sy, sz, 0.0
        write (4, *)

        rx = random(500, range=[0.0, 1000.0], seed=333)
        ry = random(500, range=[0.0, 1000.0], seed=444)
        rz = random(500, range=[0.0, 1000.0], seed=555)

        write (4, *) size(rx)

        do i = 1, size(rx, 1)
            write (4, *) rx(i), ry(i), rz(i), 1.0
        end do
        close (4)

        rr = zeros(n1, n2, n3)
        tref = zeros(n1, n2, n3)
        !$omp parallel do private(i, j, k)
        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    rr(i, j, k) = 1.0/(1.0/s0 + dot_product(g, [(i - 1)*10.0 - sz, (j - 1)*10.0 - sy, (k - 1)*10.0 - sx]))
                    tref(i, j, k) = 1.0/norm2(g)*acosh(1.0 + 0.5*rr(i, j, k)*s0*norm2(g)**2 &
                        *norm2([(i - 1)*10.0 - sz, (j - 1)*10.0 - sy, (k - 1)*10.0 - sx])**2)
                end do
            end do
        end do
        !$omp end parallel do
        call output_array(1.0/rr, './model/v3_gradient.bin')

        ! Reference solution on receivers
        t = zeros(size(rx))
        do i = 1, size(rx)
            v0 = 1.0/(1.0/s0 + dot_product(g, [rz(i) - sz, ry(i) - sy, rx(i) - sx]))
            t(i) = 1.0/norm2(g)*acosh(1.0 + 0.5*v0*s0*norm2(g)**2*norm2([rz(i) - sz, ry(i) - sy, rx(i) - sx])**2)
        end do
        call output_array(t, './tref3_gradient.bin')

        ! Output receiver locations for plotting
        open (3, file='./recr3_gradient.txt')
        do i = 1, size(rx)
            write (3, *) rz(i), ry(i), rx(i), 1.0
        end do
        close (3)
        open (3, file='./srcr3_gradient.txt')
        write (3, *) sz, sy, sx, 1.0
        close (3)

        ! Reference snapshot
        call output_array(tref, './sref3_gradient.bin')

    end block

end program test
