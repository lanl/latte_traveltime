!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

program main

    use libflit
    use librgm

    implicit none

    type(rgm2_curved) :: p
    integer :: i, ibeg, iend
    integer :: nt, nv, nm
    real, allocatable, dimension(:) :: height, slope, lwv, height2, lwv2, rmo
    integer, allocatable, dimension(:) :: nf
    real, allocatable, dimension(:, :) :: ppick, rmask
    integer :: i1, i2, l, l1, pick1, pick2
    integer :: dist, s0
    integer, allocatable, dimension(:) :: nmeq
    character(len=1024) :: dir_output

    call getpar_string('outdir', dir_output, './dataset2')
    call getpar_int('ntrain', nt, 1000)
    call getpar_int('nvalid', nv, 100)
    nm = nt + nv

    call getpar_int('ibeg', ibeg, 1)
    call getpar_int('iend', iend, nm)

    call make_directory(tidy(dir_output)//'/data_train')
    call make_directory(tidy(dir_output)//'/data_valid')
    call make_directory(tidy(dir_output)//'/target_train')
    call make_directory(tidy(dir_output)//'/target_valid')

    nf = irandom(nm, range=[1, 12], seed=111)
    height = random(nm,  range=[2.0, 12.0], seed=112)
    height2 = random(nm,  range=[10.0, 20.0], seed=113)
    lwv = zeros(nm)
    lwv2 = random(nm, range=[-0.2, 0.2], seed=115)
    slope = random(nm,  range=[-50.0, 50.0], seed=116)
    nmeq = nint(rescale(nf*1.0, [300.0, 2000.0]))
    rmo = random(nm, range=[0.1, 0.6], seed=117)

    do i = ibeg, iend

        ! Seeds of the random draws for model i; FLIT restarts its generator from
        ! the seed in each call, so every draw needs its own seed
        s0 = 200000*i

        p%n1 = 256
        p%n2 = 256
        p%nf = nf(i)
        p%refl_slope = slope(i)
        p%nl = 20
        p%fwidth = 2.0

        ! Faults with dips that change with depth, and with displacements that die
        ! out within the model; the displacements do not decay away from the faults
        p%delta_dip = [0.0, 15.0]
        p%yn_vary_disp = .false.
        p%yn_disp_decay = .false.

        if (mod(irand(range=[1, nm], seed=s0 + 1), 3) == 0) then
            p%refl_shape = 'gaussian'
            p%refl_mu2 = [0.0, p%n2 - 1.0]
            p%refl_sigma2 = [40.0, 90.0]
            p%ng = irand(range=[2, 6], seed=s0 + 2)
            p%refl_height = [0.25*height2(i), height2(i)]
            p%lwv = lwv2(i)
        else
            p%refl_shape = 'random'
            p%refl_smooth = 30
            p%refl_height = [0.0, height(i)]
            p%lwv = lwv(i)
        end if

        if (mod(irand(range=[1, nm], seed=s0 + 3), 2) == 0) then
            p%unconf = 2
            p%unconf_height = [0.05, 0.1]*p%n1
            p%unconf_z = [0.05, 0.7]
        else
            p%unconf = 0
        end if

        if (nf(i) > 8) then
            p%yn_regular_fault = .true.
            p%nf = nf(i)
            if (mod(irand(range=[1, 10], seed=s0 + 4), 2) == 0) then
                p%dip = [rand(range=[100.0, 120.0], seed=s0 + 5), rand(range=[60.0, 80.0], seed=s0 + 6)]
                p%disp = [3.0, -3.0]
            else
                p%dip = [rand(range=[60.0, 80.0], seed=s0 + 5), rand(range=[100.0, 120.0], seed=s0 + 6)]
                p%disp = [-3.0, 3.0]
            end if
        else
            p%yn_regular_fault = .false.
            p%nf = nf(i)
            p%disp = [5.0, 30.0]
            p%dip = [55.0, 125.0]
        end if

        p%yn_fault = .true.
        p%seed = i*10
        call p%generate

        where (p%fault /= 0)
            p%fault = 1.0
        end where

        ppick = zeros(p%n1, p%n2)
        l1 = 0
        do l = 1, 5*maxval(nmeq)

            if (l1 < nmeq(i)) then

                if (mod(l, 3) == 0) then
                    dist = irand(range=[0, 2], seed=s0 + 100 + 3*l)
                else
                    dist = irand(range=[0, 6], seed=s0 + 100 + 3*l)
                end if

                pick1 = irand(range=[dist + 1, p%n1 - dist], seed=s0 + 101 + 3*l)
                pick2 = irand(range=[dist + 1, p%n2 - dist], seed=s0 + 102 + 3*l)

                if (any(p%fault(pick1 - dist:pick1 + dist, pick2 - dist:pick2 + dist) == 1)) then
                    dist = 4
                    do i2 = -2*dist  - 1, 2*dist + 1
                        do i1 = -2*dist  - 1, 2*dist + 1
                            if (pick1 + i1 >= 1 .and. pick1 + i1 <= p%n1 &
                                    .and. pick2 + i2 >= 1 .and. pick2 + i2 <= p%n2) then
                                ppick(pick1 + i1, pick2 + i2) = &
                                    max(ppick(pick1 + i1, pick2 + i2), exp(-0.3*(i1**2 + i2**2)))
                            end if
                        end do
                    end do
                    l1 = l1 + 1
                end if

            end if

        end do

        rmask = random_mask_smooth(p%n1, p%n2, gs=[4.0, 4.0], mask_out=rmo(i), seed=s0 + 7)

        if (i <= nt) then
            call output_array(ppick, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_meq.bin')
            call output_array(p%fault*rmask, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0*rmask, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fdip.bin')

            call output_array(p%fault, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fdip.bin')

        else
            call output_array(ppick, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_meq.bin')
            call output_array(p%fault*rmask, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0*rmask, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fdip.bin')

            call output_array(p%fault, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fdip.bin')

        end if

        if (i <= nt) then
            print *, date_time_compact(), ' train', i - 1, l1
        else
            print *, date_time_compact(), ' valid', i - nt - 1, l1
        end if

    end do

end program main
