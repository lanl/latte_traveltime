!
! © 2024. Triad National Security, LLC. All rights reserved.
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

    type(rgm3) :: p
    integer :: i, ibeg, iend
    integer :: nt = 10, nv = 10
    integer :: nm
    real, allocatable, dimension(:) :: height, height2, slope, lwv, lwv2, rmo
    integer, allocatable, dimension(:) :: nf
    real, allocatable, dimension(:, :, :) :: ppick, rmask
    integer :: i1, i2, i3
    integer :: l1, l
    integer :: pick1, pick2, pick3
    integer :: dist
    integer, allocatable, dimension(:) :: nmeq
    character(len=1024) :: dir_output

    call getpar_string('outdir', dir_output, './dataset3')
    call getpar_int('ntrain', nt, 1000)
    call getpar_int('nvalid', nv, 100)
    nm = nt + nv

    call getpar_int('ibeg', ibeg, 1)
    call getpar_int('iend', iend, nm)

    call make_directory(tidy(dir_output)//'/data_train')
    call make_directory(tidy(dir_output)//'/data_valid')
    call make_directory(tidy(dir_output)//'/target_train')
    call make_directory(tidy(dir_output)//'/target_valid')

    nf = irandom(nm, range=[1, 12])
    height = random(nm,  range=[2.0, 20.0])
    height2 = random(nm,  range=[15.0, 25.0])
    lwv = random(nm,  range=[0.0, 0.2])
    lwv2 = random(nm,  range=[0.3, 0.7])
    slope = random(nm,  range=[-25.0, 25.0])
    p%dip = [50.0, 130.0]
    p%strike = [0.0, 180.0]
    p%rake = [0.0, 180.0]
    nmeq = nint(rescale(nf*1.0, [400.0, 5000.0]))
    rmo = random(nm, range=[0.1, 0.6])

    call getpar_int('ibeg', ibeg, 1)
    call getpar_int('iend', iend, nm)

    do i = ibeg, iend

        p%n1 = 128
        p%n2 = 180
        p%n3 = 180
        p%nf = nf(i)
        p%refl_slope = slope(i)
        p%nl = 20
        p%refl_amp = [0.1, 1.0]
        p%disp = [2.0, 30.0]
        p%fwidth = 2.0

        if (mod(irand(range=[1, nm]), 3) == 0) then
            p%refl_shape = 'gaussian'
            p%refl_mu2 = [0.0, p%n2 - 1.0]
            p%refl_mu3 = [0.0, p%n3 - 1.0]
            p%refl_sigma2 = [25.0, 50.0]
            p%refl_sigma3 = [25.0, 50.0]
            p%ng = irand(range=[2, 4])
            p%refl_height = [0.0, height2(i)]
            p%lwv = lwv2(i)
            p%secondary_refl_amp = 0.05
        else
            p%refl_shape = 'random'
            p%refl_height = [0.0, height(i)]
            p%lwv = lwv(i)
            p%secondary_refl_amp = 0.15
        end if

        if (mod(irand(range=[1, nm]), 3) == 0) then
            p%unconf = irand(range=[1, 2])
            p%unconf_z = [0.05, 0.7]
            p%unconf_amp = [0.05, 0.1]
        else
            p%unconf = 0
        end if

        if (nf(i) > 6) then
            p%yn_regular_fault = .true.
            p%nf = nf(i)
            if (mod(irand(range=[1, 10]), 2) == 0) then
                p%dip = [rand(range=[100.0, 120.0]), rand(range=[60.0, 80.0])]
                p%strike = [rand(range=[0.0, 90.0]), rand(range=[90.0, 180.0])]
                p%rake = [rand(range=[0.0, 90.0]), rand(range=[90.0, 180.0])]
                p%disp = [3.0, -3.0]
            else
                p%dip = [rand(range=[60.0, 80.0]), rand(range=[100.0, 120.0])]
                p%strike = [rand(range=[90.0, 180.0]), rand(range=[0.0, 90.0])]
                p%rake = [rand(range=[90.0, 180.0]), rand(range=[0.0, 90.0])]
                p%disp = [-3.0, 3.0]
            end if
        else
            p%yn_regular_fault = .false.
            p%nf = nf(i)
            p%disp = [2.0, 20.0]
            p%dip = [55.0, 125.0]
            p%strike = [0.0, 180.0]
        end if

        p%f0 = p%f0
        p%yn_fault = .true.
        call p%generate

        where (p%fault /= 0)
            p%fault = 1.0
        end where

        ppick = zeros(p%n1, p%n2, p%n3)
        l1 = 0
        do l = 1, 5*maxval(nmeq)

            if (l1 < nmeq(i)) then

                if (p%yn_regular_fault) then
                    dist = irand(range=[0, 2])
                else
                    dist = irand(range=[0, 5])
                end if

                pick1 = irand(range=[dist + 1, p%n1 - dist])
                pick2 = irand(range=[dist + 1, p%n2 - dist])
                pick3 = irand(range=[dist + 1, p%n3 - dist])

                if (any(p%fault(pick1 - dist:pick1 + dist, pick2 - dist:pick2 + dist, pick3 - dist:pick3 + dist) == 1)) then
                    dist = 4
                    do i3 = -2*dist - 1, 2*dist + 1
                        do i2 = -2*dist - 1, 2*dist + 1
                            do i1 = -2*dist - 1, 2*dist + 1
                                if (pick1 + i1 >= 1 .and. pick1 + i1 <= p%n1 &
                                        .and. pick2 + i2 >= 1 .and. pick2 + i2 <= p%n2 &
                                        .and. pick3 + i3 >= 1 .and. pick3 + i3 <= p%n3) then
                                    ppick(pick1 + i1, pick2 + i2, pick3 + i3) = &
                                        max(ppick(pick1 + i1, pick2 + i2, pick3 + i3), exp(-0.3*(i1**2 + i2**2 + i3**2)))
                                end if
                            end do
                        end do
                    end do
                    l1 = l1 + 1
                end if

            end if

        end do

        rmask = random_mask_smooth(p%n1, p%n2, p%n3, gs=[4.0, 4.0, 4.0], mask_out=rmo(i))

        if (i <= nt) then
            call output_array(ppick, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_meq.bin')
            call output_array(p%fault*rmask, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0*rmask, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fdip.bin')
            call output_array(p%fault_strike/180.0*rmask, tidy(dir_output)//'/data_train/'//num2str(i - 1)//'_fstrike.bin')

            call output_array(p%fault, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fdip.bin')
            call output_array(p%fault_strike/180.0, tidy(dir_output)//'/target_train/'//num2str(i - 1)//'_fstrike.bin')
        else
            call output_array(ppick, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_meq.bin')
            call output_array(p%fault*rmask, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0*rmask, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fdip.bin')
            call output_array(p%fault_strike/180.0*rmask, tidy(dir_output)//'/data_valid/'//num2str(i - nt - 1)//'_fstrike.bin')

            call output_array(p%fault, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fsem.bin')
            call output_array(p%fault_dip/180.0, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fdip.bin')
            call output_array(p%fault_strike/180.0, tidy(dir_output)//'/target_valid/'//num2str(i - nt - 1)//'_fstrike.bin')
        end if

        if (i <= nt) then
            print *, date_time_compact(), ' train', i - 1, l1
        else
            print *, date_time_compact(), ' valid', i - nt - 1, l1
        end if

    end do

end program main
