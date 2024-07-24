
program test

    use libflit

    integer :: ns, i
    real, allocatable, dimension(:) :: st0, noise
    real, allocatable, dimension(:) :: v, w
    integer :: nrx, nrz

    ns = 1200
    nrx = 50
    nrz = 30

    ! Add noise
    call make_directory('./data_noisy')
    do i = 1, ns
        st0 = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', nrx + nrz)
        noise = random(size(st0), seed=i*i)
        noise = gauss_filt(noise, 2.0)
        noise = noise - mean(noise)
        noise = noise/maxval(noise)*0.02
        st0 = st0 + noise
        call output_array(st0, './data_noisy/shot_'//num2str(i)//'_traveltime_p.bin')
        print *, i
    end do

    ! For plotting, output clean data
    call make_directory('./data_no_st0')
    st0 = random(ns, range=[0.0, 10.0], seed=789)
    do i = 1, ns
        v = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', nrx + nrz) - st0(i)
        w = load('./data_noisy/shot_'//num2str(i)//'_traveltime_p.bin', nrx + nrz) - st0(i)
        call output_array(v, './data_no_st0/shot_'//num2str(i)//'_traveltime_p.bin.clean')
        call output_array(w, './data_no_st0/shot_'//num2str(i)//'_traveltime_p.bin.noisy')
        call output_array(w - v, './data_no_st0/shot_'//num2str(i)//'_traveltime_p.bin.noise')
        print *, i
    end do

end program test
