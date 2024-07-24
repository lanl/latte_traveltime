
program test

    use libflit

    real, allocatable, dimension(:) :: st0, w, v, data, noise

    integer :: i

    ! For plotting, output clean data
    call make_directory('./data_no_st0')
    st0 = load('./model/st0.bin', 1200)
    do i = 1, 1200

        v = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', 150) - st0(i)
        w = load('./data_noisy/shot_'//num2str(i)//'_traveltime_p.bin', 150) - st0(i)
        call output_array(v, './data_no_st0/shot_'//num2str(i)//'_traveltime_p.bin.clean')
        call output_array(w, './data_no_st0/shot_'//num2str(i)//'_traveltime_p.bin.noisy')
        call output_array(w - v, './data_no_st0/shot_'//num2str(i)//'_traveltime_p.bin.noise')

        v = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', 150) - st0(i)
        w = load('./data_noisy/shot_'//num2str(i)//'_traveltime_s.bin', 150) - st0(i)
        call output_array(v, './data_no_st0/shot_'//num2str(i)//'_traveltime_s.bin.clean')
        call output_array(w, './data_no_st0/shot_'//num2str(i)//'_traveltime_s.bin.noisy')
        call output_array(w - v, './data_no_st0/shot_'//num2str(i)//'_traveltime_s.bin.noise')

        print *, i

    end do

    stop

    call make_directory('./data_noisy')

    do i = 1, 1200

        data = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', 150)
        noise = gauss_filt(random(150, seed=i + 1000), 4.0)
        noise = noise - mean(noise)
        noise = noise/maxval(abs(noise))*0.04
        data = data + noise
        call output_array(data, './data_noisy/shot_'//num2str(i)//'_traveltime_p.bin')

        data = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', 150)
        noise = gauss_filt(random(150, seed=i + 11000), 4.0)
        noise = noise - mean(noise)
        noise = noise/maxval(abs(noise))*0.05
        data = data + noise
        call output_array(data, './data_noisy/shot_'//num2str(i)//'_traveltime_s.bin')

        print *, i

    end do

end program test
