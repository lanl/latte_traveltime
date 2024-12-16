
program test

    use libflit

    real, allocatable, dimension(:) :: data, noise
    integer, allocatable, dimension(:, :) :: shot_in_rank

    integer :: i

    call mpistart

    call make_directory('./data_noisy')

    call alloc_array(shot_in_rank, [0, nrank - 1, 1, 2])
    call cut(1, 1200, nrank, shot_in_rank)

    do i = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

        data = load('./data/shot_'//num2str(i)//'_traveltime_p.bin', 150)
        noise = gauss_filt(random(150, seed=i + 1000, dist='normal'), 2.0)
        noise = noise - mean(noise)
        noise = noise/maxval(abs(noise))*0.1
        data = data + noise
        call output_array(data, './data_noisy/shot_'//num2str(i)//'_traveltime_p.bin')

        data = load('./data/shot_'//num2str(i)//'_traveltime_s.bin', 150)
        noise = gauss_filt(random(150, seed=i + 11000, dist='normal'), 2.0)
        noise = noise - mean(noise)
        noise = noise/maxval(abs(noise))*0.1
        data = data + noise
        call output_array(data, './data_noisy/shot_'//num2str(i)//'_traveltime_s.bin')

        print *, i

    end do

    call mpibarrier
    call mpiend

end program test
