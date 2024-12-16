
program main

    use libflit

    implicit none

    integer :: n1, n2, n3
    real, allocatable, dimension(:, :, :) :: w
    integer, allocatable, dimension(:) :: t1, t2, t3, s1, s2, s3

    if (command_argument_count() == 0) then
        call warn('')
        call warn(' SELECT - SELECT a subset from an array ')
        call warn('')
        call warn(' Usage:')
        call warn('   x_select <in >out [required parameters] [optional parameters]')
        call warn('')
        call warn(' Required Parameters: ')
        call warn('   n1=           number of samples in the 1st dimension ')
        call warn('')
        call warn(' Optional Parameters: ')
        call warn('   n2=           number of samples in the 2nd dimension ')
        call warn('   n3=           number of samples in the 3rd dimension ')
        call warn('   s#=beg,end    begin and end of selection along #-axis, or ')
        call warn('   s#=beg,end,interval    begin, end and interval of selection along #-axis ')
        call warn('   t#=a,b,c,...  selected indices along #-axis ')
        call warn('')
        stop
    end if

    call getpar_int('n1', n1, 0, .true.)
    call getpar_int('n2', n2, floor(get_stdin_size()/4.0d0/n1))
    call getpar_int('n3', n3, floor(get_stdin_size()/4.0d0/n1/n2))

    call getpar_nint('s1', s1, [1, n1])
    call getpar_nint('s2', s2, [1, n2])
    call getpar_nint('s3', s3, [1, n3])
    if (size(s1) == 2) then
        s1 = pad(s1, [0, 1], method=['const', 'const'], const=1)
    end if
    if (size(s2) == 2) then
        s2 = pad(s2, [0, 1], method=['const', 'const'], const=1)
    end if
    if (size(s3) == 2) then
        s3 = pad(s3, [0, 1], method=['const', 'const'], const=1)
    end if
    call getpar_nint('t1', t1, [0])
    call getpar_nint('t2', t2, [0])
    call getpar_nint('t3', t3, [0])

    call alloc_array(w, [1, n1, 1, n2, 1, n3])
    call stdin_array(w)

    if (all(t1 == 0) .and. all(t2 == 0) .and. all(t3 == 0)) then
        call stdout_array(w(s1(1):s1(2):s1(3), s2(1):s2(2):s2(3), s3(1):s3(2):s3(3)))
    else
        call stdout_array(w(t1, t2, t3))
    end if

end program main

