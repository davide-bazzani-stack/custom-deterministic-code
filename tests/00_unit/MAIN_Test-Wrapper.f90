program Test_wrapper
    use, intrinsic :: iso_fortran_env, only : error_unit, output_unit, int64, real64
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use BiCGSTAB_hybrid_suite, only : collect_suite_BiCGSTAB_hybrid
    use GaussElimination_suite, only : collect_suite_GaussElimination
    use Eigvals_Solver_suite, only : collect_suite_Eigvals_Solver

    implicit none

    integer :: stat, suite_stat, failed_suites, is, detail_unit, io_status
    integer :: dot_count
    integer(int64) :: clock_rate, suite_start, suite_end
    integer(int64) :: total_start, total_end
    real(real64) :: suite_time, total_time
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=128) :: line_prefix
    character(len=1024) :: detail_line
    character(len=8) :: outcome
    character(len=16) :: raw_time
    character(len=24) :: time_text

    stat = 0
    failed_suites = 0

    testsuites = [ &
        new_testsuite("BiCGSTAB", collect_suite_BiCGSTAB_hybrid), &
        new_testsuite("Gaussian elimination", collect_suite_GaussElimination), &
        new_testsuite("Eigenvalue spectrum", collect_suite_Eigvals_Solver) &
        ]

    call system_clock(count_rate=clock_rate)
    if (clock_rate <= 0_int64) then
        write(error_unit, '(A)') '[SYSTEM] Test timer is unavailable.'
        error stop 1
    end if
    call system_clock(total_start)

    write(output_unit, '(/,7X,A)') 'Start 1: Solver suite'
    do is = 1, size(testsuites)
        open(newunit=detail_unit, status='scratch', action='readwrite', &
             iostat=io_status)
        if (io_status /= 0) then
            write(error_unit, '(A)') '[SYSTEM] Unable to open test detail stream.'
            error stop 1
        end if

        suite_stat = 0
        call system_clock(suite_start)
        call run_testsuite(testsuites(is)%collect, detail_unit, suite_stat, &
                           parallel=.false.)
        call system_clock(suite_end)
        suite_time = real(suite_end-suite_start,real64)/real(clock_rate,real64)

        if (suite_stat == 0) then
            outcome = '[PASSED]'
        else
            outcome = '[FAILED]'
            failed_suites = failed_suites + 1
        end if
        stat = stat + suite_stat

        write(line_prefix, '(I0,A,I0,A,A)') is, '/', size(testsuites), &
            '   Testing ', trim(testsuites(is)%name)
        dot_count = max(3, 58-len_trim(line_prefix))
        write(raw_time, '(F8.2)') suite_time
        time_text = '['//trim(adjustl(raw_time))//' sec]'
        write(output_unit, '(A,1X,A,1X,A,1X,A)') trim(line_prefix), &
            repeat('.',dot_count), trim(outcome), trim(time_text)

        if (suite_stat > 0) then
            rewind(detail_unit)
            do
                read(detail_unit, '(A)', iostat=io_status) detail_line
                if (io_status /= 0) exit
                write(error_unit, '(A)') trim(detail_line)
            end do
        end if
        close(detail_unit)
    end do

    call system_clock(total_end)
    total_time = real(total_end-total_start,real64)/real(clock_rate,real64)
    write(raw_time, '(F8.2)') total_time
    write(output_unit, '(/,I0,A,I0,A,A,A)') &
        size(testsuites)-failed_suites, '/', size(testsuites), &
        ' test groups passed. Total test time: ', trim(adjustl(raw_time)), ' sec'

    if (stat > 0) then
        write(error_unit, '(A, I0, 1X, A)') &
            '[SYSTEM] ', stat, "test(s) failed!"
        error stop
    end if

end program Test_wrapper
