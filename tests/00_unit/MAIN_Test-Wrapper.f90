program Test_wrapper
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use BiCGSTAB_suite, only : collect_suite_BiCGSTAB
    use BiCGSTAB_hybrid_suite, only : collect_suite_BiCGSTAB_hybrid

    implicit none

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("1. Legacy BiCGSTAB solvers", collect_suite_BiCGSTAB), &
        new_testsuite("2. Hybrid BiCGSTAB solvers", collect_suite_BiCGSTAB_hybrid) &
        ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(A, I0, 1X, A)') &
            '[SYSTEM] ', stat, "test(s) failed!"
        error stop
    end if

end program Test_wrapper
