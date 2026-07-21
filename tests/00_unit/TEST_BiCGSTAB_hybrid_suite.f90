module BiCGSTAB_hybrid_suite
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use precision_kinds, only : prec
    use BiCGSTAB_Solver, only : banded_operator_t, bicgstab_options_t, &
        bicgstab_result_t, bicgstab_workspace_t, jacobi_preconditioner_t, &
        symmetric_gauss_seidel_preconditioner_t, ilu0_preconditioner_t, &
        user_preconditioner_t, solve_bicgstab, BICGSTAB_SUCCESS, &
        BICGSTAB_COMPONENT_FAILURE, BICGSTAB_PRECONDITIONER_FAILURE

    implicit none
    private

    type :: custom_jacobi_context_t
        real(prec), allocatable :: inverse_diagonal(:)
    end type custom_jacobi_context_t

    public :: collect_suite_BiCGSTAB_hybrid

contains

    subroutine collect_suite_BiCGSTAB_hybrid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("hybrid_classic_legacy_parity", &
                hybrid_classic_legacy_parity), &
            new_unittest("hybrid_triiso_legacy_parity", &
                hybrid_triiso_legacy_parity), &
            new_unittest("hybrid_preconditioner_variants", &
                hybrid_preconditioner_variants), &
            new_unittest("hybrid_preconditioner_actions", &
                hybrid_preconditioner_actions), &
            new_unittest("hybrid_tridiagonal_operator", &
                hybrid_tridiagonal_operator), &
            new_unittest("hybrid_heptadiagonal_operator", &
                hybrid_heptadiagonal_operator), &
            new_unittest("hybrid_eqtri_operator", &
                hybrid_eqtri_operator), &
            new_unittest("hybrid_initial_solution_and_failures", &
                hybrid_initial_solution_and_failures) &
            ]
    end subroutine collect_suite_BiCGSTAB_hybrid

    subroutine hybrid_classic_legacy_parity(error)
        use BiCGSTAB_Algorithm, only : BiCGSTAB_Classic

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 3
        integer, parameter :: n = 9
        real(prec), target :: a(n), b(n), c(n), d(n), e(n)
        real(prec) :: matrix(n,n), rhs(n), reference(n)
        real(prec) :: legacy_solution(n), hybrid_solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        integer :: legacy_info, status
        real(prec) :: tolerance

        call build_classic_fixture(a, b, c, d, e, matrix, reference)
        rhs = matmul(matrix, reference)
        legacy_solution = 0.0_prec
        hybrid_solution = 0.0_prec
        tolerance = solver_tolerance()

        call BiCGSTAB_Classic(n_x, n, a, b, c, d, e, rhs, &
            legacy_solution, tolerance, 200, legacy_info)
        call check(error, legacy_info, 0, &
            message="[SYSTEM] TEST ERROR: legacy Classic solver failed")
        if (allocated(error)) return

        call operator%bind_cartesian_pentadiagonal(1, n_x, &
            a, b, c, d, e, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: hybrid Cartesian operator binding failed")
        if (allocated(error)) return

        call configure_options(options, 200)
        call solve_bicgstab(operator, rhs, hybrid_solution, options, &
            workspace, result)
        call assert_result_success(result, "hybrid Classic solver failed", error)
        if (allocated(error)) return

        call assert_solution(hybrid_solution, legacy_solution, check_tolerance(), &
            "hybrid Classic result differs from legacy", error)
        if (allocated(error)) return
        call assert_solution(hybrid_solution, reference, check_tolerance(), &
            "hybrid Classic solution mismatch", error)
        if (allocated(error)) return
        call assert_dense_residual(matrix, hybrid_solution, rhs, check_tolerance(), &
            "hybrid Classic residual too large", error)
    end subroutine hybrid_classic_legacy_parity

    subroutine hybrid_triiso_legacy_parity(error)
        use BiCGSTAB_Algorithm, only : BiCGSTAB_TriIso

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 4
        integer, parameter :: n = 8
        real(prec), target :: a(n), b(n), c(n), d(n), e(n)
        real(prec) :: matrix(n,n), rhs(n), reference(n)
        real(prec) :: legacy_solution(n), hybrid_solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        integer :: legacy_info, status

        call build_triiso_fixture(a, b, c, d, e, matrix, reference)
        rhs = matmul(matrix, reference)
        legacy_solution = 0.0_prec
        hybrid_solution = 0.0_prec

        call BiCGSTAB_TriIso(n_x, n, a, b, c, d, e, rhs, &
            legacy_solution, solver_tolerance(), 200, legacy_info)
        call check(error, legacy_info, 0, &
            message="[SYSTEM] TEST ERROR: legacy TriIso solver failed")
        if (allocated(error)) return

        call operator%bind_triiso(n_x, a, b, c, d, e, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: hybrid TriIso operator binding failed")
        if (allocated(error)) return

        call configure_options(options, 200)
        call solve_bicgstab(operator, rhs, hybrid_solution, options, &
            workspace, result)
        call assert_result_success(result, "hybrid TriIso solver failed", error)
        if (allocated(error)) return

        call assert_solution(hybrid_solution, legacy_solution, check_tolerance(), &
            "hybrid TriIso result differs from legacy", error)
        if (allocated(error)) return
        call assert_dense_residual(matrix, hybrid_solution, rhs, check_tolerance(), &
            "hybrid TriIso residual too large", error)
    end subroutine hybrid_triiso_legacy_parity

    subroutine hybrid_preconditioner_variants(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 3
        integer, parameter :: n = 9
        real(prec), target :: a(n), b(n), c(n), d(n), e(n)
        real(prec) :: matrix(n,n), rhs(n), reference(n), solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        type(jacobi_preconditioner_t) :: jacobi
        type(symmetric_gauss_seidel_preconditioner_t) :: sgs
        type(ilu0_preconditioner_t) :: ilu0
        type(user_preconditioner_t) :: user_preconditioner
        type(custom_jacobi_context_t), target :: user_context
        integer :: status

        call build_classic_fixture(a, b, c, d, e, matrix, reference)
        rhs = matmul(matrix, reference)
        call operator%bind_cartesian_pentadiagonal(1, n_x, &
            a, b, c, d, e, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: preconditioner fixture binding failed")
        if (allocated(error)) return
        call configure_options(options, 200)

        call jacobi%setup(operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: Jacobi setup failed")
        if (allocated(error)) return
        solution = 0.0_prec
        call solve_bicgstab(operator, rhs, solution, options, workspace, result, jacobi)
        call assert_preconditioned_solution("Jacobi", result, matrix, solution, &
            rhs, reference, error)
        if (allocated(error)) return

        call sgs%setup(operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: SGS setup failed")
        if (allocated(error)) return
        solution = 0.0_prec
        call solve_bicgstab(operator, rhs, solution, options, workspace, result, sgs)
        call assert_preconditioned_solution("SGS", result, matrix, solution, &
            rhs, reference, error)
        if (allocated(error)) return

        call ilu0%setup(operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: ILU(0) setup failed")
        if (allocated(error)) return
        solution = 0.0_prec
        call solve_bicgstab(operator, rhs, solution, options, workspace, result, ilu0)
        call assert_preconditioned_solution("ILU(0)", result, matrix, solution, &
            rhs, reference, error)
        if (allocated(error)) return

        allocate(user_context%inverse_diagonal(n))
        user_context%inverse_diagonal = 1.0_prec/b
        call user_preconditioner%bind(custom_jacobi_apply, user_context)
        solution = 0.0_prec
        call solve_bicgstab(operator, rhs, solution, options, workspace, result, &
            user_preconditioner)
        call assert_preconditioned_solution("user callback", result, matrix, &
            solution, rhs, reference, error)
    end subroutine hybrid_preconditioner_variants

    subroutine hybrid_preconditioner_actions(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 6
        real(prec), target :: lower(n), diagonal(n), upper(n)
        real(prec) :: rhs(n), reference(n), output(n), expected(n), value
        type(banded_operator_t) :: operator
        type(jacobi_preconditioner_t) :: jacobi
        type(symmetric_gauss_seidel_preconditioner_t) :: sgs
        type(ilu0_preconditioner_t) :: ilu0
        integer :: row, status

        lower = [0.0_prec, -0.8_prec, -0.6_prec, -0.9_prec, -0.7_prec, -0.5_prec]
        diagonal = [4.0_prec, 4.5_prec, 5.0_prec, 4.2_prec, 4.8_prec, 5.2_prec]
        upper = [-0.4_prec, -0.7_prec, -0.5_prec, -0.8_prec, -0.6_prec, 0.0_prec]
        reference = [1.0_prec, -0.5_prec, 2.0_prec, -1.5_prec, 0.25_prec, 1.25_prec]
        rhs = diagonal*reference
        rhs(2:n) = rhs(2:n) + lower(2:n)*reference(1:n-1)
        rhs(1:n-1) = rhs(1:n-1) + upper(1:n-1)*reference(2:n)

        call operator%bind_tridiagonal(lower, diagonal, upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: preconditioner-action operator binding failed")
        if (allocated(error)) return

        call jacobi%setup(operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: Jacobi action setup failed")
        if (allocated(error)) return
        call jacobi%apply(rhs, output, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: Jacobi action failed")
        if (allocated(error)) return
        expected = rhs/diagonal
        call assert_solution(output, expected, check_tolerance(), &
            "Jacobi action mismatch", error)
        if (allocated(error)) return

        call sgs%setup(operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: SGS action setup failed")
        if (allocated(error)) return
        call sgs%apply(rhs, output, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: SGS action failed")
        if (allocated(error)) return
        expected(1) = rhs(1)/diagonal(1)
        do row = 2, n
            expected(row) = (rhs(row) - lower(row)*expected(row-1))/diagonal(row)
        end do
        expected = diagonal*expected
        expected(n) = expected(n)/diagonal(n)
        do row = n - 1, 1, -1
            value = expected(row) - upper(row)*expected(row+1)
            expected(row) = value/diagonal(row)
        end do
        call assert_solution(output, expected, check_tolerance(), &
            "SGS action mismatch", error)
        if (allocated(error)) return

        ! ILU(0) is an exact LU factorization for a tridiagonal matrix.
        call ilu0%setup(operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: ILU(0) action setup failed")
        if (allocated(error)) return
        call ilu0%apply(rhs, output, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: ILU(0) action failed")
        if (allocated(error)) return
        call assert_solution(output, reference, check_tolerance(), &
            "ILU(0) tridiagonal action mismatch", error)
    end subroutine hybrid_preconditioner_actions

    subroutine hybrid_tridiagonal_operator(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 7
        real(prec), target :: lower(n), diagonal(n), upper(n)
        real(prec) :: matrix(n,n), rhs(n), reference(n), solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        integer :: row, status

        lower = [0.0_prec, -0.8_prec, -0.7_prec, -0.9_prec, &
                 -0.6_prec, -0.75_prec, -0.85_prec]
        diagonal = [4.0_prec, 4.2_prec, 4.4_prec, 4.1_prec, &
                    4.3_prec, 4.5_prec, 4.6_prec]
        upper = [-0.5_prec, -0.6_prec, -0.55_prec, -0.7_prec, &
                 -0.65_prec, -0.8_prec, 0.0_prec]
        reference = [1.0_prec, -0.5_prec, 2.0_prec, -1.0_prec, &
                     0.75_prec, 1.5_prec, -2.0_prec]

        matrix = 0.0_prec
        do row = 1, n
            matrix(row,row) = diagonal(row)
        end do
        do row = 2, n
            matrix(row,row - 1) = lower(row)
        end do
        do row = 1, n - 1
            matrix(row,row + 1) = upper(row)
        end do
        rhs = matmul(matrix, reference)
        solution = 0.0_prec

        call operator%bind_cartesian_tridiagonal(lower, diagonal, upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: tridiagonal binding failed")
        if (allocated(error)) return
        call configure_options(options, 100)
        call solve_bicgstab(operator, rhs, solution, options, workspace, result)
        call assert_preconditioned_solution("tridiagonal", result, matrix, &
            solution, rhs, reference, error)
    end subroutine hybrid_tridiagonal_operator

    subroutine hybrid_heptadiagonal_operator(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 12
        integer, parameter :: offset_1 = 1, offset_2 = 3, offset_3 = 6
        real(prec), target :: lower_1(n), diagonal(n), upper_1(n)
        real(prec), target :: lower_2(n), upper_2(n), lower_3(n), upper_3(n)
        real(prec) :: matrix(n,n), rhs(n), reference(n), solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        integer :: row, status

        lower_1 = -0.30_prec
        upper_1 = -0.20_prec
        lower_2 = -0.25_prec
        upper_2 = -0.15_prec
        lower_3 = -0.10_prec
        upper_3 = -0.12_prec
        diagonal = 5.0_prec
        reference = [(0.25_prec*real(row,prec) - 1.0_prec, row=1,n)]

        matrix = 0.0_prec
        do row = 1, n
            matrix(row,row) = diagonal(row)
        end do
        do row = offset_1 + 1, n
            matrix(row,row-offset_1) = lower_1(row)
        end do
        do row = 1, n-offset_1
            matrix(row,row+offset_1) = upper_1(row)
        end do
        do row = offset_2 + 1, n
            matrix(row,row-offset_2) = lower_2(row)
        end do
        do row = 1, n-offset_2
            matrix(row,row+offset_2) = upper_2(row)
        end do
        do row = offset_3 + 1, n
            matrix(row,row-offset_3) = lower_3(row)
        end do
        do row = 1, n-offset_3
            matrix(row,row+offset_3) = upper_3(row)
        end do
        rhs = matmul(matrix, reference)
        solution = 0.0_prec

        call operator%bind_cartesian_heptadiagonal( &
            offset_1, offset_2, offset_3, lower_1, diagonal, upper_1, &
            lower_2, upper_2, lower_3, upper_3, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: heptadiagonal binding failed")
        if (allocated(error)) return
        call configure_options(options, 200)
        call solve_bicgstab(operator, rhs, solution, options, workspace, result)
        call assert_preconditioned_solution("heptadiagonal", result, matrix, &
            solution, rhs, reference, error)
    end subroutine hybrid_heptadiagonal_operator

    subroutine hybrid_eqtri_operator(error)
        use BiCGSTAB_Algorithm, only : BiCGSTAB_EqTri

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 3
        integer, parameter :: n_y = 4
        integer, parameter :: n = n_x*n_y
        real(prec), target :: a(n), b(n), c(n), d(n), e(n)
        real(prec) :: matrix(n,n), rhs(n), reference(n), solution(n)
        real(prec) :: legacy_solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        integer :: legacy_info, row, status

        a = -0.35_prec
        b = 4.5_prec
        c = -0.25_prec
        d = -0.15_prec
        e = -0.20_prec
        reference = [(1.0_prec - 0.1_prec*real(row,prec), row=1,n)]

        matrix = 0.0_prec
        do row = 1, n
            matrix(row,row) = b(row)
        end do
        do row = n_x + 1, n
            matrix(row,row-n_x) = a(row)
        end do
        do row = 1, n-n_x
            matrix(row,row+n_x) = c(row)
        end do
        do row = n_x + 2, n
            matrix(row,row-n_x-1) = d(row)
        end do
        do row = 1, n-n_x-1
            matrix(row,row+n_x+1) = e(row)
        end do
        rhs = matmul(matrix, reference)
        solution = 0.0_prec
        legacy_solution = 0.0_prec

        call BiCGSTAB_EqTri(n_x, n_y, a, b, c, d, e, rhs, &
            legacy_solution, solver_tolerance(), 200, legacy_info)
        call check(error, legacy_info, 0, &
            message="[SYSTEM] TEST ERROR: legacy EqTri solver failed")
        if (allocated(error)) return

        call operator%bind_eqtri(n_x, a, b, c, d, e, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: EqTri binding failed")
        if (allocated(error)) return
        call configure_options(options, 200)
        call solve_bicgstab(operator, rhs, solution, options, workspace, result)
        call assert_preconditioned_solution("EqTri", result, matrix, solution, &
            rhs, reference, error)
        if (allocated(error)) return
        call assert_solution(solution, legacy_solution, check_tolerance(), &
            "hybrid EqTri result differs from legacy", error)
    end subroutine hybrid_eqtri_operator

    subroutine hybrid_initial_solution_and_failures(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4
        real(prec), target :: lower(n), diagonal(n), upper(n)
        real(prec) :: rhs(n), reference(n), solution(n)
        type(banded_operator_t) :: operator
        type(bicgstab_options_t) :: options
        type(bicgstab_result_t) :: result
        type(bicgstab_workspace_t) :: workspace
        type(jacobi_preconditioner_t) :: jacobi
        type(user_preconditioner_t) :: unbound_preconditioner
        integer :: status

        lower = [0.0_prec, -1.0_prec, -1.0_prec, -1.0_prec]
        diagonal = 4.0_prec
        upper = [-1.0_prec, -1.0_prec, -1.0_prec, 0.0_prec]
        reference = [1.0_prec, 2.0_prec, -1.0_prec, 0.5_prec]
        rhs = diagonal*reference
        rhs(2:n) = rhs(2:n) + lower(2:n)*reference(1:n-1)
        rhs(1:n-1) = rhs(1:n-1) + upper(1:n-1)*reference(2:n)

        call operator%bind_tridiagonal(lower, diagonal, upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: edge-case operator binding failed")
        if (allocated(error)) return
        call configure_options(options, 20)

        solution = reference
        call solve_bicgstab(operator, rhs, solution, options, workspace, result)
        call assert_result_success(result, "exact initial solution was rejected", error)
        if (allocated(error)) return
        call check(error, result%iterations, 0, &
            message="[SYSTEM] TEST ERROR: exact initial solution required iterations")
        if (allocated(error)) return

        solution = 0.0_prec
        call solve_bicgstab(operator, rhs, solution, options, workspace, result, &
            unbound_preconditioner)
        call check(error, result%status, BICGSTAB_PRECONDITIONER_FAILURE, &
            message="[SYSTEM] TEST ERROR: callback preconditioner failure was not reported")
        if (allocated(error)) return

        diagonal(2) = 0.0_prec
        call operator%bind_tridiagonal(lower, diagonal, upper, status)
        call jacobi%setup(operator, status)
        call check(error, status, BICGSTAB_COMPONENT_FAILURE, &
            message="[SYSTEM] TEST ERROR: Jacobi accepted a zero diagonal pivot")
    end subroutine hybrid_initial_solution_and_failures

    subroutine build_classic_fixture(a, b, c, d, e, matrix, reference)
        real(prec), intent(out) :: a(9), b(9), c(9), d(9), e(9)
        real(prec), intent(out) :: matrix(9,9), reference(9)

        a = [0.0_prec, -0.8_prec, -0.9_prec, 0.0_prec, -0.7_prec, &
             -0.6_prec, 0.0_prec, -1.0_prec, -0.8_prec]
        b = [6.0_prec, 6.5_prec, 7.0_prec, 6.2_prec, 6.8_prec, &
             7.2_prec, 6.4_prec, 6.9_prec, 7.1_prec]
        c = [-1.0_prec, -1.1_prec, 0.0_prec, -1.2_prec, -1.0_prec, &
             0.0_prec, -0.9_prec, -1.1_prec, 0.0_prec]
        d = [0.0_prec, 0.0_prec, 0.0_prec, -0.4_prec, -0.5_prec, &
             -0.6_prec, -0.7_prec, -0.8_prec, -0.9_prec]
        e = [-0.5_prec, -0.6_prec, -0.7_prec, -0.8_prec, -0.9_prec, &
             -1.0_prec, 0.0_prec, 0.0_prec, 0.0_prec]
        reference = [1.0_prec, -2.0_prec, 0.5_prec, 3.0_prec, -1.5_prec, &
                     2.25_prec, -0.75_prec, 1.25_prec, -3.0_prec]

        matrix = 0.0_prec
        matrix(1,[1,2,4]) = [6.0_prec,-1.0_prec,-0.5_prec]
        matrix(2,[1,2,3,5]) = [-0.8_prec,6.5_prec,-1.1_prec,-0.6_prec]
        matrix(3,[2,3,6]) = [-0.9_prec,7.0_prec,-0.7_prec]
        matrix(4,[1,4,5,7]) = [-0.4_prec,6.2_prec,-1.2_prec,-0.8_prec]
        matrix(5,[2,4,5,6,8]) = [-0.5_prec,-0.7_prec,6.8_prec,-1.0_prec,-0.9_prec]
        matrix(6,[3,5,6,9]) = [-0.6_prec,-0.6_prec,7.2_prec,-1.0_prec]
        matrix(7,[4,7,8]) = [-0.7_prec,6.4_prec,-0.9_prec]
        matrix(8,[5,7,8,9]) = [-0.8_prec,-1.0_prec,6.9_prec,-1.1_prec]
        matrix(9,[6,8,9]) = [-0.9_prec,-0.8_prec,7.1_prec]
    end subroutine build_classic_fixture

    subroutine build_triiso_fixture(a, b, c, d, e, matrix, reference)
        real(prec), intent(out) :: a(8), b(8), c(8), d(8), e(8)
        real(prec), intent(out) :: matrix(8,8), reference(8)

        a = [0.0_prec,0.0_prec,-0.8_prec,0.0_prec,0.0_prec,-0.7_prec, &
             0.0_prec,-0.9_prec]
        b = [4.5_prec,5.0_prec,4.8_prec,5.2_prec,4.6_prec,5.1_prec, &
             4.9_prec,5.3_prec]
        c = [0.0_prec,-0.8_prec,0.0_prec,0.0_prec,-0.7_prec,0.0_prec, &
             -0.9_prec,0.0_prec]
        d = [0.0_prec,0.0_prec,0.0_prec,0.0_prec,-1.1_prec,-1.2_prec, &
             -1.0_prec,-1.3_prec]
        e = [-1.1_prec,-1.2_prec,-1.0_prec,-1.3_prec,0.0_prec,0.0_prec, &
             0.0_prec,0.0_prec]
        reference = [0.75_prec,1.20_prec,0.90_prec,1.30_prec,1.10_prec, &
                     0.90_prec,1.25_prec,0.95_prec]

        matrix = 0.0_prec
        matrix(1,1)=b(1); matrix(2,2)=b(2); matrix(3,3)=b(3); matrix(4,4)=b(4)
        matrix(5,5)=b(5); matrix(6,6)=b(6); matrix(7,7)=b(7); matrix(8,8)=b(8)
        matrix(1,5)=e(1); matrix(2,6)=e(2); matrix(3,7)=e(3); matrix(4,8)=e(4)
        matrix(5,1)=d(5); matrix(6,2)=d(6); matrix(7,3)=d(7); matrix(8,4)=d(8)
        matrix(2,3)=c(2); matrix(3,2)=a(3)
        matrix(5,6)=c(5); matrix(6,5)=a(6)
        matrix(7,8)=c(7); matrix(8,7)=a(8)
    end subroutine build_triiso_fixture

    subroutine configure_options(options, max_iterations)
        type(bicgstab_options_t), intent(out) :: options
        integer, intent(in) :: max_iterations

        options%max_iterations = max_iterations
        options%relative_tolerance = solver_tolerance()
        options%absolute_tolerance = 0.0_prec
        options%breakdown_tolerance = 100.0_prec*epsilon(1.0_prec)
    end subroutine configure_options

    real(prec) function solver_tolerance()
        solver_tolerance = max(1.0e-12_prec, 100.0_prec*epsilon(1.0_prec))
    end function solver_tolerance

    real(prec) function check_tolerance()
        check_tolerance = max(1.0e-10_prec, 2000.0_prec*epsilon(1.0_prec))
    end function check_tolerance

    subroutine assert_result_success(result, detail, error)
        type(bicgstab_result_t), intent(in) :: result
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        call check(error, result%status, BICGSTAB_SUCCESS, &
            message="[SYSTEM] TEST ERROR: "//detail//"; "//trim(result%reason))
    end subroutine assert_result_success

    subroutine assert_preconditioned_solution(name, result, matrix, solution, &
                                              rhs, reference, error)
        character(len=*), intent(in) :: name
        type(bicgstab_result_t), intent(in) :: result
        real(prec), intent(in) :: matrix(:,:), solution(:), rhs(:), reference(:)
        type(error_type), allocatable, intent(out) :: error

        call assert_result_success(result, name//" solver failed", error)
        if (allocated(error)) return
        call assert_solution(solution, reference, check_tolerance(), &
            name//" solution mismatch", error)
        if (allocated(error)) return
        call assert_dense_residual(matrix, solution, rhs, check_tolerance(), &
            name//" residual too large", error)
    end subroutine assert_preconditioned_solution

    subroutine assert_solution(input, reference, tolerance, detail, error)
        real(prec), intent(in) :: input(:), reference(:), tolerance
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        if (maxval(abs(input-reference)/max(1.0_prec,abs(reference))) > tolerance) then
            call check(error, .false., message="[SYSTEM] TEST ERROR: "//detail)
        end if
    end subroutine assert_solution

    subroutine assert_dense_residual(matrix, solution, rhs, tolerance, detail, error)
        real(prec), intent(in) :: matrix(:,:), solution(:), rhs(:), tolerance
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        real(prec) :: residual(size(rhs)), residual_norm, rhs_norm

        residual = matmul(matrix,solution) - rhs
        residual_norm = sqrt(dot_product(residual,residual))
        rhs_norm = max(1.0_prec, sqrt(dot_product(rhs,rhs)))
        if (residual_norm > tolerance*rhs_norm) then
            call check(error, .false., message="[SYSTEM] TEST ERROR: "//detail)
        end if
    end subroutine assert_dense_residual

    subroutine custom_jacobi_apply(context, input, output, status)
        class(*), intent(inout) :: context
        real(prec), intent(in) :: input(:)
        real(prec), intent(out) :: output(:)
        integer, intent(out) :: status

        select type (context)
        type is (custom_jacobi_context_t)
            if (.not. allocated(context%inverse_diagonal)) then
                output = 0.0_prec
                status = 1
            elseif (size(context%inverse_diagonal) /= size(input) .or. &
                    size(output) /= size(input)) then
                output = 0.0_prec
                status = 1
            else
                output = context%inverse_diagonal*input
                status = 0
            end if
        class default
            output = 0.0_prec
            status = 1
        end select
    end subroutine custom_jacobi_apply

end module BiCGSTAB_hybrid_suite
