module Eigvals_Solver_suite
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
        ieee_positive_inf, ieee_quiet_nan
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use precision_kinds, only : prec
    use Variables, only : Options_Data
    use BiCGSTAB_Solver, only : banded_operator_t, bicgstab_options_t, &
        jacobi_preconditioner_t
    use Eigen_spectrum_solver, only : eigensolver_options_t, &
        eigensolver_result_t, eigensolver_workspace_t, &
        dense_generalized_problem_t, vectorized_generalized_problem_t, &
        solve_eigenproblem, EIGEN_MODE_DOMINANT, EIGEN_MODE_FULL, &
        EIGEN_SUCCESS, EIGEN_FACTORIZATION_FAILURE, EIGEN_INVALID_INPUT, &
        EIGEN_ZERO_START_VECTOR, EIGEN_NOT_CONVERGED, &
        EIGEN_OPERATOR_FAILURE, EIGEN_INNER_SOLVER_FAILURE, EIGEN_NONFINITE_VALUE, &
        EIGEN_INFINITE_EIGENVALUE, EIGEN_SINGULAR_PENCIL

    implicit none
    private

    public :: collect_suite_Eigvals_Solver

contains

    subroutine collect_suite_Eigvals_Solver(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("eigen_spectrum_option_defaults_disabled", &
                eigen_spectrum_option_defaults_disabled), &
            new_unittest("full_one_by_one_spectrum", &
                full_one_by_one_spectrum), &
            new_unittest("full_generalized_real_spectrum", &
                full_generalized_real_spectrum), &
            new_unittest("full_complex_conjugate_spectrum", &
                full_complex_conjugate_spectrum), &
            new_unittest("full_repeated_and_zero_eigenvalues", &
                full_repeated_and_zero_eigenvalues), &
            new_unittest("repeated_complex_pair_grouping", &
                repeated_complex_pair_grouping), &
            new_unittest("close_complex_pair_grouping", &
                close_complex_pair_grouping), &
            new_unittest("dominant_dense_arnoldi", &
                dominant_dense_arnoldi), &
            new_unittest("dominant_complex_pair_preserved", &
                dominant_complex_pair_preserved), &
            new_unittest("tiny_dominant_complex_pair", &
                tiny_dominant_complex_pair), &
            new_unittest("dominant_happy_breakdown", &
                dominant_happy_breakdown), &
            new_unittest("vectorized_tridiagonal_parity", &
                vectorized_tridiagonal_parity), &
            new_unittest("vectorized_triiso_parity", &
                vectorized_triiso_parity), &
            new_unittest("vectorized_heptadiagonal_parity", &
                vectorized_heptadiagonal_parity), &
            new_unittest("vectorized_full_complex_spectrum", &
                vectorized_full_complex_spectrum), &
            new_unittest("vectorized_singular_pencil_rejected", &
                vectorized_singular_pencil_rejected), &
            new_unittest("vectorized_preconditioned_complex_pair", &
                vectorized_preconditioned_complex_pair), &
            new_unittest("extreme_scale_residual_guards", &
                extreme_scale_residual_guards), &
            new_unittest("scale_invariant_dense_problem", &
                scale_invariant_dense_problem), &
            new_unittest("eigensolver_failure_paths", &
                eigensolver_failure_paths) &
            ]
    end subroutine collect_suite_Eigvals_Solver

    subroutine eigen_spectrum_option_defaults_disabled(error)
        type(error_type), allocatable, intent(out) :: error

        type(Options_Data) :: options

        call check(error, .not. options%eigen_spectrum_analysis, &
            message="[SYSTEM] TEST ERROR: eigen spectrum analysis default is enabled")
    end subroutine eigen_spectrum_option_defaults_disabled

    subroutine full_one_by_one_spectrum(error)
        type(error_type), allocatable, intent(out) :: error

        real(prec), target :: left_matrix(1,1), right_matrix(1,1)
        complex(prec) :: expected(1)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        left_matrix(1,1) = 2.0_prec
        right_matrix(1,1) = 6.0_prec
        expected(1) = cmplx(3.0_prec,0.0_prec,kind=prec)
        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: 1x1 eigenproblem binding failed")
        if (allocated(error)) return

        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, 1, "1x1 full spectrum", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "1x1 generalized eigenvalue", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "1x1 full spectrum", error)
    end subroutine full_one_by_one_spectrum

    subroutine full_generalized_real_spectrum(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec) :: left_before(n,n), right_before(n,n), transformed(n,n)
        complex(prec) :: expected(n)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        left_matrix(1,:) = [2.0_prec, 0.0_prec, 0.0_prec, 0.0_prec]
        left_matrix(2,:) = [0.5_prec, 3.0_prec, 0.0_prec, 0.0_prec]
        left_matrix(3,:) = [0.0_prec, -0.4_prec, 4.0_prec, 0.0_prec]
        left_matrix(4,:) = [0.2_prec, 0.0_prec, 0.3_prec, 5.0_prec]

        transformed = 0.0_prec
        transformed(1,1) = 1.0_prec
        transformed(2,2) = 4.0_prec
        transformed(3,3) = -2.0_prec
        transformed(4,4) = 0.5_prec
        transformed(1,2) = 0.7_prec
        transformed(2,4) = -0.3_prec
        right_matrix = matmul(left_matrix, transformed)
        left_before = left_matrix
        right_before = right_matrix
        expected = [cmplx(4.0_prec,0.0_prec,kind=prec), &
                    cmplx(-2.0_prec,0.0_prec,kind=prec), &
                    cmplx(1.0_prec,0.0_prec,kind=prec), &
                    cmplx(0.5_prec,0.0_prec,kind=prec)]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: dense generalized binding failed")
        if (allocated(error)) return

        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, n, "full real spectrum", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "full generalized real eigenvalues", error)
        if (allocated(error)) return
        call assert_residuals(result, "full generalized real spectrum", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "full generalized real spectrum", error)
        if (allocated(error)) return
        call check(error, all(left_matrix == left_before) .and. &
            all(right_matrix == right_before), &
            message="[SYSTEM] TEST ERROR: full solver modified input matrices")
    end subroutine full_generalized_real_spectrum

    subroutine full_complex_conjugate_spectrum(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        complex(prec) :: expected(n)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        right_matrix(1,:) = [0.0_prec, -2.0_prec]
        right_matrix(2,:) = [2.0_prec,  0.0_prec]
        expected = [cmplx(0.0_prec, 2.0_prec, kind=prec), &
                    cmplx(0.0_prec,-2.0_prec, kind=prec)]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: complex-pair problem binding failed")
        if (allocated(error)) return
        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, n, "full complex spectrum", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "full complex conjugate eigenvalues", error)
        if (allocated(error)) return
        call assert_residuals(result, "full complex conjugate spectrum", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "full complex conjugate spectrum", error)
        if (allocated(error)) return
        call check(error, maxval(abs(aimag(result%eigenvectors))) > 0.1_prec, &
            message="[SYSTEM] TEST ERROR: complex eigenvector information was discarded")
    end subroutine full_complex_conjugate_spectrum

    subroutine full_repeated_and_zero_eigenvalues(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 5
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        complex(prec) :: expected(n)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,1) = 3.0_prec
        right_matrix(2,2) = 3.0_prec
        right_matrix(3,3) = -1.0_prec
        expected = [cmplx(3.0_prec,0.0_prec,kind=prec), &
                    cmplx(3.0_prec,0.0_prec,kind=prec), &
                    cmplx(-1.0_prec,0.0_prec,kind=prec), &
                    cmplx(0.0_prec,0.0_prec,kind=prec), &
                    cmplx(0.0_prec,0.0_prec,kind=prec)]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: repeated-spectrum problem binding failed")
        if (allocated(error)) return
        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, n, "repeated and zero spectrum", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "repeated and zero eigenvalues", error)
        if (allocated(error)) return
        call assert_residuals(result, "repeated and zero spectrum", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "repeated and zero spectrum", error)
    end subroutine full_repeated_and_zero_eigenvalues

    subroutine repeated_complex_pair_grouping(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4
        real(prec) :: left_matrix(n,n), right_matrix(n,n), initial_vector(n)
        complex(prec) :: expected_full(n), expected_dominant(2)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,2) = -2.0_prec
        right_matrix(2,1) =  2.0_prec
        right_matrix(3,4) = -2.0_prec
        right_matrix(4,3) =  2.0_prec
        initial_vector = [1.0_prec, 0.5_prec, 0.25_prec, 0.75_prec]
        expected_full = [cmplx(0.0_prec, 2.0_prec, kind=prec), &
                         cmplx(0.0_prec,-2.0_prec, kind=prec), &
                         cmplx(0.0_prec, 2.0_prec, kind=prec), &
                         cmplx(0.0_prec,-2.0_prec, kind=prec)]
        expected_dominant = expected_full(1:2)

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: repeated complex problem binding failed")
        if (allocated(error)) return

        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, n, &
            "full repeated complex spectrum", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected_full, &
            "grouped repeated conjugate pairs", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "full repeated complex spectrum", error)
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 2, &
            "dominant repeated complex pair", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected_dominant, &
            "dominant repeated conjugate pair", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "dominant repeated complex pair", error)
    end subroutine repeated_complex_pair_grouping

    subroutine close_complex_pair_grouping(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4
        real(prec) :: left_matrix(n,n), right_matrix(n,n), initial_vector(n)
        real(prec) :: high_frequency
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        high_frequency = 1.0_prec + 50.0_prec*epsilon(1.0_prec)
        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,2) = -high_frequency
        right_matrix(2,1) =  high_frequency
        right_matrix(3,4) = -1.0_prec
        right_matrix(4,3) =  1.0_prec
        initial_vector = [1.0_prec, 0.5_prec, 0.25_prec, 0.75_prec]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: close complex problem binding failed")
        if (allocated(error)) return

        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, n, "close complex spectrum", error)
        if (allocated(error)) return
        call check(error, result%eigenvalues(2) == conjg(result%eigenvalues(1)) .and. &
                          result%eigenvalues(4) == conjg(result%eigenvalues(3)), &
            message="[SYSTEM] TEST ERROR: close conjugate pairs were cross-matched")
        if (allocated(error)) return
        call check(error, aimag(result%eigenvalues(1)) > &
                          aimag(result%eigenvalues(3)), &
            message="[SYSTEM] TEST ERROR: close complex magnitudes were misordered")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 2, &
            "dominant close complex pair", error)
        if (allocated(error)) return
        call check(error, result%eigenvalues(2) == conjg(result%eigenvalues(1)), &
            message="[SYSTEM] TEST ERROR: dominant close pair was cross-matched")
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "dominant close complex pair", error)
    end subroutine close_complex_pair_grouping

    subroutine dominant_dense_arnoldi(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 6
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec) :: initial_vector(n), transformed(n,n)
        real(prec) :: left_copy(n,n), right_copy(n,n)
        complex(prec) :: expected(3)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: i, status

        left_matrix = 0.0_prec
        transformed = 0.0_prec
        do i = 1, n
            left_matrix(i,i) = 2.0_prec + 0.25_prec*real(i,prec)
            transformed(i,i) = real(n-i+1,prec)
            initial_vector(i) = 1.0_prec + 0.1_prec*real(i,prec)
        end do
        transformed(1,2) = 0.2_prec
        transformed(2,3) = -0.3_prec
        transformed(4,6) = 0.4_prec
        right_matrix = matmul(left_matrix, transformed)
        left_copy = left_matrix
        right_copy = right_matrix
        expected = [cmplx(6.0_prec,0.0_prec,kind=prec), &
                    cmplx(5.0_prec,0.0_prec,kind=prec), &
                    cmplx(4.0_prec,0.0_prec,kind=prec)]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: dominant dense problem binding failed")
        if (allocated(error)) return
        call configure_dominant_options(options, 3, n)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 3, "dominant dense Arnoldi", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "dominant dense eigenvalues", error)
        if (allocated(error)) return
        call assert_residuals(result, "dominant dense Arnoldi", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "dominant dense Arnoldi", error)
        if (allocated(error)) return
        call check(error, all(left_matrix == left_copy) .and. &
                          all(right_matrix == right_copy), &
            message="[SYSTEM] TEST ERROR: dominant solve modified caller matrices")
        if (allocated(error)) return
        call check(error, result%linear_solves, n, &
            message="[SYSTEM] TEST ERROR: dense LU was not reused once per Arnoldi step")
    end subroutine dominant_dense_arnoldi

    subroutine dominant_complex_pair_preserved(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 3
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec) :: initial_vector(n)
        complex(prec) :: expected(2)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,2) = -3.0_prec
        right_matrix(2,1) = 3.0_prec
        right_matrix(3,3) = 1.0_prec
        initial_vector = [1.0_prec, 0.5_prec, 0.25_prec]
        expected = [cmplx(0.0_prec, 3.0_prec, kind=prec), &
                    cmplx(0.0_prec,-3.0_prec, kind=prec)]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: dominant complex problem binding failed")
        if (allocated(error)) return
        call configure_dominant_options(options, 1, n)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 2, "dominant complex pair", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "dominant conjugate pair", error)
        if (allocated(error)) return
        call assert_residuals(result, "dominant conjugate pair", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "dominant conjugate pair", error)
    end subroutine dominant_complex_pair_preserved

    subroutine tiny_dominant_complex_pair(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4
        real(prec) :: left_matrix(n,n), right_matrix(n,n), initial_vector(n)
        real(prec) :: spectrum_scale
        complex(prec) :: scaled_values(2), expected_scaled(2)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        spectrum_scale = sqrt(tiny(1.0_prec))
        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,2) = -3.0_prec*spectrum_scale
        right_matrix(2,1) =  3.0_prec*spectrum_scale
        right_matrix(3,3) =  2.0_prec*spectrum_scale
        right_matrix(4,4) =  1.0_prec*spectrum_scale
        initial_vector = [1.0_prec, 0.5_prec, 0.25_prec, 0.75_prec]
        expected_scaled = [cmplx(0.0_prec, 3.0_prec, kind=prec), &
                           cmplx(0.0_prec,-3.0_prec, kind=prec)]

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: tiny-spectrum problem binding failed")
        if (allocated(error)) return
        call configure_dominant_options(options, 1, n)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 2, &
            "tiny dominant conjugate pair", error)
        if (allocated(error)) return

        scaled_values = result%eigenvalues/ &
            cmplx(spectrum_scale,0.0_prec,kind=prec)
        call check(error, maxval(abs(scaled_values-expected_scaled)) <= &
                          10.0_prec*test_tolerance(), &
            message="[SYSTEM] TEST ERROR: tiny dominant spectrum was misordered")
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "tiny dominant conjugate pair", error)
    end subroutine tiny_dominant_complex_pair

    subroutine dominant_happy_breakdown(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 3
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec) :: initial_vector(n)
        complex(prec) :: expected(1)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,1) = 5.0_prec
        right_matrix(2,2) = 3.0_prec
        right_matrix(3,3) = 1.0_prec
        initial_vector = [1.0_prec,0.0_prec,0.0_prec]
        expected(1) = cmplx(5.0_prec,0.0_prec,kind=prec)

        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: happy-breakdown problem binding failed")
        if (allocated(error)) return
        ! Request one extra Arnoldi step so that the exact starting
        ! eigenvector triggers a happy breakdown and the complement logic is
        ! exercised before the dominant eigenpair is selected.
        call configure_dominant_options(options, 1, 2)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 1, "happy Arnoldi breakdown", error)
        if (allocated(error)) return
        call check(error, result%arnoldi_breakdowns, 2, &
            message="[SYSTEM] TEST ERROR: happy-breakdown count is incorrect")
        if (allocated(error)) return
        call check(error, result%arnoldi_steps, 2, &
            message="[SYSTEM] TEST ERROR: complement continuation used wrong step count")
        if (allocated(error)) return
        call check(error, result%linear_solves, 2, &
            message="[SYSTEM] TEST ERROR: complement continuation used wrong solve count")
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "happy-breakdown eigenvalue", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "happy Arnoldi breakdown", error)
    end subroutine dominant_happy_breakdown

    subroutine vectorized_tridiagonal_parity(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 8
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec), target :: lower(n), diagonal(n), upper(n)
        real(prec), target :: zero_lower(n), right_diagonal(n), zero_upper(n)
        real(prec) :: initial_vector(n)
        type(banded_operator_t), target :: left_operator, right_operator
        integer :: status

        call build_tridiagonal_fixture(left_matrix, right_matrix, lower, diagonal, &
            upper, zero_lower, right_diagonal, zero_upper, initial_vector)
        call left_operator%bind_tridiagonal(lower, diagonal, upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: tridiagonal left operator binding failed")
        if (allocated(error)) return
        call right_operator%bind_tridiagonal(zero_lower, right_diagonal, zero_upper, &
                                             status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: tridiagonal right operator binding failed")
        if (allocated(error)) return

        call compare_bound_dense_and_vectorized(left_matrix, right_matrix, &
            left_operator, right_operator, initial_vector, error)
    end subroutine vectorized_tridiagonal_parity

    subroutine vectorized_triiso_parity(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 4
        integer, parameter :: n = 8
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec), target :: a(n), b(n), c(n), d(n), e(n)
        real(prec), target :: fa(n), fb(n), fc(n), fd(n), fe(n)
        real(prec) :: initial_vector(n)
        type(banded_operator_t), target :: left_operator, right_operator
        integer :: status

        call build_triiso_fixture(left_matrix, right_matrix, a, b, c, d, e, &
            fa, fb, fc, fd, fe, initial_vector)
        call left_operator%bind_triiso(n_x, a, b, c, d, e, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: TriIso left operator binding failed")
        if (allocated(error)) return
        call right_operator%bind_triiso(n_x, fa, fb, fc, fd, fe, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: TriIso right operator binding failed")
        if (allocated(error)) return

        call compare_bound_dense_and_vectorized(left_matrix, right_matrix, &
            left_operator, right_operator, initial_vector, error)
    end subroutine vectorized_triiso_parity

    subroutine vectorized_heptadiagonal_parity(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 12
        integer, parameter :: offset_1 = 1, offset_2 = 4, offset_3 = 7
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec), target :: lower_1(n), diagonal(n), upper_1(n)
        real(prec), target :: lower_2(n), upper_2(n), lower_3(n), upper_3(n)
        real(prec), target :: zeros_1(n), right_diagonal(n), zeros_2(n)
        real(prec), target :: zeros_3(n), zeros_4(n), zeros_5(n), zeros_6(n)
        real(prec) :: initial_vector(n)
        type(banded_operator_t), target :: left_operator, right_operator
        integer :: i, status

        lower_1 = -0.10_prec
        upper_1 = -0.08_prec
        lower_2 = -0.06_prec
        upper_2 = -0.05_prec
        lower_3 = -0.04_prec
        upper_3 = -0.03_prec
        diagonal = 5.0_prec
        zeros_1 = 0.0_prec
        zeros_2 = 0.0_prec
        zeros_3 = 0.0_prec
        zeros_4 = 0.0_prec
        zeros_5 = 0.0_prec
        zeros_6 = 0.0_prec

        left_matrix = 0.0_prec
        right_matrix = 0.0_prec
        do i = 1, n
            left_matrix(i,i) = diagonal(i)
            right_diagonal(i) = real(n-i+3,prec)
            right_matrix(i,i) = right_diagonal(i)
            initial_vector(i) = 0.6_prec + real(i,prec)/real(2*n,prec)
        end do
        do i = offset_1 + 1, n
            left_matrix(i,i-offset_1) = lower_1(i)
        end do
        do i = 1, n - offset_1
            left_matrix(i,i+offset_1) = upper_1(i)
        end do
        do i = offset_2 + 1, n
            left_matrix(i,i-offset_2) = lower_2(i)
        end do
        do i = 1, n - offset_2
            left_matrix(i,i+offset_2) = upper_2(i)
        end do
        do i = offset_3 + 1, n
            left_matrix(i,i-offset_3) = lower_3(i)
        end do
        do i = 1, n - offset_3
            left_matrix(i,i+offset_3) = upper_3(i)
        end do

        call left_operator%bind_heptadiagonal(offset_1, offset_2, offset_3, &
            lower_1, diagonal, upper_1, lower_2, upper_2, lower_3, upper_3, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: heptadiagonal left operator binding failed")
        if (allocated(error)) return
        call right_operator%bind_heptadiagonal(offset_1, offset_2, offset_3, &
            zeros_1, right_diagonal, zeros_2, zeros_3, zeros_4, zeros_5, &
            zeros_6, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: heptadiagonal right operator binding failed")
        if (allocated(error)) return

        call compare_bound_dense_and_vectorized(left_matrix, right_matrix, &
            left_operator, right_operator, initial_vector, error)
    end subroutine vectorized_heptadiagonal_parity

    subroutine vectorized_full_complex_spectrum(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(prec), target :: left_lower(n), left_diagonal(n), left_upper(n)
        real(prec), target :: right_lower(n), right_diagonal(n), right_upper(n)
        real(prec) :: left_matrix(n,n), right_matrix(n,n)
        complex(prec) :: expected(n)
        type(banded_operator_t), target :: left_operator, right_operator
        type(vectorized_generalized_problem_t) :: problem
        type(bicgstab_options_t) :: inner_options
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        left_lower = 0.0_prec
        left_diagonal = 1.0_prec
        left_upper = 0.0_prec
        right_lower = [0.0_prec,2.0_prec]
        right_diagonal = 0.0_prec
        right_upper = [-2.0_prec,0.0_prec]
        call build_identity(left_matrix)
        right_matrix = reshape([0.0_prec,2.0_prec,-2.0_prec,0.0_prec], [n,n])
        expected = [cmplx(0.0_prec,2.0_prec,kind=prec), &
                    cmplx(0.0_prec,-2.0_prec,kind=prec)]

        call left_operator%bind_tridiagonal(left_lower, left_diagonal, &
                                            left_upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: full-vector left binding failed")
        if (allocated(error)) return
        call right_operator%bind_tridiagonal(right_lower, right_diagonal, &
                                             right_upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: full-vector right binding failed")
        if (allocated(error)) return
        inner_options%max_iterations = 0
        call problem%bind(left_operator, right_operator, inner_options, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: full-vector problem binding failed")
        if (allocated(error)) return

        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, n, &
            "materialized vectorized full spectrum", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "materialized vectorized eigenvalues", error)
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "materialized vectorized full spectrum", error)
        if (allocated(error)) return
        call check(error, result%linear_solves, 0, &
            message="[SYSTEM] TEST ERROR: full spectrum invoked the inner solver")
    end subroutine vectorized_full_complex_spectrum

    subroutine vectorized_singular_pencil_rejected(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(prec), target :: left_lower(n), left_diagonal(n), left_upper(n)
        real(prec), target :: right_lower(n), right_diagonal(n), right_upper(n)
        real(prec) :: initial_vector(n)
        type(banded_operator_t) :: left_operator, right_operator
        type(vectorized_generalized_problem_t) :: problem
        type(bicgstab_options_t) :: inner_options
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        left_lower = 0.0_prec
        left_diagonal = [1.0_prec,0.0_prec]
        left_upper = 0.0_prec
        right_lower = 0.0_prec
        right_diagonal = 0.0_prec
        right_upper = 0.0_prec
        initial_vector = [0.0_prec,1.0_prec]

        call left_operator%bind_tridiagonal(left_lower, left_diagonal, &
                                             left_upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: singular left operator binding failed")
        if (allocated(error)) return
        call right_operator%bind_tridiagonal(right_lower, right_diagonal, &
                                              right_upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: null right operator binding failed")
        if (allocated(error)) return

        inner_options%max_iterations = 20
        inner_options%relative_tolerance = vector_tolerance()
        inner_options%absolute_tolerance = 0.0_prec
        call problem%bind(left_operator, right_operator, inner_options, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: singular vector problem binding failed")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, 1)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call check(error, result%status, EIGEN_SINGULAR_PENCIL, &
            message="[SYSTEM] TEST ERROR: common-null vector falsely converged")
        if (allocated(error)) return
        call check(error, .not. result%converged, &
            message="[SYSTEM] TEST ERROR: singular vector set convergence flag")
    end subroutine vectorized_singular_pencil_rejected

    subroutine vectorized_preconditioned_complex_pair(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(prec), target :: left_lower(n), left_diagonal(n), left_upper(n)
        real(prec), target :: right_lower(n), right_diagonal(n), right_upper(n)
        real(prec) :: left_matrix(n,n), right_matrix(n,n), initial_vector(n)
        complex(prec) :: expected(n)
        type(banded_operator_t), target :: left_operator, right_operator
        type(jacobi_preconditioner_t), target :: preconditioner
        type(vectorized_generalized_problem_t) :: problem
        type(bicgstab_options_t) :: inner_options
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        left_lower = 0.0_prec
        left_diagonal = [2.0_prec,3.0_prec]
        left_upper = 0.0_prec
        right_lower = [0.0_prec,6.0_prec]
        right_diagonal = 0.0_prec
        right_upper = [-4.0_prec,0.0_prec]
        left_matrix = reshape([2.0_prec,0.0_prec,0.0_prec,3.0_prec], [n,n])
        right_matrix = reshape([0.0_prec,6.0_prec,-4.0_prec,0.0_prec], [n,n])
        initial_vector = [1.0_prec,0.5_prec]
        expected = [cmplx(0.0_prec,2.0_prec,kind=prec), &
                    cmplx(0.0_prec,-2.0_prec,kind=prec)]

        call left_operator%bind_tridiagonal(left_lower, left_diagonal, &
                                            left_upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: preconditioned left binding failed")
        if (allocated(error)) return
        call right_operator%bind_tridiagonal(right_lower, right_diagonal, &
                                             right_upper, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: preconditioned right binding failed")
        if (allocated(error)) return
        call preconditioner%setup(left_operator, status)
        call check(error, status, 0, &
            message="[SYSTEM] TEST ERROR: Jacobi preconditioner setup failed")
        if (allocated(error)) return

        ! Exact diagonal Jacobi inversion makes every inner solve converge in
        ! one iteration; max_iterations=1 proves the supplied
        ! preconditioner is actually used.
        inner_options%max_iterations = 1
        inner_options%relative_tolerance = vector_tolerance()
        inner_options%absolute_tolerance = 0.0_prec
        call problem%bind(left_operator, right_operator, inner_options, status, &
                          preconditioner)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: preconditioned problem binding failed")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        options%residual_tolerance = vector_tolerance()
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, n, &
            "preconditioned vectorized complex pair", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "preconditioned vectorized complex eigenvalues", error, &
            vector_tolerance())
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, result, &
            "preconditioned vectorized complex pair", error, vector_tolerance())
        if (allocated(error)) return
        call check(error, result%inner_iterations, result%linear_solves, &
            message="[SYSTEM] TEST ERROR: Jacobi was not applied once per inner solve")
    end subroutine vectorized_preconditioned_complex_pair

    subroutine extreme_scale_residual_guards(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(prec) :: left_matrix(n,n), right_matrix(n,n), initial_vector(n)
        real(prec) :: matrix_scale
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        initial_vector = [1.0_prec,1.0_prec]
        call configure_dominant_options(options, 1, 1)
        options%residual_tolerance = 1.0e-6_prec

        matrix_scale = 0.90_prec*huge(1.0_prec)
        right_matrix = 0.0_prec
        right_matrix(1,1) = matrix_scale
        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: overflow-scale problem binding failed")
        if (allocated(error)) return
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_rejected_false_convergence(result, &
            "overflow-scale residual", error)
        if (allocated(error)) return

        matrix_scale = scale(tiny(1.0_prec),-8)
        call check(error, matrix_scale > 0.0_prec, &
            message="[SYSTEM] TEST ERROR: subnormal test scale is unavailable")
        if (allocated(error)) return
        right_matrix = 0.0_prec
        right_matrix(1,1) = matrix_scale
        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: underflow-scale problem binding failed")
        if (allocated(error)) return
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_rejected_false_convergence(result, &
            "underflow-scale residual", error)
    end subroutine extreme_scale_residual_guards

    subroutine scale_invariant_dense_problem(error)
        type(error_type), allocatable, intent(out) :: error

        real(prec) :: left_matrix(2,2), right_matrix(2,2), initial_vector(2)
        real(prec) :: tiny_left(1,1), tiny_right(1,1), tiny_start(1)
        real(prec) :: pencil_scale, start_scale
        complex(prec) :: expected(1)
        type(dense_generalized_problem_t) :: problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,1) = 5.0_prec
        right_matrix(2,2) = 1.0_prec
        start_scale = sqrt(tiny(1.0_prec))
        initial_vector = [start_scale,start_scale]
        expected(1) = cmplx(5.0_prec,0.0_prec,kind=prec)
        call problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: scaled-start problem binding failed")
        if (allocated(error)) return
        call configure_dominant_options(options, 1, 2)
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 1, "scaled Arnoldi start", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "scaled-start eigenvalue", error)
        if (allocated(error)) return

        start_scale = 0.90_prec*huge(1.0_prec)
        initial_vector = [start_scale,start_scale]
        call solve_eigenproblem(problem, options, workspace, result, initial_vector)
        call assert_successful_result(result, 1, "large Arnoldi start", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "large-start eigenvalue", error)
        if (allocated(error)) return

        pencil_scale = scale(tiny(1.0_prec),-8)
        tiny_left(1,1) = pencil_scale
        tiny_right(1,1) = pencil_scale
        tiny_start(1) = 1.0_prec
        expected(1) = cmplx(1.0_prec,0.0_prec,kind=prec)
        call problem%bind(tiny_left, tiny_right, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: scaled pencil binding failed")
        if (allocated(error)) return

        call configure_full_options(options)
        call solve_eigenproblem(problem, options, workspace, result)
        call assert_successful_result(result, 1, "scaled full pencil", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "scaled full-pencil eigenvalue", error)
        if (allocated(error)) return

        call configure_dominant_options(options, 1, 1)
        call solve_eigenproblem(problem, options, workspace, result, tiny_start)
        call assert_successful_result(result, 1, "scaled dominant pencil", error)
        if (allocated(error)) return
        call assert_eigenvalues(result%eigenvalues, expected, &
            "scaled dominant-pencil eigenvalue", error)
    end subroutine scale_invariant_dense_problem

    subroutine eigensolver_failure_paths(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 2
        real(prec), target :: left_matrix(n,n), right_matrix(n,n)
        real(prec), target :: lower(n), diagonal(n), upper(n)
        real(prec), target :: f_lower(n), f_diagonal(n), f_upper(n)
        real(prec) :: zero_start(n), mixed_start(n), unbound_output(n)
        real(prec) :: short_start(1)
        complex(prec) :: expected_full(n)
        type(dense_generalized_problem_t) :: dense_problem
        type(dense_generalized_problem_t) :: unbound_problem
        type(vectorized_generalized_problem_t) :: vector_problem
        type(banded_operator_t), target :: left_operator, right_operator
        type(bicgstab_options_t) :: inner_options
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: workspace
        type(eigensolver_result_t) :: result
        integer :: status

        mixed_start = 1.0_prec
        call unbound_problem%apply_left(mixed_start, unbound_output, status)
        call check(error, status, EIGEN_OPERATOR_FAILURE, &
            message="[SYSTEM] TEST ERROR: unbound dense left apply did not fail safely")
        if (allocated(error)) return
        call check(error, all(unbound_output == 0.0_prec), &
            message="[SYSTEM] TEST ERROR: unbound dense left apply returned data")
        if (allocated(error)) return
        call unbound_problem%apply_right(mixed_start, unbound_output, status)
        call check(error, status, EIGEN_OPERATOR_FAILURE, &
            message="[SYSTEM] TEST ERROR: unbound dense right apply did not fail safely")
        if (allocated(error)) return

        options = eigensolver_options_t()
        call solve_eigenproblem(unbound_problem, options, workspace, result)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: unbound dense problem was accepted")
        if (allocated(error)) return

        left_matrix = 0.0_prec
        left_matrix(1,1) = 1.0_prec
        right_matrix = 0.0_prec
        right_matrix(1,1) = 2.0_prec
        right_matrix(2,2) = 1.0_prec
        call dense_problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: singular pencil could not be represented")
        if (allocated(error)) return
        call configure_full_options(options)
        call solve_eigenproblem(dense_problem, options, workspace, result)
        call check(error, result%status, EIGEN_INFINITE_EIGENVALUE, &
            message="[SYSTEM] TEST ERROR: infinite generalized eigenvalue was accepted")
        if (allocated(error)) return
        call configure_dominant_options(options, 1, n)
        mixed_start = [1.0_prec,1.0_prec]
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_FACTORIZATION_FAILURE, &
            message="[SYSTEM] TEST ERROR: singular left solve was not rejected")
        if (allocated(error)) return

        call build_identity(left_matrix)
        right_matrix = 0.0_prec
        right_matrix(1,1) = 5.0_prec
        right_matrix(2,2) = 1.0_prec
        call dense_problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: failure-path dense binding failed")
        if (allocated(error)) return

        ! Exercise the default start-vector and full-dimension selection.
        options = eigensolver_options_t()
        call solve_eigenproblem(dense_problem, options, workspace, result)
        call assert_successful_result(result, 1, &
            "default dominant solver options", error)
        if (allocated(error)) return

        options%mode = -999
        call solve_eigenproblem(dense_problem, options, workspace, result)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: unknown eigenvalue mode was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        zero_start = 0.0_prec
        call solve_eigenproblem(dense_problem, options, workspace, result, zero_start)
        call check(error, result%status, EIGEN_ZERO_START_VECTOR, &
            message="[SYSTEM] TEST ERROR: zero Arnoldi start was accepted")
        if (allocated(error)) return

        mixed_start = [ieee_value(1.0_prec, ieee_quiet_nan), 1.0_prec]
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: nonfinite Arnoldi start was accepted")
        if (allocated(error)) return

        mixed_start = [1.0_prec, 1.0_prec]
        options%krylov_dimension = 1
        options%residual_tolerance = test_tolerance()
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_NOT_CONVERGED, &
            message="[SYSTEM] TEST ERROR: insufficient Arnoldi subspace falsely converged")
        if (allocated(error)) return

        call configure_dominant_options(options, 0, n)
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: zero eigenvalue count was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, n+1, n)
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: excessive eigenvalue count was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        options%residual_tolerance = 0.0_prec
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: zero residual tolerance was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        options%residual_tolerance = ieee_value(1.0_prec, ieee_quiet_nan)
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: nonfinite residual tolerance was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        options%residual_tolerance = ieee_value(1.0_prec, ieee_positive_inf)
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: infinite residual tolerance was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        options%breakdown_tolerance = -epsilon(1.0_prec)
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: negative breakdown tolerance was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        options%breakdown_tolerance = ieee_value(1.0_prec, ieee_quiet_nan)
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: nonfinite breakdown tolerance was accepted")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        short_start = 1.0_prec
        call solve_eigenproblem(dense_problem, options, workspace, result, short_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: wrong-sized start vector was accepted")
        if (allocated(error)) return

        options%krylov_dimension = n + 1
        call solve_eigenproblem(dense_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: invalid Krylov dimension was accepted")
        if (allocated(error)) return

        lower = 0.0_prec
        diagonal = 1.0_prec
        upper = 0.0_prec
        f_lower = 0.0_prec
        f_diagonal = [5.0_prec, 1.0_prec]
        f_upper = 0.0_prec
        call left_operator%bind_tridiagonal(lower, diagonal, upper, status)
        call right_operator%bind_tridiagonal(f_lower, f_diagonal, f_upper, status)
        inner_options%max_iterations = 0
        inner_options%relative_tolerance = 0.0_prec
        inner_options%absolute_tolerance = 0.0_prec
        call vector_problem%bind(left_operator, right_operator, inner_options, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: failure-path vector binding failed")
        if (allocated(error)) return

        call configure_dominant_options(options, 1, n)
        call solve_eigenproblem(vector_problem, options, workspace, result, mixed_start)
        call check(error, result%status, EIGEN_INNER_SOLVER_FAILURE, &
            message="[SYSTEM] TEST ERROR: inner solver failure was not propagated")
        if (allocated(error)) return
        call check(error, result%last_inner_solver_status >= 0 .and. &
                          trim(result%last_inner_solver_reason) /= &
                              'inner solver not used', &
            message="[SYSTEM] TEST ERROR: inner solver diagnostics were lost")
        if (allocated(error)) return

        options%mode = EIGEN_MODE_FULL
        call solve_eigenproblem(vector_problem, options, workspace, result)
        call assert_successful_result(result, n, &
            "full vectorized spectrum with disabled inner solver", error)
        if (allocated(error)) return
        expected_full = [cmplx(5.0_prec,0.0_prec,kind=prec), &
                         cmplx(1.0_prec,0.0_prec,kind=prec)]
        call assert_eigenvalues(result%eigenvalues, expected_full, &
            "full vectorized failure-path eigenvalues", error)
        if (allocated(error)) return

        f_diagonal(1) = ieee_value(1.0_prec, ieee_quiet_nan)
        call solve_eigenproblem(vector_problem, options, workspace, result)
        call check(error, result%status, EIGEN_NONFINITE_VALUE, &
            message="[SYSTEM] TEST ERROR: nonfinite materialized operator was accepted")
        if (allocated(error)) return

        left_matrix(1,1) = ieee_value(1.0_prec, ieee_quiet_nan)
        call dense_problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_INVALID_INPUT, &
            message="[SYSTEM] TEST ERROR: nonfinite dense matrix was accepted")
    end subroutine eigensolver_failure_paths

    subroutine compare_bound_dense_and_vectorized(left_matrix, right_matrix, &
        left_operator, right_operator, initial_vector, error)
        real(prec), target, intent(in) :: left_matrix(:,:), right_matrix(:,:)
        type(banded_operator_t), target, intent(in) :: left_operator, right_operator
        real(prec), intent(in) :: initial_vector(:)
        type(error_type), allocatable, intent(out) :: error

        type(dense_generalized_problem_t) :: dense_problem
        type(vectorized_generalized_problem_t) :: vector_problem
        type(eigensolver_options_t) :: options
        type(eigensolver_workspace_t) :: dense_workspace, vector_workspace
        type(eigensolver_result_t) :: dense_result, vector_result
        type(bicgstab_options_t) :: inner_options
        integer :: status

        call dense_problem%bind(left_matrix, right_matrix, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: parity dense problem binding failed")
        if (allocated(error)) return

        inner_options%max_iterations = 1000
        inner_options%relative_tolerance = max(1.0e-13_prec, &
                                                100.0_prec*epsilon(1.0_prec))
        inner_options%absolute_tolerance = 0.0_prec
        inner_options%breakdown_tolerance = 100.0_prec*epsilon(1.0_prec)
        call vector_problem%bind(left_operator, right_operator, inner_options, status)
        call check(error, status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: vectorized eigenproblem binding failed")
        if (allocated(error)) return

        call configure_dominant_options(options, 3, size(initial_vector))
        options%residual_tolerance = vector_tolerance()
        call solve_eigenproblem(dense_problem, options, dense_workspace, dense_result, &
                                initial_vector)
        call assert_successful_result(dense_result, 3, "dense parity reference", error)
        if (allocated(error)) return

        call solve_eigenproblem(vector_problem, options, vector_workspace, vector_result, &
                                initial_vector)
        call assert_successful_result(vector_result, 3, "vectorized parity solve", error)
        if (allocated(error)) return
        call assert_eigenvalues(vector_result%eigenvalues, dense_result%eigenvalues, &
            "dense/vectorized eigenvalue parity", error, vector_tolerance())
        if (allocated(error)) return
        call assert_residuals(vector_result, "vectorized parity solve", error, &
                              vector_tolerance())
        if (allocated(error)) return
        call assert_independent_residuals(left_matrix, right_matrix, vector_result, &
            "vectorized parity solve", error, vector_tolerance())
        if (allocated(error)) return
        call check(error, vector_result%inner_iterations > 0, &
            message="[SYSTEM] TEST ERROR: vectorized inner iterations were not recorded")
    end subroutine compare_bound_dense_and_vectorized

    subroutine build_tridiagonal_fixture(left_matrix, right_matrix, lower, diagonal, &
        upper, zero_lower, right_diagonal, zero_upper, initial_vector)
        real(prec), intent(out) :: left_matrix(:,:), right_matrix(:,:)
        real(prec), intent(out) :: lower(:), diagonal(:), upper(:)
        real(prec), intent(out) :: zero_lower(:), right_diagonal(:), zero_upper(:)
        real(prec), intent(out) :: initial_vector(:)

        integer :: i, n

        n = size(diagonal)
        lower = -0.35_prec
        upper = -0.20_prec
        lower(1) = 0.0_prec
        upper(n) = 0.0_prec
        diagonal = 4.0_prec
        zero_lower = 0.0_prec
        zero_upper = 0.0_prec
        do i = 1, n
            right_diagonal(i) = real(n-i+2,prec)
            initial_vector(i) = 0.5_prec + real(i,prec)/real(n,prec)
        end do

        left_matrix = 0.0_prec
        right_matrix = 0.0_prec
        do i = 1, n
            left_matrix(i,i) = diagonal(i)
            right_matrix(i,i) = right_diagonal(i)
        end do
        do i = 2, n
            left_matrix(i,i-1) = lower(i)
        end do
        do i = 1, n-1
            left_matrix(i,i+1) = upper(i)
        end do
    end subroutine build_tridiagonal_fixture

    subroutine build_triiso_fixture(left_matrix, right_matrix, a, b, c, d, e, &
        fa, fb, fc, fd, fe, initial_vector)
        real(prec), intent(out) :: left_matrix(8,8), right_matrix(8,8)
        real(prec), intent(out) :: a(8), b(8), c(8), d(8), e(8)
        real(prec), intent(out) :: fa(8), fb(8), fc(8), fd(8), fe(8)
        real(prec), intent(out) :: initial_vector(8)

        integer :: i

        a = [0.0_prec,0.0_prec,-0.20_prec,0.0_prec,0.0_prec,-0.15_prec, &
             0.0_prec,-0.18_prec]
        b = [4.5_prec,5.0_prec,4.8_prec,5.2_prec,4.6_prec,5.1_prec, &
             4.9_prec,5.3_prec]
        c = [0.0_prec,-0.12_prec,0.0_prec,0.0_prec,-0.14_prec,0.0_prec, &
             -0.16_prec,0.0_prec]
        d = [0.0_prec,0.0_prec,0.0_prec,0.0_prec,-0.25_prec,-0.22_prec, &
             -0.20_prec,-0.24_prec]
        e = [-0.21_prec,-0.23_prec,-0.19_prec,-0.22_prec,0.0_prec,0.0_prec, &
             0.0_prec,0.0_prec]
        fa = 0.0_prec
        fc = 0.0_prec
        fd = 0.0_prec
        fe = 0.0_prec
        do i = 1, 8
            fb(i) = real(10-i,prec)
            initial_vector(i) = 0.75_prec + 0.05_prec*real(i,prec)
        end do

        left_matrix = 0.0_prec
        right_matrix = 0.0_prec
        do i = 1, 8
            left_matrix(i,i) = b(i)
            right_matrix(i,i) = fb(i)
        end do
        left_matrix(1,5)=e(1); left_matrix(2,6)=e(2)
        left_matrix(3,7)=e(3); left_matrix(4,8)=e(4)
        left_matrix(5,1)=d(5); left_matrix(6,2)=d(6)
        left_matrix(7,3)=d(7); left_matrix(8,4)=d(8)
        left_matrix(2,3)=c(2); left_matrix(3,2)=a(3)
        left_matrix(5,6)=c(5); left_matrix(6,5)=a(6)
        left_matrix(7,8)=c(7); left_matrix(8,7)=a(8)
    end subroutine build_triiso_fixture

    subroutine configure_full_options(options)
        type(eigensolver_options_t), intent(out) :: options

        options%mode = EIGEN_MODE_FULL
        options%residual_tolerance = test_tolerance()
        options%breakdown_tolerance = 100.0_prec*epsilon(1.0_prec)
    end subroutine configure_full_options

    subroutine configure_dominant_options(options, number_of_eigenvalues, dimension)
        type(eigensolver_options_t), intent(out) :: options
        integer, intent(in) :: number_of_eigenvalues, dimension

        options%mode = EIGEN_MODE_DOMINANT
        options%number_of_eigenvalues = number_of_eigenvalues
        options%krylov_dimension = dimension
        options%residual_tolerance = test_tolerance()
        options%breakdown_tolerance = 1000.0_prec*epsilon(1.0_prec)
    end subroutine configure_dominant_options

    subroutine build_identity(matrix)
        real(prec), intent(out) :: matrix(:,:)

        integer :: i

        matrix = 0.0_prec
        do i = 1, min(size(matrix,1),size(matrix,2))
            matrix(i,i) = 1.0_prec
        end do
    end subroutine build_identity

    subroutine assert_successful_result(result, expected_count, detail, error)
        type(eigensolver_result_t), intent(in) :: result
        integer, intent(in) :: expected_count
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        integer :: index

        call check(error, result%status, EIGEN_SUCCESS, &
            message="[SYSTEM] TEST ERROR: "//detail//" failed; "//trim(result%reason))
        if (allocated(error)) return
        call check(error, result%converged, &
            message="[SYSTEM] TEST ERROR: "//detail//" did not report convergence")
        if (allocated(error)) return
        call check(error, result%returned_eigenvalues, expected_count, &
            message="[SYSTEM] TEST ERROR: "//detail//" returned wrong eigenvalue count")
        if (allocated(error)) return
        call check(error, result%converged_eigenvalues, expected_count, &
            message="[SYSTEM] TEST ERROR: "//detail//" has unconverged eigenpairs")
        if (allocated(error)) return
        call check(error, allocated(result%eigenvalues) .and. &
                          allocated(result%eigenvectors) .and. &
                          allocated(result%absolute_residuals) .and. &
                          allocated(result%relative_residuals) .and. &
                          allocated(result%ritz_residual_estimates), &
            message="[SYSTEM] TEST ERROR: "//detail//" has incomplete result arrays")
        if (allocated(error)) return
        call check(error, size(result%eigenvalues) == expected_count .and. &
                          size(result%eigenvectors,2) == expected_count .and. &
                          size(result%eigenvectors,1) > 0 .and. &
                          size(result%absolute_residuals) == expected_count .and. &
                          size(result%relative_residuals) == expected_count .and. &
                          size(result%ritz_residual_estimates) == expected_count, &
            message="[SYSTEM] TEST ERROR: "//detail//" has inconsistent result shapes")
        if (allocated(error)) return
        call check(error, all(ieee_is_finite(real(result%eigenvalues,kind=prec))) .and. &
                          all(ieee_is_finite(aimag(result%eigenvalues))) .and. &
                          all(ieee_is_finite(real(result%eigenvectors,kind=prec))) .and. &
                          all(ieee_is_finite(aimag(result%eigenvectors))) .and. &
                          all(ieee_is_finite(result%absolute_residuals)) .and. &
                          all(ieee_is_finite(result%relative_residuals)) .and. &
                          all(ieee_is_finite(result%ritz_residual_estimates)), &
            message="[SYSTEM] TEST ERROR: "//detail//" contains nonfinite results")
        if (allocated(error)) return
        call check(error, all(result%absolute_residuals >= 0.0_prec) .and. &
                          all(result%relative_residuals >= 0.0_prec) .and. &
                          all(result%ritz_residual_estimates >= 0.0_prec), &
            message="[SYSTEM] TEST ERROR: "//detail//" has negative diagnostics")
        if (allocated(error)) return
        do index = 1, expected_count
            call check(error, scaled_complex_norm(result%eigenvectors(:,index)) > &
                              0.0_prec, &
                message="[SYSTEM] TEST ERROR: "//detail//" has a zero eigenvector")
            if (allocated(error)) return
        end do
    end subroutine assert_successful_result

    subroutine assert_rejected_false_convergence(result, detail, error)
        type(eigensolver_result_t), intent(in) :: result
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        call check(error, result%status, EIGEN_NOT_CONVERGED, &
            message="[SYSTEM] TEST ERROR: "//detail//" falsely converged")
        if (allocated(error)) return
        call check(error, .not. result%converged, &
            message="[SYSTEM] TEST ERROR: "//detail//" set the convergence flag")
        if (allocated(error)) return
        call check(error, allocated(result%relative_residuals), &
            message="[SYSTEM] TEST ERROR: "//detail//" has no diagnostic residual")
        if (allocated(error)) return
        call check(error, size(result%relative_residuals) == 1, &
            message="[SYSTEM] TEST ERROR: "//detail//" has wrong residual count")
        if (allocated(error)) return
        call check(error, ieee_is_finite(result%relative_residuals(1)) .and. &
                          result%relative_residuals(1) > 0.30_prec .and. &
                          result%relative_residuals(1) < 0.55_prec, &
            message="[SYSTEM] TEST ERROR: "//detail//" residual is not O(1)")
    end subroutine assert_rejected_false_convergence

    subroutine assert_eigenvalues(actual, expected, detail, error, tolerance)
        complex(prec), intent(in) :: actual(:), expected(:)
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error
        real(prec), optional, intent(in) :: tolerance

        real(prec) :: local_tolerance

        local_tolerance = test_tolerance()
        if (present(tolerance)) local_tolerance = tolerance
        if (size(actual) /= size(expected)) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" count mismatch")
            return
        end if
        if (.not. all(ieee_is_finite(real(actual,kind=prec))) .or. &
            .not. all(ieee_is_finite(aimag(actual)))) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" contains nonfinite values")
            return
        end if
        if (maxval(abs(actual-expected)/max(1.0_prec,abs(expected))) > &
            local_tolerance) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" value mismatch")
        end if
    end subroutine assert_eigenvalues

    subroutine assert_residuals(result, detail, error, tolerance)
        type(eigensolver_result_t), intent(in) :: result
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error
        real(prec), optional, intent(in) :: tolerance

        real(prec) :: local_tolerance

        local_tolerance = test_tolerance()
        if (present(tolerance)) local_tolerance = tolerance
        if (.not. allocated(result%relative_residuals)) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" has no residuals")
            return
        end if
        if (size(result%relative_residuals) /= result%returned_eigenvalues .or. &
            .not. all(ieee_is_finite(result%relative_residuals)) .or. &
            any(result%relative_residuals < 0.0_prec)) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" has invalid residuals")
            return
        end if
        if (maxval(result%relative_residuals) > local_tolerance) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" residual too large")
        end if
    end subroutine assert_residuals

    subroutine assert_independent_residuals(left_matrix, right_matrix, result, &
                                            detail, error, tolerance)
        real(prec), intent(in) :: left_matrix(:,:), right_matrix(:,:)
        type(eigensolver_result_t), intent(in) :: result
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error
        real(prec), optional, intent(in) :: tolerance

        complex(prec) :: left_value(size(left_matrix,1))
        complex(prec) :: right_value(size(left_matrix,1))
        complex(prec) :: residual(size(left_matrix,1))
        real(prec) :: local_tolerance, relative_residual
        integer :: index

        local_tolerance = test_tolerance()
        if (present(tolerance)) local_tolerance = tolerance
        if (.not. allocated(result%eigenvalues) .or. &
            .not. allocated(result%eigenvectors)) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" has no eigenpairs")
            return
        end if
        if (size(left_matrix,2) /= size(left_matrix,1) .or. &
            size(right_matrix,1) /= size(left_matrix,1) .or. &
            size(right_matrix,2) /= size(left_matrix,1) .or. &
            size(result%eigenvectors,1) /= size(left_matrix,1) .or. &
            size(result%eigenvectors,2) /= size(result%eigenvalues)) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: "//detail//" has inconsistent dimensions")
            return
        end if
        do index = 1, size(result%eigenvalues)
            if (.not. ieee_is_finite(real(result%eigenvalues(index),kind=prec)) .or. &
                .not. ieee_is_finite(aimag(result%eigenvalues(index))) .or. &
                .not. all(ieee_is_finite(real(result%eigenvectors(:,index), &
                                                  kind=prec))) .or. &
                .not. all(ieee_is_finite(aimag(result%eigenvectors(:,index)))) .or. &
                scaled_complex_norm(result%eigenvectors(:,index)) <= 0.0_prec) then
                call check(error, .false., message= &
                    "[SYSTEM] TEST ERROR: "//detail//" has an invalid eigenpair")
                return
            end if
            left_value = matmul(cmplx(left_matrix,0.0_prec,kind=prec), &
                                result%eigenvectors(:,index))
            right_value = matmul(cmplx(right_matrix,0.0_prec,kind=prec), &
                                 result%eigenvectors(:,index))
            residual = right_value - result%eigenvalues(index)*left_value
            relative_residual = safe_test_relative_residual( &
                scaled_complex_norm(residual), scaled_complex_norm(right_value), &
                scaled_complex_norm(left_value), result%eigenvalues(index))
            if (.not. ieee_is_finite(relative_residual) .or. &
                relative_residual > local_tolerance) then
                call check(error, .false., message= &
                    "[SYSTEM] TEST ERROR: independent "//detail// &
                    " residual too large")
                return
            end if
        end do
    end subroutine assert_independent_residuals

    pure real(prec) function safe_test_relative_residual(residual_norm, &
                                                          right_norm, left_norm, &
                                                          eigenvalue)
        real(prec), intent(in) :: residual_norm, right_norm, left_norm
        complex(prec), intent(in) :: eigenvalue

        real(prec) :: denominator_scale, eigenvalue_scale, eigenvalue_shape
        real(prec) :: scaled_denominator, scaled_numerator, scaled_right

        safe_test_relative_residual = huge(1.0_prec)
        if (residual_norm == 0.0_prec) then
            safe_test_relative_residual = 0.0_prec
            return
        end if

        eigenvalue_scale = max(abs(real(eigenvalue,kind=prec)), &
                               abs(aimag(eigenvalue)))
        if (eigenvalue_scale == 0.0_prec) then
            if (right_norm == 0.0_prec) return
            safe_test_relative_residual = residual_norm/right_norm
            return
        end if
        eigenvalue_shape = sqrt( &
            (real(eigenvalue,kind=prec)/eigenvalue_scale)**2 + &
            (aimag(eigenvalue)/eigenvalue_scale)**2)

        if (eigenvalue_scale >= 1.0_prec) then
            scaled_right = right_norm/eigenvalue_scale
            denominator_scale = max(scaled_right,left_norm)
            if (denominator_scale == 0.0_prec) return
            scaled_numerator = (residual_norm/eigenvalue_scale)/denominator_scale
            scaled_denominator = scaled_right/denominator_scale + &
                eigenvalue_shape*(left_norm/denominator_scale)
        else
            denominator_scale = max(right_norm,left_norm)
            if (denominator_scale == 0.0_prec) return
            scaled_numerator = residual_norm/denominator_scale
            scaled_denominator = right_norm/denominator_scale + &
                eigenvalue_scale*eigenvalue_shape*(left_norm/denominator_scale)
        end if
        if (scaled_denominator > 0.0_prec) &
            safe_test_relative_residual = scaled_numerator/scaled_denominator
    end function safe_test_relative_residual

    pure real(prec) function scaled_complex_norm(vector)
        complex(prec), intent(in) :: vector(:)

        real(prec) :: component_value, norm_scale, sum_of_squares
        integer :: component, index

        norm_scale = 0.0_prec
        sum_of_squares = 1.0_prec
        do index = 1, size(vector)
            do component = 1, 2
                if (component == 1) then
                    component_value = abs(real(vector(index),kind=prec))
                else
                    component_value = abs(aimag(vector(index)))
                end if
                if (component_value > 0.0_prec) then
                    if (norm_scale < component_value) then
                        sum_of_squares = 1.0_prec + &
                            sum_of_squares*(norm_scale/component_value)**2
                        norm_scale = component_value
                    else
                        sum_of_squares = sum_of_squares + &
                            (component_value/norm_scale)**2
                    end if
                end if
            end do
        end do
        if (norm_scale == 0.0_prec) then
            scaled_complex_norm = 0.0_prec
        else
            scaled_complex_norm = norm_scale*sqrt(sum_of_squares)
        end if
    end function scaled_complex_norm

    pure real(prec) function test_tolerance()
        test_tolerance = max(1.0e-10_prec, 10000.0_prec*epsilon(1.0_prec))
    end function test_tolerance

    pure real(prec) function vector_tolerance()
        vector_tolerance = max(1.0e-8_prec, 100000.0_prec*epsilon(1.0_prec))
    end function vector_tolerance

end module Eigvals_Solver_suite
