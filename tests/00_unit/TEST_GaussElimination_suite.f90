module GaussElimination_suite
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use precision_kinds, only : prec
    use GaussElimination, only : GE_Expl_Main

    implicit none
    private

    public :: collect_suite_GaussElimination

contains

    subroutine collect_suite_GaussElimination(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("gauss_elimination_one_by_one", &
                gauss_elimination_one_by_one), &
            new_unittest("gauss_elimination_identity", &
                gauss_elimination_identity), &
            new_unittest("gauss_elimination_partial_pivoting", &
                gauss_elimination_partial_pivoting), &
            new_unittest("gauss_elimination_nonsymmetric", &
                gauss_elimination_nonsymmetric) &
            ]
    end subroutine collect_suite_GaussElimination

    subroutine gauss_elimination_one_by_one(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 1
        real(prec) :: matrix(n,n), matrix_before(n,n)
        real(prec) :: rhs(n), rhs_before(n), reference(n), solution(n)

        matrix(1,1) = 4.0_prec
        reference = [2.5_prec]
        rhs = matmul(matrix, reference)
        matrix_before = matrix
        rhs_before = rhs

        call GE_Expl_Main(n, matrix, rhs, solution)
        call assert_dense_solution(matrix, rhs, reference, solution, &
            "1x1 system", error)
        if (allocated(error)) return
        call assert_inputs_unchanged(matrix, matrix_before, rhs, rhs_before, &
            "1x1 system", error)
    end subroutine gauss_elimination_one_by_one

    subroutine gauss_elimination_identity(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 4
        real(prec) :: matrix(n,n), matrix_before(n,n)
        real(prec) :: rhs(n), rhs_before(n), reference(n), solution(n)
        integer :: i

        matrix = 0.0_prec
        do i = 1, n
            matrix(i,i) = 1.0_prec
        end do
        reference = [-2.0_prec, 0.5_prec, 3.0_prec, 1.25_prec]
        rhs = matmul(matrix, reference)
        matrix_before = matrix
        rhs_before = rhs

        call GE_Expl_Main(n, matrix, rhs, solution)
        call assert_dense_solution(matrix, rhs, reference, solution, &
            "identity system", error)
        if (allocated(error)) return
        call assert_inputs_unchanged(matrix, matrix_before, rhs, rhs_before, &
            "identity system", error)
    end subroutine gauss_elimination_identity

    subroutine gauss_elimination_partial_pivoting(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 3
        real(prec) :: matrix(n,n), matrix_before(n,n)
        real(prec) :: rhs(n), rhs_before(n), reference(n), solution(n)

        matrix(1,:) = [0.0_prec,  2.0_prec, -1.0_prec]
        matrix(2,:) = [1.0_prec, -2.0_prec,  3.0_prec]
        matrix(3,:) = [2.0_prec,  1.0_prec,  1.0_prec]
        reference = [1.5_prec, -2.0_prec, 0.75_prec]
        rhs = matmul(matrix, reference)
        matrix_before = matrix
        rhs_before = rhs

        call GE_Expl_Main(n, matrix, rhs, solution)
        call assert_dense_solution(matrix, rhs, reference, solution, &
            "pivot-required system", error)
        if (allocated(error)) return
        call assert_inputs_unchanged(matrix, matrix_before, rhs, rhs_before, &
            "pivot-required system", error)
    end subroutine gauss_elimination_partial_pivoting

    subroutine gauss_elimination_nonsymmetric(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n = 5
        real(prec) :: matrix(n,n), matrix_before(n,n)
        real(prec) :: rhs(n), rhs_before(n), reference(n), solution(n)

        matrix(1,:) = [ 8.0_prec, -1.0_prec,  2.0_prec,  0.0_prec,  1.0_prec]
        matrix(2,:) = [ 1.0_prec,  7.0_prec, -1.0_prec,  2.0_prec,  0.0_prec]
        matrix(3,:) = [ 0.0_prec,  2.0_prec,  9.0_prec, -1.0_prec,  1.0_prec]
        matrix(4,:) = [ 2.0_prec,  0.0_prec,  1.0_prec,  8.0_prec, -2.0_prec]
        matrix(5,:) = [-1.0_prec,  1.0_prec,  0.0_prec,  2.0_prec,  7.0_prec]
        reference = [1.0_prec, -2.0_prec, 0.5_prec, 3.0_prec, -1.5_prec]
        rhs = matmul(matrix, reference)
        matrix_before = matrix
        rhs_before = rhs

        call GE_Expl_Main(n, matrix, rhs, solution)
        call assert_dense_solution(matrix, rhs, reference, solution, &
            "nonsymmetric system", error)
        if (allocated(error)) return
        call assert_inputs_unchanged(matrix, matrix_before, rhs, rhs_before, &
            "nonsymmetric system", error)
    end subroutine gauss_elimination_nonsymmetric

    subroutine assert_dense_solution(matrix, rhs, reference, solution, detail, error)
        real(prec), intent(in) :: matrix(:,:), rhs(:), reference(:), solution(:)
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        real(prec) :: residual(size(rhs)), residual_norm, rhs_norm, tolerance

        tolerance = test_tolerance()
        if (any(abs(solution - reference) > &
            tolerance*max(1.0_prec, abs(reference)))) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: GE solution mismatch for "//detail)
            return
        end if

        residual = matmul(matrix, solution) - rhs
        residual_norm = sqrt(dot_product(residual, residual))
        rhs_norm = max(1.0_prec, sqrt(dot_product(rhs, rhs)))
        if (residual_norm > tolerance*rhs_norm) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: GE residual too large for "//detail)
        end if
    end subroutine assert_dense_solution

    subroutine assert_inputs_unchanged(matrix, matrix_before, rhs, rhs_before, detail, error)
        real(prec), intent(in) :: matrix(:,:), matrix_before(:,:)
        real(prec), intent(in) :: rhs(:), rhs_before(:)
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        if (any(matrix /= matrix_before) .or. any(rhs /= rhs_before)) then
            call check(error, .false., message= &
                "[SYSTEM] TEST ERROR: GE modified inputs for "//detail)
        end if
    end subroutine assert_inputs_unchanged

    pure real(prec) function test_tolerance()
        test_tolerance = max(1.0e-12_prec, 5000.0_prec*epsilon(1.0_prec))
    end function test_tolerance

end module GaussElimination_suite
