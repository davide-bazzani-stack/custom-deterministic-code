module BiCGSTAB_suite
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use precision_kinds, only : prec

    implicit none

    private

    public :: collect_suite_BiCGSTAB

contains

    subroutine collect_suite_BiCGSTAB(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("classic_pentadiagonal_known_solution", &
                classic_pentadiagonal_known_solution), &
            new_unittest("triiso_migration_known_solution", &
                triiso_migration_known_solution) &
            ]
    end subroutine collect_suite_BiCGSTAB

    subroutine classic_pentadiagonal_known_solution(error)
        use BiCGSTAB_Algorithm, only : BiCGSTAB_Classic

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 3
        integer, parameter :: n = 9
        integer, parameter :: max_iterations = 200
        real(prec), parameter :: solver_tolerance = &
            max(1.0e-12_prec, 100.0_prec*epsilon(1.0_prec))
        real(prec), parameter :: check_tolerance = &
            max(1.0e-10_prec, 1000.0_prec*epsilon(1.0_prec))
        real(prec), parameter :: x_reference(n) = [ &
             1.00_prec, -2.00_prec,  0.50_prec, &
             3.00_prec, -1.50_prec,  2.25_prec, &
            -0.75_prec,  1.25_prec, -3.00_prec  &
            ]

        real(prec) :: a(n), b(n), c(n), d(n), e(n)
        real(prec) :: reference_matrix(n, n)
        real(prec) :: rhs(n), x(n), residual(n)
        integer :: info

        ! Flattened 3 x 3 grid with bands at -n_x, -1, 0, +1, +n_x.
        ! Row-boundary coefficients are zero, preventing wraparound coupling.
        a = [ &
             0.0_prec, -0.8_prec, -0.9_prec, &
             0.0_prec, -0.7_prec, -0.6_prec, &
             0.0_prec, -1.0_prec, -0.8_prec  &
            ]
        b = [ &
             6.0_prec, 6.5_prec, 7.0_prec, &
             6.2_prec, 6.8_prec, 7.2_prec, &
             6.4_prec, 6.9_prec, 7.1_prec  &
            ]
        c = [ &
            -1.0_prec, -1.1_prec,  0.0_prec, &
            -1.2_prec, -1.0_prec,  0.0_prec, &
            -0.9_prec, -1.1_prec,  0.0_prec  &
            ]
        d = [ &
             0.0_prec,  0.0_prec,  0.0_prec, &
            -0.4_prec, -0.5_prec, -0.6_prec, &
            -0.7_prec, -0.8_prec, -0.9_prec  &
            ]
        e = [ &
            -0.5_prec, -0.6_prec, -0.7_prec, &
            -0.8_prec, -0.9_prec, -1.0_prec, &
             0.0_prec,  0.0_prec,  0.0_prec  &
            ]

        reference_matrix = 0.0_prec
        reference_matrix(1, [1, 2, 4]) = [6.0_prec, -1.0_prec, -0.5_prec]
        reference_matrix(2, [1, 2, 3, 5]) = [ &
            -0.8_prec, 6.5_prec, -1.1_prec, -0.6_prec]
        reference_matrix(3, [2, 3, 6]) = [-0.9_prec, 7.0_prec, -0.7_prec]
        reference_matrix(4, [1, 4, 5, 7]) = [ &
            -0.4_prec, 6.2_prec, -1.2_prec, -0.8_prec]
        reference_matrix(5, [2, 4, 5, 6, 8]) = [ &
            -0.5_prec, -0.7_prec, 6.8_prec, -1.0_prec, -0.9_prec]
        reference_matrix(6, [3, 5, 6, 9]) = [ &
            -0.6_prec, -0.6_prec, 7.2_prec, -1.0_prec]
        reference_matrix(7, [4, 7, 8]) = [-0.7_prec, 6.4_prec, -0.9_prec]
        reference_matrix(8, [5, 7, 8, 9]) = [ &
            -0.8_prec, -1.0_prec, 6.9_prec, -1.1_prec]
        reference_matrix(9, [6, 8, 9]) = [-0.9_prec, -0.8_prec, 7.1_prec]

        rhs = matmul(reference_matrix, x_reference)
        x = 0.0_prec

        call BiCGSTAB_Classic( &
            n_x, n, a, b, c, d, e, rhs, x, &
            solver_tolerance, max_iterations, info)

        call check(error, info, 0, &
            message="[SYSTEM] TEST ERROR: BiCGSTAB_Classic did not converge")
        if (allocated(error)) return

        call assert_solution(x, x_reference, check_tolerance, &
            "BiCGSTAB_Classic solution mismatch", error)
        if (allocated(error)) return

        residual = matmul(reference_matrix, x) - rhs
        call assert_residual(residual, rhs, check_tolerance, &
            "BiCGSTAB_Classic residual too large", error)
        
    end subroutine classic_pentadiagonal_known_solution

    subroutine triiso_migration_known_solution(error)
        use BiCGSTAB_Algorithm, only : BiCGSTAB_TriIso

        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: n_x = 4
        integer, parameter :: n = 8
        integer, parameter :: max_iterations = 100
        real(prec), parameter :: solver_tolerance = &
            max(1.0e-12_prec, 100.0_prec*epsilon(1.0_prec))
        real(prec), parameter :: check_tolerance = &
            max(1.0e-10_prec, 1000.0_prec*epsilon(1.0_prec))
        real(prec), parameter :: a(n) = [ &
             0.0_prec,  0.0_prec, -0.8_prec,  0.0_prec, &
             0.0_prec, -0.7_prec,  0.0_prec, -0.9_prec ]
        real(prec), parameter :: b(n) = [ &
             4.5_prec, 5.0_prec, 4.8_prec, 5.2_prec, &
             4.6_prec, 5.1_prec, 4.9_prec, 5.3_prec ]
        real(prec), parameter :: c(n) = [ &
             0.0_prec, -0.8_prec, 0.0_prec, 0.0_prec, &
            -0.7_prec,  0.0_prec, -0.9_prec, 0.0_prec ]
        real(prec), parameter :: d(n) = [ &
             0.0_prec,  0.0_prec,  0.0_prec,  0.0_prec, &
            -1.1_prec, -1.2_prec, -1.0_prec, -1.3_prec ]
        real(prec), parameter :: e(n) = [ &
            -1.1_prec, -1.2_prec, -1.0_prec, -1.3_prec, &
             0.0_prec,  0.0_prec,  0.0_prec,  0.0_prec ]
        real(prec), parameter :: x_reference(n) = [ &
            0.75_prec, 1.20_prec, 0.90_prec, 1.30_prec, &
            1.10_prec, 0.90_prec, 1.25_prec, 0.95_prec ]

        real(prec) :: reference_matrix(n, n)
        real(prec) :: rhs(n), x(n), residual(n)
        integer :: info

        ! Two packed rows with BL/BR/BL/BR and TR/TL/TR/TL orientations.
        ! The explicit dense matrix is independent of the production matvec.
        reference_matrix = 0.0_prec
        reference_matrix(1, 1) = b(1)
        reference_matrix(2, 2) = b(2)
        reference_matrix(3, 3) = b(3)
        reference_matrix(4, 4) = b(4)
        reference_matrix(5, 5) = b(5)
        reference_matrix(6, 6) = b(6)
        reference_matrix(7, 7) = b(7)
        reference_matrix(8, 8) = b(8)

        reference_matrix(1, 5) = e(1)
        reference_matrix(2, 6) = e(2)
        reference_matrix(3, 7) = e(3)
        reference_matrix(4, 8) = e(4)
        reference_matrix(5, 1) = d(5)
        reference_matrix(6, 2) = d(6)
        reference_matrix(7, 3) = d(7)
        reference_matrix(8, 4) = d(8)

        reference_matrix(2, 3) = c(2)
        reference_matrix(3, 2) = a(3)
        reference_matrix(5, 6) = c(5)
        reference_matrix(6, 5) = a(6)
        reference_matrix(7, 8) = c(7)
        reference_matrix(8, 7) = a(8)

        rhs = matmul(reference_matrix, x_reference)
        x = 0.0_prec

        call BiCGSTAB_TriIso( &
            n_x, n, a, b, c, d, e, rhs, x, &
            solver_tolerance, max_iterations, info)

        call check(error, info, 0, &
            message="[SYSTEM] TEST ERROR: BiCGSTAB_TriIso did not converge")
        if (allocated(error)) return

        call assert_solution(x, x_reference, check_tolerance, &
            "BiCGSTAB_TriIso solution mismatch", error)
        if (allocated(error)) return

        residual = matmul(reference_matrix, x) - rhs
        call assert_residual(residual, rhs, check_tolerance, &
            "BiCGSTAB_TriIso residual too large", error)
    end subroutine triiso_migration_known_solution

    subroutine assert_solution(input, reference, tolerance, detail, error)
        real(prec), intent(in) :: input(:), reference(:), tolerance
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        integer :: i

        do i = 1, size(reference)
            if (abs(input(i) - reference(i)) > &
                tolerance*max(1.0_prec, abs(reference(i)))) then
                call check(error, .false., &
                    message="[SYSTEM] TEST ERROR: "//detail)
                return
            end if
        end do
    end subroutine assert_solution

    subroutine assert_residual(residual, rhs, tolerance, detail, error)
        real(prec), intent(in) :: residual(:), rhs(:), tolerance
        character(len=*), intent(in) :: detail
        type(error_type), allocatable, intent(out) :: error

        real(prec) :: residual_norm, rhs_norm

        residual_norm = sqrt(dot_product(residual, residual))
        rhs_norm = max(1.0_prec, sqrt(dot_product(rhs, rhs)))

        if (residual_norm > tolerance*rhs_norm) then
            call check(error, .false., &
                message="[SYSTEM] TEST ERROR: "//detail)
        end if
    end subroutine assert_residual

end module BiCGSTAB_suite
