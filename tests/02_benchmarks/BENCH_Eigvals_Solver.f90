program benchmark_eigenvalue_solver
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use precision_kinds, only : prec
    use BiCGSTAB_Solver, only : banded_operator_t, bicgstab_options_t
    use Eigen_spectrum_solver, only : dense_generalized_problem_t, &
        vectorized_generalized_problem_t, eigensolver_options_t, &
        eigensolver_workspace_t, eigensolver_result_t, solve_eigenproblem, &
        EIGEN_MODE_DOMINANT, EIGEN_SUCCESS

    implicit none

    integer, parameter :: problem_sizes(3) = [128, 256, 512]
    integer, parameter :: repetitions = 7
    integer, parameter :: requested_eigenvalues = 3
    integer, parameter :: krylov_dimension = 24
    real(prec), parameter :: benchmark_residual_tolerance = &
        max(1.0e-7_prec,1000.0_prec*epsilon(1.0_prec))
    integer :: case_index

    write(*,'(A,I0,A,I0,A,I0,A,ES13.6)') '# nev=', requested_eigenvalues, &
        ',ncv=', krylov_dimension, ',repetitions=', repetitions, &
        ',residual_tolerance=', benchmark_residual_tolerance
    write(*,'(A)') 'n,dense_setup_s,vector_setup_s,dense_first_solve_s,'// &
                   'vector_first_solve_s,dense_median_s,'// &
                   'vector_median_s,dense_over_vector_speedup,parity_error,dense_residual,'// &
                   'vector_residual,vector_inner_iterations'
    do case_index = 1, size(problem_sizes)
        call run_benchmark_case(problem_sizes(case_index))
    end do

contains

    subroutine run_benchmark_case(n)
        integer, intent(in) :: n

        real(prec), allocatable, target :: left_matrix(:,:), right_matrix(:,:)
        real(prec), allocatable, target :: lower(:), diagonal(:), upper(:)
        real(prec), allocatable, target :: zero_lower(:), right_diagonal(:), zero_upper(:)
        real(prec), allocatable :: initial_vector(:)
        real(prec) :: dense_times(repetitions), vector_times(repetitions)
        real(prec) :: dense_first_solve, dense_setup
        real(prec) :: vector_first_solve, vector_setup
        real(prec) :: eigenvalue_error, speedup
        type(banded_operator_t), target :: left_operator, right_operator
        type(dense_generalized_problem_t) :: dense_problem
        type(vectorized_generalized_problem_t) :: vector_problem
        type(eigensolver_options_t) :: options
        type(bicgstab_options_t) :: inner_options
        type(eigensolver_workspace_t) :: dense_workspace, vector_workspace
        type(eigensolver_result_t) :: dense_result, vector_result
        integer :: index, repetition, status
        integer(int64) :: clock_rate, end_count, start_count

        allocate(left_matrix(n,n), right_matrix(n,n), lower(n), diagonal(n), &
                 upper(n), zero_lower(n), right_diagonal(n), zero_upper(n), &
                 initial_vector(n))

        lower = -0.10_prec
        diagonal = 4.0_prec
        upper = -0.08_prec
        lower(1) = 0.0_prec
        upper(n) = 0.0_prec
        zero_lower = 0.0_prec
        zero_upper = 0.0_prec
        left_matrix = 0.0_prec
        right_matrix = 0.0_prec
        do index = 1, n
            select case (index)
            case (1)
                right_diagonal(index) = 40.0_prec
            case (2)
                right_diagonal(index) = 20.0_prec
            case (3)
                right_diagonal(index) = 10.0_prec
            case default
                right_diagonal(index) = 1.0_prec/real(index-2,prec)
            end select
            initial_vector(index) = 1.0_prec + real(index,prec)/real(n,prec)
            left_matrix(index,index) = diagonal(index)
            right_matrix(index,index) = right_diagonal(index)
        end do
        do index = 2, n
            left_matrix(index,index-1) = lower(index)
        end do
        do index = 1, n-1
            left_matrix(index,index+1) = upper(index)
        end do

        call system_clock(count_rate=clock_rate)
        if (clock_rate <= 0_int64) error stop 'benchmark clock is unavailable'

        call system_clock(start_count)
        call dense_problem%bind(left_matrix, right_matrix, status)
        call system_clock(end_count)
        dense_setup = elapsed_seconds(start_count, end_count, clock_rate)
        if (status /= EIGEN_SUCCESS) error stop 'dense benchmark binding failed'

        call system_clock(start_count)
        call left_operator%bind_tridiagonal(lower, diagonal, upper, status)
        if (status /= 0) error stop 'left banded benchmark binding failed'
        call right_operator%bind_tridiagonal(zero_lower, right_diagonal, zero_upper, status)
        if (status /= 0) error stop 'right banded benchmark binding failed'
        inner_options%max_iterations = 1000
        inner_options%relative_tolerance = max(1.0e-12_prec, &
                                                100.0_prec*epsilon(1.0_prec))
        inner_options%absolute_tolerance = 0.0_prec
        inner_options%breakdown_tolerance = 100.0_prec*epsilon(1.0_prec)
        call vector_problem%bind(left_operator, right_operator, inner_options, status)
        call system_clock(end_count)
        vector_setup = elapsed_seconds(start_count, end_count, clock_rate)
        if (status /= EIGEN_SUCCESS) error stop 'vector benchmark binding failed'

        options%mode = EIGEN_MODE_DOMINANT
        options%number_of_eigenvalues = requested_eigenvalues
        options%krylov_dimension = min(n,krylov_dimension)
        options%residual_tolerance = benchmark_residual_tolerance
        options%breakdown_tolerance = 1000.0_prec*epsilon(1.0_prec)

        ! First calls include lazy dense LU factorization and workspace setup.
        call system_clock(start_count)
        call solve_eigenproblem(dense_problem, options, dense_workspace, dense_result, &
                                initial_vector)
        call system_clock(end_count)
        dense_first_solve = elapsed_seconds(start_count, end_count, clock_rate)
        call system_clock(start_count)
        call solve_eigenproblem(vector_problem, options, vector_workspace, vector_result, &
                                initial_vector)
        call system_clock(end_count)
        vector_first_solve = elapsed_seconds(start_count, end_count, clock_rate)
        call require_success(dense_result, 'dense warm-up')
        call require_success(vector_result, 'vector warm-up')
        call require_equivalent_results(dense_result, vector_result, options)
        call require_independent_residuals(left_matrix, right_matrix, dense_result, &
                                           options%residual_tolerance)
        call require_independent_residuals(left_matrix, right_matrix, vector_result, &
                                           options%residual_tolerance)

        do repetition = 1, repetitions
            ! Alternate execution order to reduce systematic cache and
            ! frequency-scaling bias between the two implementations.
            if (mod(repetition,2) == 1) then
                call system_clock(start_count)
                call solve_eigenproblem(dense_problem, options, dense_workspace, &
                                        dense_result, initial_vector)
                call system_clock(end_count)
                dense_times(repetition) = &
                    elapsed_seconds(start_count, end_count, clock_rate)
                call require_success(dense_result, 'dense timed solve')

                call system_clock(start_count)
                call solve_eigenproblem(vector_problem, options, vector_workspace, &
                                        vector_result, initial_vector)
                call system_clock(end_count)
                vector_times(repetition) = &
                    elapsed_seconds(start_count, end_count, clock_rate)
                call require_success(vector_result, 'vector timed solve')
            else
                call system_clock(start_count)
                call solve_eigenproblem(vector_problem, options, vector_workspace, &
                                        vector_result, initial_vector)
                call system_clock(end_count)
                vector_times(repetition) = &
                    elapsed_seconds(start_count, end_count, clock_rate)
                call require_success(vector_result, 'vector timed solve')

                call system_clock(start_count)
                call solve_eigenproblem(dense_problem, options, dense_workspace, &
                                        dense_result, initial_vector)
                call system_clock(end_count)
                dense_times(repetition) = &
                    elapsed_seconds(start_count, end_count, clock_rate)
                call require_success(dense_result, 'dense timed solve')
            end if
        end do

        call sort_times(dense_times)
        call sort_times(vector_times)
        call require_equivalent_results(dense_result, vector_result, options)
        call require_independent_residuals(left_matrix, right_matrix, dense_result, &
                                           options%residual_tolerance)
        call require_independent_residuals(left_matrix, right_matrix, vector_result, &
                                           options%residual_tolerance)
        eigenvalue_error = maxval(abs(dense_result%eigenvalues - &
                                      vector_result%eigenvalues))
        if (vector_times((repetitions+1)/2) <= 0.0_prec) &
            error stop 'benchmark timer resolution is insufficient'
        speedup = dense_times((repetitions+1)/2)/vector_times((repetitions+1)/2)

        write(*,'(I0,10(",",ES13.6),",",I0)') n, dense_setup, vector_setup, &
            dense_first_solve, vector_first_solve, dense_times((repetitions+1)/2), &
            vector_times((repetitions+1)/2), &
            speedup, eigenvalue_error, maxval(dense_result%relative_residuals), &
            maxval(vector_result%relative_residuals), vector_result%inner_iterations
    end subroutine run_benchmark_case

    subroutine require_success(result, stage)
        type(eigensolver_result_t), intent(in) :: result
        character(len=*), intent(in) :: stage

        integer :: index

        if (result%status /= EIGEN_SUCCESS .or. .not. result%converged) then
            write(*,'(A)') trim(stage)//': '//trim(result%reason)
            error stop 'eigenvalue benchmark did not converge'
        end if
        if (.not. allocated(result%eigenvalues) .or. &
            .not. allocated(result%relative_residuals)) &
            error stop 'eigenvalue benchmark result is incomplete'
        if (result%returned_eigenvalues <= 0 .or. &
            result%returned_eigenvalues /= requested_eigenvalues .or. &
            result%converged_eigenvalues /= result%returned_eigenvalues .or. &
            size(result%eigenvalues) /= result%returned_eigenvalues .or. &
            size(result%relative_residuals) /= result%returned_eigenvalues) &
            error stop 'eigenvalue benchmark result counts are inconsistent'
        if (.not. all(ieee_is_finite(real(result%eigenvalues,kind=prec))) .or. &
            .not. all(ieee_is_finite(aimag(result%eigenvalues))) .or. &
            .not. all(ieee_is_finite(result%relative_residuals))) &
            error stop 'eigenvalue benchmark result is nonfinite'
        if (.not. allocated(result%eigenvectors) .or. &
            .not. allocated(result%absolute_residuals) .or. &
            .not. allocated(result%ritz_residual_estimates)) &
            error stop 'eigenvalue benchmark diagnostics are incomplete'
        if (size(result%eigenvectors,2) /= result%returned_eigenvalues .or. &
            .not. all(ieee_is_finite(real(result%eigenvectors,kind=prec))) .or. &
            .not. all(ieee_is_finite(aimag(result%eigenvectors))) .or. &
            .not. all(ieee_is_finite(result%absolute_residuals)) .or. &
            .not. all(ieee_is_finite(result%ritz_residual_estimates)) .or. &
            any(result%absolute_residuals < 0.0_prec) .or. &
            any(result%relative_residuals < 0.0_prec) .or. &
            any(result%ritz_residual_estimates < 0.0_prec)) &
            error stop 'eigenvalue benchmark diagnostics are invalid'
        if (any([(complex_vector_norm(result%eigenvectors(:,index)) <= 0.0_prec, &
                  index=1,size(result%eigenvectors,2))])) &
            error stop 'eigenvalue benchmark returned a zero eigenvector'
    end subroutine require_success

    subroutine require_independent_residuals(left_matrix, right_matrix, result, &
                                             tolerance)
        real(prec), intent(in) :: left_matrix(:,:), right_matrix(:,:)
        type(eigensolver_result_t), intent(in) :: result
        real(prec), intent(in) :: tolerance

        complex(prec) :: left_value(size(left_matrix,1))
        complex(prec) :: residual(size(left_matrix,1))
        complex(prec) :: right_value(size(left_matrix,1))
        real(prec) :: denominator, relative_residual
        integer :: index

        if (size(left_matrix,2) /= size(left_matrix,1) .or. &
            size(right_matrix,1) /= size(left_matrix,1) .or. &
            size(right_matrix,2) /= size(left_matrix,1) .or. &
            size(result%eigenvectors,1) /= size(left_matrix,1)) &
            error stop 'independent benchmark dimensions are inconsistent'
        do index = 1, result%returned_eigenvalues
            left_value = matmul(cmplx(left_matrix,0.0_prec,kind=prec), &
                                result%eigenvectors(:,index))
            right_value = matmul(cmplx(right_matrix,0.0_prec,kind=prec), &
                                 result%eigenvectors(:,index))
            residual = right_value-result%eigenvalues(index)*left_value
            denominator = complex_vector_norm(right_value) + &
                abs(result%eigenvalues(index))*complex_vector_norm(left_value)
            relative_residual = complex_vector_norm(residual)/ &
                max(tiny(1.0_prec),denominator)
            if (.not. ieee_is_finite(relative_residual) .or. &
                relative_residual > tolerance) &
                error stop 'independent benchmark residual threshold was not met'
        end do
    end subroutine require_independent_residuals

    pure real(prec) function complex_vector_norm(vector)
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
            complex_vector_norm = 0.0_prec
        else
            complex_vector_norm = norm_scale*sqrt(sum_of_squares)
        end if
    end function complex_vector_norm

    subroutine require_equivalent_results(dense_result, vector_result, options)
        type(eigensolver_result_t), intent(in) :: dense_result, vector_result
        type(eigensolver_options_t), intent(in) :: options

        real(prec) :: normalized_parity, parity_tolerance

        if (size(dense_result%eigenvalues) /= size(vector_result%eigenvalues)) &
            error stop 'benchmark result sizes differ'
        if (dense_result%returned_eigenvalues /= vector_result%returned_eigenvalues .or. &
            dense_result%converged_eigenvalues /= &
                vector_result%converged_eigenvalues .or. &
            dense_result%arnoldi_steps /= vector_result%arnoldi_steps) &
            error stop 'benchmark methods performed unequal Arnoldi work'
        parity_tolerance = max(1.0e-9_prec, &
                               10.0_prec*options%residual_tolerance)
        normalized_parity = maxval(abs(dense_result%eigenvalues - &
            vector_result%eigenvalues)/max(1.0_prec, &
                                           abs(dense_result%eigenvalues)))
        if (.not. ieee_is_finite(normalized_parity) .or. &
            normalized_parity > parity_tolerance) &
            error stop 'dense/vectorized eigenvalues differ'
        if (maxval(dense_result%relative_residuals) > &
            options%residual_tolerance .or. &
            maxval(vector_result%relative_residuals) > &
            options%residual_tolerance) &
            error stop 'benchmark residual threshold was not met'
    end subroutine require_equivalent_results

    pure real(prec) function elapsed_seconds(start_count, end_count, clock_rate)
        integer(int64), intent(in) :: start_count, end_count, clock_rate

        elapsed_seconds = real(end_count-start_count,prec)/real(clock_rate,prec)
    end function elapsed_seconds

    subroutine sort_times(values)
        real(prec), intent(inout) :: values(:)

        real(prec) :: current
        integer :: i, position

        do i = 2, size(values)
            current = values(i)
            position = i - 1
            do while (position >= 1)
                if (values(position) <= current) exit
                values(position+1) = values(position)
                position = position - 1
            end do
            values(position+1) = current
        end do
    end subroutine sort_times

end program benchmark_eigenvalue_solver
