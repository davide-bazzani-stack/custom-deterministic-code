module Eigen_spectrum_solver
    ! Generalized real nonsymmetric eigenproblems are written as
    !
    !                    F*x = lambda*A*x.
    !
    ! Full spectra are computed from explicit dense matrices with LAPACK
    ! xGGEV.  Dominant eigenpairs are computed with an Arnoldi projection of
    ! K=A^{-1}F.  Dense and vectorized problems implement the same operator
    ! contract, so the Arnoldi iteration is independent of matrix storage.
    ! The LAPACK backend supports prec=real32 and prec=real64.  A real128
    ! selection requires a separate quadruple-precision backend.
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, &
        ieee_value
    use LAPACK_Interfaces, only : lapack_getrf, lapack_getrs, lapack_gecon, &
        lapack_geev, lapack_ggev, LAPACK_WRAPPER_ALLOCATION_FAILURE
    
    use precision_kinds, only : prec
    use BiCGSTAB_Solver, only : linear_operator_t, preconditioner_t, &
        bicgstab_options_t, bicgstab_result_t, bicgstab_workspace_t, &
        solve_bicgstab, BICGSTAB_SUCCESS, BICGSTAB_COMPONENT_SUCCESS

    implicit none

    private

    integer, parameter, public :: EIGEN_MODE_DOMINANT = 1
    integer, parameter, public :: EIGEN_MODE_FULL = 2

    integer, parameter, public :: EIGEN_SUCCESS = 0
    integer, parameter, public :: EIGEN_INVALID_INPUT = 1
    integer, parameter, public :: EIGEN_ALLOCATION_FAILURE = 2
    integer, parameter, public :: EIGEN_ZERO_START_VECTOR = 3
    integer, parameter, public :: EIGEN_FACTORIZATION_FAILURE = 4
    integer, parameter, public :: EIGEN_LAPACK_FAILURE = 5
    integer, parameter, public :: EIGEN_OPERATOR_FAILURE = 6
    integer, parameter, public :: EIGEN_INNER_SOLVER_FAILURE = 7
    integer, parameter, public :: EIGEN_NONFINITE_VALUE = 8
    integer, parameter, public :: EIGEN_NOT_CONVERGED = 9
    integer, parameter, public :: EIGEN_INFINITE_EIGENVALUE = 10
    integer, parameter, public :: EIGEN_SINGULAR_PENCIL = 11

    type, public :: eigensolver_options_t
        integer :: mode = EIGEN_MODE_DOMINANT
        integer :: number_of_eigenvalues = 1
        ! Zero selects the full available dimension.  A smaller explicit
        ! value limits dominant-mode storage and work to one Arnoldi cycle.
        ! Convergence is never claimed unless the true generalized residual
        ! passes; a future Krylov-Schur layer can add bounded restarts without
        ! changing the operator contract.
        integer :: krylov_dimension = 0
        real(prec) :: residual_tolerance = sqrt(epsilon(1.0_prec))
        real(prec) :: breakdown_tolerance = 100.0_prec*epsilon(1.0_prec)
    end type eigensolver_options_t

    type, public :: eigensolver_result_t
        logical :: converged = .false.
        integer :: status = EIGEN_INVALID_INPUT
        integer :: requested_eigenvalues = 0
        integer :: returned_eigenvalues = 0
        integer :: converged_eigenvalues = 0
        integer :: arnoldi_steps = 0
        integer :: arnoldi_breakdowns = 0
        integer :: operator_applications = 0
        integer :: linear_solves = 0
        integer :: inner_iterations = 0
        integer :: last_inner_solver_status = -1
        real(prec) :: reciprocal_condition_estimate = -1.0_prec
        character(len=24) :: method = 'not started'
        character(len=128) :: reason = 'solver not started'
        character(len=96) :: last_inner_solver_reason = 'inner solver not used'
        complex(prec), allocatable :: eigenvalues(:)
        complex(prec), allocatable :: eigenvectors(:,:)
        real(prec), allocatable :: absolute_residuals(:)
        real(prec), allocatable :: relative_residuals(:)
        real(prec), allocatable :: ritz_residual_estimates(:)
    end type eigensolver_result_t

    type, public :: eigensolver_workspace_t
        private
        real(prec), allocatable :: basis(:,:)
        real(prec), allocatable :: hessenberg(:,:)
        real(prec), allocatable :: work_vector(:)
        real(prec), allocatable :: right_hand_side(:)
        real(prec), allocatable :: residual_left_real(:)
        real(prec), allocatable :: residual_left_imaginary(:)
        real(prec), allocatable :: residual_right_real(:)
        real(prec), allocatable :: residual_right_imaginary(:)
        contains
            procedure, public :: clear => clear_eigensolver_workspace
            procedure, private :: ensure_capacity => ensure_eigensolver_workspace
    end type eigensolver_workspace_t

    type, abstract, public :: generalized_eigenproblem_t
        contains
            procedure(problem_order_interface), deferred, public :: order
            procedure(problem_valid_interface), deferred, public :: is_valid
            procedure(problem_apply_interface), deferred, public :: apply_left
            procedure(problem_apply_interface), deferred, public :: apply_right
            procedure(problem_solve_interface), deferred, public :: solve_left
            procedure(problem_reset_statistics_interface), deferred, public :: reset_statistics
            procedure(problem_get_statistics_interface), deferred, public :: get_statistics
    end type generalized_eigenproblem_t

    type, extends(generalized_eigenproblem_t), public :: dense_generalized_problem_t
        ! The matrices are copied at bind time.  This makes the cached lazy LU
        ! factorization immune to caller mutation and temporary lifetimes.
        ! One problem object is still required per concurrent solve because
        ! the reusable solve buffer and statistics are mutable.
        private
        integer :: n = 0
        logical :: ready = .false.
        logical :: factorization_ready = .false.
        real(prec), allocatable :: left_matrix(:,:)
        real(prec), allocatable :: right_matrix(:,:)
        real(prec), allocatable :: left_lu(:,:)
        real(prec), allocatable :: solve_buffer(:,:)
        integer, allocatable :: pivots(:)
        integer :: linear_solve_count = 0
        real(prec) :: factorization_scale = 1.0_prec
        real(prec) :: reciprocal_condition = -1.0_prec
        contains
            procedure, public :: bind => bind_dense_problem
            procedure, public :: clear => clear_dense_problem
            procedure, public :: condition_estimate => dense_condition_estimate
            procedure, public :: order => dense_problem_order
            procedure, public :: is_valid => dense_problem_is_valid
            procedure, public :: apply_left => apply_dense_left
            procedure, public :: apply_right => apply_dense_right
            procedure, public :: solve_left => solve_dense_left
            procedure, public :: reset_statistics => reset_dense_statistics
            procedure, public :: get_statistics => get_dense_statistics
    end type dense_generalized_problem_t

    type, extends(generalized_eigenproblem_t), public :: vectorized_generalized_problem_t
        ! Operator/preconditioner objects are copied at bind time.  Any
        ! coefficient/context pointers held inside those objects remain
        ! non-owning and must outlive the problem.  One problem object is
        ! required per concurrent Arnoldi solve.
        private
        integer :: n = 0
        logical :: ready = .false.
        class(linear_operator_t), allocatable :: left_operator
        class(linear_operator_t), allocatable :: right_operator
        class(preconditioner_t), allocatable :: preconditioner
        type(bicgstab_options_t) :: inner_options
        type(bicgstab_workspace_t) :: inner_workspace
        type(bicgstab_result_t) :: last_inner_result
        integer :: linear_solve_count = 0
        integer :: total_inner_iterations = 0
        contains
            procedure, public :: bind => bind_vectorized_problem
            procedure, public :: clear => clear_vectorized_problem
            procedure, public :: order => vectorized_problem_order
            procedure, public :: is_valid => vectorized_problem_is_valid
            procedure, public :: apply_left => apply_vectorized_left
            procedure, public :: apply_right => apply_vectorized_right
            procedure, public :: solve_left => solve_vectorized_left
            procedure, public :: reset_statistics => reset_vectorized_statistics
            procedure, public :: get_statistics => get_vectorized_statistics
    end type vectorized_generalized_problem_t

    abstract interface
        integer function problem_order_interface(this)
            import :: generalized_eigenproblem_t
            class(generalized_eigenproblem_t), intent(in) :: this
        end function problem_order_interface

        logical function problem_valid_interface(this)
            import :: generalized_eigenproblem_t
            class(generalized_eigenproblem_t), intent(in) :: this
        end function problem_valid_interface

        subroutine problem_apply_interface(this, input, output, status)
            import :: generalized_eigenproblem_t, prec
            class(generalized_eigenproblem_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status
        end subroutine problem_apply_interface

        subroutine problem_solve_interface(this, right_hand_side, solution, status)
            import :: generalized_eigenproblem_t, prec
            class(generalized_eigenproblem_t), intent(inout) :: this
            real(prec), intent(in) :: right_hand_side(:)
            real(prec), intent(out) :: solution(:)
            integer, intent(out) :: status
        end subroutine problem_solve_interface

        subroutine problem_reset_statistics_interface(this)
            import :: generalized_eigenproblem_t
            class(generalized_eigenproblem_t), intent(inout) :: this
        end subroutine problem_reset_statistics_interface

        subroutine problem_get_statistics_interface(this, linear_solves, inner_iterations)
            import :: generalized_eigenproblem_t
            class(generalized_eigenproblem_t), intent(in) :: this
            integer, intent(out) :: linear_solves, inner_iterations
        end subroutine problem_get_statistics_interface
    end interface

    public :: solve_eigenproblem

    contains

        ! Solver entry point

        subroutine solve_eigenproblem(problem, options, workspace, result, initial_vector)
            class(generalized_eigenproblem_t), intent(inout) :: problem
            type(eigensolver_options_t), intent(in) :: options
            type(eigensolver_workspace_t), intent(inout) :: workspace
            type(eigensolver_result_t), intent(out) :: result
            real(prec), optional, intent(in) :: initial_vector(:)

            integer :: n

            call initialize_result(result)

            n = problem%order()
            if (.not. common_inputs_are_valid(problem, options, n)) then
                call set_failure(result, EIGEN_INVALID_INPUT, &
                                 'invalid eigenproblem dimensions or options')
                return
            end if

            call problem%reset_statistics()

            select case (options%mode)
            case (EIGEN_MODE_FULL)
                call solve_full_spectrum(problem, options, workspace, result)

            case (EIGEN_MODE_DOMINANT)
                call solve_dominant_arnoldi(problem, options, workspace, result, &
                                            initial_vector)

            case default
                call set_failure(result, EIGEN_INVALID_INPUT, &
                                 'unknown eigenvalue solution mode')
            end select

            call capture_problem_statistics(problem, result)
        end subroutine solve_eigenproblem


        
        ! Full dense generalized spectrum

        subroutine solve_full_spectrum(problem, options, workspace, result)
            class(generalized_eigenproblem_t), intent(inout) :: problem
            type(eigensolver_options_t), intent(in) :: options
            type(eigensolver_workspace_t), intent(inout) :: workspace
            type(eigensolver_result_t), intent(inout) :: result

            real(prec), allocatable :: pencil_left(:,:), pencil_right(:,:)
            real(prec), allocatable :: alpha_real(:), alpha_imaginary(:), beta(:)
            real(prec), allocatable :: lapack_vectors(:,:)
            complex(prec), allocatable :: eigenvalues(:), eigenvectors(:,:)
            integer :: allocation_status, decode_status, index, info, n
            integer :: residual_applications

            n = problem%order()
            result%requested_eigenvalues = n

            call workspace%ensure_capacity(n, 1, allocation_status)
            if (allocation_status /= 0) then
                call set_failure(result, EIGEN_ALLOCATION_FAILURE, &
                                 'unable to allocate residual workspace')
                return
            end if

            allocate(pencil_left(n,n), pencil_right(n,n), alpha_real(n), &
                     alpha_imaginary(n), beta(n), lapack_vectors(n,n), &
                     stat=allocation_status)
            if (allocation_status /= 0) then
                call set_failure(result, EIGEN_ALLOCATION_FAILURE, &
                                 'unable to allocate full-spectrum workspace')
                return
            end if

            ! LAPACK solves pencil_left*x=lambda*pencil_right*x.  The physical
            ! convention is F*x=lambda*A*x.  Dense matrices can be copied
            ! directly; a vectorized problem is explicitly materialized one
            ! basis column at a time because a complete spectrum necessarily
            ! requires global matrix information.
            select type (dense_problem => problem)
            type is (dense_generalized_problem_t)
                result%method = 'dense LAPACK xGGEV'
                pencil_left = dense_problem%right_matrix
                pencil_right = dense_problem%left_matrix
            class default
                result%method = 'materialized LAPACK'
                do index = 1, n
                    workspace%work_vector = 0.0_prec
                    workspace%work_vector(index) = 1.0_prec
                    call problem%apply_right(workspace%work_vector, &
                                             pencil_left(:,index), decode_status)
                    if (decode_status /= EIGEN_SUCCESS) then
                        call set_failure(result, decode_status, &
                                         'unable to materialize right operator')
                        return
                    end if
                    result%operator_applications = result%operator_applications + 1
                    call problem%apply_left(workspace%work_vector, &
                                            pencil_right(:,index), decode_status)
                    if (decode_status /= EIGEN_SUCCESS) then
                        call set_failure(result, decode_status, &
                                         'unable to materialize left operator')
                        return
                    end if
                    result%operator_applications = result%operator_applications + 1
                end do
            end select
            call lapack_ggev(pencil_left, pencil_right, alpha_real, &
                             alpha_imaginary, beta, lapack_vectors, info)
            if (info /= 0) then
                if (info == LAPACK_WRAPPER_ALLOCATION_FAILURE) then
                    call set_failure(result, EIGEN_ALLOCATION_FAILURE, &
                                     'LAPACK workspace allocation failed')
                else
                    call set_failure(result, EIGEN_LAPACK_FAILURE, &
                                     'LAPACK generalized eigensolver failed')
                end if
                return
            end if

            call decode_generalized_eigensystem(alpha_real, alpha_imaginary, beta, &
                                                 lapack_vectors, eigenvalues, &
                                                 eigenvectors, decode_status)
            if (decode_status /= EIGEN_SUCCESS) then
                call set_failure(result, decode_status, &
                                 'invalid or singular generalized eigenvalue')
                return
            end if

            call sort_eigenpairs(eigenvalues, eigenvectors, decode_status)
            if (decode_status /= EIGEN_SUCCESS) then
                call set_failure(result, decode_status, &
                                 'unable to sort full-spectrum eigenpairs')
                return
            end if
            call allocate_result_diagnostics(result, n, allocation_status)
            if (allocation_status /= 0) then
                call set_failure(result, EIGEN_ALLOCATION_FAILURE, &
                                 'unable to allocate full-spectrum result')
                return
            end if

            call move_alloc(eigenvalues, result%eigenvalues)
            call move_alloc(eigenvectors, result%eigenvectors)
            result%ritz_residual_estimates = 0.0_prec

            do index = 1, n
                call normalize_and_fix_phase(result%eigenvectors(:,index), decode_status)
                if (decode_status /= EIGEN_SUCCESS) then
                    call set_failure(result, decode_status, &
                                     'LAPACK returned an invalid eigenvector')
                    return
                end if
                call evaluate_generalized_residual(problem, workspace, &
                    result%eigenvalues(index), result%eigenvectors(:,index), &
                    result%absolute_residuals(index), &
                    result%relative_residuals(index), decode_status, &
                    residual_applications)
                result%operator_applications = result%operator_applications + &
                                               residual_applications
                if (decode_status /= EIGEN_SUCCESS) then
                    call set_failure(result, decode_status, &
                                     'unable to verify a generalized eigenpair')
                    return
                end if
            end do

            result%returned_eigenvalues = n
            result%converged_eigenvalues = count( &
                result%relative_residuals <= options%residual_tolerance)
            if (result%converged_eigenvalues == n) then
                call set_success(result, 'full spectrum converged')
            else
                call set_failure(result, EIGEN_NOT_CONVERGED, &
                                 'full-spectrum residual verification failed')
            end if
        end subroutine solve_full_spectrum



        ! Dominant Arnoldi spectrum

        subroutine solve_dominant_arnoldi(problem, options, workspace, result, &
                                          initial_vector)
            class(generalized_eigenproblem_t), intent(inout) :: problem
            type(eigensolver_options_t), intent(in) :: options
            type(eigensolver_workspace_t), intent(inout) :: workspace
            type(eigensolver_result_t), intent(inout) :: result
            real(prec), optional, intent(in) :: initial_vector(:)

            complex(prec), allocatable :: projected_values(:)
            complex(prec), allocatable :: projected_vectors(:,:)
            real(prec) :: coefficient, next_norm, operator_norm, start_norm, start_scale
            integer :: allocation_status, actual_dimension, component_status
            integer :: i, j, n, number_to_return, pass, requested_dimension
            integer :: residual_applications
            logical :: complement_found, expansion_breakdown

            n = problem%order()
            requested_dimension = options%krylov_dimension
            if (requested_dimension == 0) requested_dimension = n

            result%method = 'Arnoldi projection'
            result%requested_eigenvalues = options%number_of_eigenvalues
            select type (dense_problem => problem)
            type is (dense_generalized_problem_t)
                result%reciprocal_condition_estimate = &
                    dense_problem%reciprocal_condition
            class default
                result%reciprocal_condition_estimate = -1.0_prec
            end select

            if (options%number_of_eigenvalues < 1 .or. &
                options%number_of_eigenvalues > n .or. &
                requested_dimension < options%number_of_eigenvalues .or. &
                requested_dimension > n) then
                call set_failure(result, EIGEN_INVALID_INPUT, &
                                 'invalid dominant count or Krylov dimension')
                return
            end if

            call workspace%ensure_capacity(n, requested_dimension, allocation_status)
            if (allocation_status /= 0) then
                call set_failure(result, EIGEN_ALLOCATION_FAILURE, &
                                 'unable to allocate Arnoldi workspace')
                return
            end if

            workspace%basis = 0.0_prec
            workspace%hessenberg = 0.0_prec
            workspace%work_vector = 0.0_prec
            workspace%right_hand_side = 0.0_prec

            if (present(initial_vector)) then
                if (size(initial_vector) /= n .or. &
                    .not. all(ieee_is_finite(initial_vector))) then
                    call set_failure(result, EIGEN_INVALID_INPUT, &
                                     'invalid Arnoldi initial vector')
                    return
                end if
                workspace%basis(:,1) = initial_vector
            else
                do i = 1, n
                    workspace%basis(i,1) = 1.0_prec + real(i,prec)/real(n,prec)
                end do
            end if

            start_scale = maxval(abs(workspace%basis(:,1)))
            if (.not. ieee_is_finite(start_scale) .or. start_scale == 0.0_prec) then
                call set_failure(result, EIGEN_ZERO_START_VECTOR, &
                                 'Arnoldi initial vector has zero norm')
                return
            end if
            workspace%basis(:,1) = workspace%basis(:,1)/start_scale
            start_norm = robust_real_norm(workspace%basis(:,1))
            if (.not. ieee_is_finite(start_norm) .or. start_norm == 0.0_prec) then
                call set_failure(result, EIGEN_NONFINITE_VALUE, &
                                 'unable to normalize Arnoldi initial vector')
                return
            end if
            workspace%basis(:,1) = workspace%basis(:,1)/start_norm

            actual_dimension = 0
            do j = 1, requested_dimension
                call problem%apply_right(workspace%basis(:,j), &
                                         workspace%right_hand_side, component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    call set_failure(result, component_status, &
                                     'right operator application failed')
                    return
                end if
                result%operator_applications = result%operator_applications + 1

                call problem%solve_left(workspace%right_hand_side, &
                                        workspace%work_vector, component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    call set_failure(result, component_status, &
                                     'left-hand linear solve failed')
                    return
                end if
                if (.not. all(ieee_is_finite(workspace%work_vector))) then
                    call set_failure(result, EIGEN_NONFINITE_VALUE, &
                                     'operator produced a nonfinite Arnoldi vector')
                    return
                end if

                operator_norm = robust_real_norm(workspace%work_vector)
                if (.not. ieee_is_finite(operator_norm)) then
                    call set_failure(result, EIGEN_NONFINITE_VALUE, &
                                     'operator produced a nonfinite Arnoldi norm')
                    return
                end if
                do pass = 1, 2
                    do i = 1, j
                        coefficient = dot_product(workspace%basis(:,i), &
                                                  workspace%work_vector)
                        workspace%hessenberg(i,j) = &
                            workspace%hessenberg(i,j) + coefficient
                        workspace%work_vector = workspace%work_vector - &
                            coefficient*workspace%basis(:,i)
                    end do
                end do

                next_norm = robust_real_norm(workspace%work_vector)
                if (.not. ieee_is_finite(next_norm)) then
                    call set_failure(result, EIGEN_NONFINITE_VALUE, &
                                     'orthogonalization produced a nonfinite norm')
                    return
                end if
                workspace%hessenberg(j+1,j) = next_norm
                actual_dimension = j

                expansion_breakdown = .not. (operator_norm > 0.0_prec .and. &
                    next_norm > options%breakdown_tolerance*operator_norm)
                if (expansion_breakdown) then
                    result%arnoldi_breakdowns = result%arnoldi_breakdowns + 1
                    workspace%hessenberg(j+1,j) = 0.0_prec
                end if

                if (j < requested_dimension) then
                    if (.not. expansion_breakdown) then
                        workspace%basis(:,j+1) = workspace%work_vector/next_norm
                    else
                        call build_orthogonal_complement(workspace%basis(:,1:j), &
                            workspace%basis(:,j+1), options%breakdown_tolerance, &
                            complement_found)
                        if (.not. complement_found) exit
                    end if
                end if
            end do

            result%arnoldi_steps = actual_dimension
            if (actual_dimension < options%number_of_eigenvalues) then
                call set_failure(result, EIGEN_NOT_CONVERGED, &
                                 'Arnoldi subspace is smaller than requested')
                return
            end if

            call solve_projected_problem( &
                workspace%hessenberg(1:actual_dimension,1:actual_dimension), &
                projected_values, projected_vectors, component_status)
            if (component_status /= EIGEN_SUCCESS) then
                call set_failure(result, component_status, &
                                 'projected Hessenberg eigensolve failed')
                return
            end if

            call sort_eigenpairs(projected_values, projected_vectors, component_status)
            if (component_status /= EIGEN_SUCCESS) then
                call set_failure(result, component_status, &
                                 'unable to sort Arnoldi Ritz pairs')
                return
            end if
            number_to_return = selected_eigenvalue_count(projected_values, &
                options%number_of_eigenvalues)
            call allocate_result_arrays(result, n, number_to_return, allocation_status)
            if (allocation_status /= 0) then
                call set_failure(result, EIGEN_ALLOCATION_FAILURE, &
                                 'unable to allocate Arnoldi result')
                return
            end if

            result%eigenvalues = projected_values(1:number_to_return)
            do i = 1, number_to_return
                result%eigenvectors(:,i) = cmplx( &
                    matmul(workspace%basis(:,1:actual_dimension), &
                           real(projected_vectors(:,i), kind=prec)), &
                    matmul(workspace%basis(:,1:actual_dimension), &
                           aimag(projected_vectors(:,i))), kind=prec)
                call normalize_and_fix_phase(result%eigenvectors(:,i), component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    call set_failure(result, component_status, &
                                     'invalid Arnoldi Ritz vector')
                    return
                end if

                result%ritz_residual_estimates(i) = &
                    abs(workspace%hessenberg(actual_dimension+1,actual_dimension))* &
                    abs(projected_vectors(actual_dimension,i))
                call evaluate_generalized_residual(problem, workspace, &
                    result%eigenvalues(i), result%eigenvectors(:,i), &
                    result%absolute_residuals(i), &
                    result%relative_residuals(i), component_status, &
                    residual_applications)
                result%operator_applications = result%operator_applications + &
                                               residual_applications
                if (component_status /= EIGEN_SUCCESS) then
                    call set_failure(result, component_status, &
                                     'unable to verify an Arnoldi Ritz pair')
                    return
                end if
            end do

            result%returned_eigenvalues = number_to_return
            result%converged_eigenvalues = count( &
                result%relative_residuals <= options%residual_tolerance)
            if (result%converged_eigenvalues == number_to_return) then
                call set_success(result, 'dominant eigenpairs converged')
            else
                call set_failure(result, EIGEN_NOT_CONVERGED, &
                                 'single Arnoldi cycle exhausted; enlarge subspace')
            end if
        end subroutine solve_dominant_arnoldi

        subroutine solve_projected_problem(hessenberg, eigenvalues, eigenvectors, status)
            real(prec), intent(in) :: hessenberg(:,:)
            complex(prec), allocatable, intent(out) :: eigenvalues(:)
            complex(prec), allocatable, intent(out) :: eigenvectors(:,:)
            integer, intent(out) :: status

            real(prec), allocatable :: matrix_work(:,:), values_real(:), values_imaginary(:)
            real(prec), allocatable :: lapack_vectors(:,:)
            integer :: allocation_status, info, n

            status = EIGEN_ALLOCATION_FAILURE
            n = size(hessenberg, 1)
            allocate(matrix_work(n,n), values_real(n), values_imaginary(n), &
                     lapack_vectors(n,n), stat=allocation_status)
            if (allocation_status /= 0) return

            matrix_work = hessenberg
            call lapack_geev(matrix_work, values_real, values_imaginary, &
                             lapack_vectors, info)
            if (info /= 0) then
                if (info == LAPACK_WRAPPER_ALLOCATION_FAILURE) then
                    status = EIGEN_ALLOCATION_FAILURE
                else
                    status = EIGEN_LAPACK_FAILURE
                end if
                return
            end if

            call decode_standard_eigensystem(values_real, values_imaginary, &
                                              lapack_vectors, eigenvalues, &
                                              eigenvectors, status)
        end subroutine solve_projected_problem



        ! Dense problem implementation

        subroutine bind_dense_problem(this, left_matrix, right_matrix, status)
            class(dense_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: left_matrix(:,:), right_matrix(:,:)
            integer, intent(out) :: status

            integer :: allocation_status, n

            call this%clear()
            status = EIGEN_INVALID_INPUT

            n = size(left_matrix, 1)
            if (n <= 0 .or. size(left_matrix, 2) /= n) return
            if (size(right_matrix, 1) /= n .or. size(right_matrix, 2) /= n) return
            if (.not. all(ieee_is_finite(left_matrix)) .or. &
                .not. all(ieee_is_finite(right_matrix))) return

            allocate(this%left_matrix(n,n), this%right_matrix(n,n), &
                     stat=allocation_status)
            if (allocation_status /= 0) then
                call this%clear()
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if

            this%left_matrix = left_matrix
            this%right_matrix = right_matrix
            this%n = n
            this%ready = .true.
            status = EIGEN_SUCCESS
        end subroutine bind_dense_problem

        subroutine clear_dense_problem(this)
            class(dense_generalized_problem_t), intent(inout) :: this

            if (allocated(this%left_matrix)) deallocate(this%left_matrix)
            if (allocated(this%right_matrix)) deallocate(this%right_matrix)
            call discard_dense_factorization(this)
            this%n = 0
            this%ready = .false.
            this%linear_solve_count = 0
        end subroutine clear_dense_problem

        subroutine discard_dense_factorization(this)
            class(dense_generalized_problem_t), intent(inout) :: this

            if (allocated(this%left_lu)) deallocate(this%left_lu)
            if (allocated(this%solve_buffer)) deallocate(this%solve_buffer)
            if (allocated(this%pivots)) deallocate(this%pivots)
            this%factorization_ready = .false.
            this%factorization_scale = 1.0_prec
            this%reciprocal_condition = -1.0_prec
        end subroutine discard_dense_factorization

        subroutine ensure_dense_factorization(this, status)
            class(dense_generalized_problem_t), intent(inout) :: this
            integer, intent(out) :: status

            real(prec) :: matrix_norm
            integer :: allocation_status, info

            if (this%factorization_ready) then
                status = EIGEN_SUCCESS
                return
            end if
            status = EIGEN_FACTORIZATION_FAILURE
            if (.not. this%is_valid()) return

            this%factorization_scale = maxval(abs(this%left_matrix))
            if (.not. ieee_is_finite(this%factorization_scale) .or. &
                this%factorization_scale <= 0.0_prec) return

            allocate(this%left_lu(this%n,this%n), this%pivots(this%n), &
                     this%solve_buffer(this%n,1), stat=allocation_status)
            if (allocation_status /= 0) then
                call discard_dense_factorization(this)
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if

            this%left_lu = this%left_matrix/this%factorization_scale
            matrix_norm = matrix_one_norm(this%left_lu)
            if (.not. ieee_is_finite(matrix_norm) .or. matrix_norm <= 0.0_prec) then
                call discard_dense_factorization(this)
                status = EIGEN_FACTORIZATION_FAILURE
                return
            end if
            call lapack_getrf(this%left_lu, this%pivots, info)
            if (info /= 0) then
                call discard_dense_factorization(this)
                if (info > 0) then
                    status = EIGEN_FACTORIZATION_FAILURE
                else
                    status = EIGEN_LAPACK_FAILURE
                end if
                return
            end if

            call lapack_gecon(this%left_lu, matrix_norm, &
                              this%reciprocal_condition, info)
            if (info == LAPACK_WRAPPER_ALLOCATION_FAILURE) then
                call discard_dense_factorization(this)
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if
            if (info /= 0 .or. .not. ieee_is_finite(this%reciprocal_condition) .or. &
                this%reciprocal_condition <= 0.0_prec) then
                call discard_dense_factorization(this)
                status = EIGEN_FACTORIZATION_FAILURE
                return
            end if

            this%factorization_ready = .true.
            status = EIGEN_SUCCESS
        end subroutine ensure_dense_factorization

        real(prec) function dense_condition_estimate(this)
            class(dense_generalized_problem_t), intent(in) :: this

            dense_condition_estimate = this%reciprocal_condition
        end function dense_condition_estimate

        integer function dense_problem_order(this)
            class(dense_generalized_problem_t), intent(in) :: this

            dense_problem_order = this%n
        end function dense_problem_order

        logical function dense_problem_is_valid(this)
            class(dense_generalized_problem_t), intent(in) :: this

            dense_problem_is_valid = this%ready .and. this%n > 0
            dense_problem_is_valid = dense_problem_is_valid .and. &
                allocated(this%left_matrix) .and. allocated(this%right_matrix)
            if (.not. dense_problem_is_valid) return
            dense_problem_is_valid = size(this%left_matrix,1) == this%n .and. &
                size(this%left_matrix,2) == this%n .and. &
                size(this%right_matrix,1) == this%n .and. &
                size(this%right_matrix,2) == this%n
        end function dense_problem_is_valid

        subroutine apply_dense_left(this, input, output, status)
            class(dense_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            output = 0.0_prec
            status = EIGEN_OPERATOR_FAILURE
            if (.not. this%is_valid()) return
            call apply_dense_matrix(this%left_matrix, this%n, this%ready, &
                                    input, output, status)
        end subroutine apply_dense_left

        subroutine apply_dense_right(this, input, output, status)
            class(dense_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            output = 0.0_prec
            status = EIGEN_OPERATOR_FAILURE
            if (.not. this%is_valid()) return
            call apply_dense_matrix(this%right_matrix, this%n, this%ready, &
                                    input, output, status)
        end subroutine apply_dense_right

        subroutine apply_dense_matrix(matrix, n, ready, input, output, status)
            real(prec), intent(in) :: matrix(:,:)
            integer, intent(in) :: n
            logical, intent(in) :: ready
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            output = 0.0_prec
            status = EIGEN_OPERATOR_FAILURE
            if (.not. ready) return
            if (size(matrix,1) /= n .or. size(matrix,2) /= n) return
            if (size(input) /= n .or. size(output) /= n) return
            if (.not. all(ieee_is_finite(input))) return

            output = matmul(matrix, input)
            if (.not. all(ieee_is_finite(output))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine apply_dense_matrix

        subroutine solve_dense_left(this, right_hand_side, solution, status)
            class(dense_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: right_hand_side(:)
            real(prec), intent(out) :: solution(:)
            integer, intent(out) :: status

            integer :: component_status, info

            solution = 0.0_prec
            status = EIGEN_FACTORIZATION_FAILURE
            if (.not. this%is_valid()) return
            if (size(right_hand_side) /= this%n .or. size(solution) /= this%n) return
            if (.not. all(ieee_is_finite(right_hand_side))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if

            call ensure_dense_factorization(this, component_status)
            if (component_status /= EIGEN_SUCCESS) then
                status = component_status
                return
            end if

            this%solve_buffer(:,1) = right_hand_side/this%factorization_scale
            if (.not. all(ieee_is_finite(this%solve_buffer(:,1)))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            this%linear_solve_count = this%linear_solve_count + 1
            call lapack_getrs(this%left_lu, this%pivots, this%solve_buffer, info)
            if (info /= 0) then
                status = EIGEN_LAPACK_FAILURE
                return
            end if

            solution = this%solve_buffer(:,1)
            if (.not. all(ieee_is_finite(solution))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine solve_dense_left

        subroutine reset_dense_statistics(this)
            class(dense_generalized_problem_t), intent(inout) :: this

            this%linear_solve_count = 0
        end subroutine reset_dense_statistics

        subroutine get_dense_statistics(this, linear_solves, inner_iterations)
            class(dense_generalized_problem_t), intent(in) :: this
            integer, intent(out) :: linear_solves, inner_iterations

            linear_solves = this%linear_solve_count
            inner_iterations = 0
        end subroutine get_dense_statistics



        ! Vectorized problem implementation

        subroutine bind_vectorized_problem(this, left_operator, right_operator, &
                                           inner_options, status, preconditioner)
            class(vectorized_generalized_problem_t), intent(inout) :: this
            class(linear_operator_t), intent(in) :: left_operator
            class(linear_operator_t), intent(in) :: right_operator
            type(bicgstab_options_t), intent(in) :: inner_options
            integer, intent(out) :: status
            class(preconditioner_t), optional, intent(in) :: preconditioner

            integer :: allocation_status

            call this%clear()
            status = EIGEN_INVALID_INPUT

            if (.not. left_operator%is_valid() .or. .not. right_operator%is_valid()) return
            if (left_operator%order() <= 0 .or. &
                left_operator%order() /= right_operator%order()) return
            if (inner_options%max_iterations < 0 .or. &
                inner_options%relative_tolerance < 0.0_prec .or. &
                inner_options%absolute_tolerance < 0.0_prec .or. &
                inner_options%breakdown_tolerance < 0.0_prec) return
            if (.not. ieee_is_finite(inner_options%relative_tolerance) .or. &
                .not. ieee_is_finite(inner_options%absolute_tolerance) .or. &
                .not. ieee_is_finite(inner_options%breakdown_tolerance)) return

            allocate(this%left_operator, source=left_operator, stat=allocation_status)
            if (allocation_status /= 0) then
                call this%clear()
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if
            allocate(this%right_operator, source=right_operator, stat=allocation_status)
            if (allocation_status /= 0) then
                call this%clear()
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if
            if (present(preconditioner)) then
                allocate(this%preconditioner, source=preconditioner, &
                         stat=allocation_status)
                if (allocation_status /= 0) then
                    call this%clear()
                    status = EIGEN_ALLOCATION_FAILURE
                    return
                end if
            end if
            this%n = left_operator%order()
            this%inner_options = inner_options
            this%ready = .true.
            status = EIGEN_SUCCESS
        end subroutine bind_vectorized_problem

        subroutine clear_vectorized_problem(this)
            class(vectorized_generalized_problem_t), intent(inout) :: this

            if (allocated(this%left_operator)) deallocate(this%left_operator)
            if (allocated(this%right_operator)) deallocate(this%right_operator)
            if (allocated(this%preconditioner)) deallocate(this%preconditioner)
            call this%inner_workspace%clear()
            this%n = 0
            this%ready = .false.
            this%linear_solve_count = 0
            this%total_inner_iterations = 0
        end subroutine clear_vectorized_problem

        integer function vectorized_problem_order(this)
            class(vectorized_generalized_problem_t), intent(in) :: this

            vectorized_problem_order = this%n
        end function vectorized_problem_order

        logical function vectorized_problem_is_valid(this)
            class(vectorized_generalized_problem_t), intent(in) :: this

            vectorized_problem_is_valid = this%ready .and. this%n > 0
            vectorized_problem_is_valid = vectorized_problem_is_valid .and. &
                allocated(this%left_operator) .and. allocated(this%right_operator)
            if (.not. vectorized_problem_is_valid) return
            vectorized_problem_is_valid = this%left_operator%is_valid() .and. &
                this%right_operator%is_valid()
            vectorized_problem_is_valid = vectorized_problem_is_valid .and. &
                this%left_operator%order() == this%n .and. &
                this%right_operator%order() == this%n
        end function vectorized_problem_is_valid

        subroutine apply_vectorized_left(this, input, output, status)
            class(vectorized_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            integer :: component_status

            output = 0.0_prec
            status = EIGEN_OPERATOR_FAILURE
            if (.not. this%is_valid()) return
            if (size(input) /= this%n .or. size(output) /= this%n) return

            call this%left_operator%apply(input, output, component_status)
            if (component_status /= BICGSTAB_COMPONENT_SUCCESS) return
            if (.not. all(ieee_is_finite(output))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine apply_vectorized_left

        subroutine apply_vectorized_right(this, input, output, status)
            class(vectorized_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            integer :: component_status

            output = 0.0_prec
            status = EIGEN_OPERATOR_FAILURE
            if (.not. this%is_valid()) return
            if (size(input) /= this%n .or. size(output) /= this%n) return

            call this%right_operator%apply(input, output, component_status)
            if (component_status /= BICGSTAB_COMPONENT_SUCCESS) return
            if (.not. all(ieee_is_finite(output))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine apply_vectorized_right

        subroutine solve_vectorized_left(this, right_hand_side, solution, status)
            class(vectorized_generalized_problem_t), intent(inout) :: this
            real(prec), intent(in) :: right_hand_side(:)
            real(prec), intent(out) :: solution(:)
            integer, intent(out) :: status

            solution = 0.0_prec
            status = EIGEN_INNER_SOLVER_FAILURE
            if (.not. this%is_valid()) return
            if (size(right_hand_side) /= this%n .or. size(solution) /= this%n) return
            if (.not. all(ieee_is_finite(right_hand_side))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if

            this%linear_solve_count = this%linear_solve_count + 1
            if (allocated(this%preconditioner)) then
                call solve_bicgstab(this%left_operator, right_hand_side, solution, &
                    this%inner_options, this%inner_workspace, this%last_inner_result, &
                    this%preconditioner)
            else
                call solve_bicgstab(this%left_operator, right_hand_side, solution, &
                    this%inner_options, this%inner_workspace, this%last_inner_result)
            end if
            this%total_inner_iterations = this%total_inner_iterations + &
                                          this%last_inner_result%iterations

            if (this%last_inner_result%status /= BICGSTAB_SUCCESS) return
            if (.not. all(ieee_is_finite(solution))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine solve_vectorized_left

        subroutine reset_vectorized_statistics(this)
            class(vectorized_generalized_problem_t), intent(inout) :: this

            this%linear_solve_count = 0
            this%total_inner_iterations = 0
        end subroutine reset_vectorized_statistics

        subroutine get_vectorized_statistics(this, linear_solves, inner_iterations)
            class(vectorized_generalized_problem_t), intent(in) :: this
            integer, intent(out) :: linear_solves, inner_iterations

            linear_solves = this%linear_solve_count
            inner_iterations = this%total_inner_iterations
        end subroutine get_vectorized_statistics



        ! Eigenpair decoding, ordering, and validation

        subroutine decode_standard_eigensystem(values_real, values_imaginary, &
                                               lapack_vectors, eigenvalues, &
                                               eigenvectors, status)
            real(prec), intent(in) :: values_real(:), values_imaginary(:)
            real(prec), intent(in) :: lapack_vectors(:,:)
            complex(prec), allocatable, intent(out) :: eigenvalues(:)
            complex(prec), allocatable, intent(out) :: eigenvectors(:,:)
            integer, intent(out) :: status

            integer :: allocation_status, column, n

            n = size(values_real)
            allocate(eigenvalues(n), eigenvectors(n,n), stat=allocation_status)
            if (allocation_status /= 0) then
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if

            eigenvalues = cmplx(values_real, values_imaginary, kind=prec)
            column = 1
            do while (column <= n)
                if (values_imaginary(column) == 0.0_prec) then
                    eigenvectors(:,column) = cmplx(lapack_vectors(:,column), &
                                                   0.0_prec, kind=prec)
                    column = column + 1
                elseif (values_imaginary(column) > 0.0_prec) then
                    if (column >= n) then
                        status = EIGEN_LAPACK_FAILURE
                        return
                    end if
                    if (values_imaginary(column+1) >= 0.0_prec) then
                        status = EIGEN_LAPACK_FAILURE
                        return
                    end if
                    eigenvectors(:,column) = cmplx(lapack_vectors(:,column), &
                                                   lapack_vectors(:,column+1), kind=prec)
                    eigenvalues(column+1) = conjg(eigenvalues(column))
                    eigenvectors(:,column+1) = conjg(eigenvectors(:,column))
                    column = column + 2
                else
                    status = EIGEN_LAPACK_FAILURE
                    return
                end if
            end do

            if (.not. complex_array_is_finite(eigenvalues) .or. &
                .not. all(ieee_is_finite(real(eigenvectors,kind=prec))) .or. &
                .not. all(ieee_is_finite(aimag(eigenvectors)))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine decode_standard_eigensystem

        subroutine decode_generalized_eigensystem(alpha_real, alpha_imaginary, beta, &
                                                  lapack_vectors, eigenvalues, &
                                                  eigenvectors, status)
            real(prec), intent(in) :: alpha_real(:), alpha_imaginary(:), beta(:)
            real(prec), intent(in) :: lapack_vectors(:,:)
            complex(prec), allocatable, intent(out) :: eigenvalues(:)
            complex(prec), allocatable, intent(out) :: eigenvectors(:,:)
            integer, intent(out) :: status

            real(prec) :: value_scale
            integer :: allocation_status, column, n

            n = size(alpha_real)
            allocate(eigenvalues(n), eigenvectors(n,n), stat=allocation_status)
            if (allocation_status /= 0) then
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if

            do column = 1, n
                value_scale = max(abs(alpha_real(column)), abs(alpha_imaginary(column)), &
                                  abs(beta(column)))
                if (beta(column) == 0.0_prec) then
                    if (alpha_real(column) == 0.0_prec .and. &
                        alpha_imaginary(column) == 0.0_prec) then
                        status = EIGEN_SINGULAR_PENCIL
                    else
                        status = EIGEN_INFINITE_EIGENVALUE
                    end if
                    return
                end if
                eigenvalues(column) = cmplx(alpha_real(column)/value_scale, &
                    alpha_imaginary(column)/value_scale, kind=prec)/ &
                    cmplx(beta(column)/value_scale,0.0_prec,kind=prec)
            end do

            column = 1
            do while (column <= n)
                if (alpha_imaginary(column) == 0.0_prec) then
                    eigenvectors(:,column) = cmplx(lapack_vectors(:,column), &
                                                   0.0_prec, kind=prec)
                    column = column + 1
                elseif (alpha_imaginary(column) > 0.0_prec) then
                    if (column >= n) then
                        status = EIGEN_LAPACK_FAILURE
                        return
                    end if
                    if (alpha_imaginary(column+1) >= 0.0_prec) then
                        status = EIGEN_LAPACK_FAILURE
                        return
                    end if
                    eigenvectors(:,column) = cmplx(lapack_vectors(:,column), &
                                                   lapack_vectors(:,column+1), kind=prec)
                    eigenvalues(column+1) = conjg(eigenvalues(column))
                    eigenvectors(:,column+1) = conjg(eigenvectors(:,column))
                    column = column + 2
                else
                    status = EIGEN_LAPACK_FAILURE
                    return
                end if
            end do

            if (.not. complex_array_is_finite(eigenvalues) .or. &
                .not. all(ieee_is_finite(real(eigenvectors,kind=prec))) .or. &
                .not. all(ieee_is_finite(aimag(eigenvectors)))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine decode_generalized_eigensystem

        subroutine evaluate_generalized_residual(problem, workspace, eigenvalue, &
                                                 eigenvector, absolute_residual, &
                                                 relative_residual, status, &
                                                 operator_applications)
            class(generalized_eigenproblem_t), intent(inout) :: problem
            type(eigensolver_workspace_t), intent(inout) :: workspace
            complex(prec), intent(in) :: eigenvalue, eigenvector(:)
            real(prec), intent(out) :: absolute_residual, relative_residual
            integer, intent(out) :: status
            integer, intent(out) :: operator_applications

            real(prec) :: eigenvalue_imaginary, eigenvalue_real
            real(prec) :: left_norm, right_norm
            integer :: component_status, n

            absolute_residual = huge(1.0_prec)
            relative_residual = huge(1.0_prec)
            operator_applications = 0
            status = EIGEN_NONFINITE_VALUE
            n = problem%order()
            if (size(eigenvector) /= n) return
            if (.not. ieee_is_finite(real(eigenvalue,kind=prec)) .or. &
                .not. ieee_is_finite(aimag(eigenvalue)) .or. &
                .not. complex_array_is_finite(eigenvector)) return
            if (.not. allocated(workspace%residual_left_real) .or. &
                .not. allocated(workspace%residual_left_imaginary) .or. &
                .not. allocated(workspace%residual_right_real) .or. &
                .not. allocated(workspace%residual_right_imaginary)) then
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if
            if (size(workspace%residual_left_real) /= n .or. &
                size(workspace%residual_left_imaginary) /= n .or. &
                size(workspace%residual_right_real) /= n .or. &
                size(workspace%residual_right_imaginary) /= n) then
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if

            if (any(real(eigenvector,kind=prec) /= 0.0_prec)) then
                call problem%apply_left(real(eigenvector, kind=prec), &
                                        workspace%residual_left_real, &
                                        component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    status = component_status
                    return
                end if
                operator_applications = operator_applications + 1
            else
                workspace%residual_left_real = 0.0_prec
            end if
            if (any(aimag(eigenvector) /= 0.0_prec)) then
                call problem%apply_left(aimag(eigenvector), &
                                        workspace%residual_left_imaginary, &
                                        component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    status = component_status
                    return
                end if
                operator_applications = operator_applications + 1
            else
                workspace%residual_left_imaginary = 0.0_prec
            end if
            if (any(real(eigenvector,kind=prec) /= 0.0_prec)) then
                call problem%apply_right(real(eigenvector, kind=prec), &
                                         workspace%residual_right_real, &
                                         component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    status = component_status
                    return
                end if
                operator_applications = operator_applications + 1
            else
                workspace%residual_right_real = 0.0_prec
            end if
            if (any(aimag(eigenvector) /= 0.0_prec)) then
                call problem%apply_right(aimag(eigenvector), &
                                         workspace%residual_right_imaginary, &
                                         component_status)
                if (component_status /= EIGEN_SUCCESS) then
                    status = component_status
                    return
                end if
                operator_applications = operator_applications + 1
            else
                workspace%residual_right_imaginary = 0.0_prec
            end if

            if (.not. all(ieee_is_finite(workspace%residual_left_real)) .or. &
                .not. all(ieee_is_finite(workspace%residual_left_imaginary)) .or. &
                .not. all(ieee_is_finite(workspace%residual_right_real)) .or. &
                .not. all(ieee_is_finite(workspace%residual_right_imaginary))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if

            left_norm = robust_split_complex_norm(workspace%residual_left_real, &
                                                  workspace%residual_left_imaginary)
            right_norm = robust_split_complex_norm(workspace%residual_right_real, &
                                                   workspace%residual_right_imaginary)
            if (left_norm == 0.0_prec .and. right_norm == 0.0_prec) then
                status = EIGEN_SINGULAR_PENCIL
                return
            end if
            eigenvalue_real = real(eigenvalue, kind=prec)
            eigenvalue_imaginary = aimag(eigenvalue)
            workspace%residual_right_real = workspace%residual_right_real - &
                eigenvalue_real*workspace%residual_left_real + &
                eigenvalue_imaginary*workspace%residual_left_imaginary
            workspace%residual_right_imaginary = workspace%residual_right_imaginary - &
                eigenvalue_real*workspace%residual_left_imaginary - &
                eigenvalue_imaginary*workspace%residual_left_real
            if (.not. all(ieee_is_finite(workspace%residual_right_real)) .or. &
                .not. all(ieee_is_finite(workspace%residual_right_imaginary))) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            absolute_residual = robust_split_complex_norm( &
                workspace%residual_right_real, workspace%residual_right_imaginary)
            relative_residual = safe_generalized_relative_residual( &
                absolute_residual, right_norm, left_norm, eigenvalue)

            if (.not. ieee_is_finite(absolute_residual) .or. &
                .not. ieee_is_finite(relative_residual)) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            status = EIGEN_SUCCESS
        end subroutine evaluate_generalized_residual

        subroutine normalize_and_fix_phase(vector, status)
            complex(prec), intent(inout) :: vector(:)
            integer, intent(out) :: status

            complex(prec) :: phase
            real(prec) :: vector_norm, vector_scale
            integer :: pivot_index

            if (.not. complex_array_is_finite(vector)) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            vector_scale = max(maxval(abs(real(vector,kind=prec))), &
                               maxval(abs(aimag(vector))))
            if (.not. ieee_is_finite(vector_scale) .or. vector_scale == 0.0_prec) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            vector = vector/cmplx(vector_scale,0.0_prec,kind=prec)
            vector_norm = robust_complex_norm(vector)
            if (.not. ieee_is_finite(vector_norm) .or. &
                vector_norm == 0.0_prec) then
                status = EIGEN_NONFINITE_VALUE
                return
            end if
            vector = vector/cmplx(vector_norm,0.0_prec,kind=prec)

            pivot_index = maxloc(abs(vector), dim=1)
            if (abs(vector(pivot_index)) > tiny(1.0_prec)) then
                phase = vector(pivot_index)/ &
                        cmplx(abs(vector(pivot_index)),0.0_prec,kind=prec)
                vector = vector/phase
                vector(pivot_index) = cmplx(abs(real(vector(pivot_index), kind=prec)), &
                                            0.0_prec, kind=prec)
            end if
            status = EIGEN_SUCCESS
        end subroutine normalize_and_fix_phase

        subroutine sort_eigenpairs(eigenvalues, eigenvectors, status)
            complex(prec), intent(inout) :: eigenvalues(:), eigenvectors(:,:)
            integer, intent(out) :: status

            complex(prec), allocatable :: vector_buffer(:)
            complex(prec) :: value_buffer
            integer, allocatable :: current_order(:), desired_order(:)
            logical, allocatable :: paired(:)
            real(prec) :: candidate_distance, pair_distance
            real(prec) :: pair_tolerance, pair_scale
            integer :: allocation_status, current, desired_index, original_index
            integer :: position, source_position

            if (size(eigenvectors,2) /= size(eigenvalues)) then
                status = EIGEN_INVALID_INPUT
                return
            end if
            allocate(vector_buffer(size(eigenvectors,1)), &
                     current_order(size(eigenvalues)), &
                     desired_order(size(eigenvalues)), &
                     paired(size(eigenvalues)), stat=allocation_status)
            if (allocation_status /= 0) then
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if

            do current = 1, size(eigenvalues)
                desired_order(current) = current
                current_order(current) = current
            end do

            ! Sort only integer indices first, then apply the permutation with at
            ! most n-1 full-column swaps.  This avoids the O(n^3) column traffic
            ! of direct insertion sorting.
            do current = 2, size(eigenvalues)
                original_index = desired_order(current)
                position = current - 1
                do while (position >= 1)
                    if (.not. eigenvalue_precedes(eigenvalues(original_index), &
                                                  eigenvalues(desired_order(position)))) exit
                    desired_order(position+1) = desired_order(position)
                    position = position - 1
                end do
                desired_order(position+1) = original_index
            end do

            do position = 1, size(eigenvalues)
                desired_index = desired_order(position)
                if (current_order(position) == desired_index) cycle
                source_position = findloc(current_order, desired_index, dim=1)
                if (source_position < 1) then
                    status = EIGEN_LAPACK_FAILURE
                    return
                end if

                value_buffer = eigenvalues(position)
                eigenvalues(position) = eigenvalues(source_position)
                eigenvalues(source_position) = value_buffer
                vector_buffer = eigenvectors(:,position)
                eigenvectors(:,position) = eigenvectors(:,source_position)
                eigenvectors(:,source_position) = vector_buffer
                original_index = current_order(position)
                current_order(position) = current_order(source_position)
                current_order(source_position) = original_index
            end do

            ! A real eigensolve represents every non-real eigenpair by conjugate
            ! columns.  Build a second index permutation that makes each pair
            ! adjacent without repeatedly shifting full eigenvector columns.
            paired = .false.
            desired_index = 1
            do current = 1, size(eigenvalues)
                if (paired(current)) cycle
                if (aimag(eigenvalues(current)) == 0.0_prec) then
                    desired_order(desired_index) = current
                    paired(current) = .true.
                    desired_index = desired_index + 1
                    cycle
                end if

                source_position = 0
                pair_distance = huge(1.0_prec)
                do position = 1, size(eigenvalues)
                    if (position == current .or. paired(position)) cycle
                    if (aimag(eigenvalues(current)) > 0.0_prec) then
                        if (aimag(eigenvalues(position)) >= 0.0_prec) cycle
                    else
                        if (aimag(eigenvalues(position)) <= 0.0_prec) cycle
                    end if
                    candidate_distance = abs(eigenvalues(position) - &
                                             conjg(eigenvalues(current)))
                    if (candidate_distance < pair_distance) then
                        pair_distance = candidate_distance
                        source_position = position
                    end if
                end do
                if (source_position == 0) then
                    status = EIGEN_LAPACK_FAILURE
                    return
                end if
                pair_scale = max(tiny(1.0_prec), abs(eigenvalues(current)), &
                                 abs(eigenvalues(source_position)))
                pair_tolerance = 1000.0_prec*epsilon(1.0_prec)*pair_scale
                if (pair_distance > pair_tolerance) then
                    status = EIGEN_LAPACK_FAILURE
                    return
                end if

                if (aimag(eigenvalues(current)) > 0.0_prec) then
                    desired_order(desired_index:desired_index+1) = &
                        [current, source_position]
                else
                    desired_order(desired_index:desired_index+1) = &
                        [source_position, current]
                end if
                paired(current) = .true.
                paired(source_position) = .true.
                desired_index = desired_index + 2
            end do
            if (desired_index /= size(eigenvalues) + 1) then
                status = EIGEN_LAPACK_FAILURE
                return
            end if

            do current = 1, size(eigenvalues)
                current_order(current) = current
            end do
            do position = 1, size(eigenvalues)
                desired_index = desired_order(position)
                if (current_order(position) == desired_index) cycle
                source_position = findloc(current_order, desired_index, dim=1)
                if (source_position < 1) then
                    status = EIGEN_LAPACK_FAILURE
                    return
                end if
                value_buffer = eigenvalues(position)
                eigenvalues(position) = eigenvalues(source_position)
                eigenvalues(source_position) = value_buffer
                vector_buffer = eigenvectors(:,position)
                eigenvectors(:,position) = eigenvectors(:,source_position)
                eigenvectors(:,source_position) = vector_buffer
                original_index = current_order(position)
                current_order(position) = current_order(source_position)
                current_order(source_position) = original_index
            end do
            status = EIGEN_SUCCESS
        end subroutine sort_eigenpairs

        logical function eigenvalue_precedes(candidate, current)
            complex(prec), intent(in) :: candidate, current

            real(prec) :: comparison_tolerance, candidate_magnitude, current_magnitude

            candidate_magnitude = abs(candidate)
            current_magnitude = abs(current)
            comparison_tolerance = 100.0_prec*epsilon(1.0_prec)* &
                max(tiny(1.0_prec), candidate_magnitude, current_magnitude)

            if (candidate_magnitude > current_magnitude + comparison_tolerance) then
                eigenvalue_precedes = .true.
            elseif (current_magnitude > candidate_magnitude + comparison_tolerance) then
                eigenvalue_precedes = .false.
            elseif (real(candidate,kind=prec) > &
                    real(current,kind=prec) + comparison_tolerance) then
                eigenvalue_precedes = .true.
            elseif (real(current,kind=prec) > &
                    real(candidate,kind=prec) + comparison_tolerance) then
                eigenvalue_precedes = .false.
            else
                eigenvalue_precedes = aimag(candidate) > aimag(current)
            end if
        end function eigenvalue_precedes

        integer function selected_eigenvalue_count(eigenvalues, requested_count)
            complex(prec), intent(in) :: eigenvalues(:)
            integer, intent(in) :: requested_count

            real(prec) :: pair_tolerance

            selected_eigenvalue_count = min(requested_count, size(eigenvalues))
            if (selected_eigenvalue_count >= size(eigenvalues)) return

            pair_tolerance = 1000.0_prec*epsilon(1.0_prec)* &
                max(tiny(1.0_prec), abs(eigenvalues(selected_eigenvalue_count)))
            if (aimag(eigenvalues(selected_eigenvalue_count)) > 0.0_prec) then
                if (abs(eigenvalues(selected_eigenvalue_count+1) - &
                        conjg(eigenvalues(selected_eigenvalue_count))) <= pair_tolerance) then
                    selected_eigenvalue_count = selected_eigenvalue_count + 1
                end if
            end if
        end function selected_eigenvalue_count



        ! Orthogonalization and workspace support

        subroutine build_orthogonal_complement(basis, complement, tolerance, found)
            real(prec), intent(in) :: basis(:,:)
            real(prec), intent(out) :: complement(:)
            real(prec), intent(in) :: tolerance
            logical, intent(out) :: found

            real(prec) :: coefficient, complement_norm
            integer :: coordinate, i, pass

            found = .false.
            complement = 0.0_prec
            do coordinate = 1, size(complement)
                complement = 0.0_prec
                complement(coordinate) = 1.0_prec
                do pass = 1, 2
                    do i = 1, size(basis,2)
                        coefficient = dot_product(basis(:,i), complement)
                        complement = complement - coefficient*basis(:,i)
                    end do
                end do
                complement_norm = robust_real_norm(complement)
                if (complement_norm > tolerance) then
                    complement = complement/complement_norm
                    found = .true.
                    return
                end if
            end do
        end subroutine build_orthogonal_complement

        subroutine ensure_eigensolver_workspace(this, n, krylov_dimension, status)
            class(eigensolver_workspace_t), intent(inout) :: this
            integer, intent(in) :: n, krylov_dimension
            integer, intent(out) :: status

            integer :: allocation_status

            status = EIGEN_ALLOCATION_FAILURE
            if (n <= 0 .or. krylov_dimension <= 0) return

            if (allocated(this%basis)) then
                if (size(this%basis,1) /= n .or. &
                    size(this%basis,2) /= krylov_dimension + 1) call this%clear()
            end if
            if (.not. allocated(this%basis)) then
                allocate(this%basis(n,krylov_dimension+1), &
                         this%hessenberg(krylov_dimension+1,krylov_dimension), &
                         this%work_vector(n), this%right_hand_side(n), &
                         this%residual_left_real(n), &
                         this%residual_left_imaginary(n), &
                         this%residual_right_real(n), &
                         this%residual_right_imaginary(n), &
                         stat=allocation_status)
                if (allocation_status /= 0) then
                    call this%clear()
                    return
                end if
            end if
            status = EIGEN_SUCCESS
        end subroutine ensure_eigensolver_workspace

        subroutine clear_eigensolver_workspace(this)
            class(eigensolver_workspace_t), intent(inout) :: this

            if (allocated(this%basis)) deallocate(this%basis)
            if (allocated(this%hessenberg)) deallocate(this%hessenberg)
            if (allocated(this%work_vector)) deallocate(this%work_vector)
            if (allocated(this%right_hand_side)) deallocate(this%right_hand_side)
            if (allocated(this%residual_left_real)) deallocate(this%residual_left_real)
            if (allocated(this%residual_left_imaginary)) &
                deallocate(this%residual_left_imaginary)
            if (allocated(this%residual_right_real)) deallocate(this%residual_right_real)
            if (allocated(this%residual_right_imaginary)) &
                deallocate(this%residual_right_imaginary)
        end subroutine clear_eigensolver_workspace



        ! Result and scalar support

        logical function common_inputs_are_valid(problem, options, n)
            class(generalized_eigenproblem_t), intent(in) :: problem
            type(eigensolver_options_t), intent(in) :: options
            integer, intent(in) :: n

            common_inputs_are_valid = n > 0 .and. problem%is_valid()
            common_inputs_are_valid = common_inputs_are_valid .and. &
                (options%mode == EIGEN_MODE_DOMINANT .or. &
                 options%mode == EIGEN_MODE_FULL)
            common_inputs_are_valid = common_inputs_are_valid .and. &
                options%residual_tolerance > 0.0_prec .and. &
                options%breakdown_tolerance >= 0.0_prec
            common_inputs_are_valid = common_inputs_are_valid .and. &
                ieee_is_finite(options%residual_tolerance) .and. &
                ieee_is_finite(options%breakdown_tolerance)
        end function common_inputs_are_valid

        subroutine initialize_result(result)
            type(eigensolver_result_t), intent(out) :: result

            result%converged = .false.
            result%status = EIGEN_INVALID_INPUT
            result%requested_eigenvalues = 0
            result%returned_eigenvalues = 0
            result%converged_eigenvalues = 0
            result%arnoldi_steps = 0
            result%arnoldi_breakdowns = 0
            result%operator_applications = 0
            result%linear_solves = 0
            result%inner_iterations = 0
            result%last_inner_solver_status = -1
            result%reciprocal_condition_estimate = -1.0_prec
            result%method = 'not started'
            result%reason = 'solver not started'
            result%last_inner_solver_reason = 'inner solver not used'
        end subroutine initialize_result

        subroutine allocate_result_arrays(result, vector_size, number_of_eigenvalues, status)
            type(eigensolver_result_t), intent(inout) :: result
            integer, intent(in) :: vector_size, number_of_eigenvalues
            integer, intent(out) :: status

            allocate(result%eigenvalues(number_of_eigenvalues), &
                     result%eigenvectors(vector_size,number_of_eigenvalues), &
                     result%absolute_residuals(number_of_eigenvalues), &
                     result%relative_residuals(number_of_eigenvalues), &
                     result%ritz_residual_estimates(number_of_eigenvalues), stat=status)
            if (status /= 0) then
                if (allocated(result%eigenvalues)) deallocate(result%eigenvalues)
                if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
                if (allocated(result%absolute_residuals)) deallocate(result%absolute_residuals)
                if (allocated(result%relative_residuals)) deallocate(result%relative_residuals)
                if (allocated(result%ritz_residual_estimates)) &
                    deallocate(result%ritz_residual_estimates)
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if
            result%eigenvalues = cmplx(0.0_prec, 0.0_prec, kind=prec)
            result%eigenvectors = cmplx(0.0_prec, 0.0_prec, kind=prec)
            result%absolute_residuals = huge(1.0_prec)
            result%relative_residuals = huge(1.0_prec)
            result%ritz_residual_estimates = huge(1.0_prec)
            status = EIGEN_SUCCESS
        end subroutine allocate_result_arrays

        subroutine allocate_result_diagnostics(result, number_of_eigenvalues, status)
            type(eigensolver_result_t), intent(inout) :: result
            integer, intent(in) :: number_of_eigenvalues
            integer, intent(out) :: status

            allocate(result%absolute_residuals(number_of_eigenvalues), &
                     result%relative_residuals(number_of_eigenvalues), &
                     result%ritz_residual_estimates(number_of_eigenvalues), &
                     stat=status)
            if (status /= 0) then
                if (allocated(result%absolute_residuals)) &
                    deallocate(result%absolute_residuals)
                if (allocated(result%relative_residuals)) &
                    deallocate(result%relative_residuals)
                if (allocated(result%ritz_residual_estimates)) &
                    deallocate(result%ritz_residual_estimates)
                status = EIGEN_ALLOCATION_FAILURE
                return
            end if
            result%absolute_residuals = huge(1.0_prec)
            result%relative_residuals = huge(1.0_prec)
            result%ritz_residual_estimates = huge(1.0_prec)
            status = EIGEN_SUCCESS
        end subroutine allocate_result_diagnostics

        subroutine set_success(result, reason)
            type(eigensolver_result_t), intent(inout) :: result
            character(len=*), intent(in) :: reason

            result%converged = .true.
            result%status = EIGEN_SUCCESS
            result%reason = reason
        end subroutine set_success

        subroutine set_failure(result, status, reason)
            type(eigensolver_result_t), intent(inout) :: result
            integer, intent(in) :: status
            character(len=*), intent(in) :: reason

            result%converged = .false.
            result%status = status
            result%reason = reason
        end subroutine set_failure

        subroutine capture_problem_statistics(problem, result)
            class(generalized_eigenproblem_t), intent(in) :: problem
            type(eigensolver_result_t), intent(inout) :: result

            call problem%get_statistics(result%linear_solves, result%inner_iterations)
            select type (concrete_problem => problem)
            type is (dense_generalized_problem_t)
                result%reciprocal_condition_estimate = &
                    concrete_problem%reciprocal_condition
            type is (vectorized_generalized_problem_t)
                result%reciprocal_condition_estimate = -1.0_prec
                if (result%linear_solves > 0) then
                    result%last_inner_solver_status = &
                        concrete_problem%last_inner_result%status
                    result%last_inner_solver_reason = &
                        concrete_problem%last_inner_result%reason
                end if
            class default
                result%reciprocal_condition_estimate = -1.0_prec
            end select
        end subroutine capture_problem_statistics

        real(prec) function safe_generalized_relative_residual(absolute_residual, &
                                                                right_norm, left_norm, &
                                                                eigenvalue)
            real(prec), intent(in) :: absolute_residual, right_norm, left_norm
            complex(prec), intent(in) :: eigenvalue

            real(prec) :: eigenvalue_scale, eigenvalue_shape
            real(prec) :: log_denominator, log_ratio, log_right
            real(prec) :: log_weighted_left, maximum_log, scaled_denominator
            logical :: has_right_term, has_weighted_left_term

            safe_generalized_relative_residual = huge(1.0_prec)
            if (.not. ieee_is_finite(absolute_residual) .or. &
                .not. ieee_is_finite(right_norm) .or. &
                .not. ieee_is_finite(left_norm) .or. &
                absolute_residual < 0.0_prec .or. right_norm < 0.0_prec .or. &
                left_norm < 0.0_prec .or. &
                .not. ieee_is_finite(real(eigenvalue,kind=prec)) .or. &
                .not. ieee_is_finite(aimag(eigenvalue))) return
            if (absolute_residual == 0.0_prec) then
                safe_generalized_relative_residual = 0.0_prec
                return
            end if

            has_right_term = right_norm > 0.0_prec
            if (has_right_term) log_right = log(right_norm)

            eigenvalue_scale = max(abs(real(eigenvalue,kind=prec)), &
                                   abs(aimag(eigenvalue)))
            has_weighted_left_term = eigenvalue_scale > 0.0_prec .and. &
                                     left_norm > 0.0_prec
            if (has_weighted_left_term) then
                eigenvalue_shape = sqrt( &
                    (real(eigenvalue,kind=prec)/eigenvalue_scale)**2 + &
                    (aimag(eigenvalue)/eigenvalue_scale)**2)
                log_weighted_left = log(eigenvalue_scale) + &
                                    log(eigenvalue_shape) + log(left_norm)
            end if

            if (.not. has_right_term .and. .not. has_weighted_left_term) return
            if (has_right_term .and. has_weighted_left_term) then
                maximum_log = max(log_right, log_weighted_left)
                scaled_denominator = exp(log_right-maximum_log) + &
                                     exp(log_weighted_left-maximum_log)
            elseif (has_right_term) then
                maximum_log = log_right
                scaled_denominator = 1.0_prec
            else
                maximum_log = log_weighted_left
                scaled_denominator = 1.0_prec
            end if
            log_denominator = maximum_log + log(scaled_denominator)
            log_ratio = log(absolute_residual) - log_denominator

            if (log_ratio >= log(huge(1.0_prec))) then
                safe_generalized_relative_residual = huge(1.0_prec)
            else
                ! EXP preserves representable subnormal ratios on IEEE targets;
                ! clamping at TINY would create false convergence for users who
                ! deliberately request a subnormal residual tolerance.
                safe_generalized_relative_residual = exp(log_ratio)
            end if
        end function safe_generalized_relative_residual

        real(prec) function matrix_one_norm(matrix)
            real(prec), intent(in) :: matrix(:,:)

            real(prec) :: column_sum
            integer :: column

            matrix_one_norm = 0.0_prec
            do column = 1, size(matrix,2)
                column_sum = sum(abs(matrix(:,column)))
                matrix_one_norm = max(matrix_one_norm, column_sum)
            end do
        end function matrix_one_norm

        real(prec) function robust_real_norm(vector)
            real(prec), intent(in) :: vector(:)

            real(prec) :: absolute_value, scale, sum_of_squares
            integer :: index

            if (.not. all(ieee_is_finite(vector))) then
                robust_real_norm = ieee_value(0.0_prec, ieee_quiet_nan)
                return
            end if

            scale = 0.0_prec
            sum_of_squares = 1.0_prec
            do index = 1, size(vector)
                absolute_value = abs(vector(index))
                if (absolute_value > 0.0_prec) then
                    if (scale < absolute_value) then
                        sum_of_squares = 1.0_prec + sum_of_squares*(scale/absolute_value)**2
                        scale = absolute_value
                    else
                        sum_of_squares = sum_of_squares + (absolute_value/scale)**2
                    end if
                end if
            end do
            if (scale == 0.0_prec) then
                robust_real_norm = 0.0_prec
            else
                robust_real_norm = scale*sqrt(sum_of_squares)
            end if
        end function robust_real_norm

        real(prec) function robust_complex_norm(vector)
            complex(prec), intent(in) :: vector(:)

            real(prec) :: absolute_value, scale, sum_of_squares
            integer :: index

            if (.not. complex_array_is_finite(vector)) then
                robust_complex_norm = ieee_value(0.0_prec, ieee_quiet_nan)
                return
            end if

            scale = 0.0_prec
            sum_of_squares = 1.0_prec
            do index = 1, size(vector)
                absolute_value = abs(vector(index))
                if (absolute_value > 0.0_prec) then
                    if (scale < absolute_value) then
                        sum_of_squares = 1.0_prec + sum_of_squares*(scale/absolute_value)**2
                        scale = absolute_value
                    else
                        sum_of_squares = sum_of_squares + (absolute_value/scale)**2
                    end if
                end if
            end do
            if (scale == 0.0_prec) then
                robust_complex_norm = 0.0_prec
            else
                robust_complex_norm = scale*sqrt(sum_of_squares)
            end if
        end function robust_complex_norm

        real(prec) function robust_split_complex_norm(real_part, imaginary_part)
            real(prec), intent(in) :: real_part(:), imaginary_part(:)

            real(prec) :: absolute_value, scale, sum_of_squares
            integer :: component, index

            if (size(real_part) /= size(imaginary_part) .or. &
                .not. all(ieee_is_finite(real_part)) .or. &
                .not. all(ieee_is_finite(imaginary_part))) then
                robust_split_complex_norm = ieee_value(0.0_prec, ieee_quiet_nan)
                return
            end if

            scale = 0.0_prec
            sum_of_squares = 1.0_prec
            do index = 1, size(real_part)
                do component = 1, 2
                    if (component == 1) then
                        absolute_value = abs(real_part(index))
                    else
                        absolute_value = abs(imaginary_part(index))
                    end if
                    if (absolute_value > 0.0_prec) then
                        if (scale < absolute_value) then
                            sum_of_squares = 1.0_prec + &
                                sum_of_squares*(scale/absolute_value)**2
                            scale = absolute_value
                        else
                            sum_of_squares = sum_of_squares + (absolute_value/scale)**2
                        end if
                    end if
                end do
            end do
            if (scale == 0.0_prec) then
                robust_split_complex_norm = 0.0_prec
            else
                robust_split_complex_norm = scale*sqrt(sum_of_squares)
            end if
        end function robust_split_complex_norm

        logical function complex_array_is_finite(values)
            complex(prec), intent(in) :: values(:)

            complex_array_is_finite = all(ieee_is_finite(real(values,kind=prec))) .and. &
                                      all(ieee_is_finite(aimag(values)))
        end function complex_array_is_finite

end module Eigen_spectrum_solver
