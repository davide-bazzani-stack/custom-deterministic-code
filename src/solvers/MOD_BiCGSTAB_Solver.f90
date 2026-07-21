module BiCGSTAB_Solver
    ! Matrix-free BiCGSTAB with fixed right preconditioning.  The operator
    ! bindings retain non-owning views of their coefficient arrays; those
    ! arrays must remain valid while the operator or an SGS preconditioner is
    ! in use.  Solver state is held in caller-owned workspaces so independent
    ! energy-group solves can execute concurrently.
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use precision_kinds, only : prec

        implicit none

    private

    integer, parameter, public :: BICGSTAB_SUCCESS = 0
    integer, parameter, public :: BICGSTAB_MAX_ITERATIONS = 1
    integer, parameter, public :: BICGSTAB_INVALID_INPUT = 2
    integer, parameter, public :: BICGSTAB_RHO_BREAKDOWN = 3
    integer, parameter, public :: BICGSTAB_ALPHA_BREAKDOWN = 4
    integer, parameter, public :: BICGSTAB_OMEGA_BREAKDOWN = 5
    integer, parameter, public :: BICGSTAB_OPERATOR_FAILURE = 6
    integer, parameter, public :: BICGSTAB_PRECONDITIONER_FAILURE = 7
    integer, parameter, public :: BICGSTAB_NONFINITE_VALUE = 8
    integer, parameter, public :: BICGSTAB_COMPONENT_SUCCESS = 0
    integer, parameter, public :: BICGSTAB_COMPONENT_FAILURE = 1

    integer, parameter :: APPLY_SUCCESS = BICGSTAB_COMPONENT_SUCCESS
    integer, parameter :: APPLY_FAILURE = BICGSTAB_COMPONENT_FAILURE
    integer, parameter :: MAX_BAND_PAIRS = 3

    type, public :: bicgstab_options_t
        integer :: max_iterations = 1000
        real(prec) :: relative_tolerance = sqrt(epsilon(1.0_prec))
        real(prec) :: absolute_tolerance = 0.0_prec
        real(prec) :: breakdown_tolerance = 100.0_prec*epsilon(1.0_prec)
    end type bicgstab_options_t

    type, public :: bicgstab_result_t
        logical :: converged = .false.
        integer :: status = BICGSTAB_INVALID_INPUT
        integer :: iterations = 0
        real(prec) :: initial_residual_norm = huge(1.0_prec)
        real(prec) :: residual_norm = huge(1.0_prec)
        real(prec) :: relative_residual = huge(1.0_prec)
        character(len=96) :: reason = 'solver not started'
    end type bicgstab_result_t

    type, public :: bicgstab_workspace_t
        private
        real(prec), allocatable :: r(:)
        real(prec), allocatable :: r_hat(:)
        real(prec), allocatable :: v(:)
        real(prec), allocatable :: p(:)
        real(prec), allocatable :: s(:)
        real(prec), allocatable :: t(:)
        real(prec), allocatable :: p_preconditioned(:)
        real(prec), allocatable :: s_preconditioned(:)
        contains
            procedure, public :: clear => clear_workspace
            procedure, private :: ensure_capacity => ensure_workspace_capacity
    end type bicgstab_workspace_t

    type, abstract, public :: linear_operator_t
        contains
            procedure(operator_apply_interface), deferred, public :: apply
            procedure(operator_order_interface), deferred, public :: order
            procedure(operator_valid_interface), deferred, public :: is_valid
    end type linear_operator_t

    type :: real_vector_view_t
        real(prec), pointer, contiguous :: values(:) => null()
    end type real_vector_view_t

    type, extends(linear_operator_t), public :: banded_operator_t
        ! Supports one, two, or three lower/upper diagonal pairs.  Positive
        ! offsets are strictly increasing and coefficient arrays use the full
        ! system length; entries outside each diagonal's range are ignored.
        private
        integer :: n = 0
        integer :: number_of_pairs = 0
        integer :: offsets(MAX_BAND_PAIRS) = 0
        logical :: valid_layout = .false.
        real(prec), pointer, contiguous :: diagonal(:) => null()
        type(real_vector_view_t) :: lower(MAX_BAND_PAIRS)
        type(real_vector_view_t) :: upper(MAX_BAND_PAIRS)
        contains
            procedure, public :: apply => apply_banded_operator
            procedure, public :: order => banded_operator_order
            procedure, public :: is_valid => banded_operator_is_valid
            procedure, public :: bind_tridiagonal
            procedure, public :: bind_pentadiagonal
            procedure, public :: bind_heptadiagonal
            procedure, public :: bind_cartesian_tridiagonal => bind_tridiagonal
            procedure, public :: bind_cartesian_pentadiagonal => bind_pentadiagonal
            procedure, public :: bind_cartesian_heptadiagonal => bind_heptadiagonal
            procedure, public :: bind_eqtri
            procedure, public :: bind_triiso
    end type banded_operator_t

    type, abstract, public :: preconditioner_t
        contains
            procedure(preconditioner_apply_interface), deferred, public :: apply
    end type preconditioner_t

    type, extends(preconditioner_t), public :: identity_preconditioner_t
        contains
            procedure, public :: apply => apply_identity_preconditioner
    end type identity_preconditioner_t

    type, extends(preconditioner_t), public :: jacobi_preconditioner_t
        private
        real(prec), allocatable :: inverse_diagonal(:)
        logical :: ready = .false.
        contains
            procedure, public :: setup => setup_jacobi_preconditioner
            procedure, public :: apply => apply_jacobi_preconditioner
    end type jacobi_preconditioner_t

    type, extends(preconditioner_t), public :: symmetric_gauss_seidel_preconditioner_t
        private
        type(banded_operator_t) :: operator_view
        logical :: ready = .false.
        contains
            procedure, public :: setup => setup_sgs_preconditioner
            procedure, public :: apply => apply_sgs_preconditioner
    end type symmetric_gauss_seidel_preconditioner_t

    type, extends(preconditioner_t), public :: ilu0_preconditioner_t
        private
        integer :: n = 0
        integer :: number_of_pairs = 0
        integer :: offsets(MAX_BAND_PAIRS) = 0
        real(prec), allocatable :: lower_factor(:,:)
        real(prec), allocatable :: diagonal_factor(:)
        real(prec), allocatable :: upper_factor(:,:)
        logical :: ready = .false.
        contains
            procedure, public :: setup => setup_ilu0_preconditioner
            procedure, public :: apply => apply_ilu0_preconditioner
    end type ilu0_preconditioner_t

    type, extends(preconditioner_t), public :: user_preconditioner_t
        ! Both the callback procedure and its context are non-owning views.
        ! A distinct mutable context is required for each concurrent solve.
        private
        procedure(user_preconditioner_callback), pointer, nopass :: callback => null()
        class(*), pointer :: context => null()
        contains
            procedure, public :: bind => bind_user_preconditioner
            procedure, public :: apply => apply_user_preconditioner
    end type user_preconditioner_t

    abstract interface
        subroutine operator_apply_interface(this, input, output, status)
            import :: linear_operator_t, prec
            class(linear_operator_t), intent(in) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status
        end subroutine operator_apply_interface

        integer function operator_order_interface(this)
            import :: linear_operator_t
            class(linear_operator_t), intent(in) :: this
        end function operator_order_interface

        logical function operator_valid_interface(this)
            import :: linear_operator_t
            class(linear_operator_t), intent(in) :: this
        end function operator_valid_interface

        subroutine preconditioner_apply_interface(this, input, output, status)
            import :: preconditioner_t, prec
            class(preconditioner_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status
        end subroutine preconditioner_apply_interface

        subroutine user_preconditioner_callback(context, input, output, status)
            import :: prec
            class(*), intent(inout) :: context
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status
        end subroutine user_preconditioner_callback
    end interface

    public :: solve_bicgstab
    public :: user_preconditioner_callback

    contains

        ! Solver entry point and iteration

        subroutine solve_bicgstab(operator, rhs, solution, options, workspace, result, preconditioner)
            class(linear_operator_t), intent(in) :: operator
            real(prec), intent(in) :: rhs(:)
            real(prec), intent(inout) :: solution(:)
            type(bicgstab_options_t), intent(in) :: options
            type(bicgstab_workspace_t), intent(inout) :: workspace
            type(bicgstab_result_t), intent(out) :: result
            class(preconditioner_t), optional, intent(inout) :: preconditioner

            type(identity_preconditioner_t) :: identity

            ! Omitting the optional argument selects the identity preconditioner.
            ! The identity path applies A directly and does not allocate the two
            ! additional preconditioned-vector work arrays.
            if (present(preconditioner)) then
                call bicgstab_core(operator, preconditioner, rhs, solution, &
                                   options, workspace, result)
            else
                call bicgstab_core(operator, identity, rhs, solution, &
                                   options, workspace, result)
            end if
        end subroutine solve_bicgstab

        subroutine bicgstab_core(operator, preconditioner, rhs, solution, &
                                 options, workspace, result)
            class(linear_operator_t), intent(in) :: operator
            class(preconditioner_t), intent(inout) :: preconditioner
            real(prec), intent(in) :: rhs(:)
            real(prec), intent(inout) :: solution(:)
            type(bicgstab_options_t), intent(in) :: options
            type(bicgstab_workspace_t), intent(inout) :: workspace
            type(bicgstab_result_t), intent(out) :: result

            integer :: n, iteration, operation_status, workspace_status
            real(prec) :: alpha, beta, omega, rho, rho_previous
            real(prec) :: denominator, t_dot_t, rhs_norm, residual_limit
            real(prec) :: residual_norm, shadow_norm
            logical :: identity_selected, verified

            call initialize_result(result)

            n = operator%order()
            if (.not. inputs_are_valid(operator, rhs, solution, options, n)) then
                call set_result_failure(result, BICGSTAB_INVALID_INPUT, 0, &
                                        'invalid operator, dimensions, or options')
                return
            end if

            identity_selected = .false.
            select type (preconditioner)
            type is (identity_preconditioner_t)
                identity_selected = .true.
            class default
                identity_selected = .false.
            end select

            call workspace%ensure_capacity(n, .not. identity_selected, workspace_status)
            if (workspace_status /= APPLY_SUCCESS) then
                call set_result_failure(result, BICGSTAB_INVALID_INPUT, 0, &
                                        'unable to allocate solver workspace')
                return
            end if

            call operator%apply(solution, workspace%v, operation_status)
            if (operation_status /= APPLY_SUCCESS) then
                call set_result_failure(result, BICGSTAB_OPERATOR_FAILURE, 0, &
                                        'operator application failed')
                return
            end if

            workspace%r = rhs - workspace%v
            workspace%r_hat = workspace%r
            workspace%v = 0.0_prec
            workspace%p = 0.0_prec

            rhs_norm = euclidean_norm(rhs)
            residual_norm = euclidean_norm(workspace%r)
            shadow_norm = residual_norm
            residual_limit = options%absolute_tolerance + options%relative_tolerance*rhs_norm
            result%initial_residual_norm = residual_norm
            call update_result_residual(result, residual_norm, rhs_norm)

            if (.not. ieee_is_finite(residual_norm)) then
                call set_result_failure(result, BICGSTAB_NONFINITE_VALUE, 0, &
                                        'initial residual is not finite')
                return
            end if
            if (residual_norm <= residual_limit) then
                call set_result_success(result, 0, residual_norm, rhs_norm)
                return
            end if

            rho_previous = 1.0_prec
            alpha = 1.0_prec
            omega = 1.0_prec

            do iteration = 1, options%max_iterations
                rho = dot_product(workspace%r_hat, workspace%r)
                if (scaled_dot_is_zero(rho, shadow_norm, residual_norm, &
                                       options%breakdown_tolerance)) then
                    call restart_recurrence(workspace, rho_previous, alpha, omega, residual_norm)
                    shadow_norm = residual_norm
                    rho = dot_product(workspace%r_hat, workspace%r)
                    if (abs(rho) <= tiny(1.0_prec) .or. .not. ieee_is_finite(rho)) then
                        call set_result_failure(result, BICGSTAB_RHO_BREAKDOWN, iteration - 1, &
                                                'rho breakdown after recurrence restart')
                        return
                    end if
                end if

                beta = (rho/rho_previous)*(alpha/omega)
                workspace%p = workspace%r + beta*(workspace%p - omega*workspace%v)

                call apply_effective_operator(operator, preconditioner, identity_selected, &
                    workspace%p, workspace%p_preconditioned, workspace%v, operation_status)
                if (operation_status /= BICGSTAB_SUCCESS) then
                    call set_component_failure_result(result, operation_status, iteration - 1)
                    return
                end if

                denominator = dot_product(workspace%r_hat, workspace%v)
                if (abs(denominator) <= tiny(1.0_prec) .or. &
                    .not. ieee_is_finite(denominator)) then
                    call set_result_failure(result, BICGSTAB_ALPHA_BREAKDOWN, iteration - 1, &
                                            'alpha denominator breakdown')
                    return
                end if

                alpha = rho/denominator
                workspace%s = workspace%r - alpha*workspace%v
                residual_norm = euclidean_norm(workspace%s)

                if (residual_norm <= residual_limit) then
                    if (identity_selected) then
                        solution = solution + alpha*workspace%p
                    else
                        solution = solution + alpha*workspace%p_preconditioned
                    end if
                    call verify_true_residual(operator, rhs, solution, workspace, &
                                              rhs_norm, residual_limit, iteration, result, verified)
                    if (verified) return
                    call restart_recurrence(workspace, rho_previous, alpha, omega, residual_norm)
                    shadow_norm = residual_norm
                    cycle
                end if

                call apply_effective_operator(operator, preconditioner, identity_selected, &
                    workspace%s, workspace%s_preconditioned, workspace%t, operation_status)
                if (operation_status /= BICGSTAB_SUCCESS) then
                    call set_component_failure_result(result, operation_status, iteration - 1)
                    return
                end if

                t_dot_t = dot_product(workspace%t, workspace%t)
                if (t_dot_t <= tiny(1.0_prec) .or. .not. ieee_is_finite(t_dot_t)) then
                    call set_result_failure(result, BICGSTAB_OMEGA_BREAKDOWN, iteration - 1, &
                                            'omega denominator breakdown')
                    return
                end if

                omega = dot_product(workspace%t, workspace%s)/t_dot_t
                if (abs(omega) <= options%breakdown_tolerance .or. &
                    .not. ieee_is_finite(omega)) then
                    call set_result_failure(result, BICGSTAB_OMEGA_BREAKDOWN, iteration - 1, &
                                            'omega breakdown')
                    return
                end if

                if (identity_selected) then
                    solution = solution + alpha*workspace%p + omega*workspace%s
                else
                    solution = solution + alpha*workspace%p_preconditioned + &
                               omega*workspace%s_preconditioned
                end if

                workspace%r = workspace%s - omega*workspace%t
                residual_norm = euclidean_norm(workspace%r)
                call update_result_residual(result, residual_norm, rhs_norm)
                result%iterations = iteration

                if (.not. ieee_is_finite(residual_norm)) then
                    call set_result_failure(result, BICGSTAB_NONFINITE_VALUE, iteration, &
                                            'nonfinite value generated')
                    return
                end if

                if (residual_norm <= residual_limit) then
                    call verify_true_residual(operator, rhs, solution, workspace, &
                                              rhs_norm, residual_limit, iteration, result, verified)
                    if (verified) return
                    call restart_recurrence(workspace, rho_previous, alpha, omega, residual_norm)
                    shadow_norm = residual_norm
                    cycle
                end if

                rho_previous = rho
            end do

            call verify_true_residual(operator, rhs, solution, workspace, rhs_norm, &
                                      residual_limit, options%max_iterations, result, verified)
            if (.not. verified) then
                call set_result_failure(result, BICGSTAB_MAX_ITERATIONS, &
                                        options%max_iterations, 'maximum iterations reached')
            end if
        end subroutine bicgstab_core



        ! Solver support routines

        logical function inputs_are_valid(operator, rhs, solution, options, n)
            class(linear_operator_t), intent(in) :: operator
            real(prec), intent(in) :: rhs(:), solution(:)
            type(bicgstab_options_t), intent(in) :: options
            integer, intent(in) :: n

            inputs_are_valid = operator%is_valid() .and. n > 0
            inputs_are_valid = inputs_are_valid .and. size(rhs) == n .and. size(solution) == n
            inputs_are_valid = inputs_are_valid .and. options%max_iterations >= 0
            inputs_are_valid = inputs_are_valid .and. options%relative_tolerance >= 0.0_prec
            inputs_are_valid = inputs_are_valid .and. options%absolute_tolerance >= 0.0_prec
            inputs_are_valid = inputs_are_valid .and. options%breakdown_tolerance >= 0.0_prec
            inputs_are_valid = inputs_are_valid .and. all(ieee_is_finite(rhs))
            inputs_are_valid = inputs_are_valid .and. all(ieee_is_finite(solution))
        end function inputs_are_valid

        subroutine apply_effective_operator(operator, preconditioner, identity_selected, &
                                            input, preconditioned_input, output, status)
            class(linear_operator_t), intent(in) :: operator
            class(preconditioner_t), intent(inout) :: preconditioner
            logical, intent(in) :: identity_selected
            real(prec), intent(in) :: input(:)
            real(prec), allocatable, intent(inout) :: preconditioned_input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            integer :: component_status

            if (identity_selected) then
                call operator%apply(input, output, component_status)
                if (component_status == APPLY_SUCCESS) then
                    status = BICGSTAB_SUCCESS
                else
                    status = BICGSTAB_OPERATOR_FAILURE
                end if
                return
            end if

            if (.not. allocated(preconditioned_input)) then
                output = 0.0_prec
                status = BICGSTAB_PRECONDITIONER_FAILURE
                return
            end if

            call preconditioner%apply(input, preconditioned_input, component_status)
            if (component_status /= APPLY_SUCCESS) then
                output = 0.0_prec
                status = BICGSTAB_PRECONDITIONER_FAILURE
                return
            end if

            call operator%apply(preconditioned_input, output, component_status)
            if (component_status == APPLY_SUCCESS) then
                status = BICGSTAB_SUCCESS
            else
                status = BICGSTAB_OPERATOR_FAILURE
            end if
        end subroutine apply_effective_operator

        subroutine verify_true_residual(operator, rhs, solution, workspace, rhs_norm, &
                                        residual_limit, iteration, result, verified)
            class(linear_operator_t), intent(in) :: operator
            real(prec), intent(in) :: rhs(:), solution(:)
            type(bicgstab_workspace_t), intent(inout) :: workspace
            real(prec), intent(in) :: rhs_norm, residual_limit
            integer, intent(in) :: iteration
            type(bicgstab_result_t), intent(inout) :: result
            logical, intent(out) :: verified

            integer :: status
            real(prec) :: residual_norm

            call operator%apply(solution, workspace%v, status)
            if (status /= APPLY_SUCCESS) then
                call set_result_failure(result, BICGSTAB_OPERATOR_FAILURE, iteration, &
                                        'operator failed during residual verification')
                verified = .true.
                return
            end if

            workspace%r = rhs - workspace%v
            residual_norm = euclidean_norm(workspace%r)
            call update_result_residual(result, residual_norm, rhs_norm)
            result%iterations = iteration

            if (.not. ieee_is_finite(residual_norm)) then
                call set_result_failure(result, BICGSTAB_NONFINITE_VALUE, iteration, &
                                        'true residual is not finite')
                verified = .true.
            elseif (residual_norm <= residual_limit) then
                call set_result_success(result, iteration, residual_norm, rhs_norm)
                verified = .true.
            else
                verified = .false.
            end if
        end subroutine verify_true_residual

        subroutine restart_recurrence(workspace, rho_previous, alpha, omega, residual_norm)
            type(bicgstab_workspace_t), intent(inout) :: workspace
            real(prec), intent(out) :: rho_previous, alpha, omega, residual_norm

            workspace%r_hat = workspace%r
            workspace%p = 0.0_prec
            workspace%v = 0.0_prec
            rho_previous = 1.0_prec
            alpha = 1.0_prec
            omega = 1.0_prec
            residual_norm = euclidean_norm(workspace%r)
        end subroutine restart_recurrence

        logical function scaled_dot_is_zero(value, first_norm, second_norm, tolerance)
            real(prec), intent(in) :: value, first_norm, second_norm, tolerance

            real(prec) :: scale

            scale = first_norm*second_norm
            scaled_dot_is_zero = abs(value) <= max(tiny(1.0_prec), tolerance*scale)
            scaled_dot_is_zero = scaled_dot_is_zero .or. .not. ieee_is_finite(value)
        end function scaled_dot_is_zero

        pure real(prec) function euclidean_norm(vector)
            real(prec), intent(in) :: vector(:)

            euclidean_norm = sqrt(max(0.0_prec, dot_product(vector, vector)))
        end function euclidean_norm

        subroutine initialize_result(result)
            type(bicgstab_result_t), intent(out) :: result

            result%converged = .false.
            result%status = BICGSTAB_INVALID_INPUT
            result%iterations = 0
            result%initial_residual_norm = huge(1.0_prec)
            result%residual_norm = huge(1.0_prec)
            result%relative_residual = huge(1.0_prec)
            result%reason = 'solver not started'
        end subroutine initialize_result

        subroutine set_result_success(result, iterations, residual_norm, rhs_norm)
            type(bicgstab_result_t), intent(inout) :: result
            integer, intent(in) :: iterations
            real(prec), intent(in) :: residual_norm, rhs_norm

            result%converged = .true.
            result%status = BICGSTAB_SUCCESS
            result%iterations = iterations
            result%reason = 'converged'
            call update_result_residual(result, residual_norm, rhs_norm)
        end subroutine set_result_success

        subroutine set_result_failure(result, status, iterations, reason)
            type(bicgstab_result_t), intent(inout) :: result
            integer, intent(in) :: status, iterations
            character(len=*), intent(in) :: reason

            result%converged = .false.
            result%status = status
            result%iterations = iterations
            result%reason = reason
        end subroutine set_result_failure

        subroutine set_component_failure_result(result, status, iterations)
            type(bicgstab_result_t), intent(inout) :: result
            integer, intent(in) :: status, iterations

            select case (status)
            case (BICGSTAB_PRECONDITIONER_FAILURE)
                call set_result_failure(result, status, iterations, &
                                        'preconditioner application failed')
            case (BICGSTAB_OPERATOR_FAILURE)
                call set_result_failure(result, status, iterations, &
                                        'operator application failed')
            case default
                call set_result_failure(result, BICGSTAB_INVALID_INPUT, iterations, &
                                        'unexpected component failure')
            end select
        end subroutine set_component_failure_result

        subroutine update_result_residual(result, residual_norm, rhs_norm)
            type(bicgstab_result_t), intent(inout) :: result
            real(prec), intent(in) :: residual_norm, rhs_norm

            result%residual_norm = residual_norm
            if (rhs_norm > tiny(1.0_prec)) then
                result%relative_residual = residual_norm/rhs_norm
            else
                result%relative_residual = residual_norm
            end if
        end subroutine update_result_residual



        ! Workspace management

        subroutine ensure_workspace_capacity(this, n, preconditioned, status)
            class(bicgstab_workspace_t), intent(inout) :: this
            integer, intent(in) :: n
            logical, intent(in) :: preconditioned
            integer, intent(out) :: status

            integer :: allocation_status

            status = APPLY_FAILURE
            if (n <= 0) return

            if (allocated(this%r)) then
                if (size(this%r) /= n) call this%clear()
            end if

            if (.not. allocated(this%r)) then
                allocate(this%r(n), this%r_hat(n), this%v(n), this%p(n), &
                         this%s(n), this%t(n), stat=allocation_status)
                if (allocation_status /= 0) then
                    call this%clear()
                    return
                end if
            end if

            if (preconditioned) then
                if (.not. allocated(this%p_preconditioned)) then
                    allocate(this%p_preconditioned(n), this%s_preconditioned(n), &
                             stat=allocation_status)
                    if (allocation_status /= 0) then
                        call this%clear()
                        return
                    end if
                elseif (size(this%p_preconditioned) /= n) then
                    deallocate(this%p_preconditioned, this%s_preconditioned)
                    allocate(this%p_preconditioned(n), this%s_preconditioned(n), &
                             stat=allocation_status)
                    if (allocation_status /= 0) then
                        call this%clear()
                        return
                    end if
                end if
            end if

            status = APPLY_SUCCESS
        end subroutine ensure_workspace_capacity

        subroutine clear_workspace(this)
            class(bicgstab_workspace_t), intent(inout) :: this

            if (allocated(this%r)) deallocate(this%r)
            if (allocated(this%r_hat)) deallocate(this%r_hat)
            if (allocated(this%v)) deallocate(this%v)
            if (allocated(this%p)) deallocate(this%p)
            if (allocated(this%s)) deallocate(this%s)
            if (allocated(this%t)) deallocate(this%t)
            if (allocated(this%p_preconditioned)) deallocate(this%p_preconditioned)
            if (allocated(this%s_preconditioned)) deallocate(this%s_preconditioned)
        end subroutine clear_workspace



        ! Matrix-free banded operator

        subroutine bind_tridiagonal(this, lower, diagonal, upper, status)
            class(banded_operator_t), intent(out) :: this
            real(prec), target, contiguous, intent(in) :: lower(:), diagonal(:), upper(:)
            integer, optional, intent(out) :: status

            this%n = size(diagonal)
            this%number_of_pairs = 1
            this%offsets(1) = 1
            this%diagonal => diagonal
            this%lower(1)%values => lower
            this%upper(1)%values => upper
            this%valid_layout = coefficient_sizes_are_valid(this)

            if (present(status)) status = merge(APPLY_SUCCESS, APPLY_FAILURE, this%valid_layout)
        end subroutine bind_tridiagonal

        subroutine bind_pentadiagonal(this, offset_1, offset_2, lower_1, diagonal, &
                                      upper_1, lower_2, upper_2, status)
            class(banded_operator_t), intent(out) :: this
            integer, intent(in) :: offset_1, offset_2
            real(prec), target, contiguous, intent(in) :: lower_1(:), diagonal(:), upper_1(:)
            real(prec), target, contiguous, intent(in) :: lower_2(:), upper_2(:)
            integer, optional, intent(out) :: status

            this%n = size(diagonal)
            this%number_of_pairs = 2
            this%offsets(1:2) = [offset_1, offset_2]
            this%diagonal => diagonal
            this%lower(1)%values => lower_1
            this%upper(1)%values => upper_1
            this%lower(2)%values => lower_2
            this%upper(2)%values => upper_2
            this%valid_layout = coefficient_sizes_are_valid(this)

            if (present(status)) status = merge(APPLY_SUCCESS, APPLY_FAILURE, this%valid_layout)
        end subroutine bind_pentadiagonal

        subroutine bind_heptadiagonal(this, offset_1, offset_2, offset_3, &
                                      lower_1, diagonal, upper_1, lower_2, upper_2, &
                                      lower_3, upper_3, status)
            class(banded_operator_t), intent(out) :: this
            integer, intent(in) :: offset_1, offset_2, offset_3
            real(prec), target, contiguous, intent(in) :: lower_1(:), diagonal(:), upper_1(:)
            real(prec), target, contiguous, intent(in) :: lower_2(:), upper_2(:)
            real(prec), target, contiguous, intent(in) :: lower_3(:), upper_3(:)
            integer, optional, intent(out) :: status

            this%n = size(diagonal)
            this%number_of_pairs = 3
            this%offsets = [offset_1, offset_2, offset_3]
            this%diagonal => diagonal
            this%lower(1)%values => lower_1
            this%upper(1)%values => upper_1
            this%lower(2)%values => lower_2
            this%upper(2)%values => upper_2
            this%lower(3)%values => lower_3
            this%upper(3)%values => upper_3
            this%valid_layout = coefficient_sizes_are_valid(this)

            if (present(status)) status = merge(APPLY_SUCCESS, APPLY_FAILURE, this%valid_layout)
        end subroutine bind_heptadiagonal

        subroutine bind_eqtri(this, n_x, a, b, c, d, e, status)
            class(banded_operator_t), intent(out) :: this
            integer, intent(in) :: n_x
            real(prec), target, contiguous, intent(in) :: a(:), b(:), c(:), d(:), e(:)
            integer, optional, intent(out) :: status

            integer :: local_status

            call this%bind_pentadiagonal(n_x, n_x + 1, a, b, c, d, e, local_status)
            if (present(status)) status = local_status
        end subroutine bind_eqtri

        subroutine bind_triiso(this, n_x, a, b, c, d, e, status)
            class(banded_operator_t), intent(out) :: this
            integer, intent(in) :: n_x
            real(prec), target, contiguous, intent(in) :: a(:), b(:), c(:), d(:), e(:)
            integer, optional, intent(out) :: status

            integer :: local_status

            call this%bind_pentadiagonal(1, n_x, a, b, c, d, e, local_status)
            if (present(status)) status = local_status
        end subroutine bind_triiso

        subroutine apply_banded_operator(this, input, output, status)
            class(banded_operator_t), intent(in) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            integer :: pair, offset, row

            if (.not. this%valid_layout .or. size(input) /= this%n .or. &
                size(output) /= this%n) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if

            output = this%diagonal*input

            do pair = 1, this%number_of_pairs
                offset = this%offsets(pair)
                do row = offset + 1, this%n
                    output(row) = output(row) + &
                        this%lower(pair)%values(row)*input(row - offset)
                end do
                do row = 1, this%n - offset
                    output(row) = output(row) + &
                        this%upper(pair)%values(row)*input(row + offset)
                end do
            end do

            status = APPLY_SUCCESS
        end subroutine apply_banded_operator

        integer function banded_operator_order(this)
            class(banded_operator_t), intent(in) :: this

            if (this%is_valid()) then
                banded_operator_order = this%n
            else
                banded_operator_order = 0
            end if
        end function banded_operator_order

        logical function banded_operator_is_valid(this)
            class(banded_operator_t), intent(in) :: this

            banded_operator_is_valid = this%valid_layout .and. coefficient_sizes_are_valid(this)
        end function banded_operator_is_valid

        logical function coefficient_sizes_are_valid(this)
            class(banded_operator_t), intent(in) :: this

            integer :: pair

            coefficient_sizes_are_valid = this%n > 0 .and. associated(this%diagonal)
            if (.not. coefficient_sizes_are_valid) return

            if (size(this%diagonal) /= this%n) then
                coefficient_sizes_are_valid = .false.
                return
            end if

            do pair = 1, this%number_of_pairs
                if (.not. associated(this%lower(pair)%values) .or. &
                    .not. associated(this%upper(pair)%values)) then
                    coefficient_sizes_are_valid = .false.
                    return
                end if
                if (size(this%lower(pair)%values) /= this%n .or. &
                    size(this%upper(pair)%values) /= this%n) then
                    coefficient_sizes_are_valid = .false.
                    return
                end if
                if (this%offsets(pair) <= 0 .or. this%offsets(pair) >= this%n) then
                    coefficient_sizes_are_valid = .false.
                    return
                end if
            end do

            do pair = 2, this%number_of_pairs
                if (this%offsets(pair) <= this%offsets(pair - 1)) then
                    coefficient_sizes_are_valid = .false.
                    return
                end if
            end do
        end function coefficient_sizes_are_valid



        ! Preconditioners

        subroutine apply_identity_preconditioner(this, input, output, status)
            class(identity_preconditioner_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            if (size(input) /= size(output)) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if

            output = input
            status = APPLY_SUCCESS
        end subroutine apply_identity_preconditioner

        subroutine setup_jacobi_preconditioner(this, operator, status, pivot_tolerance)
            class(jacobi_preconditioner_t), intent(out) :: this
            type(banded_operator_t), intent(in) :: operator
            integer, intent(out) :: status
            real(prec), optional, intent(in) :: pivot_tolerance

            real(prec) :: threshold

            status = APPLY_FAILURE
            if (.not. operator%is_valid()) return

            threshold = default_pivot_tolerance(operator%diagonal)
            if (present(pivot_tolerance)) threshold = max(0.0_prec, pivot_tolerance)
            if (any(abs(operator%diagonal) <= threshold)) return

            allocate(this%inverse_diagonal(operator%n))
            this%inverse_diagonal = 1.0_prec/operator%diagonal
            this%ready = .true.
            status = APPLY_SUCCESS
        end subroutine setup_jacobi_preconditioner

        subroutine apply_jacobi_preconditioner(this, input, output, status)
            class(jacobi_preconditioner_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            if (.not. this%ready .or. .not. allocated(this%inverse_diagonal)) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if
            if (size(input) /= size(this%inverse_diagonal) .or. &
                size(output) /= size(this%inverse_diagonal)) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if

            output = this%inverse_diagonal*input
            status = APPLY_SUCCESS
        end subroutine apply_jacobi_preconditioner

        subroutine setup_sgs_preconditioner(this, operator, status, pivot_tolerance)
            class(symmetric_gauss_seidel_preconditioner_t), intent(out) :: this
            type(banded_operator_t), intent(in) :: operator
            integer, intent(out) :: status
            real(prec), optional, intent(in) :: pivot_tolerance

            real(prec) :: threshold

            status = APPLY_FAILURE
            if (.not. operator%is_valid()) return

            threshold = default_pivot_tolerance(operator%diagonal)
            if (present(pivot_tolerance)) threshold = max(0.0_prec, pivot_tolerance)
            if (any(abs(operator%diagonal) <= threshold)) return

            this%operator_view = operator
            this%ready = .true.
            status = APPLY_SUCCESS
        end subroutine setup_sgs_preconditioner

        subroutine apply_sgs_preconditioner(this, input, output, status)
            class(symmetric_gauss_seidel_preconditioner_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            integer :: row, pair, offset
            real(prec) :: value

            if (.not. this%ready .or. .not. this%operator_view%is_valid()) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if
            if (size(input) /= this%operator_view%n .or. &
                size(output) /= this%operator_view%n) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if

            do row = 1, this%operator_view%n
                value = input(row)
                do pair = 1, this%operator_view%number_of_pairs
                    offset = this%operator_view%offsets(pair)
                    if (row > offset) value = value - &
                        this%operator_view%lower(pair)%values(row)*output(row - offset)
                end do
                output(row) = value/this%operator_view%diagonal(row)
            end do

            output = this%operator_view%diagonal*output

            do row = this%operator_view%n, 1, -1
                value = output(row)
                do pair = 1, this%operator_view%number_of_pairs
                    offset = this%operator_view%offsets(pair)
                    if (row <= this%operator_view%n - offset) value = value - &
                        this%operator_view%upper(pair)%values(row)*output(row + offset)
                end do
                output(row) = value/this%operator_view%diagonal(row)
            end do

            status = APPLY_SUCCESS
        end subroutine apply_sgs_preconditioner

        subroutine setup_ilu0_preconditioner(this, operator, status, pivot_tolerance)
            class(ilu0_preconditioner_t), intent(out) :: this
            type(banded_operator_t), intent(in) :: operator
            integer, intent(out) :: status
            real(prec), optional, intent(in) :: pivot_tolerance

            integer :: row, pair, previous_pair, column, previous_row
            real(prec) :: value, threshold

            status = APPLY_FAILURE
            if (.not. operator%is_valid()) return

            this%n = operator%n
            this%number_of_pairs = operator%number_of_pairs
            this%offsets = operator%offsets
            allocate(this%lower_factor(this%n, this%number_of_pairs), &
                     this%diagonal_factor(this%n), &
                     this%upper_factor(this%n, this%number_of_pairs))
            this%lower_factor = 0.0_prec
            this%diagonal_factor = 0.0_prec
            this%upper_factor = 0.0_prec

            threshold = default_pivot_tolerance(operator%diagonal)
            if (present(pivot_tolerance)) threshold = max(0.0_prec, pivot_tolerance)

            do row = 1, this%n
                do pair = this%number_of_pairs, 1, -1
                    column = row - this%offsets(pair)
                    if (column < 1) cycle

                    value = operator%lower(pair)%values(row)
                    do previous_pair = this%number_of_pairs, 1, -1
                        previous_row = row - this%offsets(previous_pair)
                        if (previous_row < 1 .or. previous_row >= column) cycle
                        value = value - this%lower_factor(row, previous_pair)* &
                            ilu_upper_value(this, previous_row, column)
                    end do

                    if (abs(this%diagonal_factor(column)) <= threshold) return
                    this%lower_factor(row, pair) = value/this%diagonal_factor(column)
                end do

                value = operator%diagonal(row)
                do previous_pair = 1, this%number_of_pairs
                    previous_row = row - this%offsets(previous_pair)
                    if (previous_row < 1) cycle
                    value = value - this%lower_factor(row, previous_pair)* &
                        ilu_upper_value(this, previous_row, row)
                end do
                if (abs(value) <= threshold .or. .not. ieee_is_finite(value)) return
                this%diagonal_factor(row) = value

                do pair = 1, this%number_of_pairs
                    column = row + this%offsets(pair)
                    if (column > this%n) cycle

                    value = operator%upper(pair)%values(row)
                    do previous_pair = 1, this%number_of_pairs
                        previous_row = row - this%offsets(previous_pair)
                        if (previous_row < 1) cycle
                        value = value - this%lower_factor(row, previous_pair)* &
                            ilu_upper_value(this, previous_row, column)
                    end do
                    this%upper_factor(row, pair) = value
                end do
            end do

            this%ready = .true.
            status = APPLY_SUCCESS
        end subroutine setup_ilu0_preconditioner

        pure real(prec) function ilu_upper_value(this, row, column)
            class(ilu0_preconditioner_t), intent(in) :: this
            integer, intent(in) :: row, column

            integer :: pair, difference

            ilu_upper_value = 0.0_prec
            if (row < 1 .or. row > this%n .or. column < row .or. column > this%n) return
            if (column == row) then
                ilu_upper_value = this%diagonal_factor(row)
                return
            end if

            difference = column - row
            do pair = 1, this%number_of_pairs
                if (difference == this%offsets(pair)) then
                    ilu_upper_value = this%upper_factor(row, pair)
                    return
                end if
            end do
        end function ilu_upper_value

        subroutine apply_ilu0_preconditioner(this, input, output, status)
            class(ilu0_preconditioner_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            integer :: row, pair, offset
            real(prec) :: value

            if (.not. this%ready .or. size(input) /= this%n .or. size(output) /= this%n) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if

            do row = 1, this%n
                value = input(row)
                do pair = 1, this%number_of_pairs
                    offset = this%offsets(pair)
                    if (row > offset) value = value - &
                        this%lower_factor(row, pair)*output(row - offset)
                end do
                output(row) = value
            end do

            do row = this%n, 1, -1
                value = output(row)
                do pair = 1, this%number_of_pairs
                    offset = this%offsets(pair)
                    if (row <= this%n - offset) value = value - &
                        this%upper_factor(row, pair)*output(row + offset)
                end do
                output(row) = value/this%diagonal_factor(row)
            end do

            status = APPLY_SUCCESS
        end subroutine apply_ilu0_preconditioner

        subroutine bind_user_preconditioner(this, callback, context)
            class(user_preconditioner_t), intent(out) :: this
            procedure(user_preconditioner_callback) :: callback
            class(*), target, intent(inout) :: context

            this%callback => callback
            this%context => context
        end subroutine bind_user_preconditioner

        subroutine apply_user_preconditioner(this, input, output, status)
            class(user_preconditioner_t), intent(inout) :: this
            real(prec), intent(in) :: input(:)
            real(prec), intent(out) :: output(:)
            integer, intent(out) :: status

            if (.not. associated(this%callback) .or. .not. associated(this%context)) then
                output = 0.0_prec
                status = APPLY_FAILURE
                return
            end if

            call this%callback(this%context, input, output, status)
        end subroutine apply_user_preconditioner

        pure real(prec) function default_pivot_tolerance(diagonal)
            real(prec), intent(in) :: diagonal(:)

            default_pivot_tolerance = max(tiny(1.0_prec), &
                100.0_prec*epsilon(1.0_prec)*maxval(abs(diagonal)))
        end function default_pivot_tolerance

end module BiCGSTAB_Solver
