module LAPACK_Interfaces
    use, intrinsic :: iso_fortran_env, only : real32, real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    implicit none

    private

    ! These checked generic wrappers intentionally support LAPACK's real32
    ! (single precision) and real64 (double precision) routines only.
    integer, parameter, public :: LAPACK_WRAPPER_ALLOCATION_FAILURE = -1000
    integer, parameter, public :: LAPACK_WRAPPER_INVALID_WORKSPACE_QUERY = -1001
    integer, parameter, public :: LAPACK_WRAPPER_INVALID_INPUT      = -1002

    interface lapack_getrf
        module procedure lapack_getrf_real32
        module procedure lapack_getrf_real64
    end interface lapack_getrf

    interface lapack_getrs
        module procedure lapack_getrs_real32
        module procedure lapack_getrs_real64
    end interface lapack_getrs

    interface lapack_gecon
        module procedure lapack_gecon_real32
        module procedure lapack_gecon_real64
    end interface lapack_gecon

    interface lapack_geev
        module procedure lapack_geev_real32
        module procedure lapack_geev_real64
    end interface lapack_geev

    interface lapack_ggev
        module procedure lapack_ggev_real32
        module procedure lapack_ggev_real64
    end interface lapack_ggev

    public :: lapack_getrf, lapack_getrs, lapack_gecon
    public :: lapack_geev, lapack_ggev

    contains

        subroutine lapack_getrf_real32(matrix, pivots, info)
            real(real32), intent(inout) :: matrix(:,:)
            integer, intent(out) :: pivots(:)
            integer, intent(out) :: info

            interface
                subroutine sgetrf(m, n, a, lda, ipiv, lapack_info)
                    import :: real32
                    integer, intent(in) :: m, n, lda
                    real(real32), intent(inout) :: a(lda,*)
                    integer, intent(out) :: ipiv(*)
                    integer, intent(out) :: lapack_info
                end subroutine sgetrf
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            if (size(matrix, 1) < 1 .or. size(matrix, 2) < 1) return
            if (size(pivots) < min(size(matrix, 1), size(matrix, 2))) return

            call sgetrf(size(matrix, 1), size(matrix, 2), matrix, &
                        size(matrix, 1), pivots, info)
        end subroutine lapack_getrf_real32

        subroutine lapack_getrf_real64(matrix, pivots, info)
            real(real64), intent(inout) :: matrix(:,:)
            integer, intent(out) :: pivots(:)
            integer, intent(out) :: info

            interface
                subroutine dgetrf(m, n, a, lda, ipiv, lapack_info)
                    import :: real64
                    integer, intent(in) :: m, n, lda
                    real(real64), intent(inout) :: a(lda,*)
                    integer, intent(out) :: ipiv(*)
                    integer, intent(out) :: lapack_info
                end subroutine dgetrf
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            if (size(matrix, 1) < 1 .or. size(matrix, 2) < 1) return
            if (size(pivots) < min(size(matrix, 1), size(matrix, 2))) return

            call dgetrf(size(matrix, 1), size(matrix, 2), matrix, &
                        size(matrix, 1), pivots, info)
        end subroutine lapack_getrf_real64

        subroutine lapack_getrs_real32(lu_matrix, pivots, right_hand_sides, info)
            real(real32), intent(in) :: lu_matrix(:,:)
            integer, intent(in) :: pivots(:)
            real(real32), intent(inout) :: right_hand_sides(:,:)
            integer, intent(out) :: info

            interface
                subroutine sgetrs(transpose_mode, n, number_of_rhs, a, lda, ipiv, &
                                  b, ldb, lapack_info)
                    import :: real32
                    character(len=1), intent(in) :: transpose_mode
                    integer, intent(in) :: n, number_of_rhs, lda, ldb
                    real(real32), intent(in) :: a(lda,*)
                    integer, intent(in) :: ipiv(*)
                    real(real32), intent(inout) :: b(ldb,*)
                    integer, intent(out) :: lapack_info
                end subroutine sgetrs
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            if (size(lu_matrix, 1) < 1) return
            if (size(lu_matrix, 2) /= size(lu_matrix, 1)) return
            if (size(pivots) < size(lu_matrix, 1)) return
            if (size(right_hand_sides, 1) /= size(lu_matrix, 1)) return
            if (size(right_hand_sides, 2) < 1) return

            call sgetrs('N', size(lu_matrix, 1), size(right_hand_sides, 2), &
                        lu_matrix, size(lu_matrix, 1), pivots, &
                        right_hand_sides, size(right_hand_sides, 1), info)
        end subroutine lapack_getrs_real32

        subroutine lapack_getrs_real64(lu_matrix, pivots, right_hand_sides, info)
            real(real64), intent(in) :: lu_matrix(:,:)
            integer, intent(in) :: pivots(:)
            real(real64), intent(inout) :: right_hand_sides(:,:)
            integer, intent(out) :: info

            interface
                subroutine dgetrs(transpose_mode, n, number_of_rhs, a, lda, ipiv, &
                                  b, ldb, lapack_info)
                    import :: real64
                    character(len=1), intent(in) :: transpose_mode
                    integer, intent(in) :: n, number_of_rhs, lda, ldb
                    real(real64), intent(in) :: a(lda,*)
                    integer, intent(in) :: ipiv(*)
                    real(real64), intent(inout) :: b(ldb,*)
                    integer, intent(out) :: lapack_info
                end subroutine dgetrs
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            if (size(lu_matrix, 1) < 1) return
            if (size(lu_matrix, 2) /= size(lu_matrix, 1)) return
            if (size(pivots) < size(lu_matrix, 1)) return
            if (size(right_hand_sides, 1) /= size(lu_matrix, 1)) return
            if (size(right_hand_sides, 2) < 1) return

            call dgetrs('N', size(lu_matrix, 1), size(right_hand_sides, 2), &
                        lu_matrix, size(lu_matrix, 1), pivots, &
                        right_hand_sides, size(right_hand_sides, 1), info)
        end subroutine lapack_getrs_real64

        subroutine lapack_gecon_real32(lu_matrix, matrix_norm, reciprocal_condition, info)
            real(real32), intent(in) :: lu_matrix(:,:)
            real(real32), intent(in) :: matrix_norm
            real(real32), intent(out) :: reciprocal_condition
            integer, intent(out) :: info

            real(real32), allocatable :: work(:)
            integer, allocatable :: integer_work(:)
            integer :: allocation_status, n

            interface
                subroutine sgecon(norm_type, n, a, lda, anorm, rcond, work, iwork, &
                                  lapack_info)
                    import :: real32
                    character(len=1), intent(in) :: norm_type
                    integer, intent(in) :: n, lda
                    real(real32), intent(in) :: a(lda,*), anorm
                    real(real32), intent(out) :: rcond, work(*)
                    integer, intent(out) :: iwork(*), lapack_info
                end subroutine sgecon
            end interface

            reciprocal_condition = 0.0_real32
            info = LAPACK_WRAPPER_INVALID_INPUT
            n = size(lu_matrix, 1)
            if (n < 1) return
            if (size(lu_matrix, 2) /= n) return
            if (.not. ieee_is_finite(matrix_norm) .or. matrix_norm < 0.0_real32) return
            if (n > (huge(n) - 3)/4) return

            allocate(work(4*n), integer_work(n), stat=allocation_status)
            if (allocation_status /= 0) then
                info = LAPACK_WRAPPER_ALLOCATION_FAILURE
                return
            end if

            call sgecon('1', n, lu_matrix, n, matrix_norm, &
                        reciprocal_condition, work, integer_work, info)
        end subroutine lapack_gecon_real32

        subroutine lapack_gecon_real64(lu_matrix, matrix_norm, reciprocal_condition, info)
            real(real64), intent(in) :: lu_matrix(:,:)
            real(real64), intent(in) :: matrix_norm
            real(real64), intent(out) :: reciprocal_condition
            integer, intent(out) :: info

            real(real64), allocatable :: work(:)
            integer, allocatable :: integer_work(:)
            integer :: allocation_status, n

            interface
                subroutine dgecon(norm_type, n, a, lda, anorm, rcond, work, iwork, &
                                  lapack_info)
                    import :: real64
                    character(len=1), intent(in) :: norm_type
                    integer, intent(in) :: n, lda
                    real(real64), intent(in) :: a(lda,*), anorm
                    real(real64), intent(out) :: rcond, work(*)
                    integer, intent(out) :: iwork(*), lapack_info
                end subroutine dgecon
            end interface

            reciprocal_condition = 0.0_real64
            info = LAPACK_WRAPPER_INVALID_INPUT
            n = size(lu_matrix, 1)
            if (n < 1) return
            if (size(lu_matrix, 2) /= n) return
            if (.not. ieee_is_finite(matrix_norm) .or. matrix_norm < 0.0_real64) return
            if (n > (huge(n) - 3)/4) return

            allocate(work(4*n), integer_work(n), stat=allocation_status)
            if (allocation_status /= 0) then
                info = LAPACK_WRAPPER_ALLOCATION_FAILURE
                return
            end if

            call dgecon('1', n, lu_matrix, n, matrix_norm, &
                        reciprocal_condition, work, integer_work, info)
        end subroutine lapack_gecon_real64

        subroutine lapack_geev_real32(matrix, eigenvalues_real, eigenvalues_imaginary, &
                                      right_eigenvectors, info)
            real(real32), intent(inout) :: matrix(:,:)
            real(real32), intent(out) :: eigenvalues_real(:), eigenvalues_imaginary(:)
            real(real32), intent(out) :: right_eigenvectors(:,:)
            integer, intent(out) :: info

            real(real32) :: left_eigenvectors(1,1), work_query(1)
            real(real32), allocatable :: work(:)
            integer :: allocation_status, lwork, n

            interface
                subroutine sgeev(job_left, job_right, n, a, lda, wr, wi, vl, ldvl, &
                                 vr, ldvr, work, lwork, lapack_info)
                    import :: real32
                    character(len=1), intent(in) :: job_left, job_right
                    integer, intent(in) :: n, lda, ldvl, ldvr, lwork
                    real(real32), intent(inout) :: a(lda,*)
                    real(real32), intent(out) :: wr(*), wi(*), vl(ldvl,*), &
                                                 vr(ldvr,*), work(*)
                    integer, intent(out) :: lapack_info
                end subroutine sgeev
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            n = size(matrix, 1)
            if (n < 1) return
            if (size(matrix, 2) /= n) return
            if (size(eigenvalues_real) /= n) return
            if (size(eigenvalues_imaginary) /= n) return
            if (size(right_eigenvectors, 1) /= n .or. &
                size(right_eigenvectors, 2) /= n) return

            lwork = -1
            call sgeev('N', 'V', n, matrix, n, eigenvalues_real, &
                       eigenvalues_imaginary, left_eigenvectors, 1, &
                       right_eigenvectors, n, work_query, lwork, info)
            if (info /= 0) return
            if (.not. ieee_is_finite(work_query(1)) .or. &
                work_query(1) < 1.0_real32 .or. &
                work_query(1) >= real(huge(lwork), real32)) then
                info = LAPACK_WRAPPER_INVALID_WORKSPACE_QUERY
                return
            end if

            lwork = max(1, ceiling(work_query(1)))
            allocate(work(lwork), stat=allocation_status)
            if (allocation_status /= 0) then
                info = LAPACK_WRAPPER_ALLOCATION_FAILURE
                return
            end if

            call sgeev('N', 'V', n, matrix, n, eigenvalues_real, &
                       eigenvalues_imaginary, left_eigenvectors, 1, &
                       right_eigenvectors, n, work, lwork, info)
        end subroutine lapack_geev_real32

        subroutine lapack_geev_real64(matrix, eigenvalues_real, eigenvalues_imaginary, &
                                      right_eigenvectors, info)
            real(real64), intent(inout) :: matrix(:,:)
            real(real64), intent(out) :: eigenvalues_real(:), eigenvalues_imaginary(:)
            real(real64), intent(out) :: right_eigenvectors(:,:)
            integer, intent(out) :: info

            real(real64) :: left_eigenvectors(1,1), work_query(1)
            real(real64), allocatable :: work(:)
            integer :: allocation_status, lwork, n

            interface
                subroutine dgeev(job_left, job_right, n, a, lda, wr, wi, vl, ldvl, &
                                 vr, ldvr, work, lwork, lapack_info)
                    import :: real64
                    character(len=1), intent(in) :: job_left, job_right
                    integer, intent(in) :: n, lda, ldvl, ldvr, lwork
                    real(real64), intent(inout) :: a(lda,*)
                    real(real64), intent(out) :: wr(*), wi(*), vl(ldvl,*), &
                                                 vr(ldvr,*), work(*)
                    integer, intent(out) :: lapack_info
                end subroutine dgeev
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            n = size(matrix, 1)
            if (n < 1) return
            if (size(matrix, 2) /= n) return
            if (size(eigenvalues_real) /= n) return
            if (size(eigenvalues_imaginary) /= n) return
            if (size(right_eigenvectors, 1) /= n .or. &
                size(right_eigenvectors, 2) /= n) return

            lwork = -1
            call dgeev('N', 'V', n, matrix, n, eigenvalues_real, &
                       eigenvalues_imaginary, left_eigenvectors, 1, &
                       right_eigenvectors, n, work_query, lwork, info)
            if (info /= 0) return
            if (.not. ieee_is_finite(work_query(1)) .or. &
                work_query(1) < 1.0_real64 .or. &
                work_query(1) >= real(huge(lwork), real64)) then
                info = LAPACK_WRAPPER_INVALID_WORKSPACE_QUERY
                return
            end if

            lwork = max(1, ceiling(work_query(1)))
            allocate(work(lwork), stat=allocation_status)
            if (allocation_status /= 0) then
                info = LAPACK_WRAPPER_ALLOCATION_FAILURE
                return
            end if

            call dgeev('N', 'V', n, matrix, n, eigenvalues_real, &
                       eigenvalues_imaginary, left_eigenvectors, 1, &
                       right_eigenvectors, n, work, lwork, info)
        end subroutine lapack_geev_real64

        subroutine lapack_ggev_real32(matrix_a, matrix_b, alpha_real, alpha_imaginary, &
                                      beta, right_eigenvectors, info)
            real(real32), intent(inout) :: matrix_a(:,:), matrix_b(:,:)
            real(real32), intent(out) :: alpha_real(:), alpha_imaginary(:), beta(:)
            real(real32), intent(out) :: right_eigenvectors(:,:)
            integer, intent(out) :: info

            real(real32) :: left_eigenvectors(1,1), work_query(1)
            real(real32), allocatable :: work(:)
            integer :: allocation_status, lwork, n

            interface
                subroutine sggev(job_left, job_right, n, a, lda, b, ldb, alphar, &
                                 alphai, beta, vl, ldvl, vr, ldvr, work, lwork, &
                                 lapack_info)
                    import :: real32
                    character(len=1), intent(in) :: job_left, job_right
                    integer, intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
                    real(real32), intent(inout) :: a(lda,*), b(ldb,*)
                    real(real32), intent(out) :: alphar(*), alphai(*), beta(*), &
                                                 vl(ldvl,*), vr(ldvr,*), work(*)
                    integer, intent(out) :: lapack_info
                end subroutine sggev
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            n = size(matrix_a, 1)
            if (n < 1) return
            if (size(matrix_a, 2) /= n) return
            if (size(matrix_b, 1) /= n .or. size(matrix_b, 2) /= n) return
            if (size(alpha_real) /= n) return
            if (size(alpha_imaginary) /= n) return
            if (size(beta) /= n) return
            if (size(right_eigenvectors, 1) /= n .or. &
                size(right_eigenvectors, 2) /= n) return

            lwork = -1
            call sggev('N', 'V', n, matrix_a, n, matrix_b, n, &
                       alpha_real, alpha_imaginary, beta, left_eigenvectors, 1, &
                       right_eigenvectors, n, work_query, lwork, info)
            if (info /= 0) return
            if (.not. ieee_is_finite(work_query(1)) .or. &
                work_query(1) < 1.0_real32 .or. &
                work_query(1) >= real(huge(lwork), real32)) then
                info = LAPACK_WRAPPER_INVALID_WORKSPACE_QUERY
                return
            end if

            lwork = max(1, ceiling(work_query(1)))
            allocate(work(lwork), stat=allocation_status)
            if (allocation_status /= 0) then
                info = LAPACK_WRAPPER_ALLOCATION_FAILURE
                return
            end if

            call sggev('N', 'V', n, matrix_a, n, matrix_b, n, &
                       alpha_real, alpha_imaginary, beta, left_eigenvectors, 1, &
                       right_eigenvectors, n, work, lwork, info)
        end subroutine lapack_ggev_real32

        subroutine lapack_ggev_real64(matrix_a, matrix_b, alpha_real, alpha_imaginary, &
                                      beta, right_eigenvectors, info)
            real(real64), intent(inout) :: matrix_a(:,:), matrix_b(:,:)
            real(real64), intent(out) :: alpha_real(:), alpha_imaginary(:), beta(:)
            real(real64), intent(out) :: right_eigenvectors(:,:)
            integer, intent(out) :: info

            real(real64) :: left_eigenvectors(1,1), work_query(1)
            real(real64), allocatable :: work(:)
            integer :: allocation_status, lwork, n

            interface
                subroutine dggev(job_left, job_right, n, a, lda, b, ldb, alphar, &
                                 alphai, beta, vl, ldvl, vr, ldvr, work, lwork, &
                                 lapack_info)
                    import :: real64
                    character(len=1), intent(in) :: job_left, job_right
                    integer, intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
                    real(real64), intent(inout) :: a(lda,*), b(ldb,*)
                    real(real64), intent(out) :: alphar(*), alphai(*), beta(*), &
                                                 vl(ldvl,*), vr(ldvr,*), work(*)
                    integer, intent(out) :: lapack_info
                end subroutine dggev
            end interface

            info = LAPACK_WRAPPER_INVALID_INPUT
            n = size(matrix_a, 1)
            if (n < 1) return
            if (size(matrix_a, 2) /= n) return
            if (size(matrix_b, 1) /= n .or. size(matrix_b, 2) /= n) return
            if (size(alpha_real) /= n) return
            if (size(alpha_imaginary) /= n) return
            if (size(beta) /= n) return
            if (size(right_eigenvectors, 1) /= n .or. &
                size(right_eigenvectors, 2) /= n) return

            lwork = -1
            call dggev('N', 'V', n, matrix_a, n, matrix_b, n, &
                       alpha_real, alpha_imaginary, beta, left_eigenvectors, 1, &
                       right_eigenvectors, n, work_query, lwork, info)
            if (info /= 0) return
            if (.not. ieee_is_finite(work_query(1)) .or. &
                work_query(1) < 1.0_real64 .or. &
                work_query(1) >= real(huge(lwork), real64)) then
                info = LAPACK_WRAPPER_INVALID_WORKSPACE_QUERY
                return
            end if

            lwork = max(1, ceiling(work_query(1)))
            allocate(work(lwork), stat=allocation_status)
            if (allocation_status /= 0) then
                info = LAPACK_WRAPPER_ALLOCATION_FAILURE
                return
            end if

            call dggev('N', 'V', n, matrix_a, n, matrix_b, n, &
                       alpha_real, alpha_imaginary, beta, left_eigenvectors, 1, &
                       right_eigenvectors, n, work, lwork, info)
        end subroutine lapack_ggev_real64

end module LAPACK_Interfaces
