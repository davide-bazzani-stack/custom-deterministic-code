module HGCMFD_Gauss_Elimination
    
    use HGCMFD_Variables, only: dprec 
        implicit none
    
    contains
    
    subroutine GE_Expl_Main(n, A_in, b_in, x) 
            implicit none
        integer, intent(in) :: n
        integer :: i, j, k
        real(dprec), intent(in) :: A_in(:,:), b_in(:)
        real(dprec), intent(out) :: x(:)
        real(dprec) :: factor
        real(dprec), allocatable :: A(:,:), b(:)
    
        allocate(A(n,n))
        allocate(b(n)) 
        A=A_in
        b=b_in
        
        do k = 1, n-1
          ! Find the pivot row
          call partial_pivot(A, b, n, k)

          ! Perform elimination
          do i = k + 1, n
            factor = A(i, k) / A(k, k)
            A(i, k:n) = A(i, k:n) - factor * A(k, k:n)
            b(i) = b(i) - factor * b(k)
          end do
        end do

        ! Back substitution
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
          x(i) = (b(i) - sum(A(i, i+1:n) * x(i+1:n))) / A(i, i)
        end do
    end subroutine GE_Expl_Main
    
    
    
    subroutine partial_pivot(A, b, n, k) 
        integer, intent(in) :: n, k
        real(dprec), intent(inout) :: A(n, n), b(n)
        integer :: max_row, i
        real(dprec) :: temp
    
        ! Find the row with the largest element in column k
        max_row = k
        do i = k + 1, n
            if (abs(A(i, k)) > abs(A(max_row, k))) then
                max_row = i
            end if
        end do
    
        ! Swap rows if needed
        if (max_row /= k) then
            A(k, :) = A(k, :) + A(max_row, :)
            A(max_row, :) = A(k, :) - A(max_row, :)
            A(k, :) = A(k, :) - A(max_row, :)
    
            temp = b(k)
            b(k) = b(max_row)
            b(max_row) = temp
        end if
    end subroutine partial_pivot
    
end module HGCMFD_Gauss_Elimination