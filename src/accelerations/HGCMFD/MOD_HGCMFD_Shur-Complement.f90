module HGCMFD_Shur_Complement
    
    use HGCMFD_Variables, only: dprec
        implicit none
    
    contains
    
        subroutine Schur_Complement_Vect(n, a_in, b_in, c_in, y_in, x) 
            integer, intent(in) :: n
            integer :: i
            real(dprec), intent(in) :: a_in(:), b_in(:), c_in(:), y_in(:)        ! Ax = y
            real(dprec), intent(out) :: x(:)
            real(dprec) :: factor, S_compl
            
            S_compl=b_in(1)
            x(1)=y_in(1)
            do i=2, n
                S_compl=S_compl-c_in(i)*a_in(i)/b_in(i)
                x(1)=x(1)-c_in(i)*y_in(i)/b_in(i)
            end do
            
            x(1)=x(1)/S_compl
            x(2:n)=(y_in(2:n)-a_in(2:n)*x(1))/b_in(2:n)
            
            return
            
        end subroutine Schur_Complement_Vect
    
end module HGCMFD_Shur_Complement