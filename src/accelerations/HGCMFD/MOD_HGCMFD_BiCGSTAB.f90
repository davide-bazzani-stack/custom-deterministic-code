module HGCMFD_BiCGSTAB
    
    use HGCMFD_Variables, only: dprec
        implicit none
    
    contains
    
        subroutine BiCGSTAB_HGCMFD(n, max_iter, tol, a, b, c, rhs, x, info) 
            integer, intent(in) :: n, max_iter
            real(dprec), intent(in) :: a(:), b(:), c(:), rhs(:), tol
            real(dprec), intent(inout) :: x(:)
            integer, intent(out) :: info
            real(dprec), allocatable :: r(:), r_old(:), v(:), p(:), s(:), t(:), x_old(:)
            real(dprec) :: alpha, beta, omega, rho, rho_old, error, norm_b
            integer :: k
                    
            allocate(r(n))
            allocate(r_old(n))
            allocate(v(n))
            allocate(p(n))
            allocate(s(n))
            allocate(t(n))
            
            allocate(x_old(n))
            
            
            call matvect_arrow(n, a, b, c, x, r)
            
            r=rhs-r		! b array in Ax=b (known)
            r_old=r
            rho_old=1.d0
            alpha=1.d0	 ! to make beta = 0 and thus p = r for the first iteration
            omega=1.d0
            v=0.d0
            p=0.d0
            
            norm_b=max(sqrt(dot_product(rhs, rhs)), 1.d0)
            
            do k=1,max_iter
                x_old=x
                
                rho=dot_product(r_old, r)
                beta=(rho/rho_old)*(alpha/omega)
                p=r+beta*(p-omega*v)   
                
                call matvect_arrow(n, a, b, c, p, v)
                            
                
                alpha=rho/dot_product(r_old, v)
                s=r-alpha*v       
                        
                error=sqrt(dot_product(s,s))
                
                if (error/norm_b<tol) then
                    x=x_old+alpha*p
                    info=0
                    return
                end if
                
                call matvect_arrow(n, a, b, c, s, t)
                omega=dot_product(t, s)/dot_product(t, t)
                x=x_old+alpha*p+omega*s

                r=s-omega*t
                error=sqrt(dot_product(r,r))
                if (error/norm_b<tol) then
                    info=0
                    return
                end if
                
                rho_old=rho
            end do
            
            info=1 ! Flag for the missed convergence
            return
        end subroutine BiCGSTAB_HGCMFD
        
        
        
        subroutine matvect_arrow(n, a, b, c, x, y) 
            integer, intent(in) :: n
            real(dprec), intent(in) :: a(:), b(:), c(:), x(:)
            real(dprec), intent(out):: y(:)
            integer :: i
            
            y(1) = b(1)*x(1)
            do i = 2, n
                y(1)=y(1)+c(i)*x(i)
                y(i)=a(i)*x(1)+b(i)*x(i)
            end do
            
        end subroutine matvect_arrow
    
end module HGCMFD_BiCGSTAB