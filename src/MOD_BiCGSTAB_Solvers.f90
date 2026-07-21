module BiCGSTAB_Algorithm 
	    implicit none
	integer, parameter :: dp_var = kind(1.d0)
    
    contains
    
    ! Made for traditional pentadiagonal matrices. If a different matrix is sought, modify the matvect_penta function
    subroutine BiCGSTAB_Classic(n_x, n, a, b, c, d, e, rhs, x, tol, max_iter, info) 
        integer, intent(in) :: n_x, n, max_iter
        real(dp_var), intent(in) :: a(:), b(:), c(:), d(:), e(:), rhs(:), tol
        real(dp_var), intent(inout) :: x(:)
        integer, intent(out) :: info
        real(dp_var), allocatable :: r(:), r_old(:), v(:), p(:), s(:), t(:), x_old(:)
        real(dp_var) :: alpha, beta, omega, rho, rho_old, error, t1, t2, t3, t4
        integer :: k
        
        allocate(r(n))
        allocate(r_old(n))
        allocate(v(n))
        allocate(p(n))
        allocate(s(n))
        allocate(t(n))
        
        allocate(x_old(n))
        

        call matvect_penta(n_x, n, a, b, c, d, e, x, r)

        r=rhs-r		! b array in Ax=b (known)
        r_old=r
        rho_old=1.d0
        alpha=1.d0	 ! to make beta = 0 and thus p = r for the first iteration
        omega=1.d0
        v=0.d0
        p=0.d0
        
        IterativeSolution: do k=1,max_iter
            
            x_old=x
            
            rho=dot_product(r_old, r)
            beta=(rho/rho_old)*(alpha/omega)
            p=r+beta*(p-omega*v)   
            
            call matvect_penta(n_x, n, a, b, c, d, e, p, v)
                        
            
            alpha=rho/dot_product(r_old, v)
            s=r-alpha*v       
                    
            error=sqrt(dot_product(s,s))
            
            !if (error<tol) then
            !    x=x+alpha*p
            !    info=0
            !    return
            !end if
            
            call matvect_penta(n_x, n, a, b, c, d, e, s, t)
            omega=dot_product(t, s)/dot_product(t, t)
            x=x+alpha*p+omega*s

            r=s-omega*t
            error=sqrt(dot_product(x-x_old,x-x_old)/dot_product(x,x))
            if (error<tol) then
                info=0
                return
            end if
            
            rho_old=rho
        end do IterativeSolution
        
        info=1 ! Flag for the missed convergence
        return
    end subroutine BiCGSTAB_Classic
        
    subroutine matvect_penta(n_x, n, a, b, c, d, e, x, y) 
        integer, intent(in) :: n_x, n
        real(dp_var), intent(in) :: a(:), b(:), c(:), d(:), e(:), x(:)
        real(dp_var), intent(out) :: y(:)
        integer :: i
        
        y=0.d0
        
        do i=2,n-1
            y(i)=y(i)+b(i)*x(i)+a(i)*x(i-1)+c(i)*x(i+1)
            if (i>n_x) then
                y(i)=y(i)+d(i)*x(i-n_x)
            end if
            if (i<=n-n_x) then
                y(i)=y(i)+e(i)*x(i+n_x)
            end if
        end do
        
        y(1)=b(1)*x(1)+c(1)*x(2)+e(1)*x(1+n_x)
        y(n)=b(n)*x(n)+a(n)*x(n-1)+d(n)*x(n-n_x)
    end subroutine matvect_penta
    
    
    
    subroutine BiCGSTAB_EqTri(n_x, n_y, a, b, c, d, e, rhs, x, tol, max_iter, info) 
        integer, intent(in) :: n_x, n_y, max_iter
        real(dp_var), intent(in) :: a(:), b(:), c(:), d(:), e(:), rhs(:), tol
        real(dp_var), intent(inout) :: x(:)
        integer, intent(out) :: info
        real(dp_var), allocatable :: r(:), r_old(:), v(:), p(:), s(:), t(:), x_old(:)
        real(dp_var) :: alpha, beta, omega, rho, rho_old, error, t1, t2, t3, t4
        integer :: k, n
        
        n=n_x*n_y
        
        allocate(r(n))
        allocate(r_old(n))
        allocate(v(n))
        allocate(p(n))
        allocate(s(n))
        allocate(t(n))
        
        allocate(x_old(n))
        
        
        call matvect_penta_EqTri(n_x, n, a, b, c, d, e, x, r)
        
        r=rhs-r		! b array in Ax=b (known)
        r_old=r
        rho_old=1.d0
        alpha=1.d0	 ! to make beta = 0 and thus p = r for the first iteration
        omega=1.d0
        v=0.d0
        p=0.d0
        
        IterativeSolution: do k=1,max_iter
            
            x_old=x
            
            rho=dot_product(r_old, r)
            beta=(rho/rho_old)*(alpha/omega)
            p=r+beta*(p-omega*v)   
            
            call matvect_penta_EqTri(n_x, n, a, b, c, d, e, p, v)
                        
            
            alpha=rho/dot_product(r_old, v)
            s=r-alpha*v       
                    
            error=sqrt(dot_product(s,s))
            
            !if (error<tol) then
            !    x=x+alpha*p
            !    info=0
            !    return
            !end if
            
            call matvect_penta_EqTri(n_x, n, a, b, c, d, e, s, t)
            omega=dot_product(t, s)/dot_product(t, t)
            x=x+alpha*p+omega*s

            r=s-omega*t
            error=sqrt(dot_product(x-x_old,x-x_old)/dot_product(x,x))
            if (error<tol) then
                info=0
                return
            end if
            
            rho_old=rho
        end do IterativeSolution
        
        info=1 ! Flag for the missed convergence
        return
    end subroutine BiCGSTAB_EqTri
    
    subroutine matvect_penta_EqTri(n_x, n, a, b, c, d, e, x, y) 
        integer, intent(in) :: n_x, n
        real(dp_var), intent(in) :: a(n), b(n), c(n), d(n), e(n), x(n)
        real(dp_var), intent(out) :: y(n)
        integer :: i
        
        y=b*x
        
        do i=1, n
            if (i>n_x) then
                y(i)=y(i)+a(i)*x(i-n_x)
            end if
            if (i>n_x+1) then
                y(i)=y(i)+d(i)*x(i-n_x-1)
            end if
            if (i<=n-n_x) then
                y(i)=y(i)+c(i)*x(i+n_x)
            end if
            if (i<=n-n_x-1) then
                y(i)=y(i)+e(i)*x(i+n_x+1)
            end if
        end do
        
    end subroutine matvect_penta_EqTri
    
    
    
    subroutine BiCGSTAB_TriIso(n_x, n, a, b, c, d, e, rhs, x, tol, max_iter, info) 
        integer, intent(in) :: n_x, n, max_iter
        real(dp_var), intent(in) :: a(:), b(:), c(:), d(:), e(:), rhs(:), tol
        real(dp_var), intent(inout) :: x(:)
        integer, intent(out) :: info
        real(dp_var), allocatable :: r(:), r_old(:), v(:), p(:), s(:), t(:), x_old(:)
        real(dp_var) :: alpha, beta, omega, rho, rho_old, error, t1, t2, t3, t4
        integer :: k
                
        allocate(r(n))
        allocate(r_old(n))
        allocate(v(n))
        allocate(p(n))
        allocate(s(n))
        allocate(t(n))
        
        allocate(x_old(n))
        
        
        call matvect_penta_TriIso(n_x, n, a, b, c, d, e, x, r)
        
        r=rhs-r		! b array in Ax=b (known)
        r_old=r
        rho_old=1.d0
        alpha=1.d0	 ! to make beta = 0 and thus p = r for the first iteration
        omega=1.d0
        v=0.d0
        p=0.d0
        
        IterativeSolution: do k=1,max_iter
            
            x_old=x
            
            rho=dot_product(r_old, r)
            beta=(rho/rho_old)*(alpha/omega)
            p=r+beta*(p-omega*v)   
            
            call matvect_penta_TriIso(n_x, n, a, b, c, d, e, p, v)
                        
            
            alpha=rho/dot_product(r_old, v)
            s=r-alpha*v       
                    
            error=sqrt(dot_product(s,s))
            
            !if (error<tol) then
            !    x=x+alpha*p
            !    info=0
            !    return
            !end if
            
            call matvect_penta_TriIso(n_x, n, a, b, c, d, e, s, t)
            omega=dot_product(t, s)/dot_product(t, t)
            x=x+alpha*p+omega*s

            r=s-omega*t
            error=sqrt(dot_product(x-x_old,x-x_old)/dot_product(x,x))
            if (error<tol) then
                info=0
                return
            end if
            
            rho_old=rho
        end do IterativeSolution
        
        info=1 ! Flag for the missed convergence
        return
    end subroutine BiCGSTAB_TriIso
    
    subroutine matvect_penta_TriIso(n_x, n, a, b, c, d, e, x, y) 
        integer, intent(in) :: n_x, n
        real(dp_var), intent(in) :: a(n), b(n), c(n), d(n), e(n), x(n)
        real(dp_var), intent(out) :: y(n)
        integer :: i
        
        y=b*x
        
        do i=1, n
            if (i>1) then
                y(i)=y(i)+a(i)*x(i-1)
            end if
            if (i>n_x) then
                y(i)=y(i)+d(i)*x(i-n_x)
            end if
            if (i<n) then
                y(i)=y(i)+c(i)*x(i+1)
            end if
            if (i<=n-n_x) then
                y(i)=y(i)+e(i)*x(i+n_x)
            end if
        end do
        
    end subroutine matvect_penta_TriIso
    
! Insert here more different versions of matvect_penta
    
end module BiCGSTAB_Algorithm 
