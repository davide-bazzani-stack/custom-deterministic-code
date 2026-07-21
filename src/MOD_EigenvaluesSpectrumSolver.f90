module Eigen_spectrum_solver
		use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
		implicit none
		
	integer, parameter :: doupr=kind(1.d0)
	
    contains
		
		
		subroutine Eigen_Routine(max_mesh_NewCMFD, n_g, NewCMFD_MigrMat, Chi_Homo, XS_nF_Homo, XS_SC_Homo, Phi, eigens_R, eigens_I) ! Adapter for the full eigenvalue problem 
			integer :: t, j, l, g
            integer, intent(in) :: max_mesh_NewCMFD, n_g
			real(doupr), intent(in) :: NewCMFD_MigrMat(:,:,:), Chi_Homo(:,:), XS_nF_Homo(:,:), XS_SC_Homo(:,:,:), Phi(:,:)
			real(doupr), allocatable, intent(out) :: eigens_R(:), eigens_I(:)
            
			real(doupr), allocatable :: Migrat_M(:,:), Out_Sc_M(:,:), Fiss_M(:,:), v1(:), Full_Phi(:,:), residuals(:)
            
            allocate( Migrat_M(max_mesh_NewCMFD*n_g, max_mesh_NewCMFD*n_g))
            allocate( Out_Sc_M(max_mesh_NewCMFD*n_g, max_mesh_NewCMFD*n_g))
            allocate(   Fiss_M(max_mesh_NewCMFD*n_g, max_mesh_NewCMFD*n_g))
            allocate(       v1(max_mesh_NewCMFD*n_g))
            allocate(residuals(max_mesh_NewCMFD*n_g))
            allocate( Full_Phi(max_mesh_NewCMFD*n_g, max_mesh_NewCMFD*n_g))
            
            allocate(eigens_R(max_mesh_NewCMFD*n_g), eigens_I(max_mesh_NewCMFD*n_g))
            
            Migrat_M=0.d0
			Out_Sc_M=0.d0
            Fiss_M=0.d0
            
            
            do t=1,n_g
                v1((t-1)*max_mesh_NewCMFD+1:t*max_mesh_NewCMFD)=Phi(:,t)
                do j=1, max_mesh_NewCMFD ! Column index
                    l=(t-1)*max_mesh_NewCMFD
                    
					Migrat_M(l+j, l+1:l+max_mesh_NewCMFD)=NewCMFD_MigrMat(j,:,t)	!NewCMFD_MigrMat((j-1)*max_mesh_NewCMFD+1:j*max_mesh_NewCMFD,t)
					Fiss_M(l+j,l+j)=Chi_Homo(j,t)*XS_nF_Homo(j,t)
                    
                    do g=1,n_g
                        if (t .NE. g) then
                            Fiss_M(l+j, (g-1)*max_mesh_NewCMFD+j)=Chi_Homo(j,t)*XS_nF_Homo(j,g)
                            Out_Sc_M(l+j, (g-1)*max_mesh_NewCMFD+j)=XS_SC_Homo(j, g, t) ! In-scattering
                        end if
                    end do
                end do
            end do
            
            Migrat_M=Migrat_M-Out_Sc_M
            
			call Arnoldi_alg(max_mesh_NewCMFD*n_g, v1, Fiss_M, Migrat_M, eigens_R, eigens_I, Full_Phi, residuals)
            
            ! Writing the Matrices and the results 
            write(4,'(A)') 'Migration Matrix M   /   Fission_Matrix   /   In-scattering'
            do j=1,max_mesh_NewCMFD*n_g
                write(4,'(100000ES21.13)') Migrat_M(j,:)+Out_Sc_M(j,:), Fiss_M(j,:), Out_Sc_M(j,:)
            end do
            write(4,*)
            write(4,*)
            write(4,*)
            write(4,'(A)') 'Combined matrix A'
            do j=1,max_mesh_NewCMFD*n_g
                write(4,'(100000ES21.13)') Migrat_M(j,:)
            end do
            write(4,*)
            write(4,*)
            write(4,*)
            write(4,'(A)') 'Eigenvalues and Eigenvectors'
            write(4,'(A)') 'Eigenvalues'
            do j=1,SIZE(eigens_R)
                write(4,'(2ES21.13,A)') eigens_R(j), eigens_I(j), 'i'
            end do
            write(4,*)
            write(4,*)
            write(4,*)
            write(4,'(A)') 'Eigenvectors'
            do j=1,max_mesh_NewCMFD*n_g
                write(4,'(100000ES21.13)') Full_Phi(j,:) ! Eigenvector elements on the column
            end do
            write(4,*)
            write(4,*)
            write(4,*)
            write(4,'(A)') 'Residuals'
            write(4,'(100000ES21.13)') residuals
            write(4,*)
            write(4,*)
            write(4,*)
            write(5,'(A)') 'Eigenvalues'
            do j=1,SIZE(eigens_R)
                write(5,'(2ES21.13,A)') eigens_R(j), eigens_I(j), 'i'
            end do
            write(5,*)
            write(5,*)
            write(5,*)
        end subroutine Eigen_Routine
	
		subroutine Arnoldi_alg(m, Phi, F, A, eigens_R, eigens_I, eigvect, res_tot) ! Arnoldi algorithm for non-symmetric large matrices. Return ALL the eigenvalues (not sorted) 
			integer, intent(in) :: m
			real(doupr), intent(in) :: Phi(:), F(:,:), A(:,:) ! A= Migration - (out-)Scattering
			real(doupr), allocatable, intent(out) :: eigens_R(:), eigens_I(:), eigvect(:,:), res_tot(:)
			real(doupr), allocatable :: H_m(:,:), H_upp_tr(:,:), V(:,:), H(:,:), z(:), w(:), Y_eigvect(:,:)
			
			integer :: i, j, n, info=0
			real(doupr) :: epsi=1.d-12
			
			n=size(Phi)
			
			allocate(V(n, m+1))
			allocate(H(m+1,m))
			allocate(w(m))
			
			V=0.d0
			H=0.d0
			w=0.d0
			
			V(:,1)=Phi/norm_l2(Phi)
			
			
			Arnoldi_loop: do j=1,m 
				z = matmul(F, V(:, j))     ! Matrix-free approach. 
                                           ! No explicit inversion of A = M - S
				
				call BiCGSTAB_Mat(n, A, z, w, 1.d-10, 1000000, info)
				if (info/=0) exit Arnoldi_loop
				
				GramShmidt_ortogon: do i=1, j
					H(i,j) = dot_product(V(:,i),w)
					w = w - H(i,j)*V(:,i)
				end do GramShmidt_ortogon
				
				H(j+1,j)=norm_l2(w) ! l_2 module of w
				
				if (H(j+1,j)<1.d-12) exit Arnoldi_loop
				
				V(:,j+1)=w/H(j+1,j)
				
            end do Arnoldi_loop
            
            j=minval([j,m])
            
			if (info==0) then
				! Eigenvalues computation 
				allocate(H_m(j,j))
				allocate(H_upp_tr(j,j))
				allocate(eigens_R(j), eigens_I(j)) ! must coincide with the 'm' INSIDE QR_main
                allocate(Y_eigvect(n,j))
                allocate(eigvect(n,j))
                allocate(res_tot(j))
				H_m=H(:j,:j)
                eigens_R=0.d0
                eigens_I=0.d0
                Y_eigvect=0.d0
				call QR_main_Complex(j, H_m, epsi, eigens_R, eigens_I, Y_eigvect)
                
                do i=1, j
					eigvect(:,i)=matmul(V(:,1:j),Y_eigvect(:,i))
					res_tot(i)=abs(H(m+1,m)*eigvect(m,i))
                end do
                
                
			else
				print*, 'The BiCGSTAB in the Arnoldi algorithm did not converge'
            end if
			
            
		end subroutine Arnoldi_alg
		
		function norm_l2(x) result(norm_x) 
			real(doupr), intent(in) :: x(:)
			real(doupr) :: norm_x
			
			norm_x=sqrt(dot_product(x,x))
		end function norm_l2
				
		subroutine BiCGSTAB_Mat(n, A, b, x, tol, max_iter, info) 
			integer, intent(in) :: n, max_iter
			real(doupr), intent(in) :: A(:,:), b(:), tol
			real(doupr), intent(inout) :: x(:)
			integer, intent(out) :: info
			real(doupr), allocatable :: r(:), r_old(:), v(:), p(:), s(:), t(:), x_old(:)
			real(doupr) :: alpha, beta, omega, rho, rho_old, error
			integer :: k
			
            allocate(r(n), r_old(n), v(n), p(n), s(n), t(n), x_old(n))
            
			r=matmul(A,x)
			
			r=b-r		! b array in Ax=b (known)
			r_old=r
			rho_old=1.d0
			alpha=1.d0	 ! to make beta = 0 and thus p = r for the first iteration
			omega=1.d0
			v=0.d0
			p=0.d0
			
			Iter_Loop: do k=1, max_iter 
				x_old=x
				
				rho=dot_product(r_old, r)
				beta=(rho/rho_old)*(alpha/omega)
				p=r+beta*(p-omega*v)
				
				v=matmul(A, p)
	      
				alpha=rho/dot_product(r_old, v)
				s=r-alpha*v       
        
				error=sqrt(dot_product(s,s))
				
				if (error<tol) then
					x=x+alpha*p
					info=0
					return
				end if
				
				t=matmul(A, s)
        
				omega=dot_product(t, s)/dot_product(t, t)
				x=x+alpha*p+omega*s

				r=s-omega*t
				error=sqrt(dot_product(x-x_old,x-x_old)/dot_product(x,x))
				if (error<tol) then
					info=0
					return
                end if
                
                if (ieee_is_nan(error)) then
					info=1
                    print*, 'Missed convergence - error = NaN'
                end if
	      
				rho_old=rho
            end do Iter_Loop
        end subroutine BiCGSTAB_Mat

        subroutine QR_main_Complex(m, H_real, epsi, eigens_R, eigens_I, Y_eigvect) ! QR Algorithm w/ shift and complex numbers handling 
			integer, intent(in) :: m
			real(doupr), intent(in) :: H_real(m,m), epsi
			real(doupr), intent(out) :: eigens_R(m), eigens_I(m), Y_eigvect(m,m)
			
			integer :: i, k, iter_max
			real(doupr) :: temp, mi
			
            complex(doupr) :: H_upp_tr(m,m), Q_acc(m,m)
            complex(doupr), allocatable :: eigens(:), H(:,:), H_temp(:,:), Q_temp(:,:), R_temp(:,:), eye(:,:)
            
			iter_max=1e6
			eigens=0.d0
            
            allocate(eigens(m), H(m,m), H_temp(m,m), Q_temp(m,m), R_temp(m,m), eye(m,m))
            
            H=cmplx(H_real, 0.d0, kind=kind(doupr)) ! To verify my need kind=kind(doupr), cmplx(H_real, kind=doupr)
			H_upp_tr=0.d0
			
			! Identity matrix building
			eye=0.d0
			do i=1,m
				eye(i,i)=(1.d0, 0.d0)
            end do	
			Q_acc=eye
            
            
			! QR loop
			H_temp=H
			QR_Iter_loop: do k=1, iter_max
				
				call QR_shift_Complex(H_temp(m-1, m-1),  H_temp(m, m-1), H_temp(m, m), mi) ! mi computation
				H_temp=H_temp-mi*eye ! Shifting of mi via identity matrix
				
				call QR_decomp_Complex(m, H_temp, Q_temp, R_temp)
				H_temp=matmul(R_temp, Q_temp)
				
				H_temp=H_temp+mi*eye ! Un-shifting of mi via identity matrix
				
                ! For eigenvectors
                Q_acc=matmul(Q_acc, Q_temp)
                
				! Check the upper-triangular layout
				temp=0.d0
				do i=2, size(H_temp, 1)
					temp=temp+sum(abs(real(H_temp(i,1:i-1))))+sum(abs(aimag(H_temp(i,1:i-1)))) ! +sum(abs(H_temp(i,1:i-1))) for the precise computation of the modulus (more expensive)
				end do
				if (temp <= epsi) exit
			end do QR_Iter_loop
			
			if (k == iter_max) print*, 'Warning: Shifted QR algorithm did not converge'
			
			! Eigenvalues extraction
			H_upp_tr=H_temp
			do i=1,size(H_upp_tr,1)
				eigens(i)=H_upp_tr(i,i)
            end do
            
            
            eigens_R=real(eigens)
            eigens_I=aimag(eigens)
            Y_eigvect=real(Q_acc)
            
        end subroutine QR_main_Complex
        
        subroutine QR_decomp_Complex(m, H_mat, Q_mat, R_mat) 
			integer, intent(in) :: m
			complex(doupr), intent(in) :: H_mat(:,:)
			complex(doupr), intent(out) :: Q_mat(:,:), R_mat(:,:)
			
			integer :: i, j
			complex(doupr) :: a, b, r, c, s, temp1,temp2
			
			complex(doupr) :: G_mat(2,2), R_temp_mat(m,m)
			
			
			R_mat=H_mat
			Q_mat=0.d0
			do i=1,size(Q_mat,1)
				Q_mat(i,i)=cmplx(1.d0,0.d0, kind=doupr) ! Identity matrix structure
			end do
			
			do i=1,m-1
				a=R_mat(i,i)
				b=R_mat(i+1,i)
				r=sqrt(conjg(a)*a+conjg(b)*b)
				if (abs(real(r))+abs(aimag(r))<1.d-14) cycle ! Skipping the cycle if r = 0.0 for robustness
				c=a/r
				s=b/r
				
				G_mat(1,1)= c
				G_mat(1,2)= s         ! G=| c   s |
				G_mat(2,1)=-conjg(s)  !   |-s*  c*|
				G_mat(2,2)= conjg(c)
				
				! Over R_mat				
				do j=1,size(R_mat,2)
					temp1=G_mat(1,1)*R_mat(i,j)+G_mat(1,2)*R_mat(i+1,j)
					temp2=G_mat(2,1)*R_mat(i,j)+G_mat(2,2)*R_mat(i+1,j)
					
					R_mat(i  ,j)=temp1
					R_mat(i+1,j)=temp2
				end do
				
				! Over Q_mat, with G transposed(!)
				do j=1,size(Q_mat,1)
					temp1=G_mat(1,1)*Q_mat(j,i)+G_mat(1,2)*Q_mat(j,i+1)
					temp2=G_mat(2,1)*Q_mat(j,i)+G_mat(2,2)*Q_mat(j,i+1)
					
					Q_mat(j  ,i)=temp1
					Q_mat(j,i+1)=temp2
				end do
			end do
			
			! Verify that
			!   Transpose(Q)*H==R
			!   H==Q*R
			
			! Debugging 
			! write(*,('')) H(i,j)-dot_product(Q(i,:), R(:,j))
				
        end subroutine QR_decomp_Complex
        
        subroutine QR_shift_Complex(a, b, c, mi) ! Works also for non-symmetric matrices 
			complex(doupr), intent(in) :: a, b, c ! sub matrix elements from the Hessemberg. 
			real(doupr), intent(out) :: mi
			
            integer :: i(1)
			real(doupr) :: roots(2), temp(2)
            
			
			! sub_M = |1   2| = |H_(m-1, m-1)   H_(m-1,   m)|
			!         |3   4|   |H_(  m, m-1)   H_(  m,   m)|  
			!
			! For robustness, take always = 3 = H_(  m, m-1)
			
			! When explicitly passing the submatrix sum_M(2,2)
			!a=sub_M(1,1)
			!b=sub_M(2,1)
			!c=sub_M(2,2)
			
			roots(1)=(a+c)/2.d0+sqrt((a-c)**2/4.d0+b**2)
			roots(2)=(a+c)/2.d0-sqrt((a-c)**2/4.d0+b**2)
			
            i=minloc(abs(roots-c))
			mi=roots(i(1))
			
		end subroutine QR_shift_Complex
        
end module Eigen_spectrum_solver
