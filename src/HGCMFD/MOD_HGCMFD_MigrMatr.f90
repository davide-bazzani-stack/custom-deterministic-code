module HGCMFD_Migration_Matrix
    
    use HGCMFD_Variables, only: dprec, GL_geom_2, GL_coeff_2
    
        implicit none
    
    contains
    
        subroutine Migration_Matr(n_g, n_tot, XS_RM, GL_Mesh, GL_Param, MM)
            integer, intent(in) :: n_g, n_tot
            real(dprec), intent(in) :: XS_RM(:)
            
            type(GL_geom_2), intent(in) :: GL_Mesh(:)
            type(GL_coeff_2), intent(in) :: GL_Param(:)
            
            real(dprec), intent(inout) :: MM(:,:)
            
            integer :: t, l, t_tot
            real(dprec) :: area_sum
            
            
            MM=0.0_dprec
            
            do t=1, n_g
                t_tot=(t-1)*n_tot
                MM(t_tot+1,t_tot+1)=XS_RM(t_tot+1)    ! BG mesh removal XS
                
                do l=2, n_tot
                    area_sum=sum(GL_Mesh(l)%A_gl(:))
                    
                    MM(t_tot+l,t_tot+l)=     XS_RM(t_tot+l)+(GL_Param(t_tot+l)%D_til_el+GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(l)%V_gl      ! Main diagonal
                    MM(t_tot+l,t_tot+1)=                   -(GL_Param(t_tot+l)%D_til_el-GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(l)%V_gl      ! First column (BG coupling)
                    MM(t_tot+1,t_tot+l)=                   -(GL_Param(t_tot+l)%D_til_el+GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(1)%V_gl      ! First Row (BG eqn.)
                    
                    MM(t_tot+1,t_tot+1)=MM(t_tot+1,t_tot+1)+(GL_Param(t_tot+l)%D_til_el-GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(1)%V_gl      ! BG diagonal contribution
                end do
                
                ! Adding the boundaries
                area_sum=sum(GL_Mesh(1)%A_gl(:))
                MM(t_tot+1,t_tot+1)=MM(t_tot+1,t_tot+1)+GL_Param(t_tot+1)%D_hat_el(1)*area_sum/GL_Mesh(1)%V_gl
            end do
            
        end subroutine Migration_Matr
        
        
        
        subroutine Migration_Vect(n_g, n_tot, XS_RM, GL_Mesh, GL_Param, a_MM, b_MM, c_MM)
            integer, intent(in) :: n_g, n_tot
            real(dprec), intent(in) :: XS_RM(:)
            
            type(GL_geom_2), intent(in) :: GL_Mesh(:)
            type(GL_coeff_2), intent(in) :: GL_Param(:)
            !type(Accel_Vars_Vect), intent(inout) :: Serv_Vect(:)
            
            real(dprec), intent(out) :: a_MM(:), b_MM(:), c_MM(:)
            
            integer :: t, l, t_tot
            real(dprec) :: area_sum
            
            
            a_MM=0.0_dprec
            b_MM=0.0_dprec
            c_MM=0.0_dprec
            
            do t=1, n_g
                t_tot=(t-1)*n_tot
                b_MM(t_tot+1)=XS_RM(t_tot+1)    ! BG mesh removal XS
                
                do l=2, n_tot
                    area_sum=sum(GL_Mesh(l)%A_gl(:))
                    
                    b_MM(t_tot+l)=XS_RM(t_tot+l)+(GL_Param(t_tot+l)%D_til_el+GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(l)%V_gl      ! Main diagonal
                    a_MM(t_tot+l)=              -(GL_Param(t_tot+l)%D_til_el-GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(l)%V_gl      ! First column (BG coupling)
                    c_MM(t_tot+l)=              -(GL_Param(t_tot+l)%D_til_el+GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(1)%V_gl      ! First Row (BG eqn.)
                    
                    b_MM(t_tot+1)= b_MM(t_tot+1)+(GL_Param(t_tot+l)%D_til_el-GL_Param(t_tot+l)%D_hat_el(1))*area_sum/GL_Mesh(1)%V_gl      ! BG diagonal contribution
                end do
                
                ! Adding the boundaries
                area_sum=sum(GL_Mesh(1)%A_gl(:))
                b_MM(t_tot+1)=b_MM(t_tot+1)+GL_Param(t_tot+1)%D_hat_el(1)*area_sum/GL_Mesh(1)%V_gl
            end do
                
        end subroutine Migration_Vect
    
    
end module HGCMFD_Migration_Matrix