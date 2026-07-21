module Element_Current
    
    use HGCMFD_Variables, only: dprec, GL_geom_2, GL_coeff_2
    
        implicit none
    
    contains
    
        subroutine J_gl_to_El(n_tot, n_g, GL_Mesh, GL_Param) 
            integer, intent(in) :: n_g, n_tot
            
            type(GL_geom_2), intent(in) :: GL_Mesh(:)
            type(GL_coeff_2), intent(inout) :: GL_Param(:)
            
            integer :: l, t, t_tot
            integer :: no_faces
            real(dprec) :: area_sum        
            
            ! Non-BG Meshes
            do l=2, n_tot ! Skipping the BG mesh, iteration "on the universe"
                
                no_faces=size(GL_Mesh(l)%A_gl(:))
                area_sum=sum(GL_Mesh(l)%A_gl(:))
                
                do t=1, n_g
                    t_tot=(t-1)*n_tot
                    
                    GL_Param(t_tot+l)%J_Net_El=dot_product(GL_Param(t_tot+l)%J_Net(1:no_faces), GL_Mesh(l)%A_gl(1:no_faces))/area_sum
                end do
            end do
            
        end subroutine J_gl_to_El
        
end module Element_Current