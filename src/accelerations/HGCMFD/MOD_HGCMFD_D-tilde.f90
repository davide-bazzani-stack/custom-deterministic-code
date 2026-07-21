module HGCMFD_D_tilde
    
    use HGCMFD_Variables, only: dprec, Options_Data_2, GL_geom_2, GL_coeff_2
        implicit none
    
    contains
    
    subroutine D_tilde(Opt_gl, GL_Mesh, D_gl, GL_Param) 
        real(dprec), intent(in) :: D_gl(:)
        
        type(Options_Data_2), intent(in) :: Opt_gl
        type(GL_geom_2), intent(in) :: GL_Mesh(:)
        type(GL_coeff_2), intent(inout) :: GL_Param(:)
        
        integer :: l, t, n, t_tot
        
        do l=1, size(GL_Param)
            GL_Param(l)%D_Til=0.0_dprec
            GL_Param(l)%D_Til_El=0.0_dprec
        end do
        
        
        do l=2, Opt_gl%n_tot ! Skipping the BG mesh, iteration "on the universe"
            
            do t=1, Opt_gl%n_g
                t_tot=(t-1)*Opt_gl%n_tot
                
                do n=1, size(GL_Param(t_tot+l)%D_Til,1)  ! Iteration on the number of sides
                    GL_Param(t_tot+l)%D_til(n,1)=D_gl(t_tot+l)/GL_Mesh(l)%dl_gl(n)*D_gl(t_tot+1)/GL_Mesh(1)%dl_gl(1)/(D_gl(t_tot+l)/GL_Mesh(l)%dl_gl(n)+D_gl(t_tot+1)/GL_Mesh(1)%dl_gl(1))
                    
                    GL_Param(t_tot+l)%D_til_el=GL_Param(t_tot+l)%D_til_el+GL_Param(t_tot+l)%D_til(n,1)*GL_Mesh(l)%A_gl(n)
                end do
                
                GL_Param(t_tot+l)%D_til_el=GL_Param(t_tot+l)%D_til_el/sum(GL_Mesh(l)%A_gl(:))
            end do
        end do
        
        
    end subroutine
    
    
end module HGCMFD_D_tilde