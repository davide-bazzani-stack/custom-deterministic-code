module HGCMFD_D_Hat
    
    use HGCMFD_Variables, only: dprec, epsi, GL_geom_2, GL_coeff_2
        implicit none
    
    contains
    
    subroutine D_Hat(n_g, n_tot, Phi_Ref, GL_Mesh, GL_Param) 
        integer, intent(in) :: n_g, n_tot
        real(dprec), intent(in) :: Phi_Ref(:)
        
        type(GL_geom_2), intent(in) :: GL_Mesh(:)
        type(GL_coeff_2), intent(inout) :: GL_Param(:)
                
        integer :: t, l, t_tot, no_faces
        real(dprec) :: area_sum, num, den
        
        ! For the BC
        no_faces=size(GL_Mesh(1)%A_gl(:))
        area_sum=sum(GL_Mesh(1)%A_gl(:))
        
        do t=1,n_g
            t_tot=(t-1)*n_tot
            do l=2, n_tot   ! Rolling on the meshes
                num=GL_Param(t_tot+l)%J_Net_el+GL_Param(t_tot+l)%D_til_el*(Phi_Ref(t_tot+1)-Phi_Ref(t_tot+l))
                den=Phi_Ref(t_tot+1)+Phi_Ref(t_tot+l)
                if (abs(den) > epsi) then
                    GL_Param(t_tot+l)%D_hat_el(1)=num/den
                else
                    GL_Param(t_tot+l)%D_hat_el(1)=0.0_dprec
                end if
                
            end do
            
            ! BC treatment ---> sum(J/Phi*A/A_tot)
            GL_Param(t_tot+1)%D_hat_el(1)=dot_product(GL_Param(t_tot+1)%J_Net(1:no_faces), GL_Mesh(1)%A_gl(1:no_faces))/(area_sum*Phi_Ref(t_tot+1))
        end do
        
        
    end subroutine D_Hat
    
    
end module HGCMFD_D_Hat