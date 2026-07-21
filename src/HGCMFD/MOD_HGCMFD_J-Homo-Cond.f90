module HGCMFD_J_Homo_Cond
    
    use HGCMFD_Variables, only: dprec, Options_Data_2, LO_geom_2, LO_coeff_2, GL_geom_2, GL_coeff_2
    
        implicit none
    
    contains
    
        
        subroutine Currents_lo_to_gl(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param) 
            type(Options_Data_2), intent(in) :: Opt_lo, Opt_gl
            type(LO_geom_2), intent(in) :: LO_Mesh(:)
            type(LO_coeff_2), intent(in) :: LO_Param(:)
            type(GL_geom_2), intent(in) :: GL_Mesh(:)
            type(GL_coeff_2), intent(inout) :: GL_Param(:)
            
            integer :: i, g, l, u, n, t, t_tot, g_tot, infc, face_ID, face_ID_GL ! In Line Face Calc.      
            real(dprec) :: Phi_S
            
            ! Initialization 
            do i=1, Opt_gl%n_tot*Opt_gl%n_g
                GL_Param(i)%J_net(:)=0.0_dprec
                GL_Param(i)%J_part(:,:)=0.0_dprec
            end do
            
            
            ! Integration 
            do l=1, Opt_lo%n_tot
                u=LO_Mesh(l)%lo_gl_Homo
                if (u>-1) then  ! ghost mesh filtering
                    
                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot
                        g=LO_Param(t_tot+l)%En_ID
                        g_tot=(g-1)*Opt_gl%n_tot
                        
                        do n=1, size(LO_Mesh(l)%Neigh_ID)     ! Iteration on the surfaces
                            
                            if (LO_Mesh(l)%J_flag(n)>0 .AND. u>1) then  !  .AND. LO_Mesh(l)%BC(n)==0
                                face_ID=LO_Mesh(l)%Sides_ID(n) ! ID Face along which the universe changes
                                face_ID_GL=findloc(GL_Mesh(u)%J_lo_gl_face, face_ID, dim=1) ! ID Face along which the universe changes in the global
                                
                                ! Net 
                                GL_Param(g_tot+u)%J_net(face_ID_GL)=GL_Param(g_tot+u)%J_net(face_ID_GL)+LO_Param(t_tot+l)%J_net(n)*LO_Mesh(l)%A_lo(n)
                                
                                ! Partial 
                                GL_Param(g_tot+u)%J_part(face_ID_GL,1)=GL_Param(g_tot+u)%J_part(face_ID_GL,1)+LO_Param(t_tot+l)%J_part(n,1)*LO_Mesh(l)%A_lo(n)
                                GL_Param(g_tot+u)%J_part(face_ID_GL,2)=GL_Param(g_tot+u)%J_part(face_ID_GL,2)+LO_Param(t_tot+l)%J_part(n,2)*LO_Mesh(l)%A_lo(n)
                            
                            elseif (LO_Mesh(l)%BC(n)>0 .AND. u==1) then 
                                face_ID=LO_Mesh(l)%Sides_ID(n) ! ID Face along which the universe changes
                                face_ID_GL=findloc(GL_Mesh(u)%J_lo_gl_face, face_ID, dim=1)
                                
                                ! Net
                                GL_Param(g_tot+u)%J_net(face_ID_GL)=GL_Param(g_tot+u)%J_net(face_ID_GL)+LO_Param(t_tot+l)%J_net(n)*LO_Mesh(l)%A_lo(n)
                                
                                ! Partial 
                                GL_Param(g_tot+u)%J_part(face_ID_GL,1)=GL_Param(g_tot+u)%J_part(face_ID_GL,1)+LO_Param(t_tot+l)%J_part(n,1)*LO_Mesh(l)%A_lo(n)
                                GL_Param(g_tot+u)%J_part(face_ID_GL,2)=GL_Param(g_tot+u)%J_part(face_ID_GL,2)+LO_Param(t_tot+l)%J_part(n,2)*LO_Mesh(l)%A_lo(n)
                                
                            end if
                        end do
                    end do
                end if
            end do
                    
            ! Homogenization & Condensation  
            do l=1, Opt_gl%n_tot
                do n=1, size(GL_Mesh(l)%A_gl)
                    do g=1, Opt_gl%n_g
                        g_tot=(g-1)*Opt_gl%n_tot
                        
                        ! Net
                        GL_Param(g_tot+l)%J_net(n)=GL_Param(g_tot+l)%J_net(n)/GL_Mesh(l)%A_gl(n)
                        
                        ! Partial
                        GL_Param(g_tot+l)%J_part(n,1)=GL_Param(g_tot+l)%J_part(n,1)/GL_Mesh(l)%A_gl(n)
                        GL_Param(g_tot+l)%J_part(n,2)=GL_Param(g_tot+l)%J_part(n,2)/GL_Mesh(l)%A_gl(n)
                    end do
                end do
            end do
            
    end subroutine Currents_lo_to_gl
    
end module HGCMFD_J_Homo_Cond