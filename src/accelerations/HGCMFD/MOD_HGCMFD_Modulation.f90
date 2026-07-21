module HGCMFD_Modulation
    
    use HGCMFD_Variables, only: dprec, Options_Data_2, LO_geom_2, LO_coeff_2
        implicit none
        
    contains
    
        
        subroutine Modulation(Opt_lo, Opt_gl, LO_Mesh, LO_Param, Phi_gl, Phi_avg_lo, Phi_lo) 
        
            type(Options_Data_2), intent(in) :: Opt_lo, Opt_gl
            type(LO_geom_2), intent(in) :: LO_Mesh(:)
            type(LO_coeff_2), intent(in) :: LO_Param(:)
            real(dprec), intent(in) :: Phi_gl(:), Phi_avg_lo(:)
            real(dprec), intent(inout) :: Phi_lo(:)
            
            integer :: l, ll, t, g, u, t_tot, t_red, g_tot
            
            do l=1, Opt_lo%n_tot
                u=LO_Mesh(l)%lo_gl_Homo ! Global spatial index relative to the total grid (w/ ghosts) 
                ll=LO_Mesh(l)%ID_Red ! Local spatial index relative to the reduced grid (w/o ghosts) 
                if (ll>0) then
                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot    ! Local energy index base relative to the total grid (w/ ghosts)
                        t_red=(t-1)*Opt_lo%n_red    ! Local energy index base relative to the reduced grid (w/o ghosts) 
                        
                        g=LO_Param(t_tot+l)%En_ID  ! Global energy index
                        g_tot=(g-1)*Opt_gl%n_tot    ! Global energy index base
                        
                        Phi_lo(t_red+ll)=Phi_lo(t_red+ll)*Phi_gl(g_tot+u)/Phi_avg_lo(g_tot+u)
                    end do
                end if
            end do
            
        end subroutine Modulation
    
end module