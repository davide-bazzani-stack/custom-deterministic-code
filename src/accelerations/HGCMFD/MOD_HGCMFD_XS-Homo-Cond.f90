module HGCMFD_XS_Homo_Cond
    
    use HGCMFD_Variables, only: dprec, Options_Data_2, XS_Data_2, LO_geom_2, LO_coeff_2, GL_geom_2
    
        implicit none
    
    contains
    
        subroutine XS_HomoCond(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, XS_lo, XS_gl, Phi_lo, Phi_gl) 
            type(Options_Data_2), intent(in) :: Opt_lo, Opt_gl
            type(LO_geom_2), allocatable, intent(in) :: LO_Mesh(:)
            type(LO_coeff_2), allocatable, intent(in) :: LO_Param(:)
            type(GL_geom_2), allocatable, intent(in) :: GL_Mesh(:)
            real(dprec), allocatable, intent(in) :: Phi_lo(:)
            real(dprec), allocatable, intent(inout) :: Phi_gl(:)
            type(XS_Data_2), intent(in) :: XS_lo
            type(XS_Data_2), intent(inout) :: XS_gl
            
            integer :: u, n, l, l_red, t, t_red, t_tot, tt, tt_red, tt_tot, g, g_red, g_tot, gg, gg_red, gg_tot, i, j
            real(dprec) :: temp
            real(dprec), allocatable :: Phi_int(:)
            type(XS_Data_2) :: XS_int
            
            ! Allocation and initialization 
            allocate(  XS_int%Tot(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Tra(Opt_gl%n_tot*Opt_gl%n_g))
            allocate( XS_int%Absr(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(XS_int%nuFis(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Fis(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Chi(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(XS_int%Scatt(Opt_gl%n_tot*Opt_gl%n_g, Opt_gl%n_g))
            allocate( XS_int%kFis(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Rem(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Dif(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(     Phi_int(Opt_gl%n_tot*Opt_gl%n_g))
            
            XS_int%Tot=  0.0_dprec
            XS_int%Tra=  0.0_dprec
            XS_int%Absr= 0.0_dprec
            XS_int%nuFis=0.0_dprec
            XS_int%Fis=  0.0_dprec
            XS_int%Chi=  0.0_dprec
            XS_int%Scatt=0.0_dprec
            XS_int%kFis= 0.0_dprec
            XS_int%Rem=  0.0_dprec
            XS_int%Dif=  0.0_dprec
            Phi_int=     0.0_dprec
            
            ! Element-wise integral 
            do l=1, Opt_lo%n_tot
                n=LO_Mesh(l)%lo_gl_Homo ! Universe
                l_red=LO_Mesh(l)%ID_Red
                if (l_red>0) then ! Ghost layer filtering
                    
                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot
                        t_red=(t-1)*Opt_lo%n_red
                        g=LO_Param(t_tot+l)%En_ID
                        g_tot=(g-1)*Opt_gl%n_tot
                        
                          XS_int%Tot(g_tot+n)=  XS_int%Tot(g_tot+n)+  XS_lo%Tot(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                          XS_int%Tra(g_tot+n)=  XS_int%Tra(g_tot+n)+  XS_lo%Tra(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                         XS_int%Absr(g_tot+n)= XS_int%Absr(g_tot+n)+ XS_lo%Absr(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                        XS_int%nuFis(g_tot+n)=XS_int%nuFis(g_tot+n)+XS_lo%nuFis(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                          XS_int%Fis(g_tot+n)=  XS_int%Fis(g_tot+n)+  XS_lo%Fis(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                         XS_int%kFis(g_tot+n)= XS_int%kFis(g_tot+n)+ XS_lo%kFis(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                             Phi_int(g_tot+n)=     Phi_int(g_tot+n)+                     Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                          
                        do tt=1, Opt_lo%n_g
                            tt_tot=(tt-1)*Opt_lo%n_tot
                            tt_red=(tt-1)*Opt_lo%n_red
                            gg=LO_Param(tt_tot+l)%En_ID
                            
                            XS_int%Scatt(g_tot+n,gg)=XS_int%Scatt(g_tot+n,gg)+XS_lo%Scatt(t_tot+l,tt)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                            XS_int%Chi(g_tot+n)=XS_int%Chi(g_tot+n)+XS_lo%Chi(t_tot+l)*XS_lo%nuFis(tt_tot+l)*Phi_lo(tt_red+l_red)*LO_Mesh(l)%V_lo
                        end do
                    end do
                end if
            end do
            
            ! Homogenization 
            do l=1, Opt_gl%n_tot
                do g=1, Opt_gl%n_g
                    g_tot=(g-1)*Opt_gl%n_tot
                    XS_gl%Tot(g_tot+l)=  XS_int%Tot(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%Tra(g_tot+l)=  XS_int%Tra(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%Absr(g_tot+l)= XS_int%Absr(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%nuFis(g_tot+l)=XS_int%nuFis(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%Fis(g_tot+l)=  XS_int%Fis(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%kFis(g_tot+l)= XS_int%kFis(g_tot+l)/Phi_int(g_tot+l)
                    
                    XS_gl%Scatt(g_tot+l,:)=XS_int%Scatt(g_tot+l,:)/Phi_int(g_tot+l)
                    
                    Phi_gl(g_tot+l)=Phi_int(g_tot+l)/GL_Mesh(l)%V_gl
                end do
            end do
            
            ! Chi 
            do l=1, Opt_gl%n_tot
                temp=0.0_dprec
                do g=1, Opt_gl%n_g
                    g_tot=(g-1)*Opt_gl%n_tot
                    temp=temp+XS_gl%nuFis(g_tot+l)*Phi_gl(g_tot+l)*GL_Mesh(l)%V_gl
                end do
                if (abs(temp)>1.d-15) then
                    do g=1, Opt_gl%n_g
                        g_tot=(g-1)*Opt_gl%n_tot
                        XS_gl%Chi(g_tot+l)=XS_int%Chi(g_tot+l)/temp
                    end do
                else
                    do g=1, Opt_gl%n_g
                        g_tot=(g-1)*Opt_gl%n_tot
                        XS_gl%Chi(g_tot+l)=0.0_dprec
                    end do
                end if
            end do
            
            ! Diffusion coefficient 
            XS_gl%Dif=1.0_dprec/3.0_dprec/XS_gl%Tra
            
            ! Removal coefficient 
            do g=1, Opt_gl%n_g
                g_tot=(g-1)*Opt_gl%n_tot
                XS_gl%Rem(g_tot+1:g_tot+Opt_gl%n_tot)=XS_gl%Absr(g_tot+1:g_tot+Opt_gl%n_tot)
                do gg=1, Opt_gl%n_g
                    if (gg .NE. g) XS_gl%Rem(g_tot+1:g_tot+Opt_gl%n_tot)=XS_gl%Rem(g_tot+1:g_tot+Opt_gl%n_tot)+XS_gl%Scatt(g_tot+1:g_tot+Opt_gl%n_tot, gg)
                end do
            end do
            
    
        end subroutine XS_HomoCond
    
end module HGCMFD_XS_Homo_Cond