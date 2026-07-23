module pCMFD_Alg
    use precision_kinds, only: prec
    use Variables
    use GaussElimination
    use Service_Fcns
        implicit none
    
    contains
        
        subroutine pCMFD_Header(Elem_CMFD, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo) 
            
            type(Figure), intent(in) :: Elem_CMFD(:)
            type(XS_Data), intent(in) :: XS_lo
            type(XS_Data), intent(inout) :: XS_gl
            type(Options_Data), intent(in) :: Opt_lo, Opt_gl
            
            type(LO_geom), allocatable, intent(in) :: LO_Mesh(:)
            type(LO_coeff), allocatable, intent(in) :: LO_Param(:)
            type(GL_geom), allocatable, intent(in) :: GL_Mesh(:)
            type(GL_coeff), allocatable, intent(inout) :: GL_Param(:)
            type(Accel_Vars_Vect), intent(inout) :: Serv_Vect
            type(Accel_Vars_Matr), intent(inout) :: Serv_Matr
            
            real(prec), allocatable, intent(inout) :: Phi_lo(:)
            
            
            integer :: i, j, l, u, t, g, t_tot, g_tot, iter_out_gl, iter_in_gl
            real(prec) :: err0, err1, err2, k_old_gl, k_gl, Phi_avg_lo(Opt_gl%n_g*Opt_gl%n_tot)
            real(prec), allocatable :: eigens_R(:), eigens_I(:), Phi_lo_temp(:)
            
            integer :: ll, t_red
            !real(prec) :: tempRR(Opt_gl%n_g*Opt_gl%n_tot)
            
            
            ! Current calculation
            call Currents_lo_to_gl_DerType(2, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param)
            
            ! Flux-Volume-weighted homogenization
            call Homogenization_DerType(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, XS_lo, XS_gl, Phi_lo, Serv_Vect%Phi)
            Serv_Matr%Phi=Serv_Vect%Phi
            
            ! Saving the homogenized flux
            Serv_Vect%Phi_avg=Serv_Vect%Phi
            Serv_Matr%Phi_avg=Serv_Matr%Phi
            
            ! Computation of the interface coefficients (global)
            call D_tilde_Build(Opt_gl, GL_Mesh, XS_gl%Dif, GL_Param)
                        
            ! Computation of the correction coefficients (global)
            call D_Hat_p_Build(Opt_gl, Serv_Vect%Phi_avg(:), GL_Mesh, GL_Param)
            call D_Hat_m_Build(Opt_gl, Serv_Vect%Phi_avg(:), GL_Mesh, GL_Param)
            
            ! Computation of the migration matrix (global)
            call MigrMat_Impl(Opt_gl, XS_gl%Rem, GL_Mesh, GL_Param, Serv_Vect%a, Serv_Vect%b, Serv_Vect%c, Serv_Vect%d, Serv_Vect%e)
            call MigrMat_Expl(Opt_gl, XS_gl%Rem, GL_Mesh, GL_Param, Serv_Matr%MM)
            
            
            ! Initialization 
            err0=1.e0_prec
            err1=1.e0_prec
            err2=1.e0_prec
            k_gl=1.e0_prec
            iter_out_gl=0
            
            ! Outer global iterations 
            do while ((err1>Opt_gl%tol1 .OR. err2>Opt_gl%tol2) .AND. iter_out_gl<Opt_gl%it_out_max) 
                iter_out_gl=iter_out_gl+1
                k_old_gl=k_gl
                Serv_Vect%Phi_old=Serv_Vect%Phi
                
                ! Fission Source building 
                Serv_Vect%S_fiss_old=0.e0_prec
                do t=1, Opt_gl%n_g
                    t_tot=(t-1)*Opt_gl%n_tot
                    Serv_Vect%S_fiss_old(1:Opt_gl%n_tot)=Serv_Vect%S_fiss_old(1:Opt_gl%n_tot)+Serv_Vect%Phi_old(t_tot+1:t_tot+Opt_gl%n_tot)*XS_gl%nuFis(t_tot+1:t_tot+Opt_gl%n_tot)
                end do
                
                ! Inner iterations 
                iter_in_gl=0
                err0=1.e0_prec
                do while (err0>Opt_gl%tol0 .AND. iter_in_gl<Opt_gl%it_in_max)
                    iter_in_gl=iter_in_gl+1
                    Serv_Vect%Phi_temp=Serv_Vect%Phi
                    
                    do t=1, Opt_gl%n_g ! Group sweeping 
                        t_tot=(t-1)*Opt_gl%n_tot
                        
                        ! Total source building
                        Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot)=XS_gl%Chi(t_tot+1:t_tot+Opt_gl%n_tot)/k_old_gl*Serv_Vect%S_fiss_old(1:Opt_gl%n_tot)
                        do g=1, Opt_gl%n_g
                            g_tot=(g-1)*Opt_gl%n_tot
                            if (t.NE.g) Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot)=Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot)+XS_gl%Scatt(g_tot+1:g_tot+Opt_gl%n_tot,t)*Serv_Vect%Phi(g_tot+1:g_tot+Opt_gl%n_tot)
                        end do
                        
                        
                        call GE_Expl_Main(Opt_gl%n_tot, Serv_Matr%MM(t_tot+1:t_tot+Opt_gl%n_tot,t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%Phi(t_tot+1:t_tot+Opt_gl%n_tot))
                        !call Schur_Complement_Vect(Opt_gl%n_tot, Serv_Vect%a(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%b(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%c(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%Phi(t_tot+1:t_tot+Opt_gl%n_tot))
                        ! Be careful in the a and c arrays!!
                        
                        
                    end do
                    
                    err0=sqrt(dot_product(Serv_Vect%Phi-Serv_Vect%Phi_temp, Serv_Vect%Phi-Serv_Vect%Phi_temp)/dot_product(Serv_Vect%Phi, Serv_Vect%Phi))
                end do
                
                ! Fission Source building 
                Serv_Vect%S_fiss=0.e0_prec
                do t=1, Opt_gl%n_g
                    t_tot=(t-1)*Opt_gl%n_tot
                    Serv_Vect%S_fiss(1:Opt_gl%n_tot)=Serv_Vect%S_fiss(1:Opt_gl%n_tot)+Serv_Vect%Phi(t_tot+1:t_tot+Opt_gl%n_tot)*XS_gl%nuFis(t_tot+1:t_tot+Opt_gl%n_tot)
                end do
                
                ! Error computation 
                k_gl=k_old_gl*dot_product(Serv_Vect%S_fiss, Serv_Vect%S_fiss)/dot_product(Serv_Vect%S_fiss, Serv_Vect%S_fiss_old)
                err1=abs(k_gl-k_old_gl)/k_old_gl
                err2=sqrt(dot_product(Serv_Vect%S_fiss-Serv_Vect%S_fiss_old, Serv_Vect%S_fiss-Serv_Vect%S_fiss_old)/dot_product(Serv_Vect%S_fiss, Serv_Vect%S_fiss))
                
            end do
            
            
            call Flux_Vol_Normalization(Opt_gl%n_g, [(GL_Mesh(i)%V_gl, i=1, Opt_gl%n_tot)], Serv_Vect%Phi)
            
            ! Modulation 
            if (product(Serv_Vect%Phi)>0) then          ! .AND. .FALSE.) then
                if (Opt_lo%flag_SOR==1) then
                    allocate(Phi_lo_temp(size(Phi_lo)))
                    Phi_lo_temp=Phi_lo
                end if
                
                call Modulation(Opt_lo, Opt_gl, LO_Mesh, LO_Param, Serv_Vect%Phi, Serv_Vect%Phi_avg, Phi_lo)
                ! SOR Acceleration
                if (Opt_lo%flag_SOR==1) then
                    call SOR_Accel(Opt_lo%w_SOR, Phi_lo_temp, Phi_lo)
                end if
                
            end if
        
        end subroutine  pCMFD_Header
        
        
        
        
        subroutine D_tilde_Build(Opt_gl, GL_Mesh, D_gl, GL_Param) 
            real(prec), intent(in) :: D_gl(:)
            
            type(Options_Data), intent(in) :: Opt_gl
            type(GL_geom), intent(in) :: GL_Mesh(:)
            type(GL_coeff), intent(inout) :: GL_Param(:)
            
            integer :: l, t, n, t_tot, nb, nb_side
            
            do l=1, size(GL_Param)
                GL_Param(l)%D_Til=0.e0_prec
            end do
            
            do l=1, Opt_gl%n_tot 
                
                do t=1, Opt_gl%n_g
                    t_tot=(t-1)*Opt_gl%n_tot
                    
                    do n=1, size(GL_Param(t_tot+l)%D_Til,1)  ! Iteration on the number of sides
                        nb=GL_Mesh(l)%Neigh_ID(n)
                        if (nb>0) then
                            nb_side=GL_Mesh(l)%Sides_ID(n)
                            
                            GL_Param(t_tot+l)%D_til(n,1)= (D_gl(t_tot+l)/GL_Mesh(l)%dl_gl(n) * D_gl(t_tot+nb)/GL_Mesh(l)%dl_gl(nb_side)) / (D_gl(t_tot+l)/GL_Mesh(l)%dl_gl(n) + D_gl(t_tot+nb)/GL_Mesh(l)%dl_gl(nb_side))!/2.e0_prec
                        else
                            GL_Param(t_tot+l)%D_til(n,1)= D_gl(t_tot+l)!/GL_Mesh(l)%dl_gl(n)
                        end if
                    end do
                    
                end do
            end do
        
        end subroutine D_tilde_Build
        
        subroutine D_Hat_p_Build(Opt_gl, Phi_Ref, GL_Mesh, GL_Param) 
            real(prec), intent(in) :: Phi_Ref(:)
            
            type(Options_Data), intent(in) :: Opt_gl
            type(GL_geom), intent(in) :: GL_Mesh(:)
            type(GL_coeff), intent(inout) :: GL_Param(:)
                    
            integer :: t, l, t_tot, no_faces, nb, n
            
            no_faces=size(GL_Mesh(1)%A_gl)
            
            do l=1, Opt_gl%n_tot   ! Rolling on the meshes
                do t=1,Opt_gl%n_g
                    t_tot=(t-1)*Opt_gl%n_tot
                    do n=1, no_faces
                        nb=GL_Mesh(l)%Neigh_ID(n)
                        
                        if (nb>0) then
                            GL_Param(t_tot+l)%D_hat(n,1)=(GL_Param(t_tot+l)%J_part(n,1)+GL_Param(t_tot+l)%D_til(n,1)*(Phi_Ref(t_tot+nb)-Phi_Ref(t_tot+l)))/Phi_Ref(t_tot+l)
                        else
                            GL_Param(t_tot+l)%D_hat(n,1)=GL_Param(t_tot+l)%J_part(n,1)/Phi_Ref(t_tot+l)
                        end if
                        
                    end do
                end do
            end do
            
        end subroutine D_Hat_p_Build
        
        subroutine D_Hat_m_Build(Opt_gl, Phi_Ref, GL_Mesh, GL_Param) 
            real(prec), intent(in) :: Phi_Ref(:)
            
            type(Options_Data), intent(in) :: Opt_gl
            type(GL_geom), intent(in) :: GL_Mesh(:)
            type(GL_coeff), intent(inout) :: GL_Param(:)
                    
            integer :: t, l, t_tot, no_faces, nb, n
            
            no_faces=size(GL_Mesh(1)%A_gl)
            
            do l=1, Opt_gl%n_tot   ! Rolling on the meshes
                do t=1,Opt_gl%n_g
                    t_tot=(t-1)*Opt_gl%n_tot
                    do n=1, no_faces
                        nb=GL_Mesh(l)%Neigh_ID(n)
                        
                        if (nb>0) then
                            GL_Param(t_tot+l)%D_hat(n,2)=(GL_Param(t_tot+l)%J_part(n,2)-GL_Param(t_tot+l)%D_til(n,1)*(Phi_Ref(t_tot+nb)-Phi_Ref(t_tot+l)))/Phi_Ref(t_tot+nb)
                        else
                            GL_Param(t_tot+l)%D_hat(n,2)=GL_Param(t_tot+l)%J_part(n,2)/Phi_Ref(t_tot+l)
                        end if
                        
                    end do
                end do
            end do
            
            
        end subroutine D_Hat_m_Build
        
        subroutine MigrMat_Impl(Opt_gl, XS_RM, GL_Mesh, Gl_Param, a_MM, b_MM, c_MM, d_MM, e_MM) 
            real(prec), intent(in) :: XS_RM(:)
            
            type(Options_data), intent(in) :: Opt_gl
            type(GL_geom), intent(in) :: GL_Mesh(:)
            type(GL_coeff), intent(in) :: GL_Param(:)
            
            real(prec), intent(out) :: a_MM(:), b_MM(:), c_MM(:), d_MM(:), e_MM(:)
            
            integer :: t, l, t_tot, i, j
            
            
            a_MM=0.e0_prec
            c_MM=0.e0_prec
            d_MM=0.e0_prec
            e_MM=0.e0_prec
            
            b_MM=XS_RM
            
            do t=1, Opt_gl%n_g
                t_tot=(t-1)*Opt_gl%n_tot
                
                do l=1, Opt_gl%n_tot
                    i=mod(l-1, Opt_gl%n_x)+1
                    j=int((l-i)/Opt_gl%n_x)+1
                    
                    if (j>1) then
                        d_MM(t_tot+l)=             -(2*GL_Param(t_tot+l)%D_til(1,1)+GL_Param(t_tot+l)%D_hat(1,2))*GL_Mesh(l)%A_gl(1)/GL_Mesh(l)%V_gl
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(2*GL_Param(t_tot+l)%D_til(1,1)+GL_Param(t_tot+l)%D_hat(1,1))*GL_Mesh(l)%A_gl(1)/GL_Mesh(l)%V_gl
                    else
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(GL_Param(t_tot+l)%D_hat(1,1)-GL_Param(t_tot+l)%D_hat(1,2))*GL_Mesh(l)%A_gl(1)/GL_Mesh(l)%V_gl
                    end if
                    
                    if (i<Opt_gl%n_x) then
                        c_MM(t_tot+l)=             -(2*GL_Param(t_tot+l)%D_til(2,1)+GL_Param(t_tot+l)%D_hat(2,2))*GL_Mesh(l)%A_gl(2)/GL_Mesh(l)%V_gl
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(2*GL_Param(t_tot+l)%D_til(2,1)+GL_Param(t_tot+l)%D_hat(2,1))*GL_Mesh(l)%A_gl(2)/GL_Mesh(l)%V_gl
                    else
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(GL_Param(t_tot+l)%D_hat(2,1)-GL_Param(t_tot+l)%D_hat(2,2))*GL_Mesh(l)%A_gl(2)/GL_Mesh(l)%V_gl
                    end if
                    
                    if (j<Opt_gl%n_y) then
                        e_MM(t_tot+l)=             -(2*GL_Param(t_tot+l)%D_til(3,1)+GL_Param(t_tot+l)%D_hat(3,2))*GL_Mesh(l)%A_gl(3)/GL_Mesh(l)%V_gl
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(2*GL_Param(t_tot+l)%D_til(3,1)+GL_Param(t_tot+l)%D_hat(3,1))*GL_Mesh(l)%A_gl(3)/GL_Mesh(l)%V_gl
                    else
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(GL_Param(t_tot+l)%D_hat(3,1)-GL_Param(t_tot+l)%D_hat(3,2))*GL_Mesh(l)%A_gl(3)/GL_Mesh(l)%V_gl
                    end if
                    
                    if (i>1) then
                        a_MM(t_tot+l)=             -(2*GL_Param(t_tot+l)%D_til(4,1)+GL_Param(t_tot+l)%D_hat(4,2))*GL_Mesh(l)%A_gl(4)/GL_Mesh(l)%V_gl
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(2*GL_Param(t_tot+l)%D_til(4,1)+GL_Param(t_tot+l)%D_hat(4,1))*GL_Mesh(l)%A_gl(4)/GL_Mesh(l)%V_gl
                    else
                        b_MM(t_tot+l)=b_MM(t_tot+l)+(GL_Param(t_tot+l)%D_hat(4,1)-GL_Param(t_tot+l)%D_hat(4,2))*GL_Mesh(l)%A_gl(4)/GL_Mesh(l)%V_gl
                    end if
                    
                end do
                
            end do
            
        end subroutine MigrMat_Impl
        
        subroutine MigrMat_Expl(Opt_gl, XS_RM, GL_Mesh, GL_Param, MM) 
            real(prec), intent(in) :: XS_RM(:)
            
            type(Options_Data), intent(in) :: Opt_gl
            type(GL_geom), intent(in) :: GL_Mesh(:)
            type(GL_coeff), intent(in) :: GL_Param(:)
            
            real(prec), intent(out) :: MM(:,:)
            
            integer :: i, j, t, l, t_tot
            
            
            MM=0.e0_prec
            
            do t=1, Opt_gl%n_g
                t_tot=(t-1)*Opt_gl%n_tot
                
                do l=1, Opt_gl%n_tot
                    MM(t_tot+l,t_tot+l)=XS_RM(t_tot+l)
                    
                    i=mod(l-1, Opt_gl%n_x)+1
                    j=int((l-i)/Opt_gl%n_x)+1
                    
                    if (j>1) then
                        MM(t_tot+l,t_tot+l-Opt_gl%n_x)=        -(2*GL_Param(t_tot+l)%D_til(1,1)+GL_Param(t_tot+l)%D_hat(1,2))*GL_Mesh(l)%A_gl(1)/GL_Mesh(l)%V_gl
                        MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(2*GL_Param(t_tot+l)%D_til(1,1)+GL_Param(t_tot+l)%D_hat(1,1))*GL_Mesh(l)%A_gl(1)/GL_Mesh(l)%V_gl
                    else
                        MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(GL_Param(t_tot+l)%D_hat(1,1)-GL_Param(t_tot+l)%D_hat(1,2))*GL_Mesh(l)%A_gl(1)/GL_Mesh(l)%V_gl
                    end if
                    
                    if (i<Opt_gl%n_x) then
                        MM(t_tot+l,t_tot+l+1)=                 -(2*GL_Param(t_tot+l)%D_til(2,1)+GL_Param(t_tot+l)%D_hat(2,2))*GL_Mesh(l)%A_gl(2)/GL_Mesh(l)%V_gl
                        MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(2*GL_Param(t_tot+l)%D_til(2,1)+GL_Param(t_tot+l)%D_hat(2,1))*GL_Mesh(l)%A_gl(2)/GL_Mesh(l)%V_gl
                    else
                        MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(GL_Param(t_tot+l)%D_hat(2,1)-GL_Param(t_tot+l)%D_hat(2,2))*GL_Mesh(l)%A_gl(2)/GL_Mesh(l)%V_gl
                    end if
                    
                    if (j<Opt_gl%n_y) then
                         MM(t_tot+l,t_tot+l+Opt_gl%n_x)=        -(2*GL_Param(t_tot+l)%D_til(3,1)+GL_Param(t_tot+l)%D_hat(3,2))*GL_Mesh(l)%A_gl(3)/GL_Mesh(l)%V_gl
                         MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(2*GL_Param(t_tot+l)%D_til(3,1)+GL_Param(t_tot+l)%D_hat(3,1))*GL_Mesh(l)%A_gl(3)/GL_Mesh(l)%V_gl
                    else
                         MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(GL_Param(t_tot+l)%D_hat(3,1)-GL_Param(t_tot+l)%D_hat(3,2))*GL_Mesh(l)%A_gl(3)/GL_Mesh(l)%V_gl
                    end if
                    
                    if (i>1) then
                        MM(t_tot+l,t_tot+l-1)=                 -(2*GL_Param(t_tot+l)%D_til(4,1)+GL_Param(t_tot+l)%D_hat(4,2))*GL_Mesh(l)%A_gl(4)/GL_Mesh(l)%V_gl
                        MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(2*GL_Param(t_tot+l)%D_til(4,1)+GL_Param(t_tot+l)%D_hat(4,1))*GL_Mesh(l)%A_gl(4)/GL_Mesh(l)%V_gl
                    else
                        MM(t_tot+l,t_tot+l)=MM(t_tot+l,t_tot+l)+(GL_Param(t_tot+l)%D_hat(4,1)-GL_Param(t_tot+l)%D_hat(4,2))*GL_Mesh(l)%A_gl(4)/GL_Mesh(l)%V_gl
                    end if
                    
                end do
            end do
            
            
        end subroutine MigrMat_Expl
    
    
        
        
        
        
end module pCMFD_Alg
