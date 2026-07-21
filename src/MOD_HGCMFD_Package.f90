!==========================================================================================================================
!   
!   HG-CMFD Module
!   
!   
!--> The following module contains the NewCMFD algorithm.
!   
!==========================================================================================================================  
module HGCMFD_Alg
    use precision_kinds, only: prec
    use Variables
    use GaussElimination
    use Service_Fcns
        implicit none
        
    contains
    
    subroutine HGCMFD_Header(Elem_CMFD, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo) 
        
        type(Figure), intent(in) :: Elem_CMFD(:)
        type(XS_Data), intent(in) :: XS_lo
        type(XS_Data), intent(inout) :: XS_gl
        type(Options_Data), intent(in) :: Opt_lo, Opt_gl
        
        type(LO_geom), allocatable, intent(in) :: LO_Mesh(:)
        type(LO_coeff), allocatable, intent(in) :: LO_Param(:)
        type(GL_geom), allocatable, intent(inout) :: GL_Mesh(:)
        type(GL_coeff), allocatable, intent(inout) :: GL_Param(:)
        type(Accel_Vars_Vect), intent(inout) :: Serv_Vect
        type(Accel_Vars_Matr), intent(inout) :: Serv_Matr
        
        real(prec), allocatable, intent(inout) :: Phi_lo(:)
        
        
        integer :: i, j, l, u, t, g, t_tot, g_tot, iter_out_gl, iter_in_gl, info
        real(prec) :: err0, err1, err2, k_old_gl, k_gl, Phi_avg_lo(Opt_gl%n_g*Opt_gl%n_tot)
        real(prec), allocatable :: eigens_R(:), eigens_I(:), Phi_lo_temp(:)
        
        integer :: ll, t_red
        real(prec) :: tempRR(Opt_gl%n_g*Opt_gl%n_tot)
        
        
        ! Current calculation
        call Currents_lo_to_gl_DerType(4, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param)

        ! Flux-Volume-weighted homogenization
        call Homogenization_DerType(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, XS_lo, XS_gl, Phi_lo, Serv_Vect%Phi)
        Serv_Matr%Phi=Serv_Vect%Phi
        
        ! Saving the homogenized flux
        Serv_Vect%Phi_avg=Serv_Vect%Phi
        Serv_Matr%Phi_avg=Serv_Matr%Phi
        
        ! Computation of the single current(global, element-wise)
        call J_gl_to_El(Opt_gl%n_tot, Opt_gl%n_g, GL_Mesh, GL_Param)
        
        
        
!! Optimal dl calculation   
!if (GL_Mesh(1)%dl_gl(1)==0.0) then
!    call Optim_dl(Opt_gl, XS_gl%Dif, GL_Mesh, GL_Param, Serv_Vect%Phi) 
!end if
        
        
        
        ! Computation of the interface coefficients (global)
        call D_tilde_Build(Opt_gl, Elem_CMFD, GL_Mesh, XS_gl%Dif, GL_Param)
        
!if (.FALSE.) then
!    call D_tilde_Mod(Opt_gl, GL_Param, Serv_Vect%Phi) 
!end if
        
        ! Computation of the correction coefficients (global)
        call D_Hat_Build(Opt_gl%n_g, Opt_gl%n_tot, Serv_Vect%Phi_avg(:), GL_Mesh, GL_Param)
        call D_Hat_Build(Opt_gl%n_g, Opt_gl%n_tot, Serv_Matr%Phi_avg(:), GL_Mesh, GL_Param)
        
        ! Computation of the migration matrix (global)
        call MigrMat_Impl(Opt_gl%n_g, Opt_gl%n_tot, XS_gl%Rem, GL_Mesh, GL_Param, Serv_Vect%a, Serv_Vect%b, Serv_Vect%c)
        call MigrMat_Expl(Opt_gl%n_g, Opt_gl%n_tot, XS_gl%Rem, GL_Mesh, GL_Param, Serv_Matr%MM)
        
        
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
                    !call BiCGSTAB_HGCMFD(Opt_gl%n_tot, Opt_gl%it_sol_max, Opt_gl%tol_solv, Serv_Vect%a(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%b(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%c(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%Phi(t_tot+1:t_tot+Opt_gl%n_tot), info)
                    !
                    !if (info==1) then 
                    !    write(*,'(A, I6.1, A)') 'ERROR! THE ', t,' GROUP IN THE HGCMFD DID NOT CONVERGE'
                    !    write(1,'(A, I6.1, A)') 'ERROR! THE ', t,' GROUP IN THE HGCMFD DID NOT CONVERGE'
                    !end if
                    
                    
                    
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
        

!! Sensitivity analysis parameters recording
!! D_tilde
!write(19, '(100000ES18.9)') [(GL_Param(i)%D_til_el(1), i=1, Opt_gl%n_tot*Opt_gl%n_g)]
!! D_hat
!write(20, '(100000ES18.9)') [(GL_Param(i)%D_hat_el(1), i=1, Opt_gl%n_tot*Opt_gl%n_g)]
!! Flux Coarse
!write(12, '(100000ES18.9)') Serv_Vect%Phi(:)
!! Flux Reference
!write(21, '(100000ES18.9)') Serv_Vect%Phi_avg(:)
!! Flux_fine
!!write(17, '(1000000ES18.9)') Phi_lo

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
        
    end subroutine  HGCMFD_Header
    
    
    
    
    
    
    subroutine D_tilde_Build(Opt_gl, Elem_det, GL_Mesh, D_gl, GL_Param) 
        real(prec), intent(in) :: D_gl(:)
        
        type(Options_Data), intent(in) :: Opt_gl
        type(Figure), intent(in) :: Elem_det(:)
        type(GL_geom), intent(in) :: GL_Mesh(:)
        type(GL_coeff), intent(inout) :: GL_Param(:)
        
        integer :: l, t, n, t_tot
        
        do l=1, size(GL_Param)
            GL_Param(l)%D_Til=0.e0_prec
            GL_Param(l)%D_Til_El=0.e0_prec
        end do
        
        
        
        do l=2, Opt_gl%n_tot ! Skipping the BG mesh, iteration "on the universe"            !maxval(lo_to_gl)
            !u=Elem_det(l)%ID ! Figure ID of the universe analyzed
            
            do t=1, Opt_gl%n_g
                t_tot=(t-1)*Opt_gl%n_tot
                
                do n=1, size(GL_Param(t_tot+l)%D_Til,1)  ! Iteration on the number of sides
                    GL_Param(t_tot+l)%D_til(n,1)=D_gl(t_tot+l)/GL_Mesh(l)%dl_gl(n)*D_gl(t_tot+1)/GL_Mesh(1)%dl_gl(1)/(D_gl(t_tot+l)/GL_Mesh(l)%dl_gl(n)+D_gl(t_tot+1)/GL_Mesh(1)%dl_gl(1))
                    
                    GL_Param(t_tot+l)%D_til_el(1)=GL_Param(t_tot+l)%D_til_el(1)+GL_Param(t_tot+l)%D_til(n,1)*GL_Mesh(l)%A_gl(n)
                end do
                
                GL_Param(t_tot+l)%D_til_el(1)=GL_Param(t_tot+l)%D_til_el(1)/sum(GL_Mesh(l)%A_gl(:))
            end do
        end do
    
    end subroutine D_tilde_Build
    
    subroutine J_gl_to_El(n_tot_gl, n_g_gl, GL_Mesh, GL_Param) 
        integer, intent(in) :: n_g_gl, n_tot_gl     !, lo_to_gl_V(:)
        
        !type(Figure), intent(in) :: Elem_det(:)
        type(GL_geom), intent(in) :: GL_Mesh(:)
        type(GL_coeff), intent(inout) :: GL_Param(:)
        
        integer :: l, t, n, u, t_tot
        integer :: no_faces
        real(prec) :: area_sum
        
        ! Non-BG Meshes
        do l=2, n_tot_gl ! Skipping the BG mesh, iteration "on the universe"    ! maxval(lo_to_gl_V)
            
            no_faces=size(GL_Mesh(l)%A_gl(:))
            area_sum=sum(GL_Mesh(l)%A_gl(:))
            
            do t=1, n_g_gl
                t_tot=(t-1)*n_tot_gl
                
                GL_Param(t_tot+l)%J_Net_El=dot_product(GL_Param(t_tot+l)%J_Net(1:no_faces), GL_Mesh(l)%A_gl(1:no_faces))/area_sum
            end do
        end do
        
    end subroutine J_gl_to_El
    
    subroutine D_Hat_Build(n_g, n_tot, Phi_Ref, GL_Mesh, GL_Param) 
        integer, intent(in) :: n_g, n_tot
        real(prec), intent(in) :: Phi_Ref(:)
        
        type(GL_geom), intent(in) :: GL_Mesh(:)
        type(GL_coeff), intent(inout) :: GL_Param(:)
                
        integer :: t, l, t_tot, no_faces
        real(prec) :: area_sum
        
        no_faces=size(GL_Mesh(1)%A_gl(:))
        area_sum=sum(GL_Mesh(1)%A_gl(:))
        
        do t=1,n_g
            t_tot=(t-1)*n_tot
            do l=2, n_tot   ! Rolling on the meshes
                GL_Param(t_tot+l)%D_hat_el(1)=(GL_Param(t_tot+l)%J_Net_el+GL_Param(t_tot+l)%D_til_el(1)*(Phi_Ref(t_tot+1)-Phi_Ref(t_tot+l)))/(Phi_Ref(t_tot+1)+Phi_Ref(t_tot+l))
            end do
            
            ! sum(J/Phi*A/A_tot)
            GL_Param(t_tot+1)%D_hat_el(1)=dot_product(GL_Param(t_tot+1)%J_Net(1:no_faces), GL_Mesh(1)%A_gl(1:no_faces))/(area_sum*Phi_Ref(t_tot+1))
        end do
        
        
    end subroutine D_Hat_Build
    
    subroutine MigrMat_Impl(n_g, n_tot, XS_RM, HGCMFD_Mesh, HGCMFD_Param, a_MM, b_MM, c_MM) 
        integer, intent(in) :: n_g, n_tot
        real(prec), intent(in) :: XS_RM(:)
        
        type(GL_geom), intent(in) :: HGCMFD_Mesh(:)
        type(GL_coeff), intent(in) :: HGCMFD_Param(:)
        !type(Accel_Vars_Vect), intent(inout) :: Serv_Vect(:)
        
        real(prec), intent(out) :: a_MM(:), b_MM(:), c_MM(:)
        
        integer :: t, l, t_tot
        real(prec) :: area_sum
        
        
        a_MM=0.e0_prec
        b_MM=0.e0_prec
        c_MM=0.e0_prec
        
        do t=1, n_g
            t_tot=(t-1)*n_tot
            b_MM(t_tot+1)=XS_RM(t_tot+1)    ! BG mesh removal XS
            
            do l=2, n_tot
                area_sum=sum(HGCMFD_Mesh(l)%A_gl(:))
                
                b_MM(t_tot+l)=XS_RM(t_tot+l)+(HGCMFD_Param(t_tot+l)%D_til_el(1)+HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(l)%V_gl      ! Main diagonal
                a_MM(t_tot+l)=              -(HGCMFD_Param(t_tot+l)%D_til_el(1)-HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(l)%V_gl      ! First column (BG coupling)
                c_MM(t_tot+l)=              -(HGCMFD_Param(t_tot+l)%D_til_el(1)+HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(1)%V_gl      ! First Row (BG eqn.)
                
                b_MM(t_tot+1)= b_MM(t_tot+1)+(HGCMFD_Param(t_tot+l)%D_til_el(1)-HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(1)%V_gl      ! BG diagonal contribution
            end do
            
            ! Adding the boundaries
            area_sum=sum(HGCMFD_Mesh(1)%A_gl(:))
            b_MM(t_tot+1)=b_MM(t_tot+1)+HGCMFD_Param(t_tot+1)%D_hat_el(1)*area_sum/HGCMFD_Mesh(1)%V_gl
        end do
        
    end subroutine MigrMat_Impl
    
    subroutine MigrMat_Expl(n_g, n_tot, XS_RM, HGCMFD_Mesh, HGCMFD_Param, MM) 
        integer, intent(in) :: n_g, n_tot
        real(prec), intent(in) :: XS_RM(:)
        
        type(GL_geom), intent(in) :: HGCMFD_Mesh(:)
        type(GL_coeff), intent(in) :: HGCMFD_Param(:)
        !type(Accel_Vars_Vect), intent(inout) :: Serv_Vect(:)
        
        real(prec), intent(out) :: MM(:,:)
        
        integer :: t, l, t_tot
        real(prec) :: area_sum
        
        
        MM=0.e0_prec
        
        do t=1, n_g
            t_tot=(t-1)*n_tot
            MM(t_tot+1,t_tot+1)=XS_RM(t_tot+1)    ! BG mesh removal XS
            
            do l=2, n_tot
                area_sum=sum(HGCMFD_Mesh(l)%A_gl(:))
                
                MM(t_tot+l,t_tot+l)=     XS_RM(t_tot+l)+(HGCMFD_Param(t_tot+l)%D_til_el(1)+HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(l)%V_gl      ! Main diagonal
                MM(t_tot+l,t_tot+1)=                   -(HGCMFD_Param(t_tot+l)%D_til_el(1)-HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(l)%V_gl      ! First column (BG coupling)
                MM(t_tot+1,t_tot+l)=                   -(HGCMFD_Param(t_tot+l)%D_til_el(1)+HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(1)%V_gl      ! First Row (BG eqn.)
                
                MM(t_tot+1,t_tot+1)=MM(t_tot+1,t_tot+1)+(HGCMFD_Param(t_tot+l)%D_til_el(1)-HGCMFD_Param(t_tot+l)%D_hat_el(1))*area_sum/HGCMFD_Mesh(1)%V_gl      ! BG diagonal contribution
            end do
            
            ! Adding the boundaries
            area_sum=sum(HGCMFD_Mesh(1)%A_gl(:))
            MM(t_tot+1,t_tot+1)=MM(t_tot+1,t_tot+1)+HGCMFD_Param(t_tot+1)%D_hat_el(1)*area_sum/HGCMFD_Mesh(1)%V_gl
        end do
        
    end subroutine MigrMat_Expl
    
    subroutine Optim_dl(Opt_gl, D_gl, GL_Mesh, GL_Param, Phi_ref) 
        real(prec), intent(in) :: D_gl(:), Phi_ref(:)
    
        type(Options_Data), intent(in) :: Opt_gl
        type(GL_geom), intent(inout) :: GL_Mesh(:)
        type(GL_coeff), intent(inout) :: GL_Param(:)
        
        integer :: l, t, n, u, t_tot
        real(prec) :: temp
        real(prec), allocatable :: dl_cand(:)
        
        allocate(dl_cand(Opt_gl%n_tot*Opt_gl%n_g))
        dl_cand(:)=1e-3_prec
        
        
        
        
        do l=2, Opt_gl%n_tot
            do t=1, Opt_gl%n_g
                ! Factor computation D_tilde_avg
                t_tot=(t-1)*Opt_gl%n_tot
                
                temp=(3e-2_prec*(Phi_Ref(t_tot+1)+Phi_ref(t_tot+l))-GL_Param(t_tot+l)%J_Net_el)/(Phi_Ref(t_tot+1)-Phi_ref(t_tot+l))
                
                dl_cand(t_tot+l)=D_gl(t_tot+1)*(1/temp-D_gl(t_tot+1)-GL_Mesh(l)%dl_gl(1)/D_gl(t_tot+l))
                
                !dl_cand(t_tot+l)=D_gl(t_tot+1)*(-(Phi_Ref(t_tot+1)-Phi_ref(t_tot+l))/GL_Param(t_tot+l)%J_Net_el-GL_Mesh(l)%dl_gl(1)/D_gl(t_tot+l))
                
            end do
        end do  
        
        do t=1, Opt_gl%n_g
            t_tot=(t-1)*Opt_gl%n_tot
            
            dl_cand(t_tot+1)=sum(dl_cand(t_tot+2:t*Opt_gl%n_tot))/(Opt_gl%n_tot-1)
        end do
                
        GL_Mesh(1)%dl_gl(:)=dl_cand(1)
    
    end subroutine Optim_dl
    
    
    
    subroutine D_tilde_Mod(Opt_gl, GL_Param, Phi_ref) 
        real(prec), intent(in) :: Phi_ref(:)
    
        type(Options_Data), intent(in) :: Opt_gl
        type(GL_coeff), intent(inout) :: GL_Param(:)
        
        integer :: l, t, n, u, t_tot
        real(prec) :: epsi
        
        
        do l=2, Opt_gl%n_tot
            do t=1, Opt_gl%n_g
                ! Factor computation D_tilde_avg
                t_tot=(t-1)*Opt_gl%n_tot
                
                epsi=(1e-3_prec*(Phi_Ref(t_tot+1)+Phi_ref(t_tot+l))-GL_Param(t_tot+l)%J_Net_el)/(GL_Param(t_tot+l)%D_til_el(1)*(Phi_Ref(t_tot+1)-Phi_ref(t_tot+l)))
                
                GL_Param(t_tot+l)%D_til_el(1)=GL_Param(t_tot+l)%D_til_el(1)*epsi
                
            end do
        end do  
        
    end subroutine D_tilde_Mod
    
    ! Solvers
    
    subroutine Schur_Complement_Vect(n, a_in, b_in, c_in, y_in, x) 
        integer, intent(in) :: n
        integer :: i
        real(prec), intent(in) :: a_in(:), b_in(:), c_in(:), y_in(:)        ! Ax = y
        real(prec), intent(out) :: x(:)
        real(prec) :: factor, S_compl
        
        S_compl=b_in(1)
        x(1)=y_in(1)
        do i=2, n
            S_compl=S_compl-c_in(i)*a_in(i)/b_in(i)
            x(1)=x(1)-c_in(i)*y_in(i)/b_in(i)
        end do
        
        x(1)=x(1)/S_compl
        x(2:n)=(y_in(2:n)-a_in(2:n)*x(1))/b_in(2:n)
        
        return
        
    end subroutine Schur_Complement_Vect
    
    subroutine BiCGSTAB_HGCMFD(n, max_iter, tol, a, b, c, rhs, x, info) 
        integer, intent(in) :: n, max_iter
        real(prec), intent(in) :: a(:), b(:), c(:), rhs(:), tol
        real(prec), intent(inout) :: x(:)
        integer, intent(out) :: info
        real(prec), allocatable :: r(:), r_old(:), v(:), p(:), s(:), t(:), x_old(:)
        real(prec) :: alpha, beta, omega, rho, rho_old, error, norm_b
        integer :: k
                
        allocate(r(n))
        allocate(r_old(n))
        allocate(v(n))
        allocate(p(n))
        allocate(s(n))
        allocate(t(n))
        
        allocate(x_old(n))
        
        
        call matvect_arrow(n, a, b, c, x, r)
        
        r=rhs-r		! b array in Ax=b (known)
        r_old=r
        rho_old=1.e0_prec
        alpha=1.e0_prec	 ! to make beta = 0 and thus p = r for the first iteration
        omega=1.e0_prec
        v=0.e0_prec
        p=0.e0_prec
        
        norm_b=max(sqrt(dot_product(rhs, rhs)), 1.e0_prec)
        
        do k=1,max_iter
            x_old=x
            
            rho=dot_product(r_old, r)
            beta=(rho/rho_old)*(alpha/omega)
            p=r+beta*(p-omega*v)   
            
            call matvect_arrow(n, a, b, c, p, v)
                        
            
            alpha=rho/dot_product(r_old, v)
            s=r-alpha*v       
                    
            error=sqrt(dot_product(s,s))
            
            if (error/norm_b<tol) then
                x=x_old+alpha*p
                info=0
                return
            end if
            
            call matvect_arrow(n, a, b, c, s, t)
            omega=dot_product(t, s)/dot_product(t, t)
            x=x_old+alpha*p+omega*s

            r=s-omega*t
            error=sqrt(dot_product(r,r))
            if (error/norm_b<tol) then
                info=0
                return
            end if
            
            rho_old=rho
        end do
        
        info=1 ! Flag for the missed convergence
        return
    end subroutine BiCGSTAB_HGCMFD
    
    subroutine matvect_arrow(n, a, b, c, x, y) 
        integer, intent(in) :: n
        real(prec), intent(in) :: a(:), b(:), c(:), x(:)
        real(prec), intent(out):: y(:)
        integer :: i
        
        y(1) = b(1)*x(1)
        do i = 2, n
            y(1)=y(1)+c(i)*x(i)
            y(i)=a(i)*x(1)+b(i)*x(i)
        end do
        
    end subroutine matvect_arrow
    
    
end module HGCMFD_Alg
