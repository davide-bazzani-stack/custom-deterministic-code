!==========================================================================================================================
!   
!   HG-CMFD Module
!   
!   
!--> The following module contains the NewCMFD algorithm.
!   
!==========================================================================================================================  
module HGCMFD_Package
    
    use HGCMFD_Variables
    
    use HGCMFD_J_Homo_Cond
    use HGCMFD_XS_Homo_Cond
    use Element_Current
    use HGCMFD_D_tilde
    use HGCMFD_D_hat
    use HGCMFD_Migration_Matrix
    use HGCMDF_Outer_Iteration
    use HGCMFD_Normalization
    use HGCMFD_Modulation
    use HGCMFD_SOR_Update
    
        implicit none
        
    contains
    
    subroutine HGCMFD_Head(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo, k) 
        
        type(XS_Data_2), intent(in) :: XS_lo
        type(XS_Data_2), intent(inout) :: XS_gl
        type(Options_Data_2), intent(in) :: Opt_lo, Opt_gl
        
        type(LO_geom_2), allocatable, intent(in) :: LO_Mesh(:)
        type(LO_coeff_2), allocatable, intent(in) :: LO_Param(:)
        type(GL_geom_2), allocatable, intent(inout) :: GL_Mesh(:)
        type(GL_coeff_2), allocatable, intent(inout) :: GL_Param(:)
        type(Accel_Vars_Vect_2), intent(inout) :: Serv_Vect
        type(Accel_Vars_Matr_2), intent(inout) :: Serv_Matr
        
        real(dprec), allocatable, intent(inout) :: Phi_lo(:)
        real(dprec), intent(inout) :: k
        
        
        integer :: i, j, l, u, t, g, t_tot, g_tot, iter_out, iter_in, info
        real(dprec) :: err0, err1, err2, Phi_avg_lo(Opt_gl%n_g*Opt_gl%n_tot)
        real(dprec), allocatable :: eigens_R(:), eigens_I(:), Phi_lo_temp(:)
        
        integer :: ll, t_red
        real(dprec) :: tempRR(Opt_gl%n_g*Opt_gl%n_tot)
        
        
        ! Current calculation
        call Currents_lo_to_gl(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param)

        ! Flux-Volume-weighted homogenization
        !call XS_HomoCond(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, XS_lo, XS_gl, Phi_lo, Serv_Vect%Phi)
        call XS_HomoCond(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, XS_lo, XS_gl, Phi_lo, Serv_Matr%Phi)
        
        ! Saving the homogenized flux
        !Serv_Vect%Phi_avg=Serv_Vect%Phi
        Serv_Matr%Phi_avg=Serv_Matr%Phi
        
        ! Computation of the single current(global, element-wise)
        call J_gl_to_El(Opt_gl%n_tot, Opt_gl%n_g, GL_Mesh, GL_Param)
        
        
        ! Computation of the interface coefficients (global)
        call D_tilde(Opt_gl, GL_Mesh, XS_gl%Dif, GL_Param)
        
        ! Computation of the correction coefficients (global)
        ! call D_Hat(Opt_gl%n_g, Opt_gl%n_tot, Serv_Vect%Phi_avg(:), GL_Mesh, GL_Param)
        call D_Hat(Opt_gl%n_g, Opt_gl%n_tot, Serv_Matr%Phi_avg(:), GL_Mesh, GL_Param)
        
        ! Computation of the migration matrix (global)
        !call Migration_Vect(Opt_gl%n_g, Opt_gl%n_tot, XS_gl%Rem, GL_Mesh, GL_Param, Serv_Vect%a, Serv_Vect%b, Serv_Vect%c)
        call Migration_Matr(Opt_gl%n_g, Opt_gl%n_tot, XS_gl%Rem, GL_Mesh, GL_Param, Serv_Matr%MM)
        
        ! Initialization 
        !call HGCMFD_Outer_Iteration_Impl(Opt_gl, XS_gl, iter_in, iter_out, k, Serv_Vect)
        call HGCMFD_Outer_Iteration_Expl(Opt_gl, XS_gl, iter_in, iter_out, k, Serv_Matr)
        
        
        ! Normalization
        call Normalize(Opt_gl%n_g, [(GL_Mesh(i)%V_gl, i=1, Opt_gl%n_tot)], Serv_Vect%Phi)
        
        
        ! Modulation 
        if (product(Serv_Vect%Phi)>0) then
            if (Opt_lo%flag_SOR==1) then
                allocate(Phi_lo_temp(size(Phi_lo)))
                Phi_lo_temp=Phi_lo
            end if
            
            !call Modulation(Opt_lo, Opt_gl, LO_Mesh, LO_Param, Serv_Vect%Phi, Serv_Vect%Phi_avg, Phi_lo)
            call Modulation(Opt_lo, Opt_gl, LO_Mesh, LO_Param, Serv_Vect%Phi, Serv_Matr%Phi_avg, Phi_lo)
            
            ! SOR Acceleration
            if (Opt_lo%flag_SOR==1) then
                call SOR_Accel(Opt_lo%w_SOR, Phi_lo_temp, Phi_lo)
            end if
        end if
        
    end subroutine  HGCMFD_Head
    
    
    
end module HGCMFD_Package