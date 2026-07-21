module HGCMDF_Outer_Iteration
    
    use HGCMFD_Variables, only: dprec, XS_Data_2, Options_Data_2, Accel_Vars_Vect_2, Accel_Vars_Matr_2
    
    use HGCMFD_S_fiss
    use HGCMFD_Inner_Iteration
    use HGCMFD_Update
    use HGCMFD_Discrepancies
    
     
        implicit none
    
    contains
        subroutine HGCMFD_Outer_Iteration_Expl(Opt_gl, XS_gl, iter_in, iter_out, k_gl, Serv_Matr)
            
            type(Options_Data_2), intent(in) :: Opt_gl
            type(XS_Data_2), intent(in) :: XS_gl
            
            integer, intent(out) :: iter_in, iter_out
            real(dprec), intent(inout) :: k_gl
            type(Accel_Vars_Matr_2), intent(inout) :: Serv_Matr
            
            real(dprec) :: err1, err2, k_old_gl
            
            err1=1.0_dprec
            err2=1.0_dprec
            iter_out=0
            
            do while ((err1>Opt_gl%tol1 .OR. err2>Opt_gl%tol2) .AND. iter_out<Opt_gl%it_out_max) 
                
                ! Update from prev. Iter.
                iter_out=iter_out+1
                k_old_gl=k_gl
                Serv_Matr%Phi_old=Serv_Matr%Phi
                
                ! Fission Source building 
                Serv_Matr%S_fiss_old=S_fiss(Opt_gl%n_tot, Opt_gl%n_g, XS_gl%nuFis, Serv_Matr%Phi_old)
                
                ! Inner iterations 
                call HGCMFD_Inner_Iteration_Expl(k_old_gl, Opt_gl, XS_gl, Serv_Matr, iter_in)
                
                ! Updated Fission Source building 
                Serv_Matr%S_fiss=S_fiss(Opt_gl%n_tot, Opt_gl%n_g, XS_gl%nuFis, Serv_Matr%Phi)
                
                
                ! Error computation 
                k_gl=k_update(k_old_gl, Serv_Matr%S_fiss_old, Serv_Matr%S_fiss)
                err1=abs_rel_err(k_gl, k_old_gl)
                err2=L2_norm_error(Serv_Matr%S_fiss-Serv_Matr%S_fiss_old, Serv_Matr%S_fiss)
                
            end do
        end subroutine HGCMFD_Outer_Iteration_Expl
        
        
        
        
            
        subroutine HGCMFD_Outer_Iteration_Impl(Opt_gl, XS_gl, iter_in, iter_out, k_gl, Serv_Vect)
            
            type(Options_Data_2), intent(in) :: Opt_gl
            type(XS_Data_2), intent(in) :: XS_gl
            
            integer, intent(out) :: iter_in, iter_out
            real(dprec), intent(inout) :: k_gl
            type(Accel_Vars_Vect_2), intent(inout) :: Serv_Vect
            
            real(dprec) :: err1, err2, k_old_gl
            
            err1=1.0_dprec
            err2=1.0_dprec
            iter_out=0
            
            do while ((err1>Opt_gl%tol1 .OR. err2>Opt_gl%tol2) .AND. iter_out<Opt_gl%it_out_max) 
                
                ! Update from prev. Iter.
                iter_out=iter_out+1
                k_old_gl=k_gl
                Serv_Vect%Phi_old=Serv_Vect%Phi
                
                ! Fission Source building 
                Serv_Vect%S_fiss_old=S_fiss(Opt_gl%n_tot, Opt_gl%n_g, XS_gl%nuFis, Serv_Vect%Phi_old)
                
                ! Inner iterations 
                call HGCMFD_Inner_Iteration_Impl(k_old_gl, Opt_gl, XS_gl, Serv_Vect, iter_in)
                
                ! Updated Fission Source building 
                Serv_Vect%S_fiss=S_fiss(Opt_gl%n_tot, Opt_gl%n_g, XS_gl%nuFis, Serv_Vect%Phi)
                
                
                ! Error computation 
                k_gl=k_update(k_old_gl, Serv_Vect%S_fiss_old, Serv_Vect%S_fiss)
                err1=abs_rel_err(k_gl, k_old_gl)
                err2=L2_norm_error(Serv_Vect%S_fiss-Serv_Vect%S_fiss_old, Serv_Vect%S_fiss)
                
            end do
        end subroutine HGCMFD_Outer_Iteration_Impl
    
end module HGCMDF_Outer_Iteration