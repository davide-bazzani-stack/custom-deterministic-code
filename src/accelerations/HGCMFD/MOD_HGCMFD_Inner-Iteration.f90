module HGCMFD_Inner_Iteration
    
    use HGCMFD_Variables, only: dprec, Options_Data_2, XS_Data_2, Accel_Vars_Vect_2, Accel_Vars_Matr_2
    use HGCMFD_S_tot
    use HGCMFD_Discrepancies
    
    use HGCMFD_Gauss_Elimination
    use HGCMFD_BiCGSTAB
    use HGCMFD_Shur_Complement
    
        implicit none
    
    contains
    
        subroutine HGCMFD_Inner_Iteration_Expl(k_old_gl, Opt_gl, XS_gl, Serv_Matr, iter)
            
            real(dprec), intent(in) :: k_old_gl
            type(Options_Data_2), intent(in) :: Opt_gl
            type(XS_Data_2), intent(in) :: XS_gl
            type(Accel_Vars_Matr_2), intent(inout) :: Serv_Matr
        
            integer, intent(out) :: iter
            integer :: t_tot, t
            real(dprec) :: err0
            
            ! Initialization
            iter=0
            err0=1.0_dprec
            
            ! Inner iterations
            do while (err0>Opt_gl%tol0 .AND. iter<Opt_gl%it_in_max)
                
                ! Update with prev. values
                iter=iter+1
                Serv_Matr%Phi_temp=Serv_Matr%Phi
                
                ! Group sweeping
                do t=1, Opt_gl%n_g 
                    t_tot=(t-1)*Opt_gl%n_tot
                    
                    ! Total source building
                    Serv_Matr%source(t_tot+1:t_tot+Opt_gl%n_tot)=S_tot(Opt_gl%n_tot, Opt_gl%n_g, t, XS_gl, k_old_gl, Serv_Matr%S_fiss_old, Serv_Matr%Q_ext(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Matr%Phi)
                    
                    ! Linear system solver
                    call GE_Expl_Main(Opt_gl%n_tot, Serv_Matr%MM(t_tot+1:t_tot+Opt_gl%n_tot,t_tot+1:t_tot+Opt_gl%n_tot), Serv_Matr%source(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Matr%Phi(t_tot+1:t_tot+Opt_gl%n_tot))
                    
                end do
                
                err0=L2_norm_error(Serv_Matr%Phi-Serv_Matr%Phi_temp, Serv_Matr%Phi)
            end do
            
        end subroutine HGCMFD_Inner_Iteration_Expl
        
        
        subroutine HGCMFD_Inner_Iteration_Impl(k_old_gl, Opt_gl, XS_gl, Serv_Vect, iter)
            
            real(dprec), intent(in) :: k_old_gl
            type(Options_Data_2), intent(in) :: Opt_gl
            type(XS_Data_2), intent(in) :: XS_gl
            type(Accel_Vars_Vect_2), intent(inout) :: Serv_Vect
        
            integer, intent(out) :: iter
            integer :: t_tot, t, info
            real(dprec) :: err0
            
            ! Initialization
            iter=0
            err0=1.0_dprec
            
            ! Inner iterations
            do while (err0>Opt_gl%tol0 .AND. iter<Opt_gl%it_in_max)
                
                ! Update with prev. values
                iter=iter+1
                Serv_Vect%Phi_temp=Serv_Vect%Phi
                
                ! Group sweeping
                do t=1, Opt_gl%n_g 
                    t_tot=(t-1)*Opt_gl%n_tot
                    
                    ! Total source building
                    Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot)=S_tot(Opt_gl%n_tot, Opt_gl%n_g, t, XS_gl, k_old_gl, Serv_Vect%S_fiss_old, Serv_Vect%Q_ext(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%Phi)
                    
                    ! Linear system solver
                    call BiCGSTAB_HGCMFD(Opt_gl%n_tot, Opt_gl%it_sol_max, Opt_gl%tol_solv, Serv_Vect%a(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%b(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%c(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%Phi(t_tot+1:t_tot+Opt_gl%n_tot), info)
                    !call Schur_Complement_Vect(Opt_gl%n_tot, Serv_Vect%a(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%b(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%c(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%source(t_tot+1:t_tot+Opt_gl%n_tot), Serv_Vect%Phi(t_tot+1:t_tot+Opt_gl%n_tot))
                    !
                    !if (info==1) then 
                    !    write(*,'(A, I6.1, A)') 'ERROR! THE ', t,' GROUP IN THE HGCMFD DID NOT CONVERGE'
                    !    write(1,'(A, I6.1, A)') 'ERROR! THE ', t,' GROUP IN THE HGCMFD DID NOT CONVERGE'
                    !end if
                    
                end do
                
                err0=L2_norm_error(Serv_Vect%Phi-Serv_Vect%Phi_temp, Serv_Vect%Phi)
            end do
            
        end subroutine HGCMFD_Inner_Iteration_Impl
        
end module HGCMFD_Inner_Iteration