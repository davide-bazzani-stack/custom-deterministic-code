module HGCMFD_S_tot
    
    use HGCMFD_Variables, only: dprec, XS_Data_2
        implicit none
    
    contains 
    
    pure function S_tot(n_tot, n_g, t, XS_gl, k_gl, S_fiss, Q_ext, Phi) result(source_tot)
        integer, intent(in) :: n_tot, n_g, t 
        type(XS_Data_2), intent(in) :: XS_gl
        real(dprec), intent(in) :: k_gl    ! Scalar
        real(dprec), intent(in) :: S_fiss(:), Phi(:) ! All-groups size
        real(dprec), intent(in) :: Q_ext(:)    ! Single group size (t group)
        real(dprec), allocatable :: source_tot(:)
        
        integer :: g_tot, g, t_tot
        
        allocate(source_tot(n_tot))
        t_tot=(t-1)*n_tot
        source_tot(1:n_tot)=XS_gl%Chi(t_tot+1:t_tot+n_tot)/k_gl*S_fiss(1:n_tot)+Q_ext(1:n_tot)
        
        do g=1, n_g
            if (t .NE. g) then  ! Exclusion of self-scattering
                g_tot=(g-1)*n_tot
                source_tot(1:n_tot)=source_tot(1:n_tot)+XS_gl%Scatt(g_tot+1:g_tot+n_tot,t)*Phi(g_tot+1:g_tot+n_tot)
            end if
        end do
        
        
    end function S_tot
    
end module HGCMFD_S_tot