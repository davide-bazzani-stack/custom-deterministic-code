module HGCMFD_S_fiss
    
    use HGCMFD_Variables, only: dprec
        implicit none
        
    contains
    
    pure function S_fiss(n_tot, n_g, nuFis, Phi) result (FSD)
    
        integer, intent(in) :: n_tot, n_g
        real(dprec), intent(in) :: nuFis(:), Phi(:)
        real(dprec), allocatable :: FSD(:)
        
        integer :: t, t_tot
    
        ! Fission Source initialization 
        allocate(FSD(n_tot))
        FSD=0.0_dprec
        
        ! Fission Source building 
        do t=1, n_g
            t_tot=(t-1)*n_tot
            FSD(1:n_tot)=FSD(1:n_tot)+nuFis(t_tot+1:t_tot+n_tot)*Phi(t_tot+1:t_tot+n_tot)
        end do
        
    end function S_fiss
end module HGCMFD_S_fiss