module HGCMFD_Normalization
    
    use HGCMFD_Variables, only: dprec
    
        implicit none
    
    contains
    
        subroutine Normalize(n_g, V, Phi) 
            integer, intent(in) :: n_g
            real(dprec), intent(in) :: V(:)
            real(dprec), intent(inout) :: Phi(:)
            
            integer :: g, g_red, n
            real(dprec) :: tmp
            
            n=size(V)
            tmp=0.0_dprec
            
            do g=1,n_g
                g_red=(g-1)*n
                tmp=tmp+sum(Phi(g_red+1:g_red+n)*V)
            end do
            tmp=tmp/sum(V)
            
            Phi=Phi/tmp
            
        end subroutine Normalize
    
    
end module HGCMFD_Normalization