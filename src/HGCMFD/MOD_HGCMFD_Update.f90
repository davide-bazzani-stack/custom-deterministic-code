module HGCMFD_Update
    
    use HGCMFD_Variables, only: dprec
    
        implicit none
    
    contains
    
    pure function k_update(k_old, S_fiss_old, S_fiss) result(k_new)
        real(dprec), intent(in) :: k_old, S_fiss_old(:), S_fiss(:)
        real(dprec) :: k_new
    
        k_new=k_old*dot_product(S_fiss, S_fiss)/dot_product(S_fiss, S_fiss_old)
    
    end function k_update
    
    
end module HGCMFD_Update