module HGCMFD_Discrepancies
    
    use HGCMFD_Variables, only: dprec
    
        implicit none
    
    contains
    
    
    
    pure function abs_rel_err(x_new, x_old) result(discr)
        real(dprec), intent(in) :: x_new, x_old
        real(dprec) :: discr
    
        discr=abs(x_new-x_old)/x_old
    
    end function abs_rel_err
    
    
    
    pure function L2_norm_error(delta, array) result(normL2)
        real(dprec), intent(in) :: delta(:), array(:)
        real(dprec) :: normL2
    
        normL2=sqrt(dot_product(delta, delta)/dot_product(array, array))
    
    end function L2_norm_error
    
end module HGCMFD_Discrepancies