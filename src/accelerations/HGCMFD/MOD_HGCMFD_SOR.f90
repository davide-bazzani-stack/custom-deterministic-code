module HGCMFD_SOR_Update
    
    use HGCMFD_Variables, only: dprec
    
        implicit none
    
    contains
        subroutine SOR_Accel(w_SOR, Var_Pre, Final_Val) 
            real(dprec), intent(in) :: w_SOR, Var_Pre(:)
            real(dprec), intent(inout) :: Final_Val(:)
            
            Final_Val=(1-w_SOR)*Var_Pre+w_SOR*Final_Val
            
        end subroutine SOR_Accel
end module HGCMFD_SOR_Update