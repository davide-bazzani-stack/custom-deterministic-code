module Intermediates
    use iso_fortran_env, only: output_unit, error_unit
    use precision_kinds, only: prec
    use IO_module, only: files, output_dir, debug_mark, file_type
    use Variables
    
        implicit none
    
    contains
    
        subroutine Reaction_Rate_Single(Phi, XS_k, V, RR) 
            
            real(prec), intent(in) :: Phi(:), XS_k(:), V(:)
            real(prec), intent(out) :: RR(:)
                        
            RR = Phi*XS_k*V
            
        end subroutine Reaction_Rate_Single
        
        subroutine Reaction_Rates_Print(Phi, Vol, Opt, XS_gl, filename) ! filename='ReactionRates_Ref' or 'ReactionRates_Acc' 
            real(prec), intent(in) :: Phi(:), Vol(:)
            type(Options_Data), intent(in) :: Opt
            type(XS_Data), intent(in) :: XS_gl
            
            character(len=17), intent(in) :: filename
            
            integer :: i, iu, t, t_tot
            real(prec), allocatable :: XS(:,:), RR(:,:)
            
            ! Size check 
            if (size(XS_gl%Tot) /= size(Phi)) then
                error stop 'ERROR: Phi and XS size do not match'
            end if
            if (size(Vol) /= Opt%n_tot) then
                error stop 'ERROR: Vol must be size n_tot'
            end if
            
            allocate(XS(8+Opt%n_g, size(Phi)))
            allocate(RR(8+Opt%n_g, size(Phi)))
            
            ! Conversion to matrix form for easier handling 
            XS(1,:)=XS_gl%Tot(:)
            XS(2,:)=XS_gl%Tra(:)
            XS(3,:)=XS_gl%Absr(:)
            XS(4,:)=XS_gl%Rem(:)
            XS(5,:)=XS_gl%Fis(:)
            XS(6,:)=XS_gl%nuFis(:)
            XS(7,:)=XS_gl%kFis(:)
            XS(8,:)=XS_gl%Chi(:)
            do i=1, Opt%n_g 
                XS(8+i,:)=XS_gl%Scatt(:,i)
            end do
            
            call files%open_single_file(trim(output_dir)//trim(debug_mark)// &
            trim(filename)//trim(file_type), iu, "replace", "write")
            write(iu, '(A)') 'Tot                  Tra                       ' //&
            'Absr                      Rem                       Fis         ' //&
            '              NuFis                     KFis                    ' //&
            '  Chi                       Sca1                      Sca2'
            do t=1, Opt%n_g
                t_tot=(t-1)*Opt%n_tot
                do i = 1, size(Vol)
                    RR(:, t_tot+i) = XS(:, t_tot+i)*Phi(t_tot+i)*Vol(i)
                    write(iu,'(10000000ES21.13)') RR(:, t_tot+i)
                end do
            end do
            call files%close_single_file(iu)
            
        end subroutine Reaction_Rates_Print
        
        
        
end module Intermediates
