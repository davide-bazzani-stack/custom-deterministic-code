module Intermediates
    use iso_fortran_env, only: output_unit, error_unit
    use precision_kinds, only: prec
    use IO_module, only: output_dir, debug_mark, file_type, output_files_t, &
                         write_file, format_builder
    use Variables
    
        implicit none
    
    type :: RR_t
        character(len = 10), allocatable :: names(:)
        real(prec), allocatable :: RR(:, :) ! n_tot * n_g array for RR storage 
        contains
            procedure :: init => Reaction_Rates_initialization
            procedure :: clear => Reaction_Rates_clear
            procedure :: assign => Reaction_Rates_assignment
            procedure :: print => Reaction_Rates_file_printing
    end type 
            
    public :: RR_t, Reaction_Rates_Single, Reaction_Rates_reset

    contains

        subroutine Reaction_Rates_initialization(this, names, XS, Phi, V)
            class(RR_t), intent(inout) :: this
            character(len = *), intent(in) :: names(:)
            real(prec), intent(in) :: XS(:, :), Phi(:), V(:)

            call this%clear()

            allocate(this%names(size(names)))
            allocate(this%RR(size(Phi), size(names)), source=0.0_prec)
            
            call this%assign(names, XS, Phi, V)

        end subroutine Reaction_Rates_initialization

        subroutine Reaction_Rates_clear(this)
            class(RR_t), intent(inout) :: this

            if (allocated(this%names)) deallocate(this%names)
            if (allocated(this%RR)) deallocate(this%RR)

        end subroutine Reaction_Rates_clear

        subroutine Reaction_Rates_reset(to_zero)
            real(prec), intent(inout) :: to_zero(:)

            to_zero = 0.0_prec
        
        end subroutine Reaction_Rates_reset

        subroutine Reaction_Rates_assignment(this, names, XS, Phi, V)
            class(RR_t), intent(inout) :: this
            character(len = *), intent(in) :: names(:)
            real(prec), intent(in) :: XS(:,:), Phi(:), V(:)
            integer :: i, t, l_0, l_1, n_RR, n_tot, n_g
            
            if (.not. allocated(this%names)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Reaction names are not allocated'
            end if
            if (.not. allocated(this%RR)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Reaction-rate storage is not allocated'
            end if

            n_RR = size(this%RR, 2)
            n_tot = size(V)

            if (n_RR <= 0) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' At least one reaction is required'
            end if
            if (n_tot <= 0) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' The volume array cannot be empty'
            end if
            if (size(Phi) <= 0) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' The flux array cannot be empty'
            end if
            if (mod(size(Phi), n_tot) /= 0) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Flux size must be a multiple of the volume size'
            end if
            if (size(this%names) /= n_RR) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Reaction names and rate columns are inconsistent'
            end if
            if (size(names) /= n_RR) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Invalid number of reaction names'
            end if
            if (size(this%RR, 1) /= size(Phi)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Flux size does not match reaction-rate storage'
            end if
            if (size(XS, 1) /= size(Phi)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Cross-section and flux sizes do not match'
            end if
            if (size(XS, 2) /= n_RR) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_assignment '// &
                          ' Cross-section and reaction counts do not match'
            end if

            n_g = size(Phi) / n_tot

            this%names = names

            do i=1, n_RR
                do t = 1, n_g
                    l_0 = (t-1)*n_tot+1 
                    l_1 = t*n_tot
                    call Reaction_Rates_Single( Phi(l_0:l_1), &
                                                XS(l_0:l_1, i), &
                                                V, &
                                                this%RR(l_0:l_1, i))
                end do
            end do

        end subroutine Reaction_Rates_assignment

        subroutine Reaction_Rates_Single(Phi, XS_i, V, RR) 
            
            real(prec), intent(in) :: Phi(:), XS_i(:), V(:)
            real(prec), intent(out) :: RR(:)

            if (size(Phi) /= size(XS_i)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_Single '// &
                           ' Flux and cross-section slice sizes do not match'
            end if
            if (size(Phi) /= size(V)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_Single'// &
                           ' Flux and volume slice sizes do not match'
            end if
            if (size(Phi) /= size(RR)) then
                error stop '[SYSTEM] MOD_DebugChecks.f90/Reaction_Rates_Single'// &
                           ' Flux and reaction-rate slice sizes do not match'
            end if

            RR = Phi*XS_i*V
            
        end subroutine Reaction_Rates_Single
        
        subroutine Reaction_Rates_file_printing(this, filename)
            class(RR_t), intent(in) :: this
            character(len = *), intent(in) :: filename
            character(len = 256) :: format_string
            type(output_files_t) :: IO_obj
            integer :: ID          
        
            call IO_obj%open_single_file(trim(output_dir) // trim(debug_mark) // &
                                         trim(filename) // trim(file_type), ID, & 
                                         'replace', 'write')
            

            format_string = format_builder(this%names)
            call write_file(ID, this%names, trim(format_string))
            call write_file(ID)
            format_string = format_builder(this%RR)
            call write_file(ID, this%RR, trim(format_string))
        
            call IO_obj%close_single_file(ID)
            
        
        
        end subroutine Reaction_Rates_file_printing


        !subroutine Reaction_Rates_Print(Phi, Vol, Opt, XS_gl, filename) ! filename='ReactionRates_Ref' or 'ReactionRates_Acc' 
        !    real(prec), intent(in) :: Phi(:), Vol(:)
        !    type(Options_Data), intent(in) :: Opt
        !    type(XS_Data), intent(in) :: XS_gl
        !    
        !    character(len=17), intent(in) :: filename
        !    
        !    integer :: i, iu, t, t_tot
        !    real(prec), allocatable :: XS(:,:), RR(:,:)
        !    
        !    ! Size check 
        !    if (size(XS_gl%Tot) /= size(Phi)) then
        !        error stop 'ERROR: Phi and XS size do not match'
        !    end if
        !    if (size(Vol) /= Opt%n_tot) then
        !        error stop 'ERROR: Vol must be size n_tot'
        !    end if
        !    
        !    allocate(XS(8+Opt%n_g, size(Phi)))
        !    allocate(RR(8+Opt%n_g, size(Phi)))
        !    
        !    ! Conversion to matrix form for easier handling 
        !    XS(1,:)=XS_gl%Tot(:)
        !    XS(2,:)=XS_gl%Tra(:)
        !    XS(3,:)=XS_gl%Absr(:)
        !    XS(4,:)=XS_gl%Rem(:)
        !    XS(5,:)=XS_gl%Fis(:)
        !    XS(6,:)=XS_gl%nuFis(:)
        !    XS(7,:)=XS_gl%kFis(:)
        !    XS(8,:)=XS_gl%Chi(:)
        !    do i=1, Opt%n_g 
        !        XS(8+i,:)=XS_gl%Scatt(:,i)
        !    end do
        !    
        !    call files%open_single_file(trim(output_dir)//trim(debug_mark)// &
        !    trim(filename)//trim(file_type), iu, "replace", "write")
        !    write(iu, '(A)') 'Tot                  Tra                       ' //&
        !    'Absr                      Rem                       Fis         ' //&
        !    '              NuFis                     KFis                    ' //&
        !    '  Chi                       Sca1                      Sca2'
        !    do t=1, Opt%n_g
        !        t_tot=(t-1)*Opt%n_tot
        !        do i = 1, size(Vol)
        !            RR(:, t_tot+i) = XS(:, t_tot+i)*Phi(t_tot+i)*Vol(i)
        !            write(iu,'(10000000ES21.13)') RR(:, t_tot+i)
        !        end do
        !    end do
        !    call files%close_single_file(iu)
        !    
        !end subroutine Reaction_Rates_Print
        
        
        
end module Intermediates
