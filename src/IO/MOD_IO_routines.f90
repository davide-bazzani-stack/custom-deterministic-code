module IO_module
    
    ! Available to outside:
    ! ---> output_dir = directory for the output
    ! ---> debug_mark = file mark for debug   
    ! ---> output_mark = file mark for output
    ! ---> file_type = file format
    !           
    ! -----> trim(output_dir)//trim(debug_mark)//trim(filename)//trim(file_type)
    !
    !
    ! ---> output_files_t = derived type skeleton for file opening / closing
    ! ---> files = global variable for file opening / closing
    ! ---> write_file = interface for line-writing
    !      call write_file(ID, values, format_string) for int / float
    !      call write_file(ID, message) for char
    ! ---> format_builder = interface for building a format, custom or standard
    !      format_var = format_builder(values, 'custom_format')

    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    use precision_kinds, only: prec
    
        implicit none

    character(len = 7), parameter, public :: output_dir = "output/"
    character(len = 7), parameter, public :: debug_mark = "[Debug]"
    character(len = 8), parameter, public :: output_mark = "[Output]"
    character(len = 4), parameter, public :: file_type = ".out"
    !logical :: output__is_enabled = .true.
    !logical :: debug__is_enabled = .false.

    private

    type, public :: output_files_t
        integer :: log = -1
        integer :: convergence = -1
        integer :: phi = -1
        integer :: eigenfunctions = -1
        integer :: eigenvalues = -1
        integer :: performance = -1
        integer :: coord_centers = -1
        integer :: coord_vertices = -1
        integer :: power = -1
        integer :: universes = -1
        integer :: XS_debug = -1
        integer :: phi_debug = -1
        integer :: J_debug = -1
        integer :: D_til_debug = -1
        integer :: D_hat_debug = -1
        integer :: D_til_iter_debug = -1
        integer :: D_hat_iter_debug = -1
        integer :: mm_debug = -1
        integer :: Phi_Gl_Avg_debug = -1
        integer :: Phi_Gl_debug = -1
        integer :: Phi_Lo_debug = -1
        integer :: phi_iter_debug = -1
        integer :: iter_error_debug = -1
        contains
            procedure :: open_IO => open_file_main
            procedure :: close_IO => close_file_main
            procedure, nopass :: open_single_file => open_single_file_proc
            procedure, nopass :: close_single_file => close_file_proc
    end type

    type(output_files_t), public :: files

    interface write_file
        module procedure write_blank_line
        module procedure write_text_proc
        module procedure write_text_vector_proc
        module procedure write_real_scalar_proc
        module procedure write_real_vector_proc
        module procedure write_real_matrix_proc
        module procedure write_integer_scalar_proc
        module procedure write_integer_vector_proc
    end interface write_file

    interface format_builder
        module procedure format_builder_text
        module procedure format_builder_text_vector
        module procedure format_builder_real_scalar
        module procedure format_builder_real_vector
        module procedure format_builder_real_matrix
        module procedure format_builder_integer_scalar
        module procedure format_builder_integer_vector
    end interface format_builder

    public :: write_file, format_builder
    
    contains

        subroutine open_file_main(this, output__is_enabled, debug__is_enabled)
            class(output_files_t), intent(inout) :: this
            logical, intent(in) :: output__is_enabled, debug__is_enabled

            call file_ID_clear(this)
            if (output__is_enabled) call open_output_files_t(this)
            if (debug__is_enabled) call open_debug_files(this)

            call write_log_header(this, 'terminal')
            if (output__is_enabled) call write_log_header(this, 'file')

        end subroutine open_file_main

        subroutine file_ID_clear(this)
            class(output_files_t), intent(inout) :: this

            this%log = -1
            this%convergence = -1
            this%phi = -1
            this%eigenfunctions = -1
            this%eigenvalues = -1
            this%coord_centers = -1
            this%coord_vertices = -1
            this%power = -1
            this%universes = -1
            this%mm_debug = -1
            this%XS_debug = -1
            this%Phi_debug = -1
            this%J_debug = -1
            this%D_til_debug = -1
            this%D_hat_debug = -1
            this%D_hat_iter_debug = -1
            this%Phi_Gl_Avg_debug = -1
            this%Phi_Gl_debug = -1
            this%Phi_Lo_debug = -1
            this%performance = -1
            this%D_til_iter_debug = -1
            this%phi_iter_debug = -1
            this%iter_error_debug = -1

        end subroutine file_ID_clear

        subroutine open_output_files_t(this)
            type(output_files_t), intent(inout) :: this
            character(len = 256) :: filename

            filename = trim(output_dir) // trim(output_mark) // 'log' // trim(file_type)
            call this%open_single_file(filename, this%log, "replace", "write")
            
            filename = trim(output_dir) // trim(output_mark) // 'k' // trim(file_type)
            call this%open_single_file(filename, this%convergence, "replace", "write")
            
            filename = trim(output_dir) // trim(output_mark) // 'Phi' // trim(file_type)
            call this%open_single_file(filename, this%phi, "replace", "write")
            
            filename = trim(output_dir) // trim(output_mark) // 'Eigfcns' // trim(file_type)
            call this%open_single_file(filename, this%eigenfunctions, "replace", "write")
            
            filename = trim(output_dir) // trim(output_mark) // 'Eigvals' // trim(file_type)
            call this%open_single_file(filename, this%eigenvalues, "replace", "write")
            
            filename = trim(output_dir) // trim(output_mark) // 'Performance_Record' // trim(file_type)
            call this%open_single_file(filename, this%performance, "replace", "write")

            filename = trim(output_dir) // trim(output_mark) // 'Coordinates_InCntrs' // trim(file_type)
            call this%open_single_file(filename, this%coord_centers, "replace", "write")
            
            filename = trim(output_dir) // trim(output_mark) // 'Coordinates_Vert' // trim(file_type)
            call this%open_single_file(filename, this%coord_vertices, "replace", "write")

            filename = trim(output_dir) // trim(output_mark) // 'Power' // trim(file_type)
            call this%open_single_file(filename, this%power, "replace", "write") 
            
            filename = trim(output_dir) // trim(output_mark) // 'Universes' // trim(file_type)
            call this%open_single_file(filename, this%universes, "replace", "write")

        end subroutine open_output_files_t
        
        subroutine open_debug_files(this)
            type(output_files_t), intent(inout) :: this
            character(len = 256) :: filename

            filename = trim(output_dir) // trim(debug_mark) // 'XS' // trim(file_type)
            call this%open_single_file(filename, this%XS_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Flux' // trim(file_type)
            call this%open_single_file(filename, this%phi_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Currents' // trim(file_type)
            call this%open_single_file(filename, this%J_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'D_tilde' // trim(file_type)
            call this%open_single_file(filename, this%D_til_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'D_hat' // trim(file_type)
            call this%open_single_file(filename, this%D_hat_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'D_til_evo' // trim(file_type)
            call this%open_single_file(filename, this%D_til_iter_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'D_hat_evo' // trim(file_type)
            call this%open_single_file(filename, this%D_hat_iter_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Migration_Matrix' // trim(file_type)
            call this%open_single_file(filename, this%mm_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Flux_GL_avg' // trim(file_type)
            call this%open_single_file(filename, this%Phi_Gl_Avg_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Flux_GL' // trim(file_type)
            call this%open_single_file(filename, this%Phi_Gl_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Flux_LO' // trim(file_type)
            call this%open_single_file(filename, this%Phi_Lo_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Iteration_Flux' // trim(file_type)
            call this%open_single_file(filename, this%phi_iter_debug, "replace", "write")

            filename = trim(output_dir) // trim(debug_mark) // 'Iteration_error_L2' // trim(file_type)
            call this%open_single_file(filename, this%iter_error_debug, "replace", "write")

        end subroutine open_debug_files

        subroutine open_single_file_proc(filename, ID, file_status, file_action)
            character(len = *), intent(in) :: filename, file_status, file_action
            integer, intent(out) :: ID
            integer :: io_status
            character(len = 512) :: io_message

            open(newunit=ID, file=filename, status=file_status, action=file_action, &
                 iostat=io_status, iomsg=io_message)

            if (io_status /= 0) then
                write(error_unit, '(A)') trim(filename)
                write(error_unit, '(A)') trim(io_message)
                write(error_unit, *)
                write(error_unit, '(A)') "[SYSTEM] Error in file opening: " 
                write(error_unit, '(A)') "[SYSTEM]   --> MOD_IO_routines.f90/" // &
                                         "open_single_file_proc/" // filename 
                
                write(output_unit, *)
                write(output_unit, '(A)') "[SYSTEM] Error in file opening: " 
                write(output_unit, '(A)') "[SYSTEM]   --> MOD_IO_routines.f90/" // &
                                         "open_single_file_proc/" // filename 
                error stop
            end if

        end subroutine open_single_file_proc

        subroutine close_file_proc(ID)
            integer, intent(inout) :: ID
            logical :: is_open
            integer :: io_status
            character(len=512) :: io_message
                
            if (ID == -1) return
                
            inquire(unit=ID, opened=is_open)
                
            if (is_open) then

                close(ID, iostat=io_status, iomsg=io_message)

                if (io_status /= 0) then
                    write(output_unit, *)
                    write(output_unit, '(A)') trim(io_message)
                    write(error_unit, *)
                    write(error_unit, '(A)') trim(io_message)

                    return
                end if
            end if

            ID = -1
        end subroutine close_file_proc

        subroutine write_log_header(this, out_dest)
            class(output_files_t), intent(in) :: this
            character(len = *), intent(in) :: out_dest

            select case(trim(out_dest))
                case('terminal')
                    write(output_unit,'(A)') "----------------------------------------------------------------------------------"
                    write(output_unit,'(A)') "                  "
                    write(output_unit,'(A)') "                     =========  \\          //  | \        / |"
                    write(output_unit,'(A)') "                    ||           \\        //   ||\\      //||"
                    write(output_unit,'(A)') "                    ||            \\      //    || \\    // ||"
                    write(output_unit,'(A)') "                    |=====         \\    //     ||  \\  //  ||"
                    write(output_unit,'(A)') "                    ||              \\  //      ||   \\//   ||"
                    write(output_unit,'(A)') "                    ||               \\//       ||          ||"
                    write(output_unit,'(A)') "                    ||                \/        ||          ||"
                    write(output_unit,'(A)') "                              "
                    write(output_unit,'(A)') "----------------------------------------------------------------------------------"
                    write(output_unit,'(A)') " "
                    write(output_unit,'(A)') "Developed by D. Bazzani"
                    write(output_unit,'(A)') " "
                    write(output_unit,'(A)') " \  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /"
                    write(output_unit,'(A)') "  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/ "
                    write(output_unit,'(A)') " "
                    write(output_unit,'(A)') " "
                    write(output_unit,'(A)') " "
                case('file')
                    write(this%log,'(A)') "----------------------------------------------------------------------------------"
                    write(this%log,'(A)') "                  "
                    write(this%log,'(A)') "                     =========  \\          //  | \        / |"
                    write(this%log,'(A)') "                    ||           \\        //   ||\\      //||"
                    write(this%log,'(A)') "                    ||            \\      //    || \\    // ||"
                    write(this%log,'(A)') "                    |=====         \\    //     ||  \\  //  ||"
                    write(this%log,'(A)') "                    ||              \\  //      ||   \\//   ||"
                    write(this%log,'(A)') "                    ||               \\//       ||          ||"
                    write(this%log,'(A)') "                    ||                \/        ||          ||"
                    write(this%log,'(A)') "                              "
                    write(this%log,'(A)') "----------------------------------------------------------------------------------"
                    write(this%log,'(A)') " "
                    write(this%log,'(A)') "Developed by D. Bazzani"
                    write(this%log,'(A)') " "
                    write(this%log,'(A)') " \  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /"
                    write(this%log,'(A)') "  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/ "
                    write(this%log,'(A)') " "
                    write(this%log,'(A)') " "
                    write(this%log,'(A)') " "
                case default
                    write(output_unit,*) 
                    write(error_unit, *) 
                    write(output_unit,'(A)') "[SYSTEM] Wrong file type selection in" // &  
                    " MOD_IO_routines.f90/write_log_header" 
                    write(error_unit,'(A)') "[SYSTEM] Wrong file type selection in" // &  
                    " MOD_IO_routines.f90/write_log_header" 
            end select
        end subroutine write_log_header

        subroutine close_file_main(this, output__is_enabled, debug__is_enabled)
	        class(output_files_t), intent(inout) :: this
            logical, intent(in) :: output__is_enabled, debug__is_enabled
    
            if (output__is_enabled) then
                call this%close_single_file(this%log)
                call this%close_single_file(this%convergence)
                call this%close_single_file(this%phi)
                call this%close_single_file(this%eigenfunctions)
                call this%close_single_file(this%eigenvalues)
                call this%close_single_file(this%performance)
                call this%close_single_file(this%coord_centers)
                call this%close_single_file(this%coord_vertices)
                call this%close_single_file(this%power)
                call this%close_single_file(this%universes)
            end if
    
            if (debug__is_enabled) then
                call this%close_single_file(this%XS_debug)
                call this%close_single_file(this%phi_debug)
                call this%close_single_file(this%J_debug)
                call this%close_single_file(this%D_til_debug)
                call this%close_single_file(this%D_hat_debug)
                call this%close_single_file(this%D_til_iter_debug)
                call this%close_single_file(this%D_hat_iter_debug)
                call this%close_single_file(this%mm_debug)
                call this%close_single_file(this%Phi_Gl_Avg_debug)
                call this%close_single_file(this%Phi_Gl_debug)
                call this%close_single_file(this%Phi_Lo_debug)
                call this%close_single_file(this%phi_iter_debug)
                call this%close_single_file(this%iter_error_debug)
            end if
        
        end subroutine close_file_main
        
        
        


        subroutine write_blank_line(ID)
            integer, intent(in) :: ID

            write(ID, *) 
            
        end subroutine write_blank_line

        subroutine write_text_proc(ID, message, format_string)
            integer, intent(in) :: ID
            character(len=*), intent(in) :: message
            character(len=*), intent(in), optional :: format_string

            if (present(format_string)) then
                write(ID, format_string) trim(message)
            else
                write(ID, '(A)') trim(message)
            end if
        end subroutine write_text_proc

        subroutine write_text_vector_proc(ID, messages, format_string)
            integer, intent(in) :: ID
            character(len=*), intent(in) :: messages(:)
            character(len=*), intent(in), optional :: format_string
            integer :: i

            if (present(format_string)) then
                write(ID, format_string) (trim(messages(i)), i=1, size(messages))
            else
                write(ID, '(*(A,1X))') (trim(messages(i)), i=1, size(messages))
            end if

        end subroutine write_text_vector_proc
        
        subroutine write_integer_scalar_proc(ID, value, format_string)
            integer, intent(in) :: ID
            integer, intent(in) :: value
            character(len=*), intent(in), optional :: format_string

            if (present(format_string)) then
                write(ID, format_string) value
            else
                write(error_unit, '(A)') "[SYSTEM] Missing format variables"
                write(output_unit, '(A)') "[SYSTEM] Missing format variables"
                error stop
            end if

        end subroutine write_integer_scalar_proc

        subroutine write_integer_vector_proc(ID, values, format_string)
            integer, intent(in) :: ID
            integer, intent(in) :: values(:)
            character(len=*), intent(in), optional :: format_string

            if (present(format_string)) then
                write(ID, format_string) values
            else
                write(error_unit, '(A)') "[SYSTEM] Missing format variables"
                write(output_unit, '(A)') "[SYSTEM] Missing format variables"
                error stop
            end if

        end subroutine write_integer_vector_proc

        subroutine write_real_scalar_proc(ID, value, format_string)
            integer, intent(in) :: ID
            real(prec), intent(in) :: value
            character(len=*), intent(in), optional :: format_string

            if (present(format_string)) then
                write(ID, format_string) value
            else
                write(error_unit, '(A)') "[SYSTEM] Missing format variables"
                write(output_unit, '(A)') "[SYSTEM] Missing format variables"
                error stop
            end if

        end subroutine write_real_scalar_proc

        subroutine write_real_vector_proc(ID, values, format_string)
            integer, intent(in) :: ID
            real(prec), intent(in) :: values(:)
            character(len=*), intent(in), optional :: format_string

            if (present(format_string)) then
                write(ID, format_string) values
            else
                write(error_unit, '(A)') "[SYSTEM] Missing format variables"
                write(output_unit, '(A)') "[SYSTEM] Missing format variables"
                error stop
            end if

        end subroutine write_real_vector_proc

        subroutine write_real_matrix_proc(ID, values, format_string)
            integer, intent(in) :: ID
            real(prec), intent(in) :: values(:,:)
            character(len = *), intent(in), optional :: format_string
            
            integer :: i
            
            do i = 1, size(values, 1)
                if (present(format_string)) then
                    write(ID, format_string) values(i,:)
                else
                    write(error_unit, '(A)') "[SYSTEM] Missing format variables"
                    write(output_unit, '(A)') "[SYSTEM] Missing format variables"
                    error stop
                end if
            end do
        end subroutine write_real_matrix_proc
        


        

        function format_builder_text(value, edit_descriptor) result(format_string)
            character(len=*), intent(in) :: value
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = scalar_format('A', edit_descriptor)

        end function format_builder_text

        function format_builder_text_vector(values, edit_descriptor) result(format_string)
            character(len=*), intent(in) :: values(:)
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = repeated_format(size(values), 'A', edit_descriptor)

        end function format_builder_text_vector

        function format_builder_integer_scalar(value, edit_descriptor) result(format_string)
            integer, intent(in) :: value
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = scalar_format('I0', edit_descriptor)

        end function format_builder_integer_scalar

        function format_builder_integer_vector(values, edit_descriptor) result(format_string)
            integer, intent(in) :: values(:)
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = repeated_format(size(values), 'I0', edit_descriptor)

        end function format_builder_integer_vector

        function format_builder_real_scalar(value, edit_descriptor) result(format_string)
            real(prec), intent(in) :: value
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = scalar_format('ES21.13', edit_descriptor)

        end function format_builder_real_scalar

        function format_builder_real_vector(values, edit_descriptor) result(format_string)
            real(prec), intent(in) :: values(:)
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = repeated_format(size(values), 'ES21.13', edit_descriptor)

        end function format_builder_real_vector

        function format_builder_real_matrix(values, edit_descriptor) result(format_string)
            real(prec), intent(in) :: values(:,:)
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string

            format_string = repeated_format(size(values, 2), 'ES21.13', edit_descriptor)

        end function format_builder_real_matrix

        function scalar_format(default_descriptor, edit_descriptor) result(format_string)
            character(len=*), intent(in) :: default_descriptor
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string
            character(len=:), allocatable :: descriptor

            descriptor = default_descriptor
            if (present(edit_descriptor)) descriptor = trim(edit_descriptor)

            if (len_trim(descriptor) == 0) then
                write(error_unit, '(A)') '[SYSTEM] IO/MOD_IO_Routines.f90/scalar_format:'// &
                                       ' The edit descriptor cannot be empty'
                write(output_unit, '(A)')'[SYSTEM] IO/MOD_IO_Routines.f90/scalar_format:'// &
                                       ' The edit descriptor cannot be empty'
                error stop 
            end if

            format_string = '(' // trim(descriptor) // ')'

        end function scalar_format

        function repeated_format(repetition_count, default_descriptor, edit_descriptor) &
            result(format_string)
            integer, intent(in) :: repetition_count
            character(len=*), intent(in) :: default_descriptor
            character(len=*), intent(in), optional :: edit_descriptor
            character(len=:), allocatable :: format_string
            character(len=:), allocatable :: descriptor
            character(len=32) :: count_string

            if (repetition_count <= 0) then
                write(error_unit, '(A)') '[SYSTEM] IO/MOD_IO_Routines.f90/repeated_format:'// &
                                       ' Cannot build a format for an empty array dimension'
                write(output_unit, '(A)')'[SYSTEM] IO/MOD_IO_Routines.f90/repeated_format:'// &
                                       ' Cannot build a format for an empty array dimension'
                error stop 
            end if

            descriptor = default_descriptor
            if (present(edit_descriptor)) descriptor = trim(edit_descriptor)

            if (len_trim(descriptor) == 0) then                
                write(error_unit, '(A)') '[SYSTEM] IO/MOD_IO_Routines.f90/repeated_format:'// &
                                       ' The edit descriptor cannot be empty'
                write(output_unit, '(A)')'[SYSTEM] IO/MOD_IO_Routines.f90/repeated_format:'// &
                                       ' The edit descriptor cannot be empty'
                error stop 
            end if

            write(count_string, '(I0)') repetition_count
            format_string = '(' // trim(count_string) // '(' // trim(descriptor) // ',1X))'

        end function repeated_format


end module IO_module
