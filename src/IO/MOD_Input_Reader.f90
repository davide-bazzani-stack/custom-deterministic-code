module Input_Reader
    use iso_fortran_env, only: output_unit, error_unit
    use precision_kinds, only: prec
    use constants, only: pi
    use Variables
    use IO_module, only : files
    use Service_Fcns

        implicit none
        
    contains
    
    
        subroutine Options_Card(flag_accel, Opt_lo, Opt_gl, dl_input, dx_input) 
            integer, intent(out) :: flag_accel   ! Flags
                
            type(Options_Data), intent(out) :: Opt_lo, Opt_gl
            type(Int_input_check) :: flag_ms, n_x, n_y, n_g, it_sol_max, it_out_max, it_in_max
            type(Flo_input_check) :: tol_solv, tol0, tol1, tol2, w_SOR, dxdl
            type(Char_input_check) :: mesh_strc
            
            real(prec), intent(out) :: dl_input, dx_input
        
            ! Service variables
            integer :: l, ios=0, iter
            character(len=256) :: file_line
            character(len=10)  :: mesh_strct, acc_type, SOR
            
            
            ! Options file reading
            open(unit=120, file='input/Options_Input.inp', status='old', action='read')
        
        
                
            ! Data extraction 
            rewind(120)   ! Restarting in case of weird file disposition
            l=0
            do while (ios==0)
                read(120,'(A)', iostat=ios) file_line
                if (ios>0) then
                    write(output_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING'
                    write(error_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING'
                    exit
                end if
            
                if (trim(file_line) == '' .or. file_line(1:1) == '#') cycle ! skip the comments and the blank lines
            
                ! Data core association 
                if (index(file_line, '% Local') > 0) then 
                    
                    ! Sentinel values 
                    Opt_lo%linear_solver='BiCGSTAB' ! Temporary hardcoded selection
                    Opt_lo%preconditioner='Identity' ! Temporary hardcoded selection
                    Opt_lo%eigen_spectrum_analysis=.false. ! Temporary hardcoded selection
                    mesh_strc%variable=''
                    flag_ms%value=-1
                    n_x%value=0
                    n_y%value=0
                    n_g%value=0
                    it_sol_max%value=0
                    it_out_max%value=0
                    it_in_max%value=0
                    w_SOR%value=0.e0_prec
                    tol_solv%value=0.e0_prec
                    tol0%value=0.e0_prec
                    tol1%value=0.e0_prec
                    tol2%value=0.e0_prec
                    mesh_strc%is_set=.FALSE.
                    flag_ms%is_set=.FALSE.
                    n_x%is_set=.FALSE.
                    n_y%is_set=.FALSE.
                    n_g%is_set=.FALSE.
                    it_sol_max%is_set=.FALSE.
                    it_out_max%is_set=.FALSE.
                    it_in_max%is_set=.FALSE.
                    w_SOR%is_set=.FALSE.
                    tol_solv%is_set=.FALSE.
                    tol0%is_set=.FALSE.
                    tol1%is_set=.FALSE.
                    tol2%is_set=.FALSE.
                    iter=0
                    
                    FVM_Reader: do 
                        iter=iter+1
                        read(120,'(A)') file_line
                        
                        ! Check the assignment
                        if (mesh_strc%is_set .AND. n_x%is_set .AND. n_y%is_set .AND. n_g%is_set .AND. it_sol_max%is_set .AND. it_out_max%is_set .AND. it_in_max%is_set .AND. w_SOR%is_set .AND. tol_solv%is_set .AND. tol0%is_set .AND. tol1%is_set .AND. tol2%is_set) then
                            exit FVM_Reader
                        end if
                        
                        ! Skip comments
                        if (trim(file_line) == '' .or. file_line(1:1) == '#') cycle ! skip the comments and the blank lines

                        if (file_line(1:4) == 'mesh') then 
                            read(120,'(A)') mesh_strct
                            mesh_strc%variable=mesh_strct
                            mesh_strc%is_set=.TRUE.
                        
                        elseif (file_line(1:2) == 'nx') then 
                            read(120,*) Opt_lo%n_x
                            n_x%value=Opt_lo%n_x
                            n_x%is_set=.TRUE.
                        
                            Opt_lo%n_x=Opt_lo%n_x+2 ! Adding the ghost layer
                        
                        elseif (file_line(1:2) == 'ny') then 
                            read(120,*) Opt_lo%n_y
                            n_y%value=Opt_lo%n_y
                            n_y%is_set=.TRUE.
                        
                            Opt_lo%n_y=Opt_lo%n_y+2 ! Adding the ghost layer
                            if (trim(mesh_strct)=='Tri_Iso' .OR. trim(mesh_strct)=='Tri_Eq') Opt_lo%n_y=2*Opt_lo%n_x
                        
                        elseif (file_line(1:2) == 'ng') then 
                            read(120,*) Opt_lo%n_g
                            n_g%value=Opt_lo%n_g
                            n_g%is_set=.TRUE.
                        
                        elseif (file_line(1:8) == 'solv_tol') then 
                            read(120,*) Opt_lo%tol_solv 
                            tol_solv%value=Opt_lo%tol_solv
                            tol_solv%is_set=.TRUE.
                        
                        elseif (file_line(1:8) == 'in_P_tol') then 
                            read(120,*) Opt_lo%tol0
                            tol0%value=Opt_lo%tol0
                            tol0%is_set=.TRUE.
                        
                        elseif (file_line(1:9) == 'out_k_tol') then 
                            read(120,*) Opt_lo%tol1
                            tol1%value=Opt_lo%tol1
                            tol1%is_set=.TRUE.
                        
                        elseif (file_line(1:9) == 'out_S_tol')then 
                            read(120,*) Opt_lo%tol2
                            tol2%value=Opt_lo%tol2
                            tol2%is_set=.TRUE.
                        
                        elseif (file_line(1:13) =='solv_max_iter') then 
                            read(120,*) Opt_lo%it_sol_max
                            it_sol_max%value=Opt_lo%it_sol_max
                            it_sol_max%is_set=.TRUE.
                        
                        elseif (file_line(1:11) == 'in_max_iter') then 
                            read(120,*) Opt_lo%it_in_max
                            it_in_max%value=Opt_lo%it_in_max
                            it_in_max%is_set=.TRUE.
                        
                        elseif (file_line(1:12) == 'out_max_iter') then 
                            read(120,*) Opt_lo%it_out_max
                            it_out_max%value=Opt_lo%it_out_max
                            it_out_max%is_set=.TRUE.
                        
                        elseif (file_line(1:3) == 'SOR') then 
                            read(120,'(A)') SOR
                        
                            if (trim(SOR)=='Y') then 
                                Opt_lo%flag_SOR=1
                                read(120,*) Opt_lo%w_SOR
                                w_SOR%value=Opt_lo%w_SOR
                                w_SOR%is_set=.TRUE.
                            
                            elseif (trim(SOR)=='N') then
                                Opt_lo%flag_SOR=0
                                Opt_lo%w_SOR=1.e0_prec
                                read(120,'(A)') file_line !# skip the line
                                w_SOR%value=1.e0_prec
                                w_SOR%is_set=.TRUE.
                            end if
                            
                        
                        elseif (file_line(1:1) == '%') then 
                                write(output_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - LOCAL not closed'
                                write(error_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - LOCAL not closed'
                                close(120)
                                exit FVM_Reader
                        end if
                        
                    ! Check on the loop
                        if (iter>1e6) then 
                            write(output_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - A Variable is not assigned'
                            write(error_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - A Variable is not assigned'
                            close(120)
                            exit FVM_Reader
                        end if
                                       
                    end do FVM_Reader

                elseif (index(file_line, '% Acceleration') > 0) then 
                    read(120,'(A)') acc_type
                    if (acc_type=='None') exit 
                    ! Sentinel values 
                    mesh_strc%variable=''
                    flag_ms%value=-1
                    n_x%value=0
                    n_y%value=0
                    n_g%value=0
                    it_sol_max%value=0
                    it_out_max%value=0
                    it_in_max%value=0
                    w_SOR%value=0.e0_prec
                    tol_solv%value=0.e0_prec
                    tol0%value=0.e0_prec
                    tol1%value=0.e0_prec
                    tol2%value=0.e0_prec
                    mesh_strc%is_set=.FALSE.
                    flag_ms%is_set=.FALSE.
                    n_x%is_set=.FALSE.
                    n_y%is_set=.FALSE.
                    n_g%is_set=.FALSE.
                    it_sol_max%is_set=.FALSE.
                    it_out_max%is_set=.FALSE.
                    it_in_max%is_set=.FALSE.
                    w_SOR%is_set=.FALSE.
                    tol_solv%is_set=.FALSE.
                    tol0%is_set=.FALSE.
                    tol1%is_set=.FALSE.
                    tol2%is_set=.FALSE.
                    if (acc_type .NE.'HGCMFD') dxdl%is_set=.TRUE. 
                    iter=0
                    
                    Accel_Reader: do 
                        iter=iter+1
                        read(120,'(A)', iostat=ios) file_line
                        
                        if (dxdl%is_set .AND. n_x%is_set .AND. n_y%is_set .AND. n_g%is_set .AND. it_sol_max%is_set .AND. it_out_max%is_set .AND. it_in_max%is_set .AND. w_SOR%is_set .AND. tol_solv%is_set .AND. tol0%is_set .AND. tol1%is_set .AND. tol2%is_set) then
                            exit Accel_Reader
                        end if
                        
                        ! Skip comments
                        if (trim(file_line) == '' .or. file_line(1:1) == '#') cycle ! skip the comments and the blank lines
                        
                        if (file_line(1:2) == 'nx') then 
                            read(120,*) Opt_gl%n_x
                            n_x%value=Opt_gl%n_x
                            n_x%is_set=.TRUE.
                        
                        elseif (file_line(1:2) == 'ny') then 
                            read(120,*) Opt_gl%n_y
                            n_y%value=Opt_gl%n_y
                            n_y%is_set=.TRUE.
                        
                        elseif (file_line(1:2) == 'ng') then 
                            read(120,*) Opt_gl%n_g
                            n_g%value=Opt_gl%n_g
                            n_g%is_set=.TRUE.
                        
                        elseif (file_line(1:8) == 'solv_tol') then 
                            read(120,*) Opt_gl%tol_solv 
                            tol_solv%value=Opt_gl%tol_solv
                            tol_solv%is_set=.TRUE.
                        
                        elseif (file_line(1:8) == 'in_P_tol') then 
                            read(120,*) Opt_gl%tol0
                            tol0%value=Opt_gl%tol0
                            tol0%is_set=.TRUE.
                        
                        elseif (file_line(1:9) == 'out_k_tol') then 
                            read(120,*) Opt_gl%tol1
                            tol1%value=Opt_gl%tol1
                            tol1%is_set=.TRUE.
                        
                        elseif (file_line(1:9) == 'out_S_tol')then 
                            read(120,*) Opt_gl%tol2
                            tol2%value=Opt_gl%tol2
                            tol2%is_set=.TRUE.
                        
                        elseif (file_line(1:13) =='solv_max_iter') then 
                            read(120,*) Opt_gl%it_sol_max
                            it_sol_max%value=Opt_gl%it_sol_max
                            it_sol_max%is_set=.TRUE.
                        
                        elseif (file_line(1:11) == 'in_max_iter') then 
                            read(120,*) Opt_gl%it_in_max
                            it_in_max%value=Opt_gl%it_in_max
                            it_in_max%is_set=.TRUE.
                        
                        elseif (file_line(1:12) == 'out_max_iter') then 
                            read(120,*) Opt_gl%it_out_max
                            it_out_max%value=Opt_gl%it_out_max
                            it_out_max%is_set=.TRUE.
                        
                        elseif (file_line(1:3) == 'SOR') then 
                            read(120,'(A)') SOR
                        
                            if (trim(SOR)=='Y') then 
                                Opt_gl%flag_SOR=1
                                read(120,*) Opt_gl%w_SOR
                                w_SOR%value=Opt_gl%w_SOR
                                w_SOR%is_set=.TRUE.
                            
                            elseif (trim(SOR)=='N') then
                                Opt_gl%flag_SOR=0
                                Opt_gl%w_SOR=1.e0_prec
                                read(120,'(A)') file_line !# skip the line
                                w_SOR%value=1.e0_prec
                                w_SOR%is_set=.TRUE.
                            end if
                                                                        
                        elseif (file_line(1:4) == 'dxdl') then 
                            read(120,*) dl_input
                            read(120,*) dx_input
                            
                            
                            dxdl%is_set=.TRUE.
                        
                        elseif (file_line(1:1) == '%') then 
                            write(output_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - LOCAL not closed'
                            write(error_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - LOCAL not closed'
                            exit Accel_Reader
                        end if
                    
                        ! Check on the loop
                        if (iter>1e6) then 
                            write(output_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - A Variable is not assigned'
                            write(error_unit,'(A)') 'ERROR IN THE OPTIONS CARD READING - A Variable is not assigned'
                            exit Accel_Reader
                        end if
                                       
                    end do Accel_Reader
                
                end if
            end do 
        
            close(120)
        
            select case (mesh_strct) 
                case ('Quad')
                    Opt_lo%flag_ms=0
                case ('Tri_Iso') 
                    Opt_lo%flag_ms=1
                case ('Tri_Eq') 
                    Opt_lo%flag_ms=2
                case default
                    mesh_strct='ERROR!'
                    write(output_unit,'(A)') " ERROR! MESH SELECTION!"
                    write(error_unit,'(A)') " ERROR! MESH SELECTION!"
            end select
        
            select case (acc_type) 
                case ('None')   ! No acceleration
                    flag_accel=0
                case ('CMFD')   ! CMFD formulation
                    flag_accel=1
                case ('pCMFD')  ! pCMFD formulation
                    flag_accel=2
                case ('lpCMFD') ! lpCMFD formulation
                    flag_accel=3
                case ('HGCMFD') ! HGCMFD formulation
                    flag_accel=4
                case ('pHGCMFD') ! HGCMFD formulation
                    flag_accel=5
                case default
                    flag_accel=-1
                    write(output_unit,'(A)') " ERROR! WRONG ACCELERATION SELECTION!"
                    write(error_unit,'(A)') " ERROR! WRONG ACCELERATION SELECTION!"
            end select
        
        
        end subroutine  Options_Card
        
        subroutine Material_Card(n_materials, n_groups, Input_Materials) 
            integer, intent(in) :: n_groups
            integer, intent(out) :: n_materials
            type(Material), allocatable, intent(out) :: Input_Materials(:)
            
            integer :: i, j, l, g, m, ios=0
            character(len=256) :: file_line, comment
            
            ! Input type structure: 
            !
            ! # Material data input
            !   n_mat = 2
            !   
            !   % mat 1
            !   XS_TO =
            !       1.2 1.4
            !   XS_TR =
            !       1.2 1.4
            !   XS_AB =
            !       1.2 1.4
            !   XS_nF =
            !       1.2 1.4
            !   XS_FS =
            !       1.2 1.4
            !   Chi = 
            !       1.2 1.4
            !   XS_kF =
            !       1.2 1.4
            !   XS_RM = 
            !       1.2 1.4
            !   XS_D =
            !       1.2 1.4
            !   XS_SC =
            !       1.2 1.4 # From g=1
            !       1.2 1.4 # From g=2 ...
            open(unit=100, file='input/Material_Input.inp', status='old', action='read')
            
            
            write(output_unit,*) 
            write(files%log,*) 
            write(output_unit,'(A)') 'Reading the material input...'
            write(files%log,'(A)') 'Reading the material input...'
            
            ! Groups and materials extraction 
            n_materials=0
            do 
                read(100,'(A)', iostat=ios) file_line
                if (ios>0) then
                    write(output_unit,'(A)') 'ERROR IN THE MATERIAL CARD READING'
                    write(error_unit,'(A)') 'ERROR IN THE MATERIAL CARD READING'
                    exit
                end if
                
                if (trim(file_line) == '' .or. file_line(1:1) == '#') cycle ! skip the comments and the blank lines
                
                if (index(file_line, 'n_mat') > 0) then ! Extracting n_materials
                    read(file_line(index(file_line, '=')+1:), '(I10,A)') n_materials, comment
                end if
                if (n_materials>0) exit
            end do
            
            ! Allocation 
            allocate(Input_Materials(n_materials)) 
            do m=1, n_materials
                allocate(Input_Materials(m)%XS_TO(n_groups))
                allocate(Input_Materials(m)%XS_TR(n_groups))
                allocate(Input_Materials(m)%XS_AB(n_groups))
                allocate(Input_Materials(m)%XS_nF(n_groups)) 
                allocate(Input_Materials(m)%XS_FS(n_groups))
                allocate(Input_Materials(m)%XS_SC(n_groups,n_groups))
                allocate(Input_Materials(m)%Chi(n_groups))
                allocate(Input_Materials(m)%XS_kF(n_groups))
                allocate(Input_Materials(m)%XS_RM(n_groups))
                allocate(Input_Materials(m)%XS_D(n_groups))
                
                Input_Materials(m)%XS_TO=0.e0_prec
                Input_Materials(m)%XS_TR=0.e0_prec
                Input_Materials(m)%XS_AB=0.e0_prec
                Input_Materials(m)%XS_nF=0.e0_prec
                Input_Materials(m)%XS_FS=0.e0_prec
                Input_Materials(m)%XS_SC=0.e0_prec
                Input_Materials(m)%Chi  =0.e0_prec
                Input_Materials(m)%XS_kF=0.e0_prec
                Input_Materials(m)%XS_RM=0.e0_prec
                Input_Materials(m)%XS_D =0.e0_prec
            end do
            
            ! Data extraction 
            rewind(100)   ! Restarting in case of weird file disposition
            l=0
            do while (ios==0)
                read(100,'(A)', iostat=ios) file_line
                if (ios>0) then
                    write(output_unit,'(A)') 'ERROR IN THE MATERIAL CARD READING'
                    write(error_unit,'(A)') 'ERROR IN THE MATERIAL CARD READING'
                    exit
                end if
                
                if (trim(file_line) == '' .or. file_line(1:1) == '#') cycle ! skip the comments and the blank lines
                
                ! Data core association
                if (index(file_line, '% mat') > 0) then
                    read(file_line(index(file_line, 'mat')+4:), *) l
                    Input_Materials(l)%ID=l
                    
                elseif (index(file_line, 'Name') > 0) then ! Extracting The material name
                    read(100, '(A)') Input_Materials(l)%name
                    Input_Materials(l)%name=trim(Input_Materials(l)%name)
                    
                elseif (index(file_line, 'XS_TO') > 0) then ! Extracting Total XS
                    read(100, *) Input_Materials(l)%XS_TO(:)
                    
                elseif (index(file_line, 'XS_TR') > 0) then ! Extracting Transport XS
                    read(100, *) Input_Materials(l)%XS_TR(:)
                    
                elseif (index(file_line, 'XS_AB') > 0) then ! Extracting Absorption XS
                    read(100, *) Input_Materials(l)%XS_AB(:)
                    
                elseif (index(file_line, 'XS_nF') > 0) then ! Extracting ni*(Fiss. XS)
                    read(100, *) Input_Materials(l)%XS_nF(:)
                    
                elseif (index(file_line, 'XS_FS') > 0) then ! Extracting Fiss. XS
                    read(100, *) Input_Materials(l)%XS_FS(:)
                    
                elseif (index(file_line, 'Chi') > 0) then ! Extracting Chi
                    read(100, *) Input_Materials(l)%Chi(:)
                    
                elseif (index(file_line, 'XS_kF') > 0) then ! Extracting kappa*(Fiss. XS) for power
                    read(100, *) Input_Materials(l)%XS_kF(:)
                    
                elseif (index(file_line, 'XS_RM') > 0) then ! Extracting Removal XS
                    read(100, *) Input_Materials(l)%XS_RM(:)
                    
                elseif (index(file_line, 'XS_D') > 0) then ! Extracting Diffusion Coeff.
                    read(100, *) Input_Materials(l)%XS_D(:)
                    
                elseif (index(file_line, 'XS_SC') > 0) then ! Extracting Scattering XS
                    do g=1,n_groups
                        read(100, *) Input_Materials(l)%XS_SC(g, :)
                    end do
                end if
            end do 
            
            ! Verification of the removal cross section presence 
            do l=1, n_materials
                if (sum(Input_Materials(l)%XS_RM(:))<1.e-12_prec) then
                    Input_Materials(l)%XS_RM=0.e0_prec
                    do g=1, n_groups
                        Input_Materials(l)%XS_RM(g)=Input_Materials(l)%XS_AB(g)+sum(Input_Materials(l)%XS_SC(g, :))-Input_Materials(l)%XS_SC(g, g)
                    end do
                end if
            end do
            
            ! Approximation of the Transport XS to the Total one in case of no input 
            do l=1, n_materials 
                if (sum(Input_Materials(l)%XS_TR(:))<1.e-12_prec) then
                    Input_Materials(l)%XS_TR(:)=Input_Materials(l)%XS_TO(:)
                end if
            end do
            
            ! Verification of the diffusion coefficient 
            do l=1, n_materials
                if (sum(Input_Materials(l)%XS_D(:))<1.e-12_prec) then
                    Input_Materials(l)%XS_D(:)=1.e0_prec/(3.e0_prec*Input_Materials(l)%XS_TR(:))
                end if
            end do
            
            close(100)
            write(output_unit,'(A)') '...Done'
            write(files%log,'(A)') '...Done'
            write(output_unit,*) 
            write(files%log,*) 
            
        end subroutine  Material_Card
    
        subroutine Geometry_Card(Opt_lo, Input_Materials, Elem_det, XS_lo, LO_Mesh, LO_Coef, dx_gh, dy_gh) 
            
            type(Material), allocatable, intent(in) :: Input_Materials(:)
            type(Options_Data), intent(inout) :: Opt_lo
            type(XS_Data), intent(out) :: XS_lo
            type(Figure), allocatable, intent(out) :: Elem_det(:)
            type(LO_geom), allocatable, intent(out) :: LO_Mesh(:)
            type(LO_coeff), allocatable, intent(out) :: LO_Coef(:)
            
            !integer, allocatable, intent(out) :: IDs_lo_Red(:)
            real(prec), intent(out) :: dx_gh, dy_gh
            
            integer :: i, j, l, n, u, t, t_tot, g, m, ll, n_tot, ios, n_elements, box_mat_ID, no_faces, iln
            real(prec) :: face_center(2)
            real(prec) :: xmin, xmax, base_width, side, height, apothem
            real(prec), allocatable :: gam_lo(:)
            character(len=256) :: file_line, Elem_type
            
            logical :: test
            logical, allocatable :: test_map(:)

            
            
            Opt_lo%n_tot=Opt_lo%n_x*Opt_lo%n_y
            allocate(LO_Mesh(Opt_lo%n_tot))
            allocate(LO_Coef(Opt_lo%n_tot*Opt_lo%n_g))
            
            ! BC dummy 
            select case(Opt_lo%flag_ms)
                case(0)
                    allocate(gam_lo(4))
                case(1)
                    allocate(gam_lo(6))
                case(2)
                    allocate(gam_lo(6))
            end select 
            gam_lo=-1.e0_prec
                        
            
            write(output_unit,*) 
            write(files%log,*) 
            write(output_unit,'(A)') 'Reading the geometry input...'
            write(files%log,'(A)') 'Reading the geometry input...'
            
            
            
            ! Elements data
            open(unit=110, file='input/Geometry_Input.inp', status='old', action='read')
            
            ! Automatic detection of the number of objects 
            n=0
            ios=0
            do while (ios==0)
                read(110,'(A)', iostat=ios) file_line
                if (ios>0) then
                    write(output_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Automatic detection routine'
                    write(error_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Automatic detection routine'
                    exit
                end if
                if (file_line(1:1) == '%') n=n+1
            end do
            n_elements=n
            
            ! Allocation and detection of the geometry details 
            allocate(Elem_det(n_elements)) 
            rewind(110)
            n=0
            Elem_det(:)%flag_BG=0
            ios=0
            do while (ios==0)
                read(110,'(A)', iostat=ios) file_line
                if (ios>0) then
                    write(output_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Element reading routine'
                    write(error_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Element reading routine'
                    exit
                end if
                
                if (trim(file_line)=='' .OR. file_line(1:1)=='#') then  ! skip the comments and the blank lines
                    cycle
                    
                elseif (file_line(1:4)=='[BC]') then 
                    select case(Opt_lo%flag_ms)
                        case(0) ! Square 
                            do i=1, 4
                                read(110,'(A)') file_line
                                select case (trim(file_line(1:2)))
                                    case('B') 
                                        read(file_line(6:),*) gam_lo(1)
                                    case('R') 
                                        read(file_line(6:),*) gam_lo(2)
                                    case('T') 
                                        read(file_line(6:),*) gam_lo(3)
                                    case('L') 
                                        read(file_line(6:),*) gam_lo(4)
                                end select
                            end do
                                
                        case(1) ! Tri Iso 
                            do i=1, 6
                                read(110,'(A)') file_line
                                select case (trim(file_line(1:2)))
                                    case('BL') 
                                        read(file_line(6:),*) gam_lo(1)
                                    case('B') 
                                        read(file_line(6:),*) gam_lo(2)
                                    case('BR') 
                                        read(file_line(6:),*) gam_lo(3)
                                    case('TR') 
                                        read(file_line(6:),*) gam_lo(4)
                                    case('T') 
                                        read(file_line(6:),*) gam_lo(5)
                                    case('TL') 
                                        read(file_line(6:),*) gam_lo(6)
                                end select
                            end do
                            
                        case(2) ! Tri Equi 
                            do i=1, 6
                                read(110,'(A)') file_line
                                select case (trim(file_line(1:2)))
                                    case('BL') 
                                        read(file_line(6:),*) gam_lo(1)
                                    case('B') 
                                        read(file_line(6:),*) gam_lo(2)
                                    case('BR') 
                                        read(file_line(6:),*) gam_lo(3)
                                    case('TR') 
                                        read(file_line(6:),*) gam_lo(4)
                                    case('T') 
                                        read(file_line(6:),*) gam_lo(5)
                                    case('TL') 
                                        read(file_line(6:),*) gam_lo(6)
                                end select
                            end do
                            
                        case default 
                            write(output_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Wrong mesh in the BC reading'
                            write(error_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Wrong mesh in the BC reading'
                    end select
                    if (any(gam_lo<0.e0_prec)) then
                        write(output_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Missing BC'
                        write(error_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Missing BC'
                        exit
                    end if
                    
                elseif (index(file_line, '%') > 0) then 
                    n=n+1
                    Elem_det(n)%ID=n
                    read(file_line(index(file_line, '%')+2:), '(A)') Elem_type
                    select case(trim(Elem_type)) 
                        
                        case ('Crfc') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,1))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                        
                        case ('Trng') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,3))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                        
                        case ('Rect') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,4))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            
                        case ('Pent') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,5))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            read(110, *) Elem_det(n)%vertices(:,5)
                        
                        case ('Hex') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,6))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            read(110, *) Elem_det(n)%vertices(:,5)
                            read(110, *) Elem_det(n)%vertices(:,6)
                        
                        case ('Ept') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,7))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            read(110, *) Elem_det(n)%vertices(:,5)
                            read(110, *) Elem_det(n)%vertices(:,6)
                            read(110, *) Elem_det(n)%vertices(:,7)
                            
                        case ('Oct') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,8))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            read(110, *) Elem_det(n)%vertices(:,5)
                            read(110, *) Elem_det(n)%vertices(:,6)
                            read(110, *) Elem_det(n)%vertices(:,7)
                            read(110, *) Elem_det(n)%vertices(:,8)
                        
                        case ('Enna') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,9))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            read(110, *) Elem_det(n)%vertices(:,5)
                            read(110, *) Elem_det(n)%vertices(:,6)
                            read(110, *) Elem_det(n)%vertices(:,7)
                            read(110, *) Elem_det(n)%vertices(:,8)
                            read(110, *) Elem_det(n)%vertices(:,9)
                            
                        case ('Deca') 
                            allocate(Elem_det(n)%centroid(2))
                            allocate(Elem_det(n)%vertices(2,10))
                            read(110, *) Elem_det(n)%centroid
                            read(110, *) Elem_det(n)%vertices(:,1)
                            read(110, *) Elem_det(n)%vertices(:,2)
                            read(110, *) Elem_det(n)%vertices(:,3)
                            read(110, *) Elem_det(n)%vertices(:,4)
                            read(110, *) Elem_det(n)%vertices(:,5)
                            read(110, *) Elem_det(n)%vertices(:,6)
                            read(110, *) Elem_det(n)%vertices(:,7)
                            read(110, *) Elem_det(n)%vertices(:,8)
                            read(110, *) Elem_det(n)%vertices(:,9)
                            read(110, *) Elem_det(n)%vertices(:,10)
                        
                        ! Orientation distinction 
                        !case ('HexX') 
                        !    allocate(Elem_det(n)%centroid(2))
                        !    allocate(Elem_det(n)%vertices(2,6))
                        !    read(110, *) Elem_det(n)%centroid
                        !    read(110, *) Elem_det(n)%vertices(:,1)
                        !    read(110, *) Elem_det(n)%vertices(:,2)
                        !    read(110, *) Elem_det(n)%vertices(:,3)
                        !    read(110, *) Elem_det(n)%vertices(:,4)
                        !    read(110, *) Elem_det(n)%vertices(:,5)
                        !    read(110, *) Elem_det(n)%vertices(:,6)
                        ! 
                        !case ('HexY') 
                        !    allocate(Elem_det(n)%centroid(2))
                        !    allocate(Elem_det(n)%vertices(2,6))
                        !    read(110, *) Elem_det(n)%centroid
                        !    read(110, *) Elem_det(n)%vertices(:,1)
                        !    read(110, *) Elem_det(n)%vertices(:,2)
                        !    read(110, *) Elem_det(n)%vertices(:,3)
                        !    read(110, *) Elem_det(n)%vertices(:,4)
                        !    read(110, *) Elem_det(n)%vertices(:,5)
                        !    read(110, *) Elem_det(n)%vertices(:,6)
                        
                        case default 
                            write(output_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Element selection'
                            write(error_unit,'(A)') 'ERROR IN THE GEOMETRY CARD READING - Element selection'
                            exit
                        end select
                        
                elseif (index(file_line, '@') > 0) then 
                    read(file_line(index(file_line, '@')+2:), '(A)') Elem_det(n)%material
                    Elem_det(n)%material=trim(Elem_det(n)%material)
                    do l=1, size(Input_Materials)
                        if (trim(Elem_det(n)%material)==trim(Input_Materials(l)%name)) then
                            Elem_det(n)%Mat_ID=Input_Materials(l)%ID
                        end if
                    end do
                    
                elseif (index(file_line, '!box') > 0) then 
                    Elem_det(n)%flag_BG=1
                    box_mat_ID=n
                end if
            end do
            
            no_faces=size(Elem_det(box_mat_ID)%vertices,2)
            xmax=maxval(Elem_det(box_mat_ID)%vertices(1,:))
            xmin=minval(Elem_det(box_mat_ID)%vertices(1,:))
            side=xmax-xmin
            apothem=side/(2.e0_prec*tan(2.e0_prec*pi/no_faces))
            
            base_width=side
            height=apothem*2.e0_prec
                
            write(output_unit,'(A)') '...Done'
            write(files%log,'(A)') '...Done'
            write(output_unit,*) 
            write(files%log,*) 
            
            write(output_unit,*) 
            write(files%log,*) 
            write(output_unit,'(A)') 'Generating the mesh...'
            write(files%log,'(A)') 'Generating the mesh...'
            
            
            ! Meshing 
            
            ! Selection of the number of faces 
            select case (Opt_lo%flag_ms)
                case (0) ! Squared-based matrix
                    no_faces=4
                    
                case (1) ! TriIso-based matrix
                    no_faces=3
                    
                case (2) ! TriEq-based matrix
                    no_faces=3
                    
            end select
            
            ! Allocate and initialize the local variables 
            do l=1, Opt_lo%n_tot
                LO_Mesh(l)%ID=l
                allocate(LO_Mesh(l)%Sides_ID(no_faces))
                allocate(LO_Mesh(l)%Sides_Neigh(no_faces))
                allocate(LO_Mesh(l)%Cent(2))
                allocate(LO_Mesh(l)%Vert(2,no_faces))
                allocate(LO_Mesh(l)%Neigh_ID(no_faces))
                allocate(LO_Mesh(l)%BC(no_faces)) ! [flag]
                allocate(LO_Mesh(l)%J_flag(no_faces))
                allocate(LO_Mesh(l)%dl_lo(no_faces))
                allocate(LO_Mesh(l)%A_lo(no_faces))
                
                LO_Mesh(l)%Sides_ID=0
                LO_Mesh(l)%Sides_Neigh=0
                LO_Mesh(l)%Cent=0.e0_prec
                LO_Mesh(l)%Vert=0.e0_prec
                LO_Mesh(l)%Neigh_ID=-1
                LO_Mesh(l)%BC=-1
                LO_Mesh(l)%J_flag=-1
                LO_Mesh(l)%dx=0.e0_prec
                LO_Mesh(l)%dy=0.e0_prec
                LO_Mesh(l)%dl_lo=0.e0_prec
                LO_Mesh(l)%A_lo=0.e0_prec
                LO_Mesh(l)%V_lo=0.e0_prec
                LO_Mesh(l)%univ=0
                
                do t=1, Opt_lo%n_g
                    t_tot=(t-1)*Opt_lo%n_tot
                    
                    allocate(LO_Coef(t_tot+l)%Gam(no_faces))
                    allocate(LO_Coef(t_tot+l)%J_Part(no_faces,2))
                    allocate(LO_Coef(t_tot+l)%J_Net(no_faces))
                    allocate(LO_Coef(t_tot+l)%D_til(no_faces))
                    
                    LO_Coef(t_tot+l)%Gam=-1.e0_prec
                    LO_Coef(t_tot+l)%J_Part=0.e0_prec
                    LO_Coef(t_tot+l)%J_Net=0.e0_prec
                    LO_Coef(t_tot+l)%D_til=0.e0_prec
                end do
            end do
            
            ! Mesher 
            select case (Opt_lo%flag_ms)   ! Allocation dependant on the mesh type 
                case (0) ! Squared-based matrix  
                    no_faces=4
                    ! To the mesher two additional lines and columns are given as ghosts
                    call  Quad_Mesher(box_mat_ID, Elem_det, Opt_lo, LO_Mesh, dx_gh, dy_gh) 
                    ! Removing the influence of the ghost layer
                    do l=1, size(Elem_det)
                        Elem_det(l)%centroid(1)=Elem_det(l)%centroid(1)+dx_gh
                        Elem_det(l)%centroid(2)=Elem_det(l)%centroid(2)+dy_gh
                        Elem_det(l)%vertices(1,:)=Elem_det(l)%vertices(1,:)+dx_gh
                        Elem_det(l)%vertices(2,:)=Elem_det(l)%vertices(2,:)+dy_gh
                    end do
                    
                case (1) ! Isoscele triangle-based matrix 
                    no_faces=3
                    ! To the mesher two additional lines and columns are given as ghosts
                    call TriIso_Mesher(box_mat_ID, Elem_det, Opt_lo, LO_Mesh, dx_gh, dy_gh)
                    ! Removing the influence of the ghost layer
                    do l=1, size(Elem_det) 
                        Elem_det(l)%centroid(1)=Elem_det(l)%centroid(1)+dx_gh
                        Elem_det(l)%centroid(2)=Elem_det(l)%centroid(2)+dy_gh!*2
                        Elem_det(l)%vertices(1,:)=Elem_det(l)%vertices(1,:)+dx_gh
                        Elem_det(l)%vertices(2,:)=Elem_det(l)%vertices(2,:)+dy_gh!*2
                    end do
                    
                !case (2) ! Equilateral triangle-based matrix 
                !    no_faces=3
                !    ! To the mesher two additional lines and columns are given as ghosts
                !    call TriEq_Mesher(n_x, n_y, dx, dy, n_tot, IDs_lo, Mesh_Ornt, base_width, height, Centr_lo, Vert_lo, dl, A_lo, V_lo, IDs_neigh_lo)
                !    ! Removing the influence of the ghost layer
                !    do l=1, size(Elem_det)
                !        Elem_det(l)%centroid(1)=Elem_det(l)%centroid(1)+dx/2.e0_prec
                !        Elem_det(l)%centroid(2)=Elem_det(l)%centroid(2)+dy
                !        Elem_det(l)%vertices(1,:)=Elem_det(l)%vertices(1,:)+dx/2.e0_prec
                !        Elem_det(l)%vertices(2,:)=Elem_det(l)%vertices(2,:)+dy
                !    end do
                    
                case default 
                    write(output_unit,'(A)') " ERROR IN THE MESH GENERATION - Wrong mesh selection!"
                    write(error_unit,'(A)') " ERROR IN THE MESH GENERATION - Wrong mesh selection!"
            end select
            
            write(output_unit,'(A)') '...Done'
            write(files%log,'(A)') '...Done'
            write(output_unit,*) 
            write(files%log,*) 
            
            write(output_unit,*) 
            write(files%log,*) 
            write(output_unit,'(A)') 'Generating the mesh...'
            write(files%log,'(A)') 'Generating the mesh...'
            
            
            ! Mesh-Problem merging 
            allocate(test_map(Opt_lo%n_tot))   ! Map of conditions
            
            ! Logic map for box 
            do l=1,Opt_lo%n_tot
                call Point_in_Polygon_Test_2D(Elem_det(box_mat_ID)%vertices, LO_Mesh(l)%cent(:), test_map(l))
                if (test_map(l)) then
                    LO_Mesh(l)%univ=1
                    LO_Mesh(l)%Mat_ID=Elem_det(box_mat_ID)%Mat_ID
                else
                    LO_Mesh(l)%univ=-1
                end if
                
                do n=1, size(Elem_det)
                    if (n .NE. box_mat_ID) then
                        call Point_in_Polygon_Test_2D(Elem_det(n)%vertices, LO_Mesh(l)%cent(:), test)
                        if (test) then
                            LO_Mesh(l)%univ=Elem_det(n)%ID
                            LO_Mesh(l)%Mat_ID=Elem_det(n)%Mat_ID
                            exit
                        end if
                    end if
                end do
            end do
            
            ! Assignment of the BC / in-mesh mask layer 
            select case(Opt_lo%flag_ms)
                case(0) ! Square mesh (n_x * n_x) 
                    do j=2, Opt_lo%n_y-1
                        do i=2, Opt_lo%n_x-1
                            l=(j-1)*Opt_lo%n_x+i
                            if (test_map(l)) then
                                LO_Mesh(l)%BC(:)=0 ! Off by default
                                do u=1, size(LO_Mesh(l)%Neigh_ID)
                                    if (.not. test_map(LO_Mesh(l)%Neigh_ID(u))) then    ! Neighbour test
                                        LO_Mesh(l)%BC(u)=1
                                    end if
                                end do
                            end if
                            
                            do u=1, no_faces
                                if (LO_Mesh(l)%BC(u)==1) then
                                    do t=1, Opt_lo%n_g
                                        t_tot=(t-1)*Opt_lo%n_tot
                                        LO_Coef(t_tot+l)%Gam(u)=gam_lo(u)
                                    end do
                                end if
                            end do
                        end do
                    end do
                    
                case(1)  ! Iso-triangular mesh (n_x * 2*n_x) 
                
                    do j=3, Opt_lo%n_y-1
                        do i=2, Opt_lo%n_x-1
                            l=(j-1)*Opt_lo%n_x+i
                            if (test_map(l)) then
                                LO_Mesh(l)%BC(:)=0 ! Off by default
                                do u=1, size(LO_Mesh(l)%Neigh_ID)
                                    if (.not. test_map(LO_Mesh(l)%Neigh_ID(u))) then    ! Neighbour test
                                        LO_Mesh(l)%BC(u)=1
                                    end if
                                end do
                            end if
                            
                            do u=1, no_faces
                                if (LO_Mesh(l)%BC(u)==1) then
                                    do t=1, Opt_lo%n_g
                                        t_tot=(t-1)*Opt_lo%n_tot
                                        LO_Coef(t_tot+l)%Gam(u)=gam_lo(u)
                                    end do
                                end if
                            end do
                        end do
                    end do
                    
                !case(2) ! TO VERIFY 
                !
                !    do j=3, Opt_lo%n_y-2
                !        do i=2, Opt_lo%n_x-1
                !            l=(j-1)*Opt_lo%n_x+i
                !            if (test_map(l)) then
                !                LO_Mesh(l)%BC(:)=0 ! Off by default
                !                do u=1, size(LO_Mesh(l)%Neigh_ID)
                !                    if (test_map(LO_Mesh(l)%Neigh_ID(u))==.FALSE.) then    ! Neighbour test
                !                        LO_Mesh(l)%BC(u)=1
                !                    end if
                !                end do
                !            end if
                !        end do
                !    end do
                
                case default
                    write(output_unit,'(A)') " ERROR IN THE BC DETECTION - Wrong mesh selection!"
                    write(error_unit,'(A)') " ERROR IN THE BC DETECTION - Wrong mesh selection!"
                    
                end select     
            
            ! Location of the current calculation gamma 
            do l=1, Opt_lo%n_tot
                i=mod(l-1,Opt_lo%n_x)+1
                j=int((l-i)/Opt_lo%n_x)+1
                
                LO_Mesh(l)%J_flag=0
                
                no_faces=size(LO_Mesh(l)%Neigh_ID)
                
                !if (flag_accel==4) then
                    do u=1, no_faces
                        iln=LO_Mesh(l)%Neigh_ID(u)
                        if (iln<1) cycle
                        if ((LO_Mesh(iln)%univ .NE. LO_Mesh(l)%univ) .AND. LO_Mesh(l)%univ>1) then
                          LO_Mesh(l)%J_flag(u)=1
                        end if
                    end do
                !elseif (flag_accel==1 .OR. flag_accel==2 .OR. flag_accel==3)
                !    do u=1, no_faces
                !        iln=LO_Mesh(l)%Neigh_ID(u)
                !        if (iln<1) cycle
                !        if ((LO_Mesh(iln)%univ .NE. LO_Mesh(l)%univ)) then
                !          LO_Mesh(l)%J_flag(u)=1
                !        end if
                !    end do
                !end if
                
            end do
       
            ! Total to Reduced map 
            !allocate(IDs_lo_Red(Opt_lo%n_tot))
            !IDs_lo_Red=-1
            n=0
            do l=1, Opt_lo%n_tot
                LO_Mesh(l)%ID_Red=-1
                if (test_map(l)) then
                    n=n+1
                    !IDs_lo_Red(l)=n
                    LO_Mesh(l)%ID_Red=n
                end if
            end do
            Opt_lo%n_red=n
        
            ! Allocation and initialization 
            allocate(XS_lo%Tot(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Tra(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Absr(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%nuFis(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Fis(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Chi(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%kFis(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Rem(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Dif(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y))
            allocate(XS_lo%Scatt(Opt_lo%n_g*Opt_lo%n_x*Opt_lo%n_y, Opt_lo%n_g))
            
               XS_lo%Tot=0.e0_prec
               XS_lo%Tra=0.e0_prec
              XS_lo%Absr=0.e0_prec
             XS_lo%nuFis=0.e0_prec
               XS_lo%Fis=0.e0_prec
               XS_lo%Chi=0.e0_prec
              XS_lo%kFis=0.e0_prec
               XS_lo%Rem=0.e0_prec
               XS_lo%Dif=0.e0_prec
             XS_lo%Scatt=0.e0_prec
            
            ! Data Reordering 
            do l=1, Opt_lo%n_tot
                i=mod(l-1,Opt_lo%n_x)+1
                j=int((l-i)/Opt_lo%n_x)+1
                
                
                m=LO_Mesh(l)%Mat_ID
                if (m<1) cycle
                
                do t=1, Opt_lo%n_g
                    ll=(t-1)*Opt_lo%n_tot+(j-1)*Opt_lo%n_x+i
                    
                    XS_lo%Tot(ll)  =Input_Materials(m)%XS_TO(t)
                    XS_lo%Tra(ll)  =Input_Materials(m)%XS_TR(t)
                    XS_lo%Absr(ll) =Input_Materials(m)%XS_AB(t)
                    XS_lo%nuFis(ll)=Input_Materials(m)%XS_nF(t)
                    XS_lo%Fis(ll)  =Input_Materials(m)%XS_FS(t)
                    XS_lo%Chi(ll)  =Input_Materials(m)%Chi(t)
                    XS_lo%kFis(ll) =Input_Materials(m)%XS_kF(t)
                    XS_lo%Rem(ll)  =Input_Materials(m)%XS_RM(t)
                    XS_lo%Dif(ll)  =Input_Materials(m)%XS_D(t)
                    do g=1, Opt_lo%n_g
                        XS_lo%Scatt(ll,g)=Input_Materials(m)%XS_SC(t,g)
                    end do
                    
                end do
            end do
            
            
            close(110)
            write(output_unit,'(A)') '...Done'
            write(files%log,'(A)') '...Done'
            write(output_unit,*) 
            write(files%log,*) 
            
            ! In-center and vertices printing 
            select case(Opt_lo%flag_ms) 
                case(0)
                    do l=1, Opt_lo%n_tot
                        write(files%coord_centers,'(2ES21.13)') LO_Mesh(l)%cent(:)
                    end do
                    do l=1, Opt_lo%n_tot
                        write(files%coord_vertices,'(8ES21.13)') LO_Mesh(l)%vert(:,1), LO_Mesh(l)%vert(:,2), LO_Mesh(l)%vert(:,3), LO_Mesh(l)%vert(:,4)
                    end do
                    
                case(1)
                    do l=1, Opt_lo%n_tot
                        write(files%coord_centers,'(2ES21.13)') LO_Mesh(l)%cent(:)
                    end do
                    do l=1, Opt_lo%n_tot
                        write(files%coord_vertices,'(8ES21.13)') LO_Mesh(l)%vert(:,1), LO_Mesh(l)%vert(:,2), LO_Mesh(l)%vert(:,3)
                    end do
                    
                case(2)
                    do l=1, Opt_lo%n_tot
                        write(files%coord_centers,'(2ES21.13)') LO_Mesh(l)%cent(:)
                    end do
                    do l=1, Opt_lo%n_tot
                        write(files%coord_vertices,'(8ES21.13)') LO_Mesh(l)%vert(:,1), LO_Mesh(l)%vert(:,2), LO_Mesh(l)%vert(:,3)
                    end do
                    
            end select
            
        
        end subroutine  Geometry_Card
        



        ! Meshers
        
        subroutine Quad_Mesher(box_mat_ID, Elem_det, Opt_lo, LO_Mesh, dx_gh, dy_gh) 
            
            integer, intent(in) :: box_mat_ID
            real(prec), intent(out) :: dx_gh, dy_gh
            type(Figure), allocatable, intent(in) :: Elem_det(:)
            type(Options_Data), intent(inout) :: Opt_lo
            type(LO_geom), allocatable, intent(inout) :: LO_Mesh(:)
            
            integer :: i, j, n, l, n_x, n_y, n_tot
            integer :: n_x_mesh, n_tot_mesh, no_faces
            real(prec) :: dx_homo, dy_homo, dx_part, dy_part, xmin, xmax, ymin, ymax, x_side, y_side
            
            
            ! Square grid box
            no_faces=size(Elem_det(box_mat_ID)%vertices,2)
            xmax=maxval(Elem_det(box_mat_ID)%vertices(1,:))
            xmin=minval(Elem_det(box_mat_ID)%vertices(1,:))
            ymax=maxval(Elem_det(box_mat_ID)%vertices(2,:))
            ymin=minval(Elem_det(box_mat_ID)%vertices(2,:))
            
            x_side=xmax-xmin
            y_side=ymax-ymin
            
            dx_homo=x_side/(Opt_lo%n_x-2)
            dy_homo=y_side/(Opt_lo%n_y-2)
            
            dx_gh=dx_homo
            dy_gh=dy_homo
            
            n_x=Opt_lo%n_x
            n_y=Opt_lo%n_y
            n_tot=Opt_lo%n_tot
            
            do l=1, n_tot
                i=mod(l-1,n_x)+1
                j=int((l-i)/n_x)+1
                
                LO_Mesh(l)%Sides_ID=[1, 2, 3, 4]
                LO_Mesh(l)%Sides_Neigh=[3, 4, 1, 2]
                LO_Mesh(l)%Cent(:)=0.e0_prec
                LO_Mesh(l)%Vert(:,:)=0.e0_prec
                LO_Mesh(l)%Neigh_ID(:)=-1
                
                ! Neighbouring cells from B, counter-clockwise 
            if   (j>1) LO_Mesh(l)%Neigh_ID(1)=l-n_x ! B
            if (i<n_x) LO_Mesh(l)%Neigh_ID(2)=l+1   ! R
            if (j<n_y) LO_Mesh(l)%Neigh_ID(3)=l+n_x ! T
            if   (i>1) LO_Mesh(l)%Neigh_ID(4)=l-1   ! L
                
                ! Mesh sides 
            LO_Mesh(l)%dx=dx_homo
            LO_Mesh(l)%dy=dy_homo
                
                ! Repr. size 
            LO_Mesh(l)%dl_lo(1)=LO_Mesh(l)%dy/2.e0_prec  ! B
            LO_Mesh(l)%dl_lo(2)=LO_Mesh(l)%dx/2.e0_prec  ! R
            LO_Mesh(l)%dl_lo(3)=LO_Mesh(l)%dy/2.e0_prec  ! T
            LO_Mesh(l)%dl_lo(4)=LO_Mesh(l)%dx/2.e0_prec  ! L
                
                ! Area 
            LO_Mesh(l)%A_lo(1)=LO_Mesh(l)%dx  ! B
            LO_Mesh(l)%A_lo(2)=LO_Mesh(l)%dy  ! R
            LO_Mesh(l)%A_lo(3)=LO_Mesh(l)%dx  ! T
            LO_Mesh(l)%A_lo(4)=LO_Mesh(l)%dy  ! L
                
                ! Volume 
            LO_Mesh(l)%V_lo=dx_homo*dy_homo
                
                ! Centroids 
                LO_Mesh(l)%Cent(1)=LO_Mesh(l)%dx/2.e0_prec
                LO_Mesh(l)%Cent(2)=LO_Mesh(l)%dy/2.e0_prec
                
                ! From BL, counter-clockwise 
                ! x 
                LO_Mesh(l)%Vert(1,1)= 0.e0_prec          ! BL
                LO_Mesh(l)%Vert(1,2)=LO_Mesh(l)%dx  ! BR
                LO_Mesh(l)%Vert(1,3)=LO_Mesh(l)%dx  ! TR
                LO_Mesh(l)%Vert(1,4)= 0.e0_prec          ! TL
                ! y 
                LO_Mesh(l)%Vert(2,1)= 0.e0_prec          ! BL
                LO_Mesh(l)%Vert(2,2)= 0.e0_prec          ! BR
                LO_Mesh(l)%Vert(2,3)=LO_Mesh(l)%dy  ! TR
                LO_Mesh(l)%Vert(2,4)=LO_Mesh(l)%dy  ! TL
                
                ! From relative to absolute 
                ! x
                if (i>1) then
                    dx_part=sum(LO_Mesh((j-1)*n_x+1:(j-1)*n_x+(i-1))%dx)
                    LO_Mesh(l)%Cent(1)=LO_Mesh(l)%Cent(1)+dx_part
                    LO_Mesh(l)%Vert(1,:)=LO_Mesh(l)%Vert(1,:)+dx_part
                end if
                ! y
                if (j>1) then
                    dy_part=sum(LO_Mesh(1:(j-1)*n_x:n_x)%dy)
                    LO_Mesh(l)%Cent(2)=LO_Mesh(l)%Cent(2)+dy_part
                    LO_Mesh(l)%Vert(2,:)=LO_Mesh(l)%Vert(2,:)+dy_part
                end if
            end do
            
            
            
        end subroutine Quad_Mesher
    
        subroutine TriIso_Mesher(box_mat_ID, Elem_det, Opt_lo, LO_Mesh, dx_gh, dy_gh) 
                ! Mesh grid visualization (problem-tailored) 
            !
            !-------------------
            !|/|\|/|\|/|\|/|\|/|
            !|\|/|\|/|\|/|\|/|\|
            !|/|\|/|\|/|\|/|\|/|
            !|\|/|\|/|\|/|\|/|\|
            !|/|\|/|\|/|\|/|\|/|
            !|\|/|\|/|\|/|\|/|\|
            !-------------------

            ! Mesh types:
            !
            ! |\
            ! |  \      ! BL mesh
            ! |____\  
            !
            !     /|
            !   /  |    ! BR mesh
            ! /____|  
            ! ______
            ! \    |
            !   \  |    ! TR mesh
            !     \|  
            ! ______
            ! |    /    
            ! |  /      ! TL mesh
            ! |/  
            integer, intent(in) :: box_mat_ID
            real(prec), intent(out) :: dx_gh, dy_gh
            type(Figure), allocatable, intent(in) :: Elem_det(:)
            type(Options_Data), intent(inout) :: Opt_lo
            type(LO_geom), allocatable, intent(inout) :: LO_Mesh(:)
            
            integer :: i, j, n, l, no_faces, n_x, n_y, n_tot
            real(prec) :: dx_homo, dy_homo, di_homo, xmin, xmax, ymin, ymax, x_side, y_side
                                                                           
            
            ! Square grid box
            no_faces=size(Elem_det(box_mat_ID)%vertices,2)
            xmax=maxval(Elem_det(box_mat_ID)%vertices(1,:))
            xmin=minval(Elem_det(box_mat_ID)%vertices(1,:))
            ymax=maxval(Elem_det(box_mat_ID)%vertices(2,:))
            ymin=minval(Elem_det(box_mat_ID)%vertices(2,:))
            
            x_side=xmax-xmin
            y_side=ymax-ymin
            
            ! Homogeneous mesh 
            dx_homo=x_side/(Opt_lo%n_x-2) ! First equal side
            dy_homo=dx_homo ! Second equal side
            di_homo=dx_homo*sqrt(2.e0_prec)
            
            dx_gh=dx_homo
            dy_gh=dx_homo
            
            n_x=Opt_lo%n_x
            n_y=Opt_lo%n_y
            n_tot=Opt_lo%n_tot
            
            
            

            
            
            ! BR THE FIRST TRIANGLE
            ! In-center & Vertices calculation (Absolute + Relative in subroutine) 
            do j=1, Opt_lo%n_y/2, 2
                do i=1, Opt_lo%n_x, 2      ! Intentional periodicity for the repeating block
                    n=(j-1)*2*Opt_lo%n_x+i
                    
                    ! BR, First triangle of the first periodic row 
                    l=n               
                    LO_Mesh(l)%dx=dx_homo
                    LO_Mesh(l)%dy=dy_homo
                    LO_Mesh(l)%ori='BR'
                    LO_Mesh(l)%Sides_ID=[1, 2, 6]
                    LO_Mesh(l)%Sides_Neigh=[3, 3, 1]
                    LO_Mesh(l)%Cent(1)=(i-1)*LO_Mesh(l)%dx
                    LO_Mesh(l)%Cent(2)=(j-1)*LO_Mesh(l)%dy
                    LO_Mesh(l)%Vert(1,:)=(i-1)*LO_Mesh(l)%dx
                    LO_Mesh(l)%Vert(2,:)=(j-1)*LO_Mesh(l)%dy
                    
                    
                    call BR_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    
                    ! BL, Second triangle of the first periodic row 
                    l=n+1            
                    if (i<n_x) then ! Additional control in case of odd n_x
                        LO_Mesh(l)%dx=dx_homo
                        LO_Mesh(l)%dy=dy_homo
                        LO_Mesh(l)%ori='BL'
                        LO_Mesh(l)%Sides_ID=[1, 5, 4]
                        LO_Mesh(l)%Sides_Neigh=[2, 1, 2]
                        LO_Mesh(l)%Cent(1)=      i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Cent(2)=  (j-1)*LO_Mesh(l)%dy
                        LO_Mesh(l)%Vert(1,:)=    i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Vert(2,:)=(j-1)*LO_Mesh(l)%dy
                        
                        call BL_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    end if
                    
                    ! TL, First triangle of the second periodic row 
                    l=n+n_x
                    LO_Mesh(l)%dx=dx_homo
                    LO_Mesh(l)%dy=dy_homo
                    LO_Mesh(l)%ori='TL'
                    LO_Mesh(l)%Sides_ID=[8, 3, 4]
                    LO_Mesh(l)%Sides_Neigh=[3, 1, 2]
                    LO_Mesh(l)%Cent(1)=(i-1)*LO_Mesh(l)%dx
                    LO_Mesh(l)%Cent(2)=(j-1)*LO_Mesh(l)%dy
                    LO_Mesh(l)%Vert(1,:)=(i-1)*LO_Mesh(l)%dx
                    LO_Mesh(l)%Vert(2,:)=(j-1)*LO_Mesh(l)%dy
                    
                    call TL_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    
                    ! TR, Second triangle of the second periodic row 
                    l=n+n_x+1      
                    if (i<n_x) then ! Additional control in case of odd n_x
                        LO_Mesh(l)%dx=dx_homo
                        LO_Mesh(l)%dy=dy_homo
                        LO_Mesh(l)%ori='TR'
                        LO_Mesh(l)%Sides_ID=[7, 2, 3]
                        LO_Mesh(l)%Sides_Neigh=[2, 3, 1]
                        LO_Mesh(l)%Cent(1)=    i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Cent(2)=(j-1)*LO_Mesh(l)%dy
                        LO_Mesh(l)%Vert(1,:)=    i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Vert(2,:)=(j-1)*LO_Mesh(l)%dy
                        
                        call TR_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    end if
                    
                    ! BL, First triangle of the third periodic row 
                    l=n+2*n_x+0   
                    if (j<n_y/2) then
                        LO_Mesh(l)%dx=dx_homo
                        LO_Mesh(l)%dy=dy_homo
                        LO_Mesh(l)%ori='BL'
                        LO_Mesh(l)%Sides_ID=[1, 5, 4]
                        LO_Mesh(l)%Sides_Neigh=[2, 1, 2]
                        LO_Mesh(l)%Cent(1)=(i-1)*LO_Mesh(l)%dx
                        LO_Mesh(l)%Cent(2)=    j*LO_Mesh(l)%dy
                        LO_Mesh(l)%Vert(1,:)=(i-1)*LO_Mesh(l)%dx
                        LO_Mesh(l)%Vert(2,:)=    j*LO_Mesh(l)%dy
                        
                        call BL_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    end if
                    
                    ! BR, Second triangle of the third periodic row 
                    l=n+2*n_x+1   
                    if (j<n_y/2 .AND. i<n_x) then
                        LO_Mesh(l)%dx=dx_homo
                        LO_Mesh(l)%dy=dy_homo
                        LO_Mesh(l)%ori='BR'
                        LO_Mesh(l)%Sides_ID=[1, 2, 6]
                        LO_Mesh(l)%Sides_Neigh=[3, 3, 1]
                        LO_Mesh(l)%Cent(1)=    i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Cent(2)=    j*LO_Mesh(l)%dy
                        LO_Mesh(l)%Vert(1,:)=  i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Vert(2,:)=  j*LO_Mesh(l)%dy
                        
                        call BR_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    end if
                    
                    ! TR, First triangle of the fourth periodic row 
                    l=n+3*n_x+0
                    if (j<n_y/2) then
                        LO_Mesh(l)%dx=dx_homo
                        LO_Mesh(l)%dy=dy_homo
                        LO_Mesh(l)%ori='TR'
                        LO_Mesh(l)%Sides_ID=[7, 2, 3]
                        LO_Mesh(l)%Sides_Neigh=[2, 3, 1]
                        LO_Mesh(l)%Cent(1)=(i-1)*LO_Mesh(l)%dx
                        LO_Mesh(l)%Cent(2)=    j*LO_Mesh(l)%dy
                        LO_Mesh(l)%Vert(1,:)=(i-1)*LO_Mesh(l)%dx
                        LO_Mesh(l)%Vert(2,:)=    j*LO_Mesh(l)%dy
                        
                        call TR_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    end if
                
                    ! TL, Second triangle of the fourth periodic row 
                    l=n+3*n_x+1
                    if (j<n_y/2 .AND. i<n_x) then
                        LO_Mesh(l)%dx=dx_homo
                        LO_Mesh(l)%dy=dy_homo     
                        LO_Mesh(l)%ori='TL'
                        LO_Mesh(l)%Sides_ID=[8, 3, 4]
                        LO_Mesh(l)%Sides_Neigh=[3, 1, 2]
                        LO_Mesh(l)%Cent(1)=    i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Cent(2)=    j*LO_Mesh(l)%dy
                        LO_Mesh(l)%Vert(1,:)=    i*LO_Mesh(l)%dx
                        LO_Mesh(l)%Vert(2,:)=    j*LO_Mesh(l)%dy
                        
                        call TL_IsoTrng(LO_Mesh(l)%dx, LO_Mesh(l)%dy, LO_Mesh(l)%Cent(:), LO_Mesh(l)%Vert(:,:))
                    end if
                
                end do
            end do
            
            
            
            
            
            
            
            ! Mesh discrimination 
            do l=1, n_tot 
                i=mod(l-1,n_x)+1
                j=int((l-i)/n_x)+1
                
                
                ! Manual in-center & Volume calculations FOR ISOSCELE TRIANGLES 
                LO_Mesh(l)%V_lo=LO_Mesh(l)%dx*LO_Mesh(l)%dy/2.e0_prec
                
                ! Mesh orientation discrimination 
                select case(LO_Mesh(l)%ori) 
                    
                    case('BL') 
                        ! Neighbouring cells from B, counter-clockwise
                        if          (j>1) LO_Mesh(l)%Neigh_ID(1)=l-n_x     ! B
                        if (l<=n_tot-n_x) LO_Mesh(l)%Neigh_ID(2)=l+n_x     ! TR
                        if          (i>1) LO_Mesh(l)%Neigh_ID(3)=l-1       ! L
                        
                        ! Repr. Size
                        LO_Mesh(l)%dl_lo(1)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec          ! B
                        LO_Mesh(l)%dl_lo(2)=LO_Mesh(l)%dx*sqrt((1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2+(1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2)           ! TR
                        LO_Mesh(l)%dl_lo(3)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec          ! L
                        
                        ! Area 
                        LO_Mesh(l)%A_lo(1)=LO_Mesh(l)%dx               ! B
                        LO_Mesh(l)%A_lo(2)=LO_Mesh(l)%dx*sqrt(2.e0_prec)    ! TR
                        LO_Mesh(l)%A_lo(3)=LO_Mesh(l)%dx               ! L
                        
                    case('BR') 
                        ! Neighbouring cells from B, counter-clockwise
                        if          (j>1) LO_Mesh(l)%Neigh_ID(1)=l-n_x     ! B
                        if        (i<n_x) LO_Mesh(l)%Neigh_ID(2)=l+1       ! R
                        if (l<=n_tot-n_x) LO_Mesh(l)%Neigh_ID(3)=l+n_x     ! TL
                        
                        ! Repr. Size
                        LO_Mesh(l)%dl_lo(1)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec     ! B
                        LO_Mesh(l)%dl_lo(2)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec     ! R
                        LO_Mesh(l)%dl_lo(3)=LO_Mesh(l)%dx*sqrt((1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2+(1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2)       ! TL
                        
                        ! Area 
                        LO_Mesh(l)%A_lo(1)=LO_Mesh(l)%dx               ! B
                        LO_Mesh(l)%A_lo(2)=LO_Mesh(l)%dx               ! R
                        LO_Mesh(l)%A_lo(3)=LO_Mesh(l)%dx*sqrt(2.e0_prec)    ! TL
                        
                    case('TR') 
                        ! Neighbouring cells from BL, counter-clockwise
                        if   (j>1) LO_Mesh(l)%Neigh_ID(1)=l-n_x     ! BL
                        if (i<n_x) LO_Mesh(l)%Neigh_ID(2)=l+1       ! R
                        if (j<n_y) LO_Mesh(l)%Neigh_ID(3)=l+n_x     ! T
                        
                        ! Repr. Size
                        LO_Mesh(l)%dl_lo(1)=LO_Mesh(l)%dx*sqrt((1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2+(1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2)       ! BL
                        LO_Mesh(l)%dl_lo(2)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec  ! R
                        LO_Mesh(l)%dl_lo(3)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec  ! T
                        
                        ! Area 
                        LO_Mesh(l)%A_lo(1)=LO_Mesh(l)%dx*sqrt(2.e0_prec)    ! BL
                        LO_Mesh(l)%A_lo(2)=LO_Mesh(l)%dx               ! R
                        LO_Mesh(l)%A_lo(3)=LO_Mesh(l)%dx               ! T
                        
                    case('TL') 
                        ! Neighbouring cells from BR, counter-clockwise
                        if   (j>1) LO_Mesh(l)%Neigh_ID(1)=l-n_x     ! BR
                        if (j<n_y) LO_Mesh(l)%Neigh_ID(2)=l+n_x     ! T
                        if   (i>1) LO_Mesh(l)%Neigh_ID(3)=l-1       ! L
                        
                        ! Repr. Size
                        LO_Mesh(l)%dl_lo(1)=LO_Mesh(l)%dx*sqrt((1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2+(1.e0_prec/3.e0_prec-1.e0_prec/2.e0_prec)**2)       ! BL
                        LO_Mesh(l)%dl_lo(2)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec  ! R
                        LO_Mesh(l)%dl_lo(3)=LO_Mesh(l)%dx*1.e0_prec/3.e0_prec  ! T
                        
                        ! Area 
                        LO_Mesh(l)%A_lo(1)=LO_Mesh(l)%dx*sqrt(2.e0_prec)    ! BR
                        LO_Mesh(l)%A_lo(2)=LO_Mesh(l)%dx               ! T
                        LO_Mesh(l)%A_lo(3)=LO_Mesh(l)%dx               ! L
                        
                    case default
                        write(output_unit,'(A)') " ERROR! Wrong Triangle orientation - Mesh generation!"
                        write(error_unit,'(A)') " ERROR! Wrong Triangle orientation - Mesh generation!"
                    end select
            
            end do
        
        end subroutine TriIso_Mesher
    
        ! Subroutines for TriIso_Mesher    
        subroutine BL_IsoTrng(dx, dy, cent, vert) 
            real(prec), intent(inout) :: cent(:), vert(:,:)
            real(prec), intent(in) :: dx, dy
            
              !cent(1)=  cent(1)+dx*1.e0_prec/3.e0_prec
              !cent(2)=  cent(2)+dy*1.e0_prec/3.e0_prec
              cent(1)=  cent(1)+dx*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec)
              cent(2)=  cent(2)+dy*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec)
            vert(1,1)=vert(1,1)+0.e0_prec
            vert(2,1)=vert(2,1)+0.e0_prec
            vert(1,2)=vert(1,2)+  dx
            vert(2,2)=vert(2,2)+0.e0_prec
            vert(1,3)=vert(1,3)+0.e0_prec
            vert(2,3)=vert(2,3)+  dy  
        end subroutine BL_IsoTrng
        
        subroutine BR_IsoTrng(dx, dy, cent, vert) 
            real(prec), intent(inout) :: cent(:), vert(:,:)
            real(prec), intent(in) :: dx, dy
            
              !cent(1)=  cent(1)+dx*2.e0_prec/3.e0_prec
              !cent(2)=  cent(2)+dy*1.e0_prec/3.e0_prec
              cent(1)=  cent(1)+dx*sqrt(2.e0_prec)/2.e0_prec          ! The in-radius is dx*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec), which means the in-center
              cent(2)=  cent(2)+dy*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec)   ! is at sqrt(2.e0_prec)/2.e0_prec if the triangle is rotated
            vert(1,1)=vert(1,1)+0.e0_prec
            vert(2,1)=vert(2,1)+0.e0_prec
            vert(1,2)=vert(1,2)+  dx
            vert(2,2)=vert(2,2)+0.e0_prec
            vert(1,3)=vert(1,3)+  dx
            vert(2,3)=vert(2,3)+  dy  
        end subroutine BR_IsoTrng
        
        subroutine TR_IsoTrng(dx, dy, cent, vert) 
            real(prec), intent(inout) :: cent(:), vert(:,:)
            real(prec), intent(in) :: dx, dy
            
              !cent(1)=  cent(1)+dx*2.e0_prec/3.e0_prec
              !cent(2)=  cent(2)+dy*2.e0_prec/3.e0_prec
              cent(1)=  cent(1)+dx*sqrt(2.e0_prec)/2.e0_prec ! The in-radius is dx*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec), which means the in-center
              cent(2)=  cent(2)+dy*sqrt(2.e0_prec)/2.e0_prec ! is at sqrt(2.e0_prec)/2.e0_prec if the triangle is rotated
            vert(1,1)=vert(1,1)+  dx
            vert(2,1)=vert(2,1)+0.e0_prec
            vert(1,2)=vert(1,2)+  dx
            vert(2,2)=vert(2,2)+  dy
            vert(1,3)=vert(1,3)+0.e0_prec
            vert(2,3)=vert(2,3)+  dy 
        end subroutine TR_IsoTrng
        
        subroutine TL_IsoTrng(dx, dy, cent, vert) 
            real(prec), intent(inout) :: cent(:), vert(:,:)
            real(prec), intent(in) :: dx, dy
            
              !cent(1)=  cent(1)+dx*1.e0_prec/3.e0_prec
              !cent(2)=  cent(2)+dy*2.e0_prec/3.e0_prec
              cent(1)=  cent(1)+dx*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec) ! The in-radius is dx*(1.e0_prec-sqrt(2.e0_prec)/2.e0_prec), which means the in-center
              cent(2)=  cent(2)+dy*sqrt(2.e0_prec)/2.e0_prec        ! is at sqrt(2.e0_prec)/2.e0_prec if the triangle is rotated
            vert(1,1)=vert(1,1)+0.e0_prec
            vert(2,1)=vert(2,1)+0.e0_prec
            vert(1,2)=vert(1,2)+  dx
            vert(2,2)=vert(2,2)+  dy
            vert(1,3)=vert(1,3)+0.e0_prec
            vert(2,3)=vert(2,3)+  dy  
        end subroutine TL_IsoTrng
    
    
    
    
    
    
    
    subroutine TriEq_Mesher(n_x, n_y, dx_m, dy_m, n_tot, ID, ori, x_max, y_max, cent, vert, dl, A, Vol, Neigh) 
        integer, intent(in) :: n_x, n_y, n_tot, ID(:)
        real(prec), intent(in) :: x_max, y_max
        
        integer, intent(out) ::  Neigh(:,:)
        real(prec), intent(out) :: cent(:,:), vert(:,:,:), dl(:,:), A(:,:), Vol(:), dx_m, dy_m
        character(len=2), allocatable, intent(inout) :: ori(:)
        
        integer :: i, j, n, l
        real(prec) :: dx, dy
        
        
        ! Homogeneous mesh 
        dx=x_max/(n_x-2)
        dy=dx*sin(pi/3.e0_prec)
        
        dx_m=dx
        dy_m=dy
        
        cent=0.e0_prec
        vert=0.e0_prec
        Neigh=-1
        
        ! Absolute Centroid & Vertices calculation 
        do j=1, n_y/2 
            do i=1, n_x
                
                ! First periodicity
                n=(j-1)*2*n_x+i
                l=ID(n)
                cent(l,1)=(j-1)*dx/2.e0_prec+(i-1)*dx
                cent(l,2)=(j-1)*dy
                vert(l,1,:)=(j-1)*dx/2.e0_prec+(i-1)*dx
                vert(l,2,:)=(j-1)*dy
                ori(l)='U'
                
                ! Second periodicity
                n=(j-1)*2*n_x+n_x+i
                l=ID(n)
                cent(l,1)=(j-1)*dx/2.e0_prec+(i-1)*dx+dx/2.e0_prec
                cent(l,2)=(j-1)*dy
                vert(l,1,:)=(j-1)*dx/2.e0_prec+(i-1)*dx+dx/2.e0_prec
                vert(l,2,:)=(j-1)*dy
                ori(l)='D'
            end do
        end do
        
        ! Relative position 
        do n=1, n_tot 
            l=ID(n)
            i=mod(l-1,n_x)+1
            j=int((l-i)/n_x)+1
            
            select case(ori(l)) 
                case('U') 
                    
                    ! Neighbouring cells from B, counter-clockwise
                    if (j>1)    Neigh(l,1)=l-n_x    ! B
                                Neigh(l,2)=l+n_x    ! TR
                    if (i>1)    Neigh(l,3)=l+n_x-1  ! TL
                    
                    ! From B, counter-clockwise
                    A(l,1)=dx   ! B
                    A(l,2)=dx   ! TR
                    A(l,3)=dx   ! TL

                    dl(l,1)=dx*sqrt(3.e0_prec)/6.e0_prec  ! B
                    dl(l,2)=dx*sqrt(3.e0_prec)/6.e0_prec  ! TR
                    dl(l,3)=dx*sqrt(3.e0_prec)/6.e0_prec  ! TL
                    
                    ! Centroids
                    cent(l,1)=cent(l,1)+dx/2.e0_prec
                    cent(l,2)=cent(l,2)+sqrt(3.e0_prec)/6.e0_prec*dx
                    
                    ! Vertices
                    vert(l,1,1)=vert(l,1,1)+0.e0_prec    ! BL - x
                    vert(l,2,1)=vert(l,2,1)+0.e0_prec    ! BL - y
                    vert(l,1,2)=vert(l,1,2)+dx      ! BR - x
                    vert(l,2,2)=vert(l,2,2)+0.e0_prec    ! BR - y
                    vert(l,1,3)=vert(l,1,3)+dx/2.e0_prec ! T  - x
                    vert(l,2,3)=vert(l,2,3)+dy      ! T  - y
                    
                case('D') 
                    
                    ! Neighbouring cells from B, counter-clockwise
                    if (i<n_x)  Neigh(l,1)=l-n_x    ! BL
                                Neigh(l,2)=l-n_x+1  ! BR
                    if (i>1)    Neigh(l,3)=l+n_x    ! T

                    dl(l,1)=dx*sqrt(3.e0_prec)/6.e0_prec  ! BL
                    dl(l,2)=dx*sqrt(3.e0_prec)/6.e0_prec  ! BR
                    dl(l,3)=dx*sqrt(3.e0_prec)/6.e0_prec  ! T
                    
                    ! Centroids
                    cent(l,1)=cent(l,1)+dx/2.e0_prec
                    cent(l,2)=cent(l,2)+(dy-sqrt(3.e0_prec)/6.e0_prec*dx)
                    
                    ! Vertices
                    vert(l,1,1)=vert(l,1,1)+dx/2.e0_prec ! B  - x
                    vert(l,2,1)=vert(l,2,1)+0.e0_prec    ! B  - y
                    vert(l,1,2)=vert(l,1,2)+dx      ! TR - x
                    vert(l,2,2)=vert(l,2,2)+dy      ! TR - y
                    vert(l,1,3)=vert(l,1,3)+0.e0_prec    ! TL - x
                    vert(l,2,3)=vert(l,2,3)+dy      ! TL - y
                    
                case default 
                    write(output_unit,'(A)') " ERROR! Wrong mesh orientation!"
                    write(error_unit,'(A)') " ERROR! Wrong mesh orientation!"
                    exit
            end select
            
            Vol(l)=dx*dy/2.e0_prec
        end do
        
    end subroutine TriEq_Mesher
    
    
    
    
    
    
    
    
    
end module  Input_Reader
