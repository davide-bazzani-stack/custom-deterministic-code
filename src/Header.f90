program Deterministic_FVM_Code
    use Variables
    use IO_module, only: files
    use Input_Reader
    use Service_Fcns
    use FVM_Quad_Alg
    use FVM_TriIso_Alg
    
        implicit none
    
    
    
    ! Options variables 
    integer :: flag_ms, flag_accel, flag_SOR_lo, flag_SOR_gl
    integer :: iter_in_max_lo, iter_out_max_lo, iter_solver_max_lo
    integer :: iter_in_max_gl, iter_out_max_gl, iter_solver_max_gl
    
    real(dp) :: tol_solv_lo, tol0_lo, tol1_lo, tol2_lo, omega_SOR_lo
    real(dp) :: tol_solv_gl, tol0_gl, tol1_gl, tol2_gl, omega_SOR_gl
    
    character(len=10)  :: mesh_strct, acc_type    
    
    type(Material), allocatable :: Input_Materials(:)
    type(Figure), allocatable :: Elem_det(:), Elem_CMFD(:)
    type(Options_Data) :: Opt_lo, Opt_gl
    
    ! Geometry & Meshing variables 
    integer :: n_x_lo, n_y_lo, n_tot_lo, n_g_lo, n_mat_lo
    real(dp) :: dx_gh, dy_gh
    integer :: n_x_gl, n_y_gl, n_tot_gl, n_g_gl
    integer, allocatable :: IDs_lo(:), IDs_neigh_lo(:,:), BC_lo(:), J_map_lo(:,:), Mat_lo(:)
    real(dp), allocatable :: Centr_lo(:,:), V_lo(:), Vert_lo(:,:,:), A_lo(:,:), dl_lo(:,:), Gam_lo(:)
    character(len=2), allocatable :: Mesh_Ornt_lo(:)
    
    ! Cross Sections 
    real(dp), allocatable :: XS_TO_lo(:), XS_TR_lo(:), XS_AB_lo(:), XS_nF_lo(:), XS_FS(:), Chi_lo(:), XS_kF(:), XS_RM_lo(:)
    real(dp), allocatable :: XS_D_lo(:), XS_SC_lo(:,:)
    
    type(XS_Data) :: XS_lo, XS_gl
    
    ! Local problem Solver-related variables 
    integer :: n_red_lo
    integer, allocatable :: IDs_lo_Red(:)
    real(dp), allocatable :: Phi_lo(:), Phi_lo_out(:), J_lo(:,:)
    
    ! Global problem Solver-related variables 
    integer, allocatable :: lo_to_gl_V(:), lo_to_gl_E(:)
    real(dp), allocatable :: A_gl(:,:), V_gl(:), dl_gl(:,:), Phi_gl(:), Phi_old_gl(:), Phi_temp_gl(:), Phi_avg_gl(:)
    real(dp), allocatable :: XS_TO_gl(:), XS_TR_gl(:), XS_AB_gl(:), XS_nF_gl(:), XS_FS_gl(:), Chi_gl(:), XS_SC_gl(:,:)
    real(dp), allocatable :: XS_kF_gl(:), XS_RM_gl(:), XS_Di_gl(:), Q_ext_gl(:), D_Til_gl(:,:), D_TiL_El_gl(:,:), D_Hat_gl(:,:)
    real(dp), allocatable :: D_Hat_BC_gl(:,:), J_Par_gl(:,:), J_Tot_gl(:,:), source_gl(:), S_gl(:), S_old_gl(:), MM_gl(:,:)
    real(dp), allocatable :: J_El_gl(:,:), J_El_BC_gl(:,:)
    
    real(dp) :: dl_input, dx_input
    
    type(GL_geom), allocatable :: HGCMFD_Mesh(:)
    type(GL_coeff), allocatable :: HGCMFD_Param(:)
    type(Accel_Vars_Vect) :: Serv_Vect 
    type(Accel_Vars_Matr) :: Serv_Matr
    
    
    
    type(LO_geom), allocatable :: LO_Mesh(:)
    type(LO_coeff), allocatable :: LO_Param(:)
    type(GL_geom), allocatable :: GL_Mesh(:)
    type(GL_coeff), allocatable :: GL_Param(:) 
    

    logical :: output__is_enabled = .true.
    logical :: debug__is_enabled = .true.
    integer :: l 
    
    
    
    
    

    
    ! Output files
    call files%open_IO(output__is_enabled, debug__is_enabled)
    
    ! Simulation options
    call Options_Card(flag_accel, Opt_lo, Opt_gl, dl_input, dx_input) 
    
    ! Material data (# of materials & XSs)
    call Material_Card(n_mat_lo, Opt_lo%n_g, Input_Materials)
    
    ! Geometry declaration & array maps
    call Geometry_Card(Opt_lo, Input_Materials, Elem_det, XS_lo, LO_Mesh, LO_Param, dx_gh, dy_gh)
    
    ! Acceleration initialization
    if (flag_accel .NE. 0) call Accel_Ini(flag_accel, Opt_lo, Opt_gl, Elem_det, Elem_CMFD, dx_gh, dy_gh, dl_input, dx_input, XS_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, Serv_Vect, Serv_Matr)
    
    
    ! FVM Solver logic
    select case (Opt_lo%flag_ms) 
        case(0)
        call FVM_Quad_Solver(Elem_CMFD, flag_accel, Opt_lo, Opt_gl, LO_Mesh, LO_Param, XS_lo, XS_gl, GL_Mesh, GL_Param, Serv_Vect, Serv_Matr, Phi_lo) 
        
        case(1)
        call FVM_TriIso_Solver(Elem_CMFD, flag_accel, Opt_lo, Opt_gl, LO_Mesh, LO_Param, XS_lo, XS_gl, GL_Mesh, GL_Param, Serv_Vect, Serv_Matr, Phi_lo)
        
        case(2)
        !call FVM_TriEq_MM()
    
    end select
    
        
    allocate(Phi_lo_out(Opt_lo%n_tot*Opt_lo%n_g))
    
    !call Output_Converter(n_g_lo, n_x_lo, IDs_lo_Red, Phi_lo, Phi_lo_out)      call Flux_Vol_Normalization(n_g_lo, V_lo, Phi_lo)
    
    call Results_Writing(Opt_lo%n_g, Opt_lo%n_x, Opt_lo%n_y, Opt_lo%n_tot, [(LO_Mesh(l)%univ, l=1, Opt_lo%n_tot)], Phi_lo, XS_lo%kFis, [(LO_Mesh(l)%V_lo, l=1, Opt_lo%n_tot)]) 
    
    
    call files%close_IO(output__is_enabled, debug__is_enabled)
    
    
end program Deterministic_FVM_Code

