module Variables
    use precision_kinds, only: prec

        implicit none
    
    private

    public :: Int_input_check, Flo_input_check, Char_input_check
    public :: Figure, Material, XS_Data, Options_Data
    public :: LO_geom, LO_coeff, GL_geom, GL_coeff
    public :: Accel_Vars_Vect, Accel_Vars_Matr

    type :: Int_input_check 
        integer :: value = 0
        logical  :: is_set = .false.
    end type
    
    type :: Flo_input_check 
        real(prec) :: value = 0.0_prec
        logical  :: is_set = .false.
    end type
    
    type :: Char_input_check 
        character(len=256) :: variable
        logical  :: is_set = .false.
    end type
    
    type :: Figure 
        integer :: ID, Fig_ID
        real(prec), allocatable :: centroid(:), vertices(:,:)
        character(len=256) :: material
        integer :: Mat_ID, universe, flag_BG
    end type
    
    type :: Material 
        integer :: ID
        character(len=256) :: name
        real(prec), allocatable :: XS_TO(:), XS_TR(:), XS_AB(:), XS_nF(:), XS_FS(:), XS_SC(:,:), Chi(:), XS_kF(:), XS_RM(:), XS_D(:)
    end type
        
    type :: XS_Data 
        real(prec), allocatable :: Tot(:), Tra(:), Absr(:), nuFis(:), Fis(:), Scatt(:,:), Chi(:), kFis(:), Rem(:), Dif(:)
    end type
    
    type :: Options_Data 
        integer :: flag_ms, n_x, n_y, n_tot, n_g, n_red
        integer :: it_sol_max, it_out_max, it_in_max
        real(prec) :: tol_solv, tol0, tol1, tol2
        character(len=16) :: preconditioner
        integer :: flag_SOR
        real(prec) :: w_SOR
    end type
    
    type :: LO_geom 
        integer :: ID, ID_Red
        integer, allocatable :: Sides_ID(:) 
        integer, allocatable :: Sides_Neigh(:) 
        real(prec), allocatable :: Cent(:), Vert(:,:)
        integer, allocatable :: Neigh_ID(:) 
        integer, allocatable :: BC(:)     ! [no_faces]
        integer, allocatable :: J_flag(:)     ! [no_faces]
        integer, allocatable :: J_lo_gl_ID(:), J_lo_gl_ID_Red(:)     ! [no_faces]
        real(prec) :: dx, dy     ! Representative sides of the mesh. dx is the side mesh on x, dy on y
        real(prec), allocatable :: dl_lo(:), A_lo(:)     ! [no_faces], [no_faces]
        real(prec) :: V_lo
        character(len=2) :: ori
        integer :: univ, lo_gl_Homo, Mat_ID
        integer, allocatable :: J_lo_gl_face(:) 
    end type
    
    type :: LO_coeff 
        integer :: En_ID, lo_gl_Cond
        real(prec), allocatable :: Gam(:)     ! [no_faces]
        real(prec), allocatable :: J_part(:,:), J_net(:)
        real(prec), allocatable :: D_til(:)
    end type
    
    type :: GL_geom 
        integer :: ID 
        integer, allocatable :: J_lo_gl_face(:), Neigh_ID(:), Sides_ID(:) 
        real(prec), allocatable :: dl_gl(:), A_gl(:)
        real(prec) :: V_gl
        integer, allocatable :: J_idx(:,:)
    end type
    
    type :: GL_coeff 
        real(prec), allocatable :: D_til(:,:), D_hat(:,:)
        real(prec), allocatable :: J_part(:,:), J_net(:)
        real(prec), allocatable :: D_til_el(:), D_hat_el(:), J_Part_El(:)
        real(prec) :: J_net_el
    end type
    
    type :: Accel_Vars_Vect 
        real(prec), allocatable :: S_fiss(:), S_fiss_old(:), source(:)
        real(prec), allocatable :: Phi(:), Phi_old(:), Phi_temp(:), Phi_avg(:)
        real(prec), allocatable :: a(:), b(:), c(:), d(:), e(:)
    end type
    
    type :: Accel_Vars_Matr 
        real(prec), allocatable :: S_fiss(:), S_fiss_old(:), source(:)
        real(prec), allocatable :: Phi(:), Phi_old(:), Phi_temp(:), Phi_avg(:), MM(:,:)
    end type
    
    
end module Variables
    
