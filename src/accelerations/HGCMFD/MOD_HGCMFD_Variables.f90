module HGCMFD_Variables
    
        implicit none
    
    
    integer, parameter :: dprec=kind(1.d0)
    
    real(dprec) :: epsi=tiny(1.d0)
    
    type :: XS_Data_2 
        real(dprec), allocatable :: Tot(:), Tra(:), Absr(:), nuFis(:), Fis(:), Scatt(:,:), Chi(:), kFis(:), Rem(:), Dif(:)
    end type
    
    type :: Options_Data_2
        integer :: flag_ms, n_x, n_y, n_tot, n_g, n_red
        integer :: it_sol_max, it_out_max, it_in_max
        real(dprec) :: tol_solv, tol0, tol1, tol2 
        integer :: flag_SOR
        real(dprec) :: w_SOR
    end type
    
    type :: LO_geom_2 
        integer :: ID, ID_Red
        integer, allocatable :: Sides_ID(:) 
        integer, allocatable :: Sides_Neigh(:) 
        real(dprec), allocatable :: Cent(:), Vert(:,:)
        integer, allocatable :: Neigh_ID(:) 
        integer, allocatable :: BC(:)     ! [no_faces]
        integer, allocatable :: J_flag(:)     ! [no_faces]
        integer, allocatable :: J_lo_gl_ID(:), J_lo_gl_ID_Red(:)     ! [no_faces]
        real(dprec) :: dx, dy     ! Representative sides of the mesh. dx is the side mesh on x, dy on y  
        real(dprec), allocatable :: dl_lo(:), A_lo(:)     ! [no_faces], [no_faces] 
        real(dprec) :: V_lo
        character(len=2) :: ori
        integer :: univ, lo_gl_Homo, Mat_ID
        integer, allocatable :: J_lo_gl_face(:) 
    end type
    
    type :: LO_coeff_2 
        integer :: En_ID, lo_gl_Cond
        real(dprec), allocatable :: Gam(:)     ! [no_faces] 
        real(dprec), allocatable :: J_part(:,:), J_net(:)
        real(dprec), allocatable :: D_til(:)
    end type
    
    type :: GL_geom_2 
        integer :: ID 
        integer, allocatable :: J_lo_gl_face(:), Neigh_ID(:), Sides_ID(:) 
        real(dprec), allocatable :: dl_gl(:), A_gl(:)
        real(dprec) :: V_gl
        integer, allocatable :: J_idx(:,:)
    end type
    
    type :: GL_coeff_2 
        real(dprec), allocatable :: D_til(:,:), D_hat(:,:)
        real(dprec), allocatable :: J_part(:,:), J_net(:)
        real(dprec), allocatable :: D_til_el, D_hat_el(:), J_Part_El(:)
        real(dprec) :: J_net_el
    end type
    
    type :: Accel_Vars_Vect_2
        real(dprec), allocatable :: S_fiss(:), S_fiss_old(:), source(:), Q_ext(:)
        real(dprec), allocatable :: Phi(:), Phi_old(:), Phi_temp(:), Phi_avg(:)
        real(dprec), allocatable :: a(:), b(:), c(:), d(:), e(:)
    end type
    
    type :: Accel_Vars_Matr_2 
        real(dprec), allocatable :: S_fiss(:), S_fiss_old(:), source(:), Q_ext(:)
        real(dprec), allocatable :: Phi(:), Phi_old(:), Phi_temp(:), Phi_avg(:), MM(:,:)
    end type
    
    
    
end module HGCMFD_Variables