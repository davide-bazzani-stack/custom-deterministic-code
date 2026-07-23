module FVM_TriIso_Alg
    use iso_fortran_env, only : output_unit, error_unit
    use precision_kinds, only: prec
    use IO_module, only : files


    use Variables
    use BiCGSTAB_Solver, only: banded_operator_t, bicgstab_options_t, &
        bicgstab_result_t, bicgstab_workspace_t, jacobi_preconditioner_t, &
        symmetric_gauss_seidel_preconditioner_t, ilu0_preconditioner_t, solve_bicgstab, &
        BICGSTAB_COMPONENT_SUCCESS
    use GaussElimination
    use CMFD_Alg
    use pCMFD_Alg
    !use lpCMFD_Alg
    use HGCMFD_Alg
    use pHGCMFD_Alg
    
        implicit none
    
    contains
    
    ! Header of the FVM solver       

    subroutine FVM_TriIso_Solver(Elem_CMFD, flag_accel, Opt_lo, Opt_gl, LO_Mesh, LO_Param, XS_lo, XS_gl, GL_Mesh, GL_Param, Serv_Vect, Serv_Matr, Phi_ext) 
        
        integer, intent(in) :: flag_accel
        
        type(Figure), intent(in) :: Elem_CMFD(:)
        type(Options_Data), intent(in) :: Opt_lo, Opt_gl
        type(XS_Data), intent(in) :: XS_lo
        type(XS_Data), intent(inout) :: XS_gl
        
        type(LO_geom), allocatable, intent(in) :: LO_Mesh(:)
        type(LO_coeff), allocatable, intent(inout) :: LO_Param(:)
        type(GL_geom), allocatable, intent(inout) :: GL_Mesh(:)
        type(GL_coeff), allocatable, intent(inout) :: GL_Param(:)
        type(Accel_Vars_Vect), intent(inout) :: Serv_Vect
        type(Accel_Vars_Matr), intent(inout) :: Serv_Matr
        
        real(prec), allocatable, intent(out) :: Phi_ext(:)
        !real(prec), allocatable:: Phi_old_ext(:)
        
        integer :: t, j, i, l, g, n, g_red, t_red, t_tot, iter_out_lo, iter_inn_lo, solver_status
        integer, allocatable :: IDs_Red_Arr(:)
        logical, allocatable :: preconditioner_is_ready(:)
        real(prec) :: k_lo, k_lo_old, err0, err1, err2, t_ini, t_fin
        real(prec), allocatable, target :: a_FVM(:), b_FVM(:), c_FVM(:), d_FVM(:), e_FVM(:)
        real(prec), allocatable :: XS_nF_Red(:), Chi_Red(:), XS_SC_Red(:,:)
        real(prec), allocatable :: J_lo(:,:), Phi_lo(:), Phi_lo_old(:), Phi_lo_temp(:), FSD(:), FSD_old(:), S_tot(:)
        type(banded_operator_t) :: migration_operator
        type(bicgstab_options_t) :: solver_options
        type(bicgstab_result_t) :: solver_result
        type(bicgstab_workspace_t) :: solver_workspace
        type(jacobi_preconditioner_t), allocatable :: jacobi_preconditioners(:)
        type(symmetric_gauss_seidel_preconditioner_t), allocatable :: sgs_preconditioners(:)
        type(ilu0_preconditioner_t), allocatable :: ilu0_preconditioners(:)

        ! Migration Matrix building
        call FVM_TriIso_MatrixBuilder(Opt_lo, LO_Mesh, LO_Param, XS_lo%Rem, XS_lo%Dif, a_FVM, b_FVM, c_FVM, d_FVM, e_FVM)


        ! Initialization 
        allocate(Phi_ext(Opt_lo%n_tot*Opt_lo%n_g))
        !allocate(Phi_old_ext(Opt_lo%n_tot*Opt_lo%n_g))
        allocate(Phi_lo_old(Opt_lo%n_red*Opt_lo%n_g))
        allocate(Phi_lo_temp(Opt_lo%n_red*Opt_lo%n_g))
        allocate(Phi_lo(Opt_lo%n_red*Opt_lo%n_g))
        allocate(FSD(Opt_lo%n_red))
        allocate(FSD_old(Opt_lo%n_red))
        allocate(S_tot(Opt_lo%n_red*Opt_lo%n_g))
        
        allocate(IDs_Red_Arr(Opt_lo%n_tot))
        
        IDs_Red_Arr=[(LO_mesh(l)%ID_Red, l=1, Opt_lo%n_tot)]
        k_lo=1.000
        Phi_lo=1.e0_prec
        Phi_ext=1.e0_prec
        err1=1.e0_prec
	    err2=1.e0_prec
        iter_out_lo=0
        solver_options%max_iterations=Opt_lo%it_sol_max
        solver_options%relative_tolerance=Opt_lo%tol_solv
        solver_options%absolute_tolerance=0.e0_prec
        allocate(preconditioner_is_ready(Opt_lo%n_g))
        preconditioner_is_ready=.false.
        select case (trim(Opt_lo%preconditioner))
            case ('Identity')
                continue
            case ('Jacobi')
                allocate(jacobi_preconditioners(Opt_lo%n_g))
            case ('SGS')
                allocate(sgs_preconditioners(Opt_lo%n_g))
            case ('ILU0')
                allocate(ilu0_preconditioners(Opt_lo%n_g))
            case default
                write(output_unit,'(2A)') 'ERROR! UNKNOWN BICGSTAB PRECONDITIONER: ', trim(Opt_lo%preconditioner)
                write(error_unit,'(2A)') 'ERROR! UNKNOWN BICGSTAB PRECONDITIONER: ', trim(Opt_lo%preconditioner)
                error stop 'Invalid BiCGSTAB preconditioner selection'
        end select
        call cpu_time(t_ini)
        
        ! Reduction for the proper source computation 
        allocate(XS_nF_Red(Opt_lo%n_g*Opt_lo%n_red))
        allocate(Chi_Red(Opt_lo%n_g*Opt_lo%n_red))
        allocate(XS_SC_Red(Opt_lo%n_g*Opt_lo%n_red,Opt_lo%n_g))
        
        n=0 
        do l=1, Opt_lo%n_tot
            if (LO_Mesh(l)%ID_Red>0) then
                n=n+1
                do t=1, Opt_lo%n_g
                    t_tot=(t-1)*Opt_lo%n_tot
                    t_red=(t-1)*Opt_lo%n_red
                    
                    XS_nF_Red(t_red+n)=XS_lo%nuFis(t_tot+l)
                    Chi_Red(t_red+n)=XS_lo%Chi(t_tot+l)
                    do g=1, Opt_lo%n_g
                        XS_SC_Red(t_red+n,g)=XS_lo%Scatt(t_tot+l,g)
                    end do
                end do
            end if
        end do
        
        
        write(output_unit,*) 
        write(files%log,*) 
        write(output_unit,'(A)') 'Start the FVM calculation (Isoscele triangle-mesh based):'
        write(files%log,'(A)') 'Start the FVM calculation (Isoscele triangle-mesh based):'
        
        
        ! FVM Outer Iterations 
        do while ((err1>Opt_lo%tol1 .OR. err2>Opt_lo%tol2) .AND. iter_out_lo<Opt_lo%it_out_max) 
            iter_out_lo=iter_out_lo+1
            k_lo_old=k_lo
            Phi_lo_old=Phi_lo
            !Phi_ext_old=Phi_ext
            
            ! Fission Source Distribution calculation 
            FSD_old=FissionSourceDistribution(Opt_lo%n_g, Opt_lo%n_red, Phi_lo_old, XS_nF_Red)
            
            ! Inner iterations (on the energy groups) 
            err0=1.e0_prec
            iter_inn_lo=0
            do while (err0>Opt_lo%tol0 .AND. iter_inn_lo<Opt_lo%it_in_max)
                Phi_lo_temp=Phi_lo
                
                ! Group sweeping 
                do t=1, Opt_lo%n_g
                    t_red=(t-1)*Opt_lo%n_red
                    S_tot(t_red+1:t_red+Opt_lo%n_red)=Total_Source_SingleGroup(t, Opt_lo%n_g, Opt_lo%n_red, k_lo_old, Chi_Red, FSD_old, XS_SC_Red, Phi_lo)
                    
                    call migration_operator%bind_triiso(Opt_lo%n_x-2, &
                        a_FVM(t_red+1:t_red+Opt_lo%n_red), b_FVM(t_red+1:t_red+Opt_lo%n_red), &
                        c_FVM(t_red+1:t_red+Opt_lo%n_red), d_FVM(t_red+1:t_red+Opt_lo%n_red), &
                        e_FVM(t_red+1:t_red+Opt_lo%n_red), solver_status)
                    if (solver_status==BICGSTAB_COMPONENT_SUCCESS) then
                        select case (trim(Opt_lo%preconditioner))
                            case ('Identity')
                                call solve_bicgstab(migration_operator, S_tot(t_red+1:t_red+Opt_lo%n_red), &
                                    Phi_lo(t_red+1:t_red+Opt_lo%n_red), solver_options, solver_workspace, solver_result)
                            case ('Jacobi')
                                if (.not. preconditioner_is_ready(t)) then
                                    call jacobi_preconditioners(t)%setup(migration_operator, solver_status)
                                    preconditioner_is_ready(t)=solver_status==BICGSTAB_COMPONENT_SUCCESS
                                end if
                                if (solver_status==BICGSTAB_COMPONENT_SUCCESS) then
                                    call solve_bicgstab(migration_operator, S_tot(t_red+1:t_red+Opt_lo%n_red), &
                                        Phi_lo(t_red+1:t_red+Opt_lo%n_red), solver_options, solver_workspace, &
                                        solver_result, jacobi_preconditioners(t))
                                end if
                            case ('SGS')
                                if (.not. preconditioner_is_ready(t)) then
                                    call sgs_preconditioners(t)%setup(migration_operator, solver_status)
                                    preconditioner_is_ready(t)=solver_status==BICGSTAB_COMPONENT_SUCCESS
                                end if
                                if (solver_status==BICGSTAB_COMPONENT_SUCCESS) then
                                    call solve_bicgstab(migration_operator, S_tot(t_red+1:t_red+Opt_lo%n_red), &
                                        Phi_lo(t_red+1:t_red+Opt_lo%n_red), solver_options, solver_workspace, &
                                        solver_result, sgs_preconditioners(t))
                                end if
                            case ('ILU0')
                                if (.not. preconditioner_is_ready(t)) then
                                    call ilu0_preconditioners(t)%setup(migration_operator, solver_status)
                                    preconditioner_is_ready(t)=solver_status==BICGSTAB_COMPONENT_SUCCESS
                                end if
                                if (solver_status==BICGSTAB_COMPONENT_SUCCESS) then
                                    call solve_bicgstab(migration_operator, S_tot(t_red+1:t_red+Opt_lo%n_red), &
                                        Phi_lo(t_red+1:t_red+Opt_lo%n_red), solver_options, solver_workspace, &
                                        solver_result, ilu0_preconditioners(t))
                                end if
                        end select
                    end if

                    if (solver_status/=BICGSTAB_COMPONENT_SUCCESS .or. .not. solver_result%converged) then
                        write(output_unit,'(A, I6.1, A)') 'ERROR! THE ', t,' GROUP DID NOT CONVERGE'
                        write(error_unit,'(A, I6.1, A)') 'ERROR! THE ', t,' GROUP DID NOT CONVERGE'
                        if (solver_status==BICGSTAB_COMPONENT_SUCCESS) then
                            write(error_unit,'(A, I0, 2A)') 'BiCGSTAB status ', solver_result%status, &
                                ': ', trim(solver_result%reason)
                        end if
                        iter_inn_lo=Opt_lo%it_in_max
                    end if
                end do
                err0=sqrt(dot_product(Phi_lo-Phi_lo_temp, Phi_lo-Phi_lo_temp)/dot_product(Phi_lo, Phi_lo))
            end do
            
            ! Inner iteration convergence check 
            if (err0>Opt_lo%tol0) then
                write(output_unit,'(A)') "MISSED INNER CONVERGENCE"
                write(error_unit,'(A)') "MISSED INNER CONVERGENCE"
            end if
            
            ! Fission Source Distribution calculation 
            FSD=FissionSourceDistribution(Opt_lo%n_g, Opt_lo%n_red, Phi_lo, XS_nF_Red)
            
            ! Error computation 
            k_lo=k_lo_old*dot_product(FSD, FSD)/dot_product(FSD, FSD_old)
            err1=abs(k_lo-k_lo_old)/k_lo_old
            err2=sqrt(dot_product(FSD-FSD_old, FSD-FSD_old)/dot_product(FSD, FSD))
            
            
            ! Log messages 
            write(output_unit,'(A)') '------------------------------------------------------------------------------------------'
            write(output_unit,*)
            write(output_unit,'(I6.1, A, F16.12, A, F16.9, A, F16.9)') iter_out_lo, '    k =  ',  k_lo, '    k_err = ',  err1, '    S_err = ',  err2
            write(output_unit,*)
            write(files%log,'(A)') '------------------------------------------------------------------------------------------'
            write(files%log,*)
            write(files%log,'(A, I6.1, A, F16.9, A)') '    At the cycle #',iter_out_lo,', of the FVM iterations, the k value is ', k_lo, '[/]'
            write(files%log,*)
            ! k value recording 
            write(files%convergence,'(I8.1,3ES21.13)') iter_out_lo, k_lo, err1, err2
            
            
            call Phi_extender(IDs_Red_Arr, Opt_lo%n_tot, Opt_lo%n_red, Opt_lo%n_g, Phi_lo, Phi_ext)

            !write(17, '(100000000ES21.13)') Phi_ext/sum(Phi_ext)
            !write(18, '(100000000ES21.13)') Phi_ext/sum(Phi_ext)

            if (err1<Opt_lo%tol1 .AND. err2<Opt_lo%tol2) exit
            
            call Local_Currents_TriIso(Opt_lo, LO_Mesh, LO_Param, XS_lo%Dif, Phi_ext)
            
            
            
            ! Accelerations 
            select case (flag_accel)
                case (0) ! FVM 
                    cycle
                    
                case (1) ! CMFD 
                    call CMFD_Header(Elem_CMFD, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo)

                case (2) ! pCMFD 
                    call pCMFD_Header(Elem_CMFD, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo)

                case (3) ! lpCMFD 
                    write(output_unit,'(A)') " ERROR! ACCELERATION NOT IMPLEMENTED YET!"
                    write(error_unit,'(A)') " ERROR! ACCELERATION NOT IMPLEMENTED YET!"
                    
                case (4) ! HG-CMFD 
                    call HGCMFD_Header(Elem_CMFD, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo)


                case (5) ! HG-CMFD 
                    call pHGCMFD_Header(Elem_CMFD, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, XS_lo, XS_gl, Serv_Vect, Serv_Matr, Phi_lo)


                case default
                    write(output_unit,'(A)') " ERROR! WRONG ACCELERATION SELECTION!"
                    write(error_unit,'(A)') " ERROR! WRONG ACCELERATION SELECTION!"
            end select
            
            !if (iter_out_lo == 1) then
            !    call Phi_extender(IDs_Red_Arr, Opt_lo%n_tot, Opt_lo%n_red, Opt_lo%n_g, Phi_lo, Phi_ext)
            !    write(17, '(100000000ES21.13)') Phi_ext/sum(Phi_ext)
            !    call Phi_extender(IDs_Red_Arr, Opt_lo%n_tot, Opt_lo%n_red, Opt_lo%n_g, Phi_lo_old, Phi_old_ext)
            !    write(18, '(100000000ES21.13)') sqrt(dot_product(Phi_ext/sum(Phi_ext) - Phi_old_ext/sum(Phi_old_ext), Phi_ext/sum(Phi_ext) - Phi_old_ext/sum(Phi_old_ext)))
            !elseif(err1<Opt_lo%tol1 .AND. err2<Opt_lo%tol2) then
            !    call Phi_extender(IDs_Red_Arr, Opt_lo%n_tot, Opt_lo%n_red, Opt_lo%n_g, Phi_lo, Phi_ext)
            !    write(17, '(100000000ES21.13)') Phi_ext/sum(Phi_ext)
            !    call Phi_extender(IDs_Red_Arr, Opt_lo%n_tot, Opt_lo%n_red, Opt_lo%n_g, Phi_lo_old, Phi_old_ext)
            !    write(18, '(100000000ES21.13)') sqrt(dot_product(Phi_ext/sum(Phi_ext) - Phi_old_ext/sum(Phi_old_ext), Phi_ext/sum(Phi_ext) - Phi_old_ext/sum(Phi_old_ext)))
            !    exit
            !end if            
            
            
        end do
        
        call cpu_time(t_fin)
        
        
        ! log.out file 
        write(files%log,'(A)') "=========================================================================================="
        write(files%log,*) 
        write(files%log,*) 
        write(files%log,*) 
        write(files%log,'(A)') "...Done. The convergence has been reached. Check the output files."
        write(files%log,*) 
        write(files%log,'(A, I6.1, A, F16.9)') "The final multiplication value, after ", iter_out_lo," iterations is k = ", k_lo
        write(files%log,"(A,ES21.13,A)") "given a total computation time of ", t_fin-t_ini," [s]"
        write(files%log,*) 
        write(files%log,*) 
        write(files%log,*)
        
        ! Display 
        write(output_unit,'(A)') "=========================================================================================="
        write(output_unit,*) 
        write(output_unit,*) 
        write(output_unit,*) 
        write(output_unit,'(A)') "...Done. The convergence has been reached. Check the output files."
        write(output_unit,*) 
        write(output_unit,'(A, I6.1, A, F16.9)') "The final multiplication value, after ", iter_out_lo," iterations is k = ", k_lo
        write(output_unit,"(A,ES21.13,A)") "given a total computation time of ", t_fin-t_ini," [s]"
        write(output_unit,*) 
        write(output_unit,*) 
        write(output_unit,*)
        
    end subroutine FVM_TriIso_Solver
    
    
    
    ! Steps for the FVM
    
    subroutine FVM_TriIso_MatrixBuilder(Opt_lo, LO_Mesh, LO_Param, XS_RM, XS_D, a_MM, b_MM, c_MM, d_MM, e_MM) 
        ! ===========================================================================================================
        ! --> Inputs:
        !
        ! n_x = # of meshes (w/ ghosts) on the x axis
        ! n_y = # of meshes (w/ ghosts) on the y axis
        ! n_g = # of energy groups
        ! n_tot = # of total meshes (w/ ghosts)
        ! IDs_lo = Array of meshes IDs
        ! IDs_neigh_lo = Array of neighbouring meshes IDs (the element (l,:) contains all the neighbours - 3 or 4)
        ! labels_lo = Array of labels for in-mesh (0), 1 BC (1), 2 BC (2), ..., n BC (n), ghost (-1) 
        ! Gam_lo = Array of gammas for the BC
        ! V_lo = Array of volumes
        ! A_lo = Array of surfaces for the currents (from BL/B, counter-clockwise)
        ! dl = Array of dx for the interface coefficients calculation (from BL/B, counter-clockwise)
        ! XS_RM = Removal XS Data  (w/ ghosts)
        ! XS_D = Diffusion Coeff. Data  (w/ ghosts)
        !       
        !
        !       
        ! --> Outputs:
        !
        ! b_MM =   main MM diagonal
        ! a_MM =  first MM  left off-diagonal
        ! c_MM =  first MM right off-diagonal 
        ! d_MM = second MM  left off-diagonal
        ! e_MM = second MM right off-diagonal
        ! IDs_lo_Red = Flag array to go back to the mesh grid w/ ghost meshes (ID_local_Reduced)
	    ! ===========================================================================================================
        		
        real(prec), intent(in) :: XS_RM(:), XS_D(:)
        real(prec), allocatable, intent(out) :: a_MM(:), b_MM(:), c_MM(:), d_MM(:), e_MM(:)
        
        type(Options_Data), intent(in) :: Opt_lo
        type(LO_geom), allocatable, intent(in) :: LO_Mesh(:)
        type(LO_coeff), allocatable, intent(in) :: LO_Param(:)
        
        integer :: i, j, l, g, g_tot, g_red, u, n, iln, face_ID  ! In-Line Neighbour calculation
        real(prec) :: bc, D_coeff, temp ! Dummy variable for the Boundary Condition contribution
        real(prec), allocatable :: a_ToRed(:), b_ToRed(:), c_ToRed(:), d_ToRed(:), e_ToRed(:)
        
        
        write(output_unit,*) 
        write(files%log,*) 
        write(output_unit,'(A)') 'Building the Migration Matrix...'
        write(files%log,'(A)') 'Building the Migration Matrix...'
                
        allocate(a_MM(Opt_lo%n_g*OPt_lo%n_red))
        allocate(b_MM(Opt_lo%n_g*OPt_lo%n_red))
        allocate(c_MM(Opt_lo%n_g*OPt_lo%n_red))
        allocate(d_MM(Opt_lo%n_g*OPt_lo%n_red))
        allocate(e_MM(Opt_lo%n_g*OPt_lo%n_red))
        a_MM=0.e0_prec
        b_MM=0.e0_prec
        c_MM=0.e0_prec
        d_MM=0.e0_prec
        e_MM=0.e0_prec
        
        n=0
        ! iln = neighbour mesh ID
        ! l = ghost-included mesh grid ID
        ! n = ghost-excluded mesh grid ID
        
		do l=1, Opt_lo%n_tot
            !i=mod(l-1,n_x)+1
            !j=int((l-i)/n_x)+1
            
            
            if (LO_Mesh(l)%ID_Red>0) then    ! Ghost filter
                n=n+1
				
                do g=1, Opt_lo%n_g    ! Energy groups loop
                    g_tot=(g-1)*Opt_lo%n_tot
                    g_red=(g-1)*Opt_lo%n_red
                    b_MM(g_red+n)=XS_RM(g_tot+l)    ! Main diagonal always present
                    
                    ! Roll on the faces 
                    do u=1,size(LO_Mesh(l)%Neigh_ID)
                        iln=LO_Mesh(l)%Neigh_ID(u)
                        if (LO_Mesh(iln)%ID_Red>0) then
                            
                            face_ID=LO_Mesh(l)%Sides_Neigh(u)
                            D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), LO_Mesh(iln)%dl_lo(face_ID), LO_Mesh(l)%dl_lo(u))
                            
                            select case(LO_Mesh(l)%Sides_ID(u))
                                case(1) 
                                    d_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! B
                                    b_MM(g_red+n)=b_MM(g_red+n)-d_MM(g_red+n)
                                case(2) 
                                    c_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! R
                                    b_MM(g_red+n)=b_MM(g_red+n)-c_MM(g_red+n)
                                case(3) 
                                    e_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! T
                                    b_MM(g_red+n)=b_MM(g_red+n)-e_MM(g_red+n)
                                case(4) 
                                    a_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! L
                                    b_MM(g_red+n)=b_MM(g_red+n)-a_MM(g_red+n)
                                case(5) 
                                    e_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! TR
                                    b_MM(g_red+n)=b_MM(g_red+n)-e_MM(g_red+n)
                                case(6) 
                                    e_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! TL
                                    b_MM(g_red+n)=b_MM(g_red+n)-e_MM(g_red+n)
                                case(7)
                                    d_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! BL
                                    b_MM(g_red+n)=b_MM(g_red+n)-d_MM(g_red+n)
                                case(8)
                                    d_MM(g_red+n)=-D_coeff*LO_Mesh(l)%A_lo(u)/LO_Mesh(l)%V_lo ! BR
                                    b_MM(g_red+n)=b_MM(g_red+n)-d_MM(g_red+n)
                            end select
                        elseif (LO_Mesh(l)%BC(u)==1 .AND. LO_Param(g_tot+l)%Gam(u)>=0.e0_prec) then
                            bc = LO_Param(g_tot+l)%Gam(u) / (1.e0_prec+LO_Param(g_tot+l)%Gam(u)*LO_Mesh(l)%dl_lo(u)/XS_D(g_tot+l)) * LO_Mesh(l)%A_lo(u) / LO_Mesh(l)%V_lo
                            b_MM(g_red+n) = b_MM(g_red+n) + bc
                        end if
                        
                    end do
                    
                    
                    
                    
                    !!Old Approach 
					!select case (ori(l))    ! Selection of the mesh type 
                    !    case('BL') ! Matrix entries (in-mesh & BC) 
					!		! Cannot be put in a loop
					!	    
                    !        
                    !        
                    !        
                    !        
					!		! Bottom Neighbour
					!		iln=IDs_neigh_lo(l, 1)
					!		if  (labels_lo(iln) .NE. -1) then 
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,2), dl(l,1))
					!			d_MM(g_red+n)=-D_coeff * A_lo(l,1)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - d_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,1)/XS_D(g_tot+l)) * A_lo(l,1)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Top-Right Neighbour
					!		iln=IDs_neigh_lo(l, 2)
					!		if (labels_lo(iln) .NE. -1) then
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,1), dl(l,2))
					!			e_MM(g_red+n)=-D_coeff * A_lo(l,2)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - e_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,2)/XS_D(g_tot+l)) * A_lo(l,2) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Left Neighbour
					!		iln=IDs_neigh_lo(l, 3)
					!		if (labels_lo(iln) .NE. -1) then
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,2), dl(l,3))
					!			a_MM(g_red+n)=-D_coeff * A_lo(l,3)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - a_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,3)/XS_D(g_tot+l)) * A_lo(l,3) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!	case('BR') ! Matrix entries (in-mesh & BC) 
					!		! Cannot be put in a loop
					!	
					!		! Bottom Neighbour
					!		iln=IDs_neigh_lo(l, 1)
					!		if  (labels_lo(iln) .NE. -1) then 
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,3), dl(l,1))
					!			d_MM(g_red+n)=-D_coeff * A_lo(l,1)/V_lo(l) ! B
					!			b_MM(g_red+n)=b_MM(g_red+n) - d_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,1)/XS_D(g_tot+l)) * A_lo(l,1)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Right Neighbour
					!		iln=IDs_neigh_lo(l, 2)
					!		if (labels_lo(iln) .NE. -1) then
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,3), dl(l,2))
					!			c_MM(g_red+n)=-D_coeff * A_lo(l,2)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - c_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,2)/XS_D(g_tot+l)) * A_lo(l,2) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Top-Left Neighbour
					!		iln=IDs_neigh_lo(l, 3)
					!		if (labels_lo(iln) .NE. -1) then
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,1), dl(l,3))
					!			e_MM(g_red+n)=-D_coeff * A_lo(l,3)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - e_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,3)/XS_D(g_tot+l)) * A_lo(l,3) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!	case('TL') ! Matrix entries (in-mesh & BC) 
					!		! Cannot be put in a loop
					!		
					!		! Bottom-Right
					!		iln=IDs_neigh_lo(l, 1)
					!		if (labels_lo(iln) .NE. -1) then ! Right neighbouring mesh
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,3), dl(l,1))
					!			d_MM(g_red+n)=-D_coeff * A_lo(l,1)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - d_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,1)/XS_D(g_tot+l)) * A_lo(l,2) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Top
					!		iln=IDs_neigh_lo(l, 2)
					!		if (labels_lo(iln) .NE. -1) then ! Right neighbouring mesh
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,1), dl(l,2))
					!			e_MM(g_red+n)=-D_coeff * A_lo(l,2)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - e_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,2)/XS_D(g_tot+l)) * A_lo(l,2) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Left
					!		iln=IDs_neigh_lo(l, 3)
					!		if (labels_lo(iln) .NE. -1) then ! Right neighbouring mesh
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,2), dl(l,3))
					!			a_MM(g_red+n)=-D_coeff * A_lo(l,3)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - a_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,3)/XS_D(g_tot+l)) * A_lo(l,3) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!	case('TR') ! Matrix entries (in-mesh & BC) 
					!		! Cannot be put in a loop
					!		
					!		! Bottom-Left
					!		iln=IDs_neigh_lo(l, 1)
					!		if (labels_lo(iln) .NE. -1) then ! Right neighbouring mesh
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,2), dl(l,1))
					!			d_MM(g_red+n)=-D_coeff * A_lo(l,1)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - d_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,1)/XS_D(g_tot+l)) * A_lo(l,1) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Right
					!		iln=IDs_neigh_lo(l, 2)
					!		if (labels_lo(iln) .NE. -1) then ! Right neighbouring mesh
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,3), dl(l,2))
					!			c_MM(g_red+n)=-D_coeff * A_lo(l,2)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - c_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,2)/XS_D(g_tot+l)) * A_lo(l,2) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
					!		
					!		! Top
					!		iln=IDs_neigh_lo(l, 3)
					!		if (labels_lo(iln) .NE. -1) then ! Right neighbouring mesh
					!			D_coeff=Int_Coeff_D(XS_D(g_tot+iln), XS_D(g_tot+l), dl(iln,1), dl(l,3))
					!			e_MM(g_red+n)=-D_coeff * A_lo(l,3)/V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) - e_MM(g_red+n)
                    !            
					!		else
					!			bc = Gam_lo(g_tot+l)/(1.e0_prec+Gam_lo(g_tot+l)*dl(l,3)/XS_D(g_tot+l)) * A_lo(l,3) / V_lo(l)
					!			b_MM(g_red+n)=b_MM(g_red+n) + bc
                    !            
					!		end if
                    !        
					!	case default 
					!		write(output_unit,'(A)') "ERROR IN THE MIGRATION MATRIX BUILDING - Wrong orientation"
					!		write(error_unit,'(A)') "ERROR IN THE MIGRATION MATRIX BUILDING - Wrong orientation"
					!
					!end select
					
                end do
            end if
        end do
        
        ! Diagonal dominance check 
        do l=1, size(a_MM)
            temp=abs(b_MM(l))-abs(a_MM(l))-abs(c_MM(l))-abs(d_MM(l))-abs(e_MM(l))
            if (temp<0) then
				write(output_unit,'(A)') "ERROR IN THE MIGRATION MATRIX BUILDING - Not diagonally dominant"
				write(error_unit,'(A)') "ERROR IN THE MIGRATION MATRIX BUILDING - Not diagonally dominant"
            end if
        end do
        
        write(output_unit,'(A)') '...Done'
        write(files%log,'(A)') '...Done'
        write(output_unit,*) 
        write(files%log,*) 
        
        
    end subroutine FVM_TriIso_MatrixBuilder
    
    pure function Int_Coeff_D(D_loc_1, D_loc_2, dl_1, dl_2) result (D_int) 
        real(prec), intent(in) :: D_loc_1, D_loc_2, dl_1, dl_2
        real(prec) :: D_int
        
        D_int=1.e0_prec/(dl_1/D_loc_1+dl_2/D_loc_2)
    end function Int_Coeff_D
        
    pure function Phi_Sup(D_loc_1, D_loc_2, dl_1, dl_2, Phi_1, Phi_2) result (Phi_s) 
        real(prec), intent(in) :: D_loc_1, D_loc_2, dl_1, dl_2, Phi_1, Phi_2
        real(prec) :: Phi_s
        
        Phi_s=(D_loc_1/dl_1*Phi_1+D_loc_2/dl_2*Phi_2)/(D_loc_1/dl_1+D_loc_2/dl_2)
    end function Phi_Sup
    
    pure function FissionSourceDistribution(n_g, n, Phi, XS_nF) result(FSD) 
        integer, intent(in) :: n_g, n
        real(prec), intent(in) :: Phi(:), XS_nF(:)
        real(prec) :: FSD(n)
        
        integer :: t
        
        FSD=0.e0_prec
        do t=1, n_g
            FSD(:)=FSD(:)+Phi((t-1)*n+1:t*n)*XS_nF((t-1)*n+1:t*n)
        end do
        
    end function FissionSourceDistribution
    
    pure function Total_Source_SingleGroup(t, n_g, n_red, k, Chi, FSD, XS_SC, Phi) result(S_tot) 
        integer, intent(in) :: t, n_g, n_red 
        real(prec), intent(in) :: k, Chi(:), FSD(:), XS_SC(:,:), Phi(:)
        
        integer :: g, t_red, g_red
        real(prec) :: S_tot(n_red)
        
        S_tot=0.e0_prec
        t_red=(t-1)*n_red
        
        S_tot(:)=S_tot(:)+Chi(t_red+1:t_red+n_red)/k*FSD(:)
        do g=1, n_g
            if (g.NE.t) then
                g_red=(g-1)*n_red
                S_tot(:)=S_tot(:)+XS_SC(g_red+1:g_red+n_red,t)*Phi(g_red+1:g_red+n_red)
            end if
        end do
    end function Total_Source_SingleGroup
    
    pure function Total_Source(n_g, n_red, k, Chi, FSD, XS_SC, Phi) result(S_tot) 
        integer, intent(in) :: n_g, n_red 
        real(prec), intent(in) :: k, Chi(:), FSD(:), XS_SC(:,:), Phi(:)
        
        integer :: t, g, t_red, g_red
        real(prec) :: S_tot(n_g*n_red)
        
        S_tot=0.e0_prec
        do t=1, n_g
            t_red=(t-1)*n_red
            S_tot(t_red+1:t_red+n_red)=S_tot(t_red+1:t_red+n_red)+Chi(t_red+1:t_red+n_red)/k*FSD(:)
            do g=1, n_g
                if (g.NE.t) then
                    g_red=(g-1)*n_red
                    S_tot(t_red+1:t_red+n_red)=S_tot(t_red+1:t_red+n_red)+XS_SC(g_red+1:g_red+n_red,t)*Phi(g_red+1:g_red+n_red)
                end if
            end do
        end do
    end function Total_Source
    
    !subroutine MM_Printing(n_g, n_x, n, a, b, c, d, e) 
    !    real(prec), intent(in) :: a(:), b(:), c(:), d(:), e(:)
    !    real(prec) :: MM(n_g*n, n_g*n)
    !    integer :: g, g_red, i
    !    integer, intent(in) :: n_g, n, n_x
    !    
    !    ! Explicit matrix building
    !    MM=0.e0_prec
    !    
    !    do g=1,n_g
    !        g_red=(g-1)*n
    !        do i=1, n
    !            MM(g_red+i, g_red+i)=b(g_red+i)
    !            if (i>1) MM(g_red+i, g_red+i-1)=a(g_red+i)
    !            if (i<n) MM(g_red+i, g_red+i+1)=c(g_red+i)
    !            if (i>n_x) MM(g_red+n, g_red+n-n_x)=d(g_red+i)
    !            if (i<n_x) MM(g_red+n, g_red+n+n_x)=e(g_red+i)
    !        end do
    !    end do
    !    
    !    do i=1, n*n_g
    !        write(16, '(1000000ES21.13)') MM(i,:)
    !    end do
    !    write(16, *) 
    !    write(16, *) 
    !    write(16, *) 
    !    
    !end subroutine MM_Printing    
    
    subroutine Phi_Extender(IDs_lo_Red, n_tot_lo, n_red, n_g_lo, Phi_red, Phi_ext) 
        integer, intent(in) :: IDs_lo_Red(:), n_tot_lo, n_red, n_g_lo
        real(prec), intent(in) :: Phi_red(:)
        
        integer :: l, t, t_tot, t_red
        real(prec), intent(out) :: Phi_ext(:)
        
        Phi_ext=0.e0_prec
        do l=1, n_tot_lo
            do t=1, n_g_lo
                t_tot=(t-1)*n_tot_lo
                t_red=(t-1)*n_red
                if (IDs_lo_Red(l)>0) Phi_ext(t_tot+l)=Phi_red(t_red+IDs_lo_Red(l))
            end do
        end do
        
    end subroutine Phi_Extender

    subroutine Local_Currents_TriIso(Opt_lo, LO_Mesh, LO_Param, XS_D, Phi_lo) 
        real(prec), intent(in) :: XS_D(:), Phi_lo(:)
        
        type(Options_Data), intent(in) :: Opt_lo
        type(LO_geom), intent(in):: LO_Mesh(:)
        type(LO_coeff), intent(inout) :: LO_Param(:)
        
        integer :: l, t, u, n, t_tot, iln_tot, face_ID
        real(prec) :: D_tilde, Phi_s
        
        n=0
        do l=1, Opt_lo%n_tot
            if (LO_Mesh(l)%ID_Red>0) then 
                n=n+1
                
                ! Roll on the E groups
                do t=1, Opt_lo%n_g
                    t_tot=(t-1)*Opt_lo%n_tot
                    LO_Param(t_tot+l)%J_net=0.e0_prec
                    LO_Param(t_tot+l)%J_part=0.e0_prec
                    
                    ! Roll on the sides
                    do u=1, size(LO_Mesh(l)%Neigh_ID)
                        iln_tot=LO_Mesh(l)%Neigh_ID(u)
                        !iln_red=LO_Mesh(l)%J_lo_gl_ID_Red(u)
                        
                        if (LO_Mesh(iln_tot)%ID_Red>0) then
                            face_ID=LO_Mesh(l)%Sides_Neigh(u)
                            D_tilde=Int_Coeff_D(XS_D(t_tot+iln_tot), XS_D(t_tot+l), LO_Mesh(iln_tot)%dl_lo(face_ID), LO_Mesh(l)%dl_lo(u))
                            
                            Phi_s=Phi_Sup(XS_D(t_tot+iln_tot), XS_D(t_tot+l), LO_Mesh(iln_tot)%dl_lo(face_ID), LO_Mesh(l)%dl_lo(u), Phi_lo(t_tot+iln_tot), Phi_lo(t_tot+l))
                            
                            LO_Param(t_tot+l)%D_til(u)=D_tilde
                            
                            ! Net
                            LO_Param(t_tot+l)%J_net(u)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
                            
                            ! Partial
                            LO_Param(t_tot+l)%J_part(u,1)=Phi_s/4.e0_prec-D_tilde/2.e0_prec*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
                            LO_Param(t_tot+l)%J_part(u,2)=Phi_s/4.e0_prec+D_tilde/2.e0_prec*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
                        
                        else
                            LO_Param(t_tot+l)%D_til(u)=XS_D(t_tot+l)
                            
                            ! Net
                            LO_Param(t_tot+l)%J_net(u)=LO_Param(t_tot+l)%Gam(u) / (1.e0_prec+LO_Param(t_tot+l)%Gam(u)*LO_Mesh(l)%dl_lo(u)/XS_D(t_tot+l)) * Phi_lo(t_tot+l)
                            
                            ! Partial
                            !Phi_s=(XS_D(t_tot+l)/LO_Mesh(l)%dl_lo(u)/2.e0_prec) / (XS_D(t_tot+l)/LO_Mesh(l)%dl_lo(u)/2.e0_prec+(LO_Param(t_tot+l)%Gam(u)-1/4.e0_prec)) * Phi_lo(t_tot+l)
                            !LO_Param(t_tot+l)%J_part(u,1)=Phi_s/4.e0_prec-D_tilde/2.e0_prec*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
                            !
                            !Phi_s=(XS_D(t_tot+l)/LO_Mesh(l)%dl_lo(u)/2.e0_prec) / (XS_D(t_tot+l)/LO_Mesh(l)%dl_lo(u)/2.e0_prec+(LO_Param(t_tot+l)%Gam(u)+1/4.e0_prec)) * Phi_lo(t_tot+l)
                            !LO_Param(t_tot+l)%J_part(u,2)=Phi_s/4.e0_prec+D_tilde/2.e0_prec*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
                            
                            Phi_s=(XS_D(t_tot+l)/LO_Mesh(l)%dl_lo(u)) / (LO_Param(t_tot+l)%Gam(u)+XS_D(t_tot+l)/LO_Mesh(l)%dl_lo(u)) * Phi_lo(t_tot+l)
                            LO_Param(t_tot+l)%J_part(u,1)=(0.5e0_prec+LO_Param(t_tot+l)%Gam(u)) * Phi_s/2.e0_prec
                            LO_Param(t_tot+l)%J_part(u,2)=(0.5e0_prec-LO_Param(t_tot+l)%Gam(u)) * Phi_s/2.e0_prec
                        end if
                    end do
                    
                end do
            end if
        end do
        
        
        
        
        
        
        
        
        
        
        
        !J_lo=0.e0_prec
        !do l=1, n_tot_lo
        !    if (J_map_lo(l)>-1) then
        !        do t=1, n_g_lo
        !            t_tot=(t-1)*n_tot_lo
        !            
        !            select case (ori(l))
        !                case('BL') 
        !                    ! Bottom Neighbour
        !                    iln_tot=IDs_neigh_lo(l, 1) ! B neighbouring mesh index
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,2), dl_lo(l,1))
        !                        J_lo(t_tot+l,1)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,1)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,1)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Top-Right Neighbour
        !                    iln_tot=IDs_neigh_lo(l, 2) 
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,1), dl_lo(l,2))
        !                        J_lo(t_tot+l,2)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,2)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,2)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Left Neighbour
        !                    iln_tot=IDs_neigh_lo(l, 3) 
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,2), dl_lo(l,3))
        !                        J_lo(t_tot+l,3)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,3)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,3)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                case('BR') 
        !                    ! Bottom Neighbour
        !                    iln_tot=IDs_neigh_lo(l, 1) 
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,3), dl_lo(l,1))
        !                        J_lo(t_tot+l,1)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,1)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,1)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Right Neighbour
        !                    iln_tot=IDs_neigh_lo(l, 2) ! R neighbouring mesh index
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,3), dl_lo(l,2))
        !                        J_lo(t_tot+l,2)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,2)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,2)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Top-Left Neighbour
        !                    iln_tot=IDs_neigh_lo(l, 3) ! L neighbouring mesh index
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,1), dl_lo(l,3))
        !                        J_lo(t_tot+l,3)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,3)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,3)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                case('TL') 
        !                    ! Bottom-Right Neighbour 
        !                    iln_tot=IDs_neigh_lo(l, 1) 
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,3), dl_lo(l,1))
        !                        J_lo(t_tot+l,1)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,1)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,1)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Right Neighbour 
        !                    iln_tot=IDs_neigh_lo(l, 2)
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,1), dl_lo(l,2))
        !                        J_lo(t_tot+l,2)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,2)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,2)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Top Neighbour 
        !                    iln_tot=IDs_neigh_lo(l, 3)
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,2), dl_lo(l,3))
        !                        J_lo(t_tot+l,3)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,3)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,3)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                case('TR') 
        !                    ! Bottom-Left Neighbour 
        !                    iln_tot=IDs_neigh_lo(l, 1) 
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,2), dl_lo(l,1))
        !                        J_lo(t_tot+l,1)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,1)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,1)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Right Neighbour 
        !                    iln_tot=IDs_neigh_lo(l, 2)
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,3), dl_lo(l,2))
        !                        J_lo(t_tot+l,2)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,2)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,2)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !                    ! Top Neighbour 
        !                    iln_tot=IDs_neigh_lo(l, 3)
        !                    if (J_map_lo(iln_tot)>-1) then
        !                        D_tilde=Int_Coeff_D(XS_D_lo(t_tot+iln_tot), XS_D_lo(t_tot+l), dl_lo(iln_tot,1), dl_lo(l,3))
        !                        J_lo(t_tot+l,3)=-D_tilde*(Phi_lo(t_tot+iln_tot)-Phi_lo(t_tot+l))
        !                    else
        !                        J_lo(t_tot+l,3)=Gam_lo(t_tot+l)/(1.e0_prec+Gam_lo(t_tot+l)*dl_lo(l,3)/XS_D_lo(t_tot+l))*Phi_lo(t_tot+l)
        !                    end if
        !                    
        !            end select
        !        end do
        !    end if
        !end do
    
    end subroutine  Local_Currents_TriIso
    
end module FVM_TriIso_Alg
    
