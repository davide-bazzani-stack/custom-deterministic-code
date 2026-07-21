module Service_Fcns 
    use Variables
    use IO_module, only : files
    use iso_fortran_env, only : output_unit, error_unit
        
        implicit none
    
    
    contains
    
    
        subroutine Flux_Vol_Normalization(n_g, V, Phi) 
            integer, intent(in) :: n_g
            real(dp), intent(in) :: V(:)
            real(dp), intent(inout) :: Phi(:)

            integer :: i, j, g, g_red, n
            real(dp) :: tmp

            n=size(V)
            tmp=0.d0

            do g=1,n_g
                g_red=(g-1)*n
                tmp=tmp+sum(Phi(g_red+1:g_red+n)*V)
            end do
            tmp=tmp/sum(V) !dble(n_mesh*n_g) ! Average flux counting all the groups

            Phi=Phi*1.d0/tmp
        end subroutine Flux_Vol_Normalization

        subroutine Flux_lo_Vol_Normalization(n_tot, n_red, n_g, labels, V, Phi) 
            integer, intent(in) :: n_tot, n_g, n_red, labels(:)
            real(dp), intent(in) :: V(:)
            real(dp), intent(inout) :: Phi(:)

            integer :: i, l, g, g_red, l_red
            real(dp) :: tmp

            tmp=0.d0
            do g=1,n_g
                l_red=0
                do l=1, n_tot
                    if (labels(l)>-1) then
                        l_red=l_red+1
                        g_red=(g-1)*n_red
                        tmp=tmp+Phi(g_red+l_red)*V(l)
                    end if
                end do
            end do
            tmp=tmp/sum(V) !dble(n_mesh*n_g) ! Average flux counting all the groups

            Phi=Phi*1.d0/tmp
        end subroutine Flux_lo_Vol_Normalization

        subroutine Output_Converter(n_g, n_x, Mesh_Position, Flux_in, Flux_out) 
            integer, intent(in) :: n_g, n_x, Mesh_Position(:)
            real(dp), intent(in) :: Flux_in(:) 

            real(dp), allocatable, intent(out) :: Flux_out(:)


            integer :: l, g, g_red, g_tot, n, n_tot, n_red

            n_tot=size(Mesh_Position)
            n_red=count(Mesh_Position>0)
            allocate(Flux_out(n_tot*n_g))

            Flux_out=0.d0
            do g=1, n_g
                g_tot=(g-1)*n_tot
                g_red=(g-1)*n_red
                do l=1, n_tot
                    if (Mesh_Position(l)>0) Flux_out(g_tot+l)=Flux_in(g_red+Mesh_Position(l))
                end do
            end do

            ! To add additional converting variablels here
        end subroutine Output_Converter

        subroutine Results_Writing(n_g, n_x, n_y, n_tot, univ, Flux, kFi_XS, V_lo) 

            use IO_module, only : files

            integer, intent(in) :: n_g, n_x, n_y, n_tot, univ(:)
            real(dp), intent(in) :: Flux(:), kFi_XS(:), V_lo(:)

            integer :: g, j, l, g_tot
            real(dp), allocatable :: Power(:)

            ! - Universes 
            do l=1, n_tot
                write(files%universes, '(I6.1)') univ(l)
                !write(files%universes, '(100000I6.1)') univ
            end do
            write(files%universes,*) 
            write(files%universes,*) 
            write(files%universes,*) 

            ! --> Flux file 
            ! - Matrix 
            write(files%phi, '(A)') "Matrix flux"
            do g=1, n_g
                g_tot=(g-1)*n_tot
                write(files%phi, '(A, I1.0)') 'Group g = ', g
                do j=1, n_y
                    write(files%phi, '(100000ES21.13)') Flux(g_tot+(j-1)*n_x+1:g_tot+j*n_x)
                end do
                write(files%phi,*) 
                write(files%phi,*) 
                write(files%phi,*) 
            end do

            ! - Single array 
            write(files%eigenfunctions, '(A)') "Dominant eigenfunction"
            write(files%eigenfunctions, *)
            do g=1, n_g
                g_tot=(g-1)*n_tot
                !write(3, '(A, I1.0)') 'Group g = ', g
                do l=1, n_tot
                    write(files%eigenfunctions, '(100000ES21.13)') Flux(g_tot+l)
                end do
                write(files%eigenfunctions,*) 
                write(files%eigenfunctions,*) 
                write(files%eigenfunctions,*) 
            end do

            ! --> Power file 
            allocate(Power(n_x*n_y))
            Power=0.d0
            do g=1, n_g
                g_tot=(g-1)*n_tot
                Power=Power+Flux(g_tot+1:g_tot+n_tot)*kFi_XS(g_tot+1:g_tot+n_tot)*V_lo
            end do

            write(files%power, '(A)') "Generated Power Profile"
            do j=1, n_y
                write(files%power, '(100000ES21.13)') Power((j-1)*n_x+1:j*n_x)
            end do
            write(files%power,*) 
            write(files%power,*) 
            write(files%power,*) 
            write(files%power,*) 
            write(files%power,*) 

            Power=Power/(sum(Power)/count(Power>0.d0))
            write(files%power, '(A)') "Relative Power Profile"
            do j=1, n_y
                write(files%power, '(100000ES21.13)') Power((j-1)*n_x+1:j*n_x)
            end do
            write(files%power,*) 
            write(files%power,*) 
            write(files%power,*) 
            write(files%power,*) 
            write(files%power,*) 
        
        
        end subroutine Results_Writing
    
        subroutine dl_gl_Auto(vertices, point, dist) 
            ! Works only for convex poligons
            ! The order of the vertices must be counter-clockwise

            integer :: i, n
            real(dp), intent(in) :: vertices(:,:), point(:)
            real(dp), allocatable :: vert(:,:)
            real(dp), intent(out) :: dist(:)

            n=size(vertices,2)

            allocate(vert(size(vertices,1), n+1))


            vert(:,1:n)=vertices(:,:)
            vert(:,n+1)=vertices(:,1)

            dist(:)=abs((vert(1,2:n+1)-vert(1,1:n))*(vert(2,1:n)-point(2))-(vert(1,1:n)-point(1))*(vert(2,2:n+1)-vert(2,1:n)))/sqrt((vert(1,2:n+1)-vert(1,1:n))**2+(vert(2,2:n+1)-vert(2,1:n))**2)

        end subroutine dl_gl_Auto
        
        subroutine dl_gl_Calc(vertices, point, distance) 
            ! Works only for convex poligons
            ! The order of the vertices must be counter-clockwise

            integer :: i, n
            real(dp), intent(in) :: vertices(:,:), point(:)
            real(dp), allocatable :: vert(:,:), dist(:)
            real(dp), intent(out) :: distance

            n=size(vertices,2)

            allocate(vert(size(vertices,1), n+1))
            allocate(dist(n))


            vert(:,1:n)=vertices(:,:)
            vert(:,n+1)=vertices(:,1)

            dist(:)=abs((vert(1,2:n+1)-vert(1,1:n))*(vert(2,1:n)-point(2))-(vert(1,1:n)-point(1))*(vert(2,2:n+1)-vert(2,1:n)))/sqrt((vert(1,2:n+1)-vert(1,1:n))**2+(vert(2,2:n+1)-vert(2,1:n))**2)

            distance=minval(dist)

            deallocate(vert)
        end subroutine dl_gl_Calc
        
        subroutine Point_in_Polygon_Test_2D(vertices, point, inside) 
            ! Works only for convex poligons
            ! The order of the vertices must be counter-clockwise
            
            integer :: i, n
            real(dp), intent(in) :: vertices(:,:), point(:)
            real(dp), allocatable :: vert(:,:), Signed_Area(:)
            logical, intent(out) :: inside
            
            n=size(vertices,2)
            inside = .FALSE.
            
            allocate(vert(size(vertices,1), n+1))
            allocate(Signed_Area(n))
            
            
            vert(:,1:n)=vertices(:,:)
            vert(:,n+1)=vertices(:,1)
            
            do i=1, n
                Signed_Area(i) = (vert(1,i+1)-vert(1,i))*(point(2)-vert(2,i))-(vert(2,i+1)-vert(2,i))*(point(1)-vert(1,i))
            end do
            
            if (abs(sum(abs(Signed_Area))-abs(sum(Signed_Area)))<1.0d-12) inside= .TRUE. 
            
            deallocate(Signed_Area)
            deallocate(vert)
            
        end subroutine Point_in_Polygon_Test_2D    
    
        subroutine Accel_Ini(flag_accel, Opt_lo, Opt_gl, Elem_det, Elem_CMFD, dx_gh, dy_gh, dl_input, dx_input, XS_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param, Serv_Vect, Serv_Matr)  

                ! Inputs
                integer, intent(in) :: flag_accel
                type(Options_Data), intent(in) :: Opt_lo
                type(Figure), allocatable, intent(in) :: Elem_det(:)
                type(LO_geom), allocatable, intent(inout) :: LO_Mesh(:)
                type(LO_coeff), allocatable, intent(inout) :: LO_Param(:)
                real(dp), intent(in) :: dx_gh, dy_gh, dl_input, dx_input

                ! Outputs
                type(Options_Data), intent(out) :: Opt_gl
                type(Figure), allocatable, intent(out) :: Elem_CMFD(:)
                type(GL_geom), allocatable, intent(out) :: GL_Mesh(:)
                type(GL_coeff), allocatable, intent(out) :: GL_Param(:)
                type(XS_Data), intent(out)  :: XS_gl
                !type() :: 

                type(Accel_Vars_Vect), intent(out) :: Serv_Vect
                type(Accel_Vars_Matr), intent(out) :: Serv_Matr






                integer :: flag_part, no_faces, g, g_tot
                real(dp) :: dl, vert_temp(2,200)
                logical :: test

                integer :: i, j, l, i_red, j_red, l_red, ii, jj, ll, t, t_tot, u, uu, n, m, nb, ios, n_elements, n_vert, box_mat_ID, lo_in_gl_meshes_x, lo_in_gl_meshes_y
                character(len=256) :: file_line, Elem_type, dummy

                ! Partial current option 
                if (flag_accel==2 .OR. flag_accel==5 ) then
                    flag_part=2
                else
                    flag_part=1
                end if

                ! File reader 
                if (flag_accel==1 .OR. flag_accel==2 .OR. flag_accel==3) then  
                    ! CMFD 
                    write(output_unit,*) 
                    write(files%log,*) 
                    write(output_unit,'(A)') 'SELECTED: CMFD Acceleration'
                    write(files%log,'(A)') 'SELECTED: CMFD Acceleration'
                    write(output_unit,'(A)') 'Reading the coarse mesh input parameters...'
                    write(files%log,'(A)') 'Reading the coarse mesh input parameters...'

                    ! Custom HGCMFD zones 
                    open(unit=130, file='CMFD_Grid.inp', status='old', action='read', iostat=ios)

                    ! Allocation and detection of the geometry details 
                    allocate(Elem_CMFD(Opt_gl%n_x*Opt_gl%n_y)) 
                    n=0
                    Elem_CMFD(:)%flag_BG=0
                    ios=0
                    do
                        read(130,'(A)', iostat=ios) file_line

                        ! In case of any error 
                        if (ios>0) then
                            write(output_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element reading routine'
                            write(error_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element reading routine'
                            exit
                        end if

                        ! End of file 
                        if (ios < 0) exit

                        ! Skip the comments and the blank lines 
                        if (trim(file_line)=='' .OR. file_line(1:1)=='#') cycle  

                        ! Figures 
                        if (index(file_line, '%') > 0) then 
                            ! Different cases 
                            n=n+1
                            Elem_CMFD(n)%ID=n
                            read(file_line(index(file_line, '%')+2:), '(A)') Elem_type

                            select case(trim(Elem_type)) 

                                    case ('Crfc') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,1))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        Elem_CMFD(n)%Fig_ID=1

                                    case ('Trng') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,3))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        Elem_CMFD(n)%Fig_ID=2

                                    case ('Rhomb') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,4))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        Elem_CMFD(n)%Fig_ID=3

                                    case ('Rect') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,4))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        Elem_CMFD(n)%Fig_ID=4

                                    case ('Pent') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,5))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        Elem_CMFD(n)%Fig_ID=5

                                    case ('Hex') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,6))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        Elem_CMFD(n)%Fig_ID=6

                                    case ('Ept') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,7))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        Elem_CMFD(n)%Fig_ID=7

                                    case ('Oct') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,8))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        Elem_CMFD(n)%Fig_ID=8

                                    case ('Enna') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,9))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        read(130, *) Elem_CMFD(n)%vertices(:,9)
                                        Elem_CMFD(n)%Fig_ID=9

                                    case ('Deca') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,10))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        read(130, *) Elem_CMFD(n)%vertices(:,9)
                                        read(130, *) Elem_CMFD(n)%vertices(:,10)
                                        Elem_CMFD(n)%Fig_ID=10

                                    case ('Custom') 
                                        !read 'n_vert' to automatically allocate!
                                        !allocate(Elem_CMFD(n)%centroid(2))
                                        !allocate(Elem_CMFD(n)%vertices(2,10))
                                        !read(130, *) Elem_CMFD(n)%centroid
                                        !read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,9)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,10)
                                        !Elem_CMFD(n)%Fig_ID=0
                                        write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                                        write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"

                                    case default 
                                        write(output_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element selection'
                                        write(error_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element selection'
                                        exit
                                end select
                        end if
                    end do   

                    ! Local ghost layer implementation 
                    do l=1, size(Elem_CMFD)
                        Elem_CMFD(l)%centroid(1)=Elem_CMFD(l)%centroid(1)+dx_gh 
                        Elem_CMFD(l)%centroid(2)=Elem_CMFD(l)%centroid(2)+dy_gh!*2
                        Elem_CMFD(l)%vertices(1,:)=Elem_CMFD(l)%vertices(1,:)+dx_gh
                        Elem_CMFD(l)%vertices(2,:)=Elem_CMFD(l)%vertices(2,:)+dy_gh!*2
                    end do

                    close(130)

                    write(output_unit,'(A)') '...Done'
                    write(files%log,'(A)') '...Done'
                    write(output_unit,*) 
                    write(files%log,*) 




                elseif (flag_accel==4 .OR. flag_accel==5 ) then 
                    ! HG-CMFD 
                    write(output_unit,*) 
                    write(files%log,*) 
                    if (flag_accel==4) then
                        write(output_unit,'(A)') 'SELECTED: HGCMFD Acceleration'
                        write(files%log,'(A)') 'SELECTED: HGCMFD Acceleration'
                    elseif (flag_accel==5) then
                        write(output_unit,'(A)') 'SELECTED: pHGCMFD-based Acceleration'
                        write(files%log,'(A)') 'SELECTED: pHGCMFD-based Acceleration'
                    end if
                    write(output_unit,'(A)') 'Reading the coarse mesh input parameters...'
                    write(files%log,'(A)') 'Reading the coarse mesh input parameters...'




                    ! Custom HGCMFD zones 
                    open(unit=130, file='HGCMFD_Grid.inp', status='old', action='read', iostat=ios)

                    if (ios .NE. 0) then 
                        ! Elem_CMFD = Elem_det in case of missing file 
                        write(output_unit,*) 
                        write(files%log,*) 
                        write(output_unit,'(A)') '[WARNING] HG-CMFD FILE NOT FOUND - Assuming HG-CMFD == Geometry elements Provided'
                        write(files%log,'(A)') '[WARNING] HG-CMFD FILE NOT FOUND - Assuming HG-CMFD == Geometry elements Provided'
                        write(output_unit,*) 
                        write(files%log,*) 
                        allocate(Elem_CMFD(size(Elem_det)))
                        Elem_CMFD=Elem_det

                    else 
                        ! Automatic detection of the number of objects 
                        n=0
                        ios=0
                        do while (ios==0)
                            read(130,'(A)', iostat=ios) file_line
                            if (ios>0) then
                                write(output_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Automatic detection routine'
                                write(error_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Automatic detection routine'
                                exit
                            end if
                            file_line=adjustl(file_line)
                            if (file_line(1:1)=='#') cycle
                            if (file_line(1:1)=='%') n=n+1
                        end do
                        n_elements=n

                        ! Allocation and detection of the geometry details 
                        allocate(Elem_CMFD(n_elements)) 
                        rewind(130)
                        n=0
                        Elem_CMFD(:)%flag_BG=0
                        ios=0
                        do
                            read(130,'(A)', iostat=ios) file_line

                            ! In case of any error 
                            if (ios>0) then
                                write(output_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element reading routine'
                                write(error_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element reading routine'
                                exit
                            end if

                            ! End of file 
                            if (ios < 0) exit

                            ! Skip the comments and the blank lines 
                            if (trim(file_line)=='' .OR. file_line(1:1)=='#') cycle  

                            if (index(file_line, '%') > 0) then
                                ! Different cases 
                                n=n+1
                                Elem_CMFD(n)%ID=n
                                read(file_line(index(file_line, '%')+2:), '(A)') Elem_type

                                select case(trim(Elem_type)) 

                                    case ('Crfc') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,1))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        Elem_CMFD(n)%Fig_ID=1

                                    case ('Trng') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,3))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        Elem_CMFD(n)%Fig_ID=2

                                    case ('Rhomb') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,4))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        Elem_CMFD(n)%Fig_ID=3

                                    case ('Rect') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,4))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        Elem_CMFD(n)%Fig_ID=4

                                    case ('Pent') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,5))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        Elem_CMFD(n)%Fig_ID=5

                                    case ('Hex') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,6))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        Elem_CMFD(n)%Fig_ID=6

                                    case ('Ept') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,7))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        Elem_CMFD(n)%Fig_ID=7

                                    case ('Oct') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,8))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        Elem_CMFD(n)%Fig_ID=8

                                    case ('Enna') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,9))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        read(130, *) Elem_CMFD(n)%vertices(:,9)
                                        Elem_CMFD(n)%Fig_ID=9

                                    case ('Deca') 
                                        allocate(Elem_CMFD(n)%centroid(2))
                                        allocate(Elem_CMFD(n)%vertices(2,10))
                                        read(130, *) Elem_CMFD(n)%centroid
                                        read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        read(130, *) Elem_CMFD(n)%vertices(:,9)
                                        read(130, *) Elem_CMFD(n)%vertices(:,10)
                                        Elem_CMFD(n)%Fig_ID=10

                                    case ('Custom') 
                                        !read 'n_vert' to automatically allocate!
                                        !allocate(Elem_CMFD(n)%centroid(2))
                                        !allocate(Elem_CMFD(n)%vertices(2,10))
                                        !read(130, *) Elem_CMFD(n)%centroid
                                        !read(130, *) Elem_CMFD(n)%vertices(:,1)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,2)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,3)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,4)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,5)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,6)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,7)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,8)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,9)
                                        !read(130, *) Elem_CMFD(n)%vertices(:,10)
                                        !Elem_CMFD(n)%Fig_ID=0
                                        write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                                        write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"

                                    case default 
                                        write(output_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element selection'
                                        write(error_unit,'(A)') 'ERROR IN THE CMFD CARD READING - Element selection'
                                        exit
                                end select

                            elseif (index(file_line, '!box') > 0) then
                                ! Box flagging 
                                Elem_CMFD(n)%flag_BG=1
                                box_mat_ID=n

                            end if
                        end do


                        ! Ghost layer implementation 
                        do l=1, size(Elem_CMFD)
                            Elem_CMFD(l)%centroid(1)=Elem_CMFD(l)%centroid(1)+dx_gh 
                            Elem_CMFD(l)%centroid(2)=Elem_CMFD(l)%centroid(2)+dy_gh!*2
                            Elem_CMFD(l)%vertices(1,:)=Elem_CMFD(l)%vertices(1,:)+dx_gh
                            Elem_CMFD(l)%vertices(2,:)=Elem_CMFD(l)%vertices(2,:)+dy_gh!*2
                        end do

                        close(130)

                    end if

                    write(output_unit,'(A)') '...Done'
                    write(files%log,'(A)') '...Done'
                    write(output_unit,*) 
                    write(files%log,*) 
                
                end if

                ! Allocation & Initialization 
                select case (flag_accel)
                    case(1)
                        Opt_gl%n_tot=Opt_gl%n_x*Opt_gl%n_y
                    case(2)
                        Opt_gl%n_tot=Opt_gl%n_x*Opt_gl%n_y
                    case(3)
                        Opt_gl%n_tot=Opt_gl%n_x*Opt_gl%n_y
                    case(4)
                        Opt_gl%n_tot=size(Elem_CMFD(:))
                    case(5)
                        Opt_gl%n_tot=size(Elem_CMFD(:))
                    case default
                        write(output_unit,'(A)') 'ERROR! WRONG ACCELERATION HANDLING'
                        write(error_unit,'(A)') 'ERROR! WRONG ACCELERATION HANDLING'
                end select

                ! Elem_CMFD to GL_Param 
                allocate(GL_Mesh(Opt_gl%n_tot))    ! number of meshes
                allocate(GL_Param(Opt_gl%n_tot*Opt_gl%n_g))    ! number of meshes and groups

                ! Mesh reordering 
                if (flag_accel==1 .OR. flag_accel==2 .OR. flag_accel==3) then  
                    ! CMFD allocation  

                    do l=1, Opt_gl%n_tot
                        no_faces=size(Elem_CMFD(l)%vertices,2)
                        allocate(GL_Mesh(l)%dl_gl(no_faces))
                        allocate(GL_Mesh(l)%A_gl(no_faces))
                        allocate(GL_Mesh(l)%Neigh_ID(no_faces))
                        allocate(GL_Mesh(l)%Sides_ID(no_faces))

                        do g=1, Opt_gl%n_g 
                            g_tot=(g-1)*Opt_gl%n_tot
                            allocate(GL_Param(g_tot+l)%D_til(no_faces, flag_part))
                            allocate(GL_Param(g_tot+l)%D_hat(no_faces, flag_part))
                            allocate(GL_Param(g_tot+l)%J_part(no_faces, 2))
                            allocate(GL_Param(g_tot+l)%J_net(no_faces))
                            allocate(GL_Param(g_tot+l)%D_til_el(flag_part))
                            allocate(GL_Param(g_tot+l)%D_hat_el(flag_part))
                            allocate(GL_Param(g_tot+l)%J_part_el(2))
                        end do
                    end do
                

                elseif (flag_accel==4 .OR. flag_accel==5) then 
                    ! Assignment of the first mesh as the BG mesh 
                    l=1
                    no_faces=size(Elem_CMFD(box_mat_ID)%vertices,2)
                    allocate(GL_Mesh(l)%dl_gl(no_faces))
                    allocate(GL_Mesh(l)%A_gl(no_faces))

                    do g=1, Opt_gl%n_g 
                        g_tot=(g-1)*Opt_gl%n_tot
                        allocate(GL_Param(g_tot+l)%D_til(no_faces, 1))
                        allocate(GL_Param(g_tot+l)%J_part(no_faces, 2))
                        allocate(GL_Param(g_tot+l)%J_net(no_faces))
                        allocate(GL_Param(g_tot+l)%D_til_el(1))
                        allocate(GL_Param(g_tot+l)%D_hat_el(flag_part))
                        allocate(GL_Param(g_tot+l)%J_part_el(2))
                    end do

                    ! Assignment of the other universes 
                    l=1 ! Skipping the first mesh 
                    do n=1, Opt_gl%n_tot
                    

                        if (n==box_mat_ID) cycle
                        l=l+1
                        no_faces=size(Elem_CMFD(n)%vertices,2)
                        allocate(GL_Mesh(l)%dl_gl(no_faces))
                        allocate(GL_Mesh(l)%A_gl(no_faces))


                        do g=1, Opt_gl%n_g
                            g_tot=(g-1)*Opt_gl%n_tot
                            allocate(GL_Param(g_tot+l)%D_til(no_faces, 1))
                            allocate(GL_Param(g_tot+l)%J_part(no_faces, 2))
                            allocate(GL_Param(g_tot+l)%J_net(no_faces))
                            allocate(GL_Param(g_tot+l)%D_til_el(1))
                            allocate(GL_Param(g_tot+l)%D_hat_el(flag_part))
                            allocate(GL_Param(g_tot+l)%J_part_el(2))
                        end do
                    end do
                end if



                ! Global XS Allocation 
                allocate(  XS_gl%Tot(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  XS_gl%Tra(Opt_gl%n_tot*Opt_gl%n_g))
                allocate( XS_gl%Absr(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(XS_gl%nuFis(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  XS_gl%Fis(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  XS_gl%Chi(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(XS_gl%Scatt(Opt_gl%n_tot*Opt_gl%n_g,Opt_gl%n_g))
                allocate( XS_gl%kFis(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  XS_gl%Rem(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  XS_gl%Dif(Opt_gl%n_tot*Opt_gl%n_g))

                ! Simulation Variables (Vectorized) 
                allocate(       Serv_Vect%Phi(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(   Serv_Vect%Phi_old(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  Serv_Vect%Phi_temp(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(   Serv_Vect%Phi_avg(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(    Serv_Vect%source(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(    Serv_Vect%S_fiss(Opt_gl%n_tot))
                allocate(Serv_Vect%S_fiss_old(Opt_gl%n_tot))
                allocate(         Serv_Vect%a(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(         Serv_Vect%b(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(         Serv_Vect%c(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(         Serv_Vect%d(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(         Serv_Vect%e(Opt_gl%n_tot*Opt_gl%n_g))

                ! Simulation Variables (Matrix) 
                allocate(       Serv_Matr%Phi(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(   Serv_Matr%Phi_old(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(  Serv_Matr%Phi_temp(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(   Serv_Matr%Phi_avg(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(    Serv_Matr%source(Opt_gl%n_tot*Opt_gl%n_g))
                allocate(    Serv_Matr%S_fiss(Opt_gl%n_tot))
                allocate(Serv_Matr%S_fiss_old(Opt_gl%n_tot))
                allocate(        Serv_Matr%MM(Opt_gl%n_tot*Opt_gl%n_g,Opt_gl%n_tot*Opt_gl%n_g))



                ! Calculations 

                ! dl 
                if (flag_accel==1 .OR. flag_accel==2 .OR. flag_accel==3) then   ! dl for CMFD 
                    do n=1, Opt_gl%n_tot ! Roll on the coarse meshes
                        call dl_gl_Auto(Elem_CMFD(n)%vertices, Elem_CMFD(n)%centroid, GL_Mesh(n)%dl_gl(:))  

                    end do

                elseif (flag_accel==4 .OR. flag_accel==5) then     ! dl for HG-CMFD 
                    do n=1, Opt_gl%n_tot ! Roll on the coarse meshes
                        if (Elem_CMFD(n)%flag_BG==1) cycle
                        !! Assuming in-center is ok
                        !call dl_gl_Calc(Elem_CMFD(n)%vertices, Elem_CMFD(n)%centroid, dl)
                        !GL_Mesh(n)%dl_gl(:)=dl
                        GL_Mesh(n)%dl_gl(:)=dx_input
                    end do

                    GL_Mesh(1)%dl_gl=dl_input
                end if

                ! Coarse meshes indexing (lo_gl_map & Elem_CMFD) 
                if (flag_accel==1 .OR. flag_accel==2 .OR. flag_accel==3) then 
                    LO_Mesh(:)%lo_gl_Homo=-1    ! for the ghost meshes


                    do l=1, Opt_lo%n_tot 
                        if (LO_Mesh(l)%ID_Red==-1) cycle

                        ! Map building
                        do n=1, Opt_gl%n_tot

                            call Point_in_Polygon_Test_2D(Elem_CMFD(n)%vertices, LO_Mesh(l)%cent(:), test)
                            if (test) then 
                                LO_Mesh(l)%lo_gl_Homo=n
                                Elem_CMFD(n)%universe=LO_Mesh(l)%lo_gl_Homo
                                exit
                            end if
                        end do

                    end do

                    do l=1, Opt_gl%n_tot
                        GL_Mesh(l)%ID=l
                    end do

                    do l=1, Opt_gl%n_tot
                        i=mod(l-1, Opt_gl%n_x)+1
                        j=int((l-i)/Opt_gl%n_x)+1

                        GL_Mesh(l)%Neigh_ID(:)=-1
                        GL_Mesh(l)%Sides_ID=[3, 4, 1, 2]

                        if (j>1) GL_Mesh(l)%Neigh_ID(1)= GL_Mesh(l-Opt_gl%n_x)%ID
                        if (i<Opt_gl%n_x) GL_Mesh(l)%Neigh_ID(2)= GL_Mesh(l+1)%ID
                        if (j<Opt_gl%n_y) GL_Mesh(l)%Neigh_ID(3)= GL_Mesh(l+Opt_gl%n_x)%ID
                        if (i>1) GL_Mesh(l)%Neigh_ID(4)= GL_Mesh(l-1)%ID


                    end do
                
                elseif (flag_accel==4 .OR. flag_accel==5) then 

                    do l=1, Opt_lo%n_tot
                        if (LO_Mesh(l)%ID_Red .NE. -1) then
                            ! Map building

                            u=1 ! Universe No. 1 is reserved
                            do n=1, Opt_gl%n_tot
                                if (n==box_mat_ID) cycle
                                u=u+1
                                call Point_in_Polygon_Test_2D(Elem_CMFD(n)%vertices, LO_Mesh(l)%cent(:), test)
                                if (test) then 
                                    LO_Mesh(l)%lo_gl_Homo=u    ! Elem_CMFD(n)%ID
                                    Elem_CMFD(n)%universe=LO_Mesh(l)%lo_gl_Homo
                                    exit
                                end if
                            end do

                            if (.not. test) then   ! Check if any of the universes was assigned
                                LO_Mesh(l)%lo_gl_Homo=1
                                Elem_CMFD(box_mat_ID)%universe=1
                            end if

                        else
                            LO_Mesh(l)%lo_gl_Homo=-1
                        end if

                    end do
                
                end if

                ! Geometric parameters calculation 
                do u=1, Opt_gl%n_tot
                    GL_Mesh(u)%A_gl=0.d0
                    GL_Mesh(u)%V_gl=0.d0
            allocate(GL_Mesh(u)%J_lo_gl_face(SIZE(GL_Mesh(u)%A_gl)))
                end do

                ! Volume 
                do l=1, Opt_lo%n_tot
                    u=LO_Mesh(l)%lo_gl_Homo

                    if (u>-1) then
                        GL_Mesh(u)%V_gl=GL_Mesh(u)%V_gl+LO_Mesh(l)%V_lo
                    end if
                end do

                ! Area 
                do u=1, Opt_gl%n_tot 
                    no_faces=size(GL_Mesh(u)%A_gl)
                    if (no_faces==1) then   ! Circumference  

                        call Distance(Elem_CMFD(u)%centroid(1), Elem_CMFD(u)%centroid(2), Elem_CMFD(u)%vertices(1,1), Elem_CMFD(u)%vertices(2,1), dl)
                        GL_Mesh(u)%A_gl(1)= atan(1.d0)*4.d0*dl**2.d0
                    else ! n-sided polygon  
                        vert_temp(:,1:no_faces)=Elem_CMFD(u)%vertices(:,1:no_faces)
                        vert_temp(:,1+no_faces)=Elem_CMFD(u)%vertices(:,1)

                        do i=1, no_faces
                            call Distance(vert_temp(1,i), vert_temp(2,i), vert_temp(1,i+1), vert_temp(2,i+1), dl)
                            GL_Mesh(u)%A_gl(i)=dl
                        end do

            ! Face ID                
            select case(Elem_CMFD(u)%Fig_ID)
            case (2)
                GL_Mesh(u)%J_lo_gl_face(1)=1
                GL_Mesh(u)%J_lo_gl_face(2)=5
                GL_Mesh(u)%J_lo_gl_face(3)=6
            case (3)    
                GL_Mesh(u)%J_lo_gl_face(1)=5
                GL_Mesh(u)%J_lo_gl_face(2)=6
                GL_Mesh(u)%J_lo_gl_face(3)=7
                GL_Mesh(u)%J_lo_gl_face(4)=8
            case (4)    
                GL_Mesh(u)%J_lo_gl_face(1)=1    
                GL_Mesh(u)%J_lo_gl_face(2)=2    
                GL_Mesh(u)%J_lo_gl_face(3)=3    
                GL_Mesh(u)%J_lo_gl_face(4)=4   
                ! Rhomboid
                !GL_Mesh(1)%J_lo_gl_face(1)=1    
                !GL_Mesh(1)%J_lo_gl_face(2)=2    
                !GL_Mesh(1)%J_lo_gl_face(3)=3    
                !GL_Mesh(1)%J_lo_gl_face(4)=4    
                !GL_Mesh(u)%J_lo_gl_face(1)=8    
                !GL_Mesh(u)%J_lo_gl_face(2)=5    
                !GL_Mesh(u)%J_lo_gl_face(3)=6    
                !GL_Mesh(u)%J_lo_gl_face(4)=7    
            case (5)    
                write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
            case (6)
                GL_Mesh(u)%J_lo_gl_face(1)=1
                GL_Mesh(u)%J_lo_gl_face(2)=8
                GL_Mesh(u)%J_lo_gl_face(3)=5
                GL_Mesh(u)%J_lo_gl_face(4)=3
                GL_Mesh(u)%J_lo_gl_face(5)=6
                GL_Mesh(u)%J_lo_gl_face(6)=7
            case (7)
                write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
            case (8)
                GL_Mesh(u)%J_lo_gl_face(1)=1
                GL_Mesh(u)%J_lo_gl_face(2)=8
                GL_Mesh(u)%J_lo_gl_face(3)=2
                GL_Mesh(u)%J_lo_gl_face(4)=5
                GL_Mesh(u)%J_lo_gl_face(5)=3
                GL_Mesh(u)%J_lo_gl_face(6)=6
                GL_Mesh(u)%J_lo_gl_face(7)=4
                GL_Mesh(u)%J_lo_gl_face(8)=7
            case (9)
                write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
            case (10)
                write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
            case (0)
                write(output_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
                write(error_unit,'(A)') " ERROR! FIGURE NOT IMPLEMENTED YET!"
            
            
            end select                

                    end if
                end do

                ! Energy groups 
                do l=1, Opt_lo%n_tot
                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot
                        if (Opt_gl%n_g==1) then 
                            LO_Param(t_tot+l)%lo_gl_Cond=1
                            LO_Param(t_tot+l)%En_ID=1
                        elseif (Opt_gl%n_g==Opt_lo%n_g) then
                            LO_Param(t_tot+l)%lo_gl_Cond=t
                            LO_Param(t_tot+l)%En_ID=t
                        else
                            write(output_unit,'(A)') 'ERROR! WRONG GROUP SELECTION IN THE MG GLOBAL PROBLEM'
                            write(error_unit,'(A)') 'ERROR! WRONG GROUP SELECTION IN THE MG GLOBAL PROBLEM'
                        end if
                    end do
                end do

                ! Currents (flags and Neigh_ID extraction) 
                no_faces=size(LO_Mesh(1)%Neigh_ID)  ! Local meshes number of faces
                do l=1, Opt_lo%n_tot
                    allocate(LO_Mesh(l)%J_lo_gl_face(no_faces))   ! Contains the face of the absolute Neighs
                    allocate(LO_Mesh(l)%J_lo_gl_ID(no_faces))   ! Contains the ID of the absolute Neighs

                    LO_Mesh(l)%J_lo_gl_face=-1
                    LO_Mesh(l)%J_lo_gl_ID=-1

                    if (LO_Mesh(l)%ID_Red>0) then
                        u=LO_Mesh(l)%lo_gl_Homo   ! Universe

                        do n=1, no_faces
                            nb=LO_Mesh(l)%Neigh_ID(n)   ! Neighbour ID
                            uu=LO_Mesh(nb)%lo_gl_Homo   ! Neigbour universe

                            if ((u .NE. uu) .AND. (LO_Mesh(nb)%ID_Red>0)) then
                                LO_Mesh(l)%J_lo_gl_face(n)=LO_Mesh(l)%Sides_Neigh(n)
                                LO_Mesh(l)%J_lo_gl_ID(n)=1
                            else
                                LO_Mesh(l)%J_lo_gl_ID(n)=0
                            end if
                        end do
                    end if
                end do
        end subroutine Accel_Ini
    
        subroutine Distance(x1, y1, x2, y2, dist) 
            real(dp), intent(in) :: x1, y1, x2, y2
            real(dp), intent(out) :: dist

            dist = sqrt((x2 - x1)**2 + (y2 - y1)**2)
        end subroutine Distance
    
        subroutine Homogenization_DerType(Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, XS_lo, XS_gl, Phi_lo, Phi_gl) 
            type(Options_Data), intent(in) :: Opt_lo, Opt_gl
            type(LO_geom), allocatable, intent(in) :: LO_Mesh(:)
            type(LO_coeff), allocatable, intent(in) :: LO_Param(:)
            type(GL_geom), allocatable, intent(in) :: GL_Mesh(:)
            real(dp), allocatable, intent(in) :: Phi_lo(:)
            real(dp), allocatable, intent(inout) :: Phi_gl(:)
            type(XS_Data), intent(in) :: XS_lo
            type(XS_Data), intent(inout) :: XS_gl

            integer :: u, n, l, l_red, t, t_red, t_tot, tt, tt_red, tt_tot, g, g_red, g_tot, gg, gg_red, gg_tot, i, j
            real(dp) :: temp
            real(dp), allocatable :: Phi_int(:)
            type(XS_Data) :: XS_int

            ! Allocation and initialization 
            allocate(  XS_int%Tot(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Tra(Opt_gl%n_tot*Opt_gl%n_g))
            allocate( XS_int%Absr(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(XS_int%nuFis(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Fis(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Chi(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(XS_int%Scatt(Opt_gl%n_tot*Opt_gl%n_g, Opt_gl%n_g))
            allocate( XS_int%kFis(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Rem(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(  XS_int%Dif(Opt_gl%n_tot*Opt_gl%n_g))
            allocate(     Phi_int(Opt_gl%n_tot*Opt_gl%n_g))

            XS_int%Tot=0.d0
            XS_int%Tra=0.d0
            XS_int%Absr=0.d0
            XS_int%nuFis=0.d0
            XS_int%Fis=0.d0
            XS_int%Chi=0.d0
            XS_int%Scatt=0.d0
            XS_int%kFis=0.d0
            XS_int%Rem=0.d0
            XS_int%Dif=0.d0
            Phi_int=0.d0

            ! Element-wise integral 
            do l=1, Opt_lo%n_tot
                n=LO_Mesh(l)%lo_gl_Homo ! Universe
                l_red=LO_Mesh(l)%ID_Red
                if (l_red>0) then ! Ghost layer filtering

                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot
                        t_red=(t-1)*Opt_lo%n_red
                        g=LO_Param(t_tot+l)%En_ID
                        g_tot=(g-1)*Opt_gl%n_tot

                          XS_int%Tot(g_tot+n)=  XS_int%Tot(g_tot+n)+  XS_lo%Tot(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                          XS_int%Tra(g_tot+n)=  XS_int%Tra(g_tot+n)+  XS_lo%Tra(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                         XS_int%Absr(g_tot+n)= XS_int%Absr(g_tot+n)+ XS_lo%Absr(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                        XS_int%nuFis(g_tot+n)=XS_int%nuFis(g_tot+n)+XS_lo%nuFis(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                          XS_int%Fis(g_tot+n)=  XS_int%Fis(g_tot+n)+  XS_lo%Fis(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                         XS_int%kFis(g_tot+n)= XS_int%kFis(g_tot+n)+ XS_lo%kFis(t_tot+l)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                             Phi_int(g_tot+n)=     Phi_int(g_tot+n)+                     Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo

                        do tt=1, Opt_lo%n_g
                            tt_tot=(tt-1)*Opt_lo%n_tot
                            tt_red=(tt-1)*Opt_lo%n_red
                            gg=LO_Param(tt_tot+l)%En_ID

                            XS_int%Scatt(g_tot+n,gg)=XS_int%Scatt(g_tot+n,gg)+XS_lo%Scatt(t_tot+l,tt)*Phi_lo(t_red+l_red)*LO_Mesh(l)%V_lo
                            XS_int%Chi(g_tot+n)=XS_int%Chi(g_tot+n)+XS_lo%Chi(t_tot+l)*XS_lo%nuFis(tt_tot+l)*Phi_lo(tt_red+l_red)*LO_Mesh(l)%V_lo
                        end do
                    end do
                end if
            end do

            ! Homogenization 
            do l=1, Opt_gl%n_tot
                do g=1, Opt_gl%n_g
                    g_tot=(g-1)*Opt_gl%n_tot
                    XS_gl%Tot(g_tot+l)=  XS_int%Tot(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%Tra(g_tot+l)=  XS_int%Tra(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%Absr(g_tot+l)= XS_int%Absr(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%nuFis(g_tot+l)=XS_int%nuFis(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%Fis(g_tot+l)=  XS_int%Fis(g_tot+l)/Phi_int(g_tot+l)
                    XS_gl%kFis(g_tot+l)= XS_int%kFis(g_tot+l)/Phi_int(g_tot+l)

                    XS_gl%Scatt(g_tot+l,:)=XS_int%Scatt(g_tot+l,:)/Phi_int(g_tot+l)

                    Phi_gl(g_tot+l)=Phi_int(g_tot+l)/GL_Mesh(l)%V_gl
                end do
            end do

            ! Chi 
            do l=1, Opt_gl%n_tot
                temp=0.d0
                do g=1, Opt_gl%n_g
                    g_tot=(g-1)*Opt_gl%n_tot
                    temp=temp+XS_gl%nuFis(g_tot+l)*Phi_gl(g_tot+l)*GL_Mesh(l)%V_gl
                end do
                if (abs(temp)>1.d-15) then
                    do g=1, Opt_gl%n_g
                        g_tot=(g-1)*Opt_gl%n_tot
                        XS_gl%Chi(g_tot+l)=XS_int%Chi(g_tot+l)/temp
                    end do
                else
                    do g=1, Opt_gl%n_g
                        g_tot=(g-1)*Opt_gl%n_tot
                        XS_gl%Chi(g_tot+l)=0.d0
                    end do
                end if
            end do

            ! Diffusion coefficient 
            XS_gl%Dif=1.d0/3.d0/XS_gl%Tra

            ! Removal coefficient 
            do g=1, Opt_gl%n_g
                g_tot=(g-1)*Opt_gl%n_tot
                XS_gl%Rem(g_tot+1:g_tot+Opt_gl%n_tot)=XS_gl%Absr(g_tot+1:g_tot+Opt_gl%n_tot)
                do gg=1, Opt_gl%n_g
                    if (gg .NE. g) XS_gl%Rem(g_tot+1:g_tot+Opt_gl%n_tot)=XS_gl%Rem(g_tot+1:g_tot+Opt_gl%n_tot)+XS_gl%Scatt(g_tot+1:g_tot+Opt_gl%n_tot, gg)
                end do
            end do

        end subroutine Homogenization_DerType
    
        subroutine Modulation(Opt_lo, Opt_gl, LO_Mesh, LO_Param, Phi_gl, Phi_avg_lo, Phi_lo) 
        
            type(Options_Data), intent(in) :: Opt_lo, Opt_gl
            type(LO_geom), intent(in) :: LO_Mesh(:)
            type(LO_coeff), intent(in) :: LO_Param(:)
            real(dp), intent(in) :: Phi_gl(:), Phi_avg_lo(:)
            real(dp), intent(inout) :: Phi_lo(:)

            integer :: l, ll, t, g, u, t_tot, t_red, g_tot

            do l=1, Opt_lo%n_tot
                u=LO_Mesh(l)%lo_gl_Homo ! Global spatial index relative to the total grid (w/ ghosts) 
                ll=LO_Mesh(l)%ID_Red ! Local spatial index relative to the reduced grid (w/o ghosts) 
                if (ll>0) then
                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot    ! Local energy index base relative to the total grid (w/ ghosts)
                        t_red=(t-1)*Opt_lo%n_red    ! Local energy index base relative to the reduced grid (w/o ghosts) 

                        g=LO_Param(t_tot+l)%En_ID  ! Global energy index
                        g_tot=(g-1)*Opt_gl%n_tot    ! Global energy index base

                        Phi_lo(t_red+ll)=Phi_lo(t_red+ll)*Phi_gl(g_tot+u)/Phi_avg_lo(g_tot+u)
                    end do
                end if
            end do

        end subroutine Modulation 
        
        subroutine Currents_lo_to_gl_DerType(flag_accel, Opt_lo, Opt_gl, LO_Mesh, LO_Param, GL_Mesh, GL_Param) 
            integer, intent(in) :: flag_accel
            type(Options_Data), intent(in) :: Opt_lo, Opt_gl
            type(LO_geom), intent(in) :: LO_Mesh(:)
            type(LO_coeff), intent(in) :: LO_Param(:)
            type(GL_geom), intent(in) :: GL_Mesh(:)
            type(GL_coeff), intent(inout) :: GL_Param(:)

            integer :: i, g, l, u, n, t, t_tot, g_tot, infc, face_ID, face_ID_GL ! In Line Face Calc.      
            real(dp) :: Phi_S

            ! Initialization 
            do i=1, Opt_gl%n_tot*Opt_gl%n_g
                GL_Param(i)%J_net(:)=0.d0
                GL_Param(i)%J_part(:,:)=0.d0
            end do


            ! Integration 
            do l=1, Opt_lo%n_tot
                u=LO_Mesh(l)%lo_gl_Homo
                if (u>-1) then  ! ghost mesh filtering

                    do t=1, Opt_lo%n_g
                        t_tot=(t-1)*Opt_lo%n_tot
                        g=LO_Param(t_tot+l)%En_ID
                        g_tot=(g-1)*Opt_gl%n_tot

                        do n=1, size(LO_Mesh(l)%Neigh_ID)     ! Iteration on the surfaces

                            if (flag_accel==4 .OR. flag_accel==5) then ! HGCMFD 
                                if (LO_Mesh(l)%J_flag(n)>0 .AND. u>1) then  !  .AND. LO_Mesh(l)%BC(n)==0
                                    face_ID=LO_Mesh(l)%Sides_ID(n) ! ID Face along which the universe changes
                                    face_ID_GL=findloc(GL_Mesh(u)%J_lo_gl_face, face_ID, dim=1) ! ID Face along which the universe changes in the global

                                    ! Net 
                                    GL_Param(g_tot+u)%J_net(face_ID_GL)=GL_Param(g_tot+u)%J_net(face_ID_GL)+LO_Param(t_tot+l)%J_net(n)*LO_Mesh(l)%A_lo(n)

                                    ! Partial 
                                    GL_Param(g_tot+u)%J_part(face_ID_GL,1)=GL_Param(g_tot+u)%J_part(face_ID_GL,1)+LO_Param(t_tot+l)%J_part(n,1)*LO_Mesh(l)%A_lo(n)
                                    GL_Param(g_tot+u)%J_part(face_ID_GL,2)=GL_Param(g_tot+u)%J_part(face_ID_GL,2)+LO_Param(t_tot+l)%J_part(n,2)*LO_Mesh(l)%A_lo(n)
                                
                                elseif (LO_Mesh(l)%BC(n)>0 .AND. u==1) then ! .AND. LO_Mesh(l)%BC(n)==1
                                    face_ID=LO_Mesh(l)%Sides_ID(n) ! ID Face along which the universe changes
                                    face_ID_GL=findloc(GL_Mesh(u)%J_lo_gl_face, face_ID, dim=1)

                                    ! Net
                                    GL_Param(g_tot+u)%J_net(face_ID_GL)=GL_Param(g_tot+u)%J_net(face_ID_GL)+LO_Param(t_tot+l)%J_net(n)*LO_Mesh(l)%A_lo(n)

                                    ! Partial 
                                    GL_Param(g_tot+u)%J_part(face_ID_GL,1)=GL_Param(g_tot+u)%J_part(face_ID_GL,1)+LO_Param(t_tot+l)%J_part(n,1)*LO_Mesh(l)%A_lo(n)
                                    GL_Param(g_tot+u)%J_part(face_ID_GL,2)=GL_Param(g_tot+u)%J_part(face_ID_GL,2)+LO_Param(t_tot+l)%J_part(n,2)*LO_Mesh(l)%A_lo(n)

                                end if
                            elseif (flag_accel==1 .OR. flag_accel==2 .OR. flag_accel==3) then !CMFD like 
                                if ((LO_Mesh(l)%J_lo_gl_ID(n)>0 .OR. LO_Mesh(l)%BC(n)>0)) then  
                                    face_ID=LO_Mesh(l)%Sides_ID(n) ! ID Face along which the universe changes
                                    face_ID_GL=findloc(GL_Mesh(u)%J_lo_gl_face, face_ID, dim=1) ! ID Face along which the universe changes in the global

                                    ! Net
                                    GL_Param(g_tot+u)%J_net(face_ID_GL)=GL_Param(g_tot+u)%J_net(face_ID_GL)+LO_Param(t_tot+l)%J_net(n)*LO_Mesh(l)%A_lo(n)


                                    ! Partial 
                                    GL_Param(g_tot+u)%J_part(face_ID_GL,1)=GL_Param(g_tot+u)%J_part(face_ID_GL,1)+LO_Param(t_tot+l)%J_part(n,1)*LO_Mesh(l)%A_lo(n)
                                    GL_Param(g_tot+u)%J_part(face_ID_GL,2)=GL_Param(g_tot+u)%J_part(face_ID_GL,2)+LO_Param(t_tot+l)%J_part(n,2)*LO_Mesh(l)%A_lo(n)

                                end if
                            end if
                        end do
                    end do
                end if
            end do

            ! Condensation & Homogenization 
            do l=1, Opt_gl%n_tot
                do n=1, size(GL_Mesh(l)%A_gl)
                    do g=1, Opt_gl%n_g
                        g_tot=(g-1)*Opt_gl%n_tot

                        ! Net
                        GL_Param(g_tot+l)%J_net(n)=GL_Param(g_tot+l)%J_net(n)/GL_Mesh(l)%A_gl(n)

                        ! Partial
                        GL_Param(g_tot+l)%J_part(n,1)=GL_Param(g_tot+l)%J_part(n,1)/GL_Mesh(l)%A_gl(n)
                        GL_Param(g_tot+l)%J_part(n,2)=GL_Param(g_tot+l)%J_part(n,2)/GL_Mesh(l)%A_gl(n)
                    end do
                end do
            end do

        end subroutine Currents_lo_to_gl_DerType
    
        subroutine SOR_Accel(w_SOR, Var_Pre, Final_Val) 
            real(dp), intent(in) :: w_SOR, Var_Pre(:)
            real(dp), intent(inout) :: Final_Val(:)

            Final_Val=(1-w_SOR)*Var_Pre+w_SOR*Final_Val
        end subroutine SOR_Accel
    
end module Service_Fcns

