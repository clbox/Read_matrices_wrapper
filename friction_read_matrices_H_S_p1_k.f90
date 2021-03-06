    subroutine friction_read_matrices_H_S_p1 &
            ( first_order_S, first_order_H, first_order_S_cmplx, first_order_H_cmplx,&
            n_spin, n_basis, ik_point, n_cells_in_hamiltonian,&
            friction_n_active_atoms, real_eigenvectors, friction_index_list)

        implicit none
        
        include 'mpif.h'
        !  ARGUMENTS
        !  imported variables


        !ADDED  BY CONNOR FOR F2PY


        logical, save :: friction_use_complex_matrices = .true.
        logical, intent(in) :: real_eigenvectors
        integer  :: n_tasks, myid
        !
        integer, intent(in) :: n_spin, n_basis, ik_point, n_cells_in_hamiltonian
        integer, intent(in) :: friction_n_active_atoms        
        character(len=50),dimension(n_spin) :: file_name

        integer, dimension(2), intent(in) :: friction_index_list
        real*8, intent(out) :: first_order_S(3, friction_n_active_atoms, 1, n_basis*(n_basis+1)/2, n_spin)
        real*8, intent(out) :: first_order_H(3, friction_n_active_atoms, 1,  n_basis*(n_basis+1)/2,1,n_spin)
        complex*16, intent(out) :: first_order_S_cmplx(3, friction_n_active_atoms, 1, n_basis*(n_basis+1)/2, n_spin)
        complex*16, intent(out) :: first_order_H_cmplx(3, friction_n_active_atoms, 1, n_basis*(n_basis+1)/2,1,n_spin)
        
        !f2py intent(out) :: first_order_S
        !f2py intent(out) :: first_order_H
        !f2py intent(out) :: first_order_S_cmplx
        !f2py intent(out) :: first_order_S_cmplx

        !  counters
        integer :: i_cart, i_atom, i_basis, j_basis, i_spin
        integer :: i_cell, i_index, i_k_task, n_k_points, i_k_point
        integer :: mpierror, numprocs
        integer status(MPI_STATUS_SIZE)

        integer, dimension(:,:),allocatable :: basis_index
        
        call MPI_Init(mpierror)
        call MPI_Comm_size(MPI_COMM_WORLD, n_tasks, mpierror)
        call MPI_Comm_rank(MPI_COMM_WORLD, myid, mpierror)

        allocate(basis_index(n_basis,n_basis))
        

        write(*,*) "Let's begin" 
        basis_index = 0
        i_index = 0
        do j_basis = 1, n_basis, 1
          do i_basis = 1, j_basis,1
            i_index = i_index + 1
            basis_index(i_basis,j_basis) = i_index
          enddo
        enddo

        n_k_points = ik_point
        !  write basis function properties
        if (friction_use_complex_matrices) then
        write(*,*) "I am using complex matrices" 
            i_k_task = 0
            !bodged it to do only the requested i_k_point, see other file for
            !'proper, multi k version
            do i_k_point = ik_point, n_k_points, 1
            if (myid.eq.MOD(i_k_point,n_tasks) .and. myid <= n_k_points) then
                write(*,*) "I am task: ", myid, "with k-point:", i_k_point 
                i_k_task = i_k_task + 1
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the overlap matrix
                    write(*,*) "Writing overlap matrix" 
                    write(*,*) friction_index_list(i_atom)
                    if (i_cart==1) then 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_x.out'
                    else if (i_cart==2) then
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_y.out'
                    else 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_k",i_k_point,"_",&
                            & friction_index_list(i_atom),'_z.out'
                    endif

                    open (50, file=file_name(1),status='old',form='unformatted')
                    if (real_eigenvectors) then
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_S(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1),&
                                i_basis=1, j_basis )
                        enddo
                    else
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_S_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1),&
                                i_basis=1, j_basis )
                        enddo
                    endif
                    close(50)
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then

                        write(*,*) "spin = 1 and writing Hamiltonian" 
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        write(*,*) "spin = 2 and writing Hamiltonian" 
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_k",i_k_point,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if
                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),status='old',form='unformatted')
                        if (real_eigenvectors) then
                            do j_basis = 1, n_basis, 1
                                read(50) (first_order_H(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1, j_basis )
                            enddo
                        else
                            do j_basis = 1, n_basis, 1
                                read(50) (first_order_H_cmplx(i_cart,i_atom,i_k_task,basis_index(i_basis,j_basis),1,i_spin),&
                                    i_basis=1, j_basis )
                            enddo
                        endif
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            endif
            enddo
        else
            write(*,*) "Real matrices?" 
            do i_cell=1, n_cells_in_hamiltonian
                do i_atom = 1, friction_n_active_atoms, 1
                do i_cart = 1, 3, 1
                    !  now, write the overlap matrix
                    if (i_cart==1) then 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_x.out'
                    else if (i_cart==2) then
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_y.out'
                    else 
                        write(file_name(1),'(A,I0.3,A,I0.3,A)') "first_order_S_N",i_cell,"_",&
                            & friction_index_list(i_atom),'_z.out'
                    endif
                    open (50, file=file_name(1),status='old',form='unformatted')
                    do j_basis = 1, n_basis, 1
                        read(50) (first_order_S(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1),&
                            i_basis=1, j_basis )
                    enddo
                    close(50)
                    
                    !  now, write the Hamiltonian matrix
                    if (n_spin.eq.1) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    else if (n_spin.eq.2) then
                        if (i_cart==1) then 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_x.out'
                        else if (i_cart==2) then
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_y.out'
                        else 
                            write(file_name(1),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_up_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                            write(file_name(2),'(A,I0.3,A,I0.3,A)') &
                                "first_order_H_dn_N",i_cell,"_",friction_index_list(i_atom),'_z.out'
                        endif
                    end if

                    do i_spin = 1, n_spin, 1
                        open (50, file=file_name(i_spin),status='old',form='unformatted')
                        do j_basis = 1, n_basis, 1
                            read(50) (first_order_H(i_cart,i_atom,i_cell,basis_index(i_basis,j_basis),1,i_spin),&
                                i_basis=1, j_basis )
                        enddo
                        close(50)
                    ! spin
                    enddo
                enddo
                enddo
            enddo
        
        endif

        deallocate(basis_index)
    call MPI_Finalize(mpierror)
    end subroutine friction_read_matrices_H_S_p1
