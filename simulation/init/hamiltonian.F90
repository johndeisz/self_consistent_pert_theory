MODULE Hamiltonian

  character*128 :: hk_file
  integer :: N_dim
  double precision :: a(3,3)
  
  integer :: N_a
  double precision, dimension(:,:), allocatable :: r_atom

  double precision, dimension(:), allocatable :: ed
  integer, dimension(:), allocatable :: N_o

CONTAINS	
  subroutine read_hamiltonian()
	
    ! MPI variables
    include 'mpif.h'
    INTEGER :: rank, ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    if (rank .eq. 0) then
       open(unit=15, file=hk_file, status='old')
       read(15,*)
       read(15,*) N_dim
       read(15,*) a(1,1), a(1,2), a(1,3)
       read(15,*) a(2,1), a(2,2), a(2,3)
       read(15,*) a(3,1), a(3,2), a(3,3)
    endif
     	
    call MPI_Bcast(N_dim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(a, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (rank .eq. 0) then
       read(15,*)
       read(15,*) 
       read(15,*) N_a
    endif
     	
    call MPI_Bcast(N_a, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    allocate (r_atom(N_a,3))

    write(6,*) 'N_a upon read = ', N_a
    
    if (rank .eq. 0) then
       do i_a = 1, N_a
          read(15,*) r_atom(i_a, 1), r_atom(i_a, 2), r_atom(i_a, 3)
       enddo
    endif
 
    call MPI_Bcast(r_atom, N_a*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  
    allocate (N_o(N_a))

    if (rank .eq. 0) then  
       read(15,*)
       read(15,*)
       nb = 0
       do i_a = 1, N_a
          read(15,*) N_o(i_a)
          write(6,*) 'N_o(i_a) = ', N_o(i_a)
          nb = nb + N_o(i_a)
       enddo
    endif
    call MPI_Bcast(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(N_o, N_a, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    allocate (ed(nb))
     	
    if (rank .eq. 0) then
       read(15,*) 
       read(15,*)
       do ib = 1, nb
          read(15,*) i_dummy_1, i_dummy_2, ed(ib)
       enddo
       write(6,*) 'orbital energies '
       ind = 0
       do i_a = 1, N_a
          do i_o = 1, N_o(i_a)
             ind = ind + 1
             write(6,*) i_a, i_O, ed(ind)
          enddo
       enddo
       write(6,*)	
       close(unit=15)
    endif
    call MPI_Bcast(ed, nb, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)  	
     	
  end subroutine read_hamiltonian

end MODULE Hamiltonian
