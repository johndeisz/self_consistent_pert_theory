MODULE Hamiltonian
  USE CONSTANTS 

  character*128 :: hk_file
  integer :: N_dim
  double precision :: a(3,3)
  
  integer :: N_a
  double precision, dimension(:,:), allocatable :: r_atom

  double precision, dimension(:), allocatable :: ed
  integer, dimension(:), allocatable :: N_o

  double complex, dimension(:,:), allocatable :: tij_temp
  double complex, dimension(:,:,:,:,:), allocatable :: tij
  
CONTAINS	
  subroutine read_hamiltonian()
	
    ! MPI variables
    implicit none
    include 'mpif.h'
    INTEGER :: rank, ierr
    INTEGER :: i_a, i_dummy_1, i_dummy_2, i1, i2, i3, jb, n_points, ip, ibp
    INTEGER :: nb, i_o, ib, ind
    double precision :: im_tij, re_tij

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
    endif
    call MPI_Bcast(ed, nb, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)  	

    allocate (tij(0:nb-1,0:nb-1,0:Nl(1)-1,0:Nl(2)-1,0:Nl(3)-1))
    tij = dcmplx(0.0d0, 0.0d0)
    allocate (tij_temp(0:nb-1, 0:nb-1))
    
    if (rank .eq. 0) then

     read(15,*) 
     read(15,*)
     read(15,*) n_points
  endif
  call MPI_Bcast(n_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  
  do ip = 1, n_points
     if (rank .eq. 0) then
        read(15,*) i1, i2, i3
        i1 = mod(i1+Nl(1),Nl(1))
        i2 = mod(i2+Nl(2),Nl(2))
        i3 = mod(i3+Nl(3),Nl(3))
        write(6,*) 'i1 = ', i1, ' i2 = ', i2, ' i3=', i3
        do jb = 1, nb*nb
           read(15,*) ib, ibp, re_tij, im_tij
           tij_temp(ib, ibp) = dcmplx(re_tij, im_tij)
        enddo
     endif
     call MPI_Bcast(i1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(i2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(i3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(tij_temp, nb*nb, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
     tij(:,:, i1, i2, i3) = tij_temp(:,:)
     
  enddo
  if (rank .eq. 0) then
     close(unit=15)
  endif

!!$     do ix = -max_x, max_x
!!$        do iy = -max_y, max_y
!!$           do iz = -max_z, max_z
!!$
!!$              read(5,*)
!!$              read(5,*)
!!$              write(6,*) 
!!$              write(6,200) ix, iy, iz
!!$              
!!$              k = + mod(iy+lly,lly)*llx + &
!!$                   mod(iz+llz,llz)*llx*lly
!!$
!!$              do ib = 0, nb-1
!!$                 do ibp = 0, nb-1
!!$                    read(5,*) id, idp, tij(ib,ibp,ix,iy,iz)
!!$                    write(6,300) ib, ibp, real(tij(ib,ibp,ix,iy,iz)), &
!!$                         aimag(tij(ib,ibp,ix,iy,iz))
!!$
!!$                 enddo
!!$              enddo
!!$
!!$           enddo
!!$        enddo
!!$     enddo

  end subroutine read_hamiltonian

end MODULE Hamiltonian
