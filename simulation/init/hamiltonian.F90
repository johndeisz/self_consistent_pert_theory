MODULE Hamiltonian
  USE CONSTANTS 

  character*128 :: hk_file, hu_file, hso_file, hfield_file
  integer :: N_dim
  double precision :: a(3,3)
  
  integer :: N_a, N_int_a, nb
  double precision, dimension(:,:), allocatable :: r_atom

  double precision, dimension(:), allocatable :: ed
  integer, dimension(:), allocatable :: N_o

  double complex, dimension(:,:,:,:,:), allocatable :: tij

  double precision, dimension(:,:,:), allocatable :: uup, uuj
  double complex, dimension(:,:,:), allocatable :: gamma0_ph
  integer :: max_orb_per_atom, max_orb_per_int_atom

  double complex, dimension(:,:,:), allocatable :: h_so

  double precision :: prfld

  double precision, dimension(:,:), allocatable :: h

  double precision :: prfld_pert
  double precision, dimension(:,:), allocatable :: h_pert
  double precision, dimension(:), allocatable :: v_pert

  
CONTAINS	
  subroutine read_hamiltonian()
	
    ! MPI variables
    implicit none
    include 'mpif.h'
    INTEGER :: rank, ierr
    INTEGER :: i_a, i_dummy_1, i_dummy_2, i1, i2, i3, jb, n_points, ip, ibp
    INTEGER :: is, isp
    INTEGER :: i_o, ib, ind, N_o_int, m_size, N_o_tmp
    double precision :: im_tij, re_tij, uuij, ujij, re_hso, im_hso, so_amp
    double precision, dimension(:,:), allocatable :: uu_temp, uj_temp
    LOGICAL :: i_terminate, int_flag, so_flag
    double complex, dimension(:,:), allocatable :: h_so_tmp
    double precision :: h_temp(1:3)
    double complex, dimension(:,:), allocatable :: tij_temp
    
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
    allocate (r_atom(0:N_a-1,3))

    if (rank .eq. 0) then
       write(6,*) 'N_a upon read = ', N_a
       do i_a = 0, N_a-1
          read(15,*) r_atom(i_a, 1), r_atom(i_a, 2), r_atom(i_a, 3)
       enddo
    endif
 
    call MPI_Bcast(r_atom, N_a*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  
    allocate (N_o(0:N_a-1))

    if (rank .eq. 0) then  
       read(15,*)
       read(15,*)
       nb = 0
       N_int_a = 0
       max_orb_per_atom = 0
       max_orb_per_int_atom = 0
       do i_a = 0, N_a-1
          read(15,*) int_flag, N_o(i_a)
          write(6,*) 'N_o(i_a) = ', N_o(i_a)
          if (int_flag) then
             N_int_a = N_int_a + 1
             if (N_o(i_a) .gt. max_orb_per_int_atom) then
                max_orb_per_int_atom = N_o(i_a)
             endif
          endif
          if (N_o(i_a) .gt. max_orb_per_atom) then
             max_orb_per_atom = N_o(i_a)
          endif
          nb = nb + N_o(i_a)
       enddo
       write(6,*) 'N_int_a = ', N_int_a
       write(6,*) 'nb = ', nb
       write(6,*) 'max_orb_per_int_atom = ', max_orb_per_int_atom
       write(6,*) 'max_orb_per_atom = ', max_orb_per_atom
    endif
    
    call MPI_Bcast(nb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(N_o, N_a, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(N_int_a, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(max_orb_per_atom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(max_orb_per_int_atom, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    allocate (ed(0:nb-1))
     	
    if (rank .eq. 0) then
       read(15,*) 
       read(15,*)
       do ib = 0, nb-1
          read(15,*) i_dummy_1, i_dummy_2, ed(ib)
       enddo
       write(6,*) 'orbital energies '
       ind = 0
       do i_a = 0, N_a-1
          do i_o = 0,  N_o(i_a)-1
             write(6,*) i_a, i_o, ed(ind)
             ind = ind + 1
          enddo
       enddo
       write(6,*)	
    endif
    call MPI_Bcast(ed, nb, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)  	

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

  !----------------------Interaction Parameters ------------------------

  allocate (uu_temp(0:max_orb_per_int_atom-1,0:max_orb_per_int_atom-1))
  allocate (uj_temp(0:max_orb_per_int_atom-1,0:max_orb_per_int_atom-1))
  
  if (rank .eq. 0) then
     open(unit=25, file=hu_file, status='old')
     read(25,*)
  endif
  
  allocate (uup(0:N_int_a-1, 0:max_orb_per_int_atom-1, 0:max_orb_per_int_atom-1))
  allocate (uuj(0:N_int_a-1, 0:max_orb_per_int_atom-1, 0:max_orb_per_int_atom-1))

  uup = 0.0d0
  uuj = 0.0d0
  
  do i_a = 0, N_int_a-1
     uu_temp = 0.0d0
     uj_temp = 0.0d0
     i_terminate = .false.
     
     if (rank .eq. 0) then

        read(25,*) N_o_int
        if ( (N_o_int .ne. N_o(i_a)) ) then
           write(6,*) 'mismatch in orbital information in hk_file and hu_file'
           i_terminate = .true.
        endif

     endif
     call MPI_Bcast(i_terminate, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     if (i_terminate) then
        if (rank .eq. 0) then
           write(6,*) 'Problem with number of interacting orbitals - stopping'
        endif
        stop
     endif

     if (rank .eq. 0) then
        do i_dummy_1 = 1, N_o_int*N_o_int
           read(25,*) i1, i2, uuij, ujij
           uu_temp(i1,i2) = uuij
           uj_temp(i1,i2) = ujij
        enddo
     endif

     m_size = max_orb_per_int_atom * max_orb_per_int_atom
     call MPI_Bcast(uu_temp, m_size, MPI_DOUBLE_PRECISION, &
          0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(uj_temp, m_size, MPI_DOUBLE_PRECISION, &
          0, MPI_COMM_WORLD, ierr)

     uup(i_a,:,:) = uu_temp
     uuj(i_a,:,:) = uj_temp

  enddo

  if (rank .eq. 0) then
     close(unit =25)
  endif

  allocate( h_so_tmp(0:2*max_orb_per_atom-1, 0:2*max_orb_per_atom-1))
  allocate( h_so(0:N_a-1, 0:2*max_orb_per_atom-1, 0:2*max_orb_per_atom-1) )
  
  if (rank .eq. 0) then
     open(unit=35, file=hso_file, status='old')
     read(35,*)
  endif
  
  do i_a = 0, N_a-1
     h_so_tmp = dcmplx(0.0d0, 0.0d0)
     
     if (rank .eq. 0) then
        read(35,*) so_flag, so_amp
        write(6,*) 'so_amp = ', so_amp
     endif
     call MPI_Bcast(so_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(so_amp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

     if (so_flag) then
        if (rank .eq. 0) then
           read(35,*) N_o_tmp
           write(6,*) 'N_o_tmp = ', N_o_tmp 
           do i_dummy_1 = 1, 2*N_o_tmp*2*N_o_tmp
              read(35,*) is, isp, re_hso, im_hso
              h_so_tmp(is,isp) = dcmplx(re_hso, im_hso)
           enddo
        endif
        m_size = (2*max_orb_per_atom)*(2*max_orb_per_atom)
        call MPI_Bcast(h_so_tmp, m_size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
        h_so_tmp = so_amp * h_so_tmp
     endif
     h_so(i_a, :, :) = h_so_tmp
  enddo
     
  if (rank .eq. 0) then
     close(unit=35)
  endif

  if (rank .eq. 0) then
     open(unit=45, file=hfield_file, status='old')
     read(45,*)
     read(45,*) prfld
     read(45,*)
  endif

  call MPI_Bcast(prfld, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  allocate ( h(0:nb-1, 3) )

  do i_a = 0, nb-1
     if (rank .eq. 0) then
        read(45,*) h_temp(1), h_temp(2), h_temp(3)
     endif
     call MPI_Bcast(h_temp, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     h(i_a,:) = h_temp
  enddo

  if (rank .eq. 0) then
     close(unit=45)
  endif
  
end subroutine read_hamiltonian

end MODULE Hamiltonian
