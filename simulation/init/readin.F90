subroutine readin(target_density, density_tol, mu, read_input, sigma_input_file, &
     write_output, sigma_output_file, max_pade_order, &
     sigma_tol, max_it, alpha, alpha_scheme )

  USE CONSTANTS
  USE hamiltonian
  IMPLICIT NONE

  include 'mpif.h'

  logical flag_nl
  double precision test_nl

  integer :: il, i_a, i_b, ib
  integer :: i_o, i_dummy_1, i_dummy_2, ind
  
! Parameters to be read and returned
  double precision :: target_density, density_tol, mu
  LOGICAL read_input, write_output
  CHARACTER*128 sigma_input_file, sigma_output_file
  INTEGER max_pade_order

  double precision :: sigma_tol, alpha
  INTEGER :: alpha_scheme, max_it

  double precision :: h_pert_amp(1:3)
  double precision :: v_pert_amp
  ! MPI variables
  INTEGER rank
  INTEGER ierr
  INTEGER icode

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  if (rank .eq. 0) then
     read(5,*) 
     read(5,*) hk_file
     read(5,*) hu_file
     read(5,*) hso_file
     read(5,*) hfield_file
     read(5,*)
     read(5,*)
     read(5,*) Nl(1), Nl(2), Nl(3)
  end if

  if (rank .eq. 0) then
     flag_nl = .false.
     do il=1,3
        print *, 'il, Nl(il) = ', il, ' ', Nl(il)
        test_nl = dlog(dfloat(Nl(il)))/dlog(2.0d0)
        test_nl = mod(test_nl, 1.0d0)
        if (abs(test_nl) .gt. 1.0d-6) then
           flag_nl = .true.
        endif
     enddo
     if (flag_nl) then
        print *, 'Each lattice dimension must equal 2^n'
        print *, 'Stopping'
     endif
  endif

  call MPI_Bcast(flag_nl, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  
  if (flag_nl) then
     call MPI_Finalize(ierr)
     stop
  endif
  call MPI_Bcast(Nl, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call read_hamiltonian()

  if (rank .eq. 0) then
     read(5,*)
     read(5,*) 
     read(5,*) t

     read(5,*) target_density, density_tol
     read(5,*) mu
        
     write(6,*) "Basic Thermodynamic Parameters" 
     write(6,*) "temperature in K = ", t  

     write(6,*) "target electron density = ", target_density
     write(6,*) "density tolerance = ", density_tol
     write(6,*) "default initial chemical potential = ", mu

     open (unit=20, file='converged_mu',status='old', iostat=icode)

     if (icode .ne. 0) then
        write (6,*) 'File converged_mu not found, using default value for mu.'
     else 
        write (6,*) 'File converged_mu found.'
        read (20,*) mu
        write (6,*) 'Using mu = ',mu
     endif
     close (20)
     write(6,*) ' '

  endif

  call MPI_Bcast(t, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
 
  call MPI_Bcast(target_density, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(density_tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mu, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  !------------------------Sigma Input and Output --------------------------
      
  if (rank .eq. 0) then
     read(5,*)
     read(5,*)
     read(5,*) read_input, sigma_input_file
     read(5,*) write_output, sigma_output_file
     read(5,*) max_pade_order

     write(6,*) "Sigma input and output"
     write(6,*) "read sigma input file = ", read_input
     write(6,*) "Sigma input file = ", sigma_input_file
     write(6,*) "write sigma input file = ", write_output
     write(6,*) "Sigma output file = ", sigma_output_file
     write(6,*) "Maximum pade order for output sigma = ", max_pade_order    
  endif
  
  call MPI_Bcast(read_input, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(write_output, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_pade_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


!----------------------------FEA parameters ------------------------------

  if (rank .eq. 0) then

     read(5,*)
     read(5,*)
     read(5,*) sigma_tol,  max_it
     read(5,*) alpha, alpha_scheme
     
     write(6,*) "Flex parameters"
     write(6,*) "Tolerance in sigma = ", sigma_tol
     write(6,*) "Maximum iterations = ", max_it 
     write(6,*) "alpha = ", alpha, " alpha scheme = ", alpha_scheme
     
  endif

  call MPI_Bcast(sigma_tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(alpha, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(alpha_scheme, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
  call MPI_Bcast(prfld, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(h, 3*nb, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  
  !-----Perturbations: applied and relaxed over initial 10 iterations -------

  allocate ( h_pert(0:nb-1,3) )
  allocate ( v_pert(0:nb-1) )
  
  if (rank .eq. 0) then
     read(5,*)
     read(5,*)
     read(5,*) prfld_pert
     read(5,*) h_pert_amp
     read(5,*) v_pert_amp

     write(6,*) "Preturbations: applied and relaxed over initial 10 iterations"
     write(6,*) "Artificial pair field = ", prfld_pert
     write(6,*) "Random orb.-dependent mag. field amp (T) ", h_pert_amp
     do ib = 0, nb-1
        h_pert(ib,1) = h_pert_amp(1) * (rand()-0.5d0)
        h_pert(ib,2) = h_pert_amp(2) * (rand()-0.5d0)
        h_pert(ib,3) = h_pert_amp(3) * (rand()-0.5d0)
     enddo
     write(6,*) "Random orbital potential amplitude  = ", v_pert_amp
     do ib = 0, nb-1
        v_pert(ib) = v_pert_amp * (rand() - 0.5d0)
     enddo
  endif

  call MPI_Bcast(prfld_pert, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(h_pert, 3*nb, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(v_pert, nb, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

300 format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')

  return
end subroutine readin
