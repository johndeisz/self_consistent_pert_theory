#include "../convert.F90"

subroutine init_environ(rank, np, starttime)

  USE CONSTANTS
  IMPLICIT NONE
       
  include 'mpif.h'

  Real starttime
  integer rank, np
  integer ierr
  logical flag_np

  double precision test_np

  call cpu_time(starttime)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

  if (rank .eq. 0) then
#if defined (FLEX)
     print *, 'Code compiled for the fluctuation exchange approximation.'
#elif defined (THIRD_ORDER)
     print *, 'Code compiled for third order perturbation theory.'
#elif defined (SECOND_ORDER)
     print *, 'Code compiled for second order perturbation theory.'
#else
     print *, 'Code compiled for first order perturbation theory.'
#endif /* defined (FEA) */

     print *
     print *, 'Running on', np, ' processes.'
     test_np = dlog(dfloat(np))/dlog(2.0d0)
     test_np = mod(test_np, 1.0d0)
     flag_np = .false.
     if (abs(test_np) .gt. 1.0d-6) then
       print *, 'Number of mpi processes must equal 2^n'
       print *, 'Stopping'
       flag_np = .true.
     endif
  endif

  call MPI_Bcast(flag_np, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  
  if (flag_np) then
     call MPI_Finalize(ierr)
     stop
  endif

  return
end subroutine init_environ
	 
