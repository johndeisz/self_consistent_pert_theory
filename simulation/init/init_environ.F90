#include "../convert.F90"

subroutine init_environ(rank, np, starttime)

  USE CONSTANTS
  IMPLICIT NONE
       
#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  Real starttime
  integer rank, np
  integer ierr

  call cpu_time(starttime)

#ifdef USE_MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#else
  rank = 0
  np = 1
#endif /* USE_MPI */

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
     print *, 'Note: number of processes must equal 2^n where n is an integer'

  endif

  return
end subroutine init_environ
	 
