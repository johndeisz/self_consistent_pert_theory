MODULE tau_epsilon_omega
  USE CONSTANTS
  double precision, dimension(:), allocatable :: tau, epsilon, omega

CONTAINS
subroutine generate_tau_eps_omega()

  USE CONSTANTS
  IMPLICIT NONE

  include 'mpif.h'

  INTEGER rank, np, ierr

  double precision delta_tau
  INTEGER l

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

  delta_tau = (1.d0 / t ) / float(m)
  allocate (tau(0:mp1))
  allocate (epsilon(0:mp1))
  allocate (omega(0:mp1))

  do l = 0, mp1
     tau(l) = float(l + rank*mp) * delta_tau  
  enddo

  if (np .eq. 1) then

     do l = 0, m / 2 - 1
        epsilon(l) = float(2 * l + 1) * pi * t
        omega(l) = float(2 * l) * pi * t
     enddo
     do l = m/2, m1
        epsilon(l) = float(2 * (l-m) + 1) * pi * t
        omega(l) = float(2 * (l-m)) * pi * t
     enddo

  else 
     
     if ( rank .lt.  np/2 ) then

        do l = 0, mp1
           epsilon(l) = float(2 * (l+rank*mp) + 1) * pi * t
           omega(l) = float(2 * (l+rank*mp) ) * pi * t
        enddo
        
     else

        do l = 0, mp1
           epsilon(l) = float(2 * (l+rank*mp - m) + 1) * pi * t
           omega(l) = float(2 * (l+rank*mp - m) ) * pi * t           
        enddo
        
     endif
  endif

  return
end subroutine generate_tau_eps_omega
end MODULE
