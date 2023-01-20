MODULE analytic_functions
  USE CONSTANTS
  USE tau_epsilon_omega
  double precision :: x(0:1,0:1), y(0:1,0:1)
  double precision, dimension(:,:,:), allocatable :: q_tau, q_mtau, r_tau
  double complex, dimension(:,:,:), allocatable :: q_epsilon
  double complex, dimension(:,:,:), allocatable :: r_omega

CONTAINS
subroutine generate_analytic_functions()
  IMPLICIT NONE
  INTEGER j, l

  !     Define the x and y arrays
  x(0,0) = 1.5d0
  y(0,0) = 1.25d0
  x(0,1) = 2.5d0
  y(0,1) = 2.25d0
  x(1,0) = 2.0d0  
  y(1,0) = 1.75d0
  x(1,1) = 3.0d0
  y(1,1) = 2.75d0

  allocate(q_tau(0:1,0:1,0:mp1))
  allocate(q_mtau(0:1,0:1,0:mp1))
  allocate(r_tau(0:1,0:1,0:mp1))
  allocate(q_epsilon(0:1,0:1,0:mp1))
  allocate(r_omega(0:1,0:1,0:mp1))

  do l = 0, mp1
     do j = 0, 1

        q_epsilon(0,j,l) = 0.5d0 * &
             ( 1.0d0 / (cmplx(0.0d0, epsilon(l)) - x(0,j) )  + &
             1.0d0 / (cmplx(0.0d0, epsilon(l)) + x(0,j) ) )
      
        q_epsilon(1,j,l) = (1.0d0 / (2.0d0 * x(1,j)) ) * &
             ( 1.0d0 / (cmplx(0.0d0, epsilon(l)) - x(1,j) ) - &
             1.0d0 / (cmplx(0.0d0, epsilon(l)) + x(1,j) ) )
      
        q_tau(0,j,l) = -0.5d0 * &
             ( exp( -x(0,j)*tau(l) ) / ( exp(-x(0,j) / t) + 1 )  + &
             exp(-x(0,j) * ( (1.0d0 / t) - tau(l))) &
             / (exp(-x(0,j) / t) + 1.0d0 ) )
      
        q_tau(1,j,l) = - ( 1.0d0 / (2.0d0 * x(1,j)) ) * &
             ( exp(-x(1,j)*tau(l) ) / ( exp( -x(1,j) / t ) + 1.0d0 ) - &
             exp( -x(1,j) * ( (1.0d0 / t) - tau(l)) ) &
             / (exp(-x(1,j) / t) + 1.0d0 ) )
      
        q_mtau(0,j,l) = 0.5d0* ( exp( -x(0,j) * (1.0d0 / t - tau(l) )) &
             / ( exp(-x(0,j) / t) + 1.0d0 )  + &
             exp(-x(0,j) * tau(l) )  / ( exp(-x(0,j) / t) + 1.0d0) )
      
        q_mtau(1,j,l) = (1.0d0 / (2.0d0*x(1,j)) ) * &
             ( exp(-x(1,j)*(1.0d0 / t - tau(l)) ) &
             / ( exp(-x(1,j) / t) + 1.0d0)  -  exp( -x(1,j)*tau(l) ) &
             / ( exp(-x(1,j) / t) + 1.0d0) )

        r_omega(0,j,l) = 0.5d0 * &
             ( 1.0d0 / ( cmplx( 0.0d0, omega(l) ) - y(0,j) ) + &
             1.0d0 / ( cmplx(0.0d0, omega(l) ) + y(0,j) ) )
    
        r_omega(1,j,l) = (1.0d0 / (2.0d0 * y(1,j) ) ) * &
             ( 1.0d0 / ( cmplx(0.0d0, omega(l)) - y(1,j) ) - &
             1.0d0 / ( cmplx(0.0d0, omega(l) ) + y(1,j) ) )
    
        r_tau(0,j,l) = -0.5d0 * ( exp( -y(0,j) * tau(l) ) &
             - exp( -y(0,j) * ( (1.0d0/ t) - tau(l) ) ) ) &
             / ( 1.0d0 - exp(-y(0,j) / t) )
    
        r_tau(1,j,l) = -( 1.0d0 / (2.0d0 * y(1,j)) ) * &
             ( exp(-y(1,j) * tau(l) ) &
             + exp(-y(1,j) * ( (1.0d0 / t) - tau(l) ) ) ) &
             / (1.0d0 - exp(-y(1,j) / t ) ) 
    
     enddo
  enddo

  return
end subroutine generate_analytic_functions
end MODULE
