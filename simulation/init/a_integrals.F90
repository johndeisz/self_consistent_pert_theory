MODULE a_integrals
  USE CONSTANTS
  USE tau_epsilon_omega
  USE analytic_functions
  double complex, dimension(:,:,:,:,:), allocatable :: a_int

CONTAINS
subroutine generate_a_integrals
  IMPLICIT NONE

  double precision a(0:1,0:1,0:1), b(0:1,0:1,0:1)
  double precision fx, ex, nby, ey
  INTEGER ia, ib, ja, jb, i, l

  double complex part1(0:mp1), part2(0:mp1), part3(0:mp1), part4(0:mp1)

  !     Define constants appearing in q and r functions.

  allocate(a_int(0:1,0:1,0:1,0:1,0:mp1))

  do i = 0,1     
     a(0,0,i) = -0.5d0
     a(1,0,i) = -0.5d0
     a(0,1,i) = 0.5d0 / x(1,i)
     a(1,1,i) = -0.5d0 / x(1,i)
        
     b(0,0,i) = -0.5d0
     b(1,0,i) = -0.5d0
     b(0,1,i) = -0.5d0 / y(1,i)
     b(1,1,i) = 0.5d0 / y(1,i)
  enddo

  !     Evaluate the integrals for the Fourier transform of  Q(-tau) R(tau)

  do ia = 0, 1
     do ib = 0, 1

        fx = 1.0d0 / ( 1.0d0 + exp( -x(ia,ib) / t) )
        ex = exp( -x(ia,ib) / t)

        do ja = 0, 1
           do jb = 0, 1
 
              nby = 1.0d0 / ( exp( -y(ja,jb) / t ) - 1.0d0 )
              ey = exp( -y(ja,jb) / t )
              
              do l = 0, mp1

                 part1(l) = (ex + ey) * a(0,ia,ib) * b(0,ja,jb) / & 
                      ( dcmplx( 0.0d0, epsilon(l) ) +  ( x(ia,ib) - y(ja,jb) ) )

                 part2(l) = (ex*ey + 1.0d0) * &
                      a(1,ia,ib) * b(0,ja,jb) / & 
                      ( dcmplx( 0.0d0, epsilon(l) ) + & 
                      ( -x(ia,ib) - y(ja,jb) ) )

                 part3(l) = -( ex*ey + 1.0d0 ) * a(0,ia,ib) * b(1,ja,jb) / & 
                      ( dcmplx( 0.0d0,epsilon(l) ) + ( x(ia,ib) + y(ja,jb) ) )
 
                 part4(l) = -(ex + ey) * a(1,ia,ib) * b(1,ja,jb) / & 
                      ( dcmplx( 0.0d0, epsilon(l) ) + ( -x(ia,ib) + y(ja,jb) ) )

                 a_int(ia,ib,ja,jb,l) = fx * nby * &
                      ( part1(l) +  part2(l) + part3(l) + part4(l) )

              enddo
              
           enddo
        enddo

     enddo
  enddo


  return
end subroutine generate_a_integrals
end MODULE
