subroutine gamma0_define()

  USE CONSTANTS
  use hamiltonian
  
  IMPLICIT NONE

  integer :: nba, msize
  integer :: na1, na2, na3, na4
  integer :: nu1, nu2, nu3, nu4
  integer :: is1, is2, is3, is4
  integer :: ind1, ind2, ind3, ind4
  integer :: ind1p, ind2p, ind3p, ind4p
  integer :: i_a, iba, jba
  
  double complex, dimension(:,:,:,:), allocatable :: gamma0

  double complex, dimension (1:3, 0:1, 0:1) :: pauli
  double precision, dimension (0:1, 0:1) :: dirac

  double precision, dimension(:), allocatable :: uum
  double precision, dimension(:,:), allocatable :: upm, ujm

  double precision :: s0, c0

  ! First, initialize the gamma0, gamma0_ph and related arrays

  msize = max_orb_per_int_atom
  allocate ( uum(0:msize-1) )
  allocate ( upm(0:msize -1, 0: msize-1) )
  allocate ( ujm(0:msize -1, 0: msize-1) )
  
  msize = 4*max_orb_per_int_atom
  allocate( gamma0(0:msize-1, 0:msize-1, 0:msize-1, 0:msize-1) )
  gamma0 = dcmplx(0.0d0, 0.0d0)

  msize = 16 * max_orb_per_int_atom * max_orb_per_int_atom
  allocate ( gamma0_ph(0:N_int_a-1, 0:msize-1, 0:msize-1 ) )
  gamma0_ph = dcmplx(0.0d0, 0.0d0)


  dirac = 0.0d0
  do is1 = 0, 1
     dirac(is1,is1) = 1.0d0
  enddo
  
  pauli(1,0,0) = 0.0d0
  pauli(1,1,0) = 1.0d0
  pauli(1,0,1) = 1.0d0
  pauli(1,1,1) = 0.0d0

  pauli(2,0,0) = 0.0d0
  pauli(2,1,0) = cmplx(0.0d0,1.0d0)
  pauli(2,0,1) = cmplx(0.0d0,-1.0d0)
  pauli(2,1,1) = 0.0d0

  pauli(3,0,0) = 1.0d0
  pauli(3,1,0) = 0.0d0
  pauli(3,0,1) = 0.0d0
  pauli(3,1,1) = -1.d0

  do i_a = 0, N_int_a-1

     nba = N_o(i_a)
     
     do iba = 0, nba-1
        uum(iba) = uup(i_a, iba, iba)
        do jba = 0, nba - 1
           if (iba .ne. jba) then
              upm(iba,jba) = uup(i_a,iba,jba)
              ujm(iba,jba) = uuj(i_a,iba,jba)
           else
              upm(iba,jba) = 0.0d0
              ujm(iba,jba) = 0.0d0
           end if
        enddo
     enddo

     do nu1 = 0, nba-1
        do nu2 = 0, nba-1
           do nu3 = 0, nba-1
              do nu4 = 0, nba-1

                 s0 = 0.0d0
                 c0 = 0.0d0

                 if (nu1 .eq. nu2) then
                    if (nu3 .eq. nu2) then
                       if (nu4 .eq. nu3) then

                          s0 = uum(nu1)
                          c0 = uum(nu1)

                       endif
                    else
                       if (nu4 .eq. nu3) then

                          s0 = ujm(nu1,nu3)
                          c0 = ujm(nu1,nu3)

                       endif
                    endif

                 else

                    if (nu1 .eq. nu3) then
                       if (nu2 .eq. nu4) then

                          s0 = ujm(nu1,nu2)
                          c0 = -ujm(nu1,nu2) + 2.0d0*upm(nu1,nu2)
                    
                       endif
                    else
                       if (nu1 .eq. nu4) then
                          if (nu2 .eq. nu3) then

                             s0 = upm(nu1,nu2)
                             c0 = 2.0*ujm(nu1,nu2) - upm(nu1,nu2)

                          endif
                       endif
                    endif
                 endif
                    
                 do is1 = 0,1
                    ind1 = 4*nu1 + is1

                    do is2 = 0,1
                       ind2 = 4*nu2 + is2

                       do is3 = 0,1
                          ind3 = 4*nu3 + is3

                          do is4 = 0,1
                             ind4 = 4*nu4 + is4

                             gamma0(ind1,ind2,ind3,ind4) = -0.5d0 * s0 * &
                                  (pauli(1,is1,is3)*pauli(1,is2,is4) + &
                                  pauli(2,is1,is3)*pauli(2,is2,is4) + &
                                  pauli(3,is1,is3)*pauli(3,is2,is4)) + &
                                  0.5d0 * c0 * dirac(is1,is3)*dirac(is2,is4)

                             
                          enddo
                       enddo
                    enddo
                 enddo

              enddo
           enddo
        enddo
     enddo

     do nu1=0,nba-1
        do nu2=0, nba-1
           do nu3=0, nba-1
              do nu4=0,nba-1
              
                 do is1 = 0,1
                    ind1 = 4*nu1 + is1
                    ind1p = ind1 + 2

                    do is2 = 0,1
                       ind2 = 4*nu2 + is2
                       ind2p = ind2 + 2

                       do is3 = 0,1
                          ind3 = 4*nu3 + is3
                          ind3p = ind3 + 2

                          do is4 = 0,1
                             ind4 = 4*nu4 + is4
                             ind4p = ind4 + 2

                             gamma0(ind1p,ind2,ind3p,ind4) = &
                                  -gamma0(ind2,ind3,ind4,ind1)

                             gamma0(ind1,ind2p,ind3,ind4p) = &
                                  -gamma0(ind1,ind4,ind3,ind2)

                             gamma0(ind1,ind2p,ind3p,ind4) = &
                                  gamma0(ind1,ind3,ind4,ind2)

                             gamma0(ind1p,ind2,ind3,ind4p) = &
                                  gamma0(ind2,ind4,ind3,ind1)

                             gamma0(ind1p,ind2p,ind3p,ind4p) = &
                                  gamma0(ind4,ind3,ind2,ind1)

                          enddo
                       enddo
                    enddo
                 enddo

              enddo
           enddo
        enddo
     enddo

     do na1 = 0, 4*nba-1
        do na2 = 0, 4*nba-1
           do na3 = 0, 4*nba-1
              do na4 = 0, 4*nba-1

                 gamma0_ph(i_a, 4*nba*na1+na3, 4*nba*na4+na2) = gamma0(na1,na2,na3,na4)

              enddo
           enddo
        enddo
     enddo

  enddo

  deallocate ( uum )
  deallocate ( ujm )
  deallocate ( upm )
  deallocate ( gamma0 )
  
  return
end subroutine gamma0_define
