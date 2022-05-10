#include "convert.F90"

program multiband_flex_dca

  USE CONSTANTS
 ! USE bare_dispersion

  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

#include "main_defs.F90"

  call init_environ(rank, size, start_time)

  ! Gather all the data from the input file
! call readin(t, flux, prfld, h, target_density, density_tol, mu, uu,
  ! up, uj, ed, tij, prfld_pert, h_pert, v_pert, h_so, read_input,
  ! sigma_input_file, write_output, sigma_output_file,
  ! max_pade_order, sigma_tol,  max_it, alpha, alpha_scheme)

  ! Put temperature in eV
! t = kb*t
  ! Multiply H by mu_B to get eV
! h = mub*h
! h_pert = mub*h_pert 

  !     Bare bandstructure and vertex
! ek_min = ek_minimum(tij, ed)
! if (rank .eq. 0) then
! write(6,*) "ek_min = ", ek_min
! endif
! call gamma0_define(rank, gamma0_ph, uu, up, uj)

#ifdef SECOND_ORDER
  !     generate the tau and epsilon matrices
! call generate_tau_eps_omega(t, tau, epsilon, omega)
  !     Create the analytic functions
! call analytic_functions(t, tau, epsilon, omega, x, y, q_tau,
  ! q_mtau, q_epsilon, r_tau, r_omega)
! call a_integrals(t, x, y, epsilon, a_int)
#endif

  !     initialize sigma
! call init_sigma1(uu, target_density, sigma1)
#ifdef SECOND_ORDER
! sigma = cmplx(0.0d0, 0.0d0)
#endif

  !     Initialize a random pair wave function.
! call init_pair_wave(psi)
! isignv(0) = 0
! isignv(1:3) = -1
! call psi_transform(psi, isignv)

! if (read_input) then
#ifdef SECOND_ORDER
!    call sigma_input(sigma_input_file, uu, sigma1, psi, sigma, mu,
  !    epsilon)
#else
!    call sigma_input(sigma_input_file, uu, sigma1, psi)
#endif
! endif

#ifdef SECOND_ORDER
! iteration = 0
! call effective_field(iteration, v_pert, h, h_pert, prfld,
  ! prfld_pert, v_pert_eff, h_eff, prfld_eff)
! call discontinuities(tij, ed, v_pert_eff, psi, h_eff, prfld_eff,
  ! mu, sigma1, h_so, delta_g_r, delta_g_k, delta_gp_r, delta_gp_k)
! call parameter_init(x, c_r, delta_g_r, delta_gp_r, y)
#endif

  !    Prepare for self-consistency loop.
! mu_old = ek_min - 1.0d0
! density_old = 0.01d0
! iteration = 1
! density_converged = .false.
! density_iteration = 0
! ft_sigma = .true.

! mu_min = ek_min - 30.0d0
! mu_max = 20.0d0

!!$  do while ( (iteration .le. max_it) .and. (.not.
  !!density_converged) )
!!$
!!$     sigma_converged = .false.
!!$
!!$     call cpu_time(last_it_time)
!!$
!!$     if (rank .eq. 0) then
!!$        write(6,*)
!!$        write(6,fmt=500,advance='NO')
!!$        write(6,*) ' Iter    Mix      %  Converged     ', $
  !!'   Pair_amp    iter. time'
!!$        write(6,fmt=500,advance='NO')
!!$#ifdef SECOND_ORDER
!!$        write(6,*) '            Sigma(e) ', ' Sigma1 '
!!$#else
!!$        write(6,*) '                          ', ' Sigma1'
!!$#endif
!!$        write(6,fmt=500,advance='NO')
!!$        write(6,*)
  !!'----------------------------------------------', $
  !!'-----------------------'
!!$     endif
!!$
!!$     do while ( (iteration .le. max_it) .and. (.not.
  !!sigma_converged) )
!!$
!!$        if (rank .eq. 0) then
!!$           write(6,fmt=600,advance='NO') iteration, alpha
!!$        endif
!!$
!!$        sigma_converged = .true.
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'start iteration = ', iteration
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        call effective_field(iteration, v_pert, h, h_pert, prfld,
!!$             prfld_pert, v_pert_eff, h_eff, prfld_eff)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'effective field'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$#ifdef SECOND_ORDER
!!$        call discontinuities(tij, ed, v_pert_eff, psi, h_eff,
  !!prfld_eff, $             mu, sigma1, h_so, delta_g_r, delta_g_k,
  !!delta_gp_r, delta_gp_k)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'discontinuities'
!!$     close(unit=9)
!!$  endif
!!$
!!$          
!!$        sigma_old = sigma

! Include these lines to supress off-diagonal terms
!        do nu1 = 0, nb-1
!         do nu2 = 0, nb-1
!          sigma(4*nu1, 4*nu2+2,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1, 4*nu2+3,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1+1, 4*nu2+2,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1+1, 4*nu2+3,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1+2, 4*nu2,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1+2, 4*nu2+1,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1+3, 4*nu2,:,:) = cmplx(0.0d0, 0.0d0)
!          sigma(4*nu1+3, 4*nu2+1,:,:) = cmplx(0.0d0, 0.0d0)
!        enddo
!       enddo

          
!!$        call dyson(rank, g, q_tau, q_epsilon, tij, ed, $
  !!v_pert_eff, psi, h_eff, prfld_eff, mu, sigma1, h_so, $
  !!sigma, epsilon, t)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'dyson'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        call calc_g_tau0_2nd(rank, g_tau0, q_tau, q_epsilon, tij,
  !!ed, $             v_pert_eff, psi, h_eff, prfld_eff, mu, sigma1,
  !!h_so, $             sigma, epsilon, t)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'calc_g_tau0'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        g_mtau = g
!!$
!!$        call g_rtau(rank, g, t, c_r, q_epsilon, q_tau, tau)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'calc_g_rtau'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$
!!$        call g_minus_tau(rank, g_mtau, t, c_r, q_epsilon, q_tau,
  !!tau)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'g_minus_tau'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        call symmetrize(rank, g, g_mtau)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'symmetrize'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        call green_parameter(rank, g, t, x, c_r, delta_g_r,
  !!delta_gp_r)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'green_parameter'
!!$     close(unit=9)
!!$  endif


!!$        call sigma_calc(rank, t, sigma, chi, g, g_mtau, c_r, $
  !!tau, epsilon, q_tau, q_epsilon, x, y, r_tau, r_omega, $
  !!a_int, gamma0_ph, overall_eigenvalue_max, $
  !!dominant_chi_eigenvector, dominant_chi_index, ft_sigma, d_r)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'sigma_calc'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        if ( (alpha_scheme .eq. 2) .and. (iteration .gt. 1) ) then
!!$
!!$           if (rank .eq. 0) then
!!$              call calc_new_alpha(alpha, delta_sigma_e0_old,
  !!sigma, sigma_old)
!!$           endif
!!$#ifdef USE_MPI
!!$           call MPI_Bcast(alpha, 1, MPI_REAL, 0, MPI_COMM_WORLD,
  !!ierr)
!!$#endif /* USE_MPI */
!!$
!!$        endif
!!$
!!$        call convergence_test(sigma_converged, rank, iteration, $
  !!sigma, sigma_old, sigma_tol, last_it_time)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'convergence test'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        if ( (iteration .gt. 1) .or. read_input) then
!!$           sigma = alpha*sigma + (1.0d0 - alpha)*sigma_old
!!$        endif
!!$
!!$        g_tau0_local = g(:,:,0,0)    
!!$#else
!!$        call calc_g_tau0(tij, ed, v_pert_eff, psi, h_eff,
  !!prfld_eff, mu, $             sigma1, h_so, t, g_tau0,
  !!g_tau0_local)
!!$#endif
!!$        call sigma_first(sigma1, sigma1_old, delta_sigma1, $
  !!delta_sigma1_old, alpha, g_tau0_local, gamma0_ph, sigma_tol, $
  !!sigma_converged, alpha_scheme, iteration)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'sigma first'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$        iteration = iteration + 1
!!$        call  pair_wave(psi, g_tau0, alpha, m_psi)
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'pair wave'
!!$     close(unit=9)
!!$  endif


!!$        call cpu_time(this_it_time)
!!$        if (rank .eq. 0) then
!!$           write(6,*) "  ", m_psi, this_it_time - last_it_time
!!$        endif
!!$        last_it_time = this_it_time
!!$
!!$#ifdef FLEX
!!$        if (rank .eq. 0) then
!!$           if (overall_eigenvalue_max .gt. 1.0d0) then
!!$              write(6,*) 'largest t-matrix eigenvalue = ', $
  !!overall_eigenvalue_max
!!$           endif
!!$        endif
!!$#endif
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'iteration completed'
!!$     close(unit=9)
!!$  endif
!!$
!!$
!!$     enddo
!!$
!!$     if (rank .eq. 0) then
!!$          
!!$        density = 0.0d0
!!$
!!$        do ib = 0, nb-1
!!$          
!!$           density = density + $
  !!(real(g_tau0_local(4*ib+0,4*ib+0)) + 1.0d0) + $
  !!(real(g_tau0_local(4*ib+1,4*ib+1)) + 1.0d0)
!!$
!!$        enddo
!!$
!!$        write(6,*)
!!$        write(6,*)
!!$        write(6,*) "mu = ", mu, "  density = ", density
!!$        write(6,*)
!!$
!!$        if (abs (density - target_density) .lt. density_tol) then
!!$
!!$           density_converged = .true.
!!$
!!$        else
!!$
!!$           if (density .gt. target_density) then
!!$              mu_max = mu
!!$           else
!!$              mu_min = mu
!!$           endif
!!$
!!$           dn_dmu = (density - density_old) / (mu - mu_old)
!!$           delta_mu = (target_density - density) / (dn_dmu)
!!$
!!$           mu_old = mu
!!$           density_old = density
!!$
!!$           if (density_iteration .eq. 0) then
!!$              mu = mu + 0.7d0 * delta_mu
!!$           else
!!$              mu = mu + delta_mu
!!$           endif
!!$
!!$           if ( (mu .gt. mu_max) .or. (mu .lt. mu_min) ) then
!!$              mu = 0.5d0 * (mu_min + mu_max)
!!$           endif
!!$
!!$        endif
!!$
!!$
!!$
!!$     endif
!!$
!!$#ifdef USE_MPI
!!$     call MPI_Bcast(mu, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
!!$     call MPI_Bcast(density_converged, 1, MPI_LOGICAL, 0, $
  !!MPI_COMM_WORLD, ierr)
!!$#endif /* USE_MPI */
!!$
!!$     density_iteration = density_iteration + 1
!!$
!!$  if (rank .eq. 0) then
!!$     open(unit=9,file='my_error_file',status='old',
  !!access='append')
!!$     write(9,*) 'density iteration completed'
!!$     close(unit=9)
!!$  endif
!!$
!!$  enddo
!!$
!!$  if (write_output) then
!!$
!!$     if (rank .eq. 0) then
!!$
!!$        write(6,*)
!!$        write(6,*) "Sigma output file written to ",
  !!sigma_output_file
!!$        write(6,*)
!!$
!!$        open(unit=40, file=sigma_output_file, status='unknown')
!!$
!!$        call sigma1_out(t, mu, flux, prfld, ed, h, tij, $
  !!uu, up, uj, psi, h_so, sigma1)
!!$
!!$     endif
!!$#ifdef SECOND_ORDER
!!$     call sigma_out(rank, max_pade_order, sigma_output_file, t,
  !!sigma)
!!$#endif /* SECOND_ORDER */
!!$
!!$     if (rank .eq. 0) then
!!$        close(unit=40)
!!$     endif
!!$
!!$  else
!!$
!!$     if (rank .eq. 0) then
!!$        write(6,*)
!!$        write(6,*) "Sigma output file not written"
!!$        write(6,*)
!!$     endif
!!$
!!$  endif
!!$
!!$
  call cpu_time(end_time)

  if (rank .eq. 0) then
     write(6,*)
     write(6,*) 'Execution time = ', end_time-start_time, ' seconds.'
  endif

#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif /* USE_MPI */

500 format('xx')
600 format('xx ', I5, '  ', f5.3,'  ')
700 format('       ',f10.8)

  stop
      
end program multiband_flex_dca
