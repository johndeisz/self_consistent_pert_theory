
  if (sigma_converged) then
     convergence_text=""
  else
     convergence_text=" NOT CONV."
  endif

  if (rank .eq. 0) then
     open(unit=25,file='converged_mu',status='unknown')
     write(25,*) mu
     close(25)

     if (.not. sigma_converged) then
        write(6,*) 'Iterations failed to converge.'
        write(6,*)
     endif

     write(6,*) 'Calculation parameters: '
     write(6,*) '-------------------------------------'
     write(6,*) 't (K) = ', t/kb
     write(6,*) 'Final mu = ', mu, convergence_text
     write(6,*) 'flux_x = ', flux(1)
     write(6,*) 'flux_y = ', flux(2)
     write(6,*) 'flux_z = ', flux(3)
     write(6,*) 'applied pair field = ', prfld
     write(6,*) 'Orbital dependent magnetic fields (T)'
     do ib = 0, nb - 1
        write(6,*) h(ib,:)/mub
     enddo
     write(6,*)
     write(6,*) 'interaction parameters'
     write(6,*) 'uu = ', uu, ' up = ', up, ' uj = ', uj
     write(6,*) 
     write(6,*) 'band structure parameters are at the top'
     write(6,*)
#ifdef SECOND_ORDER
     write(6,*) 'lcx = ', lcx,', lcy = ', lcy,', lcz =', lcz
     write(6,*) 'm = ', m
#endif
     write(6,*) 'llx = ', llx,',  lly = ',lly,',  llz = ',llz
     write(6,*)
     write(6,*)

#if defined (FLEX)
     write(6,*) 'Fluctuation exchange approximation results.'
#elif defined (THIRD_ORDER)
     write(6,*) 'Third order perturbation theory results.'
#elif defined (SECOND_ORDER)
     write(6,*) 'Second order perturbation theory results.'
#else
     write(6,*) 'First order perturbation theory results.'
#endif /* defined (FLEX) */

!!$C$$$#ifdef THIRD_ORDER
!!$C$$$        write(6,*) 'Max chi eigenvalue = ', overall_eigenvalue_max
!!$C$$$        write(6,*) 'index(0) = ', dominant_chi_index(0)
!!$C$$$        write(6,*) 'index(1) = ', dominant_chi_index(1)
!!$C$$$        write(6,*) 'index(2) = ', dominant_chi_index(2)
!!$C$$$        write(6,*) 'index(3) = ', dominant_chi_index(3)
!!$C$$$        write(6,*) 'matrix_index_00 = ', dominant_chi_eigenvector(0)
!!$C$$$        write(6,*) 'matrix_index_01 = ', dominant_chi_eigenvector(1)
!!$C$$$        write(6,*) 'matrix_index_02 = ', dominant_chi_eigenvector(2)
!!$C$$$        write(6,*) 'matrix_index_03 = ', dominant_chi_eigenvector(3)
!!$C$$$        write(6,*) 'matrix_index_04 = ', dominant_chi_eigenvector(4)
!!$C$$$        write(6,*) 'matrix_index_05 = ', dominant_chi_eigenvector(5)
!!$C$$$        write(6,*) 'matrix_index_06 = ', dominant_chi_eigenvector(6)
!!$C$$$        write(6,*) 'matrix_index_07 = ', dominant_chi_eigenvector(7)
!!$C$$$        write(6,*) 'matrix_index_08 = ', dominant_chi_eigenvector(8)
!!$C$$$        write(6,*) 'matrix_index_09 = ', dominant_chi_eigenvector(9)
!!$C$$$        write(6,*) 'matrix_index_10 = ', dominant_chi_eigenvector(10)
!!$C$$$        write(6,*) 'matrix_index_11 = ', dominant_chi_eigenvector(11)
!!$C$$$        write(6,*) 'matrix_index_12 = ', dominant_chi_eigenvector(12)
!!$C$$$        write(6,*) 'matrix_index_13 = ', dominant_chi_eigenvector(13)
!!$C$$$        write(6,*) 'matrix_index_14 = ', dominant_chi_eigenvector(14)
!!$C$$$        write(6,*) 'matrix_index_15 = ', dominant_chi_eigenvector(15)
!!$C$$$#endif /* THIRD_ORDER */

     write(6,*)
     write(6,*) 'Final density = ', density, convergence_text

     write(6,*)

     mag = 0.0d0
     magb = 0.0d0
     mag_energy = 0.0d0

     denb = 0.0d0
     pot_energy = 0.0d0

     do ib = 0, nb-1
        denb(ib) = real( g_tau0_local(4*ib+0,4*ib+0) + & 
             g_tau0_local(4*ib+1,4*ib+1) + 2.0d0 ) 
        pot_energy = pot_energy + ed(ib)*denb(ib)
     enddo

     do ib = 0, nb-1
        magb(ib,1) = real( g_tau0_local(4*ib+1,4*ib+0) + & 
             g_tau0_local(4*ib+0,4*ib+1) ) 
        mag = mag + magb(ib,1)
        mag_energy = mag_energy - & 
             mub*(gs/2.0d0)*h(ib,1)*magb(ib,1)
     enddo
     write(6,*) 'Mx_spin per atom/((g/2)*mub) = ', mag, convergence_text

     mag = 0.0d0
     do ib = 0, nb-1
        magb(ib,2) = real( cmplx(0.0d0,1.0d0)* &
             (-g_tau0_local(4*ib+1,4*ib+0) + &
             g_tau0_local(4*ib+0,4*ib+1)) ) 
        mag = mag + magb(ib,2)
        mag_energy = mag_energy - mub*(gs/2.0d0)*h(ib,2)*magb(ib,2)
     enddo
     write(6,*) 'My_spin per atom/((g/2)*mub) ', mag, convergence_text

     mag = 0.0d0
     do ib = 0, nb-1
        magb(ib,3) =  real( g_tau0_local(4*ib,4*ib) - & 
             g_tau0_local(4*ib+1,4*ib+1) )
        mag = mag + magb(ib,3)
        mag_energy = mag_energy - & 
             mub*(gs/2.0d0)*h(ib,3)*magb(ib,3)
     enddo
     write(6,*) 'Mz_spin per atom/((g/2)*mub) ', mag, convergence_text

     write(6,*) 'Density and magnetization by band'
     write(6,*) 'band  rho     mx       my      mz'
     do ib = 0, nb-1
        write(6,225) ib, denb(ib), magb(ib,1), magb(ib,2), magb(ib,3)
     enddo

     Lorb = cmplx(0.0d0, 0.0d0)
     
     call angular_matrices(Lvec)

     do ib = 0, nb-1
        do jb = 0, nb-1

           if (ib .eq. jb) then            
              Lorb = Lorb + Lvec(:,ib,ib) * (g_tau0_local(4*ib,4*ib) + &
                   g_tau0_local(4*ib+1,4*ib+1) + 2.0d0)
           else
              Lorb = Lorb + Lvec(:,ib,jb) * (g_tau0_local(4*jb,4*ib) + &
                   g_tau0_local(4*jb+1,4*ib+1))
           endif

        enddo
     enddo

     write(6,*) 'Lx/mub = ', Lorb(1)
     write(6,*) 'Ly/mub = ', Lorb(2)
     write(6,*) 'Lz/mub = ', Lorb(3)

225  format(i3,2x,d16.8,2x,d16.8,2x,d16.8,2x,d16.8)

     write(6,*)
 

  endif
 
  if (rank .eq. 0) then
     call pair_field_energy(g_tau0, prfld, psi, pair_energy)
     write(6,*) 'Pair field energy = ', pair_energy
  endif
        
  alpha = 1.0d0
  call  pair_wave(psi, g_tau0, alpha, m_psi)
  
  if (rank .eq. 0) then
     write(6,*) 'Pair amplitude = ', m_psi, convergence_text
     
     write(6,*)
     call kinetic_energy(g_tau0, tij, ed, kinetic)
     write(6,*) 'kinetic energy = ', kinetic

     so_energy = cmplx(0.0d0, 0.0d0)
     do ib =0, nb-1
        do is = 0, 1
           do jb = 0, nb-1
              do js = 0, 1

                 glocal = g_tau0_local(4*jb+js, 4*ib+is)
                 
                 if ( 4*ib+is .eq. 4*jb+js ) then
                    glocal = glocal + 1.0d0
                 endif

                 so_energy = so_energy + h_so(2*ib+is,2*jb+js)*glocal
                 
              enddo
           enddo
              
        enddo
     enddo
        
     write(6,*) 'Spin-orbit energy = ', real(so_energy), convergence_text
     write(6,*) 'Potential energy = ', pot_energy, convergence_text
     write(6,*) 'Note: magnetic energy neglects orbital terms'
     write(6,*) 'Magnetic energy = ', mag_energy, convergence_text 
     call hfock_energy(sigma1, g_tau0_local, e_hf)
     write(6,*) 'Hartree-fock energy = ', e_hf, convergence_text
  endif

  cor_energy = 0.0d0

#ifdef SECOND_ORDER
  ft_sigma = .false.
  call sigma_calc(rank, t, sigma, chi, g, g_mtau, &
       c_r, tau, epsilon, q_tau, q_epsilon, x, y, &
       r_tau, r_omega, a_int, gamma0_ph, &
       overall_eigenvalue_max,  dominant_chi_eigenvector, &
       dominant_chi_index, ft_sigma, d_r)
        
  call eval_tr_sigph_g(rank, tr_sig_g, sigma, g_mtau, t, &
       c_r, d_r, q_tau, q_mtau, r_tau, x, y)

  cor_energy = 0.25d0 * real(tr_sig_g)
      
  if (rank .eq. 0) then
     write(6,*) 'Correlation energy = ', cor_energy, convergence_text
  endif
#endif /* SECOND_ORDER */
  if (rank .eq. 0) then
     write(6,*) 'Total energy = ', kinetic + pot_energy + & 
          mag_energy + real(so_energy) + pair_energy + &
          e_hf + cor_energy, convergence_text
  endif

  if (rank .eq. 0) then
     if (m_psi .gt. 1.0d-12) then
        call analyze_psi(psi, tij)
     endif
  endif

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
