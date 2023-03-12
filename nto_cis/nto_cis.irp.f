program nto_cis
  implicit none
  BEGIN_DOC
! This program produces the natural transition orbitals (NTOs) from the state-averaged one-body transition density matrix 
! between the close-shell HF ground state and CIS excited state(s).
! The NTOs are saved into the EZFIO file, whereas the NTO excitation amplitudes are printed in the output file.
  END_DOC

  integer :: i, a, k, p


  if( elec_alpha_num .ne. elec_beta_num ) stop 'nto_cis only makes sense for systems with an even number of electrons'

  if( N_states < 2 ) stop 'nto_cis needs at least two states'

  double precision :: sum_state_average_weight_excited
  double precision :: sqrt_renormalized_state_average_weight(N_states)
  sum_state_average_weight_excited = sum( state_average_weight(2:N_states) )
  do k=2,N_states
    sqrt_renormalized_state_average_weight(k) = sqrt( state_average_weight(k) / sum_state_average_weight_excited )
  end do

  double precision :: one_e_tdm_mo_alpha(val_occ_num,vir_num)
  double precision :: one_e_tdm_mo_beta(val_occ_num,vir_num)
  one_e_tdm_mo_alpha = 0.0d0
  one_e_tdm_mo_beta  = 0.0d0

  do p = 2, N_det 
    integer          :: exc(0:2,2,2)
    integer          :: degree
    double precision :: phase
    call get_excitation(psi_det(N_int,:,1),psi_det(N_int,:,p),exc,degree,phase,N_int)
    integer :: h1,h2,p1,p2,s1,s2
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    if( degree.ne.1 ) stop 'degree different from 1. Only CIS wave functions are allowed in nto_cis'
    i = h1 - n_core_orb
    a = p1 - occ_num
    if( s1==1 ) then
      do k=2,N_states
        one_e_tdm_mo_alpha(i,a) += psi_coef(p,k) * sqrt_renormalized_state_average_weight(k)
      end do
    else if( s1==2 ) then
      do k=2,N_states
        one_e_tdm_mo_beta(i,a)  += psi_coef(p,k) * sqrt_renormalized_state_average_weight(k)
      end do
    end if
  end do

  double precision :: one_e_tdm_mo(val_occ_num,vir_num)
  one_e_tdm_mo = one_e_tdm_mo_alpha + one_e_tdm_mo_beta

  double precision :: trace_one_e_tdm_mo
  trace_one_e_tdm_mo = 0.0d0
  do a=1,vir_num
    do i=1,val_occ_num
      trace_one_e_tdm_mo += one_e_tdm_mo(i,a)**2
    end do
  end do
  trace_one_e_tdm_mo = trace_one_e_tdm_mo / 2.0d0
  write(*,*) 'Trace of the 1-body transition density matrix = ', trace_one_e_tdm_mo


! Perform a SVD decomposition of the one-body transition density matrix

  double precision :: nto_hole(val_occ_num,val_occ_num)
  double precision :: nto_particle_t(vir_num,vir_num)
  double precision :: nto_particle(vir_num,vir_num)
  double precision :: nto_amplitude(val_occ_num) ! It should be the smaller between val_occ_num and vir_num, which is almost always val_occ_num

  call svd_all(one_e_tdm_mo,val_occ_num,nto_hole,val_occ_num,nto_amplitude,nto_particle_t,vir_num,val_occ_num,vir_num)

  call dtranspose(nto_particle_t,vir_num,nto_particle,vir_num,vir_num,vir_num)

! logical :: ok
! call check_orthornormality_rotation_matrix(nto_hole,val_occ_num,ok)
! write(*,*) 'hole ', ok
! call check_orthornormality_rotation_matrix(nto_particle,vir_num,ok)
! write(*,*) 'part ', ok


! Rotate the intial molecular orbitals to the NTOs

! Rotate the occupied orbitals to the hole NTO
  call update_molecular_orbitals(mo_coef(:,n_core_orb+1:occ_num),ao_num,val_occ_num,nto_hole)

! Arrange the hole NTO in order of increasing NTO excitation amplitude
  nto_hole = 0.0d0
  do i=1,val_occ_num
    nto_hole(i,val_occ_num+1-i) = 1.0d0
  end do
  call update_molecular_orbitals(mo_coef(:,n_core_orb+1:occ_num),ao_num,val_occ_num,nto_hole)

! Rotate the virtual orbitals to the particle NTO
  call update_molecular_orbitals(mo_coef(:,occ_num+1:mo_num),ao_num,vir_num,nto_particle)

  call check_mos_orthonormality


! Save molecular orbitals to the EZFIO directory:
  call save_mos

! Resets the determinants and CI coefficients
  call save_ref_determinant


! Print NTO amplitudes
  character(len=*), parameter :: fmt = '(i4,3x,f15.8)'
  write(*,*)
  write(*,*) 'Natural transition orbital amplitudes:'
  do i=1,n_core_orb
    write(*,fmt) i, 0.0d0
  end do
  do i=n_core_orb+1,occ_num
    write(*,fmt) i, nto_amplitude(occ_num+1-i)**2 / 2.0d0
  end do
  do i=occ_num+1,occ_num+val_occ_num
    write(*,fmt) i, nto_amplitude(i-occ_num)**2 / 2.0d0
  end do
  trace_one_e_tdm_mo = sum( nto_amplitude(:)**2 / 2.0d0 )
  write(*,*)
  write(*,*) 'Trace of the 1-body transition density matrix = ', trace_one_e_tdm_mo

end
