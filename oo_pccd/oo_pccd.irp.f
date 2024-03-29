program oo_pccd
  implicit none
  BEGIN_DOC
  ! This program performs orbital optimized pair coupled cluster doubles (oo-pCCD), updating the orbitals in EZFIO at each iteration
  END_DOC

  double precision,allocatable  :: t2(:,:)
  double precision,allocatable  :: z2(:,:)
  integer :: iteration = 0
  logical :: is_converged = .false.

  provide mo_two_e_integrals_in_map

! Start of orbital optimization loop
  do while( .not.( is_converged .or. iteration==max_iter) )

    allocate(t2(r_val_occ_num,r_vir_num))
    allocate(z2(r_val_occ_num,r_vir_num))

!   Part 1: run pCCD
    call run_pCCD(t2,z2,r_val_occ_num,r_vir_num)

    if( n_frozen.gt.0 ) then

      double precision, allocatable :: t2_tmp(:,:)
      allocate( t2_tmp(r_val_occ_num,r_vir_num) )
      t2_tmp = t2
      deallocate( t2 )
      allocate( t2(r_occ_num,r_vir_num) )
      call add_zero_amplitudes_for_frozen_part(t2_tmp,t2,r_val_occ_num,r_occ_num,r_vir_num)
      deallocate( t2_tmp )

      double precision, allocatable :: z2_tmp(:,:)
      allocate( z2_tmp(r_val_occ_num,r_vir_num) )
      z2_tmp = z2
      deallocate( z2 )
      allocate( z2(r_occ_num,r_vir_num) )
      call add_zero_amplitudes_for_frozen_part(z2_tmp,z2,r_val_occ_num,r_occ_num,r_vir_num)
      deallocate( z2_tmp )

    end if

!   Part 2: compute 1-RDM and 2-RDM
    double precision,allocatable  :: r_one_e_dm_mo(:,:)
    double precision,allocatable  :: r_two_e_dm_mo(:,:,:,:)
    allocate( r_one_e_dm_mo(mo_num,mo_num) )
    allocate( r_two_e_dm_mo(mo_num,mo_num,mo_num,mo_num) )
    call compute_pccd_density_matrices(t2,z2,r_occ_num,r_vir_num,r_one_e_dm_mo,r_two_e_dm_mo,mo_num)

    deallocate( t2, z2 )

!   Part 3: orbital optimization
    call optimize_orbitals(r_one_e_dm_mo,r_two_e_dm_mo,mo_num,is_converged)

    deallocate( r_one_e_dm_mo, r_two_e_dm_mo)

    call clear_mo_map
    TOUCH mo_coef 

    iteration = iteration + 1

  end do
! End of orbital optimization loop

end program oo_pccd
