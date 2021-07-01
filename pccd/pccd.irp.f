program pccd
  implicit none
  BEGIN_DOC
  ! This program performs coupled cluster with paired doubles (pCCD)
  END_DOC

  double precision,allocatable  :: t2(:,:)
  double precision,allocatable  :: z2(:,:)

  allocate(t2(r_val_occ_num,r_vir_num))
  allocate(z2(r_val_occ_num,r_vir_num))

  call run_pccd(t2,z2,r_val_occ_num,r_vir_num)

end program pccd
