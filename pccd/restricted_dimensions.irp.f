 BEGIN_PROVIDER [ integer, r_occ_num ]
&BEGIN_PROVIDER [ integer, r_vir_num ]
&BEGIN_PROVIDER [ integer, r_val_occ_num ]
 implicit none
 BEGIN_DOC
 ! Assigns number of occupied, valence occupied, and virtual spatial orbitals
 END_DOC
 r_occ_num = elec_num / 2
 r_val_occ_num = r_occ_num - n_frozen
 r_vir_num = mo_num - r_occ_num
END_PROVIDER
