 BEGIN_PROVIDER [ integer, occ_num ]
&BEGIN_PROVIDER [ integer, val_occ_num ]
&BEGIN_PROVIDER [ integer, vir_num ]
 implicit none
 BEGIN_DOC
 ! Assigns number of occupied, valence occupied, and virtual spatial orbitals
 END_DOC
 occ_num = elec_num / 2
 val_occ_num = occ_num - n_core_orb
 vir_num = mo_num - occ_num
END_PROVIDER
