BEGIN_PROVIDER [ double precision, mo_two_e_integrals, (mo_num,mo_num,mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
! Forms the two-electron integral in the spatial orbitals basis:
! mo_two_e_integrals(p,q,r,s) = v_{pq}^{rs} = <pq|rs> (Dirac notation)
  END_DOC
  integer :: p, q, r, s

  double precision, external :: get_two_e_integral

  provide mo_two_e_integrals_in_map

  do s=1,mo_num
    do r=1,mo_num
      do q=1,mo_num
        do p=1,mo_num
          mo_two_e_integrals(p,q,r,s)=get_two_e_integral(p,q,r,s,mo_two_e_integrals_in_map)
        end do
      end do
    end do
  end do

END_PROVIDER


BEGIN_PROVIDER [ integer, mo_num_a ]
  implicit none
  BEGIN_DOC
! Number of valence occupied spatial orbitals
  END_DOC
  mo_num_a = mo_num - n_frozen
END_PROVIDER


