subroutine run_ci_mom(mo_coef_reference,swap)
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  double precision, intent(inout) :: mo_coef_reference(ao_num,mo_num)
  logical,          intent(out)   :: swap

  double precision, allocatable :: orbital_overlap(:,:)
  double precision, allocatable :: projected_overlap(:)
  double precision, allocatable :: mo_coef_old(:,:)
  integer,          allocatable :: orbs_occ(:)
  integer,          allocatable :: iorder(:)
  integer                       :: mo_num_no_core, n_occ_no_core
  integer                       :: n_doubly_occ_no_core
  integer                       :: i, j

  n_doubly_occ_no_core = n_doubly_occ - n_core_orb
  mo_num_no_core = mo_num - n_core_orb
  n_occ_no_core  = n_doubly_occ_no_core + n_singly_occ

  write(*,*) 'n_singly_occ, n_doubly_occ: ', n_singly_occ, n_doubly_occ
  write(*,*) 'n_occ_no_core, mo_num_no_core: ', n_occ_no_core, mo_num_no_core

  ! Compute the overlap between the current molecular orbitals and the initial molecular orbitals
  allocate( orbital_overlap(mo_num,mo_num) )
  call evaluate_orbital_overlap(mo_coef_reference,mo_coef,orbital_overlap)

  ! Compute the projection of the orbital overlap onto the occupied space of the reference MOs
  allocate( projected_overlap(mo_num_no_core) )
  projected_overlap = 0.0d0
  do j=1,mo_num_no_core
    do i=1,n_doubly_occ_no_core
      projected_overlap(j) = projected_overlap(j) + 4.0d0 * orbital_overlap(n_core_orb+i,n_core_orb+j)**2
    end do
  end do
  do j=1,mo_num_no_core
    do i=1,n_singly_occ
      projected_overlap(j) = projected_overlap(j) + orbital_overlap(n_doubly_occ+i,n_core_orb+j)**2
    end do
  end do
  do j=1,mo_num_no_core
    projected_overlap(j) = sqrt( projected_overlap(j) )
  end do

  deallocate( orbital_overlap )


  ! Set array of natural numbers
  allocate( iorder(mo_num_no_core) )
  call set_integer_list(iorder,mo_num_no_core)

  ! Sort the projected overlaps in decreasing order
  call dsort(projected_overlap,iorder,mo_num_no_core)

  ! Assign orbs_occ as the list of MOs (excluding the core) with the largest projected overlap
  allocate( orbs_occ(n_occ_no_core) )
  do i=1,n_occ_no_core
    orbs_occ(i) = iorder(mo_num_no_core+1-i) + n_core_orb
  end do

  ! Print the projected overlaps
  if(debug_mom) then
    write(*,*) 'projected overlap (occupied): '
    do i=1,n_occ_no_core
      write(*,*) orbs_occ(i), projected_overlap(mo_num_no_core+1-i)
    end do
    write(*,*) 'projected overlap (virtual): '
    do i=n_occ_no_core+1,mo_num_no_core
      write(*,*) iorder(mo_num_no_core+1-i) + n_core_orb, projected_overlap(mo_num_no_core+1-i)
    end do
    write(*,*)
  end if

  deallocate( projected_overlap )
  deallocate( iorder )


  ! Check if the MOs have to be swapped
  swap = .false.
  do i=n_core_orb+1,n_core_orb+n_occ_no_core
    if( .not. any(orbs_occ==i) ) then
      swap = .true.
      exit
    end if
  end do
  if(debug_mom) write(*,*) 'Swapped orbitals in MOM: ', swap


  ! Swap the MOs if required
  if( swap ) then
    allocate( mo_coef_old(ao_num,mo_num) )
    mo_coef_old = mo_coef
    integer :: count_occ, count_vir
    count_occ = 0
    count_vir = 0
    do i=n_core_orb+1,mo_num
      if( any(orbs_occ==i) ) then
        count_occ = count_occ + 1
        mo_coef(:,n_core_orb+count_occ) = mo_coef_old(:,i)
      else
        count_vir = count_vir + 1
        mo_coef(:,n_core_orb+n_occ_no_core+count_vir) = mo_coef_old(:,i)
      end if
    end do
    deallocate( mo_coef_old )
    TOUCH mo_coef
  end if

  deallocate( orbs_occ )


! Update the reference MOs for the original version of MOM. Do nothing for initial MOM (IMOM)
  if( mom_type == 'MOM' ) then
    mo_coef_reference = mo_coef
  end if

end subroutine
