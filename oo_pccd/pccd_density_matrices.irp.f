subroutine compute_pccd_density_matrices(t2,z2,nO,nV,r_one_e_dm_mo,r_two_e_dm_mo)
  implicit none
  BEGIN_DOC
! Computes the 1-RDM and 2-RDM of a pCCD wave function, in the spatial orbital basis
  END_DOC
  integer,          intent(in)  :: nO, nV
  double precision, intent(in)  :: t2(nO,nV), z2(nO,nV)
  double precision, intent(out) :: r_one_e_dm_mo(mo_num,mo_num)
  double precision, intent(out) :: r_two_e_dm_mo(mo_num,mo_num,mo_num,mo_num)
  
  call compute_pccd_one_e_dm_mo(r_one_e_dm_mo,t2,z2,nO,nV)

! call check_trace_one_e_dm(r_one_e_dm_mo,mo_num_a)

! write(*,*) 'One-body density matrix:'
! call write_ij_2d_array(r_one_e_dm_mo,mo_num,mo_num)

  call compute_pccd_two_e_dm_mo(r_two_e_dm_mo,t2,z2,nO,nV)

! call check_trace_two_e_dm(r_two_e_dm_mo,mo_num_a)

! write(*,*) 'Compute one_e_dm from the two_e_dm:'
! double precision :: one_from_two
! do q=1,mo_num
!   do p=1,mo_num
!     one_from_two = 0.d0
!     do r=1,mo_num
!       one_from_two = one_from_two + r_two_e_dm_mo(p,r,q,r)
!     end do
!     one_from_two = one_from_two / float( elec_num - 1 ) ! Check the sign of the denominator
!     write(*,*) p, q, r_one_e_dm_mo(p,q), one_from_two, r_one_e_dm_mo(p,q) - one_from_two
!   end do
! end do

! Prints total energy computed from one and two-body density matrices:
  double precision              :: energy
  call compute_core_energy(energy)
  call energy_from_density_matrices(r_one_e_dm_mo,r_two_e_dm_mo,mo_num,energy)
  write(*,*) 'Energy computed from reduced density matrices:', energy

end subroutine compute_pccd_density_matrices


! BEGIN_PROVIDER [ double precision, r_one_e_dm_mo, (mo_num,mo_num) ]
subroutine compute_pccd_one_e_dm_mo(r_one_e_dm_mo,t2,z2,nO,nV)
   implicit none
   BEGIN_DOC
   ! Computes the one-body reduced density matrix (1-RDM) of a pCCD wave function, in the spatial orbital basis
   ! r_one_e_dm_mo(p,q) = gamma_pq = sum_@ <0| (1+Z) exp(-T) c*_p@ c_q@ exp(T) |0>
   END_DOC
   integer,          intent(in)  :: nO, nV
   double precision, intent(in)  :: t2(nO,nV), z2(nO,nV)
   double precision, intent(out) :: r_one_e_dm_mo(mo_num,mo_num)

   double precision, allocatable :: xOO(:,:), xVV(:,:)
   integer                       :: i, j, a, b

   allocate( xOO(nO,nO), xVV(nV,nV) )
   xOO = 0.0d0
   xVV = 0.0d0

!  Intermediates xOO and xVV
!  These are the general expressions:
!  do j=1,nO
!    do i=1,nO
!      do a=1,nV
!        xOO(i,j) = xOO(i,j) + t2(i,a) * z2(j,a)
!      end do
!    end do
!  end do
!  do b=1,nV
!    do a=1,nV
!      do i=1,nO
!        xVV(a,b) = xVV(a,b) + t2(i,b) * z2(i,a)
!      end do
!    end do
!  end do
!  But only the diagonal terms are effectively needed for pCCD:
     do i=1,nO
       do a=1,nV
         xOO(i,i) = xOO(i,i) + t2(i,a) * z2(i,a)
       end do
     end do
     do a=1,nV
       do i=1,nO
         xVV(a,a) = xVV(a,a) + t2(i,a) * z2(i,a)
       end do
     end do

!  All off-diagonal elements are zero
   r_one_e_dm_mo = 0.0d0

!  OO diagonal
   do i=1,nO
     r_one_e_dm_mo(i,i) = 2.0d0 * ( 1.0d0 - xOO(i,i) )
   end do

!  VV diagonal
   do a=1,nV
     r_one_e_dm_mo(nO+a,nO+a) = 2.0d0 * xVV(a,a) 
   end do

end subroutine compute_pccd_one_e_dm_mo
!END_PROVIDER


!BEGIN_PROVIDER [ double precision, r_two_e_dm_mo, (mo_num,mo_num,mo_num,mo_num) ]
subroutine compute_pccd_two_e_dm_mo(r_two_e_dm_mo,t2,z2,nO,nV)
   implicit none
   BEGIN_DOC
   ! Computes the two-body reduced density matrix (2-RDM) of a pCCD wave function, in the spatial orbital basis
   ! r_two_e_dm_mo(p,q,r,s) = Gamma_pqrs = sum_@ sum_@' <0| (1+Z) exp(-T) c*_p@ c*_q@' c_s@' c_r@ exp(T) |0>
   END_DOC
   integer,          intent(in)  :: nO, nV
   double precision, intent(in)  :: t2(nO,nV), z2(nO,nV)
   double precision, intent(out) :: r_two_e_dm_mo(mo_num,mo_num,mo_num,mo_num)

   double precision, allocatable :: xOO(:,:), xVV(:,:), xOV(:,:)
   double precision              :: tmp
   integer                       :: i, j, a, b, p, q

   allocate( xOO(nO,nO), xVV(nV,nV), xOV(nO,nV) )
   xOO = 0.0d0
   xVV = 0.0d0
   xOV = 0.0d0

!  Intermediates xOO and xVV
   do j=1,nO
     do i=1,nO
       do a=1,nV
         xOO(i,j) = xOO(i,j) + t2(i,a) * z2(j,a)
       end do
     end do
   end do
   do b=1,nV
     do a=1,nV
       do i=1,nO
         xVV(a,b) = xVV(a,b) + t2(i,b) * z2(i,a)
       end do
     end do
   end do
!  Intermediates xOV
   do i=1,nO
     do a=1,nV
       do b=1,nV
         do j=1,nO
           xOV(i,a) = xOV(i,a) + t2(i,b) * t2(j,a) * z2(j,b)
         end do
       end do
     end do
   end do

   r_two_e_dm_mo = 0.0d0

! Gamma_iijj
   do i=1,nO
     do j=1,nO
       r_two_e_dm_mo(i,i,j,j) = 2.0d0 * xOO(i,j)
     end do
   end do
!  I will set the j=i case later, just after the ijij case
!  do i=1,nO
!    r_two_e_dm_mo(i,i,i,i) = 2.0d0 * ( 1.0d0 - xOO(i,i) )
!  end do
!  do i=1,nO
!    r_two_e_dm_mo(i,i,i,i) = r_two_e_dm_mo(i,i,i,i) + 2.0d0 * (1.0d0 - 2.0d0 * xOO(i,i) )
!  end do
  
! Gamma_iiaa
   do a=1,nV
     do i=1,nO
       r_two_e_dm_mo(i,i,nO+a,nO+a) = 2.0d0 * ( t2(i,a) + xOV(i,a) &
                    - 2.0d0 * t2(i,a) * ( xVV(a,a) + xOO(i,i) - t2(i,a) * z2(i,a) ) )
     end do
   end do

! Gamma_aaii
   do i=1,nO
     do a=1,nV
       r_two_e_dm_mo(nO+a,nO+a,i,i) = 2.0d0 * z2(i,a)
     end do
   end do

! Gamma_aabb
   do b=1,nV
     do a=1,nV
       r_two_e_dm_mo(nO+a,nO+a,nO+b,nO+b) = 2.0d0 * xVV(a,b)
     end do
   end do

! Gamma_ijij
   do j=1,nO
     do i=1,nO
       r_two_e_dm_mo(i,j,i,j) = 4.0d0 * ( 1.0d0 - xOO(i,i) - xOO(j,j) )
! Gamma_ijji 
! Do not need to worry about the j=i case, which will be overwritten below
       r_two_e_dm_mo(i,j,j,i) = - r_two_e_dm_mo(i,j,i,j) / 2.d0
     end do
   end do
!  do i=1,nO
!    r_two_e_dm_mo(i,i,i,i) = r_two_e_dm_mo(i,i,i,i) + 2.0d0 * ( 3.0d0 * xOO(i,i) - 1.0d0 )
!  end do

!  This is from the first expression. I use it here to avoid repetition:
   do i=1,nO
     r_two_e_dm_mo(i,i,i,i) = 2.0d0 * ( 1.0d0 - xOO(i,i) )
   end do

   do a=1,nV
     do i=1,nO
       tmp =  4.0d0 * ( xVV(a,a) - t2(i,a) * z2(i,a) )
! Gamma_iaia
       r_two_e_dm_mo(i,nO+a,i,nO+a) = tmp
! Gamma_iaai
       r_two_e_dm_mo(i,nO+a,nO+a,i) = - tmp / 2.0d0
! Gamma_aiai 
       r_two_e_dm_mo(nO+a,i,nO+a,i) = tmp
! Gamma_aiia
       r_two_e_dm_mo(nO+a,i,i,nO+a) = - tmp / 2.0d0
     end do
   end do

! Gamma_aiai
!  do i=1,nO
!    do a=1,nV
!      r_two_e_dm_mo(nO+a,i,nO+a,i) = r_two_e_dm_mo(i,nO+a,i,nO+a)
!    end do
!  end do

! Gamma_abab
!  do b=1,nV
     do a=1,nV
       r_two_e_dm_mo(nO+a,nO+a,nO+a,nO+a) = 2.0d0 * xVV(a,a)
     end do
!  end do

! Gamma_pqqp = - Gamma_pqpq / 2   (p.ne.q)
!  do q=1,mo_num
!    do p=1,mo_num
!      if( p.eq.q ) cycle
!      r_two_e_dm_mo(p,q,q,p) = - r_two_e_dm_mo(p,q,p,q) / 2.0d0
!    end do
!  end do


!open(unit=10,file='test_2rdm.dat')
!integer :: k, l
!  do l= 1, mo_num
!   do k = 1, mo_num
!     do j = 1, mo_num
!       do i = 1, mo_num
!!                write(10,*)  i, j, k, l, r_two_e_dm_mo(i,j,k,l)
!                 write(10,'(4(i2,x),e23.15)')  i, j, k, l, r_two_e_dm_mo(i,j,k,l)
!       enddo
!     enddo
!   enddo
! enddo
!close(10)

end subroutine compute_pccd_two_e_dm_mo
!END_PROVIDER


subroutine compute_core_energy(energy)
  implicit none
  BEGIN_DOC
  ! Assigns the nuclear repulsion energy
  END_DOC
  double precision, intent(out) :: energy
  energy = nuclear_repulsion
end subroutine compute_core_energy


subroutine energy_from_density_matrices(r_one_e_dm_mo,r_two_e_dm_mo,n,energy)
  implicit none
  BEGIN_DOC
  ! Computes the electronic energy from the 1-RDM and 2-RDM, written in the spatial orbital basis
  END_DOC
  integer         , intent(in)    :: n
  double precision, intent(in)    :: r_one_e_dm_mo(n,n)
  double precision, intent(in)    :: r_two_e_dm_mo(n,n,n,n)
  double precision, intent(inout) :: energy
  integer          :: p, q, r, s
  do q=1,n
    do p=1,n
      energy = energy + mo_one_e_integrals(p,q) * r_one_e_dm_mo(p,q)
    end do
  end do
  do s=1,n
    do r=1,n
      do q=1,n
        do p=1,n
          energy = energy + mo_two_e_integrals(p,q,r,s) * r_two_e_dm_mo(r,s,p,q) / 2.0d0
        end do
      end do
    end do
  end do
end subroutine energy_from_density_matrices


subroutine check_trace_one_e_dm(r_one_e_dm_mo,n)
  implicit none
  BEGIN_DOC
  ! Checkes whether the trace of the 1-RDM equals the number of electrons
  END_DOC
  integer,          intent(in) :: n
  double precision, intent(in) :: r_one_e_dm_mo(n,n)
  double precision             :: trace
  double precision, parameter  :: thresh = 1.0d-10
  integer                      :: p
  trace = 0.0d0
  do p=1,n
    trace = trace + r_one_e_dm_mo(p,p)
  end do
  if( abs( trace - float(elec_num) ) .gt. thresh ) then
    write(*,*) 'WARNING: problem with one-electron density matrix:'
    write(*,*) 'Trace of 1-RDM:      ', trace
    write(*,*) 'Number of electrons: ', elec_num
  end if
end subroutine check_trace_one_e_dm


subroutine check_trace_two_e_dm(r_two_e_dm_mo,n)
  implicit none
  BEGIN_DOC
  ! Checkes whether the trace of the 2-RDM equals N*(N-1), N being the number of electrons
  END_DOC
  integer,          intent(in) :: n
  double precision, intent(in) :: r_two_e_dm_mo(n,n,n,n)
  double precision             :: trace
  double precision, parameter  :: thresh = 1.0d-10
  integer                      :: p, q
  trace = 0.0d0
  do q=1,n
    do p=1,n
      trace = trace + r_two_e_dm_mo(p,q,p,q)
    end do
  end do
  if( abs( trace - float(elec_num*(elec_num-1)) ) .gt. thresh ) then
    write(*,*) 'WARNING: problem with two-electron density matrix:'
    write(*,*) 'Trace of 2-RDM / ( Number of electrons - 1 ): ', trace
    write(*,*) 'Number of electrons:                          ', elec_num
  end if
end subroutine check_trace_two_e_dm

