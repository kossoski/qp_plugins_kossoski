subroutine compute_r_orbrot_g(r_orbrot_g,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  implicit none
  BEGIN_DOC
  ! Computes the gradients of the orbital rotation matrix, in the spatial orbital basis, given the 1-RDM, 2-RDM, and the 1- and 2-electron integrals
  ! Implementation of eq. 25 from JCP 141 244104 (2014)
  END_DOC

  integer         , intent(in)  :: nT
  double precision, intent(out) :: r_orbrot_g(nT,nT)
  double precision, intent(in)  :: r_one_e_dm_mo(nT,nT)
  double precision, intent(in)  :: r_two_e_dm_mo(nT,nT,nT,nT)

  integer          :: p, q, r, s

! Auxiliary array to speed up the contractions (at the cost of some memory)
  double precision, allocatable  :: r_one_e_dm_mo_t(:,:)
  allocate( r_one_e_dm_mo_t(nT,nT) )
  r_one_e_dm_mo_t = transpose( r_one_e_dm_mo )

  r_orbrot_g = 0.d0

! Commented lines are from the initial code, when permutation symmetries of one and two electron integrals were not explored,


  do q=n_frozen+1,nT
    do p=n_frozen+1,nT
      if(p==q) cycle
      r_orbrot_g(p,q) = r_orbrot_g(p,q) + sum( mo_one_e_integrals(:,p) * r_one_e_dm_mo(:,q) &
                                             - mo_one_e_integrals(:,q) * r_one_e_dm_mo_t(:,p) )
!                                            - mo_one_e_integrals(:,q) * r_one_e_dm_mo(p,:) )
      do s=1,nT
        do r=1,nT
          r_orbrot_g(p,q) = r_orbrot_g(p,q) + sum( mo_two_e_integrals(:,r,s,p) * r_two_e_dm_mo(r,s,q,:) &
                                                 - mo_two_e_integrals(:,r,s,q) * r_two_e_dm_mo(p,:,r,s) )
        end do
      end do
    end do
  end do

  deallocate( r_one_e_dm_mo_t )


! Permutation:
  do q=1,nT
    do p=q+1,nT
      r_orbrot_g(p,q) = r_orbrot_g(p,q) - r_orbrot_g(q,p)
    end do
  end do


! Eliminate very small elements
  double precision :: threshold 
  threshold = 1.d-16
  call zero_small_2index(r_orbrot_g,nT,nT,threshold)

end subroutine compute_r_orbrot_g


subroutine compute_r_orbrot_h(r_orbrot_h,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  implicit none
  BEGIN_DOC
  ! Computes the Hessian matrix elements of the orbital rotation matrix, in the spatial orbital basis, given the 1-RDM, 2-RDM, and the 1- and 2-electron integrals
  ! Implementation of eq. A6 from JCP 141 244104 (2014)
  END_DOC

  integer         , intent(in)  :: nT
  double precision, intent(out) :: r_orbrot_h(nT,nT,nT,nT)
  double precision, intent(in)  :: r_one_e_dm_mo(nT,nT)
  double precision, intent(in)  :: r_two_e_dm_mo(nT,nT,nT,nT)

  integer :: p, q, r, s, t, v

! Some auxiliary arrays to speed up the contractions (at the cost of some memory)

  double precision, allocatable  :: r_one_e_dm_mo_t(:,:)
  allocate( r_one_e_dm_mo_t(nT,nT) )
  r_one_e_dm_mo_t = transpose( r_one_e_dm_mo )

  double precision, allocatable  :: r_two_e_dm_mo_t1(:,:,:,:)
  allocate( r_two_e_dm_mo_t1(nT,nT,nT,nT) )
  do s=1,nT
    do r=1,nT
      do q=1,nT
        do p=1,nT
          r_two_e_dm_mo_t1(p,q,r,s) = r_two_e_dm_mo(s,q,r,p)
        end do
      end do
    end do
  end do

  double precision, allocatable  :: r_two_e_dm_mo_t3(:,:,:,:)
  allocate( r_two_e_dm_mo_t3(nT,nT,nT,nT) )
  do s=1,nT
    do r=1,nT
      do q=1,nT
        do p=1,nT
          r_two_e_dm_mo_t3(p,q,r,s) = r_two_e_dm_mo(r,s,p,q)
        end do
      end do
    end do
  end do


  r_orbrot_h = 0.d0

! Commented lines are from the initial code, when permutation symmetries of one and two electron integrals were not explored,


  do s=n_frozen+1,nT
    do r=n_frozen+1,nT
      if(r==s) cycle
      do q=n_frozen+1,nT
        do p=n_frozen+1,nT
          if(p==q) cycle

! Third line, first term:
          do v=1,nT
            r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) + sum( mo_two_e_integrals(:,v,p,r) * r_two_e_dm_mo(:,v,q,s) &
                                                           + mo_two_e_integrals(:,v,q,s) * r_two_e_dm_mo_t3(:,v,p,r) )
!                                                          + mo_two_e_integrals(:,v,q,s) * r_two_e_dm_mo(p,r,:,v) ) )
          end do

! Third line, second term:
          do t=1,nT
            r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - sum( mo_two_e_integrals(:,p,t,s) * r_two_e_dm_mo_t1(:,t,q,r) &
                                                           + mo_two_e_integrals(:,t,s,p) * r_two_e_dm_mo_t1(:,r,q,t) &
                                                           + mo_two_e_integrals(:,q,t,r) * r_two_e_dm_mo_t1(t,:,s,p) &
                                                           + mo_two_e_integrals(:,t,r,q) * r_two_e_dm_mo(p,:,t,s) )
!           r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - sum( mo_two_e_integrals(:,p,t,s) * r_two_e_dm_mo(r,t,q,:) &
!                                                          + mo_two_e_integrals(:,t,s,p) * r_two_e_dm_mo(t,r,q,:) &
!                                                          + mo_two_e_integrals(:,q,t,r) * r_two_e_dm_mo(p,:,s,t) &
!                                                          + mo_two_e_integrals(:,t,r,q) * r_two_e_dm_mo(p,:,t,s) )
          end do

        end do
! First line, third term:
        r_orbrot_h(:,q,r,s) = r_orbrot_h(:,q,r,s) - ( mo_one_e_integrals(:,s) * r_one_e_dm_mo_t(q,r) + mo_one_e_integrals(q,r) * r_one_e_dm_mo(:,s) )  
!       r_orbrot_h(:,q,r,s) = r_orbrot_h(:,q,r,s) - ( mo_one_e_integrals(:,s) * r_one_e_dm_mo(r,q) + mo_one_e_integrals(q,r) * r_one_e_dm_mo(:,s) )  
      end do
    end do
  end do

  deallocate( r_two_e_dm_mo_t1 )
  deallocate( r_two_e_dm_mo_t3 )


  double precision, allocatable  :: r_two_e_dm_mo_t2(:,:,:,:)
  allocate( r_two_e_dm_mo_t2(nT,nT,nT,nT) )
  do s=1,nT
    do r=1,nT
      do q=1,nT
        do p=1,nT
          r_two_e_dm_mo_t2(p,q,r,s) = r_two_e_dm_mo(s,r,p,q)
        end do
      end do
    end do
  end do

! First line, first term:
  do s=n_frozen+1,nT
      do q=n_frozen+1,nT
        if(q==s) cycle
        do p=n_frozen+1,nT
          if(p==q) cycle
          r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) + sum ( mo_one_e_integrals(:,p) * r_one_e_dm_mo(:,s) &
                                                          + mo_one_e_integrals(:,s) * r_one_e_dm_mo_t(:,p) ) / 2.0d0
!                                                         + mo_one_e_integrals(:,s) * r_one_e_dm_mo(p,:) ) / 2.0d0
! Second line, first term:
          do t=1,nT
            do v=1,nT
              r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) + sum ( mo_two_e_integrals(:,v,p,t) * r_two_e_dm_mo(:,v,s,t) &
                                                              + mo_two_e_integrals(:,v,s,t) * r_two_e_dm_mo_t2(:,v,t,p) ) / 2.0d0
!                                                             + mo_two_e_integrals(:,v,s,t) * r_two_e_dm_mo(p,t,:,v) ) / 2.0d0
            end do
          end do
        end do
      end do
  end do


    do r=n_frozen+1,nT
      do q=n_frozen+1,nT
        do p=n_frozen+1,nT
          if(p==q) cycle
          if(p==r) cycle
! First line, second term:
          r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) + sum ( mo_one_e_integrals(:,r) * r_one_e_dm_mo(:,q) &
                                                          + mo_one_e_integrals(:,q) * r_one_e_dm_mo_t(:,r) ) / 2.0d0
!                                                         + mo_one_e_integrals(:,q) * r_one_e_dm_mo(r,:) ) / 2.0d0
! Second line, second term:
          do t=1,nT
            do v=1,nT
              r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) + sum ( mo_two_e_integrals(:,v,q,t) * r_two_e_dm_mo_t2(:,v,t,r) &
                                                              + mo_two_e_integrals(:,v,r,t) * r_two_e_dm_mo(:,v,q,t) ) / 2.0d0
!             r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) + sum ( mo_two_e_integrals(:,v,q,t) * r_two_e_dm_mo(r,t,:,v) &
            end do
          end do
        end do
      end do
    end do

  deallocate( r_one_e_dm_mo_t )
  deallocate( r_two_e_dm_mo_t2 )


! Permutation:
  do s=n_frozen+1,nT
    do r=s+1,nT
      do q=n_frozen+1,nT
        do p=q+1,nT
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - r_orbrot_h(p,q,s,r) - r_orbrot_h(q,p,r,s) + r_orbrot_h(q,p,s,r)
        end do
      end do
    end do
  end do


! Frozen part
  if( n_frozen.gt.0 ) then
    do q=1,nT
      do p=1,n_frozen
        if(p==q) cycle
        r_orbrot_h(p,q,p,q) = 1.0d0
        r_orbrot_h(q,p,q,p) = 1.0d0
      end do
    end do
  end if


! Eliminate very small elements
  double precision :: threshold 
  threshold = 1.d-16
  call zero_small_4index(r_orbrot_h,nT,nT,nT,nT,threshold)

end subroutine compute_r_orbrot_h


subroutine compute_r_orbrot_h_diagonal(r_orbrot_h_diag,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  implicit none
  BEGIN_DOC
  ! Computes the diagonal Hessian matrix elements of the orbital rotation matrix, in the spatial orbital basis, given the 1-RDM, 2-RDM, and the 1- and 2-electron integrals
  ! Implementation of eq. A6 from JCP 141 244104 (2014)
  END_DOC

  integer         , intent(in)  :: nT
  double precision, intent(out) :: r_orbrot_h_diag(nT,nT)
  double precision, intent(in)  :: r_one_e_dm_mo(nT,nT)
  double precision, intent(in)  :: r_two_e_dm_mo(nT,nT,nT,nT)

  integer :: p, q, r, s, t, v

! Perhaps it is not worthy introducing these arrays in this case, when only the diagonal hessian is computed 
  double precision, allocatable  :: r_one_e_dm_mo_t(:,:)
  allocate( r_one_e_dm_mo_t(nT,nT) )
  r_one_e_dm_mo_t = transpose( r_one_e_dm_mo )

  double precision, allocatable  :: r_two_e_dm_mo_t1(:,:,:,:)
  allocate( r_two_e_dm_mo_t1(nT,nT,nT,nT) )
  do s=1,nT
    do r=1,nT
      do q=1,nT
        do p=1,nT
          r_two_e_dm_mo_t1(p,q,r,s) = r_two_e_dm_mo(s,q,r,p)
        end do
      end do
    end do
  end do

  double precision, allocatable  :: r_two_e_dm_mo_t3(:,:,:,:)
  allocate( r_two_e_dm_mo_t3(nT,nT,nT,nT) )
  do s=1,nT
    do r=1,nT
      do q=1,nT
        do p=1,nT
          r_two_e_dm_mo_t3(p,q,r,s) = r_two_e_dm_mo(r,s,p,q)
        end do
      end do
    end do
  end do


  r_orbrot_h_diag = 0.d0


  do q=n_frozen+1,nT
    do p=n_frozen+1,nT

! Third line, first term:
      do v=1,nT
        r_orbrot_h_diag(p,q) = r_orbrot_h_diag(p,q) + sum( mo_two_e_integrals(:,v,p,p) * r_two_e_dm_mo(:,v,q,q) &
                                                         + mo_two_e_integrals(:,v,q,q) * r_two_e_dm_mo_t3(:,v,p,p) )
      end do
! Third line, second term:
      do t=1,nT
        r_orbrot_h_diag(p,q) = r_orbrot_h_diag(p,q) - sum( mo_two_e_integrals(:,p,t,q) * ( r_two_e_dm_mo_t1(:,t,q,p) + r_two_e_dm_mo_t1(t,:,q,p) ) &
                                                         + mo_two_e_integrals(:,t,q,p) * r_two_e_dm_mo_t1(:,p,q,t) &
                                                         + mo_two_e_integrals(:,t,p,q) * r_two_e_dm_mo(p,:,t,q) )
      end do
! First line, third term:
      r_orbrot_h_diag(p,q) = r_orbrot_h_diag(p,q) - 2.0d0 * mo_one_e_integrals(p,q) * r_one_e_dm_mo(p,q) 

    end do
  end do

  deallocate( r_two_e_dm_mo_t1 )
  deallocate( r_two_e_dm_mo_t3 )


  double precision, allocatable  :: r_two_e_dm_mo_t2(:,:,:,:)
  allocate( r_two_e_dm_mo_t2(nT,nT,nT,nT) )
  do s=1,nT
    do r=1,nT
      do q=1,nT
        do p=1,nT
          r_two_e_dm_mo_t2(p,q,r,s) = r_two_e_dm_mo(s,r,p,q)
        end do
      end do
    end do
  end do

! Second line, first and second terms:
  do p=n_frozen+1,nT
    do t=1,nT
      do v=1,nT
        r_orbrot_h_diag(p,p) = r_orbrot_h_diag(p,p) + sum( mo_two_e_integrals(:,v,p,t) * ( r_two_e_dm_mo(:,v,p,t) + r_two_e_dm_mo_t2(:,v,t,p) ) )
      end do
    end do
  end do
  deallocate( r_two_e_dm_mo_t2 )

! First line, first and second terms:
  do p=n_frozen+1,nT
    r_orbrot_h_diag(p,p) = r_orbrot_h_diag(p,p) + sum( mo_one_e_integrals(:,p) * ( r_one_e_dm_mo(:,p) + r_one_e_dm_mo_t(:,p) ) )
  end do
  deallocate( r_one_e_dm_mo_t )

! Permutation:
  do q=n_frozen+1,nT
    do p=q+1,nT
      r_orbrot_h_diag(p,q) = r_orbrot_h_diag(p,q) + r_orbrot_h_diag(q,p) 
    end do
  end do

! Frozen part
  if( n_frozen.gt.0 ) then
    do q=1,nT
      do p=1,n_frozen
        r_orbrot_h_diag(p,q) = 1.0d0
        r_orbrot_h_diag(q,p) = 1.0d0
      end do
    end do
  end if

! Eliminate very small elements
  double precision :: threshold
  threshold = 1.d-16
  call zero_small_2index(r_orbrot_h_diag,nT,nT,threshold)

end subroutine compute_r_orbrot_h_diagonal
