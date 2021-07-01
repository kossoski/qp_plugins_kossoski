subroutine compute_r_orbrot_g_pccd(r_orbrot_g,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  implicit none
  BEGIN_DOC
  ! Computes the gradients of the orbital rotation matrix, in the spatial orbital basis, given the 1-RDM, 2-RDM, and the 1- and 2-electron integrals
  ! Implementation of eq. 25 from JCP 141 244104 (2014)
  ! This subroutine is specific to a pCCD wave function
  END_DOC

  integer         , intent(in)  :: nT
  double precision, intent(out) :: r_orbrot_g(nT,nT)
  double precision, intent(in)  :: r_one_e_dm_mo(nT,nT)
  double precision, intent(in)  :: r_two_e_dm_mo(nT,nT,nT,nT)

  integer          :: p, q, r, s


  r_orbrot_g = 0.d0

! Commented lines are from the initial code, when permutation symmetries of one and two electron integrals were not explored,


  do q=n_frozen+1,nT
    do p=n_frozen+1,nT
      if(p==q) cycle
      r_orbrot_g(p,q) = r_orbrot_g(p,q) + mo_one_e_integrals(p,q) * ( r_one_e_dm_mo(q,q) - r_one_e_dm_mo(p,p) )

      r_orbrot_g(p,q) = r_orbrot_g(p,q) + 2.0d0 * ( mo_two_e_integrals(p,p,p,q) * r_two_e_dm_mo(p,p,p,p) &
                                                  - mo_two_e_integrals(p,q,q,q) * r_two_e_dm_mo(q,q,q,q) )

      do s=1,nT
        r_orbrot_g(p,q) = r_orbrot_g(p,q) + mo_two_e_integrals(s,s,p,q) * ( r_two_e_dm_mo(s,s,q,q) &
                                                                          + r_two_e_dm_mo(s,q,q,s) &
                                                                          - r_two_e_dm_mo(p,p,s,s) &
                                                                          - r_two_e_dm_mo(p,s,s,p) ) &
                                          + mo_two_e_integrals(s,p,s,q) * ( r_two_e_dm_mo(q,s,q,s) &
                                                                          - r_two_e_dm_mo(p,s,p,s) )
      end do

    end do
  end do


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

end subroutine compute_r_orbrot_g_pccd


subroutine compute_r_orbrot_h_pccd(r_orbrot_h,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  implicit none
  BEGIN_DOC
  ! Computes the Hessian matrix elements of the orbital rotation matrix, in the spatial orbital basis, given the 1-RDM, 2-RDM, and the 1- and 2-electron integrals
  ! Implementation of eq. A6 from JCP 141 244104 (2014)
  ! This subroutine is specific for the case of pCCD
  END_DOC

  integer         , intent(in)  :: nT
  double precision, intent(out) :: r_orbrot_h(nT,nT,nT,nT)
  double precision, intent(in)  :: r_one_e_dm_mo(nT,nT)
  double precision, intent(in)  :: r_two_e_dm_mo(nT,nT,nT,nT)

  integer :: p, q, r, s, t, v

! Some auxiliary arrays to speed up the contractions (at the cost of some memory)

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


  do s=n_frozen+1,nT
    do r=n_frozen+1,nT
      if(r==s) cycle
      do q=n_frozen+1,nT

        r_orbrot_h(q,r,r,s) = r_orbrot_h(q,r,r,s) - mo_one_e_integrals(q,s) * r_one_e_dm_mo(r,r) 
        r_orbrot_h(s,q,r,s) = r_orbrot_h(s,q,r,s) - mo_one_e_integrals(q,r) * r_one_e_dm_mo(s,s) 

        do t=1,nT
          r_orbrot_h(q,r,r,s) = r_orbrot_h(q,r,r,s) - mo_two_e_integrals(t,t,q,s) * r_two_e_dm_mo(t,r,r,t) 
        end do

        do t=1,nT
          r_orbrot_h(q,r,r,s) = r_orbrot_h(q,r,r,s) - mo_two_e_integrals(t,q,t,s) * r_two_e_dm_mo_t1(t,t,r,r)
        end do

        do t=1,nT
          r_orbrot_h(q,s,r,s) = r_orbrot_h(q,s,r,s) + mo_two_e_integrals(t,t,q,r) * r_two_e_dm_mo(t,t,s,s) 
        end do

        do t=1,nT
          r_orbrot_h(r,q,r,s) = r_orbrot_h(r,q,r,s) + mo_two_e_integrals(t,t,q,s) * r_two_e_dm_mo_t3(t,t,r,r)
        end do

        do t=1,nT
          r_orbrot_h(s,q,r,s) = r_orbrot_h(s,q,r,s) - mo_two_e_integrals(t,t,q,r) * r_two_e_dm_mo(s,t,t,s) 
        end do

        do t=1,nT
          r_orbrot_h(s,q,r,s) = r_orbrot_h(s,q,r,s) - mo_two_e_integrals(t,q,t,r) * r_two_e_dm_mo_t1(t,t,s,s)
        end do


        do p=n_frozen+1,nT
          if(p==q) cycle
          if(p==r) cycle
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) + mo_two_e_integrals(p,r,q,s) * r_two_e_dm_mo(p,r,p,r) &
                                                    + mo_two_e_integrals(p,q,s,r) * r_two_e_dm_mo(p,r,r,p) 
        end do

        do p=n_frozen+1,nT
          if(p==q) cycle
          if(q==r) cycle
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - mo_two_e_integrals(p,q,s,r) * ( r_two_e_dm_mo_t1(q,r,q,r) &
                                                                                    + r_two_e_dm_mo_t1(r,q,q,r) )
        end do

        do p=n_frozen+1,nT
          if(p==q) cycle
          if(q==r) cycle
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - mo_two_e_integrals(p,q,r,s) * r_two_e_dm_mo_t1(q,r,q,r) & 
                                                    - mo_two_e_integrals(p,r,q,s) * r_two_e_dm_mo(q,r,q,r) 
        end do

        do p=n_frozen+1,nT
          if(p==q) cycle
          if(p==s) cycle
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - mo_two_e_integrals(p,q,s,r) * ( r_two_e_dm_mo(p,s,s,p) &
                                                                                    + r_two_e_dm_mo(p,p,s,s) )
        end do

        do p=n_frozen+1,nT
          if(p==q) cycle
          if(p==s) cycle
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) - mo_two_e_integrals(p,r,q,s) * r_two_e_dm_mo(p,s,p,s) &
                                                    - mo_two_e_integrals(p,q,r,s) * r_two_e_dm_mo(p,p,s,s)
        end do

        do p=n_frozen+1,nT
          if(p==q) cycle
          if(q==s) cycle
          r_orbrot_h(p,q,r,s) = r_orbrot_h(p,q,r,s) + mo_two_e_integrals(p,r,q,s) * r_two_e_dm_mo(q,s,q,s) &
                                                    + mo_two_e_integrals(p,q,s,r) * r_two_e_dm_mo(s,q,q,s) 
        end do

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

  do s=n_frozen+1,nT
      do q=n_frozen+1,nT
        if(q==s) cycle
        do p=n_frozen+1,nT
          if(p==q) cycle

          r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) + mo_one_e_integrals(p,s) * ( r_one_e_dm_mo(s,s) + r_one_e_dm_mo(p,p) ) / 2.0d0

          r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) - mo_two_e_integrals(p,s,s,s) * r_two_e_dm_mo(s,s,s,s) 

          r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) - mo_two_e_integrals(p,p,p,s) * r_two_e_dm_mo(p,p,p,p) 

          do t=1,nT
            r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) + mo_two_e_integrals(t,t,p,s) * ( r_two_e_dm_mo(t,t,s,s) &
                                                      + r_two_e_dm_mo(t,s,s,t) &
                                                      + r_two_e_dm_mo_t2(t,t,p,p) &
                                                      + r_two_e_dm_mo_t2(t,p,t,p) ) / 2.0d0
          end do

          do t=1,nT
            r_orbrot_h(p,q,q,s) = r_orbrot_h(p,q,q,s) + mo_two_e_integrals(t,p,t,s) * ( r_two_e_dm_mo(s,t,s,t) &
                                                      + r_two_e_dm_mo_t2(p,t,t,p) ) / 2.0d0
          end do

        end do
      end do
  end do


    do r=n_frozen+1,nT
      do q=n_frozen+1,nT
        do p=n_frozen+1,nT
          if(p==q) cycle
          if(p==r) cycle

          r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) + mo_one_e_integrals(q,r) * ( r_one_e_dm_mo(q,q) + r_one_e_dm_mo(r,r) ) / 2.0d0

          r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) - mo_two_e_integrals(q,r,r,r) * r_two_e_dm_mo_t2(r,r,r,r)

          r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) - mo_two_e_integrals(q,q,q,r) * r_two_e_dm_mo(q,q,q,q)

          do t=1,nT
            r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) + mo_two_e_integrals(t,t,q,r) * ( r_two_e_dm_mo_t2(t,t,r,r) &
                                                                                      + r_two_e_dm_mo_t2(t,r,t,r) &
                                                                                      + r_two_e_dm_mo(t,t,q,q)    &
                                                                                      + r_two_e_dm_mo(t,q,q,t)    ) / 2.0d0
          end do

          do t=1,nT
            r_orbrot_h(p,q,r,p) = r_orbrot_h(p,q,r,p) + mo_two_e_integrals(t,q,t,r) * ( r_two_e_dm_mo_t2(r,t,t,r) &
                                                                                      + r_two_e_dm_mo(q,t,q,t)    ) / 2.0d0
          end do

        end do
      end do
    end do

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

end subroutine compute_r_orbrot_h_pccd


subroutine compute_r_orbrot_h_diagonal_pccd(r_orbrot_h_diag,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  implicit none
  BEGIN_DOC
  ! Computes the diagonal Hessian matrix elements of the orbital rotation matrix, in the spatial orbital basis, given the 1-RDM, 2-RDM, and the 1- and 2-electron integrals
  ! Implementation of eq. A6 from JCP 141 244104 (2014)
  ! This subroutine is specific for the case of pCCD
  END_DOC

  integer         , intent(in)  :: nT
  double precision, intent(out) :: r_orbrot_h_diag(nT,nT)
  double precision, intent(in)  :: r_one_e_dm_mo(nT,nT)
  double precision, intent(in)  :: r_two_e_dm_mo(nT,nT,nT,nT)

  integer :: p, q, r, s, t, v

! Perhaps it is not worthy introducing these arrays in this case, when only the diagonal Hessian is computed 

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

    end do
    r_orbrot_h_diag(q,q) = r_orbrot_h_diag(q,q) - 2.0d0 * mo_one_e_integrals(q,q) * r_one_e_dm_mo(q,q) 
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
    r_orbrot_h_diag(p,p) = r_orbrot_h_diag(p,p) + 2.0d0 * mo_one_e_integrals(p,p) * r_one_e_dm_mo(p,p)
  end do

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

end subroutine compute_r_orbrot_h_diagonal_pccd
