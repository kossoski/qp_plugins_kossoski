program c0_ref

  implicit none
  BEGIN_DOC
  ! Prints the square root of the sum of squared coefficients of the determinants in the reference space (as defined by the det_read.inp file)
  END_DOC

  integer :: i, j, k, l

  PROVIDE N_dets_read det_read

  double precision :: c0(N_states)
  c0 = 0.0d0

  do j=1,N_dets_read
    inner: do l=1,N_det
      do i=1,N_int
        if( det_read(i,1,j)/=psi_det(i,1,l) .or. det_read(i,2,j)/=psi_det(i,2,l) ) cycle inner
      end do
      write(*,*) l, ( psi_coef(l,k), k=1,N_states )
      c0(:) = c0(:) + psi_coef(l,:)**2
      exit inner
    end do inner
  end do

  c0 = sqrt( c0 )
  write(*,*) 
  write(*,*) 'c0 = [ sum_p |c_p|^2 ]^0.5, for each state:'
  write(*,*) 
  write(*,*) ( c0(k), k=1,N_states )
  write(*,*) 

! write(*,*) ( sum( psi_coef(1:N_det,k)**2 ), k=1,N_states )
! write(*,*) ( sum( psi_coef(1:10000,k)**2 ), k=1,N_states )

! do l=1,N_det
!   write(*,*) l, ( psi_coef(l,k), k=1,N_states )
! end do

end program
