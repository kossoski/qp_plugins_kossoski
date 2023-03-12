subroutine lapack_exp_antisymm_matrix(m,A,C)
  implicit none
  BEGIN_DOC
  ! Computes C = exp(A), where A is a real anti-symmetric matrix ( A^\dag = A^T = - A )
  END_DOC
  integer         , intent(in)  :: m
  double precision, intent(in)  :: A(m,m)
  double precision, intent(out) :: C(m,m)

  complex*16, allocatable       :: iA(:,:)
  double precision, allocatable :: eigvalues(:)
  complex*16, allocatable       :: eigvectors(:,:)
  complex*16                    :: img
  integer                       :: i, j

  img = dcmplx( 0.d0, 1.d0 )

  allocate( iA(m,m) )
  allocate( eigvalues(m), eigvectors(m,m) )

! Since A is real anti-symmetric, iA is hermitian and can be diagonalized with the zheevd Lapack routine
  do j=1,m
    iA(:,j) = A(:,j) * img
  end do

  call lapack_diag_hermitian_matrix(iA,m,eigvectors,eigvalues)

! Print eigenvalues:
! write(*,*) 'eigenvalues of the real anti-symmetric matrix:'
! do i=1,m
!   write(*,*) i, eigvalues(i)
! end do

! iA will be used as a temporary matrix:
! First compute U exp(-i lambda):
  do j=1,m
    do i=1,m
      iA(i,j) = eigvectors(i,j) * exp( -img * eigvalues(j) )
    end do
  end do
  deallocate( eigvalues )

! And then complete with U^\dag on the right:
! eigvectors = dconjg( transpose(eigvectors) )
  do j=1,m
    do i=1,m
      C(i,j) = real( sum( iA(i,:) * dconjg( eigvectors(j,:) ) ) )
    end do
  end do

  deallocate( iA, eigvectors )

end subroutine lapack_exp_antisymm_matrix


subroutine lapack_diag_hermitian_matrix(A,m,C,w)
  implicit none
  BEGIN_DOC
  ! Diagonalizes the complex hermitian matrix A, returning the eigenvalues in w and the eigenvectors in C
  END_DOC
  integer         , intent(in)  :: m
  complex*16      , intent(in)  :: A(m,m)
  complex*16      , intent(out) :: C(m,m)
  double precision, intent(out) :: w(m)

  complex*16, allocatable       :: work(:)
  integer                       :: lwork
  double precision, allocatable :: rwork(:)
  integer                       :: lrwork
  integer, allocatable          :: iwork(:)
  integer                       :: liwork
  integer                       :: info

  complex*16, allocatable       :: Acopy(:,:)

  allocate( Acopy(m,m) )
  Acopy = A

! First call is to obtain the optial working spaces (lwork, lrwork, and liwork)

  lwork = -1
  lrwork = -1
  liwork = -1
  allocate( work(1), rwork(1), iwork(1) )

! call zheevd('V', 'U', m, A, m, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  call zheevd('V', 'U', m, Acopy, m, w, work, lwork, rwork, lrwork, iwork, liwork, info)

  if (info < 0) then
    print *, irp_here, ': zheevd: the ',-info,'-th argument had an illegal value'
    stop 2
  endif

  lwork  = int( real(work(1)) )
  lrwork = int( rwork(1) )
  liwork = iwork(1)
  deallocate( work, rwork, iwork )
  allocate( work(lwork), rwork(lrwork), iwork(liwork) )

! Second call is for real

! call zheevd('V', 'U', m, A, m, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  call zheevd('V', 'U', m, Acopy, m, w, work, lwork, rwork, lrwork, iwork, liwork, info)

  if (info < 0) then
    print *, irp_here, ': zheevd: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0 ) then
    write(*,*) 'zheevd Failed'
    stop 1
  end if

! C = A
  C = Acopy

  deallocate( Acopy )
  deallocate( work, rwork, iwork )

end subroutine lapack_diag_hermitian_matrix


subroutine lapack_solve_linear_dsy(m,A,B,X)
  implicit none
  BEGIN_DOC
  ! Solve the linear system AX=B, for real A and B, returning the solution X
  END_DOC
  integer         , intent(in)  :: m
  double precision, intent(in)  :: A(m,m)
  double precision, intent(in)  :: B(m)
  double precision, intent(out) :: X(m)
! double precision, intent(inout) :: X(m)

  integer, allocatable          :: ipiv(:)
  double precision, allocatable :: work(:)
  integer                       :: lwork
  integer                       :: info

  double precision, allocatable :: Acopy(:,:)
  double precision, allocatable :: Bcopy(:)

  allocate( Acopy(m,m) )
  allocate( Bcopy(m) )

  Acopy(:,:) = A(:,:)
  Bcopy(:) = B(:)
! X = B

! First call is to obtain the optial working spaces (lwork)

  allocate( ipiv(m) )
  lwork = -1
  allocate( work(1) )

! call dsysv('U', m, m, A, m, ipiv, X, m, work, lwork, info)
  call dsysv('U', m, m, Acopy, m, ipiv, Bcopy, m, work, lwork, info)

  if (info < 0) then
    print *, irp_here, ': dsysv: the ',-info,'-th argument had an illegal value'
    stop 2
  endif

  lwork  = int( work(1) )
  deallocate( work )
  allocate( work(lwork) )

! Second call is for real

! call dsysv('U', m, m, A, m, ipiv, X, m, work, lwork, info)
! call dsysv('U', m, m, Acopy, m, ipiv, Bcopy, m, work, lwork, info)
  call dsysv('U', m, m, Acopy, m, ipiv, Bcopy, m, work, lwork, info)

  if (info < 0) then
    print *, irp_here, ': dsysv: the ',-info,'-th argument had an illegal value'
    stop 2
  else if( info > 0 ) then
    write(*,*) 'dsysv Failed'
    stop 1
  end if

  X(:) = Bcopy(:)

  deallocate( Acopy, Bcopy )

end subroutine lapack_solve_linear_dsy

subroutine svd_all(A,LDA,U,LDU,D,Vt,LDVt,m,n)
  implicit none
  BEGIN_DOC
  ! Exactly the same routine from utils/linear_algebra.irp.f,
  ! but replacing 'S' by 'A' as arguments of dgesvd
  !
  ! Compute A = U.D.Vt
  !
  ! LDx : leftmost dimension of x
  !
  ! Dimension of A is m x n
  !
  END_DOC

  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  double precision, intent(in)    :: A(LDA,n)
  double precision, intent(out)   :: U(LDU,min(m,n))
  double precision,intent(out)    :: Vt(LDVt,n)
  double precision,intent(out)    :: D(min(m,n))
  double precision,allocatable    :: work(:)
  integer                         :: info, lwork, i, j, k

  double precision,allocatable    :: A_tmp(:,:)
  allocate (A_tmp(LDA,n))
  do k=1,n
    do i=1,m
      A_tmp(i,k) = A(i,k)
    enddo
  enddo

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  lwork = max(int(work(1)), 10*MIN(M,N))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, info)
  deallocate(A_tmp,work)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

  do j=1,min(m,n)
    do i=1,m
      if (dabs(U(i,j)) < 1.d-14)  U(i,j) = 0.d0
    enddo
  enddo

  do j=1,n
    do i=1,n
      if (dabs(Vt(i,j)) < 1.d-14) Vt(i,j) = 0.d0
    enddo
  enddo

end
