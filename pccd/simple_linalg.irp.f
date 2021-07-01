subroutine add_shift_to_matrix(A,dim1,dim2,shift)
  implicit none
  BEGIN_DOC
  ! Adds a constant "shift" to all the elements of matrix A
  END_DOC
  integer         , intent(in)    :: dim1, dim2
  double precision, intent(inout) :: A(dim1,dim2)
  double precision, intent(in)    :: shift
  integer :: i, j
  do j=1,dim2
    do i=1,dim1
      A(i,j) = A(i,j) + shift
    end do
  end do
end subroutine add_shift_to_matrix


subroutine add_shift_to_diagonal(A,n,shift)
  implicit none
  BEGIN_DOC
  ! Adds a constant "shift" to the diagonal elements of matrix A
  END_DOC
  integer         , intent(in)    :: n
  double precision, intent(inout) :: A(n,n)
  double precision, intent(in)    :: shift
  integer :: i, j
  do i=1,n
    A(i,i) = A(i,i) + shift
  end do
end subroutine add_shift_to_diagonal


function index_lowtri_to_square(p,q,m) result(r)
  implicit none
  BEGIN_DOC
  ! Given a square matrix A(m,m), and a vector B formed from the lower diagonal elements of A (in column-major order),
  ! this function gives the index of B corresponding to the element A(p,q)
  END_DOC
  integer :: p, q, m, r
  r =  (2*m-q)*(q-1)/2 + p-q
end function index_lowtri_to_square


subroutine transform_2index_to_1index_square(A_2index,A_1index,dim1,dim2)
  implicit none
  BEGIN_DOC
  ! Transforms a matrix A_2index into a vector A_1index formed from all the elements of A_2index (in column-major order)
  END_DOC
  integer,          intent(in)  :: dim1, dim2
  double precision, intent(in)  :: A_2index(dim1,dim2)
  double precision, intent(out) :: A_1index(dim1*dim2)
  integer :: p, q, pq
  do q=1,dim2
    do p=1,dim1
      pq = (q-1)*dim1 + p
      A_1index(pq) = A_2index(p,q)
    end do
  end do
end subroutine transform_2index_to_1index_square


subroutine transform_2index_to_1index_lowtri(A_2index,A_1index,dim1)
  implicit none
  BEGIN_DOC
  ! Transforms a square matrix A_2index into a vector A_1index formed from the lower diagonal elements of A_2index (in column-major order)
  END_DOC
  integer,          intent(in)  :: dim1
  double precision, intent(in)  :: A_2index(dim1,dim1)
  double precision, intent(out) :: A_1index(dim1*(dim1-1)/2)
  integer :: p, q, pq, index_lowtri_to_square
  do q=1,dim1-1
    do p=q+1,dim1
      pq = index_lowtri_to_square(p,q,dim1)
      A_1index(pq) = A_2index(p,q)
    end do
  end do
end subroutine transform_2index_to_1index_lowtri


subroutine transform_1index_to_2index_square(A_1index,A_2index,dim1,dim2)
  implicit none
  BEGIN_DOC
  ! Transforms a vector A_1index into a square matrix A_2index (in column-major order)
  END_DOC
  integer,          intent(in)  :: dim1, dim2
  double precision, intent(in)  :: A_1index(dim1*dim2)
  double precision, intent(out) :: A_2index(dim1,dim2)
  integer :: p, q, pq
  do q=1,dim2
    do p=1,dim1
      pq = (q-1)*dim1 + p
      A_2index(p,q) = A_1index(pq)
    end do
  end do
end subroutine transform_1index_to_2index_square


subroutine transform_4index_to_2index_square(A_4index,A_2index,dim1,dim2,dim3,dim4)
  implicit none
  BEGIN_DOC
  ! Transforms all the elements of a 4-dim array A_4index into a 2-dim array A_2index (in column-major order)
  END_DOC
  integer,          intent(in)  :: dim1, dim2, dim3, dim4
  double precision, intent(in)  :: A_4index(dim1,dim2,dim3,dim4)
  double precision, intent(out) :: A_2index(dim1*dim2,dim3*dim4)
  integer :: p, q, r, s, pq, rs
  do s=1,dim4
    do r=1,dim3
      rs = (s-1)*dim3 + r
      do q=1,dim2
        do p=1,dim1
          pq = (q-1)*dim1 + p
          A_2index(pq,rs) = A_4index(p,q,r,s)
        end do
      end do
    end do
  end do
end subroutine transform_4index_to_2index_square


subroutine transform_4index_to_2index_lowtri(A_4index,A_2index,dim1,dim3)
  implicit none
  BEGIN_DOC
  ! Transforms the lower diagonal pairwise elements of a 4-dim array A_4index into a 2-dim array A_2index (in column-major order)
  END_DOC
  integer,          intent(in)  :: dim1, dim3
  double precision, intent(in)  :: A_4index(dim1,dim1,dim3,dim3)
  double precision, intent(out) :: A_2index(dim1*(dim1-1)/2,dim3*(dim3-1)/2)
  integer :: p, q, r, s, pq, rs, index_lowtri_to_square
  do s=1,dim3-1
    do r=s+1,dim3
      rs = index_lowtri_to_square(r,s,dim3)
      do q=1,dim1-1
        do p=q+1,dim1
          pq = index_lowtri_to_square(p,q,dim1)
          A_2index(pq,rs) = A_4index(p,q,r,s)
        end do
      end do
    end do
  end do
end subroutine transform_4index_to_2index_lowtri


subroutine zero_small_2index(A,dim1,dim2,threshold)
  implicit none
  BEGIN_DOC
  ! Assigns zero to matrix elements below a given threshold
  END_DOC
  integer,          intent(in)    :: dim1, dim2
  double precision, intent(inout) :: A(dim1,dim2)
  double precision, intent(in)    :: threshold
  integer :: p, q
  do q=1,dim2
    do p=1,dim1
      if( abs(A(p,q)).lt.threshold ) A(p,q) = 0.d0
    end do
  end do
end subroutine zero_small_2index


subroutine zero_small_4index(A,dim1,dim2,dim3,dim4,threshold)
  implicit none
  BEGIN_DOC
  ! Assigns zero to matrix elements below a given threshold
  END_DOC
  integer,          intent(in)    :: dim1, dim2, dim3, dim4
  double precision, intent(inout) :: A(dim1,dim2,dim3,dim4)
  double precision, intent(in)    :: threshold
  integer :: p, q, r, s
  do s=1,dim4
    do r=1,dim3
      do q=1,dim2
        do p=1,dim1
          if( abs(A(p,q,r,s)).lt.threshold ) A(p,q,r,s) = 0.d0
        end do
      end do
    end do
  end do
end subroutine zero_small_4index


subroutine form_A_square(A_lowtri,m_lowtri,A_square,m_square)
  implicit none
  BEGIN_DOC
  ! Transforms a vector A_lowtri into an antisymmetric matrix A (in column-major order)
  END_DOC
  integer,          intent(in)  :: m_lowtri
  double precision, intent(in)  :: A_lowtri(m_lowtri)
  integer,          intent(in)  :: m_square
  double precision, intent(out) :: A_square(m_square,m_square)
  integer :: p, q, r
  do q=1,m_square
    A_square(q,q) = 0.d0
  end do
  do q=1,m_square-1
    do p=q+1,m_square
      r = (2*m_square-q)*(q-1)/2 + p-q
      A_square(p,q) =   A_lowtri(r)
      A_square(q,p) = - A_lowtri(r)
    end do
  end do
end subroutine form_A_square


subroutine merge_identity_to_matrix(A_in,A_out,n_in,n_out)
  implicit none
  BEGIN_DOC
  ! Appends an indentity matrix to A_in, giving A_out:
  ! A_out = | I 0    |
  !         | 0 A_in |
  END_DOC
  integer,          intent(in)   :: n_in
  integer,          intent(in)   :: n_out
  double precision, intent(in)   :: A_in(n_in,n_in)
  double precision, intent(out)  :: A_out(n_out,n_out)
  integer :: p, q, n_diff
  n_diff = n_out - n_in
  A_out = 0.0d0   
  do p=1,n_diff
    A_out(p,p) = 1.0d0
  end do
  do q=1,n_in
    do p=1,n_in
      A_out(n_diff+p,n_diff+q) = A_in(p,q)
    end do
  end do
end subroutine merge_identity_to_matrix


subroutine check_matrix_symmetry(A,m)
  implicit none
  BEGIN_DOC
  ! Writes whether the matrix A is symmetric or not
  END_DOC
  integer,          intent(in) :: m
  double precision, intent(in) :: A(m,m)
  logical                      :: symmetric = .true.
  double precision, parameter  :: thresh = 1.d-10
  integer                      :: i, j
  outer: do j=1,m-1
    do i=j+1,m
      if( abs( A(i,j) - A(j,i) ) .gt. thresh ) then
        symmetric = .false.
        exit outer
      end if
    end do
  end do outer
  if( symmetric ) then
     write(*,*) 'Matrix is symmetric'
  else
     write(*,*) 'WARNING: matrix is not symmetric'
  end if
end subroutine check_matrix_symmetry


subroutine prepend(N,M,A,b)

  implicit none
  BEGIN_DOC
  ! Prepend the vector b of size N into the matrix A of size NxM
  END_DOC

! Input variables
  integer,intent(in)            :: N,M
  double precision,intent(in)   :: b(N)

! Output variables
! double precision,intent(out)  :: A(N,M)
  double precision,intent(inout):: A(N,M)

! Local variables
  integer                       :: i,j

! print*,'b in append'
! call matout(N,1,b)

  do i=1,N
    do j=M-1,1,-1
      A(i,j+1) = A(i,j)
    enddo
    A(i,1) = b(i)
  enddo

end subroutine prepend


subroutine append(N,M,A,b)

  implicit none
  BEGIN_DOC
  ! Append the vector b of size N into the matrix A of size NxM
  END_DOC

! Input variables
  integer,intent(in)            :: N,M
  double precision,intent(in)   :: b(N)

! Output variables
  double precision,intent(inout):: A(N,M)

! Local viaruabkes
  integer                       :: i,j

  do i=1,N
    do j=2,M
      A(i,j-1) = A(i,j)
    enddo
    A(i,M) = b(i)
  enddo

end subroutine append


subroutine build_augmented_matrix(A,b,n,C)
  implicit none
  BEGIN_DOC
  ! Builds the matrix C as:
  ! C = | A b |
  !     | b 0 |
  END_DOC
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  double precision, intent(in)  :: b(n)
  double precision, intent(out) :: C(n+1,n+1)
  integer :: i, j
  do j=1,n
    do i=1,n
      C(i,j) = A(i,j)
    end do
  end do
  do i=1,n
    C(n+1,i) = b(i)
    C(i,n+1) = b(i)
  end do
  C(n+1,n+1) = 0.0d0
end subroutine build_augmented_matrix


subroutine solve_polynomial_2nd(c0,c1,c2,root1,root2)
  implicit none
  BEGIN_DOC
  ! Solves a 2nd order polynomial equation
  END_DOC
  double precision, intent(in)  :: c0, c1, c2
  complex*16      , intent(out) :: root1, root2
  double precision :: delta
  delta = c1**2 - 4.0d0 * c2 * c0
  if( delta.eq.0.0d0 ) then
    root1 = - c1 / (2.0d0 * c2)
    root2 = root1
  elseif( delta.gt.0.0d0 ) then
    root1 = (- c1 - sqrt(delta) ) / (2.0d0 * c2 )
    root2 = (- c1 + sqrt(delta) ) / (2.0d0 * c2 )
  else
    root1 = (- c1 + sqrt(dcmplx(delta)) ) / (2.0d0 * c2 )
    root2 = dconjg( root1 )
  end if
end subroutine solve_polynomial_2nd


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

! iA will be used as a temporary matrix:
! First compute U exp(-i lambda):
  do j=1,m
    do i=1,m
      iA(i,j) = eigvectors(i,j) * exp( -img * eigvalues(j) )
    end do
  end do
  deallocate( eigvalues )

! And then complete with U^\dag on the right:
  eigvectors = dconjg( transpose(eigvectors) )
  do j=1,m
    do i=1,m
      C(i,j) = real( sum( iA(i,:) * eigvectors(:,j) ) )
    end do
  end do

! Print eigenvalues:
! write(*,*) 'eigenvalues:'
! do i=1,m
!   write(*,*) i, eigvalues(i)
! end do

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
