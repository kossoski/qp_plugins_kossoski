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


subroutine matrix_that_rotates_by_alpha_around_u(R,alpha,u)
implicit none
  BEGIN_DOC
! R = exp(A) = exp( alpha K ) = I + sin(alpha) K + (1-cos(alpha)) K**2
  END_DOC
double precision, intent(out) :: R(3,3)
double precision, intent(in)  :: alpha
double precision, intent(in)  :: u(3)

double precision :: x, y, z
double precision :: norm
double precision :: K(3,3)

norm = norm2(u)

x = u(1) / norm
y = u(2) / norm
z = u(3) / norm

K = transpose( reshape( (/ 0.d0, -z,     y, &
                           z,     0.d0, -x, &
                          -y,     x,     0.d0 /), shape(K) ) )

R = reshape( (/ 1.d0,  0.d0,  0.d0, &
                0.d0,  1.d0,  0.d0, &
                0.d0,  0.d0,  1.d0 /), shape(R) )

R = R + sin(alpha) * K + (1.0d0-cos(alpha)) * matmul(K,K)

end subroutine


subroutine set_integer_list(list,n)
  implicit none
  BEGIN_DOC
! Creates an array of integers list = (1 2 ... n)
  END_DOC
  integer, intent(in)  :: n
  integer, intent(out) :: list(n)
  integer              :: i
  do i=1,n
    list(i) = i
  end do
end subroutine


subroutine set_identity_matrix(A,n)
  implicit none
  BEGIN_DOC
! Creates an identity matrix of size A(n,n)
  END_DOC
  integer,          intent(in)  :: n
  double precision, intent(out) :: A(n,n)
  integer                       :: i
  A = 0.0d0
  do i=1,n
    A(i,i) = 1.0d0
  end do
end subroutine


subroutine update_molecular_orbitals(d,p,q,U)
  implicit none
  BEGIN_DOC
! Update the molecular orbitals d by applying a rotation matrix U
  END_DOC
  integer,          intent(in)    :: p, q
  double precision, intent(inout) :: d(p,q)
  double precision, intent(in)    :: U(q,q)
  double precision, allocatable   :: tmp(:,:)
  integer                         :: i, j
  allocate( tmp(p,q) )
  tmp = d
  do j=1,q
    do i=1,p
      d(i,j) = sum ( tmp(i,:) * U(:,j) )
    end do
  end do
  deallocate( tmp )
end subroutine update_molecular_orbitals


subroutine evaluate_orbital_overlap(mo_coef_1,mo_coef_2,orbital_overlap)

  implicit none
  BEGIN_DOC
  !
  END_DOC
  double precision, intent(in)  :: mo_coef_1(ao_num,mo_num)
  double precision, intent(in)  :: mo_coef_2(ao_num,mo_num)
  double precision, intent(out) :: orbital_overlap(mo_num,mo_num)
  double precision, allocatable :: T(:,:)

  allocate ( T(ao_num,mo_num) )

  call dgemm('N','N', ao_num, mo_num, ao_num,               &
      1.d0, ao_overlap, ao_num,                             &
      mo_coef_2, ao_num,                                    &
      0.d0, T, ao_num )

  call dgemm('T','N', mo_num, mo_num, ao_num,               &
      1.d0, mo_coef_1, ao_num,                              &
      T, ao_num,                                            &
      0.d0, orbital_overlap, mo_num)

  deallocate(T)

end subroutine


subroutine check_mos_orthonormality
  implicit none
  BEGIN_DOC
! Writes the maximum error of the expected MOs overlap
  END_DOC
  double precision, parameter :: thresh = 1.d-10
  integer                     :: i, j
  logical                     :: orthonormal = .true.
  double precision            :: max_error
  max_error = 0.0d0
  do i=1,mo_num
    max_error = max( max_error, abs( mo_overlap(i,i)-1.0d0 ) )
  end do
  do i=1,mo_num-1
    do j=i+1,mo_num
      max_error = max( max_error, abs( mo_overlap(i,j) ) )
    end do
  end do
  write(*,*) 'Maximum error in MOs overlap: ', max_error
  if( max_error .gt. thresh ) orthonormal = .false.

  if( orthonormal ) then
     write(*,*) 'MOs are orthonormal'
  else
     write(*,*) 'WARNING: large error in the orthonormality of MOs'
  end if
end subroutine check_mos_orthonormality


subroutine check_orthornormality_rotation_matrix(A,n,ok)
implicit none
  BEGIN_DOC
! Check the orthornormality of matrix A(n,n)
  END_DOC
integer,          intent(in)  :: n
double precision, intent(in)  :: A(n,n)
logical,          intent(out) :: ok

double precision              :: AAt(n,n)
double precision              :: max_error
double precision, parameter   :: thresh = 1.0d-12
integer                       :: i, j

ok = .true.

call set_identity_matrix(AAt,n)

call dgemm('N','T',n,n,n,1.0d0,A,size(A,1),A,size(A,1),-1.0d0,AAt,size(AAt,1))

max_error = 0.0d0
do j=1,n
  do i=1,n
    max_error = max( max_error, abs(AAt(i,j)) )
  enddo
enddo

if (abs(max_error) > thresh ) then
  write(*,*) 'WARNING: too large matrix element in R.R^T - 1: ', abs(max_error)
  ok = .false.
endif
write(*,*) 'Largest matrix element in R.R^T - 1: ', abs(max_error)

end subroutine

