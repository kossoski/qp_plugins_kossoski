subroutine DIIS_extrapolation(rcond,n_err,n_e,n_diis,error,e,error_in,e_inout)

  implicit none
  BEGIN_DOC
  ! Performs DIIS extrapolation
  ! This is essentially the same subroutine as implemented in Quack
  END_DOC

! include 'parameters.h'

! Output variables
  double precision,intent(out)  :: rcond

! Input variables
  integer,intent(in)            :: n_diis
  integer,intent(in)            :: n_err,n_e
  double precision,intent(in)   :: error_in(n_err)
  double precision,intent(inout):: error(n_err,n_diis),e(n_e,n_diis)
! double precision,intent(in)   :: error_in(n_err),error(n_err,n_diis),e(n_e,n_diis)
! double precision,intent(out)  :: error(n_err,n_diis),e(n_e,n_diis)

! Output variables
  double precision,intent(inout):: e_inout(n_e)

! Local variables
  double precision,allocatable  :: A(:,:),b(:),w(:)


! Memory allocation

  allocate(A(n_diis+1,n_diis+1),b(n_diis+1),w(n_diis+1))

! Update DIIS "history"

  call prepend(n_err,n_diis,error,error_in)

  call prepend(n_e,n_diis,e,e_inout)

! Build A matrix

  A(1:n_diis,1:n_diis) = matmul(transpose(error),error)

  A(1:n_diis,n_diis+1) = -1.0d0
  A(n_diis+1,1:n_diis) = -1.0d0
  A(n_diis+1,n_diis+1) = +0.0d0

! Build x matrix

  b(1:n_diis) = +0.0d0
  b(n_diis+1) = -1.0d0

! Solve linear system

  call linear_solve(n_diis+1,A,b,w,rcond)

! Extrapolate

  e_inout(:) = matmul(w(1:n_diis),transpose(e(:,1:n_diis)))

! Print weights
! write(*,*) 'DIIS coefficients:'
! integer :: i
! do i=1,n_diis
!   write(*,*) i, w(i)
! end do

end subroutine DIIS_extrapolation


subroutine linear_solve(N,A,b,x,rcond)

  implicit none
  BEGIN_DOC
  ! Solves the linear system Ax = b, returning the solution x and the condition number (rcond) of the matrix A
  END_DOC

  integer,intent(in)             :: N
  double precision,intent(in)    :: A(N,N),b(N)!,rcond
  double precision,intent(out)   :: x(N)
  double precision,intent(out)   :: rcond

  integer                        :: info,lwork
  double precision               :: ferr,berr
  integer,allocatable            :: ipiv(:),iwork(:)
  double precision,allocatable   :: AF(:,:),work(:)

  lwork = 3*N
  allocate(AF(N,N),ipiv(N),work(lwork),iwork(N))

  call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)

! if (info /= 0) then

!   print *,  info
!   stop 'error in linear_solve (dsysvx)!!'

! endif

end subroutine linear_solve
