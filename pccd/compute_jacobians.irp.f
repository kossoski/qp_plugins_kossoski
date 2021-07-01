subroutine compute_pccd_t2_jacobian(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,t2_jacobian)

  implicit none
  BEGIN_DOC
  ! Computes the Jacobian matrix of the pCCD t-amplitude residual equations, t2_jacobian(i,a,j,b) = J_{ia,jb} = \partial r_i^a / \partial t_j^b 
  END_DOC

! Input variables
  integer,         intent(in)   :: nOa, nV
  double precision,intent(in)   :: r_OOVV(nOa,nV)
  double precision,intent(in)   :: r_OVOV(nOa,nV)
  double precision,intent(in)   :: r_VVVV(nV,nV)
  double precision,intent(in)   :: r_OOOO(nOa,nOa)
  double precision,intent(in)   :: r_OVVO(nOa,nV)
  double precision,intent(in)   :: r_delta_OV(nOa,nV)
  double precision,intent(in)   :: t2(nOa,nV)

! Output variables
  double precision,intent(out)  :: t2_jacobian(nOa,nV,nOa,nV)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision,allocatable  :: yO1(:)
  double precision,allocatable  :: yV1(:)
  double precision,allocatable  :: yO2(:,:)
  double precision,allocatable  :: yV2(:,:)

  t2_jacobian = 0.0d0

  allocate( yO1(nOa) )
  yO1 = 0.0d0
  do i=1,nOa
    yO1(i) = sum( r_OOVV(i,:) * t2(i,:) )
  end do

  allocate( yV1(nV) )
  yV1 = 0.0d0
  do a=1,nV
    yV1(a) = sum( r_OOVV(:,a) * t2(:,a) )
  end do

  allocate( yO2(nOa,nOa) )
  yO2 = 0.0d0
  do j=1,nOa
    do i=1,nOa
      do b=1,nV
        yO2(i,j) = yO2(i,j) + r_OOVV(j,b) * t2(i,b)
      end do
    end do
  end do

  allocate( yV2(nV,nV) )
  yV2 = 0.0d0
  do b=1,nV
    do a=1,nV
      do j=1,nOa
        yV2(a,b) = yV2(a,b) + r_OOVV(j,b) * t2(j,a)
      end do
    end do
  end do

  do a=1,nV
    do i=1,nOa
      t2_jacobian(i,a,i,a) = r_delta_OV(i,a) &
      - 4.0d0*r_OVOV(i,a) + 2.0d0*r_OVVO(i,a) + r_VVVV(a,a) + r_OOOO(i,i) &
      - yV1(a) - yO1(i) 
    end do
  end do

  do b=1,nV
    do a=1,nV
    if( a.eq.b ) cycle
      do i=1,nOa
        t2_jacobian(i,a,i,b) = r_VVVV(a,b) - 2.0d0 * r_OOVV(i,b) * t2(i,a) + yV2(a,b)
      end do
    end do
  end do

  do j=1,nOa
    do i=1,nOa
    if( i.eq.j ) cycle
      do a=1,nV
        t2_jacobian(i,a,j,a) = r_OOOO(i,j) - 2.0d0 * r_OOVV(j,a) * t2(i,a) + yO2(i,j)
      end do
    end do
  end do

end subroutine compute_pccd_t2_jacobian


subroutine compute_pccd_t2_jacobian_diagonal(nOa,nV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2_jacobian_diagonal)

  implicit none
  BEGIN_DOC
  ! Computes the diagonal Jacobian matrix of the pCCD t-amplitude residual equations, t2_jacobian_diagonal(i,a) = J_{ia,ia} = \partial r_i^a / \partial t_i^a
  END_DOC

! Input variables
  integer,         intent(in)   :: nOa, nV
  double precision,intent(in)   :: r_OVOV(nOa,nV)
  double precision,intent(in)   :: r_VVVV(nV,nV)
  double precision,intent(in)   :: r_OOOO(nOa,nOa)
  double precision,intent(in)   :: r_OVVO(nOa,nV)
  double precision,intent(in)   :: r_delta_OV(nOa,nV)

! Output variables
  double precision,intent(out)  :: t2_jacobian_diagonal(nOa,nV)

! Local variables
  integer                       :: i,a

  do a=1,nV
    do i=1,nOa
      t2_jacobian_diagonal(i,a) = r_delta_OV(i,a) &
      - 4.0d0*r_OVOV(i,a) + 2.0d0*r_OVVO(i,a) + r_VVVV(a,a) + r_OOOO(i,i)
    end do
  end do

end subroutine compute_pccd_t2_jacobian_diagonal


subroutine compute_pccd_z2_jacobian_diagonal(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2_jacobian_diagonal)

  implicit none
  BEGIN_DOC
  ! Computes the diagonal Jacobian matrix of the pCCD z-amplitude residual equations, z2_jacobian_diagonal(i,a) = J_{ia,ia} = \partial r_i^a / \partial t_i^a
  END_DOC

! Input variables
  integer,         intent(in)   :: nOa, nV
  double precision,intent(in)   :: r_OOVV(nOa,nV)
  double precision,intent(in)   :: r_OVOV(nOa,nV)
  double precision,intent(in)   :: r_VVVV(nV,nV)
  double precision,intent(in)   :: r_OOOO(nOa,nOa)
  double precision,intent(in)   :: r_OVVO(nOa,nV)
  double precision,intent(in)   :: r_delta_OV(nOa,nV)
  double precision,intent(in)   :: t2(nOa,nV)

! Output variables
  double precision,intent(out)  :: z2_jacobian_diagonal(nOa,nV)

! Local variables
  integer                       :: i,a

  do a=1,nV
    do i=1,nOa
      z2_jacobian_diagonal(i,a) = r_delta_OV(i,a) &
      - 4.0d0*r_OVOV(i,a) + 2.0d0*r_OVVO(i,a) + r_VVVV(a,a) + r_OOOO(i,i) &
      - sum( r_OOVV(i,:)*t2(i,:) ) - sum( r_OOVV(:,a)*t2(:,a) )
    end do
  end do

end subroutine compute_pccd_z2_jacobian_diagonal


subroutine compute_pccd_z2_jacobian(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2_jacobian)

  implicit none
  BEGIN_DOC
  ! Computes the Jacobian matrix of the pCCD z-amplitude residual equations, t2_jacobian(i,a,j,b) = J_{ia,jb} = \partial r_i^a / \partial t_j^b 
  END_DOC

! Input variables
  integer,         intent(in)   :: nOa, nV
  double precision,intent(in)   :: r_OOVV(nOa,nV)
  double precision,intent(in)   :: r_OVOV(nOa,nV)
  double precision,intent(in)   :: r_VVVV(nV,nV)
  double precision,intent(in)   :: r_OOOO(nOa,nOa)
  double precision,intent(in)   :: r_OVVO(nOa,nV)
  double precision,intent(in)   :: r_delta_OV(nOa,nV)
  double precision,intent(in)   :: t2(nOa,nV)

! Output variables
  double precision,intent(out)  :: z2_jacobian(nOa,nV,nOa,nV)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision,allocatable  :: yO2(:,:)
  double precision,allocatable  :: yV2(:,:)

  z2_jacobian = 0.0d0

  allocate( yO2(nOa,nOa) )
  yO2 = 0.0d0
  do j=1,nOa
    do i=1,nOa
      do b=1,nV
        yO2(i,j) = yO2(i,j) + r_OOVV(i,b) * t2(j,b)
      end do
    end do
  end do
  allocate( yV2(nV,nV) )
  yV2 = 0.0d0
  do b=1,nV
    do a=1,nV
      do j=1,nOa
        yV2(a,b) = yV2(a,b) + r_OOVV(j,a) * t2(j,b)
      end do
    end do
  end do

  do a=1,nV
    do i=1,nOa
        z2_jacobian(i,a,i,a) = r_delta_OV(i,a) &
        - 4.0d0*r_OVOV(i,a) + 2.0d0*r_OVVO(i,a) + r_VVVV(a,a) + r_OOOO(i,i) &
        - sum( r_OOVV(i,:)*t2(i,:) ) - sum( r_OOVV(:,a)*t2(:,a) )
    end do
  end do

  do b=1,nV
    do a=1,nV
    if( a.eq.b ) cycle
      do i=1,nOa
        z2_jacobian(i,a,i,b) = r_VVVV(a,b) - 2.0d0 * r_OOVV(i,a) * t2(i,b) + yV2(a,b)
      end do
    end do
  end do

  do j=1,nOa
    do i=1,nOa
    if( i.eq.j ) cycle
      do a=1,nV
        z2_jacobian(i,a,j,a) = r_OOOO(i,j) - 2.0d0 * r_OOVV(i,a) * t2(j,a) + yO2(i,j)
      end do
    end do
  end do

end subroutine compute_pccd_z2_jacobian
