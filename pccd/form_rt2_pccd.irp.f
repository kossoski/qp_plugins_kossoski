subroutine form_rt2_pccd(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,rt2)

  implicit none
  BEGIN_DOC
  ! Computes the residuals for the pCCD t-amplitude equations
  END_DOC

! Input variables
  integer, intent(in)           :: nOa, nV
  double precision,intent(in)   :: r_OOVV(nOa,nV)
  double precision,intent(in)   :: r_OVOV(nOa,nV)
  double precision,intent(in)   :: r_VVVV(nV,nV)
  double precision,intent(in)   :: r_OOOO(nOa,nOa)
  double precision,intent(in)   :: r_OVVO(nOa,nV)
  double precision,intent(in)   :: r_delta_OV(nOa,nV)
  double precision,intent(in)   :: t2(nOa,nV)

! Output variables
  double precision,intent(out)  :: rt2(nOa,nV)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision,allocatable  :: y(:,:)

  rt2(:,:) = 0.0d0

  allocate( y(nOa,nOa) )
  y = 0.0d0
  do j=1,nOa
    do i=1,nOa
      do b=1,nV
        y(i,j) = y(i,j) + r_OOVV(j,b) * t2(i,b) 
      end do
    end do
  end do

  do a=1,nV
    do i=1,nOa

!     rt2(i,a) = r_OOVV(i,a)
!     rt2(i,a) = rt2(i,a) + r_delta_OV(i,a) * t2(i,a)
!     rt2(i,a) = rt2(i,a) - 2.0d0 * &
!                ( 2.0d0 * r_OVOV(i,a) - r_OVVO(i,a) - r_OOVV(i,a) * t2(i,a) ) * t2(i,a)
      rt2(i,a) = r_OOVV(i,a) + ( r_delta_OV(i,a) - 2.0d0 * &
                 ( 2.0d0 * r_OVOV(i,a) - r_OVVO(i,a) - r_OOVV(i,a) * t2(i,a) ) ) * t2(i,a)

!     do j=1,nOa
!         rt2(i,a) = rt2(i,a) - 2.0d0 * r_OOVV(j,a) * t2(j,a) * t2(i,a)
!     end do
!     do j=1,nOa
!         rt2(i,a) = rt2(i,a) + ( r_OOOO(i,j) + y(i,j) ) * t2(j,a)
!     end do
      do j=1,nOa
          rt2(i,a) = rt2(i,a) + ( ( r_OOOO(i,j) + y(i,j) ) - 2.0d0 * r_OOVV(j,a) * t2(i,a) ) * t2(j,a)
      end do

!     do b=1,nV
!         rt2(i,a) = rt2(i,a) - 2.0d0 * r_OOVV(i,b) * t2(i,b) * t2(i,a)
!     end do
!     do b=1,nV
!         rt2(i,a) = rt2(i,a) + r_VVVV(b,a) * t2(i,b)
!     end do
      do b=1,nV
          rt2(i,a) = rt2(i,a) + ( r_VVVV(b,a) - 2.0d0 * r_OOVV(i,b) * t2(i,a) ) * t2(i,b)
      end do

    end do
  end do

  deallocate( y )

end subroutine form_rt2_pccd
