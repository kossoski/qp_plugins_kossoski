subroutine form_delta_OV(nOa,nV,r_eO,r_eV,delta)
  implicit none
  BEGIN_DOC
  ! Computes orbital energy differences for pCCD
  END_DOC

! Input variables
  integer,         intent(in)   :: nOa, nV
  double precision,intent(in)   :: r_eO(nOa)
  double precision,intent(in)   :: r_eV(nV)
! Output variables
  double precision,intent(out)  :: delta(nOa,nV)

! Local variables
  integer                       :: i,a

  do a=1,nV
    do i=1,nOa
      delta(i,a) = 2.0d0 * ( r_eV(a) - r_eO(i) )
    enddo
  enddo

end subroutine form_delta_OV
