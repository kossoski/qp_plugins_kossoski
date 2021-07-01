subroutine write_i_1d_array(x,p)
  implicit none
  BEGIN_DOC
  ! Writes a 1-dim array
  END_DOC
  integer,           intent(in) :: p
  double precision,  intent(in) :: x(p)
  integer :: i
  do i=1,p
    write(*,*) i, x(i)
  end do
end subroutine write_i_1d_array


subroutine write_ij_2d_array(x,p,q)
  implicit none
  BEGIN_DOC
  ! Writes a 2-dim array
  END_DOC
  integer,           intent(in) :: p, q
  double precision,  intent(in) :: x(p,q)
  integer :: i, j
  do j=1,q
    do i=1,p
      write(*,*) i, j, x(i,j)
    end do
  end do
end subroutine write_ij_2d_array


subroutine read_1d_from_file(A,m,flnm)
  implicit none
  BEGIN_DOC
  ! Reads a 1-dim array from a file
  END_DOC
  integer,          intent(in)  :: m
  double precision, intent(out) :: A(m)
  character(len=*), intent(in)  :: flnm
  integer                       :: i, nunit, stat
  integer                       :: n
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    read(nunit,*) n
    if( n.ne.m ) then
      write(*,*) 'Wrong dimension in ', trim(adjustl(flnm))
      write(*,*) n, m
      stop 
    end if
    do i=1,m
      read(nunit,*) A(i)
    end do
  close(nunit)
end subroutine read_1d_from_file


subroutine write_1d_to_file(A,m,flnm)
  implicit none
  BEGIN_DOC
  ! Writes a 1-dim array to a file
  END_DOC
  integer,          intent(in)  :: m
  double precision, intent(in)  :: A(m)
  character(len=*), intent(in)  :: flnm
  integer                       :: i, nunit, stat
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    write(nunit,*) m
    do i=1,m
      write(nunit,*) A(i)
    end do
  close(nunit)
end subroutine write_1d_to_file


subroutine read_2d_from_file(A,m,flnm)
  implicit none
  BEGIN_DOC
  ! Reads a 2-dim array from a file
  END_DOC
  integer,          intent(in)  :: m
  double precision, intent(out) :: A(m,m)
  character(len=*), intent(in)  :: flnm
  integer                       :: i, j, nunit, stat
  integer                       :: n
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    read(nunit,*) n
    if( n.ne.m ) then
      write(*,*) 'Wrong dimension in ', trim(adjustl(flnm))
      write(*,*) n, m
      stop 
    end if
    do j=1,m
      do i=1,m
        read(nunit,*) A(i,j)
      end do
    end do
  close(nunit)
end subroutine read_2d_from_file


subroutine write_2d_to_file(A,m,flnm)
  implicit none
  BEGIN_DOC
  ! Writes a 2-dim array to a file
  END_DOC
  integer,          intent(in)  :: m
  double precision, intent(in)  :: A(m,m)
  character(len=*), intent(in)  :: flnm
  integer                       :: i, j, nunit, stat
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    write(nunit,*) m
    do j=1,m
      do i=1,m
        write(nunit,*) A(i,j)
      end do
    end do
  close(nunit)
end subroutine write_2d_to_file


subroutine write_double_to_file(x,flnm)
  implicit none
  BEGIN_DOC
  ! Writes a double precision variable to a file
  END_DOC
  double precision, intent(in) :: x
  character(len=*), intent(in) :: flnm
  integer                      :: nunit, stat
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    write(nunit,*) x
  close(nunit)
end subroutine write_double_to_file


subroutine read_double_from_file(x,flnm)
  implicit none
  BEGIN_DOC
  ! Reads e a double precision variable from a file
  END_DOC
  double precision, intent(out) :: x
  character(len=*), intent(in)  :: flnm
  integer                       :: nunit, stat
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    read(nunit,*) x
  close(nunit)
end subroutine read_double_from_file


subroutine write_int_to_file(x,flnm)
  implicit none
  BEGIN_DOC
  ! Writes an integer to a file
  END_DOC
  integer,          intent(in) :: x
  character(len=*), intent(in) :: flnm
  integer                      :: nunit, stat
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    write(nunit,*) x
  close(nunit)
end subroutine write_int_to_file


subroutine read_int_from_file(x,flnm)
  implicit none
  BEGIN_DOC
  ! Reads an integer from a file
  END_DOC
  integer,          intent(out) :: x
  character(len=*), intent(in)  :: flnm
  integer                       :: nunit, stat
  open(newunit=nunit,file=flnm,form='formatted',iostat=stat)
    read(nunit,*) x
  close(nunit)
end subroutine read_int_from_file

