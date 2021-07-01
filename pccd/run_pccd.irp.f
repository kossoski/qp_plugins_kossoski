subroutine run_pccd(t2,z2,nOa,nV)

  implicit none

  BEGIN_DOC
! Runs coupled-cluster with paired doubles (pCCD)
  END_DOC

  integer         ,intent(in)   :: nOa, nV
  double precision,intent(out)  :: t2(nOa,nV)
  double precision,intent(out)  :: z2(nOa,nV)

! Local variables
  double precision              :: ETHF
  double precision,allocatable  :: eHF(:)
  integer                       :: i,j,a,b,ia,p,q,r,s
  integer                       :: nO
  integer                       :: npairs
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcMP2
  double precision              :: Ecguess
  double precision              :: EpCCD, EcpCCD
  integer                       :: nunit, io, stat
  logical                       :: ex
  double precision,allocatable  :: r_eO(:), r_eV(:)
  double precision,allocatable  :: r_delta_OV(:,:)
  double precision,allocatable  :: r_OOOO(:,:)
  double precision,allocatable  :: r_OOVV(:,:)
  double precision,allocatable  :: r_OVOV(:,:)
  double precision,allocatable  :: r_VVVV(:,:)
  double precision,allocatable  :: r_OVVO(:,:)
  double precision,allocatable  :: residual(:,:)
  double precision,allocatable  :: residual_block(:)


  call check_allowed_t2_guess(t2_guess)
  call check_allowed_z2_guess(z2_guess)
  call check_allowed_t2_update_algorithm(t2_update_algorithm)
  call check_allowed_z2_update_algorithm(z2_update_algorithm)


! Number of valence occupied and virtual orbitals, and number of pairs
  nO = nOa + n_frozen
  npairs = nOa * nV
  write(*,*) 'nOa    = ', nOa
  write(*,*) 'nV     = ', nV
  write(*,*) 'npairs = ', npairs

! Hartree-Fock energy
  ETHF=hf_energy

! Hartree-Fock orbital energies
  allocate(eHF(mo_num))
  eHF(:)=fock_matrix_diag_mo(:)

! provide mo_two_e_integrals_in_map

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|         pCCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Form energy denominator
  allocate( r_eO(nOa), r_eV(nV) )
  r_eO(1:nOa) = eHF(n_frozen+1:nO)
  r_eV(:)     = eHF(nO+1:mo_num)
  allocate( r_delta_OV(nOa,nV) )
  call form_delta_OV(nOa,nV,r_eO,r_eV,r_delta_OV)
  deallocate( r_eO, r_eV )

! double precision :: shift
! shift = 1.0d-3
! call add_shift_to_matrix(r_delta_OV,nOa,nV,shift)

! Create integral batches
  allocate( r_OOVV(nOa,nV), r_OVOV(nOa,nV), &
            r_VVVV(nV,nV),  r_OOOO(nOa,nOa), &
            r_OVVO(nOa,nV) )
  do j=1,nOa
    do i=1,nOa
!     r_OOOO(i,j) = mo_two_e_integrals(i,i,j,j)
      r_OOOO(i,j) = mo_two_e_integrals(n_frozen+i,n_frozen+i,n_frozen+j,n_frozen+j)
    end do
  end do
  do a=1,nV
    do i=1,nOa
!     r_OOVV(i,a) = mo_two_e_integrals(i,i,nO+a,nO+a)
      r_OOVV(i,a) = mo_two_e_integrals(n_frozen+i,n_frozen+i,nO+a,nO+a)
    end do
  end do
  do a=1,nV
    do i=1,nOa
!     r_OVOV(i,a) = mo_two_e_integrals(i,nO+a,i,nO+a)
      r_OVOV(i,a) = mo_two_e_integrals(n_frozen+i,nO+a,n_frozen+i,nO+a)
    end do
  end do
  do a=1,nV
    do i=1,nOa
!     r_OVVO(i,a) = mo_two_e_integrals(i,nO+a,nO+a,i)
      r_OVVO(i,a) = mo_two_e_integrals(n_frozen+i,nO+a,nO+a,n_frozen+i)
    end do
  end do
  do b=1,nV
    do a=1,nV
!     r_VVVV(a,b) = mo_two_e_integrals(nO+a,nO+a,nO+b,nO+b)
      r_VVVV(a,b) = mo_two_e_integrals(nO+a,nO+a,nO+b,nO+b)
    end do
  end do


! allocate(t2(nOa,nV))

! Guess for t-amplitudes

! MP2 guess (even if the MP2 guess is not used for the t-amplitudes, this is done to compute the MP2 energy)
  call amplitude_guess_mp2(t2,r_OOVV,r_delta_OV,nOa,nV)

! MP2 correlation energy
  call compute_Ec(EcMP2,r_OOVV,t2,nOa,nV)

! if( t2_guess=='MP2' ) call amplitude_guess_mp2(t2,r_OOVV,r_delta_OV,nOa,nV)
! In this case, the t-amplitudes were just computed above

  if( t2_guess=='zero' ) t2 = 0.0d0

! TODO: There is some bug in compute_pccd_t2_jacobian_diagonal
  if( t2_guess=='residual_diagonal' ) then 
    stop 'Not yet coded for t2_guess=residual_diagonal'
    t2 = 0.0d0
    double precision, allocatable :: t2_jacobian_diagonal(:,:)
    allocate( t2_jacobian_diagonal(nOa,nV) )
    call compute_pccd_t2_jacobian_diagonal(nOa,nV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,t2_jacobian_diagonal)
    t2 = -r_OOVV / t2_jacobian_diagonal
    deallocate( t2_jacobian_diagonal )
  end if

  if( t2_guess=='residual_full' ) then

    t2 = 0.0d0

    double precision, allocatable :: t2_jacobian(:,:,:,:)
    allocate( t2_jacobian(nOa,nV,nOa,nV) )
    call compute_pccd_t2_jacobian(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,t2_jacobian)

    double precision, allocatable :: t2_jacobian_square(:,:)
    allocate( t2_jacobian_square(npairs,npairs) )
    call transform_4index_to_2index_square(t2_jacobian,t2_jacobian_square,nOa,nV,nOa,nV)
    deallocate( t2_jacobian )

    double precision, allocatable :: t2_jacobian_square_inv(:,:)
    allocate( t2_jacobian_square_inv(npairs,npairs) )
    call get_inverse(t2_jacobian_square,npairs,npairs,t2_jacobian_square_inv,npairs)
    deallocate( t2_jacobian_square )

    allocate( residual_block(npairs) )
    call transform_2index_to_1index_square(r_OOVV,residual_block,nOa,nV)

    double precision, allocatable :: t2_block(:)
    allocate( t2_block(npairs) )
    t2_block = 0.0d0
    do ia=1,npairs
      t2_block(ia) = t2_block(ia) - sum( t2_jacobian_square_inv(ia,:) * residual_block(:) )
    end do
    deallocate( t2_jacobian_square_inv, residual_block )

    call transform_1index_to_2index_square(t2_block,t2,nOa,nV)
    deallocate( t2_block )

  end if

  if( t2_guess=='read' ) then

    character(len=*), parameter :: guess_t_amplitudes_file = 'guess_t_amplitudes.dat'
    inquire(file=guess_t_amplitudes_file,exist=ex)

    if( .not. ex ) then
      write(*,*) 't2_guess == read, however could not find the file ', trim(adjustl(guess_t_amplitudes_file))
      stop
    end if
    open(newunit=nunit, file=guess_t_amplitudes_file, iostat=io, status='unknown')
      write(*,*) npairs
      do b=1,npairs
        read(nunit,*) i, a, t2(i,a)
      end do
    close(nunit)

  end if

! write(*,*) 'Guess t-amplitudes:'
! call write_ij_2d_array(t2,nOa,nV)


! Correlation energy for the guess amplitudes:
  call compute_Ec(Ecguess,r_OOVV,t2,nOa,nV)

  write(*,'(1X,A10,1X,F12.6)') 'E0          ',hf_energy
  write(*,'(1X,A10,1X,F12.6)') 'Ec(MP2)   = ',EcMP2
  write(*,'(1X,A10,1X,F12.6)') 'E (MP2)   = ',EcMP2   + ETHF
  write(*,'(1X,A10,1X,F12.6)') 'Ec(guess) = ',Ecguess
  write(*,'(1X,A10,1X,F12.6)') 'E (guess) = ',Ecguess + ETHF

! Prepare for the update of t-amplitudes:
  if( t2_update_algorithm == 'constant_diagonal' .or. &
      t2_update_algorithm == 'diagonal'          ) then
    allocate( t2_jacobian_diagonal(nOa,nV) )
    call compute_pccd_t2_jacobian_diagonal(nOa,nV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,t2_jacobian_diagonal)
  end if

! Initialization
  allocate(residual(nOa,nV))

  Conv = 1.d0
  nSCF = 0

  if ( diis_t2_amplitudes ) then
    integer :: n_diis
    integer :: max_diis
    double precision :: rcond
    double precision :: max_rcond
    n_diis    = 0
    max_diis  = 10
    max_rcond = 1.0d-15
    double precision, allocatable :: t2_diis(:,:,:)
    allocate( t2_diis(nOa,nV,max_diis) )
    t2_diis = 0.0d0
    double precision, allocatable :: error_diis(:,:,:)
    allocate( error_diis(nOa,nV,max_diis) )
    error_diis = 0.0d0
    double precision, allocatable :: error_in(:,:)
    allocate( error_in(nOa,nV) )
    error_in = 0.0d0
    double precision, allocatable :: error_in_block(:)
  end if

!------------------------------------------------------------------------
! Start of t-amplitudes loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| pCCD calculation                                 |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > cc_thresh .and. nSCF < cc_n_it_max)

!   Increment 
    nSCF = nSCF + 1

!   Compute residual
    call form_rt2_pccd(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,residual)

!   Check convergence 
    Conv = maxval(abs(residual(:,:)))
 
!   Update t-amplitudes based on the selected t2_update_algorithm:

!   Orbital energy difference in the denominator:
    if( t2_update_algorithm == 'orbital_energies' ) then
      if ( diis_t2_amplitudes ) then
        error_in(:,:) = - amplitude_damping * residual(:,:)/r_delta_OV(:,:) 
        t2(:,:) = t2(:,:) + error_in(:,:)
      else
        t2(:,:) = t2(:,:) - amplitude_damping * residual(:,:)/r_delta_OV(:,:)
      end if
    end if

!   Diagonal of the Jacobian (constant terms only):
    if( t2_update_algorithm == 'constant_diagonal' ) then
      if ( diis_t2_amplitudes ) then
        error_in(:,:) = - amplitude_damping * residual(:,:)/t2_jacobian_diagonal(:,:)
        t2(:,:) = t2(:,:) + error_in(:,:)
      else
        t2(:,:) = t2(:,:) - amplitude_damping * residual(:,:)/t2_jacobian_diagonal(:,:)
      end if
    end if

!   Diagonal of the Jacobian (constant + linear terms):
    if( t2_update_algorithm == 'diagonal' ) then
      if ( diis_t2_amplitudes ) then
        error_in(:,:) = - amplitude_damping * residual(:,:)/(t2_jacobian_diagonal(:,:)-2.0d0*r_OOVV(:,:)*t2(:,:))
        t2(:,:) = t2(:,:) + error_in(:,:)
      else
        t2(:,:) = t2(:,:) - amplitude_damping * residual(:,:)/(t2_jacobian_diagonal(:,:)-2.0d0*r_OOVV(:,:)*t2(:,:))
      end if
    end if

!   Full Jacobian:
    if( t2_update_algorithm == 'full_Newton_Raphson' ) then

      allocate( t2_jacobian(nOa,nV,nOa,nV) )
      call compute_pccd_t2_jacobian(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,t2_jacobian)

      allocate( t2_jacobian_square(npairs,npairs) )
      call transform_4index_to_2index_square(t2_jacobian,t2_jacobian_square,nOa,nV,nOa,nV)
      deallocate( t2_jacobian )

      allocate( t2_jacobian_square_inv(npairs,npairs) )
      call get_inverse(t2_jacobian_square,npairs,npairs,t2_jacobian_square_inv,npairs)
      deallocate( t2_jacobian_square )

      allocate( t2_block(npairs) )
      call transform_2index_to_1index_square(t2,t2_block,nOa,nV)

      allocate( residual_block(npairs) )
      call transform_2index_to_1index_square(residual,residual_block,nOa,nV)

      if ( diis_t2_amplitudes ) then
        allocate( error_in_block(npairs) )
        do ia=1,npairs
          error_in_block(ia) = - amplitude_damping * sum( t2_jacobian_square_inv(ia,:) * residual_block(:) )
        end do
        deallocate( t2_jacobian_square_inv, residual_block )
        do ia=1,npairs
          t2_block(ia) = t2_block(ia) + error_in_block(ia)
        end do
        call transform_1index_to_2index_square(error_in_block,error_in,nOa,nV)
        deallocate( error_in_block )
      else
        do ia=1,npairs
          t2_block(ia) = t2_block(ia) - amplitude_damping * sum( t2_jacobian_square_inv(ia,:) * residual_block(:) )
        end do
        deallocate( t2_jacobian_square_inv, residual_block )
      end if
      call transform_1index_to_2index_square(t2_block,t2,nOa,nV)
      deallocate( t2_block )

    end if

    if ( diis_t2_amplitudes ) then
      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,npairs,npairs,n_diis,error_diis,t2_diis,error_in,t2)
    end if


!   Compute correlation energy
    call compute_Ec(EcpCCD,r_OOVV,t2,nOa,nV)

    EpCCD = ETHF + EcpCCD

!   Dump results
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F12.8,1X,A1,1X)') &
      '|',nSCF,'|',EpCCD,'|',EcpCCD,'|',Conv,'|'

  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of t-amplitudes loop
!------------------------------------------------------------------------
! Did it actually converge?
  if(nSCF == cc_n_it_max) then
    call write_convergence_failed
!   Write a file called "failed", which can be a flag for a calling script
    call system('touch failed')
    stop
  else
    call save_energy(EpCCD)
  endif
  deallocate( residual )



! allocate(z2(nOa,nV))

! Guess for z-amplitudes

  if( z2_guess.eq.'t-amplitudes' ) z2 = t2

  if( z2_guess.eq.'MP2'  ) call amplitude_guess_mp2(z2,r_OOVV,r_delta_OV,nOa,nV)

  if( z2_guess.eq.'zero' ) z2 = 0.0d0

  if( z2_guess=='residual_diagonal' ) then
    double precision, allocatable :: z2_jacobian_diagonal(:,:)
    allocate( z2_jacobian_diagonal(nOa,nV) )
    call compute_pccd_z2_jacobian_diagonal(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2_jacobian_diagonal)
    z2 = -r_OOVV / z2_jacobian_diagonal
    deallocate( z2_jacobian_diagonal )
  end if

  if( z2_guess=='residual_full' ) then

    double precision, allocatable :: z2_jacobian(:,:,:,:)
    allocate( z2_jacobian(nOa,nV,nOa,nV) )
    call compute_pccd_z2_jacobian(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2_jacobian)

    double precision, allocatable :: z2_jacobian_square(:,:)
    allocate( z2_jacobian_square(npairs,npairs) )
    call transform_4index_to_2index_square(z2_jacobian,z2_jacobian_square,nOa,nV,nOa,nV)
    deallocate( z2_jacobian )

    double precision, allocatable :: z2_jacobian_square_inv(:,:)
    allocate( z2_jacobian_square_inv(npairs,npairs) )
    call get_inverse(z2_jacobian_square,npairs,npairs,z2_jacobian_square_inv,npairs)
    deallocate( z2_jacobian_square )

    double precision, allocatable :: z2_block(:)
    allocate( z2_block(npairs) )
    z2_block = 0.0d0

    allocate( residual_block(npairs) )
    call transform_2index_to_1index_square(r_OOVV,residual_block,nOa,nV)

    do ia=1,npairs
      z2_block(ia) = z2_block(ia) - sum( z2_jacobian_square_inv(ia,:) * residual_block(:) )
    end do
    deallocate( z2_jacobian_square_inv, residual_block )

    call transform_1index_to_2index_square(z2_block,z2,nOa,nV)
    deallocate( z2_block )

  end if

  if( z2_guess=='read' ) then

    character(len=*), parameter :: guess_z_amplitudes_file = 'guess_z_amplitudes.dat'
    inquire(file=guess_z_amplitudes_file,exist=ex)

    if( .not. ex ) then
      write(*,*) 'z2_guess == read, however could not find the file ', trim(adjustl(guess_z_amplitudes_file))
      stop
    end if
    open(newunit=nunit, file=guess_z_amplitudes_file, iostat=io, status='old')
      do b=1,npairs
        read(nunit,*) i, a, z2(i,a)
      end do
    close(nunit)

  end if


! write(*,*) 'Guess z-amplitudes:'
! call write_ij_2d_array(z2,nOa,nV)

!------------------------------------------------------------------------
! Start of z-amplitudes loop
!------------------------------------------------------------------------

! Prepare for the update of z-amplitudes:
  if( z2_update_algorithm == 'constant_diagonal' .or. &
      z2_update_algorithm == 'diagonal'          ) then
    allocate( z2_jacobian_diagonal(nOa,nV) )
    call compute_pccd_z2_jacobian_diagonal(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2_jacobian_diagonal)
  end if

  allocate(residual(nOa,nV))

  Conv = 1.d0
  nSCF = 0

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| z-amplitudes                                     |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > cc_thresh .and. nSCF < cc_n_it_max)

!   Increment 
    nSCF = nSCF + 1

!   Compute residual
    call form_rz2_pccd(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2,residual)

!   Check convergence 
    Conv = maxval(abs(residual(:,:)))

!   Update z-amplitudes based on the selected z2_update_algorithm:

!   Orbital energy difference in the denominator:
    if( z2_update_algorithm == 'orbital_energies' ) then
      z2(:,:) = z2(:,:) - amplitude_damping * residual(:,:)/r_delta_OV(:,:)
    end if

!   Diagonal of the Jacobian (constant terms only):
    if( z2_update_algorithm == 'constant_diagonal' ) then
      z2(:,:) = z2(:,:) - amplitude_damping * residual(:,:)/z2_jacobian_diagonal(:,:)
    end if

!   Diagonal of the Jacobian (constant + linear terms):
    if( z2_update_algorithm == 'diagonal' ) then
      z2(:,:) = z2(:,:) - amplitude_damping * residual(:,:)/(z2_jacobian_diagonal(:,:)-2.0d0*r_OOVV(:,:)*z2(:,:))
    end if

!   Full Jacobian:
    if( z2_update_algorithm == 'full_Newton_Raphson' ) then

      allocate( z2_jacobian(nOa,nV,nOa,nV) )
      call compute_pccd_z2_jacobian(nOa,nV,r_OOVV,r_OVOV,r_OVVO,r_VVVV,r_OOOO,r_delta_OV,t2,z2_jacobian)

      allocate( z2_jacobian_square(npairs,npairs) )
      call transform_4index_to_2index_square(z2_jacobian,z2_jacobian_square,nOa,nV,nOa,nV)
      deallocate( z2_jacobian )

      allocate( z2_jacobian_square_inv(npairs,npairs) )
      call get_inverse(z2_jacobian_square,npairs,npairs,z2_jacobian_square_inv,npairs)
      deallocate( z2_jacobian_square )

      allocate( z2_block(npairs) )
      call transform_2index_to_1index_square(z2,z2_block,nOa,nV)

      allocate( residual_block(npairs) )
      call transform_2index_to_1index_square(residual,residual_block,nOa,nV)

      do ia=1,npairs
        z2_block(ia) = z2_block(ia) - sum( z2_jacobian_square_inv(ia,:) * residual_block(:) )
      end do
      deallocate( z2_jacobian_square_inv, residual_block )

      call transform_1index_to_2index_square(z2_block,z2,nOa,nV)
      deallocate( z2_block )
 
    end if


!   Dump results
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F12.8,1X,A1,1X)') &
      '|',nSCF,'|',EpCCD,'|',EcpCCD,'|',Conv,'|'

  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of z-amplitudes loop
!------------------------------------------------------------------------
! Did it actually converge?
  if(nSCF == cc_n_it_max) then
    call write_convergence_failed
!   Write a file called "failed", which can be a flag for a calling script
    call system('touch failed')
    stop
  endif
  deallocate( residual )


  deallocate( r_delta_OV )
  deallocate( r_OOVV, r_OVOV, r_VVVV, r_OOOO, r_OVVO )

  if( write_final_t_amplitudes ) then
    character(len=*), parameter :: final_t_amplitudes_file = 'final_t_amplitudes.dat'
    call write_final_amplitudes(t2,final_t_amplitudes_file,nOa,nV)
  end if
  if( write_final_z_amplitudes ) then
    character(len=*), parameter :: final_z_amplitudes_file = 'final_z_amplitudes.dat'
    call write_final_amplitudes(z2,final_z_amplitudes_file,nOa,nV)
  end if

! write(*,*) 'Fock diagonal:'
! call write_i_1d_array(eHF,nbas)

  deallocate( eHF )

  write(*,*) 'EHF     ', hf_energy
  write(*,*) 'EcpCCD  ', EcpCCD
  write(*,*) 'EpCCD   ', EpCCD

! TODO:
! call system('if [ -f E_1last.dat ]; then mv E_1last.dat E_2last.dat;  fi')
! character(len=*), parameter :: E_try_file = 'E_try.dat'
! call write_double_to_file(EpCCD,E_try_file)

end subroutine run_pccd


subroutine check_allowed_t2_guess(char)
  implicit none
  BEGIN_DOC
! Stops if a wrong entry for t2_guess is given
  END_DOC
  character(len=*), intent(in) :: char
  if( char=='residual_full' ) return
  if( char=='residual_diagonal' ) return
  if( char=='MP2' ) return
  if( char=='read' ) return
  if( char=='zero' ) return
  write(*,*) 'Variable for t2_guess not allowed: ', char
  stop
end subroutine check_allowed_t2_guess


subroutine check_allowed_z2_guess(char)
  implicit none
  BEGIN_DOC
! Stops if a wrong entry for z2_guess is given
  END_DOC
  character(len=*), intent(in) :: char
  if( char=='residual_full' ) return
  if( char=='residual_diagonal' ) return
  if( char=='MP2' ) return
  if( char=='read' ) return
  if( char=='zero' ) return
  if( char=='t-amplitudes' ) return
  write(*,*) 'Variable for z2_guess not allowed: ', char
  stop
end subroutine check_allowed_z2_guess


subroutine check_allowed_t2_update_algorithm(char)
  implicit none
  BEGIN_DOC
! Stops if a wrong entry for t2_update_algorithm is given
  END_DOC
  character(len=*), intent(in) :: char
  if( char=='full_Newton_Raphson' ) return
  if( char=='diagonal' ) return
  if( char=='constant_diagonal' ) return
  if( char=='orbital_energies' ) return
  write(*,*) 'Variable for t2_update_algorithm not allowed: ', char
  stop
end subroutine check_allowed_t2_update_algorithm


subroutine check_allowed_z2_update_algorithm(char)
  implicit none
  BEGIN_DOC
! Stops if a wrong entry for z2_update_algorithm is given
  END_DOC
  character(len=*), intent(in) :: char
  if( char=='full_Newton_Raphson' ) return
  if( char=='diagonal' ) return
  if( char=='constant_diagonal' ) return
  if( char=='orbital_energies' ) return
  write(*,*) 'Variable for z2_update_algorithm not allowed: ', char
  stop
end subroutine check_allowed_z2_update_algorithm


subroutine amplitude_guess_mp2(t2,r_OOVV,r_delta_OV,nOa,nV)
  implicit none
  BEGIN_DOC
! Returns the pCCD guess amplitudes based on the MP2 amplitudes
  END_DOC
  integer,          intent(in)  :: nOa, nV
  double precision, intent(out) :: t2(nOa,nV)
  double precision, intent(in)  :: r_OOVV(nOa,nV)
  double precision, intent(in)  :: r_delta_OV(nOa,nV)
  t2(:,:) = -r_OOVV(:,:)/r_delta_OV(:,:)
end subroutine amplitude_guess_mp2


subroutine compute_Ec(Ec,r_OOVV,t2,nOa,nV)
  implicit none
  BEGIN_DOC
! Computes the correlation energy of pCCD
  END_DOC
  double precision, intent(out) :: Ec
  integer,          intent(in)  :: nOa, nV
  double precision, intent(in)  :: r_OOVV(nOa,nV)
  double precision, intent(in)  :: t2(nOa,nV)
  integer                       :: a
  Ec = 0.0d0
  do a=1,nV
    Ec = Ec + sum( r_OOVV(:,a) * t2(:,a) )
  end do
end subroutine compute_Ec


subroutine write_final_amplitudes(t2,filename,nOa,nV)
  implicit none
  BEGIN_DOC
! Write the converged pCCD amplitudes to a file
  END_DOC
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: nOa, nV
  double precision, intent(in) :: t2(nOa,nV)
  integer :: i, a
  integer :: nunit, io
  open(newunit=nunit, file=filename, iostat=io, status='unknown')
    do i=1,nOa
      do a=1,nV
        write(nunit,*) i, a, t2(i,a)
      end do
    end do
  close(nunit)
end subroutine write_final_amplitudes


subroutine write_convergence_failed
  implicit none
  BEGIN_DOC
! Convergence failed message
  END_DOC
  write(*,*)
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'                 Convergence failed                 '
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)
end subroutine write_convergence_failed


subroutine add_zero_amplitudes_for_frozen_part(t2_in,t2_out,nOa,nO,nV)
implicit none
  BEGIN_DOC
! Assign zero to the pCCD amplitudes corresponding to excitations from the core
  END_DOC
integer,          intent(in)  :: nOa, nO, nV
double precision, intent(in)  :: t2_in(nOa,nV)
double precision, intent(out) :: t2_out(nO,nV)
integer :: i, a, n_diff
n_diff = nO - nOa
do a=1,nV
  do i=1,n_diff
    t2_out(i,a) = 0.d0
  end do
  do i=1,nOa
    t2_out(n_diff+i,a) = t2_in(i,a)
  end do
end do
end subroutine add_zero_amplitudes_for_frozen_part


subroutine save_energy(E)
  implicit none
  BEGIN_DOC
! Saves the energy in |EZFIO|.
  END_DOC
  double precision, intent(in) :: E
  call ezfio_set_pccd_energy(E)
end subroutine save_energy

