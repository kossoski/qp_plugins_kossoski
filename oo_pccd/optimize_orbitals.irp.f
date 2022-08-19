subroutine optimize_orbitals(r_one_e_dm_mo,r_two_e_dm_mo,nT,is_converged)

  implicit none
  BEGIN_DOC
  ! Performs one iteration in orbital optimized pair coupled cluster doubles (oo-pCCD), updating the orbitals in EZFIO in the end
  END_DOC

  integer,         intent(in)   :: nT
  double precision,intent(in)   :: r_one_e_dm_mo(nT,nT)
  double precision,intent(in)   :: r_two_e_dm_mo(nT,nT,nT,nT)
  logical, intent(inout)        :: is_converged 

! Local variables

  integer                       :: i,j,p
  double precision,allocatable  :: r_orbrot_g(:,:)
  double precision,allocatable  :: r_orbrot_g_vector(:)
  double precision,allocatable  :: r_orbrot_h(:,:,:,:)
  double precision,allocatable  :: r_orbrot_h_square_inv(:,:)
  double precision,allocatable  :: r_orbrot_h_square(:,:)
  double precision,allocatable  :: kappa_square(:,:)
  double precision,allocatable  :: kappa_vector(:)
  double precision,allocatable  :: r_unitary_orbrot(:,:)
  logical                       :: ex

! Thresholds on convergence of optimized orbitals
  double precision :: max_grad, mean_grad, max_kappa, mean_kappa
  double precision :: max_grad_thresh, mean_grad_thresh, max_kappa_thresh, mean_kappa_thresh
  max_grad_thresh   = 1.d-5
  mean_grad_thresh  = 1.d-6
  max_kappa_thresh  = 1.d-3
  mean_kappa_thresh = 1.d-4

  call check_allowed_Hessian_update(Hessian_update)
  call check_allowed_LA_solver_orb_opt(LA_solver_orb_opt)


! Form orbital rotation gradient

  allocate( r_orbrot_g(nT,nT) )
! call compute_r_orbrot_g(r_orbrot_g,r_one_e_dm_mo,r_two_e_dm_mo,nT)
  call compute_r_orbrot_g_pccd(r_orbrot_g,r_one_e_dm_mo,r_two_e_dm_mo,nT)

! write(*,*) 'Gradients (2d):'
! call write_ij_2d_array(r_orbrot_g,nT,nT)

  integer :: nT_tri
  nT_tri = nT * (nT-1) / 2
  allocate( r_orbrot_g_vector(nT_tri) )
  call transform_2index_to_1index_lowtri(r_orbrot_g,r_orbrot_g_vector,nT)
  deallocate( r_orbrot_g )

  write(*,*) 
! write(*,*) 'Gradients (1d):'
! call write_i_1d_array(r_orbrot_g_vector,nT_tri)


! Write maximum and mean absolute gradient
  write(*,*) 
  max_grad  = maxval( abs(r_orbrot_g_vector(:)) )
  mean_grad = sum( abs(r_orbrot_g_vector(:)) ) / size(r_orbrot_g_vector)
  write(*,*) 'Maximum absolute gradient:     ', max_grad
  write(*,*) 'Mean absolute gradient:        ', mean_grad

! Start Hessian_update case

! Full Newton-Raphson:
  if( Hessian_update == 'full_Newton_Raphson' ) then

    allocate( r_orbrot_h(nT,nT,nT,nT) )
!   call compute_r_orbrot_h(r_orbrot_h,r_one_e_dm_mo,r_two_e_dm_mo,nT)
    call compute_r_orbrot_h_pccd(r_orbrot_h,r_one_e_dm_mo,r_two_e_dm_mo,nT)
    write(*,*) 'Computed orbital rotation Hessian'

    allocate( r_orbrot_h_square(nT_tri,nT_tri) )
    call transform_4index_to_2index_lowtri(r_orbrot_h,r_orbrot_h_square,nT,nT)
    write(*,*) 'Transformed orbital rotation Hessian'
    deallocate( r_orbrot_h )

  end if

! Diagonal Newton-Raphson:
  if( Hessian_update == 'diagonal' ) then

    double precision, allocatable :: r_orbrot_h_diagonal(:,:)
    allocate( r_orbrot_h_diagonal(nT,nT) )
    call compute_r_orbrot_h_diagonal(r_orbrot_h_diagonal,r_one_e_dm_mo,r_two_e_dm_mo,nT)
!   call compute_r_orbrot_h_diagonal_pccd(r_orbrot_h_diagonal,r_one_e_dm_mo,r_two_e_dm_mo,nT)
    write(*,*) 'Computed diagonal of the orbital rotation Hessian'

    double precision, allocatable :: r_orbrot_h_diagonal_square(:)
    allocate( r_orbrot_h_diagonal_square(nT_tri) )
    call transform_2index_to_1index_lowtri(r_orbrot_h_diagonal,r_orbrot_h_diagonal_square,nT)
    deallocate( r_orbrot_h_diagonal )

    allocate( r_orbrot_h_square(nT_tri,nT_tri) )
    r_orbrot_h_square = 0.0d0
    do i=1,nT_tri
      r_orbrot_h_square(i,i) = r_orbrot_h_diagonal_square(i)
    end do
    write(*,*) 'Transformed orbital rotation Hessian'
    deallocate( r_orbrot_h_diagonal_square )

  end if


! Diagonal with 1's:
  if( Hessian_update == 'identity_diagonal' ) then
    allocate( r_orbrot_h_square(nT_tri,nT_tri) )
    call set_identity_matrix(r_orbrot_h_square,nT_tri)
  end if

! End Hessian_update case

! Add diagonal mu_damping:
! call add_shift_to_diagonal(r_orbrot_h_square,nT_tri,mu_damping)


! write(*,*) 'r_orbrot_h_square:'
! call write_ij_2d_array(r_orbrot_h_square,nT_tri,nT_tri)

! Check whether r_orbrot_h_square is symmetric:
  call check_matrix_symmetry(r_orbrot_h_square,nT_tri)


! Diagonalize the Hessian matrix and write its eigenvalues
  double precision, allocatable :: hessian_eigvalues(:)
  double precision, allocatable :: hessian_eigvectors(:,:)
  allocate( hessian_eigvalues(nT_tri) )
  allocate( hessian_eigvectors(nT_tri,nT_tri) )
  call lapack_diag(hessian_eigvalues,hessian_eigvectors,r_orbrot_h_square,nT_tri,nT_tri)
! call lapack_diagd(hessian_eigvalues,hessian_eigvectors,r_orbrot_h_square,nT_tri,nT_tri)
  write(*,*) 'Diagonalized the orbital rotation Hessian'
  write(*,*) 
  write(*,*) 'Hessian eigenvalues:'
  call write_i_1d_array(hessian_eigvalues,nT_tri)
! write(*,*) 'Hessian eigenvectors:'
! call write_ij_2d_array(hessian_eigvectors,nT_tri,nT_tri)

! Avoid problems with very small eigenvalues:
  do i=1,nT_tri
    if( abs(hessian_eigvalues(i)).lt.1.0d-10 ) hessian_eigvalues(i) = 1.0d20
  end do

  allocate( kappa_vector(nT_tri) )
  kappa_vector = 0.0d0

  double precision, allocatable :: aug_r_orbrot_h_square(:,:)
  double precision, allocatable :: aug_hessian_eigvalues(:)
  double precision, allocatable :: aug_hessian_eigvectors(:,:)

  integer :: saddle_order1
  saddle_order1 = saddle_order + 1

! Start LA_solver_orb_opt case

  if( LA_solver_orb_opt=='diagonalize_augmented' ) then

    integer :: nT_tri1
    nT_tri1 = nT_tri + 1
! Build the augmented Hessian matrix
    allocate( aug_r_orbrot_h_square(nT_tri1,nT_tri1) )
    call build_augmented_matrix(r_orbrot_h_square,r_orbrot_g_vector,nT_tri,aug_r_orbrot_h_square)

! Diagonalize the augmented Hessian matrix
    allocate( aug_hessian_eigvalues(nT_tri1) )
    allocate( aug_hessian_eigvectors(nT_tri1,nT_tri1) )
    call lapack_diag(aug_hessian_eigvalues,aug_hessian_eigvectors,aug_r_orbrot_h_square,nT_tri1,nT_tri1)
    deallocate( aug_hessian_eigvectors )
    deallocate( aug_r_orbrot_h_square )
    write(*,*)
    write(*,*) 'Diagonalized augmented Hessian'

    write(*,*)
    write(*,*) 'Augmented Hessian eigenvalues:'
    call write_i_1d_array(aug_hessian_eigvalues,nT_tri1)

! Let us use kappa_vector as a buffer:
    do i=1,nT_tri
      kappa_vector(i) = sum( hessian_eigvectors(:,i) * r_orbrot_g_vector(:) )
    end do

!   write(*,*) 'Projection of gradients (2d):'
!   call write_i_1d_array(kappa_vector,nT_tri)

! Let us use r_orbrot_g_vector as a buffer:
    do i=1,nT_tri
      r_orbrot_g_vector(i) = kappa_vector(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(saddle_order1) * lambda_hessian )
    end do
! Test if the Hessian has the desired structure:
    logical :: correct_saddle_order
    correct_saddle_order = .true.
    double precision :: saddle_order_eps
    saddle_order_eps = 1.0d-10
    do i=1,saddle_order
      if( hessian_eigvalues(i).gt.(-saddle_order_eps) ) then
        correct_saddle_order = .false.
        exit
      end if
    end do
    if( hessian_eigvalues(saddle_order1).lt.saddle_order_eps ) then
      correct_saddle_order = .false.
    end if

    deallocate( hessian_eigvalues )

! And finally the actual kappa_vector:
    do i=1,nT_tri
      kappa_vector(i) = - sum( hessian_eigvectors(i,:) * r_orbrot_g_vector(:) )
    end do
    deallocate( hessian_eigvectors )

  end if



  if( LA_solver_orb_opt=='diagonalize_augmented_p' ) then

    double precision, allocatable :: kappa_vector_aux(:)
    allocate( kappa_vector_aux(nT_tri) )
! Let us use kappa_vector as a buffer:
    do i=1,nT_tri
      kappa_vector_aux(i) = sum( hessian_eigvectors(:,i) * r_orbrot_g_vector(:) )
    end do

! Negative eigevalues:

!   write(*,*) 'Projection of (negative) gradients (2d):'
!   call write_i_1d_array(kappa_vector_aux(1:saddle_order),saddle_order)

    deallocate( r_orbrot_h_square )
    allocate( r_orbrot_h_square(saddle_order,saddle_order) )
    r_orbrot_h_square = 0.0d0
    do i=1,saddle_order
      r_orbrot_h_square(i,i) = hessian_eigvalues(i)
    end do

! Build the augmented (negative) Hessian matrix
    allocate( aug_r_orbrot_h_square(saddle_order1,saddle_order1) )
    call build_augmented_matrix(r_orbrot_h_square,kappa_vector_aux(1:saddle_order),saddle_order,aug_r_orbrot_h_square)

! Diagonalize the augmented Hessian matrix
    allocate( aug_hessian_eigvalues(saddle_order1) )
    allocate( aug_hessian_eigvectors(saddle_order1,saddle_order1) )
    call lapack_diag(aug_hessian_eigvalues,aug_hessian_eigvectors,aug_r_orbrot_h_square,saddle_order1,saddle_order1)
    deallocate( aug_hessian_eigvectors )
    deallocate( aug_r_orbrot_h_square )
    write(*,*)
    write(*,*) 'Diagonalized augmented (negative) Hessian'
    write(*,*)
    write(*,*) 'Augmented (negative) Hessian eigenvalues:'
    call write_i_1d_array(aug_hessian_eigvalues,saddle_order1)

! Let us use r_orbrot_g_vector as a buffer:
    do i=1,saddle_order
      r_orbrot_g_vector(i) = kappa_vector_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(saddle_order1) )
!     r_orbrot_g_vector(i) = kappa_vector_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(saddle_order1) * lambda_hessian )
    end do
! Test if the Hessian has the desired structure:
!   logical :: correct_saddle_order
    correct_saddle_order = .true.
!   double precision :: saddle_order_eps
    saddle_order_eps = 1.0d-10
    do i=1,saddle_order
      if( hessian_eigvalues(i).gt.(-saddle_order_eps) ) then
        correct_saddle_order = .false.
        exit
      end if
    end do
    if( hessian_eigvalues(saddle_order1).lt.saddle_order_eps ) then
      correct_saddle_order = .false.
    end if

    deallocate( aug_hessian_eigvalues )

!   write(*,*) 'Projection of (negative) g (2d):'
!   call write_i_1d_array(r_orbrot_g_vector(1:saddle_order),saddle_order)

! And finally the actual kappa_vector:
    do i=1,nT_tri
      kappa_vector(i) = - sum( hessian_eigvectors(i,1:saddle_order) * r_orbrot_g_vector(1:saddle_order) )
    end do

!   write(*,*) 'Projection of (negative) kappa (2d):'
!   call write_i_1d_array(kappa_vector(1:saddle_order),saddle_order)

! Positive eigevalues:
   
    integer :: comp_saddle_order, comp_saddle_order1
    comp_saddle_order  = nT_tri - saddle_order
    comp_saddle_order1 = comp_saddle_order + 1

! Let us use kappa_vector as a buffer:
!   do i=saddle_order1,nT_tri
!     kappa_vector(i) = sum( hessian_eigvectors(:,i) * r_orbrot_g_vector(:) )
!   end do

!   write(*,*) 'Projection of (positive) gradients (2d):'
!   call write_i_1d_array(kappa_vector_aux(saddle_order1:nT_tri),comp_saddle_order)

    deallocate( r_orbrot_h_square )
    allocate( r_orbrot_h_square(comp_saddle_order,comp_saddle_order) )
    r_orbrot_h_square = 0.0d0
    do i=1,comp_saddle_order
      r_orbrot_h_square(i,i) = hessian_eigvalues(saddle_order+i)
    end do

! Build the augmented (positive) Hessian matrix
    allocate( aug_r_orbrot_h_square(comp_saddle_order1,comp_saddle_order1) )
    call build_augmented_matrix(r_orbrot_h_square,kappa_vector_aux(saddle_order1:nT_tri),comp_saddle_order,aug_r_orbrot_h_square)

! Diagonalize the augmented Hessian matrix
    allocate( aug_hessian_eigvalues(comp_saddle_order1) )
    allocate( aug_hessian_eigvectors(comp_saddle_order1,comp_saddle_order1) )

    call lapack_diag(aug_hessian_eigvalues,aug_hessian_eigvectors,aug_r_orbrot_h_square,comp_saddle_order1,comp_saddle_order1)
    deallocate( aug_hessian_eigvectors )
    deallocate( aug_r_orbrot_h_square )
    write(*,*)
    write(*,*) 'Diagonalized augmented (positive) Hessian'
    write(*,*)
    write(*,*) 'Augmented (positive) Hessian eigenvalues:'
    call write_i_1d_array(aug_hessian_eigvalues(1),1)
!   call write_i_1d_array(aug_hessian_eigvalues,comp_saddle_order1)

! Let us use r_orbrot_g_vector as a buffer:
    do i=saddle_order1,nT_tri
!     r_orbrot_g_vector(i) = kappa_vector_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(1) )
      r_orbrot_g_vector(i) = kappa_vector_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(1) * lambda_hessian )
    end do
    deallocate( aug_hessian_eigvalues )
    deallocate( kappa_vector_aux )

! And finally the actual kappa_vector:
    do i=1,nT_tri
      kappa_vector(i) = kappa_vector(i) - sum( hessian_eigvectors(i,saddle_order1:nT_tri) * r_orbrot_g_vector(saddle_order1:nT_tri) )
    end do

!   write(*,*) 'Projection of (positive) kappa (2d):'
!   call write_i_1d_array(kappa_vector(saddle_order1:nT_tri),comp_saddle_order)

    deallocate( hessian_eigvalues )
    deallocate( hessian_eigvectors )

  end if


  if( LA_solver_orb_opt=='diagonalize' ) then

! From the diagonalization:
! H-1 g = (U Z U^T)-1 g = ( U^T^-1 Z^-1 U^-1 ) g = ( U Z^-1 U^T ) g 

! Let us use kappa_vector as a buffer:
    do i=1,nT_tri
      kappa_vector(i) = sum( hessian_eigvectors(:,i) * r_orbrot_g_vector(:) )
    end do

    write(*,*) 'Projection of gradients (2d):'
    call write_i_1d_array(kappa_vector,nT_tri)
 
    do i=1,nT_tri
      if( hessian_eigvalues(i) .gt. 0.0d0 ) then
        hessian_eigvalues(i) = hessian_eigvalues(i) + mu_damping
      else
        hessian_eigvalues(i) = hessian_eigvalues(i) - mu_damping
      end if
    end do

    if( move_negative_Hessian ) then
      do i=1,min_negative_direction-1
        r_orbrot_g_vector(i) = kappa_vector(i) / hessian_eigvalues(i)
      end do
      do i=min_negative_direction,max_negative_direction
        if( hessian_eigvalues(i) .le. 0.0d0 ) then
!         r_orbrot_g_vector(i) = sign( step_negative_direction, kappa_vector(i) )
          r_orbrot_g_vector(i) = sign( - step_negative_direction / hessian_eigvalues(i), kappa_vector(i) )
        else
          r_orbrot_g_vector(i) = kappa_vector(i) / hessian_eigvalues(i)
        end if
      end do
      do i=max_negative_direction+1,nT_tri
        r_orbrot_g_vector(i) = kappa_vector(i) / hessian_eigvalues(i)
      end do
    else
! Let us use r_orbrot_g_vector as a buffer:
      do i=1,nT_tri
        r_orbrot_g_vector(i) = kappa_vector(i) / hessian_eigvalues(i)
      end do
    end if
    deallocate( hessian_eigvalues )

! And finally the actual kappa_vector:
    do i=1,nT_tri
      kappa_vector(i) = - sum( hessian_eigvectors(i,:) * r_orbrot_g_vector(:) )
    end do
    deallocate( hessian_eigvectors )

  end if


  if( LA_solver_orb_opt=='inverse' ) then

    allocate( r_orbrot_h_square_inv(nT_tri,nT_tri) )
    call get_inverse(r_orbrot_h_square,nT_tri,nT_tri,r_orbrot_h_square_inv,nT_tri)

    do p=1,nT_tri
      kappa_vector(p) = - sum( r_orbrot_h_square_inv(p,:) * r_orbrot_g_vector(:) )
    end do
    deallocate( r_orbrot_h_square_inv )

  end if


  if( LA_solver_orb_opt=='linear_system' ) then
    stop 'Not yet coded for LA_solver_orb_opt=linear_system'
! TODO: Does not work if I try to obtain kappa by solving the linear system with this subroutine
!   call lapack_solve_linear_dsy(nT_tri,r_orbrot_h_square,-r_orbrot_g_vector,kappa_vector)
  end if

! End LA_solver_orb_opt case


! Write maximum and mean kappa_vector
  write(*,*) 
  write(*,*) 'Updated k_vector'
  write(*,*) 'Norm of kappa_vector:          ', sqrt( sum( kappa_vector(:)**2 ) )
  write(*,*) 
  max_kappa = maxval( abs(kappa_vector(:)) )
  mean_kappa = sum( abs(kappa_vector(:)) ) / size(kappa_vector)
  write(*,*) 'Maximum absolute kappa_vector: ', max_kappa
  write(*,*) 'Mean absolute kappa_vector:    ', mean_kappa
  write(*,*) 
! write(*,*) 'kappa_vector:'
! call write_i_1d_array(kappa_vector,nT_tri)

  deallocate( r_orbrot_h_square )
  deallocate( r_orbrot_g_vector )


! Test for convergence
  if( max_grad .le.max_grad_thresh  .and. mean_grad .le.mean_grad_thresh .and. &
      max_kappa.le.max_kappa_thresh .and. mean_kappa.le.mean_kappa_thresh ) then
    is_converged = .true.
  end if

  if( LA_solver_orb_opt=='diagonalize_augmented_p' .or. LA_solver_orb_opt=='diagonalize_augmented_p' ) then
     if( .not.correct_saddle_order ) is_converged = .false.
  end if


  allocate( kappa_square(nT,nT) )
  call form_A_square(kappa_vector,nT_tri,kappa_square,nT)
  deallocate( kappa_vector )

! write(*,*) 'kappa_square:'
! call write_ij_2d_array(kappa_square,nT,nT)

! Form the unitary orbital rotation matrix from the antihermitian kappa matrix: U = exp( kappa )
  allocate( r_unitary_orbrot(nT,nT) )
  call lapack_exp_antisymm_matrix(nT,kappa_square,r_unitary_orbrot)
  write(*,*) 'Formed the rotation matrix U = exp( kappa )'
  deallocate( kappa_square )

! Apply the orbital rotation U to the molecular orbitals: d' = d U
  call update_molecular_orbitals(mo_coef,ao_num,nT,r_unitary_orbrot)
  deallocate( r_unitary_orbrot )

! Make sure the orbitals are orthonormal after the rotation:
  call check_mos_orthonormality

! After one (quasi-)Newton-Raphson step, save molecular orbitals to the EZFIO directory:
  call save_mos
  write(*,*) 'Updated the molecular orbitals'

end subroutine optimize_orbitals


subroutine check_allowed_Hessian_update(char)
  implicit none
  BEGIN_DOC
! Stops if a wrong entry for Hessian_update is given
  END_DOC
  character(len=*), intent(in) :: char
  if( char=='full_Newton_Raphson' ) return
  if( char=='diagonal' ) return
! if( char=='SR1' ) return
  if( char=='identity_diagonal' ) return
  write(*,*) 'Variable for Hessian_update not allowed: ', char
  stop
end subroutine check_allowed_Hessian_update


subroutine check_allowed_LA_solver_orb_opt(char)
  implicit none
  BEGIN_DOC
! Stops if a wrong entry for LA_solver_orb_opt is given
  END_DOC
  character(len=*), intent(in) :: char
  if( char=='diagonalize' ) return
  if( char=='diagonalize_augmented' ) return
  if( char=='diagonalize_augmented' ) return
  if( char=='diagonalize_augmented_p' ) return
  if( char=='inverse' ) return
  if( char=='linear_system' ) return
  write(*,*) 'Variable for LA_solver_orb_opt not allowed: ', char
  stop
end subroutine check_allowed_LA_solver_orb_opt

