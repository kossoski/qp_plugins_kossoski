subroutine optimize_orbitals_ci

  implicit none
  BEGIN_DOC
  ! Optimize the orbtitals of a general CI wave function
  END_DOC

! Local variables
  integer :: iteration = 1
  logical :: is_converged

  integer                       :: i,j,p
  double precision,allocatable  :: r_orbrot_g(:,:)
  double precision,allocatable  :: r_orbrot_g_vector(:)
  double precision,allocatable  :: r_orbrot_h(:,:,:,:)
  double precision,allocatable  :: r_orbrot_h_square_inv(:,:)
  double precision,allocatable  :: r_orbrot_h_square(:,:)
  double precision,allocatable  :: kappa_square(:,:)
  double precision,allocatable  :: kappa_vector(:)
  double precision,allocatable  :: r_unitary_orbrot(:,:)
  double precision :: prev_criterion, criterion
  logical                       :: ex

! Thresholds on convergence of optimized orbitals
  double precision :: max_grad, mean_grad, max_kappa, mean_kappa
  double precision :: norm_grad
  double precision :: max_grad_thresh, mean_grad_thresh, max_kappa_thresh, mean_kappa_thresh
  max_grad_thresh   = 1.d-5
  mean_grad_thresh  = 1.d-6
  max_kappa_thresh  = 1.d-3
  mean_kappa_thresh = 1.d-4

  call check_allowed_Hessian_update(Hessian_update)
  call check_allowed_LA_solver_orb_opt(LA_solver_orb_opt)

! provide mo_two_e_integrals_in_map
  PROVIDE mo_two_e_integrals_in_map ci_energy psi_det psi_coef

  ! MOM related variables
  double precision, allocatable :: mo_coef_reference(:,:)
  logical                       :: swap

  allocate( mo_coef_reference(ao_num,mo_num) )
  mo_coef_reference = mo_coef

  integer :: nT
  nT = dim_list_act_orb
  integer :: nT_tri
  nT_tri = nT * (nT-1) / 2

  ! Renormalization of the weights of the states
  call state_weight_normalization

  ! Compute the criterion before the loop
  call state_average_energy(prev_criterion)

  is_converged = .false.
! Start of orbital optimization loop
  do while( .not.( is_converged .or. iteration.gt.max_iter) )


! Form orbital rotation gradient

  allocate( r_orbrot_g_vector(nT_tri) )
! call gradient_list_opt(tmp_n, m, tmp_list, v_grad, max_elem_grad, norm_grad)
  call gradient_list_opt(nT_tri, nT, list_act, r_orbrot_g_vector, max_grad, norm_grad)

! allocate( r_orbrot_g(nT,nT) )
! call compute_r_orbrot_g(r_orbrot_g,r_one_e_dm_mo,r_two_e_dm_mo,nT)
! call compute_r_orbrot_g_pccd(r_orbrot_g,r_one_e_dm_mo,r_two_e_dm_mo,nT)

!! write(*,*) 'Gradients (2d):'
!! call write_ij_2d_array(r_orbrot_g,nT,nT)

! call transform_2index_to_1index_lowtri(r_orbrot_g,r_orbrot_g_vector,nT)
! deallocate( r_orbrot_g )

  write(*,*) 
! write(*,*) 'Gradients (1d):'
! call write_i_1d_array(r_orbrot_g_vector,nT_tri)


! Write maximum and mean absolute gradient
  write(*,*) 
! max_grad  = maxval( abs(r_orbrot_g_vector(:)) )
  mean_grad = sum( abs(r_orbrot_g_vector(:)) ) / size(r_orbrot_g_vector)
  write(*,*) 'Maximum absolute gradient:     ', max_grad
  write(*,*) 'Mean absolute gradient:        ', mean_grad
  write(*,*) 'Norm absolute gradient:        ', norm_grad

! Start Hessian_update case

! Full Newton-Raphson:
  if( Hessian_update == 'full_Newton_Raphson' ) then

    allocate( r_orbrot_h_square(nT_tri,nT_tri) )
    allocate( r_orbrot_h(nT,nT,nT,nT) )
!   call hessian_list_opt(tmp_n, m, tmp_list, H, h_f)
    call hessian_list_opt(nT_tri, nT, list_act, r_orbrot_h_square, r_orbrot_h)
    deallocate( r_orbrot_h )

!   allocate( r_orbrot_h(nT,nT,nT,nT) )
!   call compute_r_orbrot_h(r_orbrot_h,r_one_e_dm_mo,r_two_e_dm_mo,nT)
!   call compute_r_orbrot_h_pccd(r_orbrot_h,r_one_e_dm_mo,r_two_e_dm_mo,nT)
!   write(*,*) 'Computed orbital rotation Hessian'

!   allocate( r_orbrot_h_square(nT_tri,nT_tri) )
!   call transform_4index_to_2index_lowtri(r_orbrot_h,r_orbrot_h_square,nT,nT)
!   write(*,*) 'Transformed orbital rotation Hessian'
!   deallocate( r_orbrot_h )

  end if

! Diagonal Newton-Raphson:
  if( Hessian_update == 'diagonal' ) then

    allocate( r_orbrot_h_square(nT_tri,nT_tri) )
    allocate( r_orbrot_h(nT,nT,nT,nT) )
    call diag_hessian_list_opt(nT_tri, nT, list_act, r_orbrot_h_square, r_orbrot_h)
    deallocate( r_orbrot_h )

!   double precision, allocatable :: r_orbrot_h_diagonal(:,:)
!   allocate( r_orbrot_h_diagonal(nT,nT) )
!   call compute_r_orbrot_h_diagonal(r_orbrot_h_diagonal,r_one_e_dm_mo,r_two_e_dm_mo,nT)
!!  call compute_r_orbrot_h_diagonal_pccd(r_orbrot_h_diagonal,r_one_e_dm_mo,r_two_e_dm_mo,nT)
!   write(*,*) 'Computed diagonal of the orbital rotation Hessian'

!   double precision, allocatable :: r_orbrot_h_diagonal_square(:)
!   allocate( r_orbrot_h_diagonal_square(nT_tri) )
!   call transform_2index_to_1index_lowtri(r_orbrot_h_diagonal,r_orbrot_h_diagonal_square,nT)
!   deallocate( r_orbrot_h_diagonal )

!   allocate( r_orbrot_h_square(nT_tri,nT_tri) )
!   r_orbrot_h_square = 0.0d0
!   do i=1,nT_tri
!     r_orbrot_h_square(i,i) = r_orbrot_h_diagonal_square(i)
!   end do
!   write(*,*) 'Transformed orbital rotation Hessian'
!   deallocate( r_orbrot_h_diagonal_square )

  end if


! Diagonal with 1's:
  if( Hessian_update == 'identity_diagonal' ) then
    allocate( r_orbrot_h_square(nT_tri,nT_tri) )
    r_orbrot_h_square = 0.0d0
    do i=1,nT_tri
      r_orbrot_h_square(i,i) = 1.0d0
    end do
  end if

! End Hessian_update case

! Add diagonal mu_damping:
! call add_shift_to_diagonal(r_orbrot_h_square,nT_tri,mu_damping)


!! write(*,*) 'r_orbrot_h_square:'
!! call write_ij_2d_array(r_orbrot_h_square,nT_tri,nT_tri)

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
! do i=1,nT_tri
!   if( abs(hessian_eigvalues(i)).lt.1.0d-10 ) hessian_eigvalues(i) = 1.0d20
!   if( abs(hessian_eigvalues(i)).lt.1.0d-10 ) hessian_eigvalues(i) = 1.0d80
! end do

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
! H is real symmetric => H = U Z U^-1 = U Z U^T, or U^-1 = U^T
! H-1 g = (U Z U^T)-1 g = ( U^T^-1 Z^-1 U^-1 ) g = ( U Z^-1 U^T ) g 

! Let us use kappa_vector as a buffer:
!   do i=1,nT_tri
!     kappa_vector(i) = sum( hessian_eigvectors(:,i) * r_orbrot_g_vector(:) )
!   end do
!
    double precision, allocatable :: tmp_wtg(:)
    allocate( tmp_wtg(nT_tri) )
    tmp_wtg = 0.0d0
    do i=1,nT_tri
      tmp_wtg(i) = sum( hessian_eigvectors(:,i) * r_orbrot_g_vector(:) )
    end do

!   write(*,*) 'Projection of gradients (2d):'
!!  call write_i_1d_array(kappa_vector,nT_tri)
!   call write_i_1d_array(tmp_wtg,nT_tri)
 
!   do i=1,nT_tri
!     if( hessian_eigvalues(i) .gt. 0.0d0 ) then
!       hessian_eigvalues(i) = hessian_eigvalues(i) + mu_damping
!     else
!       hessian_eigvalues(i) = hessian_eigvalues(i) - mu_damping
!     end if
!   end do

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
!       r_orbrot_g_vector(i) = kappa_vector(i) / hessian_eigvalues(i)
        r_orbrot_g_vector(i) = tmp_wtg(i) / hessian_eigvalues(i)
      end do
    end if
!   deallocate( hessian_eigvalues )

! And finally the actual kappa_vector:
!   do i=1,nT_tri
!     kappa_vector(i) = - sum( hessian_eigvectors(i,:) * r_orbrot_g_vector(:) )
!   end do
    kappa_vector = 0.0d0
    do i=1,nT_tri
      if( hessian_eigvalues(i) < -1.0d-14 ) then
        do j=1,nT_tri
          kappa_vector(j) = kappa_vector(j) + tmp_wtg(i) * hessian_eigvectors(j,i) / ( hessian_eigvalues(i) - mu_damping )
        end do
      else
        do j=1,nT_tri
          kappa_vector(j) = kappa_vector(j) - tmp_wtg(i) * hessian_eigvectors(j,i) / ( hessian_eigvalues(i) + mu_damping )
        end do
      end if
    end do
    deallocate( hessian_eigvectors )
    deallocate( hessian_eigvalues )

    deallocate( tmp_wtg ) 

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

  write(*,*) 'kappa_vector:'
  call write_i_1d_array(kappa_vector,nT_tri)

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
  ! 1D tmp -> 2D tmp 
  call vec_to_mat_v2(nT_tri, nT, kappa_vector, kappa_square)
! call form_A_square(kappa_vector,nT_tri,kappa_square,nT)

  deallocate( kappa_vector )

!     write(*,*) 'kappa_square'
!     do i=1,20
!       write(*,'(20(f10.5,x))') (kappa_square(j,i), j=1,20)
!     end do

! Damp rotation between orbitals orb_a, orb_b and all the others, by damp_rotation
  integer :: orb_a, orb_b
  double precision :: damp_rotation
  orb_a = 4
  orb_b = 5
  damp_rotation = 0.001d0
! do i=1,mo_num
!   kappa_square(i,orb_a) = kappa_square(i,orb_a) * damp_rotation
!   kappa_square(orb_a,i) = kappa_square(i,orb_a)
!   kappa_square(i,orb_b) = kappa_square(i,orb_b) * damp_rotation
!   kappa_square(orb_b,i) = kappa_square(i,orb_b)
! end do

! write(*,*) 'kappa_square:'
! call write_ij_2d_array(kappa_square,nT,nT)

! Form the unitary orbital rotation matrix from the antihermitian kappa matrix: U = exp( kappa )
  allocate( r_unitary_orbrot(nT,nT) )
  call lapack_exp_antisymm_matrix(nT,kappa_square,r_unitary_orbrot)
  write(*,*) 'Formed the rotation matrix U = exp( kappa )'
  deallocate( kappa_square )

!     write(*,*) 'rotation matrix'
!     do i=1,20
!       write(*,'(20(f10.5,x))') (r_unitary_orbrot(j,i), j=1,20)
!     end do

! Check the rotation matrix r_unitary_orbrot
  logical :: rotation_matrix_ok
  call check_rotation_matrix(r_unitary_orbrot,nT,rotation_matrix_ok)
  if(.not.rotation_matrix_ok) then
!   if( mu_damping == 0.0d0 ) then
!     mu_damping = 1.0d0
!   else
!     mu_damping = mu_damping * 2.0d0
!   end if
    deallocate( r_unitary_orbrot )
    cycle
  end if

! Final rotation matrix
  double precision, allocatable :: R(:,:)
  allocate( R(mo_num,mo_num) )

! r_unitary_orbrot to R, active space to full space
  call sub_to_full_rotation_matrix(nT, list_act, r_unitary_orbrot, R)
  deallocate( r_unitary_orbrot )


! Apply the orbital rotation U to the molecular orbitals: d' = d U
! call update_molecular_orbitals(mo_coef,ao_num,nT,r_unitary_orbrot)
  call update_molecular_orbitals(mo_coef,ao_num,mo_num,R)
  deallocate( R )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   call run_ci_mom(mo_coef_reference,swap)

  ! Swap orbitals
!   if( swap ) then

!     call clear_mo_map
!     TOUCH mo_coef psi_det psi_coef ci_energy two_e_dm_mo

      ! Diagonalization of the hamiltonian
!     FREE ci_energy! To enforce the recomputation
!     call diagonalize_ci
!     call save_wavefunction_unsorted

      ! Energy obtained after the diagonalization of the CI matrix
!     call state_average_energy(prev_criterion)

!   end if

! Update the reference MOs for the original version of MOM. Do nothing for initial MOM (IMOM)
!   if( mom_type == 'MOM' ) then
!     mo_coef_reference = mo_coef
!   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FREE ci_energy! To enforce the recomputation

  call clear_mo_map
  TOUCH mo_coef
! TOUCH mo_coef psi_det psi_coef ci_energy two_e_dm_mo
! call state_average_energy(criterion)

  ! Diagonalization of the hamiltonian
! FREE ci_energy! To enforce the recomputation
  call diagonalize_ci
  call save_wavefunction_unsorted

  ! Energy obtained after the diagonalization of the CI matrix
! call state_average_energy(prev_criterion)

! Make sure the orbitals are orthonormal after the rotation:
  call check_mos_orthonormality

! After one (quasi-)Newton-Raphson step, save molecular orbitals to the EZFIO directory:
  call save_mos
  write(*,*) 'Updated the molecular orbitals'

  write(*,*) '#############################'
  write(*,*) 'Done iteration ', iteration
  write(*,*) '#############################'
  iteration = iteration + 1

  end do
! End of orbital optimization loop

end subroutine optimize_orbitals_ci


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


