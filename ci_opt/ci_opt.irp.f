program ci_opt

  implicit none
  BEGIN_DOC
  ! Optimize the orbitals of a general CI wave function
  END_DOC

! Local variables
  integer :: iteration = 0
  logical :: is_converged

  double precision,allocatable  :: gradient(:)
  double precision,allocatable  :: hessian(:,:)
  double precision,allocatable  :: kappa(:)
  double precision,allocatable  :: rotation_matrix(:,:)
  double precision              :: prev_criterion, criterion
  integer                       :: nT, nT_tri, nT_tri1
  integer                       :: i, j

! Thresholds on convergence of optimized orbitals
  double precision :: max_grad, mean_grad, max_kappa, mean_kappa
  double precision :: max_grad_thresh, mean_grad_thresh, max_kappa_thresh, mean_kappa_thresh

  read_wf = .true.
  TOUCH read_wf

  max_grad_thresh   = 1.d-5
  mean_grad_thresh  = 1.d-6
  max_kappa_thresh  = 1.d-3
  mean_kappa_thresh = 1.d-4

  call check_allowed_Hessian_update(Hessian_update)
  call check_allowed_LA_solver_orb_opt(LA_solver_orb_opt)

  provide mo_two_e_integrals_in_map
! provide mo_two_e_integrals_in_map ci_energy psi_det psi_coef

! Diagonalization of the hamiltonian
  FREE ci_energy
  call diagonalize_ci
  call save_wavefunction_unsorted

! MOM related variables
  double precision, allocatable :: mo_coef_reference(:,:)
  logical                       :: swap
  allocate( mo_coef_reference(ao_num,mo_num) )
  mo_coef_reference = mo_coef

  nT = dim_list_act_orb
  nT_tri = nT * (nT-1) / 2

! Renormalization of the weights of the states
! call state_weight_normalization

! Compute the criterion before the loop
  call state_average_energy(prev_criterion)

  is_converged = .false.

! Start of orbital optimization loop
  do while( .not.( is_converged .or. iteration==max_iter) )


! Computes orbital rotation gradient
  allocate( gradient(nT_tri) )
  double precision :: norm_grad
  call gradient_list_opt(nT_tri, nT, list_act, gradient, max_grad, norm_grad)
  max_grad  = maxval( abs(gradient(:)) )
  mean_grad = sum( abs(gradient(:)) ) / size(gradient)

! Prints the gradient
  if( debug_orbital_optimization ) then
    write(*,*) 
    write(*,*) 'Gradients (1d):'
    call write_i_1d_array(gradient,nT_tri)
  end if


! Computes orbital rotation Hessian

  allocate( hessian(nT_tri,nT_tri) )

! Full Newton-Raphson:
  double precision,allocatable  :: hessian_4d(:,:,:,:)
  if( Hessian_update == 'full_Newton_Raphson' ) then
    allocate( hessian_4d(nT,nT,nT,nT) )
    call hessian_list_opt(nT_tri, nT, list_act, hessian, hessian_4d)
    deallocate( hessian_4d )
  end if

! Diagonal Newton-Raphson:
  if( Hessian_update == 'diagonal' ) then
    allocate( hessian_4d(nT,nT,nT,nT) )
    call diag_hessian_list_opt(nT_tri, nT, list_act, hessian, hessian_4d)
    deallocate( hessian_4d )
  end if

! Identity:
  if( Hessian_update == 'identity_diagonal' ) then
    call set_identity_matrix(hessian,nT_tri)
  end if

! Prints the Hessian:
! if( debug_orbital_optimization ) then
! write(*,*) 
! write(*,*) 'Hessian (2d)'
! call write_ij_2d_array(hessian,nT_tri,nT_tri)
! end if

! Checks whether the Hessian matrix is symmetric:
  call check_matrix_symmetry(hessian,nT_tri)

! Diagonalizes the Hessian matrix and writes its eigenvalues
  double precision, allocatable :: hessian_eigvalues(:)
  double precision, allocatable :: hessian_eigvectors(:,:)
  allocate( hessian_eigvalues(nT_tri) )
  allocate( hessian_eigvectors(nT_tri,nT_tri) )
  call lapack_diag(hessian_eigvalues,hessian_eigvectors,hessian,nT_tri,nT_tri)
! call lapack_diagd(hessian_eigvalues,hessian_eigvectors,hessian,nT_tri,nT_tri)

  integer :: negatives
  call count_numbers_below_thresh(hessian_eigvalues,nT_tri,-1.d-7,negatives)
  write(*,*) 
  write(*,*) 'Number of negative Hessian eigenvalues below -1.0d-7: ', negatives
  write(*,*) 

  do i=1,nT_tri
    if( hessian_eigvalues(i).gt.-1.d-7 ) exit
    write(*,*) 'Largest elements of the eigenvalue ', i
    write(*,*) 'Eigenvalue = ', hessian_eigvalues(i)
    call print_largest_eigenvectors(abs(hessian_eigvectors(:,i)),nT_tri)
    write(*,*) 
  end do

! Prints the Hessian eigenvalue:
  if( debug_orbital_optimization ) then
    write(*,*) 
    write(*,*) 'Hessian eigenvalues:'
    call write_i_1d_array(hessian_eigvalues,nT_tri)
!   write(*,*) 'Hessian eigenvectors:'
!   call write_ij_2d_array(hessian_eigvectors,nT_tri,nT_tri)
  end if


  allocate( kappa(nT_tri) )
  kappa = 0.0d0

  double precision, allocatable :: aug_hessian(:,:)
  double precision, allocatable :: aug_hessian_eigvalues(:)
  double precision, allocatable :: aug_hessian_eigvectors(:,:)

  integer :: saddle_order1
  saddle_order1 = saddle_order + 1

! Start LA_solver_orb_opt case

  if( LA_solver_orb_opt=='diagonalize_augmented' ) then

    nT_tri1 = nT_tri + 1
! Build the augmented Hessian matrix
    allocate( aug_hessian(nT_tri1,nT_tri1) )
    call build_augmented_matrix(hessian,gradient,nT_tri,aug_hessian)

! Diagonalize the augmented Hessian matrix
    allocate( aug_hessian_eigvalues(nT_tri1) )
    allocate( aug_hessian_eigvectors(nT_tri1,nT_tri1) )
    call lapack_diag(aug_hessian_eigvalues,aug_hessian_eigvectors,aug_hessian,nT_tri1,nT_tri1)
    deallocate( aug_hessian_eigvectors )
    deallocate( aug_hessian )

    if( debug_orbital_optimization ) then
      write(*,*)
      write(*,*) 'Augmented Hessian eigenvalues:'
      call write_i_1d_array(aug_hessian_eigvalues,nT_tri1)
    end if

! Let us use kappa as a buffer:
    do i=1,nT_tri
      kappa(i) = sum( hessian_eigvectors(:,i) * gradient(:) )
    end do

!   write(*,*) 'Projection of gradients (2d):'
!   call write_i_1d_array(kappa,nT_tri)

! Let us use gradient as a buffer:
    do i=1,nT_tri
      gradient(i) = kappa(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(saddle_order1) * lambda_hessian )
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

! And finally the actual kappa:
    do i=1,nT_tri
      kappa(i) = - sum( hessian_eigvectors(i,:) * gradient(:) )
    end do
    deallocate( hessian_eigvectors )

  end if



  if( LA_solver_orb_opt=='diagonalize_augmented_p' ) then

    double precision, allocatable :: kappa_aux(:)
    allocate( kappa_aux(nT_tri) )
! Let us use kappa as a buffer:
    do i=1,nT_tri
      kappa_aux(i) = sum( hessian_eigvectors(:,i) * gradient(:) )
    end do

! Negative eigevalues:

!   write(*,*) 'Projection of (negative) gradients (2d):'
!   call write_i_1d_array(kappa_aux(1:saddle_order),saddle_order)

    deallocate( hessian )
    allocate( hessian(saddle_order,saddle_order) )
    hessian = 0.0d0
    do i=1,saddle_order
      hessian(i,i) = hessian_eigvalues(i)
    end do

! Build the augmented (negative) Hessian matrix
    allocate( aug_hessian(saddle_order1,saddle_order1) )
    call build_augmented_matrix(hessian,kappa_aux(1:saddle_order),saddle_order,aug_hessian)

! Diagonalize the augmented Hessian matrix
    allocate( aug_hessian_eigvalues(saddle_order1) )
    allocate( aug_hessian_eigvectors(saddle_order1,saddle_order1) )
    call lapack_diag(aug_hessian_eigvalues,aug_hessian_eigvectors,aug_hessian,saddle_order1,saddle_order1)
    deallocate( aug_hessian_eigvectors )
    deallocate( aug_hessian )
    if( debug_orbital_optimization ) then
      write(*,*)
      write(*,*) 'Augmented (negative) Hessian eigenvalues:'
      call write_i_1d_array(aug_hessian_eigvalues,saddle_order1)
    end if

! Let us use gradient as a buffer:
    do i=1,saddle_order
      gradient(i) = kappa_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(saddle_order1) )
!     gradient(i) = kappa_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(saddle_order1) * lambda_hessian )
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

    if( debug_orbital_optimization ) then
      write(*,*) 'Projection of the (negative) gradient:'
      call write_i_1d_array(gradient(1:saddle_order),saddle_order)
    end if

! And finally the actual kappa:
    do i=1,nT_tri
      kappa(i) = - sum( hessian_eigvectors(i,1:saddle_order) * gradient(1:saddle_order) )
    end do

    if( debug_orbital_optimization ) then
      write(*,*) 'Projection of the (negative) kappa:'
      call write_i_1d_array(kappa(1:saddle_order),saddle_order)
    end if

! Positive eigevalues:
   
    integer :: comp_saddle_order, comp_saddle_order1
    comp_saddle_order  = nT_tri - saddle_order
    comp_saddle_order1 = comp_saddle_order + 1

! Let us use kappa as a buffer:
!   do i=saddle_order1,nT_tri
!     kappa(i) = sum( hessian_eigvectors(:,i) * gradient(:) )
!   end do

    if( debug_orbital_optimization ) then
      write(*,*) 'Projection of the (positive) gradients:'
      call write_i_1d_array(kappa_aux(saddle_order1:nT_tri),comp_saddle_order)
    end if

    deallocate( hessian )
    allocate( hessian(comp_saddle_order,comp_saddle_order) )
    hessian = 0.0d0
    do i=1,comp_saddle_order
      hessian(i,i) = hessian_eigvalues(saddle_order+i)
    end do

! Build the augmented (positive) Hessian matrix
    allocate( aug_hessian(comp_saddle_order1,comp_saddle_order1) )
    call build_augmented_matrix(hessian,kappa_aux(saddle_order1:nT_tri),comp_saddle_order,aug_hessian)

! Diagonalize the augmented Hessian matrix
    allocate( aug_hessian_eigvalues(comp_saddle_order1) )
    allocate( aug_hessian_eigvectors(comp_saddle_order1,comp_saddle_order1) )

    call lapack_diag(aug_hessian_eigvalues,aug_hessian_eigvectors,aug_hessian,comp_saddle_order1,comp_saddle_order1)
    deallocate( aug_hessian_eigvectors )
    deallocate( aug_hessian )
    if( debug_orbital_optimization ) then
      write(*,*)
      write(*,*) 'Augmented (positive) Hessian eigenvalues:'
      call write_i_1d_array(aug_hessian_eigvalues(1),1)
!     call write_i_1d_array(aug_hessian_eigvalues,comp_saddle_order1)
    end if

! Let us use gradient as a buffer:
    do i=saddle_order1,nT_tri
!     gradient(i) = kappa_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(1) )
      gradient(i) = kappa_aux(i) / ( hessian_eigvalues(i) - aug_hessian_eigvalues(1) * lambda_hessian )
    end do
    deallocate( aug_hessian_eigvalues )
    deallocate( kappa_aux )

! And finally the actual kappa:
    do i=1,nT_tri
      kappa(i) = kappa(i) - sum( hessian_eigvectors(i,saddle_order1:nT_tri) * gradient(saddle_order1:nT_tri) )
    end do

    if( debug_orbital_optimization ) then
      write(*,*) 'Projection of the (positive) kappa:'
      call write_i_1d_array(kappa(saddle_order1:nT_tri),comp_saddle_order)
    end if

    deallocate( hessian_eigvalues )
    deallocate( hessian_eigvectors )

  end if


  if( LA_solver_orb_opt=='diagonalize' ) then

! From the diagonalization:
! H is real symmetric => H = U Z U^-1 = U Z U^T, or U^-1 = U^T
! H-1 g = (U Z U^T)-1 g = ( U^T^-1 Z^-1 U^-1 ) g = ( U Z^-1 U^T ) g 

! Let us use kappa as a buffer:
    do i=1,nT_tri
      kappa(i) = sum( hessian_eigvectors(:,i) * gradient(:) )
    end do
 
    do i=1,nT_tri
      if( hessian_eigvalues(i) < -1.0d-14 ) then
!       hessian_eigvalues(i) = - hessian_eigvalues(i) + mu_damping
        hessian_eigvalues(i) = + hessian_eigvalues(i) - mu_damping
      else
        hessian_eigvalues(i) = + hessian_eigvalues(i) + mu_damping
      end if
    end do

    if( debug_orbital_optimization ) then
      write(*,*) 'Projection of gradients:'
      call write_i_1d_array(kappa,nT_tri)
    end if
 
    if( move_negative_Hessian ) then
      do i=1,min_negative_direction-1
        gradient(i) = kappa(i) / hessian_eigvalues(i)
      end do
      do i=min_negative_direction,max_negative_direction
        if( hessian_eigvalues(i) .le. 0.0d0 ) then
!         gradient(i) = sign( step_negative_direction, kappa(i) )
          gradient(i) = sign( - step_negative_direction / hessian_eigvalues(i), kappa(i) )
        else
          gradient(i) = kappa(i) / hessian_eigvalues(i)
        end if
      end do
      do i=max_negative_direction+1,nT_tri
        gradient(i) = kappa(i) / hessian_eigvalues(i)
      end do
    else
! Let us use gradient as a buffer:
      do i=1,nT_tri
        gradient(i) = kappa(i) / hessian_eigvalues(i)
      end do
    end if
!   deallocate( hessian_eigvalues )

! And finally the actual kappa:
    do i=1,nT_tri
      kappa(i) = - sum( hessian_eigvectors(i,:) * gradient(:) )
    end do

    deallocate( hessian_eigvectors )
    deallocate( hessian_eigvalues )

  end if


  if( LA_solver_orb_opt=='inverse' ) then

    double precision,allocatable  :: hessian_inv(:,:)
    allocate( hessian_inv(nT_tri,nT_tri) )
    call get_inverse(hessian,nT_tri,nT_tri,hessian_inv,nT_tri)

    do i=1,nT_tri
      kappa(i) = - sum( hessian_inv(i,:) * gradient(:) )
    end do
    deallocate( hessian_inv )

  end if


  if( LA_solver_orb_opt=='linear_system' ) then
    stop 'Not yet coded for LA_solver_orb_opt=linear_system'
! TODO: Does not work if I try to obtain kappa by solving the linear system with this subroutine
!   call lapack_solve_linear_dsy(nT_tri,hessian,-gradient,kappa)
  end if

! End LA_solver_orb_opt case

  deallocate( hessian )
  deallocate( gradient )


  max_kappa  = maxval( abs(kappa(:)) )
  mean_kappa = sum( abs(kappa(:)) ) / size(kappa)
! Prints the displacement vector kappa
  if( debug_orbital_optimization ) then
    write(*,*) 'kappa:'
    call write_i_1d_array(kappa,nT_tri)
  end if

! Prints information on the gradient and kappa
  write(*,*) 
  write(*,*) 'Maximum absolute gradient: ', max_grad
  write(*,*) 'Mean    absolute gradient: ', mean_grad
  write(*,*) 'Maximum absolute kappa:    ', max_kappa
  write(*,*) 'Mean    absolute kappa:    ', mean_kappa
  write(*,*) 

! Test for convergence
  if( max_grad .le.max_grad_thresh  .and. mean_grad .le.mean_grad_thresh .and. &
      max_kappa.le.max_kappa_thresh .and. mean_kappa.le.mean_kappa_thresh ) then
    is_converged = .true.
  end if

  if( LA_solver_orb_opt=='diagonalize_augmented_p' .or. LA_solver_orb_opt=='diagonalize_augmented_p' ) then
     if( .not.correct_saddle_order ) is_converged = .false.
  end if


  double precision,allocatable  :: kappa_2d(:,:)
  allocate( kappa_2d(nT,nT) )
  call vec_to_mat_v2(nT_tri, nT, kappa, kappa_2d)

  deallocate( kappa )

! Form the unitary orbital rotation matrix from the antihermitian kappa matrix: U = exp( kappa )
  allocate( rotation_matrix(nT,nT) )
  call lapack_exp_antisymm_matrix(nT,kappa_2d,rotation_matrix)
  deallocate( kappa_2d )
  write(*,*) 'Formed the rotation matrix U = exp( kappa )'

! Check the orthonormality of the rotation matrix
  logical :: rotation_matrix_ok
  call check_orthornormality_rotation_matrix(rotation_matrix,nT,rotation_matrix_ok)

! Full rotation matrix
  double precision, allocatable :: R(:,:)
  allocate( R(mo_num,mo_num) )
! rotation_matrix to R, active space to full space
  call sub_to_full_rotation_matrix(nT, list_act, rotation_matrix, R)
  deallocate( rotation_matrix )

! Apply the orbital rotation U to the molecular orbitals: mo_coef' = mo_coef R
  call update_molecular_orbitals(mo_coef,ao_num,mo_num,R)
  deallocate( R )


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
  call state_average_energy(prev_criterion)

! Make sure the orbitals are orthonormal after the rotation:
  call check_mos_orthonormality

! After one (quasi-)Newton-Raphson step, save molecular orbitals to the EZFIO directory:
  call save_mos
! write(*,*) 'Updated the molecular orbitals'

  iteration = iteration + 1
  write(*,*) 
  write(*,*) '#############################'
  write(*,*) 'Done iteration ', iteration
  write(*,*) '#############################'
  write(*,*) 

  end do
! End of orbital optimization loop

!end subroutine optimize_orbitals_ci


end program

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


subroutine print_n_lowest_matrx_elements(A,n,n_lowest)
  implicit none
  BEGIN_DOC
! 
  END_DOC
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n,n)
  integer,          intent(in)  :: n_lowest
  double precision, allocatable :: A_copy(:,:)
! double precision, allocatable :: A_1d(:)
  integer         , allocatable :: iorder(:)
! double precision              :: A_1d(n*n)
! integer                       :: iorder(n*n)
  integer                       :: n_tri
  integer                       :: i, j, k, l, p, q, r, s, pq, rs, ij

! n_tri = n * (n-1) / 2
  n_tri = n * n

  ! Set array of natural numbers
  allocate( iorder(n_tri) )
  call set_integer_list(iorder,n_tri)
! call set_integer_list(iorder,n*n)

  allocate( A_copy(n,n) )
  A_copy = A
! allocate( A_1d(n_tri) )
! A_1d = reshape( A_copy, shape(A_1d) )
! A_1d = reshape( A_copy, (/1/) )

  ! Sort A_1d in increasing order
! call dsort(A_1d,iorder,n_tri)
! call dsort(A_1d,iorder,n*n)

  call dsort(A_copy,iorder,n_tri)
  write(*,*) 'Smallest matrix elements of the Hessian:'
  do k=1,n_lowest
    pq = mod(k-1,n) + 1
    rs = (k-1)/n + 1
    l = iorder(k)
    i = mod(l-1,n) + 1
    j = (l-1)/n + 1
    call vec_to_mat_index(j,r,s)
    call vec_to_mat_index(i,p,q)
    write(*,'(5(i5,x),e16.8)') k, p, q, r, s, A_copy(pq,rs)
!   write(*,'(5(i5,x),e16.8)') k, p, q, r, s, A(pq,rs)
!   write(*,'(7(i5,x),e16.8)') k, i, j, p, q, r, s, A(pq,rs)
!   write(*,*) k, pq, rs, A(pq,rs)
  end do

end subroutine

subroutine print_largest_eigenvectors(A,n)
  implicit none
  BEGIN_DOC
! 
  END_DOC
  integer,          intent(in)  :: n
  double precision, intent(in)  :: A(n)
  double precision, allocatable :: A_copy(:)
  double precision              :: thresh = 0.2d0
  integer,          allocatable :: iorder(:)
  integer                       :: k, l, i, j

! Set array of natural numbers
  allocate( iorder(n) )
  call set_integer_list(iorder,n)

  allocate( A_copy(n) )
  A_copy = A

! Sort A_copy in increasing order
  call dsort(A_copy,iorder,n)

  do k=n,1,-1
    if( A_copy(k).lt.thresh ) exit
    l = iorder(k)
    call vec_to_mat_index(l,i,j)
!   write(*,'(2(i5,x),e16.8)') i, j, A_copy(k)
    write(*,'(2(i5,x),e16.8)') list_act(i), list_act(j), A_copy(k)
  end do

end subroutine
