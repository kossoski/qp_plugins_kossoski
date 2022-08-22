!
! This subroutine is adapted from the orbital optimization subroutine developed by Yann Damour: https://github.com/Ydrnan/qp_plugins_damour
!

! Subroutine : run_orb_opt_trust



! Subroutine to optimize the MOs using a trust region algorithm:
! - choice of the method
! - initialization
! - optimization until convergence

! The optimization use the trust region algorithm, the different parts
! are explained in the corresponding subroutine files.

! qp_edit:
! | thresh_opt_max_elem_grad |
! | optimization_max_nb_iter |
! | optimization_method      |

! Provided:
! | mo_num                         | integer          | number of MOs                  |
! | ao_num                         | integer          | number of AOs                  |
! | N_states                       | integer          | number of states               |
! | ci_energy(N_states)            | double precision | CI energies                    |
! | state_average_weight(N_states) | double precision | Weight of the different states |

! Variables:
! | m                         | integer          | number of active MOs                               |
! | tmp_n                     | integer          | m*(m-1)/2, number of MO parameters                 |
! | v_grad(tmp_n)             | double precision | gradient                                           |
! | H(tmp_n,tmp_n)            | double precision | hessian (2D)                                       |
! | h_f(m,m,m,m)              | double precision | hessian (4D)                                       |
! | e_val(m)                  | double precision | eigenvalues of the hessian                         |
! | w(m,m)                    | double precision | eigenvectors of the hessian                        |
! | x(m)                      | double precision | step given by the trust region                     |
! | m_x(m,m)                  | double precision | step given by the trust region after               |
! | tmp_R(m,m)                | double precision | rotation matrix for active MOs                     |
! | R(mo_num,mo_num)          | double precision | full rotation matrix                               |
! | prev_mos(ao_num,mo_num)   | double precision | previous MOs (before the rotation)                 |
! | new_mos(ao_num,mo_num)    | double precision | new MOs (after the roration)                       |
! | delta                     | double precision | radius of the trust region                         |
! | rho                       | double precision | agreement between the model and the exact function |
! | max_elem                  | double precision | maximum element in the gradient                    |
! | i                         | integer          | index                                              |
! | tmp_i,tmp_j               | integer          | indexes in the subspace containing only            |
! |                           |                  | the active MOs                                     |
! | converged                 | logical          | convergence of the algorithm                       |
! | cancel_step               | logical          | if the step must be cancelled                      |
! | nb_iter                   | integer          | number of iterations (accepted)                    |
! | nb_diag                   | integer          | number of diagonalizations of the CI matrix        |
! | nb_cancel                 | integer          | number of cancelled steps for the actual iteration |
! | nb_cancel_tot             | integer          | total number of cancel steps                       |
! | info                      | integer          | if 0 ok, else problem in the diagonalization of    |
! |                           |                  | the hessian with the Lapack routine                |
! | criterion                 | double precision | energy at a given step                             |
! | prev_criterion            | double precision | energy before the rotation                         |
! | criterion_model           | double precision | estimated energy after the rotation using          |
! |                           |                  | a Taylor series                                    |
! | must_exit                 | logical          | To exit the trust region algorithm when            |
! |                           |                  | criterion - criterion_model is too small           |
! | enforce_step_cancellation | logical          | To force the cancellation of the step if the       |
! |                           |                  | error in the rotation matrix is too large          |

subroutine run_ci_mom_orb_opt_trust_v2

  include 'constants.h'

  implicit none

  ! Variables

  double precision, allocatable :: R(:,:)
  double precision, allocatable :: H(:,:),h_f(:,:,:,:)
  double precision, allocatable :: v_grad(:)
  double precision, allocatable :: prev_mos(:,:),new_mos(:,:)
  integer                       :: info
  integer                       :: n
  integer                       :: i,j,p,q,k
  double precision              :: max_elem_grad, delta, rho, norm_grad, normalization_factor
  logical                       :: cancel_step
  integer                       :: nb_iter, nb_diag, nb_cancel, nb_cancel_tot, nb_sub_iter
  double precision              :: t1, t2, t3
  double precision              :: prev_criterion, criterion, criterion_model
  logical                       :: not_converged, must_exit, enforce_step_cancellation
  integer                       :: m, tmp_n, tmp_i, tmp_j, tmp_k
  integer,allocatable           :: tmp_list(:)
  double precision, allocatable :: tmp_m_x(:,:),tmp_R(:,:), tmp_x(:), W(:,:), e_val(:)

  PROVIDE mo_two_e_integrals_in_map ci_energy psi_det psi_coef

! Method
!    There are three different methods : 
!    - the "full" hessian, which uses all the elements of the hessian
!      matrix"
!    - the "diagonal" hessian, which uses only the diagonal elements of the
!      hessian
!    - without the hessian (hessian = identity matrix)

!Display the method
 print*, 'Method :', optimization_method
if (optimization_method == 'full') then 
  print*, 'Full hessian'
elseif (optimization_method == 'diag') then
  print*,'Diagonal hessian'
elseif (optimization_method == 'none') then
  print*,'No hessian'
else
  print*,'Unknown optimization_method, please select full, diag or none'
  call abort
endif
print*, 'Absolute value of the hessian:', absolute_eig

! Allocation

allocate(R(mo_num,mo_num))  ! rotation matrix
allocate(prev_mos(ao_num,mo_num), new_mos(ao_num,mo_num)) ! old and new MOs

! Definition of m and tmp_n
m = dim_list_act_orb
tmp_n = m*(m-1)/2

allocate(tmp_list(m))
allocate(tmp_R(m,m), tmp_m_x(m,m), tmp_x(tmp_n))
allocate(H(tmp_n,tmp_n), h_f(m,m,m,m), v_grad(tmp_n),W(tmp_n,tmp_n),e_val(tmp_n))

! Algorithm

! Here is the main algorithm of the optimization:
! - First of all we initialize some parameters and we compute the
!   criterion (the ci energy) before doing any MO rotations
! - We compute the gradient and the hessian for the active MOs
! - We diagonalize the hessian
! - We compute a step and loop to reduce the radius of the
!   trust region (and the size of the step by the way) until the step is
!   accepted 
! - We repeat the process until the convergence 
!   NB: the convergence criterion can be changed

! Loop until the convergence of the optimization
  ! call diagonalize_ci 

  !### Initialization ###
  nb_iter = 0
  rho = 0.5d0
  not_converged = .True.
  tmp_list = list_act ! Optimization of the active MOs
  nb_cancel_tot = 0

  ! Renormalization of the weights of the states
  call state_weight_normalization

  !if (normalized_st_av_weight) then
  !  ! To put the same weight everywhere
  !  do i = 1, N_states
  !    state_average_weight(i) = 1d0/DBLE(N_states)
  !  enddo
  !  TOUCH state_average_weight
  !elseif (normalized_st_weight) then
  !  ! To normalize the weights
  !  normalization_factor = 0d0
  !  do i = 1, N_states
  !    normalization_factor = normalization_factor + state_average_weight(i)**2
  !  enddo
  !  normalization_factor = 1d0 / dsqrt(normalization_factor)
  !  
  !  do i = 1, N_states
  !    state_average_weight(i) = state_average_weight(i) * normalization_factor
  !  enddo
  !  TOUCH state_average_weight
  !else
  !  ! Nothing
  !endif
  !print*, 'Number of states:', N_states
  !print*, 'State average weights:'
  !print*, state_average_weight(:)

  ! Compute the criterion before the loop
  call state_average_energy(prev_criterion)

  ! MOM related variables
  double precision, allocatable :: mo_coef_reference(:,:)
  logical                       :: swap 

  allocate( mo_coef_reference(ao_num,mo_num) )
  mo_coef_reference = mo_coef

! Identify which orbitals belong to the singly occupied and doubly occupied spaces
! TODO
! do i=1,n_det
!   psi_det(1:N_int,1,i)
!   psi_det(1:N_int,2,i)
! end do
! n_singly_occ = 1
! n_doubly_occ = 2

  do while (not_converged)
    print*,''
    print*,'******************'
    print*,'Iteration', nb_iter
    print*,'******************'
    print*,''

    ! Gradient
    call gradient_list_opt(tmp_n, m, tmp_list, v_grad, max_elem_grad, norm_grad)
    
    ! Hessian
    if (optimization_method == 'full') then
      ! Full hessian
      call hessian_list_opt(tmp_n, m, tmp_list, H, h_f)
    elseif (optimization_method == 'diag') then
      ! Diagonal hessian 
      call diag_hessian_list_opt(tmp_n, m, tmp_list, H, h_f)
    else
      ! Identity matrix 
      H = 0d0
      do tmp_i = 1, tmp_n
        H(tmp_i,tmp_i) = 1d0
      enddo
    endif
 
    ! Diagonalization of the hessian 
    call diagonalization_hessian(tmp_n, H, e_val, w)

    ! Init before the internal loop
    cancel_step = .True. ! To enter in the loop just after 
    nb_cancel = 0
    nb_sub_iter = 0

    ! Loop to reduce the trust radius until the criterion decreases and rho >= thresh_rho
    do while (cancel_step)
      print*,'-----------------------------'
      print*,'Iteration:', nb_iter
      print*,'Sub iteration:', nb_sub_iter
      print*,'-----------------------------'

      ! Hessian,gradient,Criterion -> x 
      call trust_region_step_w_expected_e(tmp_n,H,W,e_val,v_grad,prev_criterion,rho,nb_iter,delta,criterion_model,tmp_x,must_exit) 

      if (must_exit) then
        print*,'step_in_trust_region sends the message : Exit'
        exit
      endif

      ! 1D tmp -> 2D tmp 
      call vec_to_mat_v2(tmp_n, m, tmp_x, tmp_m_x)

      ! Rotation matrix for the active MOs
      call rotation_matrix(tmp_m_x, m, tmp_R, m, m, info, enforce_step_cancellation)

      ! Security to ensure an unitary transformation
      if (enforce_step_cancellation) then
        print*, 'Forces the step cancellation, too large error in the rotation matrix'
        rho = 0d0
        cycle
      endif

      ! tmp_R to R, subspace to full space
      call sub_to_full_rotation_matrix(m, tmp_list, tmp_R, R)
    
      ! MO rotations
      call apply_mo_rotation(R, prev_mos)   

      ! Update of the energy before the diagonalization of the hamiltonian
      call clear_mo_map
      TOUCH mo_coef psi_det psi_coef ci_energy two_e_dm_mo 
      call state_average_energy(criterion)

      ! Criterion -> step accepted or rejected 
      call trust_region_is_step_cancelled(nb_iter, prev_criterion, criterion, criterion_model, rho, cancel_step)

      ! Cancellation of the step if necessary
      if (cancel_step) then
        mo_coef = prev_mos
        call save_mos()
        nb_cancel = nb_cancel + 1
        nb_cancel_tot = nb_cancel_tot + 1
      else
        ! Diagonalization of the hamiltonian
        FREE ci_energy! To enforce the recomputation
        call diagonalize_ci
        call save_wavefunction_unsorted

        ! Energy obtained after the diagonalization of the CI matrix
        call state_average_energy(prev_criterion)
      endif

      ! Internal loop exit conditions
      if (nb_cancel > nb_cancel_max) then
        print*,'###################################'
        print*,'nb_cancel >', nb_cancel_max,', exit'
        print*,'###################################'
        must_exit = .True.
        exit
      endif
 
      if (nb_cancel_tot > nb_cancel_tot_max) then
        print*,'###########################################'
        print*,'nb_cancel_tot >', nb_cancel_tot_max,', exit'
        print*,'###########################################'
        must_exit = .True.
        exit
      endif

      nb_sub_iter = nb_sub_iter + 1
    enddo

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

    call save_mos() !### depend of the time for 1 iteration

    ! To exit the external loop if must_exit = .True.
    if (must_exit) then
      exit
    endif 

    ! Step accepted, nb iteration + 1
    nb_iter = nb_iter + 1

    ! External loop exit conditions
    if (DABS(max_elem_grad) < thresh_opt_max_elem_grad) then
      print*,'Converged: DABS(max_elem_grad) < thresh_opt_max_elem_grad'
      not_converged = .False.
    endif
    if (nb_iter >= optimization_max_nb_iter) then
      print*,'Not converged: nb_iter >= optimization_max_nb_iter'
      not_converged = .False. 
    endif
    if (prev_criterion - criterion < thresh_opt_energy_gain .and. DABS(max_elem_grad) < thresh_opt_max_elem_grad * 10d0) then
      print*,'Converged: prev_criterion - criterion < thresh_opt_energy_gain .and. DABS(max_elem_grad) < thresh_opt_max_elem_grad * 10d0'
      not_converged = .False.
    endif

    if (.not. not_converged) then
      print*,'#############################'
      print*,'   End of the optimization'
      print*,'#############################'
    endif
  enddo

!  call diagonalize_ci

! Deallocation, end

deallocate(v_grad,H,R,W,e_val)
  deallocate(h_f,prev_mos,new_mos)

end
