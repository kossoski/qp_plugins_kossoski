program ci_mom
  read_wf = .true.
  TOUCH read_wf
! call run_ci_mom_orb_opt_trust_v2
  call optimize_orbitals_ci
end
