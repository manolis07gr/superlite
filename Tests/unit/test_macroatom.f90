!-- Macro-atom consistency test (nltehtrmod): in a photoionized cell the
!-- energy the macro-atom emits in H-alpha and H-beta per recombination
!-- activation must reproduce the line emissivities of the statistical-
!-- equilibrium populations (Lucy 2002), and the energy must be conserved
!-- (outcome probabilities sum to 1).
program test_macroatom
  use randommod
  use nltehmod
  use nltehtrmod
  use gridmod, only:grd_tempinv
  implicit none
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
  integer,parameter :: nlev = 30, ng = 6000, ndraw = 2000000
  real*8 :: wlg(ng+1), jg(ng), x(0:nlev), ne, beta(nlev*(nlev-1)/2), jb(nlev*(nlev-1)/2)
  logical :: ljb(nlev*(nlev-1)/2)
  real*8 :: t, nh, q, w, t2, jha, jhb, rmc, rse, err
  integer :: i, k, ierr, kind, idx, kha, khb, nha, nhb, nfail
  real*8 :: wlc
  type(rnd_t) :: st
  nfail = 0
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  call nlteh_tr_init(nlev, 1, 1, ng, wlg)
!-- photoionized low-density cell, Lyman lines thick (as test_nlteh case B)
  t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-9
  jg = 0d0
  do k = 1, ng
     if(nht_nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nht_nug(k), t2)
  enddo
  jb = 0d0; ljb = .false.
  x = 0d0; x(0) = 0.99d0; x(1) = 0.01d0
  call nlteh_solve_cell(t, t2, nh, 0d0, q, q, ng, nht_nug, nht_dnug, jg, jb, ljb, 0d0, x, ne, beta, ierr)
  if(ierr /= 0) stop 'solve failed'
  nht_x(:,1) = x; nht_nh(1) = nh; nht_ne(1) = ne
  nht_jg(:,1) = jg; nht_jb(:,1) = 0d0; nht_ljb(:,1) = .false.; nht_have_mc = .true.
  allocate(grd_tempinv(1)); grd_tempinv = 1d0/t
  call nlteh_tables(1, [1d0/t], [t2], [0d0], [1d15, 1.001d15], [q*1d15, q*1.001d15])
!-- probabilities are cumulative and end at 1
  do k = 1, nht_nst
     if(abs(nht_acum(nht_nout,k,1) - 1d0) > 1d-12 .or. any(nht_acum(2:,k,1) < nht_acum(:nht_nout-1,k,1) - 1d-14)) then
        write(*,*) 'bad cumulative table for state', k
        nfail = nfail + 1
     endif
  enddo
!-- recombination activations
  call rnd_init(1, 3); st = rnd_states(1)
  kha = nh_lin(2,3); khb = nh_lin(2,4)
  nha = 0; nhb = 0
  do i = 1, ndraw
     call nlteh_activate(1, 0, st, kind, idx, wlc)
     if(kind == em_line .and. idx == kha) nha = nha + 1
     if(kind == em_line .and. idx == khb) nhb = nhb + 1
  enddo
  rmc = dble(nha)/dble(nhb)
  err = rmc*sqrt(1d0/nha + 1d0/nhb)
  jha = x(3)*nh_a(kha)*nht_beta(kha,1)*nh_nu(kha)
  jhb = x(4)*nh_a(khb)*nht_beta(khb,1)*nh_nu(khb)
  rse = jha/jhb
  write(*,'(a,f8.4,a,f7.4,a,f8.4)') '  H-alpha/H-beta energy from the macro-atom:', rmc, ' +-', err, &
       '   from the populations:', rse
  write(*,'(a,f7.4,a,f7.4)') '  fraction of recombination energy in H-alpha:', dble(nha)/ndraw, &
       '   H-beta:', dble(nhb)/ndraw
  if(abs(rmc - rse) > 4d0*err + 0.005d0*rse) nfail = nfail + 1
  if(nfail > 0) then
     write(*,'(a,i3)') 'test_macroatom: FAILED checks:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_macroatom: PASS'
end program test_macroatom
