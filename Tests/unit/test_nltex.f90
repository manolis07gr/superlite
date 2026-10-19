!-- Tests of the multi-species NLTE treatment (nlteatommod, nltexmod,
!-- nltextrmod):
!--  1. built-in hydrogen reproduces the H-only treatment (nltehmod):
!--     recombination coefficients, photoionization rates, populations in the
!--     LTE limit and in a Case B nebula
!--  2. with He model atoms (directory given as the first argument; skipped
!--     without): LTE limit at high density, He I 5876/4471 in a Case B
!--     nebula, and the macro-atom: cumulative tables, and the 5876/4471
!--     energy ratio of recombination activations against the populations
program test_nltex
  use nltehmod
  use nlteatommod
  use nltexmod
  use nltextrmod
  use nltedismod
  use randommod
  use gridmod, only:grd_tempinv
  implicit none
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
  real*8,parameter :: pi_ = 3.141592653589793d0
  integer,parameter :: nlev = 30, ng = 6000
  real*8 :: wlg(ng+1), nug(ng), dnug(ng), jg(ng)
  real*8 :: xo(0:nlev), betao(nlev*(nlev-1)/2), jbo(nlev*(nlev-1)/2), neo
  logical :: ljbo(nlev*(nlev-1)/2)
  real*8,allocatable :: x(:), xl(:), beta(:), jb(:), nti(:), ntx(:), nel(:)
  type(nxrt) :: rr
  real*8,allocatable :: wv(:)
  logical,allocatable :: ljb(:)
  real*8 :: gamo(nlev), alpo(nlev), alpso(nlev), cupo(nlev), c3bo(nlev), qupo(nlev,nlev), qdno(nlev,nlev)
  real*8 :: jlo(nlev*(nlev-1)/2)
  type(nxion) :: ions(3)
  real*8 :: t, nh, w, t2, q, ne, dev, d, tl(3), nl(3)
  integer :: i, k, n, ierr, nfail, nargs
  character(256) :: adir
  nfail = 0
!-- group grid as in the A4 test (1 A - 30000 A, logarithmic)
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i = 1, ng
     nug(i) = 0.5d0*(cc/wlg(i) + cc/wlg(i+1))
     dnug(i) = cc/wlg(i) - cc/wlg(i+1)
  enddo
!
!-- 1. built-in hydrogen
  call nlteh_init(nlev)
  call nxa_builtin_h(nlev, ions(1))
  call nxa_setup(1, ions(1:1))
  call nxa_groups(ng, wlg)
  allocate(x(nxa_nst), xl(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst))
  allocate(nel(nxa_nel))
!-- recombination coefficients against a fine integration in ln u (the
!-- integrand of high levels at high T falls steeply above the edge)
  dev = 0d0
  d = 0d0
  tl = [5d3, 1d4, 3d4]
  do i = 1, 3
     t = tl(i)
     do n = 1, nlev, 4
        w = alpha_ref(n, t)
        dev = max(dev, abs(nxa_alpha(n, t)/w - 1d0))
        d = max(d, abs(nlteh_alpha(n, t)/w - 1d0))
     enddo
  enddo
  write(*,'(a,es10.3,a,es10.3,a)') '  max |alpha/alpha_ref - 1| =', dev, '   (nltehmod:', d, ')'
  if(dev > 3d-3) nfail = nfail + 1
!-- photoionization rates in a diluted 5e4 K field
  w = 1d-14; t2 = 5d4
  do k = 1, ng
     jg(k) = w*bplanck(nug(k), t2)
  enddo
  jbo = 0d0; ljbo = .false.
  call nlteh_cell_rates(1d4, t2, ng, nug, dnug, jg, jbo, ljbo, 0d0, gamo, alpo, alpso, cupo, c3bo, qupo, qdno, jlo)
  call nx_el_rates(1, 1d4, t2, ng, nug, dnug, jg, rr)
  dev = maxval(abs(rr%gam(1:nlev)/gamo - 1d0))
  d = maxval(abs(rr%cq(1,2:nlev)/qupo(1,2:nlev) - 1d0))
  write(*,'(a,es10.3,a,es10.3)') '  max |Gamma/Gamma(nltehmod) - 1| =', dev, '   collisions:', d
  if(dev > 1d-2 .or. d > 1d-2) nfail = nfail + 1
!-- LTE limit and Case B: populations against nlteh_solve_cell
  tl = [5d3, 1d4, 2d4]; nl = [1d17, 1d17, 1d18]
  do i = 1, 4
     if(i <= 3) then
        t = tl(i); nh = nl(i); t2 = t; q = 1d0/1.7d6
        do k = 1, ng
           jg(k) = bplanck(nug(k), t)
        enddo
        jbo = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]
        ljbo = .true.
     else
!-- Case B: nH = 1e4, T = 1e4 K, ionizing diluted 5e4 K blackbody only
        t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-9
        jg = 0d0
        do k = 1, ng
           if(nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nug(k), t2)
        enddo
        jbo = 0d0; ljbo = .true.
     endif
     xo = 0d0; xo(0) = 0.5d0; xo(1) = 0.5d0
     call nlteh_solve_cell(t, t2, nh, 0d0, q, q, ng, nug, dnug, jg, jbo, ljbo, 0d0, xo, neo, betao, ierr)
     x = 0d0; x(1) = 0.5d0; x(nxa_nst) = 0.5d0
     jb = jbo; ljb = ljbo; nti = 0d0; ntx = 0d0; nel = nh
     call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
     dev = max(abs(x(nxa_nst)/xo(0) - 1d0), maxval(abs(x(1:nlev)/xo(1:nlev) - 1d0)))
     write(*,'(a,i2,a,es9.2,a,f7.0,a,i2,a,es10.3,a,es10.3)') '  H case', i, ' nH=', nh, ' T=', t, '  ierr=', &
          ierr, '  x(HII)=', x(nxa_nst), '  max |x/x(nltehmod) - 1| =', dev
     if(i <= 3) then
        if(ierr /= 0 .or. dev > 1d-2) nfail = nfail + 1
     else
!-- Case B Balmer decrement (Storey & Hummer 1995: 2.86)
        k = ions(1)%lin(2,3); n = ions(1)%lin(2,4)
        w = (x(3)*nxa_la(k)*beta(k)*nxa_lnu(k))/(x(4)*nxa_la(n)*beta(n)*nxa_lnu(n))
        write(*,'(a,f7.3,a)') '  H Case B H-alpha/H-beta =', w, '   (Storey & Hummer 1995: 2.86)'
        if(ierr /= 0 .or. abs(w/2.86d0 - 1d0) > 0.05d0) nfail = nfail + 1
     endif
  enddo
  call pc_tests
  call td_tests
  call ali_tests
  call ldb_tests
  call dis_tests
  deallocate(x, xl, beta, jb, ljb, nti, ntx, nel)
!
!-- 2. He model atoms
  nargs = command_argument_count()
  if(nargs < 1) then
     write(*,'(a)') '  (no model-atom directory given: He tests skipped)'
  else
     call get_command_argument(1, adir)
     call he_tests(trim(adir))
     call all_lte(trim(adir))
  endif
  if(nfail > 0) then
     write(*,'(a,i3)') 'test_nltex: FAILED checks:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_nltex: PASS'
contains
!
  subroutine td_tests
!-- time-dependent statistical equilibrium (built-in H, Case B nebula):
!-- the steady state is a solution for any dt; dt -> 0 keeps the previous
!-- populations; one backward-Euler step of recombination in the dark,
!-- x = x0/(1 + alpha nH x dt), against its root
    real*8 :: xss(nxa_nst), xp(nxa_nst), at, a, xa, d1, d2, d3
    integer :: lv
    t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-9
    jg = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nug(k), t2)
    enddo
    jb = 0d0; ljb = .true.; nti = 0d0; ntx = 0d0; nel = nh
    xss = 0d0; xss(1) = 0.5d0; xss(nxa_nst) = 0.5d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xss, ne, beta, ierr)
    x = 0d0; x(1) = 0.3d0; x(nxa_nst) = 0.7d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr, &
         xprev=xss, dtinv=1d-6)
    d1 = maxval(abs(x/xss - 1d0), mask=xss>1d-20)
!-- dt = 1e-3 s from an ionized state: frozen
    xp = 0d0; xp(1) = 1d-3; xp(nxa_nst) = 1d0 - 1d-3
    x = xp
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr, &
         xprev=xp, dtinv=1d3)
    d2 = abs(x(nxa_nst)/xp(nxa_nst) - 1d0)
!-- recombination in the dark from x(H II) = 1 - 1e-3 over dt = 1e8 s
    jg = 0d0
    at = 0d0
    do lv = 1, nlev
       at = at + nxa_alpha(lv, t)
    enddo
    a = at*nh*1d8
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr, &
         xprev=xp, dtinv=1d-8)
!-- x = xp/(1 + a x): a x^2 + x - xp = 0
    xa = (-1d0 + sqrt(1d0 + 4d0*a*xp(nxa_nst)))/(2d0*a)
    d3 = abs(x(nxa_nst)/xa - 1d0)
    write(*,'(a,es9.2,a,es9.2,a,f7.4,a,f7.4,a)') '  time-dependent SE: steady state kept to', d1, &
         ', dt -> 0 frozen to', d2, ', recombination in the dark x(H II) =', x(nxa_nst), &
         ' (backward Euler', xa, ')'
    if(d1 > 1d-6 .or. d2 > 1d-5 .or. d3 > 1d-2) nfail = nfail + 1
  end subroutine td_tests
!
  subroutine ali_tests
!-- accelerated lambda iteration of the bound-free transitions (built-in H,
!-- 1e4 K, nH = 1e10, diluted 5e4 K field above the Lyman edge): at the
!-- populations with which J was made the corrections cancel, so the
!-- solution must not move; with J made of the cell's own Lyman-continuum
!-- emission (J = Lambda* S + the rest) the iteration with the corrections
!-- must reach the fixed point of the plain iteration in fewer steps
    real*8 :: x0(nxa_nst), kap(ng), lst(ng), x1(nxa_nst), wv(nxa_nlv), d1
    real*8 :: jext(ng), xa(nxa_nst), xb(nxa_nst), da, db, sl
    integer :: lv, ig, itn, na, nb
    t = 1d4; nh = 1d10; w = 1d-14; t2 = 5d4; q = 1d-9
    jext = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jext(k) = w*bplanck(nug(k), t2)
    enddo
    jb = 0d0; ljb = .true.; nti = 0d0; ntx = 0d0; nel = nh
    jg = jext
    x0 = 0d0; x0(1) = 0.5d0; x0(nxa_nst) = 0.5d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x0, ne, beta, ierr)
!-- Lyman-continuum opacity of x0 and a local operator of 0.9
    call nxa_weights(t, wv)
    kap = 0d0
    do lv = 1, nlev
       do ig = 1, nxa_g1(lv)
          kap(ig) = kap(ig) + nh*x0(nxa_lvs(lv))*wv(lv)*nxa_sgm(nxa_sgo(lv)+ig)
       enddo
    enddo
    lst = 0.9d0
    if(allocated(nx_ali_lst)) deallocate(nx_ali_lst, nx_ali_kap)
    allocate(nx_ali_lst(ng), nx_ali_kap(ng))
    nx_ali_lst = lst; nx_ali_kap = kap
    nx_ali_on = .true.
    x1 = x0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x1, ne, beta, ierr)
    nx_ali_on = .false.
    d1 = maxval(abs(x1/x0 - 1d0), mask=x0>1d-20)
!-- lambda iteration: J = jext + Lambda* S_Lyc(x) in the Lyman continuum,
!-- S = eta/kappa of the ground state (Milne), plain and accelerated
    na = 0; nb = 0
    xa = 0d0; xa(1) = 0.999d0; xa(nxa_nst) = 1d-3
    xb = xa
    da = 1d0; db = 1d0
    do itn = 1, 200
       x1 = xa
       call lyc_j(xa, jg)
       call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xa, ne, beta, ierr)
       da = abs(xa(nxa_nst)/x1(nxa_nst) - 1d0)
       if(da < 1d-6) then
          na = itn
          exit
       endif
    enddo
    do itn = 1, 200
       x1 = xb
       call lyc_j(xb, jg)
       call set_ali(xb)
       call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xb, ne, beta, ierr)
       nx_ali_on = .false.
       db = abs(xb(nxa_nst)/x1(nxa_nst) - 1d0)
       if(db < 1d-6) then
          nb = itn
          exit
       endif
    enddo
    sl = abs(xb(nxa_nst)/xa(nxa_nst) - 1d0)
    write(*,'(a,es9.2,a,i4,a,i4,a,f8.5,a,f8.5,a,es9.2)') '  ALI: fixed point kept to', d1, &
         '; lambda iteration converged in', na, ' (plain) and', nb, ' (ALI) steps, x(H II) =', &
         xa(nxa_nst), ' and', xb(nxa_nst), ', difference', sl
    if(d1 > 1d-6 .or. na == 0 .or. nb == 0 .or. nb >= na .or. sl > 1d-4) nfail = nfail + 1
  end subroutine ali_tests
!
  subroutine pc_tests
!-- continuum destruction of the trapped line photons (nx_pcont):
!--  a) P(r) against a quadrature of (2/sqrt(pi)) int_0^inf r/(1 + r e^{x^2}) dx;
!--  b) LTE (J = B in the groups and the lines): the populations do not move
!--     when a continuum takes the trapped photons (P S_c = P S_line);
!--  c) Case B nebula with a continuum that takes every trapped photon at
!--     the lines (P -> 1) where J = 0: the lines are as if thin (the same
!--     solution as with a velocity gradient that makes beta = 1)
    real*8 :: rr(14), pq, xx, dx, dvp, xa(nxa_nst), xb(nxa_nst), d1, d2
    integer :: i2, j2
    rr = [1d-14, 1d-10, 3.3d-7, 1d-6, 1d-3, 0.05d0, 1d-1, 1d0, 2.6d0, 1d1, 3.1d2, 1d3, 1d8, 1d12]
    dvp = 0d0
    do i2 = 1, 14
       pq = 0d0
       dx = 1d-4
       do j2 = 1, 200000
          xx = (j2 - 0.5d0)*dx
          pq = pq + rr(i2)/(1d0 + rr(i2)*exp(min(xx*xx, 700d0)))*dx
       enddo
       pq = 2d0/sqrt(pi_)*pq
       dvp = max(dvp, abs(nx_pcont(rr(i2))/pq - 1d0))
    enddo
!-- b) LTE at 1e4 K, n_H = 1e17 (as H case 2), continuum of 1e-5 /cm in every group
    t = 1d4; nh = 1d17; t2 = t; q = 1d0/1.7d6
    do k = 1, ng
       jg(k) = bplanck(nug(k), t)
    enddo
    jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.; nti = 0d0; ntx = 0d0; nel = nh
    xa = 0d0; xa(1) = 0.5d0; xa(nxa_nst) = 0.5d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xa, ne, beta, ierr)
    allocate(nx_kcg(ng)); nx_kcg = 1d-5; nx_pc_vtb = 1d6; nx_pc_on = .true.
    xb = 0d0; xb(1) = 0.5d0; xb(nxa_nst) = 0.5d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xb, ne, beta, ierr)
    d1 = maxval(abs(xb/xa - 1d0), mask=xa > 1d-30)
!-- c) Case B (as H case 4) with P -> 1 at every line, against beta -> 1
    t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-9
    jg = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nug(k), t2)
    enddo
    jb = 0d0; ljb = .true.; nel = nh
    nx_kcg = 1d10
    xa = 0d0; xa(1) = 0.5d0; xa(nxa_nst) = 0.5d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xa, ne, beta, ierr)
    nx_pc_on = .false.
    xb = 0d0; xb(1) = 0.5d0; xb(nxa_nst) = 0.5d0
    call nx_solve_cell(t, t2, nel, 0d0, 1d10, 1d10, ng, nug, dnug, jg, jb, ljb, nti, ntx, xb, ne, beta, ierr)
    d2 = maxval(abs(xb/xa - 1d0), mask=xa > 1d-30)
    deallocate(nx_kcg); nx_pc_vtb = 0d0
    write(*,'(a,es9.2,a,es9.2,a,es9.2,a,f8.5)') '  continuum destruction: P(r) against quadrature to', dvp, &
         '; LTE kept to', d1, '; P -> 1 against thin lines to', d2, ', x(H II)', xa(nxa_nst)
    if(dvp > 1d-3 .or. d1 > 1d-6 .or. d2 > 1d-4) nfail = nfail + 1
  end subroutine pc_tests
!
  subroutine ldb_tests
!-- local detailed balance in thick groups (nxt_ldb, nxt_ali_cell):
!--  a) the share Lambda* of its own emission that a slab keeps, against
!--     1 - (1 - 2 E_3(tau))/(2 tau) with E_3 by quadrature;
!--  b) J of a one-cell grid: Lambda* S + J_ext in the thick groups (S +
!--     J_ext with DDMC), the DDMC lump tally shared as the source
!--     spectrum, the others kept;
!--  c) the toy cell of ali_tests with a cell width that makes its Lyman
!--     continuum thick: J = jext + Lambda* S (Lambda* of the slab), with a
!--     Monte Carlo estimate that misses the cell's own emission in the
!--     groups above nxt_ldbtau (J = jext there, as when few packets are
!--     emitted in them) the gas is too neutral (case A); with local
!--     detailed balance it reaches the fixed point of the full J in a few
!--     steps
    real*8,parameter :: sb_ = 5.670374419d-5
    real*8 :: tt(4), e3, sv, ds, dev1, dev2, u, sg, jw, w0, dsl
    real*8 :: vol1(1), dr1(1), ti1(1)
    real*8,allocatable :: plk(:,:), spk1(:,:)
    real*4,allocatable :: cap1(:,:), cape1(:,:)
    integer*2,allocatable :: gl1(:,:)
    integer :: nl1(1), i2, j2, ig, na, nb
    real*8 :: xa(nxa_nst), xb(nxa_nst), xc(nxa_nst), x1(nxa_nst), kap(ng), sl(ng), dra, da, db
!-- a) slab share
    tt = [0.1d0, 2d0, 5d0, 30d0]
    dev1 = 0d0
    do i2 = 1, 4
       e3 = 0d0
       ds = 60d0/2d6
       do j2 = 1, 2000000
          sv = (j2 - 0.5d0)*ds
          e3 = e3 + exp(-sv)/(1d0 + sv/tt(i2))**3*ds
       enddo
       e3 = e3*exp(-tt(i2))/tt(i2)
       dev1 = max(dev1, abs((1d0 - nxt_lsldb(tt(i2)))/((1d0 - 2d0*e3)/(2d0*tt(i2))) - 1d0))
    enddo
!-- b) one-cell grid: groups 1-10 thick (tau = 100), lump = groups 1-20
    allocate(plk(ng,1), spk1(ng,1), cap1(ng,1), cape1(ng,1), gl1(ng,1))
    if(allocated(nxt_jg)) deallocate(nxt_jg)
    if(allocated(nxt_jxsum)) deallocate(nxt_jxsum)
    if(allocated(nxt_jlsum)) deallocate(nxt_jlsum)
    if(allocated(nxt_dnug)) deallocate(nxt_dnug)
    if(allocated(grd_emnbf)) deallocate(grd_emnbf)
    allocate(nxt_jg(ng,1), nxt_jxsum(ng,1), nxt_jlsum(1), nxt_dnug(ng), grd_emnbf(ng,1))
    nxt_dnug = dnug
    vol1 = 1d45; dr1 = 1d14; ti1 = 1d-4
    plk(:,1) = 1d0/ng
    spk1(:,1) = [(1d0 + 0.01d0*ig, ig=1,ng)]
    spk1 = spk1/sum(spk1)
    cap1 = 1e-16; cap1(1:10,1) = 1e-12
    cape1 = 0.5*cap1
    grd_emnbf(:,1) = 0.3*cap1(:,1)*real(plk(:,1))
    nxt_jg = 1d0
    nxt_jxsum(:,1) = 1d20*[(dble(ig), ig=1,ng)]
    nxt_jlsum = 7d33
    nl1 = 20
    gl1(:,1) = [(int(ig,2), ig=1,ng)]
    nxt_ldbtau = 30d0
    call nxt_ldb(1, ng, vol1, dr1, ti1, cap1, cape1, plk, spk1, nl1, gl1)
    u = sb_/(pi_*ti1(1)**4)
    w0 = sum(spk1(1:20,1))
    dev2 = 0d0
    do ig = 1, 30
       if(ig <= 10) then
          sg = (plk(ig,1)*cape1(ig,1) + grd_emnbf(ig,1))/cap1(ig,1)*u/dnug(ig)
          jw = nxt_lsldb(dble(cap1(ig,1))*dr1(1))*sg + nxt_jxsum(ig,1)/(4d0*pi_*vol1(1)*dnug(ig))
       elseif(ig <= 20) then
          jw = 1d0 + nxt_jlsum(1)*spk1(ig,1)/w0/(4d0*pi_*vol1(1)*dnug(ig))
       else
          jw = 1d0
       endif
       dev2 = max(dev2, abs(nxt_jg(ig,1)/jw - 1d0))
    enddo
!-- (with DDMC the thick groups are in the diffusion limit: Lambda* = 1)
    nxt_ldbone = .true.
    nxt_jg = 1d0
    call nxt_ldb(1, ng, vol1, dr1, ti1, cap1, cape1, plk, spk1, nl1, gl1)
    do ig = 1, 10
       sg = (plk(ig,1)*cape1(ig,1) + grd_emnbf(ig,1))/cap1(ig,1)*u/dnug(ig)
       jw = sg + nxt_jxsum(ig,1)/(4d0*pi_*vol1(1)*dnug(ig))
       dev2 = max(dev2, abs(nxt_jg(ig,1)/jw - 1d0))
    enddo
!-- (and the net emission of a line inside a thick group is part of its
!-- source function: Lyman beta in a thick group, S = (kappa S_c + eta_l)/
!-- (kappa + kappa_l) from the populations and the escape probability)
    allocate(nxt_x(nxa_nst,1), nxt_nel(nxa_nel,1), nxt_w(nxa_nlv,1), nxt_beta(nxa_nlin,1), nxt_pc(nxa_nlin,1), &
         nxt_jgc(ng,1))
    nxt_pc = 0d0
    nxt_x = 0d0
    nxt_x(1,1) = 1d-3; nxt_x(2,1) = 1d-9; nxt_x(3,1) = 2d-10; nxt_x(nxa_nst,1) = 1d0 - 1.0000012d-3
    nxt_nel = 1d9
    nxt_w = 1d0
    nxt_beta = 0d0
    do k = 1, nxa_nlin
       if(nxa_lsl(k)==1 .and. nxa_lsu(k)==3) exit
    enddo
    nxt_beta(k,1) = 0.3d0
    ig = nxa_lgn(k)
    cap1(ig,1) = 1e-12
    cape1(ig,1) = 0.5*cap1(ig,1)
    grd_emnbf(ig,1) = 0.3*cap1(ig,1)*real(plk(ig,1))
    nxt_jg = 1d0
    call nxt_ldb(1, ng, vol1, dr1, ti1, cap1, cape1, plk, spk1, nl1, gl1)
    sg = (plk(ig,1)*cape1(ig,1) + grd_emnbf(ig,1))/cap1(ig,1)*u/dnug(ig)
    da = 0.3d0*hh*nxa_lnu(k)/(4d0*pi_*dnug(ig))
    db = da*(1d9*1d-3*nxa_lblu(k) - 1d9*2d-10*nxa_lbul(k))
    da = da*1d9*2d-10*nxa_la(k)
    jw = (cap1(ig,1)*sg + da)/(cap1(ig,1) + db) + nxt_jxsum(ig,1)/(4d0*pi_*vol1(1)*dnug(ig))
    dsl = nxt_jg(ig,1)/jw - 1d0
!-- (and the lines see the continuum alone: nxt_jgc = S_c + J_ext)
    dsl = max(abs(dsl), abs(nxt_jgc(ig,1)/(sg + nxt_jxsum(ig,1)/(4d0*pi_*vol1(1)*dnug(ig))) - 1d0))
!-- (J moves from S_c towards the line's source function da/db)
    w0 = nxt_jg(ig,1) - nxt_jxsum(ig,1)/(4d0*pi_*vol1(1)*dnug(ig))
    write(*,'(a,es10.3,a,es10.3,a,es10.3)') '  line inside a thick DDMC group: J/(S_c + line terms) - 1 and J_c/S_c - 1 to', dsl, &
         '; (J - J_ext)/S_c =', w0/sg,', S_line/S_c =', da/db/sg
    if(abs(dsl) > 1d-12 .or. (w0 - sg)*(da/db - sg) <= 0d0) nfail = nfail + 1
    deallocate(nxt_x, nxt_nel, nxt_w, nxt_beta, nxt_pc, nxt_jgc)
    nxt_ldbone = .false.
    nxt_ldbtau = 0d0
    deallocate(plk, spk1, cap1, cape1, gl1)
    deallocate(nxt_jg, nxt_jxsum, nxt_jlsum, nxt_dnug, grd_emnbf)
!-- c) toy cell, Lyman continuum thick (cell width dra)
    t = 1d4; nh = 1d10; w = 1d-14; t2 = 5d4; q = 1d-9; dra = 1d13
    jb = 0d0; ljb = .true.; nti = 0d0; ntx = 0d0; nel = nh
!-- (case A: the Monte Carlo estimate without the cell's own emission in
!-- the thick groups)
    nxt_ldbtau = 30d0
    xa = 0d0; xa(1) = 0.5d0; xa(nxa_nst) = 0.5d0
    do i2 = 1, 5000
       x1 = xa
       call lyc_sk(xa, sl, kap, jg)
       do ig = 1, ng
          if(kap(ig) > 0d0 .and. kap(ig)*dra < nxt_ldbtau) jg(ig) = jg(ig) + nxt_lsldb(kap(ig)*dra)*sl(ig)
       enddo
       call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xa, ne, beta, ierr)
       if(abs(xa(nxa_nst)/x1(nxa_nst) - 1d0) < 1d-9) exit
    enddo
!-- (reference: plain lambda iteration of J = jext + Lambda*(tau) S)
    xb = 0d0; xb(1) = 0.5d0; xb(nxa_nst) = 0.5d0
    na = 0
    do i2 = 1, 5000
       x1 = xb
       call lyc_sk(xb, sl, kap, jg)
       do ig = 1, ng
          if(kap(ig) > 0d0) jg(ig) = jg(ig) + nxt_lsldb(kap(ig)*dra)*sl(ig)
       enddo
       call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xb, ne, beta, ierr)
       da = abs(xb(nxa_nst)/x1(nxa_nst) - 1d0)
       if(da < 1d-9) then
          na = i2
          exit
       endif
    enddo
!-- (local detailed balance: J of nxt_ldb, own cycle out with nxt_ali_cell;
!-- the thin groups with the full J of a well-sampled estimate)
    xc = 0d0; xc(1) = 0.5d0; xc(nxa_nst) = 0.5d0
    nb = 0
    do i2 = 1, 200
       x1 = xc
       call lyc_sk(xc, sl, kap, jg)
       do ig = 1, ng
          if(kap(ig) > 0d0) jg(ig) = jg(ig) + nxt_lsldb(kap(ig)*dra)*sl(ig)
       enddo
       call nxt_ali_cell(ng, real(kap,4), 1d0, dra, lali=.false.)
       call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, xc, ne, beta, ierr)
       call nxt_ali_off
       db = abs(xc(nxa_nst)/x1(nxa_nst) - 1d0)
       if(db < 1d-9) then
          nb = i2
          exit
       endif
    enddo
    nxt_ldbtau = 0d0
    write(*,'(a,es9.2,a,es9.2)') '  LDB: slab share against E_3 quadrature to', dev1, &
         '; J of the one-cell grid to', dev2
    write(*,'(a,f9.6,a,f9.6,a,i5,a,f9.6,a,i3,a)') '  LDB toy cell: x(H II) case A', xa(nxa_nst), &
         ', reference', xb(nxa_nst), ' (', na, ' steps), local detailed balance', xc(nxa_nst), ' (', nb, ' steps)'
    if(dev1 > 1d-4 .or. dev2 > 1d-6 .or. na == 0 .or. nb == 0 .or. nb >= na) nfail = nfail + 1
    if(abs(xc(nxa_nst)/xb(nxa_nst) - 1d0) > 1d-5 .or. abs(xa(nxa_nst)/xb(nxa_nst) - 1d0) < 1d-2) nfail = nfail + 1
  end subroutine ldb_tests
!
  subroutine dis_tests
!-- level dissolution (nltedismod):
!--  a) 1 - w of hydrogen and He II levels against CMFGEN's OCCUPATION_PROB
!--     (values computed independently from its formula, n_e = N_ion);
!--  b) the LTE limit of hydrogen at nH = 1e15 cm^-3, T = 1e4 K in a
!--     Planck field, where the levels above n = 10 dissolve: the
!--     statistical equilibrium with the dissolved LTE populations, shares
!--     and continua below the edges gives the dissolved Saha-Boltzmann
!--     populations (detailed balance of all the dissolution terms)
    real*8 :: p(3), dref(6), dval, dmx
    integer :: j
    real*8,parameter :: nn(5) = [8d0, 10d0, 15d0, 20d0, 30d0]
    nxd_on = .true.
    call nxd_setup(ng, wlg)
!-- a) n_e = 1e13, 1e15 cm^-3, T = 1e4 K; the last value He II n = 20
    dmx = 0d0
    dref = [2.501514d-04, 1.255859d-03, 2.771890d-02, 3.422844d-01, 9.925372d-01, 9.158527d-03]
    call nxd_cell(1d13, 1d13, 1d4, p)
    do j = 1, 5
       dval = 1d0 - nxd_wocc(nn(j), 1d0, p)
       dmx = max(dmx, abs(dval/dref(j) - 1d0))
    enddo
    dmx = max(dmx, abs((1d0 - nxd_wocc(20d0, 2d0, p))/dref(6) - 1d0))
    dref = [2.932838d-02, 2.030115d-01, 9.803182d-01, 9.996923d-01, 9.999992d-01, 8.764411d-01]
    call nxd_cell(1d15, 1d15, 1d4, p)
    do j = 1, 5
       dval = 1d0 - nxd_wocc(nn(j), 1d0, p)
       dmx = max(dmx, abs(dval/dref(j) - 1d0))
    enddo
    dmx = max(dmx, abs((1d0 - nxd_wocc(20d0, 2d0, p))/dref(6) - 1d0))
    write(*,'(a,es10.3)') '  dissolution: max |(1-w)/(1-w)_CMFGEN - 1| =', dmx
    if(dmx > 1d-5) nfail = nfail + 1
!-- b) LTE limit with dissolution
    t = 1d4; nh = 1d15
    do k = 1, ng
       jg(k) = bplanck(nug(k), t)
    enddo
    jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.
    nti = 0d0; ntx = 0d0; nel = nh
    nx_dis = 0d0
    call nx_lte(t, nh, x)
    ne = nh*x(nxa_nst)
!-- (the cell's parameters follow the solution: three passes)
    do j = 1, 3
       call nxd_cellx(ne, nel, x, t, nx_dis)
       call nx_solve_cell(t, t, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
    enddo
    call nxd_cellx(ne, nel, x, t, nx_dis)
    call nx_lte(t, ne, xl)
    dev = maxval(abs(x(1:nxa_nst)/xl(1:nxa_nst) - 1d0), mask=xl(1:nxa_nst) > 1d-30)
    write(*,'(a,es10.3,a,es10.3,a,es10.3,a,i2)') '  dissolution LTE limit: n_e =', ne, &
         '  1-w(n=12) =', 1d0 - nxd_lvw(12, nx_dis), '  max |x/x_LTE - 1| =', dev, '  ierr=', ierr
    if(ierr /= 0 .or. dev > 1d-2 .or. nxd_lvw(12, nx_dis) > 0.99d0) nfail = nfail + 1
    nx_dis = 0d0
    nxd_on = .false.
  end subroutine dis_tests
!
  subroutine lyc_sk(xx, ss, kk, jj)
!-- toy cell: source function and opacity of the ground state's Lyman
!-- continuum (Milne, populations xx at t), and jext
    real*8, intent(in) :: xx(nxa_nst)
    real*8, intent(out) :: ss(ng), kk(ng), jj(ng)
    real*8 :: wv(nxa_nlv), nstar, nu, phit
    integer :: ig
    call nxa_weights(t, wv)
    ss = 0d0; kk = 0d0; jj = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jj(k) = w*bplanck(nug(k), t2)
    enddo
    phit = nxa_lvphi(1)*nxa_sahac*t**(-1.5d0)
    nstar = nh*xx(nxa_nst)*nh*xx(nxa_nst)*phit
    do ig = 1, nxa_g1(1)
       nu = max(nug(ig), nxa_lvnu0(1))
       kk(ig) = nh*xx(nxa_lvs(1))*wv(1)*nxa_sgm(nxa_sgo(1)+ig)
       if(kk(ig) > 0d0) ss(ig) = 2d0*hh*nu**3/cc**2*nstar*exp(-hh*(nu-nxa_lvnu0(1))/(kb*t))* &
            nxa_sgm(nxa_sgo(1)+ig)/kk(ig)
    enddo
  end subroutine lyc_sk
!
  subroutine lyc_j(xx, jj)
!-- J = jext + Lambda* S of the ground state's Lyman continuum (Lambda* =
!-- 0.999), S from the Milne relation with the populations xx at t
    real*8, intent(in) :: xx(nxa_nst)
    real*8, intent(out) :: jj(ng)
    real*8 :: wv(nxa_nlv), nstar, kp, et, nu, phit
    integer :: ig
    call nxa_weights(t, wv)
    jj = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jj(k) = w*bplanck(nug(k), t2)
    enddo
!-- n1*: Saha with the actual protons and electrons (ne from x)
    phit = nxa_lvphi(1)*nxa_sahac*t**(-1.5d0)
    nstar = nh*xx(nxa_nst)*nh*xx(nxa_nst)*phit
    do ig = 1, nxa_g1(1)
       nu = max(nug(ig), nxa_lvnu0(1))
       kp = nh*xx(nxa_lvs(1))*wv(1)*nxa_sgm(nxa_sgo(1)+ig)
       et = 2d0*hh*nu**3/cc**2*nstar*exp(-hh*(nu-nxa_lvnu0(1))/(kb*t))*nxa_sgm(nxa_sgo(1)+ig)
       if(kp > 0d0) jj(ig) = jj(ig) + 0.999d0*et/kp
    enddo
  end subroutine lyc_j
!
  subroutine set_ali(xx)
!-- ALI data of the toy cell: Lambda* = 0.999, kappa of the ground state
    real*8, intent(in) :: xx(nxa_nst)
    real*8 :: wv(nxa_nlv)
    integer :: ig
    call nxa_weights(t, wv)
    if(allocated(nx_ali_lst)) deallocate(nx_ali_lst, nx_ali_kap)
    allocate(nx_ali_lst(ng), nx_ali_kap(ng))
    nx_ali_lst = 0.999d0
    nx_ali_kap = 0d0
    do ig = 1, nxa_g1(1)
       nx_ali_kap(ig) = nh*xx(nxa_lvs(1))*wv(1)*nxa_sgm(nxa_sgo(1)+ig)
    enddo
    nx_ali_on = .true.
  end subroutine set_ali
!
  subroutine he_tests(dir)
    character(*), intent(in) :: dir
    character(256) :: files(2)
    integer :: ihe, i2, o, k5876, k4471, ndraw, n58, n44, kind, idx, j, s, iel, ns
    real*8 :: yhe, r, rmc, rse, err, wlc, j58, j44, wls(6)
    type(rnd_t) :: st
    files(1) = dir//'/he1.atom'
    files(2) = dir//'/he2.atom'
    call nxt_init(1, 1, ng, wlg, 2, files, nlev)
    allocate(x(nxa_nst), xl(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst))
    allocate(nel(nxa_nel))
    ihe = 0
    do iel = 1, nxa_nel
       if(nxa_el(iel)%z == 2) ihe = iel
    enddo
    write(*,'(a,i3,a,i5,a,i6,a,i5)') '  He model atoms: elements', nxa_nel, '  states', nxa_nst, &
         '  lines', nxa_nlin, '  levels', nxa_nlv
!-- LTE limit at high density
    yhe = 0.1d0
    do j = 1, 2
       t = 1.2d4; if(j == 2) t = 3d4
       nh = 1d17
       do k = 1, ng
          jg(k) = bplanck(nug(k), t)
       enddo
       jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.
       nel = 0d0; nel(1) = nh; nel(ihe) = yhe*nh
       call nx_lte(t, nh, x)
       nti = 0d0; ntx = 0d0
       call nx_solve_cell(t, t, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
       call nx_lte(t, ne, xl)
       o = nxa_el(ihe)%xoff
       dev = 0d0
       do s = 1, nxa_el(ihe)%ns
          if(xl(o+s) > 1d-12) dev = max(dev, abs(x(o+s)/xl(o+s) - 1d0))
       enddo
       write(*,'(a,f7.0,a,i2,a,3es10.3,a,es10.3)') '  He LTE test T=', t, '  ierr=', ierr, '  He I,II,III:', &
            ionfr(o, 1), ionfr(o, 2), x(o+nxa_el(ihe)%ns), '  max |n/n_LTE - 1| =', dev
       if(ierr /= 0 .or. dev > 1d-2) nfail = nfail + 1
    enddo
!-- Case B nebula: He-ionizing diluted 5e4 K blackbody, resonance lines thick
    t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-12
    jg = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nug(k), t2)
    enddo
    jb = 0d0; ljb = .true.
    nel = 0d0; nel(1) = nh; nel(ihe) = yhe*nh
    x = 0d0
    do iel = 1, nxa_nel
       x(nxa_el(iel)%xoff+nxa_el(iel)%ns) = 0.9d0
       x(nxa_el(iel)%xoff+1) = 0.1d0
    enddo
    nti = 0d0; ntx = 0d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
    o = nxa_el(ihe)%xoff
    k5876 = lnear(5875.6d-8); k4471 = lnear(4471.5d-8)
    j58 = popu(k5876)*nxa_la(k5876)*beta(k5876)*nxa_lnu(k5876)
    j44 = popu(k4471)*nxa_la(k4471)*beta(k4471)*nxa_lnu(k4471)
    r = j58/j44
    write(*,'(a,i2,a,3es10.3,a,es10.3)') '  He Case B: ierr=', ierr, '  He I,II,III:', ionfr(o,1), ionfr(o,2), &
         x(o+nxa_el(ihe)%ns), '  ne=', ne
    write(*,'(a,f7.3,a)') '  He I 5876/4471 =', r, '   (Case B at 1e4 K, low density: ~2.7-2.8)'
    wls = [584.33d-8, 10830.3d-8, 3888.6d-8, 5875.6d-8, 4471.5d-8, 6678.2d-8]
    allocate(wv(nxa_nlv))
    call nxa_weights(t, wv)
    do j = 1, 6
       k = lnear(wls(j))
       write(*,'(a,f9.1,a,es10.3,a,es10.3,a,es10.3)') '    He I', nxa_lwl(k)*1d8, ' tau=', nx_kq(k, nel(ihe), &
            x(o+1:o+nxa_el(ihe)%ns), wv(nxa_el(ihe)%lvoff+1:nxa_el(ihe)%lvoff+nxa_el(ihe)%nlv))/q, &
            ' beta=', beta(k), ' x_u=', popu(k)
    enddo
    deallocate(wv)
    if(ierr /= 0 .or. r < 2.4d0 .or. r > 3.2d0) nfail = nfail + 1
!-- macro-atom against the populations (Lucy 2002): exact when the
!-- activations carry all energy input; recombination activations do so only
!-- without collisional excitation (which the thermal pool drives), so the
!-- comparison is made with the electron-impact rates off
    nxa_colfac = 0d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
    j58 = popu(k5876)*nxa_la(k5876)*beta(k5876)*nxa_lnu(k5876)
    j44 = popu(k4471)*nxa_la(k4471)*beta(k4471)*nxa_lnu(k4471)
!-- macro-atom tables of this cell
    nxt_x(:,1) = x; nxt_nel(:,1) = nel; nxt_ne(1) = ne
    nxt_jg(:,1) = jg; nxt_jb(:,1) = jb; nxt_ljb(:,1) = ljb; nxt_have_mc = .true.
    allocate(grd_tempinv(1)); grd_tempinv = 1d0/t
    call nxt_tables(1, [1d0/t], [t2], [0d0], [1d15, 1.001d15], [q*1d15, q*1.001d15])
    do iel = 1, nxa_nel
       ns = nxa_el(iel)%ns
       do s = 1, ns
          j = nxt_ecoff(iel) + (s-1)*(ns+1)
          if(abs(nxt_ecum(j+ns+1,1) - 1d0) > 1d-12 .or. &
               any(nxt_ecum(j+2:j+ns+1,1) < nxt_ecum(j+1:j+ns,1) - 1d-14)) then
             write(*,*) 'bad exit table: element', iel, ' state', s
             nfail = nfail + 1
          endif
       enddo
    enddo
    if(abs(nxt_kcum(nxt_nko,1) - 1d0) > 1d-12) nfail = nfail + 1
!-- recombination activations of He I: macro-atom in the He II ground state
    i2 = nxa_ion(nxa_el(ihe)%ion(2))%soff + 1
    call rnd_init(1, 7); st = rnd_states(1)
    ndraw = 2000000
    n58 = 0; n44 = 0
    do j = 1, ndraw
       call nxt_activate(1, ihe, i2, st, kind, idx, wlc)
       if(kind == em_line .and. idx == k5876) n58 = n58 + 1
       if(kind == em_line .and. idx == k4471) n44 = n44 + 1
    enddo
!-- packets carry equal energy: the count ratio is the energy ratio
    rmc = dble(n58)/dble(max(n44,1))
    err = rmc*sqrt(1d0/max(n58,1) + 1d0/max(n44,1))
    rse = j58/j44
    write(*,'(a,f8.4,a,f7.4,a,f8.4)') '  He I 5876/4471 energy from the macro-atom:', rmc, ' +-', err, &
         '   from the populations:', rse
    if(abs(rmc - rse) > 4d0*err + 0.02d0*rse) nfail = nfail + 1
!-- the same for H-alpha/H-beta (H recombination)
    k5876 = ions_lin(1, 2, 3); k4471 = ions_lin(1, 2, 4)
    n58 = 0; n44 = 0
    do j = 1, ndraw
       call nxt_activate(1, 1, nxa_el(1)%ns, st, kind, idx, wlc)
       if(kind == em_line .and. idx == k5876) n58 = n58 + 1
       if(kind == em_line .and. idx == k4471) n44 = n44 + 1
    enddo
    rmc = dble(n58)/dble(max(n44,1))
    err = rmc*sqrt(1d0/max(n58,1) + 1d0/max(n44,1))
    rse = (x(3)*nxa_la(k5876)*nxt_beta(k5876,1)*nxa_lnu(k5876))/(x(4)*nxa_la(k4471)*nxt_beta(k4471,1)*nxa_lnu(k4471))
    write(*,'(a,f8.4,a,f7.4,a,f8.4)') '  H-alpha/H-beta energy from the macro-atom:', rmc, ' +-', err, &
         '   from the populations:', rse
    if(abs(rmc - rse) > 4d0*err + 0.02d0*rse) nfail = nfail + 1
    nxa_colfac = 1d0
  end subroutine he_tests
!
  integer function ions_lin(iel, l, u)
!-- global index of the line (l,u) of the first ion of element iel
    integer, intent(in) :: iel, l, u
    integer :: ii
    ii = nxa_el(iel)%ion(1)
    ions_lin = nxa_ion(ii)%loff + nxa_ion(ii)%lin(l,u)
  end function ions_lin
!
  subroutine all_lte(dir)
!-- every model atom of the directory (list in atoms.txt, if present): the
!-- populations of all elements tend to Saha-Boltzmann at high density in a
!-- Planck field (detailed balance of the converted atomic data)
    character(*), intent(in) :: dir
    character(256) :: files(80), line
    integer :: nf, ios, iel, s, o, j, nbad, ncx
    real*8 :: devel, dmax
    integer :: cxf(4096), cxt(4096)
    real*8 :: cxr(4096), fnet(2000), fabs(2000)
    nf = 0
    open(9, file=dir//'/atoms.txt', status='old', iostat=ios)
    if(ios /= 0) then
       write(*,'(a)') '  (no atoms.txt in the model-atom directory: all-element LTE test skipped)'
       return
    endif
    do
       read(9,'(a)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       nf = nf + 1
       files(nf) = dir//'/'//trim(line)
    enddo
    close(9)
    call nxt_init(1, 1, ng, wlg, nf, files(:nf), nlev)
    call nxa_read_cx(dir//'/chgexch.dat', .true.)
    if(allocated(x)) deallocate(x, xl, beta, jb, ljb, nti, ntx, nel)
    allocate(x(nxa_nst), xl(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst))
    allocate(nel(nxa_nel))
    write(*,'(a,i3,a,i3,a,i5,a,i7)') '  all atoms:', nxa_nel, ' elements', nxa_nion, ' ions', nxa_nst, &
         ' states', nxa_nlin
    nbad = 0
!-- (1e5 K for the highest ions: He III, C V, N V-VI, O V-VII)
    do j = 1, 3
       t = 8d3; if(j == 2) t = 2.5d4; if(j == 3) t = 1d5
       nh = 1d17
       do k = 1, ng
          jg(k) = bplanck(nug(k), t)
       enddo
       jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.
       do iel = 1, nxa_nel
          nel(iel) = nh*1d-4
       enddo
       nel(1) = nh
       call nx_lte(t, nh, x)
       nti = 0d0; ntx = 0d0
       call nx_solve_cell(t, t, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
       call nx_lte(t, ne, xl)
!-- charge exchange: no net flow into any state at LTE (detailed balance)
       if(nxa_ncx > 0) then
          dmax = 0d0
          allocate(wv(nxa_nlv))
          call nxa_weights(t, wv)
          do iel = 1, nxa_nel
             o = nxa_el(iel)%xoff
             call nx_cx_rates(iel, t, nel, xl, wv, ncx, cxf, cxt, cxr)
             fnet = 0d0; fabs = 0d0
             do k = 1, ncx
                fnet(cxt(k)) = fnet(cxt(k)) + xl(o+cxf(k))*cxr(k)
                fnet(cxf(k)) = fnet(cxf(k)) - xl(o+cxf(k))*cxr(k)
                fabs(cxf(k)) = fabs(cxf(k)) + xl(o+cxf(k))*cxr(k)
             enddo
             do k = 1, nxa_el(iel)%ns
                if(fabs(k) > 0d0) dmax = max(dmax, abs(fnet(k))/fabs(k))
             enddo
          enddo
          deallocate(wv)
          write(*,'(a,f7.0,a,i3,a,es9.2)') '  charge exchange at LTE, T=', t, ':', nxa_ncx, &
               ' reactions, max |net/out| per state =', dmax
          if(dmax > 1d-8) nbad = nbad + 1
       endif
       do iel = 1, nxa_nel
          o = nxa_el(iel)%xoff
          devel = 0d0
          do s = 1, nxa_el(iel)%ns
             if(xl(o+s) > 1d-10) devel = max(devel, abs(x(o+s)/xl(o+s) - 1d0))
          enddo
          if(devel > 2d-2 .or. ierr /= 0) then
             nbad = nbad + 1
             write(*,'(a,f7.0,a,i3,a,i2,a,es10.3)') '    T=', t, '  Z=', nxa_el(iel)%z, '  ierr=', ierr, &
                  '  max |n/n_LTE - 1| =', devel
          endif
       enddo
    enddo
    write(*,'(a,i3)') '  all-element LTE test, elements off LTE:', nbad
    if(nbad > 0) nfail = nfail + 1
!-- the same with level dissolution (nltedismod) at 25 000 K: the dissolved
!-- shares within the super-levels and the states' mean w of the
!-- collisional de-excitations keep every element at its dissolved LTE
    nxd_on = .true.
    call nxd_setup(ng, wlg)
    t = 2.5d4; nh = 1d17
    do k = 1, ng
       jg(k) = bplanck(nug(k), t)
    enddo
    jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.
    nx_dis = 0d0
    call nx_lte(t, nh, x)
    ne = nh
    do j = 1, 3
       call nxd_cellx(ne, nel, x, t, nx_dis)
       nti = 0d0; ntx = 0d0
       call nx_solve_cell(t, t, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
    enddo
    call nxd_cellx(ne, nel, x, t, nx_dis)
    call nx_lte(t, ne, xl)
    nbad = 0
    dmax = 0d0
    do iel = 1, nxa_nel
       o = nxa_el(iel)%xoff
       devel = 0d0
       do s = 1, nxa_el(iel)%ns
          if(xl(o+s) > 1d-10) devel = max(devel, abs(x(o+s)/xl(o+s) - 1d0))
       enddo
       dmax = max(dmax, devel)
       if(devel > 2d-2 .or. ierr /= 0) then
          nbad = nbad + 1
          write(*,'(a,i3,a,i2,a,es10.3)') '    dissolution: Z=', nxa_el(iel)%z, '  ierr=', ierr, &
               '  max |n/n_LTE - 1| =', devel
       endif
    enddo
    write(*,'(a,es10.3,a,es10.3,a,i3)') '  all-element LTE test with dissolution (n_e =', ne, &
         '): max |n/n_LTE - 1| =', dmax, ', elements off LTE:', nbad
    if(nbad > 0) nfail = nfail + 1
    nx_dis = 0d0
    nxd_on = .false.
    call bf_pick_test
  end subroutine all_lte
!
  subroutine bf_pick_test
!-- bound-free absorption: the level choice (nxt_bfpick) against the
!-- weights n_l <sigma_l>_g at a few groups, all atoms at LTE (25 000 K),
!-- with the default lists and with lists of 4 levels (most choices then
!-- walk the other levels), and with level dissolution just below the
!-- Lyman, Balmer and Paschen edges (the dissolved continua, bfsig);
!-- chi-square over the levels expected at least 10 times, the others
!-- lumped
    real*8,allocatable :: p(:)
    integer,allocatable :: cnt(:)
    real*8 :: wls(4), chi2, e, ptot, elump, olump, tb, nhb, neb, pd(3), dm, dh, dhs, dsp
    integer :: k, iw, ig, lv, iel, j, ndr, nbin, kb, lsel, nbad, kbs(3)
    type(rnd_t) :: st
    nbad = 0
    tb = 2.5d4; nhb = 1d17
    do k = 1, ng
       jg(k) = bplanck(nug(k), tb)
    enddo
    jb = [(bplanck(nxa_lnu(k), tb), k=1,nxa_nlin)]; ljb = .true.
    nel = nhb*1d-4; nel(1) = nhb
    call nx_lte(tb, nhb, x)
    nti = 0d0; ntx = 0d0
    call nx_solve_cell(tb, tb, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, neb, beta, ierr)
    nxt_x(:,1) = x; nxt_nel(:,1) = nel; nxt_ne(1) = neb
    nxt_jg(:,1) = jg; nxt_jb(:,1) = jb; nxt_ljb(:,1) = ljb; nxt_have_mc = .true.
    if(.not.allocated(grd_tempinv)) allocate(grd_tempinv(1))
    grd_tempinv = 1d0/tb
    allocate(p(nxa_nlv), cnt(nxa_nlv))
    kbs = [256, 4, 256]
    if(.not.allocated(rnd_states)) call rnd_init(1, 11)
    st = rnd_states(1)
    do kb = 1, 3
       nxt_bfk = kbs(kb)
       wls = [3000d-8, 900d-8, 500d-8, 200d-8]
       pd = 0d0
       if(kb == 3) then
          nxd_on = .true.
          call nxd_setup(ng, wlg)
          call nxd_cellx(neb, nel, x, tb, pd)
          wls = [3700d-8, 915d-8, 8400d-8, 3000d-8]
       endif
       call nxt_tables(1, [1d0/tb], [tb], [0d0], [1d15, 1.001d15], [1d15/1.7d6, 1.001d15/1.7d6])
       ndr = 200000; if(kb == 2) ndr = 50000
       do iw = 1, 4
          ig = minloc(abs(wlg(1:ng) - wls(iw)), 1)
          do lv = 1, nxa_nlv
             p(lv) = 0d0
             iel = nxa_lvel(lv)
             if(ig <= nxa_g1(lv)) p(lv) = nxa_sgm(nxa_sgo(lv)+ig)
             if(kb == 3 .and. nxd_ga(lv) > 0) then
                if(ig >= nxd_ga(lv) .and. ig <= nxd_gb(lv)) then
                   call nxd_gint(lv, ig, pd, 1d4, nxd_numin(nxa_lvnu0(lv), nxd_z(lv), pd), dm, dh, dhs, dsp)
                   p(lv) = p(lv) + dm
                endif
             endif
             p(lv) = nxt_nel(iel,1)*nxt_x(nxa_el(iel)%xoff+nxa_lvs(lv),1)*nxt_w(lv,1)*p(lv)
          enddo
          ptot = sum(p)
          if(ptot <= 0d0) cycle
          p = p/ptot
          cnt = 0
          do j = 1, ndr
             call nxt_bfpick(1, ig, st, lsel)
             if(lsel < 1) then
                nbad = nbad + 1
                exit
             endif
             cnt(lsel) = cnt(lsel) + 1
          enddo
          chi2 = 0d0; nbin = 0; elump = 0d0; olump = 0d0
          do lv = 1, nxa_nlv
             if(p(lv) <= 0d0 .and. cnt(lv) > 0) nbad = nbad + 1
             e = p(lv)*ndr
             if(e >= 10d0) then
                chi2 = chi2 + (cnt(lv) - e)**2/e
                nbin = nbin + 1
             else
                elump = elump + e
                olump = olump + cnt(lv)
             endif
          enddo
          if(elump >= 10d0) then
             chi2 = chi2 + (olump - elump)**2/elump
             nbin = nbin + 1
          endif
          write(*,'(a,i4,a,f7.1,a,i3,a,f8.1,a,i4,a,es9.2)') '  bound-free level choice: list', nxt_bfk, ', ', &
               wlg(ig)*1d8, ' A:', nbin, ' bins, chi2 =', chi2, ' for', max(nbin-1,1), &
               '; mean share of walks', nxt_bfstat(2)/max(nxt_bfstat(1),1d0)
          if(nbin > 1) then
             if(chi2 > nbin - 1 + 5d0*sqrt(2d0*(nbin - 1))) nbad = nbad + 1
          endif
       enddo
    enddo
    nxt_bfk = 256
    nxd_on = .false.
    deallocate(p, cnt)
    if(nbad > 0) then
       write(*,'(a,i4)') '  bound-free level choice: failed checks', nbad
       nfail = nfail + 1
    endif
  end subroutine bf_pick_test
!
  real*8 function ionfr(o, j)
    integer, intent(in) :: o, j
    integer :: ii
    ionfr = 0d0
    do ii = 1, nxa_nion
       if(nxa_el(nxa_ion(ii)%el)%xoff /= o) cycle
       if(nxa_ion(ii)%stage /= j) cycle
       ionfr = sum(x(o+nxa_ion(ii)%soff+1:o+nxa_ion(ii)%soff+nxa_ion(ii)%nsl))
    enddo
  end function ionfr
!
  integer function lnear(wl)
!-- He line closest to wl [cm]
    real*8, intent(in) :: wl
    integer :: kk
    real*8 :: dmin
    dmin = 1d99
    lnear = 0
    do kk = 1, nxa_nlin
       if(nxa_el(nxa_lel(kk))%z /= 2) cycle
       if(abs(nxa_lwl(kk) - wl) < dmin) then
          dmin = abs(nxa_lwl(kk) - wl)
          lnear = kk
       endif
    enddo
  end function lnear
!
  real*8 function alpha_ref(nn, tt)
!-- Milne integral of nlteh_sigma in ln u (trapezoid, 200000 points to y = 60)
    integer, intent(in) :: nn
    real*8, intent(in) :: tt
    integer :: jj
    integer, parameter :: np = 200000
    real*8 :: nu0, ykt, lu, umax, u, y, f, fo, dl, ss
    nu0 = nh_nu0(nn)
    ykt = hh*nu0/(kb*tt)
    umax = 1d0 + 60d0/ykt
    dl = log(umax)/np
    ss = 0d0
    fo = 0d0
    do jj = 0, np
       lu = jj*dl
       u = exp(lu)
       y = (u - 1d0)*ykt
       f = nlteh_sigma(nn, u*nu0)*(u*nu0)**2*exp(-y)*u*nu0     !d nu = nu d ln u
       if(jj > 0) ss = ss + 0.5d0*(f + fo)*dl
       fo = f
    enddo
    alpha_ref = 4d0*pi_*dble(nn*nn)*4.1412d-16*tt**(-1.5d0)*2d0/cc**2*ss
  end function alpha_ref
!
  real*8 function popu(kk)
!-- population (fraction of the element) of the upper level of line kk
    integer, intent(in) :: kk
    popu = x(nxa_el(nxa_lel(kk))%xoff + nxa_lsu(kk))
  end function popu
end program test_nltex
