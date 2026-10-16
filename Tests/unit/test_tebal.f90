!-- Unit test for tebalmod (electron temperature from the thermal balance).
!--  1. A Planck radiation field J = B(T0) with LTE opacities and the NLTE
!--     elements must give T = T0 (detailed balance of every heating/cooling
!--     pair: LTE opacities, photoionization/recombination, collisions,
!--     Compton with induced scattering), from above and below, at low and
!--     high density: built-in H, and H with He I and He II if a model-atom
!--     directory is given as the first argument.
!--  2. Pure H photoionized by a dilute 4e4 K blackbody at low density, no other
!--     coolants: T_e is set by photoionization heating against recombination
!--     and collisional cooling; printed with its terms (expected ~1-3e4 K).
!--  3. Fast-electron energy shares (Shull & van Steenberg 1985): heat, H
!--     ionization and excitation sum to 1, heat grows and ionization and
!--     excitation fall with the ionized fraction, all heat when fully ionized,
!--     the published neutral-gas values at x = 1e-4, and the He shares.
program test_tebal
  use nltehmod, only:nlteh_init,nlteh_svs,bplanck
  use nlteatommod
  use nltexmod
  use nltextrmod, only:nxt_shares,nxt_init,nxt_elh,nxt_elhe
  use tebalmod
  implicit none
  real*8,parameter :: cc = 2.99792458d10
  integer,parameter :: nlev = 30, ng = 3000
  real*8 :: wlg(ng+1), nug(ng), dnug(ng), jg(ng), tnew, nh, t0, ne, terms(te_nterm)
  real*4 :: capth(ng)
  real*8,allocatable :: x(:), jb(:), beta(:), nti(:), ntx(:), nel(:)
  logical,allocatable :: ljb(:)
  real*8, parameter :: t0s(2) = [7d3, 2d4], nhs(2) = [1d8, 1d13], tst(2) = [0.6d0, 1.7d0]
  integer :: i, k, j, m, ierr, nfail, ipass
  real*8 :: f(3), fo(3), g(5)
  character(256) :: adir, files(2)
  type(nxion) :: ions(1)
  nfail = 0
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i = 1, ng
     nug(i) = 0.5d0*(cc/wlg(i) + cc/wlg(i+1))
     dnug(i) = cc/wlg(i) - cc/wlg(i+1)
     capth(i) = real(1d-14*(nug(i)/1d15)**(-2), 4)      !free-free-like LTE opacity
  enddo
  te_damp = 1d0
  te_maxfac = 100d0
  do ipass = 1, 2
     if(ipass == 1) then
!-- built-in hydrogen
        call nlteh_init(nlev)
        call nxa_builtin_h(nlev, ions(1))
        call nxa_setup(1, ions)
        call nxa_groups(ng, wlg)
     else
        if(command_argument_count() < 1) then
           write(6,'(a)') '  (no model-atom directory given: He test skipped)'
           exit
        endif
        call get_command_argument(1, adir)
        files(1) = trim(adir)//'/he1.atom'
        files(2) = trim(adir)//'/he2.atom'
        call nxt_init(1, 1, ng, wlg, 2, files, nlev)
        write(6,'(a,i3,a)') ' with He:', nxa_nion, ' ions'
     endif
     if(allocated(x)) deallocate(x, jb, ljb, beta, nti, ntx, nel)
     allocate(x(nxa_nst), jb(nxa_nlin), ljb(nxa_nlin), beta(nxa_nlin), nti(nxa_nst), ntx(nxa_nst), nel(nxa_nel))
     jb = 0d0; ljb = .false.; nti = 0d0; ntx = 0d0
!-- 1. Planck field
     do i = 1, 2
        t0 = t0s(i)
        do k = 1, ng
           jg(k) = bplanck(nug(k), t0)
        enddo
        do j = 1, 2
           nh = nhs(j)
           nel = 0d0
           nel(1) = nh
           if(nxa_nel > 1) nel(2) = 0.1d0*nh
           call nx_lte(t0, nh, x)
           do m = 1, 2
              call te_balance(t0*tst(m), t0, nel, 0d0, 1d-7, 1d-7, ng, nug, dnug, jg, jb, ljb, nti, ntx, &
                   capth, capth, 0d0, x, tnew, ierr)
              write(6,'(a,f8.0,a,es8.1,a,f8.0,a,f8.1,a,i2)') 'Planck T0=', t0, ' nH=', nh, '  start ', t0*tst(m), &
                   ' -> T=', tnew, '  ierr', ierr
              if(abs(tnew/t0 - 1d0) > 5d-3 .or. ierr /= 0) nfail = nfail + 1
           enddo
        enddo
     enddo
  enddo
!-- 2. dilute 4e4 K blackbody, W = 1e-14 (a hot star at 0.1-1 pc), n_H = 100,
!-- hydrogen only
  call nxa_builtin_h(nlev, ions(1))
  call nxa_setup(1, ions)
  call nxa_groups(ng, wlg)
  deallocate(x, jb, ljb, beta, nti, ntx, nel)
  allocate(x(nxa_nst), jb(nxa_nlin), ljb(nxa_nlin), beta(nxa_nlin), nti(nxa_nst), ntx(nxa_nst), nel(nxa_nel))
  jb = 0d0; ljb = .false.; nti = 0d0; ntx = 0d0
  do k = 1, ng
     jg(k) = 1d-14*bplanck(nug(k), 4d4)
  enddo
  capth = 0.
  nh = 1d2
  nel = nh
  x = 0d0; x(nxa_nst) = 0.99d0; x(1) = 0.01d0
  call te_balance(1d4, 4d4, nel, 0d0, 1d-15, 1d-15, ng, nug, dnug, jg, jb, ljb, nti, ntx, capth, capth, 0d0, &
       x, tnew, ierr)
  write(6,'(a,f9.0,a,i2)') 'pure H nebula, 4e4 K star: T_e =', tnew, '  ierr', ierr
  call nx_solve_cell(tnew, 4d4, nel, 0d0, 1d-15, 1d-15, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
  call te_terms(tnew, ng, nug, dnug, jg, capth, capth, nel, x, ne, 0d0, terms)
  write(6,'(a,6es11.3)') '  terms (LTE, PI heat, rec cool, coll, Compton, ext):', terms
  write(6,'(a,es10.3)') '  x(H II) =', x(nxa_nst)
  if(tnew < 1d4 .or. tnew > 4d4 .or. ierr /= 0) nfail = nfail + 1
!-- 3. fast-electron shares
  fo = [0d0, 1d0, 1d0]
  do k = -4, 0
     call nlteh_svs(10d0**k, f(1), f(2), f(3))
     write(6,'(a,es8.1,a,3f7.3)') 'SvS85 x =', 10d0**k, '  heat, ion, exc:', f
     if(abs(sum(f) - 1d0) > 1d-12) nfail = nfail + 1
     if(k > -4 .and. (f(1) < fo(1) .or. f(2) > fo(2) .or. f(3) > fo(3))) nfail = nfail + 1
     fo = f
  enddo
  if(abs(f(1) - 1d0) > 1d-12) nfail = nfail + 1
  call nlteh_svs(1d-4, f(1), f(2), f(3))
  if(abs(f(2) - 0.375d0) > 0.01d0 .or. abs(f(3) - 0.419d0) > 0.01d0) nfail = nfail + 1
!-- with He at the SvS abundance: shares sum to 1, He ionization 0.0554 at x -> 0
  call nxt_shares(1d-4, 1d0, g(1), g(2), g(3), g(4), g(5))
  write(6,'(a,5f7.3)') 'SvS85 with He, x = 1e-4: heat, H ion, H exc, He ion, He exc:', g
  if(abs(sum(g) - 1d0) > 1d-12 .or. abs(g(4) - 0.0554d0) > 0.003d0) nfail = nfail + 1
  if(nfail == 0) then
     write(6,*) 'test_tebal: PASS'
  else
     write(6,*) 'test_tebal: FAIL', nfail
     stop 1
  endif
  if(nxt_elh < 0 .or. nxt_elhe < 0) stop 'unreachable'
end program test_tebal
