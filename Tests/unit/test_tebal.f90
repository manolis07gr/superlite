!-- Unit test for tebalmod (electron temperature from the thermal balance).
!--  1. A Planck radiation field J = B(T0) with LTE opacities and NLTE H must
!--     give T = T0 (detailed balance of every heating/cooling pair: LTE
!--     opacities, H photoionization/recombination, H collisions, Compton with
!--     induced scattering), from above and below, at low and high density.
!--  2. Pure H photoionized by a dilute 4e4 K blackbody at low density, no other
!--     coolants: T_e is set by photoionization heating against recombination
!--     and collisional cooling; printed with its terms (expected ~1-3e4 K).
!--  3. Fast-electron energy shares (Shull & van Steenberg 1985): heat, H
!--     ionization and excitation sum to 1, heat grows and ionization and
!--     excitation fall with the ionized fraction, all heat when fully ionized,
!--     and the published neutral-gas values at x = 1e-4.
program test_tebal
  use nltehmod
  use tebalmod
  implicit none
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16, pi = 3.141592653589793d0
  integer,parameter :: nlev = 30, ng = 3000
  real*8 :: wlg(ng+1), nug(ng), dnug(ng), jg(ng), jb(1000), x(0:nlev), tnew, nh, t0, ne, beta(1000), terms(6)
  real*4 :: capth(ng)
  logical :: ljb(1000)
  real*8, parameter :: t0s(2) = [7d3, 2d4], nhs(2) = [1d8, 1d13], tst(2) = [0.6d0, 1.7d0]
  integer :: i, k, j, m, ierr, nfail
  real*8 :: f(3), fo(3)
  nfail = 0
  call nlteh_init(nlev)
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i = 1, ng
     nug(i) = 0.5d0*(cc/wlg(i) + cc/wlg(i+1))
     dnug(i) = cc/wlg(i) - cc/wlg(i+1)
     capth(i) = real(1d-14*(nug(i)/1d15)**(-2), 4)      !free-free-like LTE opacity
  enddo
  call te_init(ng, nug)
  te_damp = 1d0
  te_maxfac = 100d0
  jb = 0d0
  ljb = .false.
!-- 1. Planck field
  do i = 1, 2
     t0 = t0s(i)
     do k = 1, ng
        jg(k) = bnu(nug(k), t0)
     enddo
     do j = 1, 2
        nh = nhs(j)
        call nlteh_lte(t0, nh, x)
        do m = 1, 2
           call te_balance(t0*tst(m), t0, nh, 0d0, 1d-7, 1d-7, ng, nug, dnug, jg, jb, ljb, 0d0, 0d0, capth, &
                0d0, x, tnew, ierr)
           write(6,'(a,f8.0,a,es8.1,a,f8.0,a,f8.1,a,i2)') 'Planck T0=', t0, ' nH=', nh, '  start ', t0*tst(m), &
                ' -> T=', tnew, '  ierr', ierr
           if(abs(tnew/t0 - 1d0) > 5d-3 .or. ierr /= 0) nfail = nfail + 1
        enddo
     enddo
  enddo
!-- 2. dilute 4e4 K blackbody, W = 1e-14 (a hot star at 0.1-1 pc), n_H = 100
  do k = 1, ng
     jg(k) = 1d-14*bnu(nug(k), 4d4)
  enddo
  capth = 0.
  nh = 1d2
  x = 0d0; x(0) = 0.99d0; x(1) = 0.01d0
  call te_balance(1d4, 4d4, nh, 0d0, 1d-15, 1d-15, ng, nug, dnug, jg, jb, ljb, 0d0, 0d0, capth, 0d0, x, tnew, ierr)
  write(6,'(a,f9.0,a,i2)') 'pure H nebula, 4e4 K star: T_e =', tnew, '  ierr', ierr
  call nlteh_solve_cell(tnew, 4d4, nh, 0d0, 1d-15, 1d-15, ng, nug, dnug, jg, jb, ljb, 0d0, x, ne, beta, ierr)
  call te_terms(tnew, ng, nug, dnug, jg, capth, nh, x, ne, 0d0, terms)
  write(6,'(a,6es11.3)') '  terms (LTE, PI heat, rec cool, coll, Compton, ext):', terms
  write(6,'(a,es10.3)') '  x(H II) =', x(0)
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
  if(nfail == 0) then
     write(6,*) 'test_tebal: PASS'
  else
     write(6,*) 'test_tebal: FAIL', nfail
     stop 1
  endif
contains
  real*8 function bnu(nu, t)
    real*8, intent(in) :: nu, t
    real*8 :: e
    e = hh*nu/(kb*t)
    if(e > 700d0) then
       bnu = 0d0
    else
       bnu = 2d0*hh*nu**3/cc**2/(exp(e) - 1d0)
    endif
  end function bnu
end program test_tebal
