!-- Unit test for bbsobmod (Sobolev expansion opacity and line thermalization).
!--   * F(K) = <|Q|(1-exp(-K/|Q|))> against a brute-force angle integral for
!--     random accelerating, decelerating and mixed gradients;
!--   * thin limit F = K (smeared opacity), thick limit F = <|Q|>;
!--   * homologous flow: F = (1-exp(-K t))/t exactly (Eastman & Pinto 1993);
!--   * van Regemorter P(y): table values, continuity at the table ends;
!--   * eps for H-alpha at 1e4 K, n_e = 1e10: C_32/A_32 = 20.6 lam^3 n_e P/sqrt(T).
program test_bbsob
  use bbsobmod
  implicit none
  real*8 :: q(bbs_nq), w(bbs_nq), qavg, qmn, qmx, a, b, kk, f, fref, err, maxerr, t, x, y, p1, p2, eps
  integer :: i, k, nfail
  real*8 :: u
  nfail = 0
!-- brute force over random gradients
  maxerr = 0d0
  do i = 1, 2000
     call random_number(u); a = (u - 0.3d0)*2d-7            ! dv/dr, both signs
     call random_number(u); b = 1d-9 + u*1d-7                ! v/r > 0
     call random_number(u); kk = 10d0**(-12d0 + 8d0*u)       ! K from thin to thick
     call bbs_quad(a, b, 0d0, q, w, qavg, qmn, qmx)
     f = bbs_fq(kk, q, w, qavg, qmn, qmx)
     fref = brute(kk, a, b)
     err = abs(f/fref - 1d0)
     maxerr = max(maxerr, err)
  enddo
  write(6,'(a,es10.3)') 'F(K) vs brute-force angle integral, max rel. error: ', maxerr
  if(maxerr > 2d-2) nfail = nfail + 1
!-- limits
  a = 3d-8; b = 1d-8
  call bbs_quad(a, b, 0d0, q, w, qavg, qmn, qmx)
  f = bbs_fq(1d-20, q, w, qavg, qmn, qmx)
  write(6,'(a,es10.3)') 'thin limit F/K - 1: ', f/1d-20 - 1d0
  if(abs(f/1d-20 - 1d0) > 1d-6) nfail = nfail + 1
  f = bbs_fq(1d0, q, w, qavg, qmn, qmx)
  write(6,'(a,es10.3)') 'thick limit F/<|Q|> - 1: ', f/qavg - 1d0
  if(abs(f/qavg - 1d0) > 1d-12) nfail = nfail + 1
!-- homologous
  t = 20d0*86400d0
  maxerr = 0d0
  do k = -4, 4
     kk = 10d0**k/t
     call bbs_quad(1d0/t, 1d0/t, 0d0, q, w, qavg, qmn, qmx)
     f = bbs_fq(kk, q, w, qavg, qmn, qmx)
     fref = (1d0 - exp(-kk*t))/t
     maxerr = max(maxerr, abs(f/fref - 1d0))
  enddo
  write(6,'(a,es10.3)') 'homologous F vs (1-exp(-K t))/t, max rel. error: ', maxerr
  if(maxerr > 1d-10) nfail = nfail + 1
!-- van Regemorter P(y)
  p1 = bbs_pvr(1d0, .false.); p2 = bbs_pvr(1d0, .true.)
  write(6,'(a,2f8.3)') 'P(1) neutral, ion (0.100, 0.214): ', p1, p2
  if(abs(p1 - .100d0) > 1d-6 .or. abs(p2 - .214d0) > 1d-6) nfail = nfail + 1
  write(6,'(a,4f8.4)') 'P continuity at y = 0.01, 10 (neutral): ', bbs_pvr(.0099999d0, .false.), &
       bbs_pvr(.0100001d0, .false.), bbs_pvr(9.99999d0, .false.), bbs_pvr(10.00001d0, .false.)
  if(abs(bbs_pvr(.0099999d0,.false.) - bbs_pvr(.0100001d0,.false.)) > 1d-4) nfail = nfail + 1
  if(abs(bbs_pvr(9.99999d0,.false.) - bbs_pvr(10.00001d0,.false.)) > 1d-4) nfail = nfail + 1
!-- eps for H-alpha
  x = 6.5628d-5; y = 6.6260693d-27*2.99792458d10/(x*1.3806505d-16*1d4)
  p1 = bbs_pvr(y, .false.)
  eps = bbs_epsl(x, 1d10, 1d4, 1, 1d0)
  t = 20.6d0*x**3*1d10*p1/1d2               ! X = C_32/A_32
  write(6,'(a,es10.3,a,es10.3,a,f6.3)') 'H-alpha 1e4 K, n_e 1e10: X = ', t, '  eps(beta=1) = ', eps, '  P = ', p1
  if(abs(eps - t/(t + 1d0)) > 1d-12) nfail = nfail + 1
!-- C_32 = X*A_32 is of the order of the hydrogenic H-alpha rate (~1e-7 cm^3/s)
  write(6,'(a,es10.3,a)') '  -> C_32/n_e = ', t*4.41d7/1d10, ' cm^3/s (hydrogenic H-alpha rates at 1e4 K: ~1e-7)'
  if(nfail == 0) then
     write(6,*) 'test_bbsob: PASS'
  else
     write(6,*) 'test_bbsob: FAIL', nfail
     stop 1
  endif
contains
  real*8 function brute(kk, a, b)
    real*8, intent(in) :: kk, a, b
    integer, parameter :: m = 200000
    integer :: j
    real*8 :: mu, qq, s
    s = 0d0
    do j = 1, m
       mu = (j - .5d0)/m
       qq = max(abs(mu*mu*a + (1d0 - mu*mu)*b), 1d-3*max(abs(a), abs(b)))
       s = s + qq*(1d0 - exp(-kk/qq))
    enddo
    brute = s/m
  end function brute
end program test_bbsob
