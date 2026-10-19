!-- Electron scattering of line photons trapped in a Sobolev resonance
!-- (nltexmod: nx_ytrap, nx_es_table, nx_bes, nx_wwin, nx_bew):
!--  1. y(tau) = tau G(tau) against a direct quadrature of the double
!--     integral of nx_ytrap's header at several depths
!--  2. nx_bes in a homologous cell (Q independent of the direction) against
!--     (1-beta)(1-exp(-sigma v_D y(K/q)/(K (1-beta)) f)), f = exp(-K wq/q) the
!--     escape from the Lorentz wings (f = 1 and f < 1)
!--  3. a Monte Carlo of one resonance (one direction; complete
!--     redistribution over the Doppler profile, frequency drift, electron
!--     scattering): the share of the photons scattered by an electron on a
!--     flight that the line would have ended, before they escape, against
!--     bes/(beta+bes), from the Sobolev regime (sigma L = 0.03) to the one
!--     where electrons scatter a photon before it drifts a Doppler width
!--     (sigma L = 1, 10)
!--  4. the window of the transport's far-side absorption, W(t) =
!--     int_{-5}^{5} (1 - exp(-t erfc(-x)/2)) dx (nx_wwin), against a
!--     quadrature, below, inside and above its table, and the expectation
!--     nx_bew = sigma v_D W(K/q)/K in a homologous cell
program test_estrap
  use nltexmod
  use randommod
  implicit none
  real*8,parameter :: pi = 3.141592653589793d0
  real*8 :: taus(5), yref, y, dev, q, kq, sv, beta, bes, ref, etab(0:nx_nbt), btab(0:nx_nbt)
  real*8 :: ewtab(0:nx_nbt), ts(6), w
  real*8 :: tau, sl, p, pmc
  integer :: i, nfail
  type(rnd_t) :: st
  nfail = 0
!
!-- 1. the table against the quadrature
  taus = [0.1d0, 3d0, 100d0, 1d4, 1d9]
  dev = 0d0
  do i = 1, 5
     yref = taus(i)*gquad(taus(i))
     y = nx_ytrap(taus(i))
     dev = max(dev, abs(y/yref - 1d0))
     write(6,'(a,es9.2,a,f9.5,a,f9.5)') ' tau', taus(i), ': y table', y, ', quadrature', yref
  enddo
  if(dev > 1d-2) then
     write(6,*) 'FAIL: y(tau) table, max rel. dev.', dev
     nfail = nfail + 1
  endif
!
!-- 2. homologous cell: <y(K/|Q|)> = y(K/q)
  q = 1d-9
  call nx_beta_table(q, q, btab)
  call nx_es_table(q, q, etab)
  dev = 0d0
  do i = 1, 5
     kq = taus(i)*q
     sv = 0.05d0*q
     beta = nx_beta(kq, q, q, btab)
     bes = nx_bes(kq, sv, 0d0, beta, q, q, etab)
     ref = (1d0-beta)*(1d0 - exp(-sv*nx_ytrap(taus(i))/(kq*(1d0-beta))))
     dev = max(dev, abs(bes/ref - 1d0))
     bes = nx_bes(kq, sv, 0.3d0/taus(i), beta, q, q, etab)
     ref = (1d0-beta)*(1d0 - exp(-sv*nx_ytrap(taus(i))/(kq*(1d0-beta))*exp(-0.3d0)))
     dev = max(dev, abs(bes/ref - 1d0))
  enddo
  write(6,'(a,es10.3)') ' nx_bes, homologous cell: max rel. dev.', dev
  if(dev > 1d-2) then
     write(6,*) 'FAIL: nx_bes'
     nfail = nfail + 1
  endif
!
!-- 3. Monte Carlo of one resonance
  call rnd_init(1, 11)
  st = rnd_states(1)
  do i = 1, 4
     tau = merge(30d0, 1d3, i <= 3)
     sl = merge(0.03d0, merge(1d0, 10d0, i == 2 .or. i == 4), i == 1)
     pmc = mc_trap(tau, sl, 40000)
     beta = (1d0 - exp(-tau))/tau
     bes = sl*nx_ytrap(tau)/tau
     bes = (1d0-beta)*(1d0 - exp(-bes/(1d0-beta)))
     p = bes/(beta + bes)
     write(6,'(a,es8.1,a,f6.2,a,f7.4,a,f7.4)') ' trapped e.s. before escape, tau', tau, ' sigma L', sl, &
          ': Monte Carlo', pmc, ', bes/(beta+bes)', p
     if(abs(pmc/p - 1d0) > 0.08d0) then
        write(6,*) 'FAIL: Monte Carlo of the resonance'
        nfail = nfail + 1
     endif
  enddo
!
!-- 4. the window W(t) and the expectation nx_bew
  ts = [3d-3, 0.5d0, 2d2, 1d8, 1d11, 1d13]
  dev = 0d0
  do i = 1, 6
     w = nx_wwin(ts(i))
     ref = wquad(ts(i))
     dev = max(dev, abs(w/ref - 1d0))
     write(6,'(a,es9.2,a,f9.5,a,f9.5)') ' t', ts(i), ': W', w, ', quadrature', ref
  enddo
  call nx_es_table(q, q, etab, ewtab)
  do i = 2, 5
     kq = ts(i)*q
     sv = 0.05d0*q
     dev = max(dev, abs(nx_bew(kq, sv, q, q, ewtab)/(sv*wquad(ts(i))/kq) - 1d0))
  enddo
  write(6,'(a,es10.3)') ' window of the far-side absorption and nx_bew: max rel. dev.', dev
  if(dev > 1d-2) then
     write(6,*) 'FAIL: nx_wwin or nx_bew'
     nfail = nfail + 1
  endif
!
  if(nfail == 0) then
     write(6,*) 'test_estrap: PASS'
  else
     write(6,*) 'test_estrap: FAIL', nfail
     stop 1
  endif
contains
  real*8 function phi(x)
    real*8,intent(in) :: x
    phi = exp(-x*x)/sqrt(pi)
  end function phi
  real*8 function cphi(x)
!-- Phi(x) = int_{-inf}^x phi
    real*8,intent(in) :: x
    cphi = 0.5d0*erfc(-x)
  end function cphi
  real*8 function gquad(t)
!-- G(t) by quadrature: x0 on a uniform grid, the inner integral on a
!-- logarithmic grid in d = x0 - x
    real*8,intent(in) :: t
    integer,parameter :: nx = 1401, nd = 1200
    real*8 :: x0, dd(nd), f(nd), s, ph0
    integer :: i, j
    dd(1) = 0d0
    do j = 2, nd
       dd(j) = 10d0**(-14d0 + 15.6d0*(j-2)/dble(nd-2))
    enddo
    gquad = 0d0
    do i = 1, nx
       x0 = -7d0 + 14d0*(i-1)/dble(nx-1)
       ph0 = phi(x0)
       do j = 1, nd
          f(j) = exp(-t*0.5d0*(erfc(-x0) - erfc(-(x0-dd(j)))))*(1d0 - exp(-t*cphi(x0-dd(j))))
       enddo
       s = 0d0
       do j = 2, nd
          s = s + 0.5d0*(f(j)+f(j-1))*(dd(j)-dd(j-1))
       enddo
       gquad = gquad + merge(0.5d0, 1d0, i==1 .or. i==nx)*ph0*s
    enddo
    gquad = gquad*14d0/dble(nx-1)
  end function gquad
  real*8 function wquad(t)
!-- W(t) by the trapezoidal rule
    real*8,intent(in) :: t
    integer,parameter :: nx = 20001
    real*8 :: x, f
    integer :: i
    wquad = 0d0
    do i = 1, nx
       x = -5d0 + 10d0*(i-1)/dble(nx-1)
       f = 1d0 - exp(-0.5d0*t*erfc(-x))
       wquad = wquad + merge(0.5d0, 1d0, i==1 .or. i==nx)*f
    enddo
    wquad = wquad*10d0/dble(nx-1)
  end function wquad
  real*8 function mc_trap(t, sl, n)
!-- photons emitted over the profile; flights end at a line absorption
!-- (re-emission over the profile) or, if the line cannot absorb the photon
!-- any more, in escape; an electron scattering (distance -ln r/sl in
!-- Sobolev lengths) before the line absorption counts as trapped
    real*8,intent(in) :: t, sl
    integer,intent(in) :: n
    integer :: i, ntr, k
    real*8 :: x0, r, tgt, xa, xb, xm, sa, sc
    ntr = 0
    do i = 1, n
       x0 = gauss()
       do
          call rnd_r(r, st)
          tgt = cphi(x0) + log(max(r,1d-300))/t
          if(tgt <= 0d0) exit                          !escapes the line
!-- the absorption point: Phi(x0 - sa) = tgt (bisection)
          xa = -12d0
          xb = x0
          do k = 1, 60
             xm = 0.5d0*(xa + xb)
             if(cphi(xm) < tgt) then
                xa = xm
             else
                xb = xm
             endif
          enddo
          sa = x0 - 0.5d0*(xa + xb)
          call rnd_r(r, st)
          sc = -log(max(r,1d-300))/sl
          if(sc < sa) then
             ntr = ntr + 1
             exit
          endif
          x0 = gauss()
       enddo
    enddo
    mc_trap = dble(ntr)/dble(n)
  end function mc_trap
  real*8 function gauss()
!-- deviate of the Doppler profile exp(-x^2)/sqrt(pi) (Box-Muller)
    real*8 :: u1, u2
    call rnd_r(u1, st)
    call rnd_r(u2, st)
    gauss = sqrt(-log(max(u1,1d-300)))*cos(2d0*pi*u2)
  end function gauss
end program test_estrap
