!-- Unit tests for nltehmod (n-resolved NLTE hydrogen)
!--  1. atomic data against exact hydrogenic values
!--  2. total recombination coefficients against Hui & Gnedin (1997)
!--  3. J = B(T) everywhere must give exact Saha-Boltzmann populations at any
!--     density (detailed balance of every rate pair)
!--  4. Case B: photoionized low-density gas with optically thick Lyman lines
!--     must give the Storey & Hummer (1995) Balmer decrement, 2.86, to 5%:
!--     the l-bundled atom is exact only where collisions mix the l states
!--     (n_e >~ 1e8 cm^-3 for the levels that matter); at 1e4 cm^-3 bundling
!--     lowers H-alpha/H-beta by ~4% (2.74 with 20 or 30 levels)
program test_nlteh
  use nltehmod
  implicit none
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16, pi = 3.141592653589793d0
  integer,parameter :: nlev = 30, ng = 6000
  real*8 :: wlg(ng+1), nug(ng), dnug(ng), jg(ng), jb(1000), x(0:nlev), xl(0:nlev), ne, beta(1000)
  real*8 :: t, nh, aA, aB, hg_a, hg_b, lam, r, dev, w, t2, q
  real*8, parameter :: tr(3) = [5d3, 1d4, 2d4], tl(3) = [6d3, 1d4, 3d4], nl(3) = [1d6, 1d10, 1d14]
  logical :: ljb(1000)
  integer :: i, k, ierr, nfail, n, kha, khb
  nfail = 0
  call nlteh_init(nlev)
!-- 1. atomic data
  call check('f(1-2)', nh_f(1), 0.4162d0, 0.005d0)
  k = minloc(abs(nh_wl - 6563d-8), 1)
  call check('f(2-3)', nh_f(k), 0.6407d0, 0.005d0)
  call check('A(3-2) [1/s]', nh_a(k), 4.410d7, 0.01d0)
  call check('lambda(H-alpha, air) [A]', nh_wl(k)*1d8, 6562.80d0, 2d-5)
  call check('sigma_1 at threshold [cm^2]', nlteh_sigma(1, nh_nu0(1)), 6.30d-18, 0.01d0)
  call check('sigma_2 at threshold [cm^2]', nlteh_sigma(2, nh_nu0(2)), 1.39d-17, 0.02d0)
!-- 2. recombination: Hui & Gnedin (1997) case A/B fits
  do i = 1, 3
     t = tr(i)
     lam = 2d0*157807d0/t
     hg_a = 1.269d-13*lam**1.503d0/(1d0 + (lam/0.522d0)**0.470d0)**1.923d0
     hg_b = 2.753d-14*lam**1.5d0/(1d0 + (lam/2.740d0)**0.407d0)**2.242d0
     aA = 0d0
     do n = 1, nlev
        aA = aA + nlteh_alpha(n, t)
     enddo
     aB = aA - nlteh_alpha(1, t)
     write(*,'(a,f7.0,a,2es11.3,a,2es11.3)') '  T=', t, '  alpha_A, alpha_B (n<=30):', aA, aB, &
          '   Hui&Gnedin:', hg_a, hg_b
!-- levels above n=30 add a few per cent; Johnson Gaunt factors are good to ~2%
     call check('alpha_B(n<=30)/HG', aB/hg_b, 0.96d0, 0.04d0)
  enddo
!-- group grid as in the A4 test (1 A - 30000 A, logarithmic)
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i = 1, ng
     nug(i) = 0.5d0*(cc/wlg(i) + cc/wlg(i+1))
     dnug(i) = cc/wlg(i) - cc/wlg(i+1)
  enddo
!-- 3. LTE limit at low and high density
  do i = 1, 3
     t = tl(i)
     nh = nl(i)
     do k = 1, ng
        jg(k) = bnu(nug(k), t)
     enddo
     do k = 1, nh_nlin
        jb(k) = bnu(nh_nu(k), t)
     enddo
     ljb = .true.
     x = 0d0; x(0) = 0.5d0; x(1) = 0.5d0
     call nlteh_solve_cell(t, t, nh, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, 0d0, x, ne, beta, ierr)
     call nlteh_lte(t, ne, xl)
     dev = maxval(abs(x(0:nlev)/xl(0:nlev) - 1d0))
     write(*,'(a,es9.2,a,f7.0,a,i2,a,es10.3,a,es10.3)') '  LTE test nH=', nh, ' T=', t, '  ierr=', ierr, &
          '  x(HII)=', x(0), '  max |n/n_LTE - 1| =', dev
     if(ierr /= 0 .or. dev > 0.005d0) nfail = nfail + 1
  enddo
!-- 4. Case B: nH = 1e4, T = 1e4 K, ionizing diluted 5e4 K blackbody only, near-static (Lyman thick)
  t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4
  jg = 0d0
  do k = 1, ng
     if(nug(k) > nh_nu0(1)) jg(k) = w*bnu(nug(k), t2)
  enddo
  jb = 0d0; ljb = .true.
  x = 0d0; x(0) = 0.99d0; x(1) = 0.01d0
  q = 1d-9       !velocity gradient [1/s]: Lyman lines very thick, Balmer thin
  call nlteh_solve_cell(t, t2, nh, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, 0d0, x, ne, beta, ierr)
  kha = minloc(abs(nh_wl - 6563d-8), 1); khb = minloc(abs(nh_wl - 4861d-8), 1)
  r = (x(3)*nh_a(kha)*beta(kha)/nh_wl(kha))/(x(4)*nh_a(khb)*beta(khb)/nh_wl(khb))
  write(*,'(a,i2,a,f8.5,a,es9.2,a,es9.2,a,2f8.5)') '  Case B test: ierr=', ierr, '  x(HII)=', x(0), &
       '  beta(Ly-a)=', beta(1), '  ne=', ne, '  beta(Ha,Hb)=', beta(kha), beta(khb)
  if(beta(1) > 1d-3 .or. beta(kha) < 0.99d0) write(*,'(a)') '  (not Case B conditions)'
  call check('Case B H-alpha/H-beta (Storey & Hummer 2.86)', r, 2.86d0, 0.05d0)
  if(ierr /= 0) nfail = nfail + 1
!
  if(nfail > 0) then
     write(*,'(a,i3)') 'test_nlteh: FAILED checks:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_nlteh: PASS'
contains
  function bnu(nu, tt) result(b)
    real*8, intent(in) :: nu, tt
    real*8 :: b, xx
    xx = hh*nu/(kb*tt)
    if(xx > 700d0) then
       b = 0d0
    else
       b = 2d0*hh*nu**3/cc**2/(exp(xx) - 1d0)
    endif
  end function bnu
  subroutine check(name, val, ref, tol)
    character(*), intent(in) :: name
    real*8, intent(in) :: val, ref, tol
    logical :: ok
    ok = abs(val/ref - 1d0) <= tol
    write(*,'(2x,a,a,es12.5,a,es12.5,a,l1)') name, ': ', val, '  expected ', ref, '  ok=', ok
    if(.not. ok) nfail = nfail + 1
  end subroutine check
end program test_nlteh
