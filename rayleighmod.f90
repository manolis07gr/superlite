!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module rayleighmod
!-----------------------------------------------------------------------
! Rayleigh scattering by H I in the ground state, in units of the Thomson
! cross-section, with the resonances of the Lyman lines.  Below the Lyman
! limit the Kramers-Heisenberg amplitude
!   A(x) = sum_n f_n x^2/(x_n^2 - x^2) + Int_1^inf (df/dx') x^2/(x'^2 - x^2) dx',
! x = nu/nu_1 (nu_1 the Lyman limit), sigma = sigma_T A^2 (Nussbaumer,
! Schmid & Vogel 1989, A&A 211, L27; Lee & Kim 2004, MNRAS 347, 802):
! hydrogenic oscillator strengths f_n of 1s - np (n = 2..ray_nmax, x_n =
! 1 - 1/n^2) and the oscillator-strength density of the exact hydrogenic
! 1s photoionization cross-section (the sums give 0.565 and 0.435).  Far
! from the resonances this is Rayleigh's sigma = 1.265 sigma_T x^4 (as
! Dalgarno's fit, Kurucz's HRAYOP); near Lyman alpha it is the Lorentz wing
! of the line, which grows as 1/(nu - nu_n)^2.  There the cross-section is
! capped at 1e4 f_n^2 sigma_T of the nearest resonance (n <= 4), within
! about 1500 km/s of Lyman alpha 1732 sigma_T, as CMFGEN's RAYLEIGH_SCAT
! (the lines themselves carry their cores).  Above the Lyman limit
! sigma = sigma_T (1 + 1.66/sqrt(x)) (Gavrila 1967), where photoionization
! dominates.  The scattering is coherent in the comoving frame, as the
! other line scattering of the group opacity (grd_capsc).
!-----------------------------------------------------------------------
  implicit none
  private
  real*8,parameter :: pi = 3.141592653589793d0
!-- Lyman limit of H I [Hz] (109678.764 cm^-1); pi e^2/(m_e c) [cm^2 Hz]
  real*8,parameter :: ray_nu1 = 2.99792458d10*109678.764d0
  real*8,parameter :: fsob = 0.02654008d0
!-- (2^9 pi^2/3) alpha a_0^2 [cm^2]: the 1s photoionization cross-section
!-- is a0s x^-4 exp(-4 atan(e)/e)/(1 - exp(-2 pi/e)), e = sqrt(x-1)
  real*8,parameter :: a0s = 512d0*pi**2/3d0/137.035999d0*(0.529177211d-8)**2
  integer,parameter :: ray_nmax = 200
!-- the continuum term of the amplitude tabulated in x on [0, xtab]
  integer,parameter :: ntab = 2000
  real*8,parameter :: xtab = 0.9999d0
  real*8 :: fosc(2:ray_nmax),xres(2:ray_nmax)
  real*8 :: actab(0:ntab)
  logical :: ray_ready = .false.
!-- group averages (ray_groups), cached for the group boundaries in use
  real*8,allocatable :: ray_gsig(:)
  real*8,allocatable :: ray_gnu(:)
!
  public :: ray_sigma,ray_average,ray_groups,ray_init
!
contains
!
!
subroutine ray_init
!-- oscillator strengths, resonance positions and the continuum term
  integer :: n,i,k
  integer,parameter :: nq = 4000
  real*8 :: x,t,xp,s,w,dt
  do n=2,ray_nmax
!-- f(1s-np) = 2^8 n^5 (n-1)^(2n-4)/(3 (n+1)^(2n+4))
     fosc(n) = exp(8d0*log(2d0) + 5d0*log(dble(n)) + (2*n-4)*log(dble(n-1)) - &
          log(3d0) - (2*n+4)*log(dble(n+1)))
     xres(n) = 1d0 - 1d0/dble(n)**2
  enddo
!-- Int_1^inf (df/dx') x^2/(x'^2 - x^2) dx' with x' = 1/t, by Simpson's rule
!-- in t on (0, 1] (the integrand vanishes as t^3.5 at t = 0)
  dt = 1d0/nq
  do i=0,ntab
     x = xtab*dble(i)/ntab
     s = 0d0
     do k=1,nq
        t = k*dt
        xp = 1d0/t
        w = merge(1d0,merge(4d0,2d0,mod(k,2)==1),k==nq)
        s = s + w*dfdx(xp)*x*x/(xp*xp - x*x)/(t*t)
     enddo
     actab(i) = s*dt/3d0
  enddo
  ray_ready = .true.
end subroutine ray_init
!
!
pure function dfdx(x) result(df)
!-- oscillator-strength density of the 1s continuum per unit x = nu/nu_1
  real*8,intent(in) :: x
  real*8 :: df,e
  e = sqrt(max(x-1d0,1d-30))
  df = a0s/x**4*exp(-4d0*atan(e)/e)/(1d0 - exp(-2d0*pi/e))*ray_nu1/fsob
end function dfdx
!
!
function ray_sigma(nu) result(sig)
!-- Rayleigh cross-section of H I (1s) at frequency nu [Hz] / sigma_T
  real*8,intent(in) :: nu
  real*8 :: sig,x,a,d,dmin,w
  integer :: n,ires,i
  if(.not.ray_ready) call ray_init
  x = nu/ray_nu1
  if(x >= 1d0) then
     sig = 1d0 + 1.66d0/sqrt(x)
     return
  endif
  a = 0d0
  dmin = huge(1d0)
  ires = 2
  do n=2,ray_nmax
     d = (xres(n)/x)**2 - 1d0
     if(abs(d) < dmin) then
        dmin = abs(d)
        ires = n
     endif
!-- (keep the term finite at the resonance; the cap below applies there)
     d = sign(max(abs(d),1d-9),d)
     a = a + fosc(n)/d
  enddo
  w = min(x,xtab)/xtab*ntab
  i = min(int(w),ntab-1)
  w = w - i
  a = a + (1d0-w)*actab(i) + w*actab(i+1)
  sig = min(a*a,1d4*fosc(min(ires,4))**2)
end function ray_sigma
!
!
function ray_average(nu1,nu2) result(sig)
!-- mean of ray_sigma over [nu1, nu2] (midpoint rule in nu)
  real*8,intent(in) :: nu1,nu2
  real*8 :: sig
  integer,parameter :: ns = 64
  integer :: k
  sig = 0d0
  do k=1,ns
     sig = sig + ray_sigma(nu1 + (k-0.5d0)*(nu2-nu1)/ns)
  enddo
  sig = sig/ns
end function ray_average
!
!
subroutine ray_groups(nb,nub,sig)
!-- ray_average over the nb-1 intervals of the boundaries nub [Hz]; the
!-- result is cached while the boundaries stay the same
  integer,intent(in) :: nb
  real*8,intent(in) :: nub(nb)
  real*8,intent(out) :: sig(nb-1)
  integer :: k
  if(allocated(ray_gnu)) then
     if(size(ray_gnu)==nb) then
        if(all(ray_gnu==nub)) then
           sig = ray_gsig
           return
        endif
     endif
     deallocate(ray_gnu,ray_gsig)
  endif
  allocate(ray_gnu(nb),ray_gsig(nb-1))
  ray_gnu = nub
  do k=1,nb-1
     ray_gsig(k) = ray_average(min(nub(k),nub(k+1)),max(nub(k),nub(k+1)))
  enddo
  sig = ray_gsig
end subroutine ray_groups
!
end module rayleighmod
