!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module bbsobmod
!-----------------------------------------------------------------------
! Bound-bound lines of the LTE elements in the Sobolev approximation for
! arbitrary 1D velocity fields (expansion opacity), and the probability
! that a line interaction thermalizes the photon.
!
! A line with frequency-integrated opacity k (cm^-1 s^-1) has the Sobolev
! optical depth tau(mu) = K/|Q(mu)|, K = k*lambda (s^-1), Q(mu) =
! mu^2 dv/dr + (1-mu^2) v/r.  Its direction-averaged contribution to the
! opacity of a group of width dnu is (1/(lambda dnu)) <|Q| (1-exp(-K/|Q|))>
! (Eastman & Pinto 1993 for v = r/t, here with the local gradient): the
! smeared opacity k/dnu for thin lines, saturating at <|Q|>/(lambda dnu)
! for thick ones, where the smeared opacity overestimates it by tau.
!
! A line interaction ends in collisional de-excitation (the photon joins
! the thermal pool) with probability eps = X/(X+beta) (two-level atom with
! Sobolev escape), X = C_ul/A_ul = 20.6 lambda^3 n_e P(y)/sqrt(T) (van
! Regemorter 1962, lambda in cm) and beta = <(1-exp(-tau))/tau>; otherwise
! the photon is re-emitted in the line (coherent scattering in the comoving
! frame).  bbs_eps >= 0 replaces eps by a constant (1: pure absorption,
! the released code).
!-----------------------------------------------------------------------
  use physconstmod
  implicit none
  private
  integer,parameter,public :: bbs_nq = 16      !directions (midpoints of mu in [0,1])
  logical,public :: bbs_sobolev = .false.      !expansion opacity (else smeared)
  real*8,public :: bbs_eps = 1d0               !<0: van Regemorter, else constant
  logical,public :: bbs_on = .false.           !line scattering opacity in use
!-- line scattering part of the bound-bound opacity
  real*4,allocatable,public :: grd_capsc(:,:)  !(ng,ncell)
  real*4,allocatable,public :: gas_capsc(:,:)  !(ng,gas_ncell)
!-- the share of Rayleigh scattering by H I in it
  real*4,allocatable,public :: grd_fray(:,:)   !(ng,ncell)
  real*4,allocatable,public :: gas_fray(:,:)   !(ng,gas_ncell)
  public :: bbs_init,bbs_quad,bbs_fq,bbs_epsl,bbs_pvr
!-- van Regemorter (1962) effective Gaunt factor at y = dE/kT
  real*8,parameter :: ytab(10) = [.01d0,.02d0,.04d0,.1d0,.2d0,.4d0,1d0,2d0,4d0,10d0]
  real*8,parameter :: pneu(10) = [1.16d0,.956d0,.758d0,.493d0,.331d0,.209d0,.100d0,.063d0,.040d0,.023d0]
  real*8,parameter :: pion(10) = [1.16d0,.977d0,.788d0,.554d0,.403d0,.290d0,.214d0,.201d0,.200d0,.200d0]

contains

subroutine bbs_init(ng,ncell,gncell)
  integer,intent(in) :: ng,ncell,gncell
  bbs_on = bbs_sobolev .or. bbs_eps/=1d0
  allocate(grd_capsc(ng,ncell),gas_capsc(ng,gncell))
  grd_capsc = 0.
  gas_capsc = 0.
  allocate(grd_fray(ng,ncell),gas_fray(ng,gncell))
  grd_fray = 0.
  gas_fray = 0.
end subroutine bbs_init


pure subroutine bbs_quad(a,b,qabs,q,w,qavg,qmn,qmx)
!-- |Q(mu)| at bbs_nq directions mu in [0,1] for dv/dr = a, v/r = b, floored
!-- as in the transport and at qabs (e.g. a thermal speed over r), with
!-- quadrature weights (sum 1): two 8-point Gauss-Legendre panels, split at
!-- the direction where Q changes sign (compressing flows)
  real*8,intent(in) :: a,b,qabs
  real*8,intent(out) :: q(bbs_nq),w(bbs_nq),qavg,qmn,qmx
  real*8,parameter :: t8(4) = [.1834346424956498d0,.5255324099163290d0, &
       .7966664774136267d0,.9602898564975363d0]
  real*8,parameter :: w8(4) = [.3626837833783620d0,.3137066458778873d0, &
       .2223810344533745d0,.1012285362903763d0]
  real*8,parameter :: tt(8) = [-t8(4),-t8(3),-t8(2),-t8(1),t8(1),t8(2),t8(3),t8(4)]
  real*8,parameter :: ww(8) = [w8(4),w8(3),w8(2),w8(1),w8(1),w8(2),w8(3),w8(4)]
  integer :: k,j,i
  real*8 :: mu,qfl,ms,lo(2),hi(2),c,h
  qfl = max(1d-3*max(abs(a),abs(b)),qabs)
!-- panels
  ms = .5d0
  if(a*b<0d0) ms = sqrt(b/(b-a))
  lo = [0d0,ms]
  hi = [ms,1d0]
  k = 0
  do j=1,2
     c = .5d0*(lo(j)+hi(j))
     h = .5d0*(hi(j)-lo(j))
     do i=1,8
        k = k + 1
        mu = c + h*tt(i)
        q(k) = max(abs(mu*mu*a+(1d0-mu*mu)*b),qfl)
        w(k) = h*ww(i)
     enddo
  enddo
  qavg = sum(w*q)
  qmn = minval(q)
  qmx = maxval(q)
end subroutine bbs_quad


pure function bbs_fq(kk,q,w,qavg,qmn,qmx) result(f)
!-- <|Q| (1-exp(-K/|Q|))> over directions; beta = f/K
  real*8,intent(in) :: kk,q(bbs_nq),w(bbs_nq),qavg,qmn,qmx
  real*8 :: f,t
  integer :: k
  if(kk<=1d-6*qmn) then
     f = kk
  elseif(kk>=40d0*qmx) then
     f = qavg
  else
     f = 0d0
     do k=1,bbs_nq
        t = kk/q(k)
        if(t<1d-6) then
           f = f + w(k)*kk*(1d0-.5d0*t)
        else
           f = f + w(k)*q(k)*(1d0-exp(-t))
        endif
     enddo
  endif
end function bbs_fq


pure function bbs_pvr(y,lion) result(p)
!-- van Regemorter P(y): tabulated, log-interpolated; Bethe limit below
!-- y = 0.01, threshold limits (0.066/sqrt(y) neutrals, 0.2 ions) above 10
  real*8,intent(in) :: y
  logical,intent(in) :: lion
  real*8 :: p,s
  integer :: k
  if(y<=ytab(1)) then
     p = pneu(1) + .276d0*log(ytab(1)/max(y,1d-30))
  elseif(y>=ytab(10)) then
     if(lion) then
        p = pion(10)
     else
        p = pneu(10)*sqrt(ytab(10)/y)
     endif
  else
     do k=1,9
        if(y<ytab(k+1)) exit
     enddo
     s = log(y/ytab(k))/log(ytab(k+1)/ytab(k))
     if(lion) then
        p = (1d0-s)*pion(k) + s*pion(k+1)
     else
        p = (1d0-s)*pneu(k) + s*pneu(k+1)
     endif
  endif
end function bbs_pvr


pure function bbs_epsl(wl,ne,temp,ii,beta) result(eps)
!-- thermalization probability of a line interaction (wl in cm, ne in
!-- cm^-3, ii: ionization stage, 1 = neutral)
  real*8,intent(in) :: wl,ne,temp,beta
  integer,intent(in) :: ii
  real*8 :: eps,y,x
  y = pc_h*pc_c/(wl*pc_kb*temp)
  x = 20.6d0*wl**3*ne*bbs_pvr(y,ii>1)/sqrt(temp)
  eps = x/(x+max(beta,1d-300))
end function bbs_epsl

end module bbsobmod
! vim: fdm=marker
