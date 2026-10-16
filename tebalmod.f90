!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module tebalmod
!-----------------------------------------------------------------------
! Electron temperature from the thermal balance of the gas, with the
! radiation field of the previous Monte Carlo iteration:
!
!   f(T) = Lambda(T) - Gamma(T) = 0,
!
! per unit volume [erg/s/cm^3]:
!  * opacities of the LTE elements (free-free, bound-free, thermalized
!    share of the lines), kappa_g from the last opacity calculation:
!      4 pi sum_g kappa_g (B_g(T) - J_g) dnu_g
!  * NLTE hydrogen, populations re-solved at each trial T:
!      photoionization heating
!        4 pi n_H sum_n sum_g (x_n - x_n* e^{-h nu/kT}) sigma_n J_g (1 - nu_n/nu) dnu_g
!      recombination cooling  n_e n_p sum_n rc_n(T)
!      collisions  n_e sum_{l<u} (n_l q_lu - n_u q_ul) (E_u - E_l)
!                + n_e sum_n (n_n c_n - n_e n_p c3b_n) chi_n
!  * Compton  sigma_T n_e/(m_e c^2) 4 pi sum_g J_g dnu_g (4kT - h nu_g (1 + n_g)),
!    n_g = c^2 J_g/(2 h nu_g^3) the photon occupation number (induced
!    scattering; the balance is exact for a Planck field at T)
!  * heating not carried by the radiation field: energy injected in the
!    cell outside shocks ('lumg' source)
!  * photons above te_ntemin make fast photoelectrons: only the heat share
!    fh (Shull & van Steenberg 1985) of their absorption heats the gas, the
!    rest ionizes and excites H (rates nonth, nonexc in the populations)
!
! f is bracketed in log T and the root found by bisection/secant steps
! (Illinois).  The new temperature is mixed geometrically with the old one
! (te_damp) and its change per iteration limited to a factor te_maxfac.
!-----------------------------------------------------------------------
  use nltehmod
  implicit none
  private
  real*8,parameter :: pi = 3.141592653589793d0
  real*8,parameter :: hh = 6.6260693d-27, cc = 2.99792458d10, kb = 1.3806505d-16
  real*8,parameter :: sigt = 6.6524587d-25, mec2 = 9.1093826d-28*2.99792458d10**2
  real*8,parameter :: sahac = 4.1412d-16       !(h^2/(2 pi m_e k))^1.5 for g_n/g_p = n^2
  real*8,public :: te_damp = 0.5d0              !weight of the new temperature
  real*8,public :: te_maxfac = 1.5d0            !maximum change per iteration
  real*8,public :: te_tmin = 1d3, te_tmax = 1d7
  real*8,public :: te_ntemin = 0d0              !fast photoelectrons above [erg], 0: off
  real*8,allocatable :: sigtab(:,:)             !(nlev,ng) sigma_n at group centres
  integer :: tng = 0
  public :: te_init,te_balance,te_terms

contains

subroutine te_init(ng,nug)
!-- photoionization cross sections at the group centres
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng)
  integer :: n,ig
  tng = ng
  allocate(sigtab(nh_nlev,ng))
  do ig=1,ng
     do n=1,nh_nlev
        sigtab(n,ig) = nlteh_sigma(n,nug(ig))
     enddo
  enddo
end subroutine te_init


pure function rc_coef(n,temp) result(rc)
!-- spontaneous recombination cooling coefficient to level n [erg cm^3/s]:
!-- kinetic energy h(nu - nu_n) = y kT carried off per recombination,
!-- Milne relation as nlteh_alpha
  integer,intent(in) :: n
  real*8,intent(in) :: temp
  real*8 :: rc
  integer,parameter :: ni = 200
  real*8,parameter :: ymax = 40d0
  real*8 :: s,nu,y,w,h
  integer :: i
  h = ymax/ni
  s = 0d0
  do i=1,ni
     y = i*h
     w = merge(1d0, merge(4d0, 2d0, mod(i,2)==1), i==ni)
     nu = nh_nu0(n) + y*kb*temp/hh
     s = s + w*nlteh_sigma(n,nu)*2d0*nu*nu/cc**2*exp(-y)*y
  enddo
  s = s*h/3d0
  rc = 4d0*pi*dble(n*n)*sahac*temp**(-1.5d0)*(kb*temp/hh)*(kb*temp)*s
end function rc_coef


subroutine te_terms(temp,ng,nug,dnug,jg,capth,nh,x,ne,gam_ext,terms)
!-----------------------------------------------------------------------
! Heating and cooling terms at temperature temp for fixed populations.
! terms(1) LTE opacities, net cooling      terms(2) H photoionization heating
! terms(3) H recombination cooling         terms(4) H collisional net cooling
! terms(5) Compton net cooling             terms(6) external heating
! f = t1 - t2 + t3 + t4 + t5 - t6
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,nug(ng),dnug(ng),jg(ng),nh,x(0:nh_nlev),ne,gam_ext
  real*4,intent(in) :: capth(ng)
  integer,intent(in) :: ng
  real*8,intent(out) :: terms(6)
  integer :: ig,n,l,u
  real*8 :: hkt,bg,s,np,saha,xs,qup(nh_nlev,nh_nlev),qdn(nh_nlev,nh_nlev),cup,c3b,eu
  real*8 :: fh,fi,fe,jh
!
  hkt = hh/(kb*temp)
  np = x(0)*nh
  terms = 0d0
  fh = 1d0
  if(te_ntemin > 0d0) call nlteh_svs(x(0),fh,fi,fe)
!-- LTE opacities and Compton
  do ig=1,ng
     if(hkt*nug(ig) < 700d0) then
        bg = 2d0*hh*nug(ig)**3/cc**2/(exp(hkt*nug(ig)) - 1d0)
     else
        bg = 0d0
     endif
     jh = jg(ig)
     if(te_ntemin > 0d0 .and. hh*nug(ig) > te_ntemin) jh = fh*jg(ig)
     terms(1) = terms(1) + capth(ig)*(bg - jh)*dnug(ig)
     terms(5) = terms(5) + jg(ig)*dnug(ig)*(4d0*kb*temp - &
          hh*nug(ig)*(1d0 + cc**2*jg(ig)/(2d0*hh*nug(ig)**3)))
  enddo
  terms(1) = 4d0*pi*terms(1)
  terms(5) = 4d0*pi*sigt*ne/mec2*terms(5)
!-- H photoionization heating (net of stimulated recombination)
  if(nh > 0d0) then
     saha = ne*x(0)*sahac*temp**(-1.5d0)
     s = 0d0
     do ig=1,ng
        do n=1,nh_nlev
           if(nug(ig) < nh_nu0(n)) cycle
           xs = saha*n*n*exp(-hkt*(nug(ig) - nh_nu0(n)))
           jh = jg(ig)
           if(te_ntemin > 0d0 .and. hh*nug(ig) > te_ntemin) jh = fh*jg(ig)
           s = s + (x(n) - xs)*sigtab(n,ig)*jh*(1d0 - nh_nu0(n)/nug(ig))*dnug(ig)
        enddo
     enddo
     terms(2) = 4d0*pi*nh*s
!-- recombination cooling
     s = 0d0
     do n=1,nh_nlev
        s = s + rc_coef(n,temp)
     enddo
     terms(3) = ne*np*s
!-- collisions: excitation/de-excitation and ionization/three-body recombination
     call nlteh_qcoll(temp,qup,qdn)
     s = 0d0
     do l=1,nh_nlev-1
        do u=l+1,nh_nlev
           eu = nh_eion(l) - nh_eion(u)
           s = s + (x(l)*qup(l,u) - x(u)*qdn(u,l))*eu
        enddo
     enddo
     do n=1,nh_nlev
        call nlteh_cion(n,temp,cup,c3b)
        s = s + (x(n)*cup - x(0)*ne*c3b)*nh_eion(n)
     enddo
     terms(4) = ne*nh*s
  endif
  terms(6) = gam_ext
end subroutine te_terms


subroutine te_balance(told,trad,nh,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nonth,nonexc,capth, &
     gam_ext,x,tnew,ierr)
!-----------------------------------------------------------------------
! Solve f(T) = 0 with the H populations re-solved at each trial T; x is
! the initial guess for the populations (unchanged on return).  tnew is
! the damped and limited new temperature; ierr = 1 if f has no sign
! change in [te_tmin, te_tmax] (tnew moves toward the bound).
!-----------------------------------------------------------------------
  real*8,intent(in) :: told,trad,nh,neoth,qb,qv,nonth,nonexc,gam_ext
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng),jb(nh_nlin)
  logical,intent(in) :: ljb(nh_nlin)
  real*4,intent(in) :: capth(ng)
  real*8,intent(in) :: x(0:nh_nlev)
  real*8,intent(out) :: tnew
  integer,intent(out) :: ierr
  integer :: it,side
  real*8 :: a,b,c,fa,fb,fc,tsol,lo,hi
!
  ierr = 0
!-- bracket in log T around the old temperature
  a = log(told)
  fa = fval(a)
  b = a
  fb = fa
  lo = log(te_tmin)
  hi = log(te_tmax)
  do it=1,40
     if(fa < 0d0) then
        b = min(a + 0.3d0, hi)
     else
        b = max(a - 0.3d0, lo)
     endif
     fb = fval(b)
     if(fa*fb <= 0d0) exit
     if(b <= lo .or. b >= hi) exit
     a = b
     fa = fb
  enddo
  if(fa*fb > 0d0) then
!-- no root: the bound
     ierr = 1
     tsol = exp(b)
  else
!-- Illinois (regula falsi) iterations
     side = 0
     do it=1,60
        c = (a*fb - b*fa)/(fb - fa)
        fc = fval(c)
        if(abs(b - a) < 1d-4 .or. fc == 0d0) exit
        if(fc*fb < 0d0) then
           a = b
           fa = fb
           b = c
           fb = fc
           side = 0
        else
           b = c
           fb = fc
           if(side == -1) fa = 0.5d0*fa
           side = -1
        endif
     enddo
     tsol = exp(c)
  endif
!-- damping and limiter
  tnew = told**(1d0-te_damp)*tsol**te_damp
  tnew = min(max(tnew, told/te_maxfac), told*te_maxfac)
  tnew = min(max(tnew, te_tmin), te_tmax)
contains
  real*8 function fval(lt)
    real*8,intent(in) :: lt
    real*8 :: t,xx(0:nh_nlev),ne,beta(nh_nlin),terms(6)
    integer :: ie
    t = exp(lt)
    xx = x
    ne = neoth + xx(0)*nh
    if(nh > 0d0) then
       call nlteh_solve_cell(t,trad,nh,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nonth,xx,ne,beta,ie,nonexc)
       if(ie /= 0) then
          xx = x
          ne = neoth + xx(0)*nh
       endif
    endif
    call te_terms(t,ng,nug,dnug,jg,capth,nh,xx,ne,gam_ext,terms)
    fval = terms(1) - terms(2) + terms(3) + terms(4) + terms(5) - terms(6)
  end function fval
end subroutine te_balance

end module tebalmod
! vim: fdm=marker
