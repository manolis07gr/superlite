!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nltedismod
!-----------------------------------------------------------------------
! Level dissolution by the plasma microfield (Hummer & Mihalas 1988), as
! CMFGEN with DO_DIS (OCCUPATION_PROB, COMP_LEV_DIS_BLK, GENOPAETA,
! QUAD_MULTI), with the fit of Hubeny, Hummer & Lanz (1994, A&A 282,
! 151) to the occupation probability by charged perturbers.  Per cell
! (N_ion: all ions, n_e: electrons [cm^-3]):
!   b = (8.589e14 N_ion^(1/3)/n_e)^1.5,  a = 0.0009 n_e^(1/6)/sqrt(T),
!   x = (1 + a)^3.15,
! and per level of effective quantum number n = sqrt(R Z^2/nu0) (R the
! Rydberg frequency, Z the charge seen by the outer electron, nu0 the
! edge):
!   K = 16 n/(3 (n+1)^2) (1 for n <= 3),  beta = (K Z^3/n^4)^1.5 b,
!   f = 1.091 (x + 4 (Z-1) a) beta^2/(7.782 + x beta),  w = f/(1 + f);
! w = 1 for n <= 2Z.  The LTE population of a level carries w (the n* of
! the rate equations), and the bound-free continuum of every level goes
! on below its edge into the dissolved upper levels: at nu < nu0 the
! cross-section is its threshold value times 1 - w(n*), with
! n* = sqrt(R Z^2/(nu0 - nu)) the effective quantum number of the upper
! level, down to 0.8 nu0 and where 1 - w > 1e-4 (CMFGEN's QUAD_MULTI and
! PHOT_DIS_PARAMETER).  The lines are not changed (as CMFGEN).  At the
! densities of a IIn wind (n_e = 1e9 - 1e11) only levels above n = 35 -
! 60 dissolve; in a dense shell (1e13) those above n = 20, and the Balmer
! continuum reaches 3700 Å.
!-----------------------------------------------------------------------
  use nlteatommod
  implicit none
  private
!
  real*8,parameter :: ryd = 3.289395d15     !Rydberg frequency [Hz] (CMFGEN)
  real*8,parameter :: dmin = 1d-4           !smallest dissolved fraction (PHOT_DIS_PARAMETER)
  real*8,parameter :: fmin = 0.8d0          !lowest frequency of a dissolved continuum [nu0]
  real*8,parameter :: hh = 6.6260693d-27, kb = 1.3806505d-16, cc = 2.99792458d10
  integer,parameter :: nzd = 10             !charges with a cached critical n per cell
!
  logical,public :: nxd_on = .false.
!-- per level: charge, threshold cross-section [cm^2], groups (increasing
!-- wavelength) that the dissolved continuum can reach (0: none)
  real*8,allocatable,public :: nxd_z(:),nxd_sth(:)   !(nlv)
  integer,allocatable,public :: nxd_ga(:),nxd_gb(:)  !(nlv)
  real*8,allocatable :: nxd_nub(:)                   !(ng+1) group boundaries [Hz], decreasing
  integer :: nxd_ng = 0
!
  public :: nxd_setup,nxd_cell,nxd_cellx,nxd_wocc,nxd_lvw,nxd_dfrac,nxd_numin,nxd_gint,nxd_sample,nxd_weights,nxd_nion
!
contains
!
!
subroutine nxd_setup(ng,wlg)
!-----------------------------------------------------------------------
! Static data of the dissolved continua on the group grid wlg (group
! boundaries [cm], increasing): per level the threshold cross-section and
! the groups between 0.8 nu0 and nu0 (with nu > 1e14 Hz, CMFGEN's
! MIN_DIS).
!-----------------------------------------------------------------------
  integer,intent(in) :: ng
  real*8,intent(in) :: wlg(ng+1)
  integer :: lv,i,k,ig
  real*8 :: nu0,nulo
!
  if(allocated(nxd_z)) deallocate(nxd_z,nxd_sth,nxd_ga,nxd_gb,nxd_nub)
  allocate(nxd_z(nxa_nlv),nxd_sth(nxa_nlv),nxd_ga(nxa_nlv),nxd_gb(nxa_nlv),nxd_nub(ng+1))
  nxd_ng = ng
  nxd_nub = cc/wlg
  do lv=1,nxa_nlv
     i = nxa_lvion(lv)
     k = lv - nxa_ion(i)%lvoff
     nu0 = nxa_lvnu0(lv)
     nxd_z(lv) = nxa_ion(i)%zion
     nxd_sth(lv) = nxa_sigma(nxa_ion(i),k,nu0*(1d0+1d-12))
     nxd_ga(lv) = 0
     nxd_gb(lv) = 0
     nulo = max(fmin*nu0,1d14)
     if(nxd_sth(lv) <= 0d0 .or. nulo >= nu0) cycle
!-- groups with a part in (nulo, nu0): from the group containing nu0 to
!-- the one containing nulo
     do ig=1,ng
        if(nxd_nub(ig+1) < nu0 .and. nxd_nub(ig) > nulo) then
           if(nxd_ga(lv)==0) nxd_ga(lv) = ig
           nxd_gb(lv) = ig
        endif
     enddo
  enddo
end subroutine nxd_setup
!
!
pure subroutine nxd_cell(ne,nion,temp,p)
!-- parameters of a cell (b, a, x; p(1) = 0: no dissolution) for its
!-- electron and ion densities [cm^-3] and temperature [K]
  real*8,intent(in) :: ne,nion,temp
  real*8,intent(out) :: p(3)
  p = 0d0
  if(.not.nxd_on .or. ne <= 0d0 .or. nion <= 0d0 .or. temp <= 0d0) return
  p(1) = (8.589d14*nion**(1d0/3d0)/ne)**1.5d0
  p(2) = 0.0009d0*ne**0.1667d0/sqrt(temp)
  p(3) = (1d0 + p(2))**3.15d0
end subroutine nxd_cell
!
!
pure subroutine nxd_cellx(ne,nel,x,temp,p)
!-- parameters of a cell (nxd_cell) from its electron density, the NLTE
!-- element densities nel and state fractions x: the ions of the NLTE
!-- elements, and the electrons of the others as singly charged ions
  real*8,intent(in) :: ne,nel(:),x(:),temp
  real*8,intent(out) :: p(3)
  real*8 :: fe
  integer :: iel,s,o
  fe = 0d0
  do iel=1,nxa_nel
     o = nxa_el(iel)%xoff
     do s=1,nxa_el(iel)%ns
        fe = fe + nel(iel)*x(o+s)*nxa_el(iel)%chg(s)
     enddo
  enddo
  call nxd_cell(ne,nxd_nion(nel,x,max(ne-fe,0d0)),temp,p)
end subroutine nxd_cellx
!
!
pure function nxd_nion(nel,x,neoth) result(ni)
!-- number density of all ions [cm^-3]: the states of the NLTE elements
!-- with charge (element densities nel, state fractions x), and the
!-- electrons of the other elements neoth as singly charged ions
  real*8,intent(in) :: nel(:),x(:),neoth
  real*8 :: ni
  integer :: iel,s,o
  ni = max(neoth,0d0)
  do iel=1,nxa_nel
     o = nxa_el(iel)%xoff
     do s=1,nxa_el(iel)%ns
        if(nxa_el(iel)%chg(s) > 0d0) ni = ni + nel(iel)*x(o+s)
     enddo
  enddo
end function nxd_nion
!
!
pure subroutine nxd_weights(temp,p,w)
!-- shares of all levels in their states (nxa_weights) in a cell with the
!-- parameters p: in proportion to the dissolved LTE populations w g e^{-E/kT}
  real*8,intent(in) :: temp,p(3)
  real*8,intent(out) :: w(nxa_nlv)
  real*8 :: wo(nxa_nlv),ws(nxa_nst)
  integer :: lv,g
  call nxa_weights(temp,w)
  if(p(1) <= 0d0) return
  ws = 0d0
  do lv=1,nxa_nlv
     wo(lv) = nxd_lvw(lv,p)
     g = nxa_el(nxa_lvel(lv))%xoff + nxa_lvs(lv)
     ws(g) = ws(g) + wo(lv)*w(lv)
  enddo
  do lv=1,nxa_nlv
     g = nxa_el(nxa_lvel(lv))%xoff + nxa_lvs(lv)
     if(ws(g) > 0d0) w(lv) = w(lv)*wo(lv)/ws(g)
  enddo
end subroutine nxd_weights
!
!
pure function nxd_wocc(neff,z,p) result(w)
!-- occupation probability of a level of effective quantum number neff
!-- (charge z) in a cell with parameters p
  real*8,intent(in) :: neff,z,p(3)
  real*8 :: w
  w = 1d0
  if(p(1) <= 0d0 .or. neff <= 2.01d0*z) return
  w = 1d0 - dis(neff,z,p)
  if(w > 0.99999999d0) w = 1d0
end function nxd_wocc
!
!
pure function nxd_lvw(lv,p) result(w)
!-- occupation probability of global level lv
  integer,intent(in) :: lv
  real*8,intent(in) :: p(3)
  real*8 :: w,t
  w = 1d0
  if(p(1) <= 0d0) return
  t = ryd*nxd_z(lv)**2/nxa_lvnu0(lv)
  if(t <= 0d0 .or. t >= 961d0) return
  w = nxd_wocc(sqrt(t),nxd_z(lv),p)
end function nxd_lvw
!
!
pure function dis(neff,z,p) result(d)
!-- 1 - w = 1/(1 + f) for neff > 2z
  real*8,intent(in) :: neff,z,p(3)
  real*8 :: d,k,beta,y,f
  k = min(1d0,16d0*neff/(1d0+neff)**2/3d0)
  beta = (k*z**3/neff**4)**1.5d0*p(1)
  y = 1.091d0*(p(3) + 4d0*(z-1d0)*p(2))
  f = y*beta*beta/(7.782d0 + p(3)*beta)
  d = 1d0/(1d0 + f)
end function dis
!
!
pure function nxd_dfrac(nu,nu0,z,p) result(d)
!-- dissolved fraction 1 - w(n*) of the upper level that a photon of
!-- frequency nu below the edge nu0 reaches (0 above the edge, below
!-- 0.8 nu0 or 1e14 Hz, for n* <= 2z, and below 1e-4)
  real*8,intent(in) :: nu,nu0,z,p(3)
  real*8 :: d,neff
  d = 0d0
  if(p(1) <= 0d0 .or. nu >= nu0 .or. nu <= max(fmin*nu0,1d14)) return
  neff = sqrt(ryd*z*z/(nu0-nu))
  if(neff <= 2d0*z) return
  d = dis(neff,z,p)
  if(d <= dmin) d = 0d0
end function nxd_dfrac
!
!
pure function nxd_numin(nu0,z,p) result(numin)
!-- lowest frequency of the dissolved continuum below the edge nu0: where
!-- 1 - w(n*) falls to 1e-4, at least 0.8 nu0 and 1e14 Hz.  f = 1/1e-4 - 1
!-- at the critical beta of the quadratic y beta^2 - f x beta - 7.782 f = 0,
!-- K z^3/n^4 = (beta/b)^(2/3) =: q there: n = (z^3/q)^(1/4) where K = 1
!-- (n <= 3), else n^3 (1+n)^2 = 16 z^3/(3 q) by Newton's method
  real*8,intent(in) :: nu0,z,p(3)
  real*8 :: numin,y,ff,bs,q,n,c,g,dg
  integer :: it
  numin = nu0
  if(p(1) <= 0d0) return
  y = 1.091d0*(p(3) + 4d0*(z-1d0)*p(2))
  ff = 1d0/dmin - 1d0
  bs = (ff*p(3) + sqrt((ff*p(3))**2 + 4d0*y*ff*7.782d0))/(2d0*y)
  q = (bs/p(1))**(2d0/3d0)
  n = (z**3/q)**0.25d0
  if(n > 3d0) then
     c = 16d0*z**3/(3d0*q)
     n = c**0.2d0
     do it=1,8
        g = n**3*(1d0+n)**2 - c
        dg = 3d0*n**2*(1d0+n)**2 + 2d0*n**3*(1d0+n)
        n = n - g/dg
     enddo
  endif
  n = max(n,2d0*z)
  numin = max(nu0 - ryd*z*z/(n*n),fmin*nu0,1d14)
end function nxd_numin
!
!
pure subroutine nxd_gint(lv,ig,p,temp,nuc,dm,dh,dhs,dsp)
!-----------------------------------------------------------------------
! Integrals of the dissolved cross-section sigma_d = sigma_th (1 - w(n*))
! of level lv over the part of group ig between nuc (nxd_numin) and the
! edge, with 4-point Gauss-Legendre on it:
!   dm  = (1/dnu_g) Int sigma_d dnu                (group mean: opacity)
!   dh  = Int sigma_d/(h nu) dnu                   (photoionization per J)
!   dhs = Int sigma_d/(h nu) e^{h(nu0-nu)/kT} dnu  (stimulated recombination
!                                                   per J and LTE factor)
!   dsp = Int sigma_d nu^2 e^{h(nu0-nu)/kT} dnu    (spontaneous: alpha_d =
!                                                   4 pi phi T^-1.5 sahac (2/c^2) dsp)
!-----------------------------------------------------------------------
  integer,intent(in) :: lv,ig
  real*8,intent(in) :: p(3),temp,nuc
  real*8,intent(out) :: dm,dh,dhs,dsp
  real*8,parameter :: xg(4) = [-0.8611363115940526d0,-0.3399810435848563d0, &
       0.3399810435848563d0,0.8611363115940526d0]
  real*8,parameter :: wg(4) = [0.3478548451374538d0,0.6521451548625461d0, &
       0.6521451548625461d0,0.3478548451374538d0]
  real*8 :: nu0,a,b,nu,d,hkt,ex
  integer :: j
  dm = 0d0; dh = 0d0; dhs = 0d0; dsp = 0d0
  nu0 = nxa_lvnu0(lv)
  a = max(nxd_nub(ig+1),nuc)
  b = min(nxd_nub(ig),nu0)
  if(b <= a) return
  hkt = hh/(kb*temp)
  do j=1,4
     nu = 0.5d0*(a+b) + 0.5d0*(b-a)*xg(j)
     d = nxd_dfrac(nu,nu0,nxd_z(lv),p)*nxd_sth(lv)*0.5d0*(b-a)*wg(j)
     if(d <= 0d0) cycle
     ex = exp(min(hkt*(nu0-nu),700d0))
     dm = dm + d
     dh = dh + d/(hh*nu)
     dhs = dhs + d/(hh*nu)*ex
     dsp = dsp + d*nu*nu*ex
  enddo
  dm = dm/(nxd_nub(ig) - nxd_nub(ig+1))
end subroutine nxd_gint
!
!
pure subroutine nxd_sample(lv,p,temp,nuc,r1,r2,nu,ok)
!-- a frequency of the dissolved recombination emission of level lv
!-- (energy spectrum sigma_d nu^3 e^{h(nu0-nu)/kT} on [nuc, nu0]) by
!-- rejection from the uniform distribution (ok = f: rejected, draw again)
  integer,intent(in) :: lv
  real*8,intent(in) :: p(3),temp,nuc,r1,r2
  real*8,intent(out) :: nu
  logical,intent(out) :: ok
  real*8 :: nu0,f,fmax,hkt
  nu0 = nxa_lvnu0(lv)
  nu = nuc + r1*(nu0 - nuc)
  hkt = hh/(kb*temp)
  fmax = nu0**3*exp(min(hkt*(nu0-nuc),700d0))
  f = nxd_dfrac(nu,nu0,nxd_z(lv),p)*nu**3*exp(min(hkt*(nu0-nu),700d0))
  ok = r2*fmax <= f
end subroutine nxd_sample
!
end module nltedismod
