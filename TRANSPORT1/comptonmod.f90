!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module comptonmod
!-----------------------------------------------------------------------
! Thomson scattering off thermal electrons (non-coherent electron
! scattering) in the comoving frame of the gas.
!
! The electron velocity is drawn from a Maxwellian at the electron
! temperature, weighted by the relative flux factor (1 - n.beta).  The photon
! is boosted into the electron rest frame, scattered there (isotropic or
! dipole, elastic apart from Compton recoil) and boosted back.  To first
! order this broadens a monochromatic line by
!   dnu/nu = (n_out - n_in).beta_e,  rms = sqrt(2(1-cos theta) kT/(m_e c^2)),
! about 550 km/s per scattering at 10^4 K, and to second order it gives the
! Kompaneets mean shift (4kT - h nu)/(m_e c^2) per scattering.  Multiple
! scattering in an ionized CSM builds the broad, roughly exponential wings
! of SN IIn lines (Chugai 2001; Huang & Chevalier 2018).
!
! The Maxwellian is non-relativistic; for kT/(m_e c^2) > 0.05 (T > 3e8 K)
! the electron speed distribution should be Maxwell-Juttner.
!-----------------------------------------------------------------------
  use randommod
  implicit none
  private
  public :: thermal_thomson

  real*8,parameter :: kb = 1.3806505d-16, me = 9.1093826d-28, cl = 2.99792458d10
  real*8,parameter :: hpl = 6.6260693d-27
  real*8,parameter :: pi2 = 6.283185307179586d0

contains

pure subroutine gauss2(z1,z2,state)
!-- two independent standard normal deviates (Box-Muller)
  real*8,intent(out) :: z1,z2
  type(rnd_t),intent(inout) :: state
  real*8 :: u1,u2,r
  call rnd_r(u1,state)
  call rnd_r(u2,state)
  u1 = max(u1,1d-300)
  r = sqrt(-2d0*log(u1))
  z1 = r*cos(pi2*u2)
  z2 = r*sin(pi2*u2)
end subroutine gauss2


pure subroutine iso_dir(n,state)
  real*8,intent(out) :: n(3)
  type(rnd_t),intent(inout) :: state
  real*8 :: u1,u2,st
  call rnd_r(u1,state)
  call rnd_r(u2,state)
  n(3) = 1d0 - 2d0*u1
  st = sqrt(max(0d0,1d0-n(3)**2))
  n(1) = st*cos(pi2*u2)
  n(2) = st*sin(pi2*u2)
end subroutine iso_dir


pure subroutine dipole_dir(n0,n,state)
!-- direction with p(cos t) ~ 1 + cos^2 t about the unit vector n0
  real*8,intent(in) :: n0(3)
  real*8,intent(out) :: n(3)
  type(rnd_t),intent(inout) :: state
  real*8 :: u,ct,st,ph,e1(3),e2(3),a
  do
     call rnd_r(u,state)
     ct = 1d0 - 2d0*u
     call rnd_r(u,state)
     if(2d0*u <= 1d0 + ct*ct) exit
  enddo
  call rnd_r(u,state)
  ph = pi2*u
  st = sqrt(max(0d0,1d0-ct*ct))
!-- orthonormal basis around n0
  if(abs(n0(3)) < 0.9d0) then
     e1 = [-n0(2), n0(1), 0d0]
  else
     e1 = [0d0, -n0(3), n0(2)]
  endif
  a = sqrt(sum(e1*e1))
  e1 = e1/a
  e2 = [n0(2)*e1(3)-n0(3)*e1(2), n0(3)*e1(1)-n0(1)*e1(3), n0(1)*e1(2)-n0(2)*e1(1)]
  n = ct*n0 + st*(cos(ph)*e1 + sin(ph)*e2)
end subroutine dipole_dir


pure subroutine thermal_thomson(mu_in,wl_in,tempe,dipole,state,mu_out,ratio)
!-----------------------------------------------------------------------
! mu_in   comoving direction cosine (relative to the radial direction)
! wl_in   comoving wavelength [cm] (for the recoil term)
! tempe   electron temperature [K]
! dipole  .true.: dipole phase function in the electron frame, else isotropic
! mu_out  comoving direction cosine after scattering
! ratio   nu_out/nu_in in the comoving frame (= wl_in/wl_out)
!-----------------------------------------------------------------------
  real*8,intent(in) :: mu_in,wl_in,tempe
  logical,intent(in) :: dipole
  type(rnd_t),intent(inout) :: state
  real*8,intent(out) :: mu_out,ratio
!
  real*8 :: sig,b(3),bb,gam,nin(3),ne(3),nout(3),no(3),z(4),u
  real*8 :: dop_in,dop_out,nb,x_e,cth
!
  nin = [sqrt(max(0d0,1d0-mu_in*mu_in)), 0d0, mu_in]
  sig = sqrt(kb*max(tempe,0d0)/(me*cl*cl))
!
!-- electron velocity: Maxwellian weighted by the flux factor (1 - n.beta)
  do
     call gauss2(z(1),z(2),state)
     call gauss2(z(3),z(4),state)
     b = sig*z(1:3)
     bb = sum(b*b)
     if(bb >= 0.99d0) cycle
     nb = sum(nin*b)
     call rnd_r(u,state)
     if(u*(1d0+sqrt(bb)) <= 1d0 - nb) exit
  enddo
  gam = 1d0/sqrt(1d0-bb)
!
!-- into the electron rest frame
  dop_in = gam*(1d0 - nb)
  if(bb > 0d0) then
     ne = (nin + ((gam-1d0)*nb/bb - gam)*b)/dop_in
  else
     ne = nin
  endif
  ne = ne/sqrt(sum(ne*ne))
!
!-- scatter (elastic in the electron frame apart from recoil)
  if(dipole) then
     call dipole_dir(ne,no,state)
  else
     call iso_dir(no,state)
  endif
  cth = sum(ne*no)
  x_e = hpl*cl/wl_in*dop_in/(me*cl*cl)       !h nu_e/(m_e c^2)
!
!-- back to the comoving frame
  nb = sum(no*b)
  dop_out = gam*(1d0 + nb)
  if(bb > 0d0) then
     nout = (no + ((gam-1d0)*nb/bb + gam)*b)/dop_out
  else
     nout = no
  endif
  mu_out = max(-1d0,min(1d0,nout(3)/sqrt(sum(nout*nout))))
  ratio = dop_in*dop_out/(1d0 + x_e*(1d0-cth))
end subroutine thermal_thomson

end module comptonmod
! vim: fdm=marker
