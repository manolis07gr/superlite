!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module dopplermod
!-----------------------------------------------------------------------
! Comoving-frame group crossings of a photon in free flight through a
! cell with a velocity field that is linear in radius, v(r) = a + b*r.
!
! Along a straight path from (r0,mu0):
!   r(s)^2 = (s + r0*mu0)^2 + p^2,   p^2 = r0^2*(1 - mu0^2)
!   mu(s)  = (s + r0*mu0)/r(s)
!   g(s)   = mu(s)*v(r(s)) = b*(s + r0*mu0) + a*(s + r0*mu0)/r(s)
! The lab-frame wavelength wl is constant in free flight; to first order in
! v/c the comoving wavelength is wl/(1 - g(s)/c).  A photon therefore leaves
! its comoving group [wlo,whi) where g(s) reaches
!   g_red  = c*(1 - wl/whi)   (redshift)  or   g_blue = c*(1 - wl/wlo)   (blueshift).
! g'(s) = b + a*p^2/r(s)^3 vanishes only where r(s) = (-a*p^2/b)^(1/3), so g
! has at most two stationary points along the path and is monotonic between
! them.  This covers accelerating, decelerating and non-monotonic flows
! (a and b of any sign), where photons can blueshift, and redshift and
! blueshift again along one path.
!-----------------------------------------------------------------------
  implicit none
  private
  public :: doppler_crossing, dop_g

contains

pure function dop_g(s,r0,mu0,a,b) result(g)
  real*8,intent(in) :: s,r0,mu0,a,b
  real*8 :: g
  real*8 :: q,r
  q = s + r0*mu0
  r = sqrt(q*q + r0*r0*max(0d0,1d0-mu0*mu0))
  if(r>0d0) then
     g = b*q + a*q/r
  else
     g = 0d0
  endif
end function dop_g


pure subroutine doppler_crossing(r0,mu0,xl,xr,vl,vr,wl,wlo,whi,smax,cspeed,s,idir)
  real*8,intent(in) :: r0,mu0      !position, lab-frame direction cosine
  real*8,intent(in) :: xl,xr,vl,vr !cell edges and edge velocities
  real*8,intent(in) :: wl          !lab-frame wavelength
  real*8,intent(in) :: wlo,whi     !comoving group boundaries, wlo < whi
  real*8,intent(in) :: smax        !path length to the cell boundary
  real*8,intent(in) :: cspeed
  real*8,intent(out) :: s          !distance to the crossing (huge if none)
  integer,intent(out) :: idir      !+1 red (through whi), -1 blue (through wlo), 0 none
!
  real*8 :: a,b,p2,gred,gblue,re,h,sa,sb,ga,gb
  real*8 :: sbrk(4)
  integer :: nbrk,i,j
!
  s = huge(s)
  idir = 0
  if(smax<=0d0) return
!
  b = (vr - vl)/(xr - xl)
  a = vl - b*xl
  p2 = r0*r0*max(0d0,1d0-mu0*mu0)
  gred = cspeed*(1d0 - wl/whi)
  gblue = cspeed*(1d0 - wl/wlo)
!
!-- breakpoints: path start, stationary points of g, path end
  nbrk = 1
  sbrk(1) = 0d0
  if(p2>0d0 .and. a*b<0d0) then
     re = (-a*p2/b)**(1d0/3d0)
     h = re*re - p2
     if(h>0d0) then
        h = sqrt(h)
        do j=-1,1,2
           sa = -r0*mu0 + j*h
           if(sa>0d0 .and. sa<smax) then
              nbrk = nbrk + 1
              sbrk(nbrk) = sa
           endif
        enddo
     endif
  endif
  nbrk = nbrk + 1
  sbrk(nbrk) = smax
!
!-- monotonic segments in path order
  do i=1,nbrk-1
     sa = sbrk(i)
     sb = sbrk(i+1)
     if(sb<=sa) cycle
     ga = dop_g(sa,r0,mu0,a,b)
     gb = dop_g(sb,r0,mu0,a,b)
     if(gb>ga) then
!-- comoving wavelength increasing: redshift through whi
        if(ga>=gred .and. i==1) then
           s = 0d0
           idir = 1
           return
        endif
        if(ga<gred .and. gb>=gred) then
           s = root(sa,sb,ga,gb,gred)
           idir = 1
           return
        endif
     elseif(gb<ga) then
!-- comoving wavelength decreasing: blueshift through wlo
        if(ga<=gblue .and. i==1) then
           s = 0d0
           idir = -1
           return
        endif
        if(ga>gblue .and. gb<=gblue) then
           s = root(sa,sb,ga,gb,gblue)
           idir = -1
           return
        endif
     endif
  enddo

contains

  pure function root(s1,s2,g1,g2,gt) result(sr)
!-- g monotonic on [s1,s2], g1 and g2 on opposite sides of gt
    real*8,intent(in) :: s1,s2,g1,g2,gt
    real*8 :: sr
    real*8 :: lo,hi,flo,fhi,fm,dg,q,r,sn
    integer :: it
    lo = s1
    hi = s2
    flo = g1 - gt
    fhi = g2 - gt
    sr = lo - flo*(hi-lo)/(fhi-flo)
    do it=1,100
!-- Newton step, kept inside the bracket
       q = sr + r0*mu0
       r = sqrt(q*q + p2)
       fm = b*q + a*q/r - gt
       if(fm==0d0) return
       if((fm<0d0) .eqv. (flo<0d0)) then
          lo = sr
          flo = fm
       else
          hi = sr
          fhi = fm
       endif
       dg = b + a*p2/(r*r*r)
       if(dg/=0d0) then
          sn = sr - fm/dg
       else
          sn = 0.5d0*(lo+hi)
       endif
       if(sn<=min(lo,hi) .or. sn>=max(lo,hi)) sn = 0.5d0*(lo+hi)
       if(abs(sn-sr)<=1d-13*max(abs(sr),abs(s2-s1))) then
          sr = sn
          return
       endif
       sr = sn
    enddo
  end function root

end subroutine doppler_crossing

end module dopplermod
! vim: fdm=marker
