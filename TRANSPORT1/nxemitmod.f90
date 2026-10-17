!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nxemitmod
!-----------------------------------------------------------------------
! Re-emission of a packet after a macro-atom or k-packet process of the
! NLTE elements, shared by the IMC (transport11) and DDMC (diffusion11)
! steps and the interior source: the comoving wavelength and group of the
! emission, the packet placed in the lab frame (on the line it was emitted
! in, or into DDMC), the escape direction of a line photon, and the
! electron scattering of line photons trapped in their resonance zone.
!-----------------------------------------------------------------------
  use randommod
  use miscmod
  use gridmod
  use groupmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  use nltextrmod
  use bbsobmod, only:grd_capsc
  use comptonmod, only:thermal_thomson
  use nlteatommod, only:nxa_lwl
  implicit none
  private
  real*8,parameter :: cinv = 1d0/pc_c
  public :: nx_emitwl,nx_place,nx_reemit,line_mu,nx_linees,nx_ilem
!
contains
!
!
  pure subroutine nx_emitwl(kind,idx,wlcmf,ic,rndstate,wlc,ig,lout)
!-----------------------------------------------------------------------
! Comoving wavelength wlc and group ig of an emission of kind em_line,
! em_lines, em_bf, em_2ph (at wlcmf) or em_thermal (the emissivity of the
! other opacities, uniform in frequency within the group).  lout: outside
! the group grid.
!-----------------------------------------------------------------------
    integer,intent(in) :: kind,idx,ic
    real*8,intent(in) :: wlcmf
    type(rnd_t),intent(inout) :: rndstate
    real*8,intent(out) :: wlc
    integer,intent(out) :: ig
    logical,intent(out) :: lout
    real*8 :: rr
    lout = .false.
    if(kind==em_thermal) then
       call rnd_r(rr,rndstate)
       ig = emitgroup(rr,ic)
       call rnd_r(rr,rndstate)
       wlc = 1d0/((1d0-rr)*grp_wlinv(ig)+rr*grp_wlinv(ig+1))
    else
       wlc = wlcmf
       if(kind==em_line) then
          ig = nxt_lgrp(idx)
       else
          ig = 0
       endif
       if(ig<1) ig = binsrch(wlc,grp_wl,grp_ng+1,.false.)
       lout = wlc<grp_wl(1) .or. wlc>=grp_wl(grp_ng+1)
    endif
  end subroutine nx_emitwl
!
!
  pure subroutine nx_place(kind,idx,wlc,ig,lout,wl,e,e0,mu,vx,x,ix,ic,ptcl2,totevelo)
!-----------------------------------------------------------------------
! Packet of comoving wavelength wlc in group ig after a re-emission at x
! (velocity vx); mu (lab) and the lab energy are already set for the new
! direction.  Outside the group grid (e.g. mid-IR fine-structure lines)
! the medium is not described: the photon escapes (flux tally: outside
! the bins unless the flux grid is wider).
!-----------------------------------------------------------------------
    integer,intent(in) :: kind,idx,ig,ix,ic
    logical,intent(in) :: lout
    real*8,intent(in) :: wlc,mu,vx,x
    real*8,intent(inout) :: wl,e,e0,totevelo
    type(packet2),intent(inout) :: ptcl2
    real*8 :: qn,dvdrl
    if(lout) then
       wl = wlc*(1d0-vx*mu*cinv)
       ptcl2%stat = 'flux'
       return
    endif
    if(kind==em_line .and. nxt_isrt(idx)>0) then
!-- the packet sits on the (active) line it was emitted in
       ptcl2%iline = nxt_isrt(idx)
       dvdrl = (grd_vxarr(ix+1)-grd_vxarr(ix))/(grd_xarr(ix+1)-grd_xarr(ix))
       qn = ((mu-vx*cinv)/(1d0-mu*vx*cinv))**2
       qn = qn*dvdrl + (1d0-qn)*vx/x
       ptcl2%ilside = merge(1,-1,qn>=0d0)
    endif
    wl = wlc*(1d0-vx*mu*cinv)
!-- DDMC in the new group
    if(((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*(grd_xarr(ix+1)-grd_xarr(ix)) >= trn_tauddmc) &
         .and. .not.in_puretran) then
       ptcl2%itype = 2
       ptcl2%iline = 0
       totevelo = totevelo+e*vx*mu*cinv
       e = e*(1d0-vx*mu*cinv)
       e0 = e0*(1d0-vx*mu*cinv)
       wl = wlc
    endif
  end subroutine nx_place
!
!
  pure subroutine nx_reemit(kind,idx,wlcmf,ig,wl,e,e0,mu,vx,x,ix,ic,ptcl2,totevelo,rndstate)
!-----------------------------------------------------------------------
! Re-emit the packet after a macro-atom/k-packet process with comoving
! wavelength wlcmf (kind em_line/em_bf/em_2ph) or the thermal emissivity
! (em_thermal).  mu (lab) and the lab energy are already set for the new
! direction.
!-----------------------------------------------------------------------
    integer,intent(in) :: kind,idx,ix,ic
    real*8,intent(in) :: wlcmf,mu,vx,x
    integer,intent(inout) :: ig
    real*8,intent(inout) :: wl,e,e0,totevelo
    type(packet2),intent(inout) :: ptcl2
    type(rnd_t),intent(inout) :: rndstate
    real*8 :: wlc
    logical :: lout
    call nx_emitwl(kind,idx,wlcmf,ic,rndstate,wlc,ig,lout)
    call nx_place(kind,idx,wlc,ig,lout,wl,e,e0,mu,vx,x,ix,ic,ptcl2,totevelo)
  end subroutine nx_reemit
!
!
  pure subroutine nx_linees(kind,idx,ic,rndstate,mucmf,wlcmf)
!-----------------------------------------------------------------------
! Where the trapped-photon electron scattering is in the rates and the
! macro-atom (nxt_useb; otherwise only in the formal integral): a photon of
! line idx leaving the macro-atom in cell ic was scattered by a thermal
! electron while trapped in the resonance zone with probability
! bes/(beta+bes) (nxt_bes).  It then leaves the line centre with the
! electron's Doppler shift (comoving wavelength wlcmf/ratio) in a new
! comoving direction mucmf, kind em_lines; a shift to the blue brings it
! back to the line further along its path.  Otherwise kind stays em_line
! (mucmf is not set).
!-----------------------------------------------------------------------
    integer,intent(inout) :: kind
    integer,intent(in) :: idx,ic
    type(rnd_t),intent(inout) :: rndstate
    real*8,intent(inout) :: mucmf,wlcmf
    real*8 :: r,b,ratio
    if(kind/=em_line) return
    if(.not.nxt_useb(idx,ic)) return
    b = nxt_bes(idx,ic)
    if(b <= 0d0) return
    call rnd_r(r,rndstate)
    if(r*(nxt_beta(idx,ic)+b) >= b) return
!-- the photon's direction in the zone is isotropic
    call rnd_r(r,rndstate)
    call thermal_thomson(1d0-2d0*r,wlcmf,1d0/grd_tempinv(ic),trn_thm_aniso, &
         rndstate,mucmf,ratio)
    wlcmf = wlcmf/ratio
    kind = em_lines
  end subroutine nx_linees
!
!
  pure function nx_ilem(kind,idx,wlcmf) result(ilem)
!-- packet%ilem after an emission: the line of a line photon while it is on
!-- the red side of the line (comoving wavelength wlcmf)
    integer,intent(in) :: kind,idx
    real*8,intent(in) :: wlcmf
    integer :: ilem
    ilem = 0
    if(kind==em_line) then
       ilem = idx
    elseif(kind==em_lines) then
       if(wlcmf > nxa_lwl(idx)) ilem = idx
    endif
  end function nx_ilem
!
!
  pure subroutine line_mu(tau0,a,b,rndstate,mucmf)
!-----------------------------------------------------------------------
! Comoving direction of a photon leaving a Sobolev resonance in a flow
! with dv/dr = a and v/r = b.  After any number of local scatterings the
! photon escapes in direction mu with probability density proportional to
! beta(mu) = (1-exp(-tau(mu)))/tau(mu), tau(mu) = tau0/|Q(mu)|,
! Q(mu) = mu**2*a + (1-mu**2)*b (rejection sampling; isotropic when
! a = b, as in homologous flows, or for optically thin lines).
!-----------------------------------------------------------------------
    real*8,intent(in) :: tau0,a,b
    type(rnd_t),intent(inout) :: rndstate
    real*8,intent(out) :: mucmf
    real*8 :: qmax,qlo,q,fmax,rr
    integer :: itry
    qmax = max(abs(a),abs(b))
    call rnd_r(rr,rndstate)
    mucmf = 1d0-2d0*rr
    if(tau0<=0d0 .or. qmax<=0d0) return
    qlo = 1d-3*qmax
!-- q*(1-exp(-tau0/q)) grows with q: maximal at qmax
    fmax = qmax*(1d0-exp(-tau0/qmax))
    do itry=1,10000
       q = max(abs(mucmf**2*a+(1d0-mucmf**2)*b),qlo)
       call rnd_r(rr,rndstate)
       if(rr*fmax <= q*(1d0-exp(-tau0/q))) return
       call rnd_r(rr,rndstate)
       mucmf = 1d0-2d0*rr
    enddo
  end subroutine line_mu
!
end module nxemitmod
! vim: fdm=marker
