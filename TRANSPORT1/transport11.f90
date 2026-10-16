!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
pure subroutine transport11(ptcl,ptcl2,vx,vy,vz,rndstate,&
  eraddens,eamp,jrad,totevelo,ierr,jbk,jbval)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  use dopplermod
  use comptonmod
  use nlteatommod, only:nxa_lnu,nxa_lwl
  use nltextrmod
  use bbsobmod, only:bbs_on,grd_capsc
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(inout) :: vx,vy,vz
  real*8,intent(out) :: eraddens, eamp
  real*8,intent(out) :: jrad
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
  integer,intent(out) :: jbk      !NLTE line crossed (resonance estimator), 0 if none
  real*8,intent(out) :: jbval     !its contribution e*c/(nu*|Q|)
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c

  logical :: lout
  integer :: ixnext
  real*8 :: r1, r2
  real*8 :: db, dcol, dthm, ddop, dsc
  real*8 :: darr(5)
  real*8 :: elabfact,elabfactold
  real*8 :: xold, muold
  real*8 :: help
!-- distance out of physical reach
  real*8 :: far
  real*8 :: emitlump
!-- interpolate velocity
  real*8 :: vhelp,vhelp1,vhelp2,help1,help2
  real*8 :: dummy

  integer,pointer :: ix,ic,ig
  integer,parameter :: iy=1, iz=1
  real*8,pointer :: x, mu, e, e0, wl, d
!-- anisotropic Thomson scattering
  real*8 :: mutemp,muhelp,munew1,munew2,B,C
!-- Doppler distance calculations
  real*8 :: wlhelp!,vxhelp,labfacthelp
  ! logical :: lsuplum
  real*8 :: shelp,xhelp,epstol,phi,dphi
  integer :: it
  integer :: idop !direction of the comoving group crossing
!-- thermal electron scattering
  real*8 :: mucmf,wlcmf,nuratio
!-- NLTE hydrogen: line resonances and macro-atom
  integer :: jlo,jhi,jres,kres,kind,idx,ilin0,ilsd0
  real*8 :: wlo_e,whi_e,qmu,tau,qmin
  logical :: lremit

!-- statement function
  integer :: l
  real*8 :: dx
  real*8 :: dvx, dvdr
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dvx(l) = grd_vxarr(l+1) - grd_vxarr(l)
  dvdr(l) = dvx(l) / dx(l)

  ix => ptcl2%ix
  ic => ptcl2%ic
  ig => ptcl2%ig
  d => ptcl2%dist
  x => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

  dummy = vy
  dummy = vz

!-- no error by default
  ierr = 0
!-- init
  eraddens = 0d0
  eamp = 0d0
  jrad = 0d0
  jbk = 0
  jbval = 0d0
  lremit = .false.
!-- the packet sits on an NLTE line only right after a resonance event
  ilin0 = ptcl2%iline
  ilsd0 = ptcl2%ilside
  ptcl2%iline = 0
  ptcl2%ilside = 0
  jlo = 0
  jhi = 0
  jres = 0
!
!-- calculating initial transformation factors
  elabfact = 1d0 - mu*vx*cinv
  elabfactold = elabfact

!-- distance longer than logest possible photon path
  far = 2d0*abs(sqrt(grd_xarr(ix+1)**2-grd_xarr(ix)**2))

!
!-- boundary distances
  if(ix==1 .or. mu>=-sqrt(1d0-(grd_xarr(ix)/x)**2)) then
!-- outer boundary
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
     ixnext = ix+1
  else
!-- inner boundary
     db = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
     ixnext = ix-1
  endif
!-- sanity check
  if(db/=db) then
!    stop 'transport11: db/=db'
     ierr = 1
     return
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = -log(r1)/(elabfact*grd_sig(ic))
  else
     dthm = far
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = far
  else
     call rnd_r(r1,rndstate)
     dcol = -log(r1)/(elabfact*grd_cap(ig,ic))
  endif
!
!-- line scattering distance (lines of the LTE elements)
  dsc = far
  if(bbs_on) then
     if(grd_capsc(ig,ic)>0.) then
        call rnd_r(r1,rndstate)
        dsc = -log(r1)/(elabfact*grd_capsc(ig,ic))
     endif
  endif
!
!-- Doppler shift distance
  if(trn_dopexact) then
!-- exact distance to the first comoving group boundary or NLTE line along the
!-- path, for the cell's linear v(r) of any slope (decelerating and
!-- non-monotonic flows have blueshifts and repeated resonances)
     wlo_e = grp_wl(ig)
     whi_e = grp_wl(ig+1)
     if(nxt_on) then
        call nxt_neighbors(wl/elabfact,ilin0,ilsd0,jlo,jhi)
        if(jlo>0) then
           if(nxt_wlact(jlo)>wlo_e) then
              wlo_e = nxt_wlact(jlo)
           else
              jlo = 0
           endif
        endif
        if(jhi>0) then
           if(nxt_wlact(jhi)<whi_e) then
              whi_e = nxt_wlact(jhi)
           else
              jhi = 0
           endif
        endif
     endif
     call doppler_crossing(x,mu,grd_xarr(ix),grd_xarr(ix+1),grd_vxarr(ix), &
          grd_vxarr(ix+1),wl,wlo_e,whi_e,db,pc_c,ddop,idop)
     if(idop==0) ddop = far
     if(idop>0 .and. jhi>0) jres = jhi
     if(idop<0 .and. jlo>0) jres = jlo
     goto 66
  endif
!-- legacy: \tilde{v} frame
  vhelp = grd_vxarr(ix) - dvdr(ix)*grd_xarr(ix)
  muhelp = (mu - vx*cinv) / (1 - mu*vx*cinv)
  wlhelp = (wl/elabfact) / (1 + (muhelp*(vx-vhelp)*cinv))
  ! lsuplum = .false.
  if(abs(vx - vhelp)/pc_c.ge.0.1d0) then  ! check for superluminal flow
     ! lsuplum = .true.
     ! ddop = far
     goto 67
  ! elseif(dvdr(ix)>0.and.ig<grp_ng) then !redshift
  !    ddop = pc_c*wlhelp/(wl*dvdr(ix))*(elabfact - wl*grp_wlinv(ig+1))
  ! elseif(dvdr(ix)<0.and.ig>1) then ! blueshift
  !    ddop = pc_c*wlhelp/(wl*dvdr(ix))*(elabfact - wl*grp_wlinv(ig))
  ! else
  !    ddop = far
  ! endif
!-- SuperNu method modified for non-homologous outflows
  elseif(dvdr(ix)>0.and.ig<grp_ng) then !redshift
     ddop = 3d0*pc_c/grd_divv(ix)*(elabfact - wl*grp_wlinv(ig+1))
  elseif(dvdr(ix)<0.and.ig>1) then ! blueshift
     ddop = 3d0*pc_c/grd_divv(ix)*(elabfact - wl*grp_wlinv(ig))
  else
     ddop = far
  endif

  if(ddop<0d0) then
     ddop = far
  endif
  goto 66

!-- Newton-Raphson Method (Wagle et al. 2023, ApJ, 953, 132)
67 continue
  epstol = 1d-3 ! tolerance for convergence
  it = 0
  shelp = db ! initial value (s_0)
  !-- new x & mu
  xhelp = sqrt(x**2 + shelp**2 + 2d0*shelp*x*mu)
  muhelp = (mu*x + shelp) / xhelp
  !-- phi
  phi = dvdr(ix)*mu*x - pc_c*(1-wl*grp_wlinv(ig+1))+&
        vhelp*muhelp + dvdr(ix)*shelp
  !-- derivative of phi
  dphi = dvdr(ix) + vhelp*(1 - muhelp**2)/xhelp
  !-- next value, s_k
  ddop = shelp - phi / dphi
  !-- iterate to find ddop
  do while((abs(ddop-shelp).gt.ddop*epstol).and.it.lt.100)
    it = it + 1
    shelp = ddop ! store previous value
    !-- new x & mu
    xhelp = sqrt(x**2 + shelp**2 + 2d0*shelp*x*mu)
    muhelp = (mu*x + shelp) / xhelp
    !-- phi
    if(dvdr(ix)>0.and.ig<grp_ng) then
      phi = dvdr(ix)*mu*x - pc_c*(1-wl*grp_wlinv(ig+1))+&
            vhelp*muhelp + dvdr(ix)*shelp
    elseif(dvdr(ix)<0.and.ig>1) then
      phi = dvdr(ix)*mu*x - pc_c*(1-wl*grp_wlinv(ig))+&
            vhelp*muhelp + dvdr(ix)*shelp
    else
      ddop = far
      exit
    endif
    !-- derivative of phi
    dphi = dvdr(ix) + vhelp*(1 - muhelp**2)/xhelp
    !-- update ddop
    ddop = shelp - phi / dphi
  enddo

  if(ddop.le.0d0.or.it.eq.100) then
     ddop = far
  endif

66 continue
!
!--finding minimum distance
  darr = [dcol,dthm,ddop,db,dsc]
  ptcl2%idist = minloc(darr,dim=1)
  d = minval(darr)
  if(any(darr/=darr)) then
     ierr = 3
     return
  endif
  if(d<0d0) then
     ierr = 4
     return
  endif

!-- updating position, angle
  xold = x
  x = sqrt(x**2 + d**2 + 2d0*d*x*mu)
!-- interpolate velocity
  vhelp1 = grd_vxarr(ix)
  vhelp2 = grd_vxarr(ix+1)
  help1 = grd_xarr(ix)
  help2 = grd_xarr(ix+1)
  vx = (vhelp1*(help2-x)+vhelp2*(x-help1))/dx(ix)
  muold = mu
  if(x==0d0) then
     mu = 1d0
  else
     mu = (xold*mu+d)/x
  endif

!-- tallying energy densities
  eraddens = e*elabfact**2*d*cinv

!-- updating raditation intensity
  jrad = pc_c*eraddens
!
!-- updating transformation factors
  elabfact = 1d0 - mu*vx*cinv

!
!-- instantiate new particle at inner boundary
  if(x < grd_xarr(1)) then
     if(trn_innerbc=='absb') then
!-- absorbed by the boundary
        ptcl2%stat = 'dead'
        return
     endif
     x = grd_xarr(1)
     vx = grd_vxarr(1)
     ix = 1
     ic = grd_icell(ix,iy,iz)
     call rnd_r(r1,rndstate)
     if(trn_innerbc=='thrm') then
!-- thermalized: Lambertian direction and the boundary blackbody spectrum,
!-- comoving energy conserved
        e = e*elabfact
        e0 = e0*elabfact
        mucmf = sqrt(r1)
        mu = (mucmf+vx*cinv)/(1d0+vx*mucmf*cinv)
        help = trn_tinner
        if(help<=0d0) help = grd_radtemp(1)
        call planck_wl(help,rndstate,wlcmf)
        wlcmf = min(max(wlcmf,grp_wl(1)),grp_wl(grp_ng+1)*(1d0-1d-12))
        ig = binsrch(wlcmf,grp_wl,grp_ng+1,.false.)
        wl = wlcmf*(1d0-mu*vx*cinv)
        e = e/(1d0-mu*vx*cinv)
        e0 = e0/(1d0-mu*vx*cinv)
     else
     !-- forward scattering 0 <= mu <=1
        mu = 1d0 - r1
     endif
     return
  endif
!

!-- common manipulations for collisions
  nuratio = 1d0
  if(d==dthm.and.trn_thm_thermal) then
!-- scattering off a thermal electron: new comoving direction and frequency
     mucmf = (mu-vx*cinv)/(1d0-mu*vx*cinv)
     wlcmf = wl/elabfact
     call thermal_thomson(mucmf,wlcmf,1d0/grd_tempinv(ic),trn_thm_aniso, &
          rndstate,mu,nuratio)
     mu=(mu+vx*cinv)/(1d0+vx*mu*cinv)
  elseif(d==dthm.or.d==dcol) then
     call rnd_r(r1,rndstate)
!-- anisotropic scattering (experimental)
     if(d==dthm.and.trn_thm_aniso) then
       mutemp = ((4d0*r1-2d0+SQRT(1d0+16d0*(r1-0.5)**2d0))**(1d0/3d0) - &
             (SQRT(1d0+16d0*(r1-0.5)**2d0)-(4d0-2d0*r1))**(1d0/3d0))
       muhelp = mutemp ! need to verify
       B = 2*mutemp*muhelp
       C = mutemp**2 + muhelp**2 - 1
       munew1 = (-B + SQRT(B**2 - 4*C))/2
       munew2 = (-B - SQRT(B**2 - 4*C))/2
       !-- randomizing the angle sampling
       if(munew1.le.1.and.munew1.ge.-1.and.&
         munew2.le.1.and.munew2.ge.-1) then
         if(r1<0.5) then
           mu = munew1
         else
           mu = munew2
         endif
       elseif(munew1.le.1.and.munew1.ge.-1) then
         mu = munew1
       elseif(munew2.le.1.and.munew2.ge.-1) then
         mu = munew2
       else
         mu = 1d0-2d0*r1
       endif
!-- isotropic scattering
     else
       mu = 1d0-2d0*r1
     endif
!-- checking velocity dependence
     mu=(mu+vx*cinv)/(1d0+vx*mu*cinv)
!-- outer domain boundary check
  elseif(d==db) then
     lout = mu>=0d0.and.ix==grd_nx
     if(lout) then
!-- ending particle
        ptcl2%stat = 'flux'
        return
     endif
  endif

!
!-- Thomson scatter
  if(d == dthm) then
!-- lab wavelength (comoving frequency changed by nuratio for thermal electrons)
     wl = wl*(1d0-mu*vx*cinv)/(elabfact*nuratio)
     help = elabfact*nuratio/(1d0-mu*vx*cinv)
!-- velocity effects accounting (bulk flow and electron thermal motion)
     totevelo=totevelo+e*(1d0-help)
!-- energy weight
     e = e*help
     e0 = e0*help
     if(nuratio/=1d0) then
!-- new comoving group
        ig = binsrch(wl/(1d0-mu*vx*cinv),grp_wl,grp_ng+1,.false.)
        if(((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dx(ix) >= trn_tauddmc) &
             .and. .not.in_puretran) then
           ptcl2%itype = 2
           totevelo = totevelo+e*vx*mu*cinv
           e = e*(1d0-vx*mu*cinv)
           e0 = e0*(1d0-vx*mu*cinv)
           wl = wl/(1d0-vx*mu*cinv)
        endif
     endif

!
!-- line scattering: isotropic and coherent in the comoving frame
  elseif(d == dsc) then
     call rnd_r(r1,rndstate)
     mucmf = 1d0-2d0*r1
     mu = (mucmf+vx*cinv)/(1d0+vx*mucmf*cinv)
     wl = wl*(1d0-mu*vx*cinv)/elabfact
     help = elabfact/(1d0-mu*vx*cinv)
     totevelo = totevelo+e*(1d0-help)
     e = e*help
     e0 = e0*help

!
!-- effective scattering
  elseif(d == dcol) then
     if(nxt_on) then
!-- NLTE bound-free absorption activates the macro-atom, the rest the thermal
!-- pool; X-rays make fast photoelectrons (heat, H/He ionization and excitation)
        call rnd_r(r1,rndstate)
        if(nxt_ntemin>0d0 .and. pc_h*pc_c*elabfact/wl>nxt_ntemin) then
           call nxt_nonthermal(ic,rndstate,kind,idx,wlcmf)
        elseif(r1*grd_cap(ig,ic) < grd_capnbf(ig,ic)) then
           call nxt_bfabs(ic,ig,pc_c*elabfact/wl,rndstate,kind,idx,wlcmf)
        else
           call nxt_kpacket(ic,rndstate,kind,idx,wlcmf)
        endif
        lremit = .true.
!-- a line photon leaves its Sobolev resonance along directions of large |Q(mu)|
        if(kind==em_line) then
           call line_mu(nxt_kq(idx,ic),dvdr(ix),vx/x,rndstate,mucmf)
           mu = (mucmf+vx*cinv)/(1d0+vx*mucmf*cinv)
        endif
     endif
     call rnd_r(r1,rndstate)
!-- transforming to lab
     help = elabfact/(1d0-mu*vx*cinv)
!-- velocity effects accounting
     totevelo = totevelo+e*(1d0-help)
!-- energy weight
     e = e*help
     e0 = e0*help
!-- redistributing wavelength
     if(in_nlte) then !NLTE
        emitlump = grd_opaclump(9,ic)/grd_capemitgrey(ic)
     else !LTE
        emitlump = grd_opaclump(8,ic)/grd_capgrey(ic)
     endif
     if(nxt_on) then
!-- outcome chosen above
     elseif(grp_ng==1) then
        !-- do nothing
     else
        call rnd_r(r1,rndstate)
        ig = emitgroup(r1,ic)
!-- checking if DDMC in new group
        if(((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dx(ix) >= trn_tauddmc) &
            .and. .not.in_puretran) ptcl2%itype = 2
     endif
!-- checking for DDMC in new group
     if(ptcl2%itype==2) then
!-- transforming to cmf
!-- velocity effects accounting
        totevelo = totevelo+e*vx*mu*cinv
!-- energy weight
        e = e*(1d0-vx*mu*cinv)
        e0 = e0*(1d0-vx*mu*cinv)
        !wl = 0d0 !workaround ifort 13.1.3 bug
     elseif(.not.lremit) then
!-- uniformly in new group
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- converting comoving wavelength to lab frame wavelength
        wl = wl*(1d0-vx*mu*cinv)
     endif
     if(lremit) then
        call nx_reemit(kind,idx,wlcmf,ig,wl,e,e0,mu,vx,x,ix,ic,ptcl2,totevelo,rndstate)
        ptcl%em = nxt_emcode(kind,idx)
     endif

!
!-- outer radial bound
  elseif(d==db .and. ixnext>ix) then

     l = grd_icell(ix+1,iy,iz)
     if((grd_sig(l)+grd_cap(ig,l)+grd_capsc(ig,l))*dx(ix+1)<trn_tauddmc &
          .or. in_puretran) then
!-- IMC in adjacent cell
        x = grd_xarr(ix+1)
        vx = grd_vxarr(ix+1)
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        mu = (mu-vx*cinv)/(1d0-vx*mu*cinv)
        help= (grd_cap(ig,l)+grd_capsc(ig,l)+grd_sig(l))*dx(ix+1)
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5*abs(mu))) then
           ptcl2%itype = 2
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-elabfact)
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
!-- update
           ix = ix+1
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
           mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
           x = grd_xarr(ix+1)
           vx = grd_vxarr(ix+1)
        endif
     endif

!
!-- inner radial bound
  elseif(d==db) then

     l = grd_icell(ix-1,iy,iz)
     if((grd_sig(l)+grd_cap(ig,l)+grd_capsc(ig,l))*dx(ix-1) &
           < trn_tauddmc .or. &
          in_puretran) then
!-- IMC in adjacent cell
        x = grd_xarr(ix)
        vx = grd_vxarr(ix)
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
!-- transforming x-cosine to cmf
        mu = (mu-vx*cinv)/(1d0-vx*mu*cinv)
!-- amplification factor
        if(.not.trn_noampfact .and. mu<0d0) then
            help = 1d0/abs(mu)
            help = min(100d0, help) !-- truncate singularity
!-- velocity effects accounting
            totevelo = totevelo-e*2d0*(0.55d0*help-1.25d0*abs(mu))*vx*cinv
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*vx*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0 = e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*vx*cinv)
              e = e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*vx*cinv)
        endif
        help = (grd_cap(ig,l)+grd_capsc(ig,l)+grd_sig(l))*dx(ix-1)
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl2%itype = 2
!-- velocity effects accounting
           totevelo = totevelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
           ix = ix-1
           ic = grd_icell(ix,iy,iz)
        else
!-- resampling direction
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
!-- transforming mu to lab
           mu=(mu+vx*cinv)/(1d0+vx*mu*cinv)
           x = grd_xarr(ix)
           vx = grd_vxarr(ix)
        endif
     endif

!
!-- Doppler shift
  elseif(d == ddop .and. jres>0) then
!-- NLTE line resonance (Sobolev): projected velocity gradient at this point
     kres = nxt_act(jres)
     mucmf = (mu-vx*cinv)/(1d0-mu*vx*cinv)
     qmu = mucmf**2*dvdr(ix) + (1d0-mucmf**2)*vx/x
     qmin = 1d-3*max(abs(dvdr(ix)),abs(vx/x))
     qmu = max(abs(qmu),qmin)
     tau = nxt_kq(kres,ic)/qmu
!-- resonance estimator of the incident intensity
     jbk = kres
     jbval = e*elabfact*pc_c/(nxa_lnu(kres)*qmu)
     call rnd_r(r1,rndstate)
     if(tau>0d0 .and. r1<1d0-exp(-tau)) then
!-- absorbed: macro-atom in the upper level
        if(nxt_tau0>0d0) then
!-- benchmark: pure scattering in the same line
           kind = em_line
           idx = kres
           wlcmf = nxa_lwl(kres)
        else
           call nxt_activate_line(ic,kres,rndstate,kind,idx,wlcmf)
        endif
!-- line photons escape along directions of large |Q(mu)|, continuum isotropically
        if(kind==em_line) then
           call line_mu(nxt_kq(idx,ic),dvdr(ix),vx/x,rndstate,mucmf)
        else
           call rnd_r(r1,rndstate)
           mucmf = 1d0-2d0*r1
        endif
        mu = (mucmf+vx*cinv)/(1d0+vx*mucmf*cinv)
        help = elabfact/(1d0-mu*vx*cinv)
        totevelo = totevelo+e*(1d0-help)
        e = e*help
        e0 = e0*help
        call nx_reemit(kind,idx,wlcmf,ig,wl,e,e0,mu,vx,x,ix,ic,ptcl2,totevelo,rndstate)
        ptcl%em = nxt_emcode(kind,idx)
     else
!-- passes the resonance: now on the far side of the line
        ptcl2%iline = jres
        ptcl2%ilside = idop
     endif
  elseif(d == ddop) then
     if(trn_dopexact) then
!-- the lab-frame wavelength does not change in free flight; only the
!-- comoving group does
        if(idop>0) then
           if(ig<grp_ng) then
              ig = ig+1
           else
!-- resampling wavelength in highest group
              call rnd_r(r1,rndstate)
              wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
              wl = wl*elabfact
           endif
        else
           if(ig>1) then
              ig = ig-1
           else
!-- resampling wavelength in lowest group
              call rnd_r(r1,rndstate)
              wl=1d0/(r1*grp_wlinv(1) + (1d0-r1)*grp_wlinv(2))
              wl = wl*elabfact
           endif
        endif
     elseif(dvdr(ix).gt.0) then ! redshit
        if(ig<grp_ng) then
!-- shifting group
          ig = ig+1
          wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfactold
        else
!-- resampling wavelength in highest group
          call rnd_r(r1,rndstate)
          wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
          wl = wl*elabfact
        endif
     elseif(dvdr(ix)<0) then ! blueshift
       if(ig>1) then
!-- shifting group
          wl = (grp_wl(ig)+1d-6*(grp_wl(ig-1)-grp_wl(ig)))*elabfactold
          ig = ig-1
        else
!-- resampling wavelength in lowest group
          call rnd_r(r1,rndstate)
          wl=1d0/(r1*grp_wlinv(1) + (1d0-r1)*grp_wlinv(2))
          wl = wl*elabfact
        endif
     endif
!-- check if ddmc region
     if (((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dx(ix) &
           >= trn_tauddmc) .and. .not.in_puretran) then
        ptcl2%itype = 2
!-- velocity effects accounting
        totevelo=totevelo+e*vx*mu*cinv
!
        e = e*elabfact
        e0 = e0*elabfact
        wl = wl/elabfact
     endif
  else
!    stop 'transport11: invalid distance'
     ierr = 17
     return
  endif

  return

contains

  pure subroutine planck_wl(temp,rs,wlp)
!-- wavelength drawn from a blackbody (Carter & Cashwell 1975: x = h nu/kT)
    real*8,intent(in) :: temp
    type(rnd_t),intent(inout) :: rs
    real*8,intent(out) :: wlp
    real*8 :: r0,s,u(4)
    integer :: l2,j2
    call rnd_r(r0,rs)
    r0 = r0*pc_pi**4/90d0
    s = 0d0
    do l2=1,10000
       s = s + 1d0/dble(l2)**4
       if(s >= r0) exit
    enddo
    do j2=1,4
       call rnd_r(u(j2),rs)
    enddo
    wlp = pc_h*pc_c/(pc_kb*temp*(-log(max(u(1)*u(2)*u(3)*u(4),1d-300))/dble(l2)))
  end subroutine planck_wl

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
    real*8 :: wlc,qn,rr,dvdrl
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
       if(wlc<grp_wl(1) .or. wlc>=grp_wl(grp_ng+1)) then
!-- outside the group grid (e.g. mid-IR fine-structure lines): the medium
!-- is not described there, the photon escapes (flux tally: outside the
!-- bins unless the flux grid is wider)
          wl = wlc*(1d0-vx*mu*cinv)
          ptcl2%stat = 'flux'
          return
       elseif(kind==em_line .and. nxt_isrt(idx)>0) then
!-- the packet sits on the (active) line it was emitted in
          ptcl2%iline = nxt_isrt(idx)
          dvdrl = (grd_vxarr(ix+1)-grd_vxarr(ix))/(grd_xarr(ix+1)-grd_xarr(ix))
          qn = ((mu-vx*cinv)/(1d0-mu*vx*cinv))**2
          qn = qn*dvdrl + (1d0-qn)*vx/x
          ptcl2%ilside = merge(1,-1,qn>=0d0)
       endif
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
  end subroutine nx_reemit

! !-- updating the co-moving group for superluminal flow for \tilde{v} frame
!   if(lsuplum.and.ptcl2%itype == 1) then
! !-- update velocity coordinate
!     vhelp1 = grd_vxarr(ix)
!     vhelp2 = grd_vxarr(ix+1)
!     help1 = grd_xarr(ix)
!     help2 = grd_xarr(ix+1)
!     vxhelp = (vhelp1*(help2-x)+vhelp2*(x-help1))/(help2-help1)
!     labfacthelp = 1d0-vxhelp*mu/pc_c
!     !-- comoving wavelength
!     wlhelp = wl/labfacthelp
!     !-- comoving group index
!     ig = binsrch(wlhelp,grp_wl,grp_ng+1,.false.)
! !-- check for ddmc again for new ig
!     if (((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dx(ix) &
!            >= trn_tauddmc) .and. .not.in_puretran) then
!       ptcl2%itype = 2
! !-- velocity effects accounting
!       e = e*labfacthelp
!       e0 = e0*labfacthelp
!       wl = wlhelp
!     endif
!   endif

end subroutine transport11
! vim: fdm=marker
