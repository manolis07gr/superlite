!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
pure subroutine transport11(ptcl,ptcl2,vx,vy,vz,rndstate,&
  eraddens,eamp,jrad,totevelo,ierr,jbk,jbval,stl,stk,swl)

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
  use nxemitmod
  use bbsobmod, only:bbs_on,grd_capsc,grd_fray
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
  real*8,intent(out) :: stl(4)    !flight in a slow cell (nxt_dopmc): comoving wavelength, comoving energy times length (0: none), projected velocity gradient, comoving length
  integer,intent(out) :: stk      !the packet's own line on that flight (nxt_static_tally)
  real*8,intent(out) :: swl(3)    !the flight for the window estimator of the line incident intensity (nxt_jbwin): comoving wavelength at its start and end, comoving energy times length (0: none)
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
  real*8 :: db, dcol, dthm, ddop, dsc, dlin, ddrift
  real*8 :: qsl, mus
  real*8 :: darr(7)
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
!-- Doppler-profile lines in slow cells
  logical :: lslow
  real*8 :: kapl,wlc0,qmx

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
  stl = 0d0
  stk = 0
  swl = 0d0
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

!-- (nxt_esrates: a photon just scattered out of a line core by an electron
!-- lands in the line profiles as after a thermal electron scattering)
  if(ptcl2%estest) then
     ptcl2%estest = .false.
     if(nxt_on .and. nxt_esrates) then
        mucmf = (mu-vx*cinv)/(1d0-mu*vx*cinv)
        qmu = mucmf**2*dvdr(ix) + (1d0-mucmf**2)*vx/x
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        call nxt_esfar(ic,wl/elabfact,qmu,vx,1d0/grd_tempinv(ic),grd_sig(ic),r1,r2,kres,ptcl2%kskip)
        if(kres>0) then
           ptcl2%kfa = kres
           ptcl2%efa = e*elabfact
           d = 0d0
           call nx_lineabs(kres,ic,ix,elabfact,dvdr(ix),x,vx,wl,ig,e,e0,mu,totevelo, &
                ptcl,ptcl2,rndstate)
           return
        endif
     endif
  endif

!-- distance longer than logest possible photon path
  far = 2d0*abs(sqrt(grd_xarr(ix+1)**2-grd_xarr(ix)**2))

!
!-- boundary distances (the innermost cell of a grid with a hollow core,
!-- grd_xarr(1) > 0, has the inner boundary as a surface: before, its
!-- inward rays crossed the core as part of the cell, with the cell's
!-- opacity and their path in its J)
  if((ix==1 .and. grd_xarr(1)<=0d0) .or. &
       mu>=-sqrt(max(0d0,1d0-(grd_xarr(ix)/x)**2))) then
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
!-- Doppler-profile NLTE lines in a slow cell (nxt_dopmc): static opacity
!-- of the lines near the comoving wavelength (from Z = nxt_dopz) but the
!-- packet's own, and the flight ends after a drift of half a Doppler width
!-- to evaluate it again
  lslow = .false.
  dlin = far
  ddrift = far
  wlc0 = wl/elabfact
  if(nxt_on .and. nxt_dopmc>0d0) then
     if(nxt_slow(ic) .and. nxt_anystat) then
        lslow = .true.
        call nxt_static_line(ic,wlc0,ptcl%ilem,-1d0,kapl,kres)
        if(kapl>0d0) then
           call rnd_r(r1,rndstate)
           dlin = -log(r1)/(elabfact*kapl)
        endif
        qmx = max(abs(dvdr(ix)),abs(vx/max(x,1d-300)))
        if(qmx>0d0) ddrift = 0.5d0*nxt_vdlmin(ic)/qmx
!-- projected velocity gradient of the flight: the side of a line its
!-- comoving frequency drifts in from (nxt_static_tally)
        mus = (mu-vx*cinv)/(1d0-mu*vx*cinv)
        qsl = mus**2*dvdr(ix) + (1d0-mus**2)*vx/max(x,1d-300)
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
!-- (a slow cell keeps the resonances of the lines without static profiles)
        if(lslow) call nxt_sobolev_neighbors(jlo,jhi)
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
  darr = [dcol,dthm,ddop,db,dsc,dlin,ddrift]
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
!-- (in_rayleigh = 'lymes', benchmark: a scattering by H I is redistributed
!-- by the thermal electrons, as CMFGEN's iterations with incoherent electron
!-- scattering treat their Rayleigh opacity: the event becomes a Thomson
!-- scattering, by its share of the scattering opacity of the group)
  if(trn_raythm .and. d==dsc) then
     call rnd_r(r1,rndstate)
     if(r1 < grd_fray(ig,ic)) dthm = d
  endif

!-- the flight for the incident intensity of the lines of a slow cell
  if(lslow) then
     stl(1) = wlc0
     stl(2) = e*elabfact*d
     stl(3) = qsl
     stl(4) = elabfact*d
     stk = ptcl%ilem
  endif
!
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
!-- the flight for the window estimator of the line incident intensity
  if(nxt_on .and. nxt_jbwin>0d0) then
     swl(1) = wlc0
     swl(2) = wl/elabfact
     swl(3) = jrad
  endif

!
!-- instantiate new particle at inner boundary
  if(x < grd_xarr(1) .or. (d==db .and. ixnext<1)) then
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
!-- diffuse reflection: Lambertian direction, as the surface source emits
!-- (a uniform mu, as before, put an isotropic field's reflected packets on
!-- near-tangent rays and raised the J of a thin innermost cell)
        mu = sqrt(r1)
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
!-- (nxt_esrates: scattered into the profile of a line, where the point
!-- resonance would absorb it only if it drifted through the centre before
!-- the next scattering: absorbed in the profile, as the electron scattering
!-- of the rates brings photons into the cores; the lines ahead whose whole
!-- profile was tested are then passed without absorption, ptcl2%kskip)
     ptcl2%kskip = 0
     if(nuratio/=1d0 .and. nxt_on .and. nxt_esrates) then
        mucmf = (mu-vx*cinv)/(1d0-mu*vx*cinv)
        qmu = mucmf**2*dvdr(ix) + (1d0-mucmf**2)*vx/x
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        call nxt_esfar(ic,wl/(1d0-mu*vx*cinv),qmu,vx,1d0/grd_tempinv(ic),grd_sig(ic),r1,r2,kres, &
             ptcl2%kskip)
        if(kres>0) then
           ptcl2%kfa = kres
           ptcl2%efa = e*(1d0-mu*vx*cinv)
           call nx_lineabs(kres,ic,ix,1d0-mu*vx*cinv,dvdr(ix),x,vx,wl,ig,e,e0,mu,totevelo, &
                ptcl,ptcl2,rndstate)
           return
        endif
     endif
     if(nuratio/=1d0) then
!-- (a line photon scattered to the blue of its line will cross it again)
        if(ptcl%ilem>0) then
           if(wl/(1d0-mu*vx*cinv) < nxa_lwl(ptcl%ilem)) ptcl%ilem = 0
        endif
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
!-- a line photon leaves its Sobolev resonance along directions of large
!-- |Q(mu)|, or scattered by an electron in the resonance zone
        if(kind==em_line) then
           call nx_linees(kind,idx,ic,rndstate,mucmf,wlcmf)
           if(kind==em_line) call line_mu(nxt_kq(idx,ic),dvdr(ix),vx/x,rndstate,mucmf)
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
        ptcl%icem = 0
        if(nxt_emloc(ptcl%em)) ptcl%icem = ic
        ptcl%ilem = nx_ilem(kind,idx,wlcmf)
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
!-- absorbed in the static profile of a line (slow cell): the line by its
!-- share of the opacity at the comoving wavelength of the flight, the
!-- macro-atom in its upper level; a line photon leaves along the
!-- directions of its Sobolev escape, as from a resonance
  elseif(d == dlin) then
     call rnd_r(r1,rndstate)
     call nxt_static_line(ic,wlc0,ptcl%ilem,r1,kapl,kres)
     mucmf = (mu-vx*cinv)/(1d0-mu*vx*cinv)
     if(kres>0) then
        call nxt_activate_line(ic,kres,rndstate,kind,idx,wlcmf)
     else
        call nxt_kpacket(ic,rndstate,kind,idx,wlcmf)
     endif
     if(kind==em_line) call nx_linees(kind,idx,ic,rndstate,mucmf,wlcmf)
     if(kind==em_line) then
        call line_mu(nxt_kq(idx,ic),dvdr(ix),vx/x,rndstate,mucmf)
     elseif(kind/=em_lines) then
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
     ptcl%icem = 0
     if(nxt_emloc(ptcl%em)) ptcl%icem = ic
     ptcl%ilem = nx_ilem(kind,idx,wlcmf)
!
!-- (end of a drift in a slow cell: the static line opacity is evaluated
!-- again on the next flight)
  elseif(d == ddrift) then
     continue
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
!-- (nxt_esrates: the line's profile was tested for the whole path at the
!-- packet's last electron scattering, nxt_esfar)
     if(any(ptcl2%kskip==kres)) then
        where(ptcl2%kskip==kres) ptcl2%kskip = 0
        jbk = 0
        jbval = 0d0
        tau = 0d0
     endif
     call rnd_r(r1,rndstate)
     if(tau>0d0 .and. r1<1d0-exp(-tau)) then
!-- absorbed: macro-atom in the upper level
        ptcl2%kab = kres
        ptcl2%eab = e*elabfact
        if(nxt_tau0>0d0) then
!-- benchmark: pure scattering in the same line
           kind = em_line
           idx = kres
           wlcmf = nxa_lwl(kres)
        else
           call nxt_activate_line(ic,kres,rndstate,kind,idx,wlcmf)
        endif
!-- line photons escape along directions of large |Q(mu)| (or are
!-- scattered by an electron in the resonance zone), continuum isotropically
        if(kind==em_line) call nx_linees(kind,idx,ic,rndstate,mucmf,wlcmf)
        if(kind==em_line) then
           call line_mu(nxt_kq(idx,ic),dvdr(ix),vx/x,rndstate,mucmf)
        elseif(kind/=em_lines) then
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
        ptcl%icem = 0
        if(nxt_emloc(ptcl%em)) ptcl%icem = ic
        ptcl%ilem = nx_ilem(kind,idx,wlcmf)
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
