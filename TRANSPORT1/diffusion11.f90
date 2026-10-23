!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
pure subroutine diffusion11(ptcl,ptcl2,vx,vy,vz,cache,rndstate,&
  eraddens,jrad,totevelo,ierr)

  use randommod
  use miscmod
  use groupmod
  use gridmod
  use bbsobmod, only:grd_capsc
  use physconstmod
  use particlemod
  use transportmod
  use fluxmod
  use inputparmod
  use nltextrmod, only:nxt_on,nxt_ntemin,nxt_kq,grd_capnbf,nxt_bfabs,nxt_kpacket, &
       nxt_nonthermal,nxt_emcode,nxt_emloc,em_line,em_lines
  use nxemitmod, only:nx_emitwl,nx_place,line_mu,nx_linees,nx_ilem
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(grp_t_cache),target,intent(inout) :: cache
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(inout) :: vx,vy,vz
  real*8,intent(out) :: eraddens
  real*8,intent(out) :: jrad
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
!
  integer :: iig, iiig
  logical :: lhelp
  real*8 :: r1, r2
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, pa, pdop
!-- lumped quantities -----------------------------------------

  real*8 :: emitlump, caplump, capemitlump, doplump
  real*8 :: specig
  real*8 :: mfphelp, ppl, ppr
  real*8 :: opacleak(2)
  real*8 :: probleak(2)
  real*8 :: resopacleak, resdopleak
  integer :: glump
  real*8,pointer :: capgreyinv
  real*8,pointer :: capemitgreyinv
  real*8,pointer :: speclump
  real*8 :: dist, help, wtb
  real*8 :: vhelp1,vhelp2,help1,help2
  real*8 :: dummy
!-- NLTE macro-atoms: absorption and re-emission
  integer :: kind,idx,ilo,ihi,imid
  real*8 :: nu,wlcmf,wlc,mucmf
  logical :: lout,lddmc

  integer,pointer :: ix, ic, ig
  integer,parameter :: iy=1, iz=1
  real*8,pointer :: x, mu, e, e0, wl
!-- statement function
  integer :: l
  real*8 :: dx,dx3
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx3(l) = grd_xarr(l+1)**3 - grd_xarr(l)**3

  ix => ptcl2%ix
  ic => ptcl2%ic
  ig => ptcl2%ig
  x => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!-- (a DDMC packet has no frequency history)
  ptcl%ilem = 0
  capgreyinv => cache%capgreyinv
  capemitgreyinv => cache%capemitgreyinv
  speclump => cache%speclump
!
  dummy = vy
  dummy = vz

!-- no error by default
  ierr = 0
!-- (nxt_esrates: line profiles are tested on IMC steps only)
  ptcl2%kskip = 0
  ptcl2%estest = .false.
!-- init
  eraddens = 0d0
  jrad = 0d0
!
  dist = dx(ix)

!
!-- update cache: the lumping of the cell (ddmc_lumps, once per iteration)
  if(ic/=cache%ic) then
     cache%ic = ic
!-- (mean absorption opacity over the DDMC source spectrum, ddmc_spectra)
     capgreyinv = max(1d0/grd_capgreyk(ic),0d0) !catch nans
     capemitgreyinv = max(1d0/grd_capemitgrey(ic),0d0)
     speclump = grd_opaclump(7,ic)
     cache%nlump = grd_nlump(ic)
     cache%emitlump = grd_vlump(1,ic)
     cache%caplump = grd_vlump(2,ic)
     cache%capemitlump = grd_vlump(3,ic)
     cache%doplump = grd_vlump(4,ic)
  endif !cache%ic /= ic

!
!-- in lump?
  if(grd_cap(ig,ic)*dist >= trn_taulump) then
     glump = cache%nlump
  else
     glump = 0
  endif
!
!-- sanity check
  if((grd_sig(ic) + grd_cap(ig,ic)+grd_capsc(ig,ic))*dist < trn_tauddmc) then
     ierr = 100
     return
  endif
!
!-- retrieve from cache
  if(glump>0) then
     emitlump = cache%emitlump
     caplump = cache%caplump
     capemitlump = cache%capemitlump
     doplump = cache%doplump
  else
!-- outside the lump
     if(in_nlte .and. .not.nxt_on) then !NLTE (H only)
       emitlump = specint0(grd_tempinv(ic),ig)* &
                  capemitgreyinv*grd_capemit(ig,ic)
     elseif(nxt_on) then ! NLTE macro-atoms: source spectrum (ddmc_spectra)
       emitlump = grd_spk(ig,ic)*capgreyinv*grd_cap(ig,ic)
     else ! LTE
       emitlump = specint0(grd_tempinv(ic),ig)* &
                  capgreyinv*grd_cap(ig,ic)
     endif
     caplump = grd_cap(ig,ic)
     capemitlump = grd_capemit(ig,ic)
     if(ig/=grp_ng.and.grd_divv(ix).ge.0) then ! redshift
       doplump = dopspeccalc(grd_tempinv(ic),ig)*grd_divv(ix) &
                   /(grd_planck(ig,ic)*3*pc_c)
     elseif(ig/=1.and.grd_divv(ix).lt.0) then ! blueshift
       doplump = -1*dopspeccalc(grd_tempinv(ic),ig-1)*grd_divv(ix) &
                   /(grd_planck(ig,ic)*3*pc_c)
     else
       doplump = 0d0
     endif
  endif
!
!-- calculate lumped values
  if(glump>0) then
!-- leakage opacities
     opacleak = grd_opaclump(1:2,ic)
!-- calculating unlumped values
  else

!-- inward
     if(ix/=1) l = grd_icell(ix-1,iy,iz)
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,l)+grd_capsc(ig,l)+ &
        grd_sig(l))*dx(ix-1)<trn_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_capsc(ig,ic)+grd_sig(ic))*dx(ix)
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)= 1.5d0*ppl*(grd_xarr(ix))**2/dx3(ix)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l)+grd_capsc(ig,l))*dx(ix-1))
        opacleak(1)=2.0d0*(grd_xarr(ix))**2/(mfphelp*dx3(ix))
     endif
!
!-- outward
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+grd_capsc(ig,l)+ &
           grd_sig(l))*dx(ix+1)<trn_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_capsc(ig,ic)+grd_sig(ic))*dx(ix)
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)=1.5d0*ppr*(grd_xarr(ix+1))**2/dx3(ix)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l)+grd_capsc(ig,l))*dx(ix+1))
        opacleak(2)=2.0d0*(grd_xarr(ix+1))**2/(mfphelp*dx3(ix))
     endif
  endif
!
!-------------------------------------------------------------
!

!-- calculate time to event (NLTE macro-atoms: every absorption is an
!-- event, unless all groups are lumped)
  if(nxt_on .and. glump<grp_ng) then
     denom = sum(opacleak) + caplump + doplump
  else
     denom = sum(opacleak) + (1d0-emitlump)*caplump + doplump
  endif
  denom = 1d0/denom

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)*denom*cinv)
  ddmct = tau

!
!-- calculating energy depostion and density
  eraddens  = e*ddmct
!

!-- updating radiation intensity
  jrad = pc_c*eraddens


!-- perform event
  call rnd_r(r1,rndstate)

!-- leak probability
  probleak = opacleak*denom

!-- absorption probability
  pa = 0d0

!-- Doppler shift
  pdop = doplump*denom

!-- doppler shift
  if (r1>=pa .and. r1<pa+pdop) then

     if(glump==0) then
        iiig = ig
     else
!-- sample group
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/doplump
        if(grd_divv(ix).ge.0) then ! redshift
          do iig=1,glump
             iiig = grd_glump(iig,ic)
             if(iiig == grp_ng) cycle
             if(grd_cap(iiig+1,ic)*dist >= trn_taulump) cycle
             specig = 0d0
             if(grd_planck(iiig,ic)>0d0) specig = grd_spk(iiig,ic)/grd_planck(iiig,ic)
             resdopleak = dopspeccalc(grd_tempinv(ic),iiig)&
                          *grd_divv(ix)/(3*pc_c)*specig
             denom2 = denom2+resdopleak*speclump*help
             if(denom2>r1) exit
          enddo
         else ! blueshift
           do iig=glump,1,-1
             iiig = grd_glump(iig,ic)
             if(iiig == 1) cycle
             if(grd_cap(iiig-1,ic)*dist >= trn_taulump) cycle
             specig = 0d0
             if(grd_planck(iiig,ic)>0d0) specig = grd_spk(iiig,ic)/grd_planck(iiig,ic)
             resdopleak = dopspeccalc(grd_tempinv(ic),iiig-1)&
                          *grd_divv(ix)/(3*pc_c)*specig
             denom2 = denom2-resdopleak*speclump*help
             if(denom2>r1) exit
           enddo
         endif
     endif

!-- reshift/blueshift particle in this group
     if(grd_divv(ix).ge.0) then ! redshift
       ig = iiig+1
       wl = grp_wl(ig)
       ig = min(ig,grp_ng)
     else ! blueshift
       ig = iiig-1
       wl = grp_wl(ig+1)
       ig = max(ig,1)
     endif

!-- method changes to IMC
     if((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dist < trn_tauddmc) then
        ptcl2%itype = 1
!-- direction sampled isotropically
        call rnd_r(r1,rndstate)
        mu = 1d0-2d0*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- determine velocity coordinates (currently uses linear interpolation)
!-- x
        vhelp1 = grd_vxarr(ix)
        vhelp2 = grd_vxarr(ix+1)
        help1 = grd_xarr(ix)
        help2 = grd_xarr(ix+1)
        vx = (vhelp1*(help2-x)+vhelp2*(x-help1))/dx(ix)
!
!-- velocity effects accounting
        mu = (mu+vx*cinv)/(1.0+vx*mu*cinv)
        wl = wl*(1.0-vx*mu*cinv)
        help = 1d0/(1.0-vx*mu*cinv)
        totevelo = totevelo+e*(1d0 - help)
        e = e*help
        e0 = e0*help
     endif

!-- left leakage sample
  elseif (r1>=pa+pdop .and. r1<pa+pdop+probleak(1)) then
     ptcl2%idist = -3

!-- checking if at inner bound
     if(ix==1) then
!       stop 'diffusion11: non-physical inward leakage'
!        ierr = 101
        return

!-- sample adjacent group (assumes aligned ig bounds)
     else

        l = grd_icell(ix-1,iy,iz)
        if(glump==0) then
           iiig = ig
        else
!-- sample group
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(1)
           do iig=1,glump
              iiig = grd_glump(iig,ic)
              specig = grd_spk(iiig,ic)
              !specig = specint0(grd_tempinv(ic),iiig)
!-- calculating resolved leakage opacities
              if((grd_cap(iiig,l)+grd_capsc(iiig,l)+ &
                   grd_sig(l))*dx(ix-1)<trn_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,ic)+grd_capsc(iiig,ic)+grd_sig(ic))*dx(ix)
                 ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*ppl*(grd_xarr(ix))**2/dx3(ix)
              else
!-- IMC interface
                 mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)+grd_capsc(iiig,ic))*dx(ix)+&
                      (grd_sig(l)+grd_cap(iiig,l)+grd_capsc(iiig,l))*dx(ix-1))
                 resopacleak = 2.0d0*(grd_xarr(ix))**2/ &
                      (mfphelp*dx3(ix))
              endif
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
        endif
!
!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l)+grd_capsc(iiig,l))*dx(ix-1) < trn_tauddmc) then
           ptcl2%itype = 1
!-- sampling wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!-- location set right bound of left cell
           x = grd_xarr(ix)
           vx = grd_vxarr(ix)
!-- particle angle sampled from isotropic b.c. inward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
!-- doppler and aberration corrections
           mu = (mu+vx*cinv)/(1.0+vx*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1.0-vx*mu*cinv)
           totevelo = totevelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1.0-vx*mu*cinv)
        endif
!
!-- update particle
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
        ig = iiig

     endif


!-- right leakage sample
  elseif (r1>=pa+pdop+probleak(1) .and. r1<pa+pdop+sum(probleak)) then
     ptcl2%idist = -4
!
!-- checking if at outer bound
     if(ix==grd_nx) then
        ptcl2%stat = 'flux'
!-- outbound luminosity tally
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = max(r1,r2)
        if(glump==0) then
        else
!-- sample group
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig=grd_glump(iig,ic)
              specig = grd_spk(iiig,ic)
              !specig = specint0(grd_tempinv(ic),iiig)
!-- calculating resolved leakage opacities
              mfphelp = (grd_cap(iiig,ic)+grd_capsc(iiig,ic)+grd_sig(ic))*dx(ix)
              ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 1.5d0*ppr*(grd_xarr(ix+1))**2/dx3(ix)
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
           ig = iiig
        endif
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(ig+1) + (1d0-r1)*grp_wlinv(ig))
!-- position
        x=grd_xarr(grd_nx+1)
        vx = grd_vxarr(grd_nx+1)
!-- changing from comoving frame to observer frame
        mu = (mu+vx*cinv)/(1d0+vx*mu*cinv)
        mu = min(mu,1d0)
        help = 1d0/(1d0-mu*vx*cinv)
!-- velocity effects accounting
        totevelo = totevelo+e*(1d0 - help)
        wl = wl/help
        e = e*help
        e0 = e0*help
        return
!
!
     else

        l = grd_icell(ix+1,iy,iz)
!-- sample adjacent group (assumes aligned ig bounds)
        if(glump==0) then
           iiig = ig
        else
!-- sample group
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig = grd_glump(iig,ic)
              specig = grd_spk(iiig,ic)
!-- calculating resolved leakage opacities
              if((grd_cap(iiig,l)+grd_capsc(iiig,l)+ &
                   grd_sig(l))*dx(ix+1)<trn_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,ic)+grd_capsc(iiig,ic)+grd_sig(ic))*dx(ix)
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*ppr*(grd_xarr(ix+1))**2/dx3(ix)
              else
!-- IMC interface
                 mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)+grd_capsc(iiig,ic))*dx(ix)+&
                      (grd_sig(l)+grd_cap(iiig,l)+grd_capsc(iiig,l))*dx(ix+1))
                 resopacleak = 2.0d0*(grd_xarr(ix+1))**2/ &
                      (mfphelp*dx3(ix))
              endif
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
        endif

!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l)+grd_capsc(iiig,l))*dx(ix+1)< trn_tauddmc) then
!
           ptcl2%itype = 1
!-- sampling wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!-- location set left bound of right cell
           x = grd_xarr(ix+1)
           vx = grd_vxarr(ix+1)
!-- particle angle sampled from isotropic b.c. outward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
!
!-- doppler and aberration corrections
           mu = (mu+vx*cinv)/(1.0+vx*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1.0-vx*mu*cinv)
           totevelo = totevelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1.0-vx*mu*cinv)
        endif
!
!-- update particle
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
        ig = iiig
!
     endif!


!-- effective scattering sample
  else
     ptcl2%idist = -2
     if(nxt_on) then
!-- NLTE: an absorption in the lump (group by its rate kappa J, J ~ S:
!-- grd_spk) or in the group, followed by the macro-atom or the thermal
!-- pool as on IMC steps (transport11), so that the energy of the thick
!-- continua cascades through the lines.  A re-emission in a DDMC group of
!-- the cell keeps the packet in DDMC (in the lump: no change); otherwise it
!-- leaves as an IMC packet from a random point of the cell.
        iiig = ig
        if(glump>0) then
           call rnd_r(r1,rndstate)
           r1 = r1*grd_lumpcum(glump,ic)
!-- (bisection: the first lump group whose cumulative rate exceeds r1)
           ilo = 1
           ihi = glump
           do while(ilo<ihi)
              imid = (ilo+ihi)/2
              if(grd_lumpcum(imid,ic) > r1) then
                 ihi = imid
              else
                 ilo = imid + 1
              endif
           enddo
           iiig = grd_glump(ilo,ic)
        endif
!-- comoving frequency uniform in the group
        call rnd_r(r1,rndstate)
        nu = pc_c*((1d0-r1)*grp_wlinv(iiig+1) + r1*grp_wlinv(iiig))
        call rnd_r(r1,rndstate)
        wtb = 1d0
        if(nxt_ntemin>0d0 .and. pc_h*nu>nxt_ntemin) then
           call nxt_nonthermal(ic,rndstate,kind,idx,wlcmf,wtb)
        elseif(r1*grd_cap(iiig,ic) < grd_capnbf(iiig,ic)) then
           call nxt_bfabs(ic,iiig,nu,rndstate,kind,idx,wlcmf,wtb)
        else
           call nxt_kpacket(ic,rndstate,kind,idx,wlcmf,wtb)
        endif
!-- (nxt_fbias: the weight of a biased emission)
        e = e*wtb
        e0 = e0*wtb
        ptcl%em = nxt_emcode(kind,idx)
        ptcl%icem = 0
        if(nxt_emloc(ptcl%em)) ptcl%icem = ic
!-- (a line photon scattered by an electron in its resonance zone: new
!-- comoving wavelength and direction)
        if(kind==em_line) call nx_linees(kind,idx,ic,rndstate,mucmf,wlcmf)
!-- (diagnostics: the line emitted)
        if(kind==em_line .or. kind==em_lines) then
           ptcl2%kem = idx
           ptcl2%cem = ic
           ptcl2%eem = e
        endif
        call nx_emitwl(kind,idx,wlcmf,ic,rndstate,wlc,iiig,lout)
        lddmc = .false.
        if(.not.lout) lddmc = (grd_sig(ic)+grd_cap(iiig,ic)+grd_capsc(iiig,ic))*dist >= trn_tauddmc
        if(lddmc) then
!-- DDMC in the new group
           ig = iiig
        else
!-- IMC (or escape outside the group grid): position uniform in the cell,
!-- direction isotropic, a line photon's along large |Q(mu)|
           ptcl2%itype = 1
           call rnd_r(r1,rndstate)
           x = (r1*grd_xarr(ix+1)**3 + (1d0-r1)*grd_xarr(ix)**3)**(1d0/3d0)
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           vx = (grd_vxarr(ix)*(grd_xarr(ix+1)-x) + grd_vxarr(ix+1)*(x-grd_xarr(ix)))/dist
           if(kind==em_line) then
              call line_mu(nxt_kq(idx,ic),(grd_vxarr(ix+1)-grd_vxarr(ix))/dist,vx/x,rndstate,mucmf)
           elseif(kind/=em_lines) then
              call rnd_r(r1,rndstate)
              mucmf = 1d0-2d0*r1
           endif
           mu = (mucmf+vx*cinv)/(1d0+vx*mucmf*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-vx*mu*cinv)
           totevelo = totevelo+e*(1d0-help)
           e = e*help
           e0 = e0*help
           ig = iiig
           ptcl%ilem = nx_ilem(kind,idx,wlc)
           call nx_place(kind,idx,wlc,ig,lout,wl,e,e0,mu,vx,x,ix,ic,ptcl2,totevelo)
           if(ptcl2%stat=='flux') return
        endif
     else
!-- re-emitted with the cell's emission spectrum
        ptcl%icem = ic

        if(glump==grp_ng) then
!       stop 'diffusion11: effective scattering with glump==ng'
           ierr = 102
           return
        endif

        if(glump==0) then
!-- sample group
           call rnd_r(r1,rndstate)
           iiig = emitgroup(r1,ic)
        else
!
           call rnd_r(r1,rndstate)
           denom2 = 1d0/(1d0-emitlump)
           denom3 = 0d0
           do iig=grp_ng,glump+1,-1
              iiig = grd_glump(iig,ic)
              if(in_nlte) then !NLTE (H only)
                help = grd_planck(iiig,ic)*grd_capemit(iiig,ic)&
                       *capemitgreyinv
              else
                help = grd_spk(iiig,ic)*grd_cap(iiig,ic)*capgreyinv
              endif
              denom3 = denom3 + help*denom2
              if(denom3>r1) exit
           enddo
        endif
        ig = iiig
        if((grd_sig(ic)+grd_cap(ig,ic)+grd_capsc(ig,ic))*dist < trn_tauddmc) then
           ptcl2%itype = 1
!-- sample wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))
!-- direction sampled isotropically
           call rnd_r(r1,rndstate)
           mu = 1d0-2d0*r1
!-- position sampled uniformly
           call rnd_r(r1,rndstate)
           x = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
!-- determine velocity coordinates (currently uses linear interpolation)
!-- x
           vhelp1 = grd_vxarr(ix)
           vhelp2 = grd_vxarr(ix+1)
           help1 = grd_xarr(ix)
           help2 = grd_xarr(ix+1)
           vx = (vhelp1*(help2-x)+vhelp2*(x-help1))/dx(ix)
!
!-- doppler and aberration corrections
           mu = (mu+vx*cinv)/(1.0+vx*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-vx*mu*cinv)
           totevelo = totevelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1.0-vx*mu*cinv)
        endif

     endif
  endif

!-- eliminate particles crossing inner boundary
  if (x  < grd_xarr(1)) then
    ptcl2%stat = 'dead'
    return
  endif

!-- Doppler shift energy weights and wavelengths ! added - gaw
  totevelo = totevelo + e*(1d0-exp(-grd_divv(ix)*ddmct/3))
  e = e*exp(-grd_divv(ix)*ddmct/3)
  e0 = e0*exp(-grd_divv(ix)*ddmct/3)

end subroutine diffusion11
! vim: fdm=marker
