!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine ddmc_lumps
!-----------------------------------------------------------------------
! DDMC group lumping of every cell, once per iteration for diffusion11
! (before: at each step that entered a new cell).  The groups with
! absorption depth >= trn_taulump and total depth >= trn_tauddmc over the
! cell form the lump, grd_glump(1:grd_nlump) in ascending order; the other
! groups follow from the end in descending order.  grd_vlump holds the
! lumped emission share, absorption and emission opacities, and Doppler
! opacity.
!-----------------------------------------------------------------------
  use gridmod
  use groupmod
  use transportmod
  use inputparmod
  use bbsobmod, only:grd_capsc
  use nltextrmod, only:nxt_on
  implicit none
  integer :: ix,ic,iig,glump,gunlump
  real*8 :: dist,speclump,capgreyinv,capemitgreyinv,emitlump
!
  do ix=1,grd_nx
     ic = grd_icell(ix,1,1)
     dist = grd_xarr(ix+1) - grd_xarr(ix)
     capgreyinv = max(1d0/grd_capgreyk(ic),0d0) !catch nans
     capemitgreyinv = max(1d0/grd_capemitgrey(ic),0d0)
     glump = 0
     gunlump = grp_ng
     grd_glump(:,ic) = 0
     speclump = grd_opaclump(7,ic)
     if(speclump/=0d0) then
        do iig=1,grp_ng
           if(grd_cap(iig,ic)*dist >= trn_taulump .and. &
                (grd_sig(ic) + grd_cap(iig,ic)+grd_capsc(iig,ic))*dist >= trn_tauddmc) then
              glump = glump + 1
              grd_glump(glump,ic) = int(iig,2)
           else
              grd_glump(gunlump,ic) = int(iig,2)
              gunlump = gunlump - 1
           endif
        enddo
     endif
     grd_nlump(ic) = glump
!-- cumulative absorption rate over the lump (NLTE: absorbing group)
     emitlump = 0d0
     do iig=1,glump
        emitlump = emitlump + grd_spk(grd_glump(iig,ic),ic)*grd_cap(grd_glump(iig,ic),ic)
        grd_lumpcum(iig,ic) = real(emitlump)
     enddo
     if(glump==grp_ng) then
        grd_vlump(:,ic) = [1d0, grd_capgrey(ic), grd_capemitgrey(ic), 0d0]
     else
!-- (with the NLTE macro-atoms: the emission share of the absorption
!-- opacities in the source spectrum; diffusion11 then re-emits through the
!-- macro-atom instead)
        if(in_nlte .and. .not.nxt_on) then
           emitlump = grd_opaclump(9,ic)*capemitgreyinv
        else
           emitlump = grd_opaclump(8,ic)*capgreyinv
        endif
        grd_vlump(:,ic) = [min(emitlump,1d0), grd_opaclump(8,ic)*speclump, &
             grd_opaclump(9,ic)*speclump, grd_opaclump(11,ic)*speclump]
     endif
  enddo
end subroutine ddmc_lumps
! vim: fdm=marker
!
!
subroutine ddmc_spectra
!-----------------------------------------------------------------------
! Spectral weights of the radiation in the DDMC cells, once per iteration
! for all cells: the Planck group integrals at T_e (grd_planck) and, with
! the NLTE macro-atoms, the source function S = (kappa_LTE B + eta)/kappa
! of the absorption opacity kappa (grd_spk, normalized; eta: the
! recombination emissivity of the NLTE levels, grd_emnbf), with the mean
! absorption opacity over it (grd_capgreyk).  DDMC samples leakage and
! re-emission with these weights: in LTE S = B and nothing changes; in a
! thick NLTE continuum S = B/b (the Lyman continuum of a wind, b1 << 1),
! where the Planck weights at T_e would thermalize the radiation.
!-----------------------------------------------------------------------
  use gridmod
  use groupmod
  use nltextrmod, only:nxt_on,grd_emnbf
  implicit none
  integer :: l,ig
  real*8 :: s(grp_ng),tot
!
  do l=1,grd_ncell
     call specintv(grd_tempinv(l),grp_ng,grd_planck(:,l))
     if(.not.nxt_on) then
        grd_spk(:,l) = grd_planck(:,l)
        grd_capgreyk(l) = grd_capgrey(l)
        cycle
     endif
     do ig=1,grp_ng
        if(grd_cap(ig,l)>0.) then
           s(ig) = (grd_planck(ig,l)*grd_capemit(ig,l) + grd_emnbf(ig,l))/grd_cap(ig,l)
        else
           s(ig) = grd_planck(ig,l)
        endif
     enddo
     tot = sum(s)
     if(tot>0d0) then
        grd_spk(:,l) = s/tot
     else
        grd_spk(:,l) = grd_planck(:,l)
     endif
     grd_capgreyk(l) = sum(grd_spk(:,l)*grd_cap(:,l))
  enddo
end subroutine ddmc_spectra
! vim: fdm=marker
