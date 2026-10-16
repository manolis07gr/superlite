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
     capgreyinv = max(1d0/grd_capgrey(ic),0d0) !catch nans
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
     if(glump==grp_ng) then
        grd_vlump(:,ic) = [1d0, grd_capgrey(ic), grd_capemitgrey(ic), 0d0]
     else
!-- (with the NLTE macro-atoms: Kirchhoff emission with all absorption
!-- opacities; the macro-atom is not followed in DDMC)
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
