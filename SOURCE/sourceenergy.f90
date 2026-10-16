!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
subroutine sourceenergy

  use gridmod
  use sourcemod
  use totalsmod
  use physconstmod
  use inputparmod
  use groupmod
  use nltextrmod, only:nxt_on,nxt_dtinv,nxt_release,nxt_erel,nxt_prel
  use xraysrcmod, only:xs_on,xs_pow,xs_px
  implicit none
!##################################################
!This subroutine computes the distribution of source particles
!A fraction of the source particle number src_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################
  integer :: i
  real*8 :: help
  ! real*8,dimension(grd_ncell) :: P_adia
  ! real*8,dimension(grp_ng) :: wlm,dwl

! Calculating fictitious emission energy per cell
!-- thermal source
  grd_emit = 0.d0
  if(in_srctype=='lumg') then
!-- energy injected where the luminosity of the input model grows outward
!-- (shock dissipation, cooling, diffusion from the interior); grd_lum(i) is
!-- L at the inner edge of cell i and in_L_bol the emergent luminosity.  The
!-- luminosity entering at the inner boundary is a boundary source.
     do i=1,grd_ncell-1
        grd_emit(i) = max(0d0, grd_lum(i+1) - grd_lum(i))
     enddo
     grd_emit(grd_ncell) = max(0d0, in_L_bol - grd_lum(grd_ncell))
!-- a difference of in_L_bol and the input's outer L below 1e-4 is rounding
!-- (it would heat the outermost, often near-empty cell)
     if(grd_emit(grd_ncell) < 1d-4*in_L_bol) grd_emit(grd_ncell) = 0d0
     tot_esurf = max(0d0, grd_lum(1))
     help = sum(grd_emit) + tot_esurf
     if(help<=0d0) stop 'sourceenergy: lumg source needs a lum column'
     grd_emit = grd_emit*in_L_bol/help
     tot_esurf = tot_esurf*in_L_bol/help
!-- time-dependent NLTE statistical equilibrium: the power released by the
!-- recombining gas (nxt_release) is emitted by macro-atom activation in the
!-- states that empty.  It is part of the injected energy (the hydro L(r)
!-- contains the energy the gas releases), added where it is larger; stored
!-- energy (ionizing gas) reduces the emission.  Shock cells keep theirs.
     if(nxt_on .and. nxt_dtinv>0d0) then
        call nxt_release(grd_ncell,grd_vol*grd_fcl)
        do i=1,grd_ncell
           if(src_tsh(i)>0d0) then
              nxt_prel(i) = 0d0
           elseif(nxt_erel(i)>=0d0) then
              grd_emit(i) = max(grd_emit(i),nxt_erel(i))
              if(grd_emit(i)>0d0) nxt_prel(i) = nxt_erel(i)/grd_emit(i)
           else
              grd_emit(i) = max(0d0,grd_emit(i) + nxt_erel(i))
              nxt_prel(i) = 0d0
           endif
        enddo
     endif
!-- X-ray emission of wind shocks (xraysrcmod): mechanical energy beyond
!-- the injected energy
     if(xs_on) then
        do i=1,grd_ncell
           help = grd_emit(i)
           grd_emit(i) = grd_emit(i) + xs_pow(i)
           xs_px(i) = 0d0
           if(grd_emit(i)>0d0) then
              xs_px(i) = xs_pow(i)/grd_emit(i)
              if(nxt_on .and. nxt_dtinv>0d0) nxt_prel(i) = nxt_prel(i)*help/grd_emit(i)
           endif
        enddo
     endif
     tot_sthermal = sum(grd_emit) + tot_esurf
     return
  endif
!-- surface source
  grd_emit = pc_pi*pc_c*grd_rcell**2*pc_acoef*grd_radtemp**4
!-- volume source (experimental)
  ! grd_emit = grd_vol*pc_acoef*grd_radtemp**4
!-- luminosity source (experimental)
  ! grd_emit = grd_lum !CAUTION:grd_emit can be >=0, will result in sampling problems
!-- emissivity source (experimental)
  ! do i=1,grp_ng
  !   dwl(i) = grp_wl(i+1) - grp_wl(i)
  !   wlm(i) = 0.5d0*(grp_wl(i+1) + grp_wl(i))
  ! enddo !i
  ! do i=1,grd_ncell-1
  !  grd_emit(i) = 4*pc_pi*grd_vol(i)*sum(grd_emiss(:,i)*dwl*pc_c/wlm**2)
!--  Fick's law (experimental)
  ! do i=1,grd_ncell-1
  !   grd_emit(i) = (-1.d0*pc_pi16*grd_rcell(i)**2/(3.d0*grd_capross(i)))* &
  !                 pc_acoef*pc_c*grd_radtemp(i)**3* &
  !                 (grd_radtemp(i+1)-grd_radtemp(i)) / &
  !                 (grd_rcell(i+1)-grd_rcell(i))
  ! enddo
!-- plus adiabatic compression/expansion source (experimental)
  ! P_adia = grd_natom * pc_kb / grd_tempinv
  ! grd_emit = grd_emit - P_adia*grd_divv

!-- renormalize source energy to the user input bolometric luminosity
  grd_emit = grd_emit / sum(grd_emit) * in_L_bol
  tot_sthermal = sum(grd_emit)

!-- inner boundary surface luminosity source (experimental)
  ! if(in_L_inner>0d0) tot_esurf = in_L_inner

end subroutine sourceenergy
! vim: fdm=marker
