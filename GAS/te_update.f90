!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine te_update(it)
!-----------------------------------------------------------------------
! Electron temperature of this rank's gas cells from the thermal balance
! with the radiation field of the previous Monte Carlo iteration
! (tebalmod).  Non-radiative heating: energy injected outside shocks by the
! 'lumg' source (shock cells radiate theirs with the shock spectrum).
!-----------------------------------------------------------------------
  use gasmod
  use gridmod
  use sourcemod, only:src_tsh
  use inputparmod
  use nltehmod, only:nh_nlev
  use nltehtrmod
  use tebalmod
  use mpimod, only:lmpi0
  implicit none
  integer,intent(in) :: it
  integer :: i,l,ierr,nfail
  real*8 :: nh,neoth,qb,qv,gext,tnew,dmax,dmin,edep,nonth,nonexc
!
  if(.not.nht_have_mc) return
  nfail = 0
  dmax = 1d0
  dmin = 1d0
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
     l = gas_icell1 + i - 1
     nh = nht_nhloc(i)
     neoth = nht_neothloc(i)*gas_natom(i)/gas_vol(i)
     qb = (grd_vxarr(l+1)-grd_vxarr(l))/(grd_xarr(l+1)-grd_xarr(l))
     qv = (grd_vxarr(l+1)+grd_vxarr(l))/(grd_xarr(l+1)+grd_xarr(l))
     gext = 0d0
     if(in_srctype=='lumg' .and. src_tsh(l)<=0d0) gext = grd_emit(l)/grd_vol(l)
!-- fast photoelectrons: non-thermal ionization and excitation of H
     edep = nlteh_ntdep(nht_ng,nht_nug,nht_dnug,nht_jg(:,l),gas_cap(:,i))
     call nlteh_ntrates(nh,nht_xloc(:,i),edep,nonth,nonexc)
     call te_balance(gas_temp(i),gas_radtemp(i),nh,neoth,qb,qv,nht_ng,nht_nug,nht_dnug, &
          nht_jg(:,l),nht_jb(:,l),nht_ljb(:,l),nonth,nonexc,gas_capemit(:,i),gext,nht_xloc(:,i), &
          tnew,ierr)
     if(ierr/=0) nfail = nfail + 1
     dmax = max(dmax,tnew/gas_temp(i))
     dmin = min(dmin,tnew/gas_temp(i))
     gas_temp(i) = tnew
  enddo
  if(lmpi0) write(6,'(1x,a,i3,a,f6.3,a,f6.3,a,i4)') 'T_e (it',it,'): change factor', &
       dmin,' ..',dmax,' (rank 0), cells without root:',nfail
end subroutine te_update
! vim: fdm=marker
