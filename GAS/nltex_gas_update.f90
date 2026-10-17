!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine nltex_gas_update(it)
!-----------------------------------------------------------------------
! NLTE elements in this rank's gas cells: statistical equilibrium with the
! previous iteration's Monte Carlo estimators, then the electron density
! (NLTE elements from their solution, the others from the EOS).
!-----------------------------------------------------------------------
  use gasmod
  use gridmod
  use nlteatommod, only:nxa_nel,nxa_el
  use nltextrmod
  use mpimod, only:lmpi0
  use inputparmod, only:in_nltex_ali,in_nltex_ldb
  implicit none
  integer,intent(in) :: it
  integer :: i,l,iel
  real*8 :: nelv(nxa_nel,gas_ncell),neoth(gas_ncell),qb(gas_ncell),qv(gas_ncell)
  real*8 :: edep(gas_ncell),dr(gas_ncell),vel(gas_ncell)
!
  do i=1,gas_ncell
     l = gas_icell1 + i - 1
     if(gas_mass(i)<=0d0) then
        nelv(:,i) = 0d0
        neoth(i) = 0d0
     else
        do iel=1,nxa_nel
           nelv(iel,i) = gas_natom(i)*gas_natomfr(nxa_el(iel)%z,i)/gas_volcl(i)
        enddo
        neoth(i) = nxt_neothloc(i)*gas_natom(i)/gas_volcl(i)
     endif
     dr(i) = grd_xarr(l+1)-grd_xarr(l)
     vel(i) = 0.5d0*(grd_vxarr(l)+grd_vxarr(l+1))
     qb(i) = (grd_vxarr(l+1)-grd_vxarr(l))/(grd_xarr(l+1)-grd_xarr(l))
     qv(i) = (grd_vxarr(l+1)+grd_vxarr(l))/(grd_xarr(l+1)+grd_xarr(l))
!-- microclumping: the Sobolev depths of the clump populations with the
!-- volume-mean opacity, K f/Q
     qb(i) = qb(i)/gas_fcl(i)
     qv(i) = qv(i)/gas_fcl(i)
!-- power absorbed from X-rays (fast photoelectrons), with the opacity of the
!-- previous iteration
     edep(i) = 0d0
     if(it>1 .and. nxt_have_mc .and. gas_mass(i)>0d0) &
          edep(i) = nxt_ntdep(nxt_ng,nxt_nug,nxt_dnug,nxt_jg(:,l),gas_cap(:,i))
  enddo
  if(in_nltex_ali .or. in_nltex_ldb>0d0) then
     call nxt_gas_solve(it,gas_ncell,gas_temp,gas_radtemp,nelv,neoth,qb,qv,gas_icell1,edep, &
          kap=gas_cap,fcl=gas_fcl,dr=dr,vel=vel)
  else
     call nxt_gas_solve(it,gas_ncell,gas_temp,gas_radtemp,nelv,neoth,qb,qv,gas_icell1,edep,vel=vel)
  endif
  if(lmpi0 .and. it>1 .and. nxt_have_mc) write(6,'(1x,a,i3,a,es9.2,a,es9.2,a)') 'NLTE (it',it, &
       '): largest change of a state above 1e-4 per cell, median',nxt_dxstat(1),', max',nxt_dxstat(2),' (rank 0)'
!-- electrons per atom
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
     gas_nelec(i) = nxt_neloc(i)*gas_volcl(i)/gas_natom(i)
  enddo
end subroutine nltex_gas_update
! vim: fdm=marker
