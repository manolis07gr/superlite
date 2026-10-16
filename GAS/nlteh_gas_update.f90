!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine nlteh_gas_update(it)
!-----------------------------------------------------------------------
! NLTE hydrogen in this rank's gas cells: statistical equilibrium with the
! previous iteration's Monte Carlo estimators, then the electron density
! (H from the NLTE solution, other elements from the EOS).
!-----------------------------------------------------------------------
  use gasmod
  use gridmod
  use nltehtrmod
  implicit none
  integer,intent(in) :: it
  integer :: i,l
  real*8 :: nh(gas_ncell),neoth(gas_ncell),qb(gas_ncell),qv(gas_ncell)
!
  do i=1,gas_ncell
     l = gas_icell1 + i - 1
     if(gas_mass(i)<=0d0) then
        nh(i) = 0d0
        neoth(i) = 0d0
     else
        nh(i) = gas_natom(i)*gas_natomfr(1,i)/gas_vol(i)
        neoth(i) = nht_neothloc(i)*gas_natom(i)/gas_vol(i)
     endif
     qb(i) = (grd_vxarr(l+1)-grd_vxarr(l))/(grd_xarr(l+1)-grd_xarr(l))
     qv(i) = (grd_vxarr(l+1)+grd_vxarr(l))/(grd_xarr(l+1)+grd_xarr(l))
  enddo
  call nlteh_gas_solve(it,gas_ncell,gas_temp,gas_radtemp,nh,neoth,qb,qv,gas_icell1)
!-- electrons per atom
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
     gas_nelec(i) = nht_neloc(i)*gas_vol(i)/gas_natom(i)
  enddo
end subroutine nlteh_gas_update
! vim: fdm=marker
