!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine nltex_gas_update(it)
!-----------------------------------------------------------------------
! NLTE elements in this rank's gas cells: statistical equilibrium with the
! previous iteration's Monte Carlo estimators, then the electron density
! (NLTE elements from their solution, the others from the EOS).
!-----------------------------------------------------------------------
  use gasmod
  use gridmod
  use nlteatommod, only:nxa_nel,nxa_el,nxa_nlin
  use nltextrmod
  use mpimod, only:lmpi0,impi
  use inputparmod, only:in_nltex_ali,in_nltex_ldb,in_nt_lumg,in_srctype
  use sourcemod, only:src_tsh
  use nltentmod, only:nxn_on,nxn_sloc
  implicit none
  integer,intent(in) :: it
  integer :: i,l,iel
  character(32) :: fname
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
!-- and the share in_nt_lumg of the injected energy (radioactive decay)
     if(nxn_on .and. in_nt_lumg>0d0 .and. in_srctype=='lumg' .and. gas_mass(i)>0d0) then
        if(src_tsh(l)<=0d0) edep(i) = edep(i) + in_nt_lumg*grd_emit(l)/gas_volcl(i)
     endif
  enddo
  if(in_nltex_ali .or. in_nltex_ldb>0d0) then
     call nxt_gas_solve(it,gas_ncell,gas_temp,gas_radtemp,nelv,neoth,qb,qv,gas_icell1,edep, &
          kap=gas_cap,fcl=gas_fcl,dr=dr,vel=vel)
  else
     call nxt_gas_solve(it,gas_ncell,gas_temp,gas_radtemp,nelv,neoth,qb,qv,gas_icell1,edep,vel=vel)
  endif
  if(lmpi0 .and. it>1 .and. nxt_have_mc) write(6,'(1x,a,i3,a,es9.2,a,es9.2,a)') 'NLTE (it',it, &
       '): largest change of a state above 1e-4 per cell, median',nxt_dxstat(1),', max',nxt_dxstat(2),' (rank 0)'
  if(lmpi0 .and. nxt_nactl>0) write(6,'(1x,a,i8,a,i8,a,es8.1,a)') 'NLTE lines: transport steps at',nxt_nactl, &
       ' of',nxa_nlin,' lines (Sobolev depth above',nxt_taumin,' somewhere)'
  if(lmpi0 .and. it>1 .and. nxt_have_mc .and. nxt_dmpsec) write(6,'(1x,a,f6.3,a,f6.3)') &
       'NLTE damping (in_nltex_dmpsec): weight of the new solution, median',nxt_dmpstat(1),', smallest',nxt_dmpstat(2)
!-- (the history of the adaptive damping, for a continued run)
  if(nxt_dmpsec) call nxt_write_dmp(impi)
!-- (fast electrons: the shares of their energy in each cell)
  if(nxn_on) then
     write(fname,'(a,i4.4)') 'output.nltex_nt_',impi
     open(42,file=fname,status='replace')
     write(42,'(a,i3)') '# fast electrons (Spencer-Fano): shares of their energy, iteration',it
     write(42,'(a)') '#  cell   edep[erg/s/cm^3]   heat      ionization  excitation'
     do i=1,gas_ncell
        write(42,'(i7,es14.4,3f12.5)') gas_icell1+i-1,edep(i),nxn_sloc(:,i)
     enddo
     close(42)
  endif
!-- electrons per atom
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
     gas_nelec(i) = nxt_neloc(i)*gas_volcl(i)/gas_natom(i)
  enddo
end subroutine nltex_gas_update
! vim: fdm=marker
