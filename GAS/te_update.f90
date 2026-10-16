!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine te_update(it)
!-----------------------------------------------------------------------
! Electron temperature of this rank's gas cells from the thermal balance
! with the radiation field of the previous Monte Carlo iteration
! (tebalmod).  Non-radiative heating: energy injected outside shocks by the
! 'lumg' source (shock cells radiate theirs with the shock spectrum).
! At the last iteration each rank writes the heating and cooling terms of
! its cells at the new temperature to output.tebal_<rank>.
!-----------------------------------------------------------------------
  use gasmod
  use gridmod
  use sourcemod, only:src_tsh
  use inputparmod
  use nlteatommod, only:nxa_nel,nxa_el,nxa_nst,nxa_nlin
  use nltexmod, only:nx_solve_cell,nx_jline
  use nltextrmod
  use tebalmod
  use xraysrcmod, only:xs_on,xs_px
  use mpimod, only:lmpi0,impi
  implicit none
  integer,intent(in) :: it
  integer :: i,l,ierr,nfail,iel,k
  real*8 :: neoth,qb,qv,gext,tnew,dmax,dmin,edep,ne,tsol
  real*8 :: nel(nxa_nel),x(nxa_nst),nti(nxa_nst),ntx(nxa_nst),beta(nxa_nlin),terms(te_nterm)
  real*8 :: jl(nxa_nlin)
  real*4 :: capabs(nxt_ng)
!-- per cell: previous proposed change of ln T and damping weight
  real*8,allocatable,save :: dprev(:),dloc(:)
  logical :: lout
  character(24) :: fname
!
  if(.not.nxt_have_mc) return
  if(.not.allocated(dprev)) then
     allocate(dprev(gas_ncell),dloc(gas_ncell))
     dprev = 0d0
     dloc = te_damp
  endif
  lout = it==in_niter
  if(lout) then
     write(fname,'(a,i4.4)') 'output.tebal_',impi
     open(41,file=fname,status='replace')
     write(41,'(a)') '# thermal balance per cell at the new T_e [erg/s/cm^3]: '// &
          'LTE-opacity net cooling, NLTE photoionization heating, NLTE recombination cooling, '// &
          'NLTE collisional net cooling, Compton net cooling, external heating'
     write(41,'(a)') '#  cell        T_e       T_old     lte_net     pi_heat     rec_cool    coll_net    compton_net ext_heat'
  endif
  nfail = 0
  dmax = 1d0
  dmin = 1d0
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
     l = gas_icell1 + i - 1
     do iel=1,nxa_nel
        nel(iel) = gas_natom(i)*gas_natomfr(nxa_el(iel)%z,i)/gas_volcl(i)
     enddo
     neoth = nxt_neothloc(i)*gas_natom(i)/gas_volcl(i)
     qb = (grd_vxarr(l+1)-grd_vxarr(l))/(grd_xarr(l+1)-grd_xarr(l))/gas_fcl(i)
     qv = (grd_vxarr(l+1)+grd_vxarr(l))/(grd_xarr(l+1)+grd_xarr(l))/gas_fcl(i)
!-- (per unit clump volume: microclumping)
     gext = 0d0
     if(in_srctype=='lumg' .and. src_tsh(l)<=0d0) gext = grd_emit(l)/gas_volcl(i)
!-- (time-dependent NLTE: the share released by the recombining gas is
!-- emitted by the macro-atom, not by the thermal pool)
     if(nxt_dtinv>0d0) then
        if(allocated(nxt_prel)) gext = gext*(1d0 - nxt_prel(l))
     endif
!-- (and the X-rays of wind shocks are mechanical, not thermal)
     if(xs_on) then
        if(allocated(xs_px)) gext = gext - grd_emit(l)*xs_px(l)/gas_volcl(i)
     endif
     gext = max(gext,0d0)
!-- fast photoelectrons: non-thermal ionization and excitation of H and He
     edep = nxt_ntdep(nxt_ng,nxt_nug,nxt_dnug,nxt_jg(:,l),gas_cap(:,i))
     call nxt_ntrates(nel,nxt_xloc(:,i),edep,nti,ntx)
!-- absorption by the LTE elements (NLTE bound-free is in the NLTE terms)
     capabs = gas_cap(:,i) - gas_capnbf(:,i)
!-- (accelerated lambda iteration of the bound-free transitions)
     if(in_nltex_ali) call nxt_ali_cell(nxt_ng,gas_cap(:,i),gas_fcl(i), &
          grd_xarr(l+1)-grd_xarr(l))
     if(nxt_dtinv>0d0) then
        call te_balance(gas_temp(i),gas_radtemp(i),nel,neoth,qb,qv,nxt_ng,nxt_nug,nxt_dnug, &
             nxt_jg(:,l),nxt_jb(:,l),nxt_ljb(:,l),nti,ntx,capabs,gas_capemit(:,i),gext, &
             nxt_xloc(:,i),tnew,ierr,tsol,xprev=nxt_xprev(:,l),dtinv=nxt_dtinv)
     else
        call te_balance(gas_temp(i),gas_radtemp(i),nel,neoth,qb,qv,nxt_ng,nxt_nug,nxt_dnug, &
             nxt_jg(:,l),nxt_jb(:,l),nxt_ljb(:,l),nti,ntx,capabs,gas_capemit(:,i),gext, &
             nxt_xloc(:,i),tnew,ierr,tsol)
     endif
     call te_step(gas_temp(i),tsol,dprev(i),dloc(i),tnew)
     if(ierr/=0) nfail = nfail + 1
     if(lout) then
        x = nxt_xloc(:,i)
        if(nxt_dtinv>0d0) then
           call nx_solve_cell(tnew,gas_radtemp(i),nel,neoth,qb,qv,nxt_ng,nxt_nug,nxt_dnug, &
                nxt_jg(:,l),nxt_jb(:,l),nxt_ljb(:,l),nti,ntx,x,ne,beta,ierr, &
                xprev=nxt_xprev(:,l),dtinv=nxt_dtinv)
        else
           call nx_solve_cell(tnew,gas_radtemp(i),nel,neoth,qb,qv,nxt_ng,nxt_nug,nxt_dnug, &
                nxt_jg(:,l),nxt_jb(:,l),nxt_ljb(:,l),nti,ntx,x,ne,beta,ierr)
        endif
        do k=1,nxa_nlin
           jl(k) = nx_jline(k,nxt_ng,nxt_nug,nxt_jg(:,l),nxt_jb(k,l),nxt_ljb(k,l))
        enddo
        call te_terms(tnew,nxt_ng,nxt_nug,nxt_dnug,nxt_jg(:,l),capabs,gas_capemit(:,i),nel,x,ne, &
             gext,terms,beta=beta,jl=jl)
        write(41,'(i7,2f12.1,6es12.4)') l,tnew,gas_temp(i),terms
     endif
     call nxt_ali_off
     dmax = max(dmax,tnew/gas_temp(i))
     dmin = min(dmin,tnew/gas_temp(i))
     gas_temp(i) = tnew
  enddo
  if(lout) close(41)
  if(lmpi0) write(6,'(1x,a,i3,a,f6.3,a,f6.3,a,i4)') 'T_e (it',it,'): change factor', &
       dmin,' ..',dmax,' (rank 0), cells without root:',nfail
end subroutine te_update
! vim: fdm=marker
!
!
subroutine nebular_wtr
!-----------------------------------------------------------------------
! Radiation temperature and dilution factor of the Monte Carlo radiation
! field of the previous iteration in this rank's cells (Mazzali & Lucy
! 1993): T_R from the mean frequency, nu_mean = 3.832 kT_R/h, and
! W = pi J/(sigma T_R^4).  The LTE elements are excited at T_R
! (gas_tion) and ionized in the nebular approximation.  Without estimators
! (first iteration) T_R is the input temperature and W = 1.
!-----------------------------------------------------------------------
  use gasmod
  use nltextrmod
  implicit none
  real*8,parameter :: hh = 6.6260693d-27, kb = 1.3806505d-16, sb = 5.670400d-5, pi = 3.141592653589793d0
  integer :: i,l
  real*8 :: j0,j1
  if(.not.nxt_have_mc) return
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
     l = gas_icell1 + i - 1
     j0 = sum(nxt_jg(:,l)*nxt_dnug)
     j1 = sum(nxt_nug*nxt_jg(:,l)*nxt_dnug)
     if(j0<=0d0) cycle
     gas_tion(i) = hh*(j1/j0)/(3.832d0*kb)
     gas_wdil(i) = min(1d0, pi*j0/(sb*gas_tion(i)**4))
  enddo
end subroutine nebular_wtr
! vim: fdm=marker
