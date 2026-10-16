!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nltextrmod
!-----------------------------------------------------------------------
! NLTE elements in the radiation transport: cell state, Monte Carlo
! estimators and Lucy's (2002, A&A 384, 725; 2003, A&A 403, 261) macro-atom
! for all NLTE elements (model atoms: nlteatommod; statistical
! equilibrium: nltexmod).
!
! The lines of the NLTE ions are not part of the group opacities.  A packet
! meets a line where its comoving wavelength equals the line wavelength
! (exact crossing, so a non-monotonic flow can have several resonance points
! per line) and interacts there with probability 1 - exp(-tau),
! tau = K/|Q(mu)| (Sobolev, Q(mu) = mu^2 dv/dr + (1-mu^2) v/r).  Every
! crossing is tallied in the resonance estimator Jb.
!
! An absorbed line packet activates the macro-atom of the line's element in
! the upper level; bound-free absorption at nu by level l activates the
! continuum state of l with probability nu_l/nu and the thermal pool
! otherwise.  The macro-atom jumps between the states of the element (all
! levels of its NLTE ions and the top ion) with probabilities set by the
! rates, weighted by the state energies (Lucy), and de-activates by emitting
! a line, recombination-continuum or two-photon packet, or by collisional
! de-excitation and three-body recombination into the thermal pool
! (k-packet).  Stimulated emission and stimulated recombination are negative
! absorption: they reduce the upward rates.  A k-packet is emitted thermally
! (emissivity of the other opacities), as recombination cooling, or
! collisionally excites or ionizes an NLTE element, in proportion to the
! cooling rates.  Packets keep their comoving energy, so the scheme
! conserves energy exactly.
!
! The random walk is an absorbing Markov chain.  Per element and cell the
! fundamental matrix N = (I-Q)^-1 of its states gives, for each start state,
! the probabilities to de-activate from each state (or to enter the thermal
! pool); the emission channel is drawn from the de-activating state's
! channels.  The thermal pool's outcome distribution includes all returns
! from the elements (it is proportional to the non-returning flows), so an
! activation takes at most three draws, with the statistics of the full
! walk and without the long walks of collision-dominated gas.
!-----------------------------------------------------------------------
  use nlteatommod
  use nltexmod
  use nltehmod, only:bplanck,nlteh_twophoton_y,nlteh_svs
  use randommod
  implicit none
  private
!
  real*8,parameter :: pi = 3.141592653589793d0, cc = 2.99792458d10
  real*8,parameter :: hh = 6.6260693d-27, kb = 1.3806505d-16, arad = 7.5657d-15
!
  logical,public :: nxt_on = .false.
  real*8,public :: nxt_damp = 0.5d0       !weight of the new solution between iterations
  integer,public :: nxt_njbmin = 30       !crossings needed to use the resonance estimator
  logical,public :: nxt_have_mc = .false. !Monte Carlo estimators available
  real*8,public :: nxt_tau0 = 0d0         !benchmark: pure-scattering H-alpha with this Sobolev depth
  real*8,public :: nxt_ntemin = 0d0       !photons above this energy [erg] make fast photoelectrons (0: off)
  real*8,public :: nxt_taumin = 1d-4      !active lines: Sobolev depth above this somewhere
!
  integer,public :: nxt_ncell = 0, nxt_ng = 0
!-- all grid cells
  real*8,allocatable,public :: nxt_x(:,:)     !(nst,ncell) populations (fractions of each element)
  real*8,allocatable,public :: nxt_nel(:,:)   !(nel,ncell) element number densities [cm^-3]
  real*8,allocatable,public :: nxt_ne(:)      !(ncell) electron density [cm^-3]
  real*8,allocatable,public :: nxt_kq(:,:)    !(nlin,ncell) Sobolev strength K [1/s]
  real*8,allocatable,public :: nxt_beta(:,:)  !(nlin,ncell) escape probability
  real*8,allocatable,public :: nxt_jg(:,:)    !(ng,ncell) MC mean intensity [erg/s/cm^2/Hz/sr]
  real*8,allocatable,public :: nxt_jb(:,:)    !(nlin,ncell) resonance mean intensity
  logical,allocatable,public :: nxt_ljb(:,:)  !(nlin,ncell) jb from enough crossings
  real*8,allocatable,public :: nxt_jbsum(:,:) !(nlin,ncell) tally: sum e c/(nu |Q|)
  real*8,allocatable,public :: nxt_jbcnt(:,:) !(nlin,ncell) tally: crossings
  real*4,allocatable,public :: grd_capnbf(:,:) !(ng,ncell) NLTE bound-free part of grd_cap
!-- this rank's gas cells
  real*8,allocatable,public :: nxt_xloc(:,:)  !(nst,gas_ncell)
  real*8,allocatable,public :: nxt_neloc(:)   !(gas_ncell)
  real*8,allocatable,public :: nxt_nelloc(:,:) !(nel,gas_ncell)
  real*8,allocatable,public :: nxt_neothloc(:) !(gas_ncell) electrons per atom from the other elements
  real*4,allocatable,public :: gas_capnbf(:,:) !(ng,gas_ncell)
!-- lines and groups
  integer,allocatable,public :: nxt_lgrp(:)   !(nlin) group of the line
  integer,allocatable,public :: nxt_isrt(:)   !(nlin) position of line k in the active list, 0 if inactive
  integer,public :: nxt_nact = 0
  integer,allocatable,public :: nxt_act(:)    !(nlin) active lines sorted by wavelength
  real*8,allocatable,public :: nxt_wlact(:)   !(nlin) their wavelengths
  real*8,allocatable,public :: nxt_nug(:),nxt_dnug(:) !(ng) group centre frequency and width
  integer,public :: nxt_kha = 0               !H-alpha (benchmark)
!-- macro-atom tables
  integer,allocatable,public :: nxt_ecoff(:)  !(nel) offset of the element's exit tables
  integer,public :: nxt_nec = 0
  real*8,allocatable,public :: nxt_ecum(:,:)  !(nec,ncell) per start state: cumulative over exit state, then k
  integer,public :: nxt_nch = 0               !emission channels
  integer,allocatable,public :: nxt_chst(:)   !(nst+1) first channel of each (global) state
  integer,allocatable,public :: nxt_chk(:),nxt_chi(:) !(nch) kind (em_*), index (line, level, two-photon)
  real*8,allocatable,public :: nxt_chcum(:,:) !(nch,ncell) cumulative within each state
  integer,public :: nxt_nko = 0               !k-packet outcomes: thermal, bf cooling per level, exit per state
  real*8,allocatable,public :: nxt_kcum(:,:)  !(nko,ncell)
  integer,public :: nxt_n2p = 0               !two-photon transitions
  integer,allocatable,public :: nxt_2pion(:),nxt_2pk(:)
  real*8,allocatable,public :: nxt_2pnu(:)
  real*8,allocatable,public :: nxt_cmax(:)    !(table points) running max of sigma u^3 (bf emission sampling)
  integer,allocatable,public :: nxt_cmoff(:)  !(nion) offset of each ion in nxt_cmax
!-- non-thermal: target states (global), 0 if absent
  integer,public :: nxt_sh1 = 0, nxt_sh2 = 0, nxt_shc = 0      !H I ground, its excitation target, continuum
  integer,public :: nxt_she1 = 0, nxt_she2 = 0, nxt_shec = 0   !He I ground, excitation target, continuum
  integer,public :: nxt_elh = 0, nxt_elhe = 0
  logical,public :: nxt_isel(0:99) = .false. !element (Z) in NLTE
  integer,public :: nxt_elz(0:99) = 0         !NLTE element index of Z, 0 if none
!-- NLTE ion of (Z, stage); -1: the top ion of an NLTE element; 0: none
  integer,public :: nxt_ionzs(0:99,30) = 0
  real*8,allocatable,public :: nxt_etop(:)    !(nion) highest level of the model atom [cm^-1]
!-- emission kinds returned by the macro-atom
  integer,parameter,public :: em_line = 1, em_bf = 2, em_thermal = 3, em_2ph = 4
!
  public :: nxt_init,nxt_gas_solve,nxt_store_jg,nxt_tally_reset,nxt_tally_finish
  public :: nxt_tables,nxt_neighbors,nxt_activate,nxt_activate_line,nxt_kpacket
  public :: nxt_bfabs,nxt_nonthermal,nxt_ntdep,nxt_ntrates,nxt_write,nxt_bffreq
  public :: nxt_shares,nxt_dealloc
!
contains
!
!
subroutine nxt_init(ncell,gncell,ng,wlg,nfile,files,hnlev)
!-----------------------------------------------------------------------
! Read the model atoms (files; built-in hydrogen with hnlev levels if
! hnlev > 0) and set up the cell arrays and macro-atom structure.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell,gncell,ng,nfile,hnlev
  real*8,intent(in) :: wlg(ng+1)       !group boundaries [cm]
  character(*),intent(in) :: files(nfile)
  type(nxion),allocatable :: ions(:)
  integer :: nion,i,k,ig,j,iel,s,lv,n,nst,ic,is
  integer,allocatable :: cnt(:)
!
  call nxt_dealloc
  nion = nfile + merge(1,0,hnlev>0)
  allocate(ions(nion))
  i = 0
  if(hnlev>0) then
     i = 1
     call nxa_builtin_h(hnlev,ions(1))
  endif
  do k=1,nfile
     i = i + 1
     call nxa_read(files(k),ions(i))
  enddo
  call nxa_setup(nion,ions)
  deallocate(ions)
  call nxa_groups(ng,wlg)
!
  nxt_on = .true.
  nxt_ncell = ncell
  nxt_ng = ng
  nst = nxa_nst
  allocate(nxt_x(nst,ncell),nxt_nel(nxa_nel,ncell),nxt_ne(ncell))
  allocate(nxt_kq(nxa_nlin,ncell),nxt_beta(nxa_nlin,ncell))
  allocate(nxt_jg(ng,ncell),nxt_jb(nxa_nlin,ncell),nxt_ljb(nxa_nlin,ncell))
  allocate(nxt_jbsum(nxa_nlin,ncell),nxt_jbcnt(nxa_nlin,ncell))
  allocate(grd_capnbf(ng,ncell))
  allocate(nxt_xloc(nst,gncell),nxt_neloc(gncell),nxt_nelloc(nxa_nel,gncell),nxt_neothloc(gncell))
  allocate(gas_capnbf(ng,gncell))
  allocate(nxt_lgrp(nxa_nlin),nxt_isrt(nxa_nlin),nxt_nug(ng),nxt_dnug(ng))
  allocate(nxt_act(nxa_nlin),nxt_wlact(nxa_nlin))
  nxt_x = 0d0; nxt_nel = 0d0; nxt_ne = 0d0; nxt_kq = 0d0; nxt_beta = 1d0
  nxt_jg = 0d0; nxt_jb = 0d0; nxt_ljb = .false.; nxt_jbsum = 0d0; nxt_jbcnt = 0d0
  grd_capnbf = 0.; gas_capnbf = 0.
  nxt_xloc = 0d0; nxt_neloc = 0d0; nxt_nelloc = 0d0; nxt_neothloc = 0d0
!
  do ig=1,ng
     nxt_nug(ig) = 0.5d0*(cc/wlg(ig) + cc/wlg(ig+1))
     nxt_dnug(ig) = cc/wlg(ig) - cc/wlg(ig+1)
  enddo
  do k=1,nxa_nlin
     nxt_lgrp(k) = 0
     do ig=1,ng
        if(nxa_lwl(k) >= wlg(ig) .and. nxa_lwl(k) < wlg(ig+1)) then
           nxt_lgrp(k) = ig
           exit
        endif
     enddo
  enddo
!-- all lines active until the first tables are built
  nxt_nact = nxa_nlin
  do j=1,nxa_nlin
     nxt_act(j) = nxa_lsrt(j)
     nxt_wlact(j) = nxa_lwl(nxa_lsrt(j))
     nxt_isrt(nxa_lsrt(j)) = j
  enddo
!
!-- exit tables
  allocate(nxt_ecoff(nxa_nel))
  nxt_nec = 0
  do iel=1,nxa_nel
     nxt_ecoff(iel) = nxt_nec
     nxt_nec = nxt_nec + nxa_el(iel)%ns*(nxa_el(iel)%ns+1)
  enddo
  allocate(nxt_ecum(nxt_nec,ncell))
  nxt_ecum = 1d0
!-- two-photon list
  nxt_n2p = sum(nxa_ion(:)%ntwo)
  allocate(nxt_2pion(nxt_n2p),nxt_2pk(nxt_n2p),nxt_2pnu(nxt_n2p))
  n = 0
  do i=1,nxa_nion
     do k=1,nxa_ion(i)%ntwo
        n = n + 1
        nxt_2pion(n) = i
        nxt_2pk(n) = k
        nxt_2pnu(n) = (nxa_ion(i)%e(nxa_ion(i)%tu(k)) - nxa_ion(i)%e(nxa_ion(i)%tl(k)))/hh
     enddo
  enddo
!-- emission channels of each global state: lines (upper state), two-photon
!-- (upper state), recombination (continuum state)
  allocate(cnt(nst+1))
  cnt = 0
  do k=1,nxa_nlin
     is = nxa_el(nxa_lel(k))%xoff + nxa_lsu(k)
     cnt(is) = cnt(is) + 1
  enddo
  do n=1,nxt_n2p
     i = nxt_2pion(n)
     is = nxa_el(nxa_ion(i)%el)%xoff + nxa_ion(i)%soff + nxa_ion(i)%tu(nxt_2pk(n))
     cnt(is) = cnt(is) + 1
  enddo
  do lv=1,nxa_nlv
     is = nxa_el(nxa_lvel(lv))%xoff + nxa_lvc(lv)
     cnt(is) = cnt(is) + 1
  enddo
  allocate(nxt_chst(nst+1))
  nxt_chst(1) = 1
  do s=1,nst
     nxt_chst(s+1) = nxt_chst(s) + cnt(s)
  enddo
  nxt_nch = nxt_chst(nst+1) - 1
  allocate(nxt_chk(nxt_nch),nxt_chi(nxt_nch),nxt_chcum(nxt_nch,ncell))
  nxt_chcum = 1d0
  cnt = 0
  do k=1,nxa_nlin
     is = nxa_el(nxa_lel(k))%xoff + nxa_lsu(k)
     j = nxt_chst(is) + cnt(is)
     cnt(is) = cnt(is) + 1
     nxt_chk(j) = em_line
     nxt_chi(j) = k
  enddo
  do n=1,nxt_n2p
     i = nxt_2pion(n)
     is = nxa_el(nxa_ion(i)%el)%xoff + nxa_ion(i)%soff + nxa_ion(i)%tu(nxt_2pk(n))
     j = nxt_chst(is) + cnt(is)
     cnt(is) = cnt(is) + 1
     nxt_chk(j) = em_2ph
     nxt_chi(j) = n
  enddo
  do lv=1,nxa_nlv
     is = nxa_el(nxa_lvel(lv))%xoff + nxa_lvc(lv)
     j = nxt_chst(is) + cnt(is)
     cnt(is) = cnt(is) + 1
     nxt_chk(j) = em_bf
     nxt_chi(j) = lv
  enddo
!-- k-packet outcomes: thermal, bf cooling per level, exit per global state
  nxt_nko = 1 + nxa_nlv + nst
  allocate(nxt_kcum(nxt_nko,ncell))
  nxt_kcum = 1d0
!-- bf emission sampling: running maximum of sigma u^3 along each table
  allocate(nxt_cmax(sum([(size(nxa_ion(i)%ps),i=1,nxa_nion)])),nxt_cmoff(nxa_nion))
  n = 0
  do i=1,nxa_nion
     nxt_cmoff(i) = n
     do k=1,nxa_ion(i)%nlev
        j = nxa_ion(i)%pst(k)
        do ic=j,j+nxa_ion(i)%pn(k)-1
           n = n + 1
           nxt_cmax(n) = nxa_ion(i)%ps(ic)*nxa_ion(i)%pu(ic)**3
           if(ic>j) nxt_cmax(n) = max(nxt_cmax(n),nxt_cmax(n-1))
        enddo
     enddo
  enddo
!-- H and He for the non-thermal electrons and the benchmark
  nxt_isel = .false.
  nxt_elz = 0
  nxt_ionzs = 0
  allocate(nxt_etop(nxa_nion))
  do iel=1,nxa_nel
     nxt_isel(nxa_el(iel)%z) = .true.
     nxt_elz(nxa_el(iel)%z) = iel
     do j=1,nxa_el(iel)%nion
        i = nxa_el(iel)%ion(j)
        nxt_ionzs(nxa_el(iel)%z,nxa_ion(i)%stage) = i
        nxt_etop(i) = maxval(nxa_ion(i)%e)/(hh*cc)
     enddo
     i = nxa_el(iel)%ion(nxa_el(iel)%nion)
     if(nxa_ion(i)%stage < 30) nxt_ionzs(nxa_el(iel)%z,nxa_ion(i)%stage+1) = -1
  enddo
  do iel=1,nxa_nel
     i = nxa_el(iel)%ion(1)
     if(nxa_ion(i)%stage/=1) cycle
     if(nxa_el(iel)%z==1) then
        nxt_elh = iel
        call nt_targets(iel,nxt_sh1,nxt_sh2,nxt_shc)
        if(nxa_ion(i)%nlev>=3) then
           if(nxa_ion(i)%lin(2,3)>0) nxt_kha = nxa_ion(i)%loff + nxa_ion(i)%lin(2,3)
        endif
     elseif(nxa_el(iel)%z==2) then
        nxt_elhe = iel
        call nt_targets(iel,nxt_she1,nxt_she2,nxt_shec)
     endif
  enddo
contains
  subroutine nt_targets(je,s1,s2,sc)
!-- global states: ground of the neutral, upper level of its strongest line, continuum
    integer,intent(in) :: je
    integer,intent(out) :: s1,s2,sc
    integer :: ii,kk,ku
    ii = nxa_el(je)%ion(1)
    s1 = nxa_el(je)%xoff + 1
    ku = 0
    do kk=1,nxa_ion(ii)%nlin
       if(nxa_ion(ii)%ll(kk)/=1) cycle
       if(ku==0) then
          ku = kk
       elseif(nxa_ion(ii)%f(kk) > nxa_ion(ii)%f(ku)) then
          ku = kk
       endif
    enddo
    s2 = 0
    if(ku>0) s2 = nxa_el(je)%xoff + nxa_ion(ii)%lu(ku)
    sc = nxa_el(je)%xoff + nxa_lvc(nxa_ion(ii)%lvoff+1)
  end subroutine nt_targets
end subroutine nxt_init
!
!
subroutine nxt_dealloc
!-- free the cell arrays and tables (set up again by nxt_init)
  if(allocated(nxt_x)) deallocate(nxt_x,nxt_nel,nxt_ne,nxt_kq,nxt_beta,nxt_jg,nxt_jb,nxt_ljb, &
       nxt_jbsum,nxt_jbcnt,grd_capnbf,nxt_xloc,nxt_neloc,nxt_nelloc,nxt_neothloc,gas_capnbf, &
       nxt_lgrp,nxt_isrt,nxt_nug,nxt_dnug,nxt_act,nxt_wlact)
  if(allocated(nxt_ecoff)) deallocate(nxt_ecoff,nxt_ecum)
  if(allocated(nxt_2pion)) deallocate(nxt_2pion,nxt_2pk,nxt_2pnu)
  if(allocated(nxt_chst)) deallocate(nxt_chst,nxt_chk,nxt_chi,nxt_chcum)
  if(allocated(nxt_kcum)) deallocate(nxt_kcum)
  if(allocated(nxt_cmax)) deallocate(nxt_cmax,nxt_cmoff)
  if(allocated(nxt_etop)) deallocate(nxt_etop)
  nxt_elh = 0; nxt_elhe = 0; nxt_kha = 0
  nxt_sh1 = 0; nxt_sh2 = 0; nxt_shc = 0; nxt_she1 = 0; nxt_she2 = 0; nxt_shec = 0
end subroutine nxt_dealloc
!
!
pure subroutine nxt_shares(x,yhe,fh,fih,feh,fihe,fehe)
!-----------------------------------------------------------------------
! Shares of the energy of a fast electron (E >~ 100 eV) in H/He gas with
! ionized fraction x (Shull & van Steenberg 1985, for n_He/n_H = 1/12):
! heat fh, H ionization fih and excitation feh, He I ionization fihe and
! excitation fehe.  The He shares scale with yhe = (n_He/n_H)/(1/12)
! (0 when He is not in NLTE: its share then heats the gas).
!-----------------------------------------------------------------------
  real*8,intent(in) :: x,yhe
  real*8,intent(out) :: fh,fih,feh,fihe,fehe
  real*8 :: xx,s
  xx = min(max(x,1d-4),1d0)
  call nlteh_svs(xx,fh,fih,feh)
  fihe = yhe*0.0554d0*(1d0-xx**0.4614d0)**1.6660d0
  fehe = yhe*0.0246d0*(1d0-xx**0.4049d0)**1.6594d0
  s = fih + feh + fihe + fehe
  if(s > 0.99d0) then
     fih = fih*0.99d0/s; feh = feh*0.99d0/s
     fihe = fihe*0.99d0/s; fehe = fehe*0.99d0/s
  endif
  fh = 1d0 - fih - feh - fihe - fehe
end subroutine nxt_shares
!
!
pure function nxt_yhe(nel) result(y)
!-- He/H relative to Shull & van Steenberg's 1/12 (0 without NLTE He)
  real*8,intent(in) :: nel(nxa_nel)
  real*8 :: y
  y = 0d0
  if(nxt_elhe>0 .and. nxt_elh>0) then
     if(nel(nxt_elh) > 0d0) y = (nel(nxt_elhe)/nel(nxt_elh))*12d0
  endif
end function nxt_yhe
!
!
pure subroutine nxt_ntrates(nel,x,edep,nti,ntx)
!-----------------------------------------------------------------------
! Non-thermal ionization (nti) and excitation (ntx) rates per particle of
! the ground states of H I and He I [1/s] from the power edep [erg/s/cm^3]
! deposited by fast electrons.
!-----------------------------------------------------------------------
  real*8,intent(in) :: nel(nxa_nel),x(nxa_nst),edep
  real*8,intent(out) :: nti(nxa_nst),ntx(nxa_nst)
  real*8 :: fh,fih,feh,fihe,fehe,n1,xion,eex
  integer :: i
  nti = 0d0
  ntx = 0d0
  if(edep <= 0d0 .or. nxt_elh==0) return
  xion = x(nxa_el(nxt_elh)%xoff+nxa_el(nxt_elh)%ns)
  call nxt_shares(xion,nxt_yhe(nel),fh,fih,feh,fihe,fehe)
  i = nxa_el(nxt_elh)%ion(1)
  n1 = nel(nxt_elh)*x(nxt_sh1)
  if(n1 > 1d-30) then
     nti(nxt_sh1) = fih*edep/(nxa_ion(i)%chi*n1)
     if(nxt_sh2>0) then
        eex = nxa_el(nxt_elh)%eps(nxt_sh2-nxa_el(nxt_elh)%xoff)
        ntx(nxt_sh1) = feh*edep/(eex*n1)
     endif
  endif
  if(nxt_elhe>0) then
     i = nxa_el(nxt_elhe)%ion(1)
     n1 = nel(nxt_elhe)*x(nxt_she1)
     if(n1 > 1d-30) then
        nti(nxt_she1) = fihe*edep/(nxa_ion(i)%chi*n1)
        if(nxt_she2>0) then
           eex = nxa_el(nxt_elhe)%eps(nxt_she2-nxa_el(nxt_elhe)%xoff)
           ntx(nxt_she1) = fehe*edep/(eex*n1)
        endif
     endif
  endif
end subroutine nxt_ntrates
!
!
subroutine nxt_gas_solve(it,gncell,temp,trad,nelv,neoth,qb,qv,icell1,edep)
!-----------------------------------------------------------------------
! Statistical equilibrium of the NLTE elements in this rank's gas cells.
! The first iteration uses J = B(trad); later ones the Monte Carlo
! estimators of the previous iteration, and the new populations are mixed
! geometrically with the old ones (weight nxt_damp).  edep: power absorbed
! from photons above nxt_ntemin [erg/s/cm^3], deposited by fast electrons.
!-----------------------------------------------------------------------
  integer,intent(in) :: it,gncell,icell1
  real*8,intent(in) :: temp(gncell),trad(gncell),nelv(nxa_nel,gncell),neoth(gncell)
  real*8,intent(in) :: edep(gncell),qb(gncell),qv(gncell)
  integer :: i,l,ig,ierr,nfail
  real*8 :: x(nxa_nst),xold(nxa_nst),ne,beta(nxa_nlin),jg(nxt_ng),jb(nxa_nlin)
  real*8 :: nti(nxa_nst),ntx(nxa_nst),ntot
  logical :: ljb(nxa_nlin)
!
  nfail = 0
  do i=1,gncell
     l = icell1 + i - 1
     ntot = sum(nelv(:,i))
     nxt_nelloc(:,i) = nelv(:,i)
     if(ntot <= 0d0) then
        nxt_xloc(:,i) = 0d0
        nxt_neloc(i) = neoth(i)
        cycle
     endif
     if(it==1 .or. .not.nxt_have_mc) then
        do ig=1,nxt_ng
           jg(ig) = bplanck(nxt_nug(ig),trad(i))
        enddo
        jb = 0d0
        ljb = .false.
        call nx_lte(temp(i),max(neoth(i),1d-3*ntot),x)
     else
        jg = nxt_jg(:,l)
        jb = nxt_jb(:,l)
        ljb = nxt_ljb(:,l)
        x = nxt_xloc(:,i)
     endif
     xold = x
     call nxt_ntrates(nelv(:,i),x,edep(i),nti,ntx)
     call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
          nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr)
     if(ierr==1) then
        nfail = nfail + 1
        x = xold
     elseif(it>1 .and. nxt_have_mc .and. nxt_damp<1d0) then
        x = max(xold,1d-300)**(1d0-nxt_damp)*x**nxt_damp
        call renorm(x)
     endif
     nxt_xloc(:,i) = x
     nxt_neloc(i) = neoth(i) + free_el(nelv(:,i),x)
  enddo
  if(nfail>0) write(0,*) 'nxt_gas_solve: no solution in cells:',nfail
contains
  subroutine renorm(xx)
    real*8,intent(inout) :: xx(nxa_nst)
    integer :: je,jo,ns
    do je=1,nxa_nel
       jo = nxa_el(je)%xoff
       ns = nxa_el(je)%ns
       xx(jo+1:jo+ns) = xx(jo+1:jo+ns)/sum(xx(jo+1:jo+ns))
    enddo
  end subroutine renorm
end subroutine nxt_gas_solve
!
!
pure function free_el(nel,x) result(ne)
!-- free electrons from the NLTE elements [cm^-3]
  real*8,intent(in) :: nel(nxa_nel),x(nxa_nst)
  real*8 :: ne
  integer :: je,jo
  ne = 0d0
  do je=1,nxa_nel
     jo = nxa_el(je)%xoff
     ne = ne + nel(je)*sum(nxa_el(je)%chg*x(jo+1:jo+nxa_el(je)%ns))
  enddo
end function free_el
!
!
subroutine nxt_store_jg(ncell,ng,jrad,vol)
!-- MC path-length tally (sum of e*d per group) to mean intensity per Hz per sr
  integer,intent(in) :: ncell,ng
  real*8,intent(in) :: jrad(ng,ncell),vol(ncell)
  integer :: i
  do i=1,ncell
     nxt_jg(:,i) = jrad(:,i)/(4d0*pi*vol(i)*nxt_dnug)
  enddo
end subroutine nxt_store_jg
!
!
subroutine nxt_tally_reset
  nxt_jbsum = 0d0
  nxt_jbcnt = 0d0
end subroutine nxt_tally_reset
!
!
subroutine nxt_tally_finish(ncell,vol)
!-- resonance tally (sum of e c/(nu |Q|)) to mean intensity incident on the line
  integer,intent(in) :: ncell
  real*8,intent(in) :: vol(ncell)
  integer :: i
  do i=1,ncell
     nxt_jb(:,i) = nxt_jbsum(:,i)/(4d0*pi*vol(i))
     nxt_ljb(:,i) = nxt_jbcnt(:,i) >= nxt_njbmin
  enddo
  nxt_have_mc = .true.
end subroutine nxt_tally_finish
!
!
subroutine nxt_tables(ncell,tempinv,radtemp,capemitgrey,xarr,vxarr)
!-----------------------------------------------------------------------
! Sobolev strengths, escape probabilities, macro-atom and k-packet tables
! for all cells, from the gathered populations.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: tempinv(ncell),radtemp(ncell),capemitgrey(ncell)
  real*8,intent(in) :: xarr(ncell+1),vxarr(ncell+1)
  integer :: ic,iel,ig,k,j,s,nfail,info,lv,l0,l1,ns,o,ko,ch0,ch1
  real*8 :: t,ne,qb,qv,ww,wth
  real*8 :: jg(nxt_ng),jb(nxa_nlin),jl(nxa_nlin),btab(0:nx_nbt)
  logical :: ljb(nxa_nlin)
  real*8 :: kw(nxt_nko)            !k-packet: non-returning flows
  real*8,allocatable :: cq(:,:),ci(:),c3(:),gam(:),gst(:),ar(:)
  real*8,allocatable :: q(:,:),nn(:,:),tot(:),emt(:),kex(:),kin(:),chw(:)
  integer,allocatable :: ipiv(:)
!
  nfail = 0
  do ic=1,ncell
     t = 1d0/tempinv(ic)
     ne = nxt_ne(ic)
     qb = (vxarr(ic+1) - vxarr(ic))/(xarr(ic+1) - xarr(ic))
     qv = (vxarr(ic+1) + vxarr(ic))/(xarr(ic+1) + xarr(ic))
     if(nxt_have_mc) then
        jg = nxt_jg(:,ic)
        jb = nxt_jb(:,ic)
        ljb = nxt_ljb(:,ic)
     else
        do ig=1,nxt_ng
           jg(ig) = bplanck(nxt_nug(ig),radtemp(ic))
        enddo
        jb = 0d0
        ljb = .false.
     endif
     call nx_beta_table(qb,qv,btab)
     do k=1,nxa_nlin
        jl(k) = nx_jline(k,nxt_ng,nxt_nug,jg,jb(k),ljb(k))
        o = nxa_el(nxa_lel(k))%xoff
        nxt_kq(k,ic) = nx_kq(k,nxt_nel(nxa_lel(k),ic),nxt_x(o+1:o+nxa_el(nxa_lel(k))%ns,ic))
        nxt_beta(k,ic) = nx_beta(nxt_kq(k,ic),qb,qv,btab)
     enddo
     if(nxt_tau0>0d0 .and. nxt_kha>0) then
!-- benchmark: only H-alpha, tau = tau0 for a gradient v/r
        nxt_kq(:,ic) = 0d0
        nxt_kq(nxt_kha,ic) = nxt_tau0*qv
     endif
!
     kw = 0d0
!-- thermal emission of the other opacities
     wth = arad*cc*t**4*capemitgrey(ic)
     kw(1) = wth
     do iel=1,nxa_nel
        ns = nxa_el(iel)%ns
        o = nxa_el(iel)%xoff
        allocate(cq(ns,ns),ci(ns),c3(ns),gam(ns),gst(ns),ar(ns))
        allocate(q(ns,ns),nn(ns,ns),tot(ns),emt(ns),kex(ns),kin(ns),ipiv(ns))
        call nx_el_rates(iel,t,radtemp(ic),nxt_ng,nxt_nug,nxt_dnug,jg,cq,ci,c3,gam,gst,ar)
        l0 = nxa_ion(nxa_el(iel)%ion(1))%loff + 1
        l1 = nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%loff + &
             nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%nlin
        ch0 = nxt_chst(o+1)
        ch1 = nxt_chst(o+ns+1) - 1
        allocate(chw(ch0:ch1))
        call el_chain(iel,ne,nxt_nel(iel,ic),nxt_x(o+1:o+ns,ic),cq,ci,c3,gam,gst,ar, &
             jl(l0:l1),nxt_beta(l0:l1,ic),q,tot,emt,kex,kin,chw,ch0)
!-- recombination cooling per level of this element
        do lv=1,nxa_nlv
           if(nxa_lvel(lv)/=iel) cycle
           kw(1+lv) = nxt_nel(iel,ic)*nxt_x(o+nxa_lvc(lv),ic)*ne*nxa_rcool(lv,t)
        enddo
!-- fundamental matrix N = (I-Q)^-1
        nn = 0d0
        do s=1,ns
           nn(s,s) = 1d0
        enddo
        q = -q
        do s=1,ns
           q(s,s) = q(s,s) + 1d0
        enddo
        call dgesv(ns,ns,q,ns,ipiv,nn,ns,info)
        if(info/=0) then
           nfail = nfail + 1
           nn = 0d0
           do s=1,ns
              nn(s,s) = 1d0
           enddo
        endif
        nn = max(nn,0d0)
!-- exit distributions per start state: emission from t, then the k-packet
        do s=1,ns
           j = nxt_ecoff(iel) + (s-1)*(ns+1)
           ww = 0d0
           do k=1,ns
              ww = ww + nn(s,k)*emt(k)
              nxt_ecum(j+k,ic) = ww
           enddo
           ww = ww + sum(nn(s,:)*kex)
           nxt_ecum(j+ns+1,ic) = ww
           if(ww > 0d0) then
              nxt_ecum(j+1:j+ns+1,ic) = nxt_ecum(j+1:j+ns+1,ic)/ww
           else
              nxt_ecum(j+1:j+ns+1,ic) = 1d0
           endif
           nxt_ecum(j+ns+1,ic) = 1d0
        enddo
!-- k-packet exits through this element's states: sum_s kin(s) N(s,t) emt(t)
        ko = 1 + nxa_nlv + o
        do k=1,ns
           kw(ko+k) = sum(kin*nn(:,k))*emt(k)
        enddo
!-- channel cumulatives
        do s=1,ns
           j = o + s
           ww = 0d0
           do k=nxt_chst(j),nxt_chst(j+1)-1
              ww = ww + chw(k)
              nxt_chcum(k,ic) = ww
           enddo
           if(ww > 0d0) then
              nxt_chcum(nxt_chst(j):nxt_chst(j+1)-1,ic) = nxt_chcum(nxt_chst(j):nxt_chst(j+1)-1,ic)/ww
           endif
           if(nxt_chst(j+1)>nxt_chst(j)) nxt_chcum(nxt_chst(j+1)-1,ic) = 1d0
        enddo
        deallocate(cq,ci,c3,gam,gst,ar,q,nn,tot,emt,kex,kin,ipiv,chw)
     enddo
!-- k-packet outcome distribution
     ww = 0d0
     do k=1,nxt_nko
        ww = ww + kw(k)
        nxt_kcum(k,ic) = ww
     enddo
     if(ww > 0d0) then
        nxt_kcum(:,ic) = nxt_kcum(:,ic)/ww
     else
        nxt_kcum(:,ic) = 1d0
     endif
     nxt_kcum(nxt_nko,ic) = 1d0
  enddo
  if(nfail>0) write(0,*) 'nxt_tables: singular macro-atom chain in cells:',nfail
!
!-- active lines: tau > nxt_taumin for the smallest projected gradient
  nxt_nact = 0
  nxt_isrt = 0
  do j=1,nxa_nlin
     k = nxa_lsrt(j)
     ww = 0d0
     do ic=1,ncell
        qb = (vxarr(ic+1) - vxarr(ic))/(xarr(ic+1) - xarr(ic))
        qv = (vxarr(ic+1) + vxarr(ic))/(xarr(ic+1) + xarr(ic))
        ww = max(ww, nxt_kq(k,ic)/max(min(abs(qb),abs(qv)),1d-3*max(abs(qb),abs(qv))))
     enddo
     if(ww > nxt_taumin) then
        nxt_nact = nxt_nact + 1
        nxt_act(nxt_nact) = k
        nxt_wlact(nxt_nact) = nxa_lwl(k)
        nxt_isrt(k) = nxt_nact
     endif
  enddo
end subroutine nxt_tables
!
!
subroutine el_chain(iel,ne,nel,x,cq,ci,c3,gam,gst,ar,jl,beta,q,tot,emt,kex,kin,chw,ch0)
!-----------------------------------------------------------------------
! Macro-atom transition probabilities of element iel in cell ic, from the
! energy-weighted rates (Lucy 2002): per state s with energy eps_s
!  internal jump down to l:  (R_sl + C_sl) eps_l
!  internal jump up to u:    (R_su + C_su) eps_s
!  emission (to l):          R_sl (eps_s - eps_l)
!  into the k-packet:        C_sl (eps_s - eps_l)
! normalized by their sum.  Out: q(s,s') jump probabilities, emt(s) and
! kex(s) total emission and k-packet probabilities of s, chw (channels of
! the element, from ch0) emission probabilities, kin(s) the k-packet
! entries into s (collisional cooling rates per volume), tot (scratch).
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,ch0
  real*8,intent(in) :: ne,nel,x(:),cq(:,:),ci(:),c3(:),gam(:),gst(:),ar(:),jl(:),beta(:)
  real*8,intent(out) :: q(:,:),tot(:),emt(:),kex(:),kin(:),chw(ch0:)
  integer :: ns,j,i,n,o,k,kg,kl,l,u,s,sc,lv,ch,g
  real*8 :: eps(nxa_el(iel)%ns),r,gl,gu,cor
!
  ns = nxa_el(iel)%ns
  eps = nxa_el(iel)%eps
  g = nxa_el(iel)%xoff
  q = 0d0
  kex = 0d0
  kin = 0d0
  chw = 0d0
  do j=1,nxa_el(iel)%nion
     i = nxa_el(iel)%ion(j)
     n = nxa_ion(i)%nlev
     o = nxa_ion(i)%soff
!-- collisions within the ion
     do l=1,n
        do u=1,n
           if(u==l) cycle
           r = ne*cq(o+l,o+u)
           if(r<=0d0) cycle
           if(eps(o+u)>eps(o+l)) then
!-- excitation: jump up; the k-packet excites it at the cooling rate
              q(o+l,o+u) = q(o+l,o+u) + r*eps(o+l)
              kin(o+u) = kin(o+u) + nel*x(o+l)*r*(eps(o+u) - eps(o+l))
           else
!-- de-excitation: jump down and into the k-packet
              q(o+l,o+u) = q(o+l,o+u) + r*eps(o+u)
              kex(o+l) = kex(o+l) + r*(eps(o+l) - eps(o+u))
           endif
        enddo
     enddo
!-- lines: emission and jump down, net absorption up
     do k=1,nxa_ion(i)%nlin
        kg = nxa_ion(i)%loff + k
        kl = kg - nxa_ion(nxa_el(iel)%ion(1))%loff
        l = o + nxa_ion(i)%ll(k)
        u = o + nxa_ion(i)%lu(k)
        r = nxa_la(kg)*beta(kl)
        q(u,l) = q(u,l) + r*eps(l)
        ch = find_ch(g+u,em_line,kg)
        chw(ch) = chw(ch) + r*(eps(u) - eps(l))
        gl = nxa_ion(i)%g(nxa_ion(i)%ll(k))
        gu = nxa_ion(i)%g(nxa_ion(i)%lu(k))
        cor = 1d0
        if(x(l) > 0d0) cor = max(0d0,1d0 - gl*x(u)/(gu*x(l)))
        q(l,u) = q(l,u) + beta(kl)*nxa_lblu(kg)*jl(kl)*cor*eps(l)
     enddo
!-- two-photon
     do k=1,nxa_ion(i)%ntwo
        l = o + nxa_ion(i)%tl(k)
        u = o + nxa_ion(i)%tu(k)
        r = nxa_ion(i)%ta(k)
        q(u,l) = q(u,l) + r*eps(l)
        ch = find_ch(g+u,em_2ph,two_index(i,k))
        chw(ch) = chw(ch) + r*(eps(u) - eps(l))
     enddo
!-- bound-free with the continuum state
     if(j<nxa_el(iel)%nion) then
        sc = nxa_ion(nxa_el(iel)%ion(j+1))%soff + 1
     else
        sc = ns
     endif
     do k=1,n
        s = o + k
        lv = nxa_ion(i)%lvoff + k
!-- up: photoionization less stimulated recombination, collisional ionization
        r = gam(s)
        if(x(s) > 0d0) r = max(0d0, r - ne*gst(s)*x(sc)/x(s))
        q(s,sc) = q(s,sc) + (r + ne*ci(s))*eps(s)
        kin(sc) = kin(sc) + nel*x(s)*ne*ci(s)*(eps(sc) - eps(s))
!-- down: radiative and three-body recombination
        q(sc,s) = q(sc,s) + ne*(ar(s) + ne*c3(s))*eps(s)
        ch = find_ch(g+sc,em_bf,lv)
        chw(ch) = chw(ch) + ne*ar(s)*(eps(sc) - eps(s))
        kex(sc) = kex(sc) + ne*ne*c3(s)*(eps(sc) - eps(s))
     enddo
  enddo
!-- normalize per state
  do s=1,ns
     r = 0d0
     do ch=nxt_chst(g+s),nxt_chst(g+s+1)-1
        r = r + chw(ch)
     enddo
     emt(s) = r
     tot(s) = sum(q(s,:)) + emt(s) + kex(s)
     if(tot(s) <= 0d0) then
!-- a state without any channel goes into the k-packet
        q(s,:) = 0d0
        emt(s) = 0d0
        kex(s) = 1d0
        do ch=nxt_chst(g+s),nxt_chst(g+s+1)-1
           chw(ch) = 0d0
        enddo
        cycle
     endif
     q(s,:) = q(s,:)/tot(s)
     emt(s) = emt(s)/tot(s)
     kex(s) = kex(s)/tot(s)
  enddo
!-- (channel weights stay unnormalized: their cumulative is normalized per state)
contains
  integer function find_ch(gs,kind,idx)
    integer,intent(in) :: gs,kind,idx
    integer :: c
    find_ch = 0
    do c=nxt_chst(gs),nxt_chst(gs+1)-1
       if(nxt_chk(c)==kind .and. nxt_chi(c)==idx) then
          find_ch = c
          return
       endif
    enddo
    stop 'el_chain: channel not found'
  end function find_ch
  integer function two_index(ii,kk)
    integer,intent(in) :: ii,kk
    integer :: m
    do m=1,nxt_n2p
       if(nxt_2pion(m)==ii .and. nxt_2pk(m)==kk) then
          two_index = m
          return
       endif
    enddo
    two_index = 0
    stop 'el_chain: two-photon transition not found'
  end function two_index
end subroutine el_chain
!
!
pure subroutine nxt_neighbors(wlc,iline,ilside,jlo,jhi)
!-----------------------------------------------------------------------
! Indices in the active line list of the lines just below (jlo) and above
! (jhi) the comoving wavelength wlc (0 if none).  A packet sitting exactly on
! line iline (active index) is on its red side if ilside>0, blue side otherwise.
!-----------------------------------------------------------------------
  real*8,intent(in) :: wlc
  integer,intent(in) :: iline,ilside
  integer,intent(out) :: jlo,jhi
  integer :: a,b,m
  if(iline>0) then
     if(ilside>0) then
        jlo = iline
        jhi = iline + 1
     else
        jlo = iline - 1
        jhi = iline
     endif
  else
     a = 0
     b = nxt_nact + 1
     do while(b-a>1)
        m = (a+b)/2
        if(nxt_wlact(m) <= wlc) then
           a = m
        else
           b = m
        endif
     enddo
     jlo = a
     jhi = a + 1
  endif
  if(jhi>nxt_nact) jhi = 0
  if(jlo<1) jlo = 0
end subroutine nxt_neighbors
!
!
pure subroutine pick(cum,r,i)
!-- first index with cum(i) >= r
  real*8,intent(in) :: cum(:),r
  integer,intent(out) :: i
  integer :: a,b,m
  a = 0
  b = size(cum)
  do while(b-a>1)
     m = (a+b)/2
     if(cum(m) >= r) then
        b = m
     else
        a = m
     endif
  enddo
  i = b
end subroutine pick
!
!
pure subroutine nxt_bffreq(lv,t,state,nu)
!-----------------------------------------------------------------------
! Frequency of a recombination photon to global level lv at temperature t:
! energy spectrum j_nu ~ sigma(nu) nu^3 e^{-h nu/kT}, by rejection from
! nu = nu0 + y kT/h, y ~ e^{-y}, with the running maximum of sigma u^3.
!-----------------------------------------------------------------------
  integer,intent(in) :: lv
  real*8,intent(in) :: t
  type(rnd_t),intent(inout) :: state
  real*8,intent(out) :: nu
  real*8 :: r,y,u,f,fmax,nu0
  integer :: i,k,it,j0,j1,jc,a,b,m,ic0
  i = nxa_lvion(lv)
  k = lv - nxa_ion(i)%lvoff
  nu0 = nxa_lvnu0(lv)
  j0 = nxa_ion(i)%pst(k)
  j1 = j0 + nxa_ion(i)%pn(k) - 1
!-- offset of this level's running maximum
  ic0 = nxt_cmoff(i) + (j0 - 1)
  do it=1,1000
     call rnd_r(r,state)
     y = -log(max(r,1d-300))
     u = 1d0 + y*kb*t/(hh*nu0)
     nu = nu0*u
     f = nxa_sigma(nxa_ion(i),k,nu)*u**3
!-- envelope: running maximum up to u (beyond the table: the last value, sigma u^3 const)
     if(u >= nxa_ion(i)%pu(j1)) then
        jc = j1
     else
        a = j0
        b = j1
        do while(b-a>1)
           m = (a+b)/2
           if(nxa_ion(i)%pu(m) <= u) then
              a = m
           else
              b = m
           endif
        enddo
        jc = b
     endif
     fmax = nxt_cmax(ic0 + jc - j0 + 1)
     call rnd_r(r,state)
     if(r*fmax <= f) exit
  enddo
end subroutine nxt_bffreq
!
!
pure subroutine emit_state(ic,gs,state,kind,idx,wlc)
!-- emission channel of the de-activating global state gs
  use gridmod, only:grd_tempinv
  integer,intent(in) :: ic,gs
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: r,nu
  integer :: c
  if(nxt_chst(gs+1)<=nxt_chst(gs)) then
!-- no emission channel (only reached through a degenerate table)
     kind = em_thermal
     idx = 0
     wlc = 0d0
     return
  endif
  call rnd_r(r,state)
  call pick(nxt_chcum(nxt_chst(gs):nxt_chst(gs+1)-1,ic),max(r,1d-300),c)
  c = nxt_chst(gs) + c - 1
  kind = nxt_chk(c)
  idx = nxt_chi(c)
  select case(kind)
  case(em_line)
     wlc = nxa_lwl(idx)
  case(em_bf)
     call nxt_bffreq(idx,1d0/grd_tempinv(ic),state,nu)
     wlc = cc/nu
  case default
     call rnd_r(r,state)
     wlc = cc/(nlteh_twophoton_y(r)*nxt_2pnu(idx))
  end select
end subroutine emit_state
!
!
pure subroutine nxt_kpacket(ic,state,kind,idx,wlc)
!-- a k-packet (thermal energy) in cell ic
  use gridmod, only:grd_tempinv
  integer,intent(in) :: ic
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: r,nu
  integer :: e
  call rnd_r(r,state)
  call pick(nxt_kcum(:,ic),max(r,1d-300),e)
  idx = 0
  wlc = 0d0
  if(e==1) then
     kind = em_thermal
  elseif(e<=1+nxa_nlv) then
     kind = em_bf
     idx = e - 1
     call nxt_bffreq(idx,1d0/grd_tempinv(ic),state,nu)
     wlc = cc/nu
  else
     call emit_state(ic,e-1-nxa_nlv,state,kind,idx,wlc)
  endif
end subroutine nxt_kpacket
!
!
pure subroutine nxt_activate(ic,iel,s,state,kind,idx,wlc)
!-- macro-atom of element iel activated in its state s
  integer,intent(in) :: ic,iel,s
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: r
  integer :: e,j,ns
  ns = nxa_el(iel)%ns
  j = nxt_ecoff(iel) + (s-1)*(ns+1)
  call rnd_r(r,state)
  call pick(nxt_ecum(j+1:j+ns+1,ic),max(r,1d-300),e)
  if(e>ns) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
  else
     call emit_state(ic,nxa_el(iel)%xoff+e,state,kind,idx,wlc)
  endif
end subroutine nxt_activate
!
!
pure subroutine nxt_activate_line(ic,k,state,kind,idx,wlc)
!-- packet absorbed in global line k: macro-atom in the upper level
  integer,intent(in) :: ic,k
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  call nxt_activate(ic,nxa_lel(k),nxa_lsu(k),state,kind,idx,wlc)
end subroutine nxt_activate_line
!
!
pure subroutine nxt_nonthermal(ic,state,kind,idx,wlc)
!-- energy of a fast photoelectron in cell ic: heat (k-packet), ionization
!-- (macro-atom in the continuum) or excitation of H I and He I with the
!-- Shull & van Steenberg (1985) shares
  integer,intent(in) :: ic
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: fh,fih,feh,fihe,fehe,r,xion
  if(nxt_elh==0) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
     return
  endif
  xion = nxt_x(nxa_el(nxt_elh)%xoff+nxa_el(nxt_elh)%ns,ic)
  call nxt_shares(xion,nxt_yhe(nxt_nel(:,ic)),fh,fih,feh,fihe,fehe)
  call rnd_r(r,state)
  if(r < fh) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
  elseif(r < fh + fih) then
     call nxt_activate(ic,nxt_elh,nxt_shc-nxa_el(nxt_elh)%xoff,state,kind,idx,wlc)
  elseif(r < fh + fih + feh .and. nxt_sh2>0) then
     call nxt_activate(ic,nxt_elh,nxt_sh2-nxa_el(nxt_elh)%xoff,state,kind,idx,wlc)
  elseif(r < fh + fih + feh + fihe .and. nxt_elhe>0) then
     call nxt_activate(ic,nxt_elhe,nxt_shec-nxa_el(nxt_elhe)%xoff,state,kind,idx,wlc)
  elseif(nxt_elhe>0 .and. nxt_she2>0) then
     call nxt_activate(ic,nxt_elhe,nxt_she2-nxa_el(nxt_elhe)%xoff,state,kind,idx,wlc)
  else
     call nxt_kpacket(ic,state,kind,idx,wlc)
  endif
end subroutine nxt_nonthermal
!
!
pure function nxt_ntdep(ng,nug,dnug,jg,cap) result(edep)
!-- power absorbed from photons above nxt_ntemin [erg/s/cm^3]
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng)
  real*4,intent(in) :: cap(ng)
  real*8 :: edep
  integer :: ig
  edep = 0d0
  if(nxt_ntemin <= 0d0) return
  do ig=1,ng
     if(hh*nug(ig) > nxt_ntemin) edep = edep + cap(ig)*jg(ig)*dnug(ig)
  enddo
  edep = 4d0*pi*edep
end function nxt_ntdep
!
!
pure subroutine nxt_bfabs(ic,ig,nu,state,kind,idx,wlc)
!-----------------------------------------------------------------------
! NLTE bound-free absorption of a packet of comoving frequency nu (group
! ig) in cell ic: choose the level by n_l <sigma_l>_g, then activate its
! continuum state with probability nu_l/nu (the ionization energy) or the
! thermal pool (the photo-electron's kinetic energy).
!-----------------------------------------------------------------------
  integer,intent(in) :: ic,ig
  real*8,intent(in) :: nu
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: s,r,w
  integer :: lv,iel,lsel
  s = 0d0
  do lv=1,nxa_nlv
     if(ig > nxa_g1(lv)) cycle
     iel = nxa_lvel(lv)
     s = s + nxt_nel(iel,ic)*nxt_x(nxa_el(iel)%xoff+nxa_lvs(lv),ic)*nxa_sgm(ig,lv)
  enddo
  if(s<=0d0) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
     return
  endif
  call rnd_r(r,state)
  r = r*s
  w = 0d0
  lsel = 0
  do lv=1,nxa_nlv
     if(ig > nxa_g1(lv)) cycle
     iel = nxa_lvel(lv)
     w = w + nxt_nel(iel,ic)*nxt_x(nxa_el(iel)%xoff+nxa_lvs(lv),ic)*nxa_sgm(ig,lv)
     lsel = lv
     if(w >= r) exit
  enddo
  call rnd_r(r,state)
  if(r < nxa_lvnu0(lsel)/nu) then
     call nxt_activate(ic,nxa_lvel(lsel),nxa_lvc(lsel),state,kind,idx,wlc)
  else
     call nxt_kpacket(ic,state,kind,idx,wlc)
  endif
end subroutine nxt_bfabs
!
!
subroutine nxt_write(ncell,rcell,vcell,tempinv,radtemp)
!-----------------------------------------------------------------------
! Per-cell NLTE diagnostics: output.nltex (ion fractions of each element,
! departure coefficients of the first levels of each ion, line depths of
! selected lines) and, for hydrogen, output.nlteh in the format of the
! H-only treatment.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: rcell(ncell),vcell(ncell),tempinv(ncell),radtemp(ncell)
  integer :: ic,iel,j,i,o,k,kla,kha,khb,ih
  real*8 :: t,xl(nxa_nst),qv,b(3),fr(nxa_maxion+1)
  character(16) :: hd(nxa_maxion+1)
  open(4,file='output.nltex',status='replace')
  write(4,'(a)') '# NLTE elements per cell: ion fractions (last = top ion), electron density'
  do iel=1,nxa_nel
     do j=1,nxa_el(iel)%nion
        hd(j) = nxa_ion(nxa_el(iel)%ion(j))%name
     enddo
     write(hd(nxa_el(iel)%nion+1),'(a,i0)') 'top_Z',nxa_el(iel)%z
     write(4,'(a,i3,a,99a12)') '# element Z =',nxa_el(iel)%z,':  r[cm] T[K] ne n_el ', &
          (trim(hd(j)),j=1,nxa_el(iel)%nion+1)
     do ic=1,ncell
        o = nxa_el(iel)%xoff
        do j=1,nxa_el(iel)%nion
           i = nxa_el(iel)%ion(j)
           fr(j) = sum(nxt_x(o+nxa_ion(i)%soff+1:o+nxa_ion(i)%soff+nxa_ion(i)%nlev,ic))
        enddo
        fr(nxa_el(iel)%nion+1) = nxt_x(o+nxa_el(iel)%ns,ic)
!-- below 1e-99 the e12.4 format drops the exponent letter
        where(fr<1d-99) fr = 0d0
        write(4,'(1p,99e12.4)') rcell(ic),1d0/tempinv(ic),nxt_ne(ic),nxt_nel(iel,ic), &
             (fr(j),j=1,nxa_el(iel)%nion+1)
     enddo
  enddo
  close(4)
!-- hydrogen in the format of the H-only treatment
  if(nxt_elh==0) return
  iel = nxt_elh
  ih = nxa_el(iel)%ion(1)
  o = nxa_el(iel)%xoff
  kla = nxa_ion(ih)%loff + nxa_ion(ih)%lin(1,2)
  kha = nxa_ion(ih)%loff + nxa_ion(ih)%lin(2,3)
  khb = nxa_ion(ih)%loff + nxa_ion(ih)%lin(2,4)
  open(4,file='output.nlteh',status='replace')
  write(4,'(a)') '# NLTE hydrogen per cell: b_n = departure coefficients relative to the'// &
       ' actual ion density (b = 1: LTE); tau = K/(v/r) (homologous equivalent)'
  write(4,'(a1,a11,19a12)') '#','r[cm]','v[cm/s]','T[K]','Trad[K]','nH','ne','x_HII','x_1','x_2','x_3', &
       'b_1','b_2','b_3','tau_Lya','tau_Ha','tau_Hb','beta_Ha','Jb_Ha/B','ljb_Ha'
  do ic=1,ncell
     t = 1d0/tempinv(ic)
     call nx_lte(t,nxt_ne(ic),xl)
     do k=1,3
        b(k) = (nxt_x(o+k,ic)/max(nxt_x(o+nxa_el(iel)%ns,ic),1d-300))/ &
             (xl(o+k)/max(xl(o+nxa_el(iel)%ns),1d-300))
     enddo
     qv = max(vcell(ic)/rcell(ic),1d-30)
     write(4,'(1p,19e12.4)') rcell(ic),vcell(ic),t,radtemp(ic),nxt_nel(iel,ic),nxt_ne(ic), &
          nxt_x(o+nxa_el(iel)%ns,ic),nxt_x(o+1,ic),nxt_x(o+2,ic),nxt_x(o+3,ic),b(1),b(2),b(3), &
          nxt_kq(kla,ic)/qv,nxt_kq(kha,ic)/qv,nxt_kq(khb,ic)/qv,nxt_beta(kha,ic), &
          nxt_jb(kha,ic)/max(bplanck(nxa_lnu(kha),t),1d-300),merge(1d0,0d0,nxt_ljb(kha,ic))
  enddo
  close(4)
end subroutine nxt_write
!
end module nltextrmod
! vim: fdm=marker
