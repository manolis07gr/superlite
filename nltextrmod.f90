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
! Super-levels (states of several levels, Boltzmann at T_e within them):
! the macro-atom jumps between states with the rates summed over their
! levels; the energy of a state is the mean of its levels.  The energy
! emitted or given to the electrons by a transition is that of the level
! pair (line h nu, ionization energy, collisional energy exchange), the
! jump keeps the rest: weight (rate x energy of the state - emitted).  A
! line within a super-level is thermal (source function B(T_e)): absorption
! goes to the k-packet, which emits it at the rate n_u A beta h nu.
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
  integer,public :: nxt_nko = 0               !k-packet outcomes: thermal, bf cooling per level, exit per
                                              !state, lines within super-levels
  real*8,allocatable,public :: nxt_kcum(:,:)  !(nko,ncell)
  integer,allocatable,public :: nxt_lch(:)    !(nlin) emission channel of the line (0: within a super-level)
  integer,allocatable,public :: nxt_bch(:)    !(nlv) recombination channel of the level
  integer,allocatable,public :: nxt_2pch(:)   !(n2p) channel of the two-photon transition
  integer,public :: nxt_nin = 0               !lines within super-levels
  integer,allocatable,public :: nxt_inl(:)    !(nin) their global index
  real*8,allocatable,public :: nxt_w(:,:)     !(nlv,ncell) share of each level in its state (nxa_weights)
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
!-- time-dependent statistical equilibrium (in_nltex_prev): 1/dt to the
!-- previous model (0: steady state), its populations at the same mass
!-- coordinate, and per cell the power released by the change of the
!-- ionization and excitation energy, sum n_el eps_s (x_prev - x)/dt V
!-- [erg/s] (negative: stored), and the share of the cell's 'lumg' energy
!-- emitted by macro-atom activation in the states that empty
  real*8,public :: nxt_dtinv = 0d0
  real*8,allocatable,public :: nxt_xprev(:,:) !(nst,ncell)
  real*8,allocatable,public :: nxt_erel(:)    !(ncell)
  real*8,allocatable,public :: nxt_prel(:)    !(ncell)
  real*8,allocatable :: nxt_rcum(:,:)         !(nst,ncell) cumulative weight of the activation states
!-- microclumping: volume filling factor per cell (not allocated: smooth).
!-- The populations are those of the clumps; the Sobolev strengths K of
!-- the transport and of the escape probabilities are volume means, f K
  real*8,allocatable,public :: nxt_fcl(:)    !(ncell)
!-- emission kinds returned by the macro-atom
  integer,parameter,public :: em_line = 1, em_bf = 2, em_thermal = 3, em_2ph = 4
!
  public :: nxt_init,nxt_gas_solve,nxt_store_jg,nxt_tally_reset,nxt_tally_finish
  public :: nxt_tables,nxt_neighbors,nxt_activate,nxt_activate_line,nxt_kpacket
  public :: nxt_bfabs,nxt_nonthermal,nxt_ntdep,nxt_ntrates,nxt_write,nxt_bffreq
  public :: nxt_shares,nxt_dealloc,nxt_emcode
  public :: nxt_write_state,nxt_read_prev,nxt_release,nxt_activate_rel
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
  allocate(nxt_x(nst,ncell),nxt_nel(nxa_nel,ncell),nxt_ne(ncell),nxt_w(nxa_nlv,ncell))
  nxt_w = 1d0
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
     if(nxa_lwl(k) < wlg(1) .or. nxa_lwl(k) >= wlg(ng+1)) cycle
     i = 1
     j = ng + 1
     do while(j-i>1)
        n = (i+j)/2
        if(wlg(n) <= nxa_lwl(k)) then
           i = n
        else
           j = n
        endif
     enddo
     nxt_lgrp(k) = i
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
!-- emission channels of each global state: lines (upper state; not those
!-- within a super-level), two-photon (upper state), recombination
!-- (continuum state)
  allocate(cnt(nst+1),nxt_lch(nxa_nlin),nxt_bch(nxa_nlv),nxt_2pch(nxt_n2p))
  cnt = 0
  nxt_nin = 0
  do k=1,nxa_nlin
     if(nxa_lsl(k)==nxa_lsu(k)) then
        nxt_nin = nxt_nin + 1
        cycle
     endif
     is = nxa_el(nxa_lel(k))%xoff + nxa_lsu(k)
     cnt(is) = cnt(is) + 1
  enddo
  do n=1,nxt_n2p
     i = nxt_2pion(n)
     if(nxa_ion(i)%sl(nxa_ion(i)%tl(nxt_2pk(n)))==nxa_ion(i)%sl(nxa_ion(i)%tu(nxt_2pk(n)))) &
          stop 'nxt_init: two-photon transition within a super-level'
     is = nxa_el(nxa_ion(i)%el)%xoff + nxa_ion(i)%soff + nxa_ion(i)%sl(nxa_ion(i)%tu(nxt_2pk(n)))
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
  allocate(nxt_chk(nxt_nch),nxt_chi(nxt_nch),nxt_chcum(nxt_nch,ncell),nxt_inl(nxt_nin))
  nxt_chcum = 1d0
  cnt = 0
  nxt_nin = 0
  do k=1,nxa_nlin
     if(nxa_lsl(k)==nxa_lsu(k)) then
        nxt_nin = nxt_nin + 1
        nxt_inl(nxt_nin) = k
        nxt_lch(k) = 0
        cycle
     endif
     is = nxa_el(nxa_lel(k))%xoff + nxa_lsu(k)
     j = nxt_chst(is) + cnt(is)
     cnt(is) = cnt(is) + 1
     nxt_chk(j) = em_line
     nxt_chi(j) = k
     nxt_lch(k) = j
  enddo
  do n=1,nxt_n2p
     i = nxt_2pion(n)
     is = nxa_el(nxa_ion(i)%el)%xoff + nxa_ion(i)%soff + nxa_ion(i)%sl(nxa_ion(i)%tu(nxt_2pk(n)))
     j = nxt_chst(is) + cnt(is)
     cnt(is) = cnt(is) + 1
     nxt_chk(j) = em_2ph
     nxt_chi(j) = n
     nxt_2pch(n) = j
  enddo
  do lv=1,nxa_nlv
     is = nxa_el(nxa_lvel(lv))%xoff + nxa_lvc(lv)
     j = nxt_chst(is) + cnt(is)
     cnt(is) = cnt(is) + 1
     nxt_chk(j) = em_bf
     nxt_chi(j) = lv
     nxt_bch(lv) = j
  enddo
!-- k-packet outcomes: thermal, bf cooling per level, exit per global state,
!-- lines within super-levels
  nxt_nko = 1 + nxa_nlv + nst + nxt_nin
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
    if(ku>0) s2 = nxa_el(je)%xoff + nxa_ion(ii)%soff + nxa_ion(ii)%sl(nxa_ion(ii)%lu(ku))
    sc = nxa_el(je)%xoff + nxa_lvc(nxa_ion(ii)%lvoff+1)
  end subroutine nt_targets
end subroutine nxt_init
!
!
subroutine nxt_dealloc
!-- free the cell arrays and tables (set up again by nxt_init)
  if(allocated(nxt_x)) deallocate(nxt_x,nxt_nel,nxt_ne,nxt_kq,nxt_beta,nxt_jg,nxt_jb,nxt_ljb, &
       nxt_jbsum,nxt_jbcnt,grd_capnbf,nxt_xloc,nxt_neloc,nxt_nelloc,nxt_neothloc,gas_capnbf, &
       nxt_lgrp,nxt_isrt,nxt_nug,nxt_dnug,nxt_act,nxt_wlact,nxt_w)
  if(allocated(nxt_lch)) deallocate(nxt_lch,nxt_bch,nxt_2pch)
  if(allocated(nxt_inl)) deallocate(nxt_inl)
  nxt_nin = 0
  if(allocated(nxt_ecoff)) deallocate(nxt_ecoff,nxt_ecum)
  if(allocated(nxt_2pion)) deallocate(nxt_2pion,nxt_2pk,nxt_2pnu)
  if(allocated(nxt_chst)) deallocate(nxt_chst,nxt_chk,nxt_chi,nxt_chcum)
  if(allocated(nxt_kcum)) deallocate(nxt_kcum)
  if(allocated(nxt_cmax)) deallocate(nxt_cmax,nxt_cmoff)
  if(allocated(nxt_etop)) deallocate(nxt_etop)
  if(allocated(nxt_xprev)) deallocate(nxt_xprev)
  if(allocated(nxt_fcl)) deallocate(nxt_fcl)
  if(allocated(nxt_erel)) deallocate(nxt_erel,nxt_prel,nxt_rcum)
  nxt_dtinv = 0d0
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
     if(nxt_dtinv>0d0) then
        call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
             nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,xprev=nxt_xprev(:,l),dtinv=nxt_dtinv)
     else
        call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
             nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr)
     endif
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
  integer :: ic,iel,ig,k,j,s,nfail,info,lv,l0,l1,ns,o,ko,ch0,ch1,lo,m
  real*8 :: t,ne,qb,qv,ww,wth
  real*8 :: jg(nxt_ng),jb(nxa_nlin),jl(nxa_nlin),btab(0:nx_nbt)
  logical :: ljb(nxa_nlin)
  real*8 :: kw(nxt_nko)            !k-packet: non-returning flows
  type(nxrt) :: r
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
     call nxa_weights(t,nxt_w(:,ic))
     do k=1,nxa_nlin
        jl(k) = nx_jline(k,nxt_ng,nxt_nug,jg,jb(k),ljb(k))
        iel = nxa_lel(k)
        o = nxa_el(iel)%xoff
        lo = nxa_el(iel)%lvoff
        nxt_kq(k,ic) = nx_kq(k,nxt_nel(iel,ic),nxt_x(o+1:o+nxa_el(iel)%ns,ic), &
             nxt_w(lo+1:lo+nxa_el(iel)%nlv,ic))
        if(allocated(nxt_fcl)) nxt_kq(k,ic) = nxt_fcl(ic)*nxt_kq(k,ic)
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
        lo = nxa_el(iel)%lvoff
        allocate(q(ns,ns),nn(ns,ns),tot(ns),emt(ns),kex(ns),kin(ns),ipiv(ns))
        call nx_el_rates(iel,t,radtemp(ic),nxt_ng,nxt_nug,nxt_dnug,jg,r)
        l0 = nxa_ion(nxa_el(iel)%ion(1))%loff + 1
        l1 = nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%loff + &
             nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%nlin
        ch0 = nxt_chst(o+1)
        ch1 = nxt_chst(o+ns+1) - 1
        allocate(chw(ch0:ch1))
        call el_chain(iel,ne,nxt_nel(iel,ic),nxt_x(o+1:o+ns,ic),r, &
             jl(l0:l1),nxt_beta(l0:l1,ic),q,tot,emt,kex,kin,chw,ch0)
!-- recombination cooling per level of this element (from level 1 of the next ion)
        do k=1,nxa_el(iel)%nlv
           lv = lo + k
           kw(1+lv) = nxt_nel(iel,ic)*nxt_x(o+nxa_lvc(lv),ic)*r%wc(k)*ne*nxa_rcool(lv,t)
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
        deallocate(q,nn,tot,emt,kex,kin,ipiv,chw)
     enddo
!-- lines within super-levels: thermal emission n_u A beta h nu
     ko = 1 + nxa_nlv + nxa_nst
     do m=1,nxt_nin
        k = nxt_inl(m)
        iel = nxa_lel(k)
        kw(ko+m) = nxt_nel(iel,ic)*nxt_x(nxa_el(iel)%xoff+nxa_lsu(k),ic)*nxt_w(nxa_luv(k),ic)* &
             nxa_la(k)*nxt_beta(k,ic)*hh*nxa_lnu(k)
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
subroutine el_chain(iel,ne,nel,x,r,jl,beta,q,tot,emt,kex,kin,chw,ch0)
!-----------------------------------------------------------------------
! Macro-atom transition probabilities of element iel in cell ic, from the
! energy-weighted rates (Lucy 2002), per state s with energy eps_s (r%eps):
!  a transition out of s with rate R that emits (or gives the electrons)
!  the energy E: emission (k-packet) weight R E, internal jump R (eps_s - E)
!  (Lucy: E = eps_s - eps_l for a downward jump to l, 0 upward), with E of
!  the level pair: line h nu, ionization energy h nu_0, the collisional
!  exchange (r%edn).
! normalized by their sum.  Out: q(s,s') jump probabilities, emt(s) and
! kex(s) total emission and k-packet probabilities of s, chw (channels of
! the element, from ch0) emission probabilities, kin(s) the k-packet
! entries into s (collisional cooling rates per volume), tot (scratch).
! Lines within a super-level are not part of the macro-atom.
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,ch0
  real*8,intent(in) :: ne,nel,x(:),jl(:),beta(:)
  type(nxrt),intent(in) :: r
  real*8,intent(out) :: q(:,:),tot(:),emt(:),kex(:),kin(:),chw(ch0:)
  integer :: ns,j,i,o,k,kg,kl,l,u,s,sc,lv,ch,lv0,k0,lo,sl
  real*8 :: eps(nxa_el(iel)%ns),rr,cor,hnu,chi,nl,nu_
!
  ns = nxa_el(iel)%ns
  eps = r%eps
  lv0 = nxa_el(iel)%lvoff
  k0 = nxa_ion(nxa_el(iel)%ion(1))%loff
  q = 0d0
  kex = 0d0
  kin = 0d0
  chw = 0d0
!-- collisions between the states
  do l=1,ns
     do u=1,ns
        if(u==l) cycle
        rr = ne*r%cq(l,u)
        if(rr<=0d0) cycle
        q(l,u) = q(l,u) + max(rr*eps(l) - ne*r%edn(l,u),0d0)
        kex(l) = kex(l) + ne*r%edn(l,u)
        kin(u) = kin(u) + nel*x(l)*ne*r%eup(l,u)
     enddo
  enddo
  do j=1,nxa_el(iel)%nion
     i = nxa_el(iel)%ion(j)
     o = nxa_ion(i)%soff
     lo = nxa_ion(i)%lvoff - lv0
!-- lines: emission and jump down, net absorption up
     do k=1,nxa_ion(i)%nlin
        kg = nxa_ion(i)%loff + k
        kl = kg - k0
        l = nxa_lsl(kg)
        u = nxa_lsu(kg)
        if(l==u) cycle
        hnu = hh*nxa_lnu(kg)
        rr = nxa_la(kg)*beta(kl)*r%w(nxa_luv(kg)-lv0)
        q(u,l) = q(u,l) + max(rr*(eps(u) - hnu),0d0)
        ch = nxt_lch(kg)
        chw(ch) = chw(ch) + rr*hnu
        nl = x(l)*r%w(nxa_llv(kg)-lv0)
        cor = 1d0
        if(nl > 0d0) cor = max(0d0,1d0 - nxa_ion(i)%g(nxa_ion(i)%ll(k))*x(u)*r%w(nxa_luv(kg)-lv0)/ &
             (nxa_ion(i)%g(nxa_ion(i)%lu(k))*nl))
        q(l,u) = q(l,u) + r%w(nxa_llv(kg)-lv0)*beta(kl)*nxa_lblu(kg)*jl(kl)*cor*eps(l)
     enddo
!-- two-photon
     do k=1,nxa_ion(i)%ntwo
        l = o + nxa_ion(i)%sl(nxa_ion(i)%tl(k))
        u = o + nxa_ion(i)%sl(nxa_ion(i)%tu(k))
        hnu = nxa_ion(i)%e(nxa_ion(i)%tu(k)) - nxa_ion(i)%e(nxa_ion(i)%tl(k))
        rr = nxa_ion(i)%ta(k)*r%w(lo+nxa_ion(i)%tu(k))
        q(u,l) = q(u,l) + max(rr*(eps(u) - hnu),0d0)
        ch = nxt_2pch(two_index(i,k))
        chw(ch) = chw(ch) + rr*hnu
     enddo
!-- bound-free with the continuum state (level 1 of the next ion)
     sc = nxa_lvc(nxa_ion(i)%lvoff+1)
     do s=o+1,o+nxa_ion(i)%nsl
!-- up: photoionization less stimulated recombination, collisional ionization
        rr = r%gam(s)
        if(x(s) > 0d0) rr = max(0d0, rr - ne*r%gst(s)*x(sc)/x(s))
        q(s,sc) = q(s,sc) + (rr + ne*r%ci(s))*eps(s)
     enddo
     do k=1,nxa_ion(i)%nlev
        lv = nxa_ion(i)%lvoff + k
        s = nxa_lvs(lv)
        chi = hh*nxa_ion(i)%nu0(k)
        kin(sc) = kin(sc) + nel*x(s)*r%w(lo+k)*ne*r%cil(lo+k)*chi
!-- down: radiative and three-body recombination into level k
        rr = ne*r%wc(lo+k)*r%arl(lo+k)
        nu_ = ne*ne*r%wc(lo+k)*r%c3l(lo+k)
        q(sc,s) = q(sc,s) + max((rr + nu_)*(eps(sc) - chi),0d0)
        ch = nxt_bch(lv)
        chw(ch) = chw(ch) + rr*chi
        kex(sc) = kex(sc) + nu_*chi
     enddo
  enddo
!-- jumps within a state (super-level) do not change it
  do s=1,ns
     q(s,s) = 0d0
  enddo
!-- normalize per state
  do s=1,ns
     rr = 0d0
     do ch=nxt_chst(nxa_el(iel)%xoff+s),nxt_chst(nxa_el(iel)%xoff+s+1)-1
        rr = rr + chw(ch)
     enddo
     emt(s) = rr
     tot(s) = sum(q(s,:)) + emt(s) + kex(s)
     if(tot(s) <= 0d0) then
!-- a state without any channel goes into the k-packet
        q(s,:) = 0d0
        emt(s) = 0d0
        kex(s) = 1d0
        do ch=nxt_chst(nxa_el(iel)%xoff+s),nxt_chst(nxa_el(iel)%xoff+s+1)-1
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
pure function nxt_emcode(kind,idx) result(code)
!-- emission process for the diagnostics (packet%em): 2 thermal emissivity,
!-- 3 bound-free, 4 line, 5 two-photon, + 100 Z of the element
  integer,intent(in) :: kind,idx
  integer :: code
  code = 2
  if(idx <= 0) return
  select case(kind)
  case(em_line)
     code = 4 + 100*nxa_el(nxa_lel(idx))%z
  case(em_bf)
     code = 3 + 100*nxa_el(nxa_lvel(idx))%z
  case(em_2ph)
     code = 5 + 100*nxa_el(nxa_ion(nxt_2pion(idx))%el)%z
  end select
end function nxt_emcode
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
  elseif(e<=1+nxa_nlv+nxa_nst) then
     call emit_state(ic,e-1-nxa_nlv,state,kind,idx,wlc)
  else
!-- a line within a super-level
     kind = em_line
     idx = nxt_inl(e-1-nxa_nlv-nxa_nst)
     wlc = nxa_lwl(idx)
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
!-- packet absorbed in global line k: macro-atom in the upper level (its
!-- state); a line within a super-level is thermal (k-packet)
  integer,intent(in) :: ic,k
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  if(nxt_lch(k)==0) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
  else
     call nxt_activate(ic,nxa_lel(k),nxa_lsu(k),state,kind,idx,wlc)
  endif
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
     s = s + nxt_nel(iel,ic)*nxt_x(nxa_el(iel)%xoff+nxa_lvs(lv),ic)*nxt_w(lv,ic)*nxa_sgm(nxa_sgo(lv)+ig)
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
     w = w + nxt_nel(iel,ic)*nxt_x(nxa_el(iel)%xoff+nxa_lvs(lv),ic)*nxt_w(lv,ic)*nxa_sgm(nxa_sgo(lv)+ig)
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
  real*8 :: t,xl(nxa_nst),qv,b(3),fr(nxa_maxion+1),b8(8),trw(2)
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
           fr(j) = sum(nxt_x(o+nxa_ion(i)%soff+1:o+nxa_ion(i)%soff+nxa_ion(i)%nsl,ic))
        enddo
        fr(nxa_el(iel)%nion+1) = nxt_x(o+nxa_el(iel)%ns,ic)
!-- below 1e-99 the e12.4 format drops the exponent letter
        where(fr<1d-99) fr = 0d0
        write(4,'(1p,99e12.4)') rcell(ic),1d0/tempinv(ic),nxt_ne(ic),nxt_nel(iel,ic), &
             (fr(j),j=1,nxa_el(iel)%nion+1)
     enddo
  enddo
  close(4)
!-- Monte Carlo mean intensity in bands: 4 pi Int J dnu [erg/s/cm^2]
  open(4,file='output.nltexj',status='replace')
  write(4,'(a)') '# Monte Carlo radiation field per cell: 4 pi Int J dnu [erg/s/cm^2] in bands, '// &
       'T_R from the mean frequency and dilution W (Mazzali & Lucy 1993)'
  write(4,'(a1,a11,11a12)') '#','r[cm]','T[K]','<228A','228-504','504-912','912-2000', &
       '2000-4000','4000-1e4','1e4-3e4','>3e4A','T_R[K]','W'
  do ic=1,ncell
     call bands(nxt_jg(:,ic),b8)
     write(4,'(1p,12e12.4)') rcell(ic),1d0/tempinv(ic),b8,trw(1),trw(2)
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
contains
  subroutine bands(jg,bo)
    real*8,intent(in) :: jg(nxt_ng)
    real*8,intent(out) :: bo(8)
    real*8,parameter :: edge(9) = [0d0,228d-8,504d-8,912d-8,2000d-8,4000d-8,1d-4,3d-4,1d99]
    real*8 :: wl,j0,j1
    integer :: ig,ib
    bo = 0d0
    j0 = 0d0
    j1 = 0d0
    do ig=1,nxt_ng
       wl = cc/nxt_nug(ig)
       do ib=1,8
          if(wl>=edge(ib) .and. wl<edge(ib+1)) bo(ib) = bo(ib) + 4d0*pi*jg(ig)*nxt_dnug(ig)
       enddo
       j0 = j0 + jg(ig)*nxt_dnug(ig)
       j1 = j1 + nxt_nug(ig)*jg(ig)*nxt_dnug(ig)
    enddo
    trw = 0d0
    if(j0>0d0) then
       trw(1) = hh*(j1/j0)/(3.832d0*kb)
       trw(2) = pi*j0/(5.670400d-5*trw(1)**4)
    endif
  end subroutine bands
end subroutine nxt_write
!
!
subroutine nxt_write_state(ncell,mass,vol,tempinv)
!-----------------------------------------------------------------------
! output.nltex_state: the populations of all NLTE states per cell, with
! the exterior mass coordinate of the cell centre, for the time-dependent
! statistical equilibrium of a later model (in_nltex_prev).
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: mass(ncell),vol(ncell),tempinv(ncell)
  integer :: ic,iel
  real*8 :: mext(ncell)
  mext(ncell) = 0.5d0*mass(ncell)
  do ic=ncell-1,1,-1
     mext(ic) = mext(ic+1) + 0.5d0*(mass(ic+1) + mass(ic))
  enddo
  open(4,file='output.nltex_state',status='replace')
  write(4,'(a)') '# SuperLite NLTE state populations per cell (in_nltex_prev of a later model)'
  write(4,'(a,3i9)') 'nel nst ncell',nxa_nel,nxa_nst,ncell
  write(4,'(a,99(1x,i3,1x,i5))') 'elements Z ns',(nxa_el(iel)%z,nxa_el(iel)%ns,iel=1,nxa_nel)
  write(4,'(a)') '# per cell: exterior mass of the cell centre [g], mass [g], density [g/cm^3], '// &
       'T_e [K], n_e [cm^-3]; then the fractions of the states of each element'
  do ic=1,ncell
     write(4,'(1p,5es16.7e3)') mext(ic),mass(ic),mass(ic)/max(vol(ic),1d-300),1d0/tempinv(ic),nxt_ne(ic)
     write(4,'(1p,8es16.7e3)') max(nxt_x(:,ic),0d0)
  enddo
  close(4)
end subroutine nxt_write_state
!
!
subroutine nxt_read_prev(fname,ncell,mass,dt,lprint)
!-----------------------------------------------------------------------
! Populations of the previous model (its output.nltex_state) at the
! exterior mass coordinate of each cell centre, linear in that coordinate
! between the previous cell centres (the end values beyond), for the
! time-dependent statistical equilibrium with dt [s] since that model.
! The exterior mass is Lagrangian when both models end at the same mass
! shell (the outer edge of the CSM); the model atoms must be the same.
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname
  integer,intent(in) :: ncell
  real*8,intent(in) :: mass(ncell),dt
  logical,intent(in) :: lprint
  integer :: iu,ios,nel,nst,npc,ic,iel,a,b,m,k,o
  integer :: zs(200)
  character(4096) :: line
  real*8,allocatable :: mp(:),xp(:,:)
  real*8 :: mext(ncell),rd(5),w
!
  if(dt<=0d0) stop 'nxt_read_prev: in_nltex_dtprev must be > 0'
  open(newunit=iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) stop 'nxt_read_prev: cannot open in_nltex_prev'
  call nextrec
  if(index(line,'nel nst ncell')/=1) stop 'nxt_read_prev: not an output.nltex_state file'
  read(line(14:),*) nel,nst,npc
  if(nel/=nxa_nel .or. nst/=nxa_nst) stop 'nxt_read_prev: other model atoms than this model'
  call nextrec
  if(nel>100) stop "nxt_read_prev: too many elements"
  read(line(14:),*) zs(:2*nel)
  do iel=1,nel
     if(zs(2*iel-1)/=nxa_el(iel)%z .or. zs(2*iel)/=nxa_el(iel)%ns) &
          stop 'nxt_read_prev: other model atoms than this model'
  enddo
  allocate(mp(npc),xp(nst,npc))
  do ic=1,npc
     call nextrec
     read(line,*) rd
     read(iu,*) xp(:,ic)
     mp(ic) = rd(1)
  enddo
  close(iu)
!-- exterior mass of the cell centres of this model
  mext(ncell) = 0.5d0*mass(ncell)
  do ic=ncell-1,1,-1
     mext(ic) = mext(ic+1) + 0.5d0*(mass(ic+1) + mass(ic))
  enddo
  if(allocated(nxt_xprev)) deallocate(nxt_xprev)
  allocate(nxt_xprev(nst,ncell))
!-- mp decreases outward (index npc: the outermost cell)
  do ic=1,ncell
     if(mext(ic) >= mp(1)) then
        nxt_xprev(:,ic) = xp(:,1)
     elseif(mext(ic) <= mp(npc)) then
        nxt_xprev(:,ic) = xp(:,npc)
     else
        a = 1
        b = npc
        do while(b-a>1)
           m = (a+b)/2
           if(mp(m) >= mext(ic)) then
              a = m
           else
              b = m
           endif
        enddo
        w = (mp(a) - mext(ic))/(mp(a) - mp(b))
        nxt_xprev(:,ic) = (1d0-w)*xp(:,a) + w*xp(:,b)
     endif
     do iel=1,nxa_nel
        o = nxa_el(iel)%xoff
        k = nxa_el(iel)%ns
        w = sum(nxt_xprev(o+1:o+k,ic))
        if(w>0d0) nxt_xprev(o+1:o+k,ic) = nxt_xprev(o+1:o+k,ic)/w
     enddo
  enddo
  nxt_dtinv = 1d0/dt
  if(lprint) write(6,'(1x,a,a,a,es10.3,a,i5,a,es10.3,a,es10.3,a)') 'NLTE time-dependent: previous model ', &
       trim(fname),', dt =',dt,' s;',npc,' cells, exterior mass',mp(1),' g (this model',mext(1),' g)'
  deallocate(mp,xp)
contains
  subroutine nextrec
!-- next record that is not a comment
    do
       read(iu,'(a)',iostat=ios) line
       if(ios/=0) stop 'nxt_read_prev: file truncated'
       if(line(1:1)/='#' .and. len_trim(line)>0) exit
    enddo
  end subroutine nextrec
end subroutine nxt_read_prev
!
!
subroutine nxt_release(ncell,vol)
!-----------------------------------------------------------------------
! Power released by the change of the ionization and excitation energy of
! the NLTE elements since the previous model (time-dependent statistical
! equilibrium), per cell [erg/s]:
!   nxt_erel = V/dt sum_el n_el sum_s eps_s (x_prev,s - x_s)
! (negative: energy stored), and the cumulative activation weights of the
! states that empty, n_el eps_s (x_prev,s - x_s) > 0.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: vol(ncell)
  integer :: ic,iel,o,s
  real*8 :: e,c,w
  if(.not.allocated(nxt_erel)) allocate(nxt_erel(ncell),nxt_prel(ncell),nxt_rcum(nxa_nst,ncell))
  nxt_erel = 0d0
  nxt_prel = 0d0
  nxt_rcum = 0d0
  if(nxt_dtinv<=0d0 .or. .not.allocated(nxt_xprev)) return
  do ic=1,ncell
     e = 0d0
     c = 0d0
     do iel=1,nxa_nel
        o = nxa_el(iel)%xoff
        do s=1,nxa_el(iel)%ns
           w = nxt_nel(iel,ic)*nxa_el(iel)%eps(s)*(nxt_xprev(o+s,ic) - nxt_x(o+s,ic))
           e = e + w
           c = c + max(w,0d0)
           nxt_rcum(o+s,ic) = c
        enddo
     enddo
     nxt_erel(ic) = e*nxt_dtinv*vol(ic)
  enddo
end subroutine nxt_release
!
!
pure subroutine nxt_activate_rel(ic,state,kind,idx,wlc)
!-- packet of the energy released by the recombining gas (nxt_release):
!-- macro-atom activated in a state that empties, with its share
  integer,intent(in) :: ic
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: r
  integer :: g,iel
  if(nxt_rcum(nxa_nst,ic) <= 0d0) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
     return
  endif
  call rnd_r(r,state)
  call pick(nxt_rcum(:,ic),max(r,1d-300)*nxt_rcum(nxa_nst,ic),g)
  do iel=nxa_nel,1,-1
     if(g > nxa_el(iel)%xoff) exit
  enddo
  call nxt_activate(ic,iel,g-nxa_el(iel)%xoff,state,kind,idx,wlc)
end subroutine nxt_activate_rel
!
end module nltextrmod
! vim: fdm=marker
