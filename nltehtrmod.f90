!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nltehtrmod
!-----------------------------------------------------------------------
! NLTE hydrogen in the radiation transport: cell state, Monte Carlo
! estimators and Lucy's (2002, A&A 384, 725; 2003, A&A 403, 261) macro-atom.
!
! The H lines are not part of the group opacities.  A packet meets an H
! line where its comoving wavelength equals the line wavelength (found with
! the same exact crossing as group boundaries, so a non-monotonic flow can
! have several resonance points per line) and interacts there with
! probability 1 - exp(-tau), tau = K/|Q(mu)| (Sobolev, with the projected
! velocity gradient Q(mu) = mu^2 dv/dr + (1-mu^2) v/r).  Every crossing is
! tallied in the resonance estimator Jb of the incident intensity.
!
! An absorbed line packet activates the macro-atom in the line's upper
! level; H bound-free absorption at nu activates the continuum with
! probability nu_n/nu and the thermal pool otherwise.  The macro-atom jumps
! between levels with probabilities set by the statistical-equilibrium rates
! and de-activates by emitting a line, recombination-continuum or two-photon
! packet, or by collisional de-excitation into the thermal pool (k-packet).
! A k-packet is emitted thermally (non-H emissivity at the gas temperature)
! or excites/ionizes H, or is radiated as recombination cooling, in
! proportion to the cooling rates.  Packets keep their comoving energy
! throughout, so the scheme conserves energy exactly.
!
! The random walk through the macro-atom and the thermal pool is an
! absorbing Markov chain (transient states: k-packet, bound levels,
! continuum; absorbing: the emission channels).  Its absorption
! probabilities B = (I - Q)^-1 R are computed once per cell and iteration,
! so an activation is a single draw: identical statistics to the walk,
! without the very long walks of collision-dominated (near-LTE) gas.
!-----------------------------------------------------------------------
  use nltehmod
  use randommod
  implicit none
  private
!
  real*8,parameter :: pi = 3.141592653589793d0, cc = 2.99792458d10
  real*8,parameter :: hh = 6.6260693d-27, kb = 1.3806505d-16, arad = 7.5657d-15
!
  logical,public :: nht_on = .false.
  real*8,public :: nht_damp = 0.5d0       !weight of the new solution between iterations
  integer,public :: nht_njbmin = 30       !crossings needed to use the resonance estimator
  logical,public :: nht_have_mc = .false. !Monte Carlo estimators available
  real*8,public :: nht_tau0 = 0d0         !benchmark: pure-scattering H-alpha with this Sobolev depth
  real*8,public :: nht_ntemin = 0d0       !photons above this energy [erg] make fast photoelectrons (0: off)
!
  integer,public :: nht_ncell = 0, nht_ng = 0
  integer,public :: nht_nst = 0           !transient states: k-packet, levels 1..nlev, continuum
  integer,public :: nht_nout = 0          !outcomes: lines, bf to each level, two-photon, thermal
!-- all grid cells
  real*8,allocatable,public :: nht_x(:,:)     !(0:nlev,ncell) fractions of n_H, 0 = H II
  real*8,allocatable,public :: nht_nh(:)      !(ncell) H number density [cm^-3]
  real*8,allocatable,public :: nht_ne(:)      !(ncell) electron density [cm^-3]
  real*8,allocatable,public :: nht_kq(:,:)    !(nlin,ncell) Sobolev strength K [1/s]
  real*8,allocatable,public :: nht_beta(:,:)  !(nlin,ncell) escape probability
  real*8,allocatable,public :: nht_jg(:,:)    !(ng,ncell) MC mean intensity [erg/s/cm^2/Hz/sr]
  real*8,allocatable,public :: nht_jb(:,:)    !(nlin,ncell) resonance mean intensity
  logical,allocatable,public :: nht_ljb(:,:)  !(nlin,ncell) jb from enough crossings
  real*8,allocatable,public :: nht_jbsum(:,:) !(nlin,ncell) tally: sum e c/(nu |Q|)
  real*8,allocatable,public :: nht_jbcnt(:,:) !(nlin,ncell) tally: crossings
  real*8,allocatable,public :: nht_acum(:,:,:) !(nout,nst,ncell) cumulative outcome probabilities per start state
  real*4,allocatable,public :: grd_caphbf(:,:) !(ng,ncell) H bound-free part of grd_cap
!-- this rank's gas cells
  real*8,allocatable,public :: nht_xloc(:,:)  !(0:nlev,gas_ncell)
  real*8,allocatable,public :: nht_neloc(:),nht_nhloc(:) !(gas_ncell)
  real*8,allocatable,public :: nht_neothloc(:) !(gas_ncell) electrons per atom from other elements
  real*4,allocatable,public :: gas_caphbf(:,:) !(ng,gas_ncell)
!-- lines and groups
  integer,allocatable,public :: nht_lgrp(:)   !(nlin) group of the line
  integer,allocatable,public :: nht_isrt(:)   !(nlin) position of line k in the active list, 0 if inactive
!-- active lines: Sobolev depth above nht_taumin somewhere in the model; only
!-- these are resonances in the transport (thin lines use the group J)
  real*8,public :: nht_taumin = 1d-4
  integer,public :: nht_nact = 0
  integer,allocatable,public :: nht_act(:)    !(nlin) active line indices sorted by wavelength
  real*8,allocatable,public :: nht_wlact(:)   !(nlin) their wavelengths
  real*8,allocatable,public :: nht_nug(:),nht_dnug(:) !(ng) group centre frequency and width
!-- emission kinds returned by the macro-atom
  integer,parameter,public :: em_line = 1, em_bf = 2, em_thermal = 3, em_2ph = 4
!
  public :: nlteh_tr_init,nlteh_gas_solve,nlteh_store_jg,nlteh_tally_reset
  public :: nlteh_tally_finish,nlteh_tables,nlteh_neighbors,nlteh_activate
  public :: nlteh_kpacket,nlteh_bfabs,nlteh_write,nlteh_outcome
  public :: nlteh_nonthermal,nlteh_ntdep,nlteh_ntrates
!
contains
!
!
subroutine nlteh_tr_init(nlev,ncell,gncell,ng,wlg)
!-----------------------------------------------------------------------
  integer,intent(in) :: nlev,ncell,gncell,ng
  real*8,intent(in) :: wlg(ng+1)  !group boundaries [cm]
  integer :: k,ig,j
!
  call nlteh_init(nlev)
  nht_on = .true.
  nht_ncell = ncell
  nht_ng = ng
  nht_nst = nlev + 2
  nht_nout = nh_nlin + nlev + 2
  allocate(nht_x(0:nlev,ncell),nht_nh(ncell),nht_ne(ncell))
  allocate(nht_kq(nh_nlin,ncell),nht_beta(nh_nlin,ncell))
  allocate(nht_jg(ng,ncell),nht_jb(nh_nlin,ncell),nht_ljb(nh_nlin,ncell))
  allocate(nht_jbsum(nh_nlin,ncell),nht_jbcnt(nh_nlin,ncell))
  allocate(nht_acum(nht_nout,nht_nst,ncell))
  allocate(grd_caphbf(ng,ncell))
  allocate(nht_xloc(0:nlev,gncell),nht_neloc(gncell),nht_nhloc(gncell),nht_neothloc(gncell))
  allocate(gas_caphbf(ng,gncell))
  allocate(nht_lgrp(nh_nlin),nht_isrt(nh_nlin),nht_nug(ng),nht_dnug(ng))
  allocate(nht_act(nh_nlin),nht_wlact(nh_nlin))
  nht_x = 0d0; nht_nh = 0d0; nht_ne = 0d0; nht_kq = 0d0; nht_beta = 1d0
  nht_jg = 0d0; nht_jb = 0d0; nht_ljb = .false.; nht_jbsum = 0d0; nht_jbcnt = 0d0
  nht_acum = 1d0; grd_caphbf = 0.; gas_caphbf = 0.
  nht_xloc = 0d0; nht_neloc = 0d0; nht_nhloc = 0d0; nht_neothloc = 0d0
!
  do ig=1,ng
     nht_nug(ig) = 0.5d0*(cc/wlg(ig) + cc/wlg(ig+1))
     nht_dnug(ig) = cc/wlg(ig) - cc/wlg(ig+1)
  enddo
  do k=1,nh_nlin
     nht_lgrp(k) = 0
     do ig=1,ng
        if(nh_wl(k) >= wlg(ig) .and. nh_wl(k) < wlg(ig+1)) then
           nht_lgrp(k) = ig
           exit
        endif
     enddo
  enddo
!-- all lines active until the first tables are built
  nht_nact = nh_nlin
  do j=1,nh_nlin
     nht_act(j) = nh_srt(j)
     nht_wlact(j) = nh_wlsrt(j)
     nht_isrt(nh_srt(j)) = j
  enddo
end subroutine nlteh_tr_init
!
!
subroutine nlteh_gas_solve(it,gncell,temp,trad,nh,neoth,qb,qv,icell1,edep)
!-----------------------------------------------------------------------
! Statistical equilibrium of H in this rank's gas cells.  The first
! iteration uses J = B(trad); later ones the Monte Carlo estimators of the
! previous iteration, and the new populations are mixed geometrically with
! the old ones (weight nht_damp).  edep: power absorbed from photons above
! nht_ntemin [erg/s/cm^3], deposited by fast electrons (non-thermal ionization
! and excitation of H).
!-----------------------------------------------------------------------
  integer,intent(in) :: it,gncell,icell1
  real*8,intent(in) :: temp(gncell),trad(gncell),nh(gncell),neoth(gncell)
  real*8,intent(in) :: edep(gncell)
  real*8,intent(in) :: qb(gncell),qv(gncell)
  integer :: i,l,ig,ierr,nfail
  real*8 :: x(0:nh_nlev),xold(0:nh_nlev),ne,beta(nh_nlin),jg(nht_ng),jb(nh_nlin)
  real*8 :: nonth,nonexc
  logical :: ljb(nh_nlin)
!
  nfail = 0
  do i=1,gncell
     l = icell1 + i - 1
     if(nh(i) <= 1d-10*max(neoth(i),1d0)) then
!-- no hydrogen
        nht_xloc(:,i) = 0d0
        nht_xloc(0,i) = 1d0
        nht_nhloc(i) = nh(i)
        nht_neloc(i) = neoth(i) + nh(i)
        cycle
     endif
     if(it==1 .or. .not.nht_have_mc) then
        do ig=1,nht_ng
           jg(ig) = bplanck(nht_nug(ig),trad(i))
        enddo
        jb = 0d0
        ljb = .false.
        call nlteh_lte(temp(i),max(neoth(i),1d-3*nh(i)),x)
        x(0) = max(x(0),1d-6)
        x = x/sum(x)
     else
        jg = nht_jg(:,l)
        jb = nht_jb(:,l)
        ljb = nht_ljb(:,l)
        x = nht_xloc(:,i)
     endif
     xold = x
     call nlteh_ntrates(nh(i),x,edep(i),nonth,nonexc)
     call nlteh_solve_cell(temp(i),trad(i),nh(i),neoth(i),qb(i),qv(i),nht_ng, &
          nht_nug,nht_dnug,jg,jb,ljb,nonth,x,ne,beta,ierr,nonexc)
     if(ierr/=0) then
        nfail = nfail + 1
        x = xold
     elseif(it>1 .and. nht_have_mc .and. nht_damp<1d0) then
        x = xold**(1d0-nht_damp)*x**nht_damp
        x = x/sum(x)
     endif
     nht_xloc(:,i) = x
     nht_nhloc(i) = nh(i)
     nht_neloc(i) = neoth(i) + x(0)*nh(i)
  enddo
  if(nfail>0) write(0,*) 'nlteh_gas_solve: no solution in cells:',nfail
end subroutine nlteh_gas_solve
!
!
subroutine nlteh_store_jg(ncell,ng,jrad,vol)
!-- MC path-length tally (sum of e*d per group) to mean intensity per Hz per sr
  integer,intent(in) :: ncell,ng
  real*8,intent(in) :: jrad(ng,ncell),vol(ncell)
  integer :: i
  do i=1,ncell
     nht_jg(:,i) = jrad(:,i)/(4d0*pi*vol(i)*nht_dnug)
  enddo
end subroutine nlteh_store_jg
!
!
subroutine nlteh_tally_reset
  nht_jbsum = 0d0
  nht_jbcnt = 0d0
end subroutine nlteh_tally_reset
!
!
subroutine nlteh_tally_finish(ncell,vol)
!-- resonance tally (sum of e c/(nu |Q|)) to mean intensity incident on the line
  integer,intent(in) :: ncell
  real*8,intent(in) :: vol(ncell)
  integer :: i
  do i=1,ncell
     nht_jb(:,i) = nht_jbsum(:,i)/(4d0*pi*vol(i))
     nht_ljb(:,i) = nht_jbcnt(:,i) >= nht_njbmin
  enddo
  nht_have_mc = .true.
end subroutine nlteh_tally_finish
!
!
subroutine nlteh_tables(ncell,tempinv,radtemp,capemitgrey,xarr,vxarr)
!-----------------------------------------------------------------------
! Sobolev strengths, escape probabilities, macro-atom and k-packet tables
! for all cells, from the gathered populations.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: tempinv(ncell),radtemp(ncell),capemitgrey(ncell)
  real*8,intent(in) :: xarr(ncell+1),vxarr(ncell+1)
  integer :: ic,n,l,u,k,nl,j,js,nst,nout,nlin,info,nfail
  real*8 :: t,nh,ne,x(0:nh_nlev),qb,qv,rr,eps(0:nh_nlev),s
  real*8 :: q(nht_nst,nht_nst),rm(nht_nst,nht_nout)
  integer :: ipiv(nht_nst)
  real*8 :: gam(nh_nlev),alp(nh_nlev),alps(nh_nlev),cup(nh_nlev),c3b(nh_nlev)
  real*8 :: qup(nh_nlev,nh_nlev),qdn(nh_nlev,nh_nlev),jl(nh_nlin),jb(nh_nlin)
  logical :: ljb(nh_nlin)
  real*8 :: jg(nht_ng)
  integer :: ig
!
  nl = nh_nlev
  nlin = nh_nlin
  nst = nht_nst
  nout = nht_nout
  nfail = 0
  eps(0) = nh_chi
  eps(1:nl) = nh_eexc
  do ic=1,ncell
     t = 1d0/tempinv(ic)
     nh = nht_nh(ic)
     ne = nht_ne(ic)
     x = nht_x(:,ic)
     qb = (vxarr(ic+1) - vxarr(ic))/(xarr(ic+1) - xarr(ic))
     qv = (vxarr(ic+1) + vxarr(ic))/(xarr(ic+1) + xarr(ic))
     if(nht_have_mc) then
        jg = nht_jg(:,ic)
        jb = nht_jb(:,ic)
        ljb = nht_ljb(:,ic)
     else
        do ig=1,nht_ng
           jg(ig) = bplanck(nht_nug(ig),radtemp(ic))
        enddo
        jb = 0d0
        ljb = .false.
     endif
     call nlteh_cell_rates(t,radtemp(ic),nht_ng,nht_nug,nht_dnug,jg,jb,ljb,0d0, &
          gam,alp,alps,cup,c3b,qup,qdn,jl)
     do k=1,nh_nlin
        nht_kq(k,ic) = nlteh_kq(k,nh,x)
        nht_beta(k,ic) = nlteh_beta(nht_kq(k,ic),qb,qv)
     enddo
     if(nht_tau0>0d0) then
!-- benchmark: only H-alpha, tau = tau0 for a gradient v/r
        nht_kq(:,ic) = 0d0
        nht_kq(nh_lin(2,3),ic) = nht_tau0*qv
     endif
!
!-- transition probabilities of the macro-atom/thermal-pool chain
!-- state 1: k-packet, 1+n: level n, nl+2: continuum
!-- outcome k: line k, nlin+l: bf to l, nlin+nl+1: two-photon, nlin+nl+2: thermal
     q = 0d0
     rm = 0d0
!-- bound levels
     do n=1,nl
        js = 1 + n
        do l=1,n-1
           k = nh_lin(l,n)
           rr = nh_a(k)*nht_beta(k,ic)
           rm(js,k) = rm(js,k) + rr*(eps(n) - eps(l))           !line emission
           q(js,1+l) = q(js,1+l) + (rr + ne*qdn(n,l))*eps(l)     !internal jump down
           q(js,1) = q(js,1) + ne*qdn(n,l)*(eps(n) - eps(l))     !collisional -> k
        enddo
        if(n==2) rm(js,nlin+nl+1) = 0.25d0*a2s*(eps(2) - eps(1)) !two-photon
        do u=n+1,nl
           k = nh_lin(n,u)
           q(js,1+u) = q(js,1+u) + (nht_beta(k,ic)*nh_blu(k)*jl(k) + ne*qup(n,u))*eps(n)
        enddo
        q(js,nl+2) = (gam(n) + ne*cup(n))*eps(n)               !to the continuum
     enddo
!-- continuum
     js = nl + 2
     do l=1,nl
        rm(js,nlin+l) = ne*alp(l)*(eps(0) - eps(l))             !recombination emission
        q(js,1+l) = q(js,1+l) + ne*(alp(l) + alps(l) + ne*c3b(l))*eps(l)
        q(js,1) = q(js,1) + ne*ne*c3b(l)*(eps(0) - eps(l))      !three-body -> k
     enddo
!-- k-packet: cooling rates per volume
     rm(1,nlin+nl+2) = arad*cc*t**4*capemitgrey(ic)             !thermal, non-H
     do k=1,nlin
        l = nh_ll(k); u = nh_lu(k)
        q(1,1+u) = q(1,1+u) + nh*x(l)*ne*qup(l,u)*(eps(u) - eps(l))
     enddo
     do n=1,nl
        q(1,nl+2) = q(1,nl+2) + nh*x(n)*ne*cup(n)*nh_eion(n)
        rm(1,nlin+n) = rm(1,nlin+n) + nh*x(0)*ne*alp(n)*kb*t
     enddo
!-- normalize rows; a state without any channel is emitted thermally
     do js=1,nst
        s = sum(q(js,:)) + sum(rm(js,:))
        if(s <= 0d0) then
           q(js,:) = 0d0
           rm(js,:) = 0d0
           rm(js,nlin+nl+2) = 1d0
        else
           q(js,:) = q(js,:)/s
           rm(js,:) = rm(js,:)/s
        endif
     enddo
!-- absorption probabilities: (I - Q) B = R
     do js=1,nst
        q(js,:) = -q(js,:)
        q(js,js) = q(js,js) + 1d0
     enddo
     call dgesv(nst,nout,q,nst,ipiv,rm,nst,info)
     if(info/=0) then
        rm = 0d0
        rm(:,nlin+nl+2) = 1d0
        nfail = nfail + 1
     endif
     do js=1,nst
        rm(js,:) = max(rm(js,:),0d0)
        s = sum(rm(js,:))
        if(s<=0d0) then
           rm(js,:) = 0d0
           rm(js,nlin+nl+2) = 1d0
           s = 1d0
        endif
        nht_acum(1,js,ic) = rm(js,1)/s
        do k=2,nout
           nht_acum(k,js,ic) = nht_acum(k-1,js,ic) + rm(js,k)/s
        enddo
        nht_acum(nout,js,ic) = 1d0
     enddo
  enddo
  if(nfail>0) write(0,*) 'nlteh_tables: singular macro-atom chain in cells:',nfail
!
!-- active lines: tau > nht_taumin for the smallest projected gradient
  nht_nact = 0
  nht_isrt = 0
  do j=1,nh_nlin
     k = nh_srt(j)
     s = 0d0
     do ic=1,ncell
        qb = (vxarr(ic+1) - vxarr(ic))/(xarr(ic+1) - xarr(ic))
        qv = (vxarr(ic+1) + vxarr(ic))/(xarr(ic+1) + xarr(ic))
        s = max(s, nht_kq(k,ic)/max(min(abs(qb),abs(qv)),1d-3*max(abs(qb),abs(qv))))
     enddo
     if(s > nht_taumin) then
        nht_nact = nht_nact + 1
        nht_act(nht_nact) = k
        nht_wlact(nht_nact) = nh_wl(k)
        nht_isrt(k) = nht_nact
     endif
  enddo
end subroutine nlteh_tables
!
!
pure subroutine nlteh_neighbors(wlc,iline,ilside,jlo,jhi)
!-----------------------------------------------------------------------
! Indices in the active line list of the H lines just below (jlo) and above
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
!-- largest index with nh_wlsrt <= wlc
     a = 0
     b = nht_nact + 1
     do while(b-a>1)
        m = (a+b)/2
        if(nht_wlact(m) <= wlc) then
           a = m
        else
           b = m
        endif
     enddo
     jlo = a
     jhi = a + 1
  endif
  if(jhi>nht_nact) jhi = 0
  if(jlo<1) jlo = 0
end subroutine nlteh_neighbors
!
!
pure subroutine bf_freq(n,t,state,nu)
!-- frequency of a recombination photon to level n: j_nu ~ g_bf e^{-h nu/kT}
  integer,intent(in) :: n
  real*8,intent(in) :: t
  type(rnd_t),intent(inout) :: state
  real*8,intent(out) :: nu
  real*8 :: r,gmax,g
  integer :: it
  gmax = max(nlteh_gbf(n,1d0), nlteh_gbf(n,1d30), nlteh_gbf(n,2d0), nlteh_gbf(n,4d0))
  do it=1,1000
     call rnd_r(r,state)
     nu = nh_nu0(n) - log(max(r,1d-300))*kb*t/hh
     g = nlteh_gbf(n,nu/nh_nu0(n))
     call rnd_r(r,state)
     if(r*gmax <= g) exit
  enddo
end subroutine bf_freq
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
pure subroutine nlteh_outcome(ic,js,state,kind,idx,wlc)
!-----------------------------------------------------------------------
! Draw the emission of a packet entering the macro-atom/thermal-pool chain
! of cell ic in state js (1: k-packet, 1+n: level n, nlev+2: continuum).
! kind: em_line (idx = line, wlc = its wavelength), em_bf (idx = level, wlc
! sampled from its recombination spectrum), em_2ph (wlc sampled) or
! em_thermal (the caller samples the thermal emissivity).
!-----------------------------------------------------------------------
  use gridmod, only:grd_tempinv
  integer,intent(in) :: ic,js
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: r,nu
  integer :: e
  call rnd_r(r,state)
  call pick(nht_acum(:,js,ic),max(r,1d-300),e)
  idx = 0
  wlc = 0d0
  if(e<=nh_nlin) then
     kind = em_line
     idx = e
     wlc = nh_wl(e)
  elseif(e<=nh_nlin+nh_nlev) then
     kind = em_bf
     idx = e - nh_nlin
     call bf_freq(idx,1d0/grd_tempinv(ic),state,nu)
     wlc = cc/nu
  elseif(e==nh_nlin+nh_nlev+1) then
     kind = em_2ph
     call rnd_r(r,state)
     nu = nlteh_twophoton_y(r)*nh_nu(nh_lin(1,2))
     wlc = cc/nu
  else
     kind = em_thermal
  endif
end subroutine nlteh_outcome
!
!
pure subroutine nlteh_activate(ic,ilev,state,kind,idx,wlc)
!-- macro-atom activated in level ilev (0 = continuum)
  integer,intent(in) :: ic,ilev
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  if(ilev==0) then
     call nlteh_outcome(ic,nh_nlev+2,state,kind,idx,wlc)
  else
     call nlteh_outcome(ic,1+ilev,state,kind,idx,wlc)
  endif
end subroutine nlteh_activate
!
!
pure subroutine nlteh_nonthermal(ic,state,kind,idx,wlc)
!-- energy of a fast photoelectron in cell ic: heat (k-packet), H ionization
!-- (macro-atom in the continuum) or H excitation (macro-atom in n = 2) with
!-- the Shull & van Steenberg (1985) shares
  integer,intent(in) :: ic
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: fh,fi,fe,r
  call nlteh_svs(nht_x(0,ic),fh,fi,fe)
  call rnd_r(r,state)
  if(r < fh) then
     call nlteh_kpacket(ic,state,kind,idx,wlc)
  elseif(r < fh + fi) then
     call nlteh_activate(ic,0,state,kind,idx,wlc)
  else
     call nlteh_activate(ic,2,state,kind,idx,wlc)
  endif
end subroutine nlteh_nonthermal
!
!
pure function nlteh_ntdep(ng,nug,dnug,jg,cap) result(edep)
!-- power absorbed from photons above nht_ntemin [erg/s/cm^3]
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng)
  real*4,intent(in) :: cap(ng)
  real*8 :: edep
  integer :: ig
  edep = 0d0
  if(nht_ntemin <= 0d0) return
  do ig=1,ng
     if(hh*nug(ig) > nht_ntemin) edep = edep + cap(ig)*jg(ig)*dnug(ig)
  enddo
  edep = 4d0*pi*edep
end function nlteh_ntdep
!
!
pure subroutine nlteh_ntrates(nh,x,edep,nonth,nonexc)
!-- non-thermal ionization and excitation (1 -> 2) rates per ground-state H
!-- atom from the deposited power: chi_H and 3/4 chi_H per event
  real*8,intent(in) :: nh,x(0:nh_nlev),edep
  real*8,intent(out) :: nonth,nonexc
  real*8 :: fh,fi,fe,n1
  nonth = 0d0
  nonexc = 0d0
  n1 = nh*x(1)
  if(edep <= 0d0 .or. n1 <= 1d-30) return
  call nlteh_svs(x(0),fh,fi,fe)
  nonth = fi*edep/(nh_chi*n1)
  nonexc = fe*edep/(0.75d0*nh_chi*n1)
end subroutine nlteh_ntrates
!
!
pure subroutine nlteh_kpacket(ic,state,kind,idx,wlc)
!-- a k-packet (thermal energy) in cell ic
  integer,intent(in) :: ic
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  call nlteh_outcome(ic,1,state,kind,idx,wlc)
end subroutine nlteh_kpacket
!
!
pure subroutine nlteh_bfabs(ic,nu,state,kind,idx,wlc)
!-----------------------------------------------------------------------
! H bound-free absorption of a packet of comoving frequency nu in cell ic:
! choose the level by n_n sigma_n(nu), then activate the continuum with
! probability nu_n/nu (the ionization energy) or the thermal pool (the
! photo-electron's kinetic energy).
!-----------------------------------------------------------------------
  integer,intent(in) :: ic
  real*8,intent(in) :: nu
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: w(nh_nlev),s,r
  integer :: n
  do n=1,nh_nlev
     w(n) = nht_x(n,ic)*nlteh_sigma(n,nu)
  enddo
  s = sum(w)
  if(s<=0d0) then
     call nlteh_kpacket(ic,state,kind,idx,wlc)
     return
  endif
  call rnd_r(r,state)
  r = r*s
  s = 0d0
  do n=1,nh_nlev
     s = s + w(n)
     if(s >= r) exit
  enddo
  n = min(n,nh_nlev)
  call rnd_r(r,state)
  if(r < nh_nu0(n)/nu) then
     call nlteh_activate(ic,0,state,kind,idx,wlc)
  else
     call nlteh_kpacket(ic,state,kind,idx,wlc)
  endif
end subroutine nlteh_bfabs
!
!
subroutine nlteh_write(ncell,rcell,vcell,tempinv,radtemp)
!-- per-cell NLTE H diagnostics
  integer,intent(in) :: ncell
  real*8,intent(in) :: rcell(ncell),vcell(ncell),tempinv(ncell),radtemp(ncell)
  integer :: ic,kha,khb,kla,n
  real*8 :: t,b(0:3),xl(0:nh_nlev),qv
  kla = nh_lin(1,2); kha = nh_lin(2,3); khb = nh_lin(2,4)
  open(4,file='output.nlteh',status='replace')
  write(4,'(a)') '# NLTE hydrogen per cell: b_n = departure coefficients relative to the'// &
       ' actual ion density (b = 1: LTE); tau = K/(v/r) (homologous equivalent)'
  write(4,'(a1,a11,19a12)') '#','r[cm]','v[cm/s]','T[K]','Trad[K]','nH','ne','x_HII','x_1','x_2','x_3', &
       'b_1','b_2','b_3','tau_Lya','tau_Ha','tau_Hb','beta_Ha','Jb_Ha/B','ljb_Ha'
  do ic=1,ncell
     t = 1d0/tempinv(ic)
     call nlteh_lte(t,nht_ne(ic),xl)
!-- departure coefficient: (n_n/n_+)/(n_n/n_+)_LTE at the actual ne
     do n=1,3
        b(n) = (nht_x(n,ic)/max(nht_x(0,ic),1d-300))/(xl(n)/max(xl(0),1d-300))
     enddo
     qv = max(vcell(ic)/rcell(ic),1d-30)
     write(4,'(1p,19e12.4)') rcell(ic),vcell(ic),t,radtemp(ic),nht_nh(ic),nht_ne(ic), &
          nht_x(0,ic),nht_x(1,ic),nht_x(2,ic),nht_x(3,ic),b(1),b(2),b(3), &
          nht_kq(kla,ic)/qv,nht_kq(kha,ic)/qv,nht_kq(khb,ic)/qv,nht_beta(kha,ic), &
          nht_jb(kha,ic)/max(bplanck(nh_nu(kha),t),1d-300),merge(1d0,0d0,nht_ljb(kha,ic))
  enddo
  close(4)
end subroutine nlteh_write
!
end module nltehtrmod
! vim: fdm=marker
