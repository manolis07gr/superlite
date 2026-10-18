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
! otherwise, on IMC steps and in DDMC groups alike (diffusion11).  The
! macro-atom jumps between the states of the element (all levels of its
! NLTE ions and the top ion) with probabilities set by the
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
  use nltedismod
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
  real*8,allocatable,public :: nxt_bes(:,:)   !(nlin,ncell) electron scattering of trapped photons (nx_bes)
  real*8,allocatable,public :: nxt_jg(:,:)    !(ng,ncell) MC mean intensity [erg/s/cm^2/Hz/sr]
  real*8,allocatable,public :: nxt_jb(:,:)    !(nlin,ncell) resonance mean intensity
  logical,allocatable,public :: nxt_ljb(:,:)  !(nlin,ncell) jb set (MC; nxt_tally_finish)
  real*8,allocatable,public :: nxt_jbsum(:,:) !(nlin,ncell) tally: sum e c/(nu |Q|)
  real*8,allocatable,public :: nxt_jbcnt(:,:) !(nlin,ncell) tally: crossings
!-- incident intensity of the Sobolev lines from the flights in a window next
!-- to each line (nxt_window_tally; 0: off): its width [cm/s], the tallies
!-- (sum of the contributions and of their squares), and the widths to the
!-- blue and red of each active line (up to the next active line)
  real*8,public :: nxt_jbwin = 0d0
  real*8,allocatable,public :: nxt_jwsum(:,:),nxt_jwsq(:,:) !(nlin,ncell)
  real*8,allocatable :: nxt_wwb(:),nxt_wwr(:) !(nlin) by active index
!-- (diagnostics, output.nlteh: Halpha emitted by the macro-atom and absorbed
!-- at its resonances in each cell on the last Monte Carlo step [erg/s])
  real*8,allocatable,public :: nxt_haem(:),nxt_haab(:) !(ncell)
!-- (largest Sobolev depth of each active line over the cells; the windows end
!-- at the next line with at least nxt_tauwin: the weaker ones absorb little)
  real*8,allocatable :: nxt_tmact(:)          !(nlin) by active index
  real*8,parameter :: nxt_tauwin = 0.1d0
!-- radiation of the optically thick groups (nxt_ldb): tallies of the packets
!-- on IMC steps outside the cell of their last thermal or bound-free
!-- emission (radiation from elsewhere), and of the DDMC packets in the
!-- cell's group lump; groups with absorption depth across the cell above
!-- nxt_ldbtau (0: none) are in local detailed balance
  real*8,public :: nxt_ldbtau = 0d0
  logical,public :: nxt_alion = .true.        !accelerated lambda iteration in the other groups
!-- (with DDMC, groups that thick are DDMC groups: their radiation is the
!-- diffusion limit, J = S to second order in 1/tau, with no Monte Carlo
!-- tally of the radiation from elsewhere: Lambda* = 1)
  logical,public :: nxt_ldbone = .false.
  real*8,allocatable,public :: nxt_jxsum(:,:) !(ng,ncell)
!-- tally of the photons of NLTE lines on IMC steps in the line's group with
!-- their comoving frequency history (packet%ilem): not incident on the line
  real*8,allocatable,public :: nxt_jlosum(:,:) !(ng,ncell)
  real*8,allocatable,public :: nxt_jlsum(:)   !(ncell)
!-- (diagnostics, output.nltex_ldb: per cell, over the groups in local
!-- detailed balance at 228-504 and 504-912 A, 4 pi Int dnu of Lambda* S,
!-- J_ext and the Monte Carlo estimate)
  real*8,allocatable :: nxt_ldbd(:,:,:)       !(3,2,ncell)
  real*4,allocatable,public :: grd_capnbf(:,:) !(ng,ncell) NLTE bound-free part of grd_cap
!-- NLTE recombination emission per group, in the units of the Planck group
!-- integrals times opacity (int eta dnu/(sigma T^4/pi)): DDMC source spectrum
  real*4,allocatable,public :: grd_emnbf(:,:) !(ng,ncell)
!-- this rank's gas cells
  real*8,allocatable,public :: nxt_xloc(:,:)  !(nst,gas_ncell)
  real*8,allocatable,public :: nxt_neloc(:)   !(gas_ncell)
  real*8,allocatable,public :: nxt_nelloc(:,:) !(nel,gas_ncell)
  real*8,allocatable,public :: nxt_neothloc(:) !(gas_ncell) electrons per atom from the other elements
  real*4,allocatable,public :: gas_capnbf(:,:) !(ng,gas_ncell)
  real*4,allocatable,public :: gas_emnbf(:,:) !(ng,gas_ncell)
!-- lines and groups
  integer,allocatable,public :: nxt_lgrp(:)   !(nlin) group of the line
  integer,allocatable,public :: nxt_isrt(:)   !(nlin) position of line k in the active list, 0 if inactive
  integer,public :: nxt_nact = 0
  integer,allocatable,public :: nxt_act(:)    !(nlin) active lines sorted by wavelength
  real*8,allocatable,public :: nxt_wlact(:)   !(nlin) their wavelengths
!-- lookup of the active lines by ln(wavelength) for nxt_neighbors:
!-- nxt_wcnt(k) lines at or below exp(nxt_lw0 + k nxt_dlw)
  integer,parameter :: nxt_nwt = 65536
  integer,allocatable :: nxt_wcnt(:)          !(0:nwt)
  real*8 :: nxt_lw0 = 0d0, nxt_dlw = 1d0
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
!-- (in_nltex_bbext) the line-list lines of NLTE ions beyond the model atoms,
!-- and of the top ions, in the Sobolev expansion opacity
  logical,public :: nxt_bbext = .true.
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
!-- warm start (in_nltex_init): populations, T_e and n_e of a saved state
!-- (output.nltex_state) at this model's cells; the first iteration keeps
!-- them instead of solving with J = B(T_R), and makes its tables with the
!-- saved Monte Carlo estimators (output.nltex_jmc) when there are any
  logical,public :: nxt_warm = .false.
!-- (in_nltex_hold) the saved populations in every iteration (diagnostic)
  logical,public :: nxt_hold = .false.
  real*8,allocatable,public :: nxt_xinit(:,:) !(nst,ncell)
  real*8,allocatable,public :: nxt_tinit(:),nxt_neinit(:) !(ncell)
!-- microclumping: volume filling factor per cell (not allocated: smooth).
!-- The populations are those of the clumps; the Sobolev strengths K of
!-- the transport and of the escape probabilities are volume means, f K
  real*8,allocatable,public :: nxt_fcl(:)    !(ncell)
!-- bound-free absorption: per block of groups and cell the nxt_bfk levels
!-- of largest bound n_l w_l max_{g in block} sigma_l(g) of their weight
!-- (decreasing), their weight at each group, and the sum of the bounds of
!-- all other levels (nxt_bfabs)
  integer,parameter :: nxt_nbb = 256
  integer,public :: nxt_bfk = 256            !levels per list
  integer :: nxt_bsz = 0                     !groups per block
  integer :: nxt_nwm = 0                     !64-bit words of a level mask
  real*4,allocatable :: nxt_smaxb(:,:)      !(nlv,nbb) largest sigma of the level in the block
!-- level dissolution (nltedismod): per cell its parameters, and per level
!-- the lowest frequency of its dissolved continuum and the share of its
!-- radiative recombination that goes into it (nxt_tables)
  real*8,allocatable :: nxt_disp(:,:)       !(3,ncell)
  real*8,allocatable :: nxt_dnuc(:,:)       !(nlv,ncell)
  real*4,allocatable :: nxt_fdis(:,:)       !(nlv,ncell)
  integer,allocatable :: nxt_bfl(:,:,:)     !(bfk,nbb,ncell) the levels
  integer,allocatable :: nxt_bfn(:,:)       !(nbb,ncell) their number
  integer*8,allocatable :: nxt_bfm(:,:,:)   !(nwm,nbb,ncell) mask of the listed levels
  real*8,allocatable :: nxt_bfr(:,:)        !(nbb,ncell) bound of the others
  real*8,allocatable :: nxt_bfs(:,:)        !(ng,ncell) weight of the listed levels
!-- (groups with bound-free absorption, sum over them of the share of
!-- absorptions that walk the other levels, groups where it exceeds 0.1)
  real*8,public :: nxt_bfstat(3) = 0d0
!-- convergence of the populations: median and largest over this rank's
!-- cells of the cell's largest relative change of a state above 1e-4
!-- (after the damping; nxt_gas_solve)
  real*8,public :: nxt_dxstat(2) = 0d0
!-- electron scattering of line photons trapped in their Sobolev resonance
!-- zones (in_nltex_vturb >= 0): microturbulent velocity [cm/s] at v = 0
!-- and at the outer edge nxt_vmax (linear in v) for the Doppler widths.
!-- The probability per emission (nxt_bes) gives the scattered emission of
!-- the formal integral, n_u A bes h nu (formalmod).  With nxt_esrates also
!-- the rates and the macro-atom have beta+bes, and a line de-activation
!-- emits the photon scattered by a thermal electron with probability
!-- bes/(beta+bes) (em_lines, nxemitmod)
  logical,public :: nxt_estrap = .false.
  real*8,public :: nxt_vturb(2) = 0d0, nxt_vmax = 0d0
!-- (and in the statistical equilibrium, the macro-atom and the transport;
!-- off: only in the formal integral, as CMFGEN adds the incoherent
!-- scattering to the observer spectrum of iterations with coherent
!-- scattering)
  logical,public :: nxt_esrates = .false.
  logical,public :: nxt_esdamp = .true.   !Lorentz-wing escape factor (nx_wq)
!-- (with nxt_esrates, the photons that electrons scatter into the core of a
!-- Sobolev line on the side they drift away to are absorbed in the
!-- transport, nxt_esfar; the rates take that absorption from its tally,
!-- the comoving energy E per line and cell, which nxt_tally_finish turns
!-- into the beta-weighted incident intensity E c/(4 pi V nu K) with the K
!-- of the transport: in the rates it absorbs E K/K_transport)
  real*8,allocatable,public :: nxt_efsum(:,:) !(nlin,ncell)
!-- Doppler-profile lines in the Monte Carlo transport of slow cells
!-- (in_nltex_dopmc > 0): where the velocity changes across a cell by less
!-- than nxt_dopmc Doppler widths, a Sobolev resonance would lie further
!-- away than the cell; there a packet sees the active lines within
!-- nxt_dopx Doppler widths of its comoving frequency as static profiles,
!-- K/(sqrt(pi) v_D) exp(-x^2) per unit length (the overlapping line
!-- forest of a slow inner wind), except the line it was emitted in, whose
!-- trapping stays with the Sobolev escape of the macro-atom; their
!-- incident intensity is the path-length estimator over a window of
!-- nxt_dopw Doppler widths outside the profile, on the side the comoving
!-- frequency drifts in from (nxt_static_tally).  Only the lines of the
!-- elements from Z = nxt_dopz (nxt_lstat): the H and He lines keep their
!-- Sobolev resonances there.  (In a dense wind an electron scatters a
!-- photon about as often as it drifts across a Doppler width, and a static
!-- profile absorbs every photon scattered into the core of a thick line,
!-- while the macro-atom's Sobolev escape treats each absorption as a
!-- whole trapping episode without escape by electron scattering: the
!-- H/He benchmark lost a fifth of Halpha and two thirds of Hbeta.)
  real*8,public :: nxt_dopmc = 0d0
  real*8,parameter,public :: nxt_dopx = 4d0, nxt_dopw = 2d0
  integer,public :: nxt_dopz = 3
  logical,allocatable,public :: nxt_lstat(:)  !(nlin) static profile in slow cells
  logical,public :: nxt_anystat = .false.     !an active line with a static profile
!-- (nxt_dopes: the static lines of the slow cells escape their trapping
!-- also by the electron scattering of the trapped photons, and absorb the
!-- photons scattered into their cores, in the rates and the macro-atom,
!-- as with nxt_esrates for all lines: beta + bes, (beta Jb + bes Jes)/
!-- (beta + bes), and the scattered emission em_lines.  The transport
!-- scatters photons into the static cores explicitly, so without it the
!-- macro-atom converts every one of them as a whole trapping episode.)
  logical,public :: nxt_dopes = .false.
!-- nearest active line at or below (nxt_sjlo) and at or above (nxt_sjhi)
!-- an active index that keeps its Sobolev resonance in slow cells (0: none)
  integer,allocatable :: nxt_sjlo(:),nxt_sjhi(:) !(0:nlin+1)
  logical,allocatable,public :: nxt_slow(:)   !(ncell)
  real*8,allocatable,public :: nxt_vdl(:,:)   !(nel,ncell) Doppler widths [cm/s]
  real*8,allocatable,public :: nxt_vdlmax(:),nxt_vdlmin(:) !(ncell)
!-- emission kinds returned by the macro-atom (em_lines: a line photon
!-- scattered by an electron in the resonance zone, nxemitmod)
  integer,parameter,public :: em_line = 1, em_bf = 2, em_thermal = 3, em_2ph = 4, em_lines = 5
!
  public :: nxt_init,nxt_gas_solve,nxt_store_jg,nxt_tally_reset,nxt_tally_finish
  public :: nxt_tables,nxt_neighbors,nxt_activate,nxt_activate_line,nxt_kpacket
  public :: nxt_bfabs,nxt_bfpick,nxt_nonthermal,nxt_ntdep,nxt_ntrates,nxt_write,nxt_bffreq
  public :: nxt_shares,nxt_dealloc,nxt_emcode,nxt_emloc
  public :: nxt_write_state,nxt_read_prev,nxt_read_init,nxt_release,nxt_activate_rel
  public :: nxt_ali_cell,nxt_ali_off,nxt_ldb,nxt_lsldb,nxt_vtb,nxt_vturbv
  public :: nxt_static_line,nxt_static_tally,nxt_sobolev_neighbors,nxt_act_lookup,nxt_window_tally,nxt_esfar
  public :: nxt_useb,nxt_es_cell
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
  allocate(nxt_kq(nxa_nlin,ncell),nxt_beta(nxa_nlin,ncell),nxt_bes(nxa_nlin,ncell))
  allocate(nxt_jg(ng,ncell),nxt_jb(nxa_nlin,ncell),nxt_ljb(nxa_nlin,ncell))
  allocate(nxt_jbsum(nxa_nlin,ncell),nxt_jbcnt(nxa_nlin,ncell))
  allocate(nxt_jwsum(nxa_nlin,ncell),nxt_jwsq(nxa_nlin,ncell),nxt_haem(ncell),nxt_haab(ncell))
  nxt_haem = 0d0; nxt_haab = 0d0
  allocate(nxt_efsum(nxa_nlin,ncell))
  nxt_efsum = 0d0
  allocate(nxt_jxsum(ng,ncell),nxt_jlsum(ncell),nxt_jlosum(ng,ncell))
  nxt_jxsum = 0d0; nxt_jlsum = 0d0; nxt_jlosum = 0d0; nxt_jwsum = 0d0; nxt_jwsq = 0d0
  allocate(grd_capnbf(ng,ncell),grd_emnbf(ng,ncell))
  allocate(nxt_xloc(nst,gncell),nxt_neloc(gncell),nxt_nelloc(nxa_nel,gncell),nxt_neothloc(gncell))
  allocate(gas_capnbf(ng,gncell),gas_emnbf(ng,gncell))
  allocate(nxt_lgrp(nxa_nlin),nxt_isrt(nxa_nlin),nxt_nug(ng),nxt_dnug(ng))
  allocate(nxt_act(nxa_nlin),nxt_wlact(nxa_nlin),nxt_tmact(nxa_nlin))
  nxt_tmact = 1d99
  nxt_x = 0d0; nxt_nel = 0d0; nxt_ne = 0d0; nxt_kq = 0d0; nxt_beta = 1d0; nxt_bes = 0d0
  nxt_jg = 0d0; nxt_jb = 0d0; nxt_ljb = .false.; nxt_jbsum = 0d0; nxt_jbcnt = 0d0
  grd_capnbf = 0.; gas_capnbf = 0.; grd_emnbf = 0.; gas_emnbf = 0.
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
  call nxt_act_lookup
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
  if(allocated(nxt_x)) deallocate(nxt_x,nxt_nel,nxt_ne,nxt_kq,nxt_beta,nxt_bes,nxt_jg,nxt_jb,nxt_ljb, &
       nxt_jbsum,nxt_jbcnt,nxt_jxsum,nxt_jlsum,nxt_jlosum,grd_capnbf,nxt_xloc,nxt_neloc,nxt_nelloc,nxt_neothloc,gas_capnbf, &
       grd_emnbf,gas_emnbf,nxt_jwsum,nxt_jwsq,nxt_haem,nxt_haab,nxt_efsum, &
       nxt_lgrp,nxt_isrt,nxt_nug,nxt_dnug,nxt_act,nxt_wlact,nxt_w,nxt_tmact)
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
  if(allocated(nxt_bfl)) deallocate(nxt_bfl,nxt_bfn,nxt_bfm,nxt_bfr,nxt_bfs,nxt_smaxb)
  if(allocated(nxt_wcnt)) deallocate(nxt_wcnt)
  if(allocated(nxt_lstat)) deallocate(nxt_lstat,nxt_sjlo,nxt_sjhi)
  if(allocated(nxt_wwb)) deallocate(nxt_wwb,nxt_wwr)
  if(allocated(nxt_erel)) deallocate(nxt_erel,nxt_prel,nxt_rcum)
  if(allocated(nxt_xinit)) deallocate(nxt_xinit,nxt_tinit,nxt_neinit)
  nxt_warm = .false.
  nxt_hold = .false.
  nxt_dtinv = 0d0
  nxt_elh = 0; nxt_elhe = 0; nxt_kha = 0
  if(allocated(nxt_ldbd)) deallocate(nxt_ldbd)
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
subroutine nxt_gas_solve(it,gncell,temp,trad,nelv,neoth,qb,qv,icell1,edep,kap,fcl,dr,vel)
!-----------------------------------------------------------------------
! Statistical equilibrium of the NLTE elements in this rank's gas cells.
! The first iteration uses J = B(trad); later ones the Monte Carlo
! estimators of the previous iteration, and the new populations are mixed
! geometrically with the old ones (weight nxt_damp).  edep: power absorbed
! from photons above nxt_ntemin [erg/s/cm^3], deposited by fast electrons.
! kap, fcl, dr (optional): absorption opacity of the groups with which J
! was made, volume filling factor and width of the cells, for the
! accelerated lambda iteration of the bound-free transitions (nxt_ali_cell).
! vel (optional): flow velocity of the cells, for the microturbulence of
! the electron scattering of trapped line photons (nxt_vtb).
!-----------------------------------------------------------------------
  integer,intent(in) :: it,gncell,icell1
  real*8,intent(in) :: temp(gncell),trad(gncell),nelv(nxa_nel,gncell),neoth(gncell)
  real*8,intent(in) :: edep(gncell),qb(gncell),qv(gncell)
  real*4,intent(in),optional :: kap(:,:)
  real*8,intent(in),optional :: fcl(gncell),dr(gncell),vel(gncell)
  integer :: i,l,ig,ierr,nfail,k,nc,nnan
  real*8 :: x(nxa_nst),xold(nxa_nst),ne,beta(nxa_nlin),jg(nxt_ng),jb(nxa_nlin)
  real*8 :: nti(nxa_nst),ntx(nxa_nst),ntot,dxc(gncell),w,efc(nxa_nlin)
  logical :: ljb(nxa_nlin),lfar
!
  nfail = 0
  nnan = 0
  nc = 0
  do i=1,gncell
     l = icell1 + i - 1
     ntot = sum(nelv(:,i))
     nxt_nelloc(:,i) = nelv(:,i)
     if(ntot <= 0d0) then
        nxt_xloc(:,i) = 0d0
        nxt_neloc(i) = neoth(i)
        cycle
     endif
!-- (warm start: the saved state, kept in the first iteration, or in all
!-- with nxt_hold)
     if(nxt_warm .and. (it==1 .or. .not.nxt_have_mc .or. nxt_hold)) then
        nxt_neloc(i) = neoth(i) + free_el(nelv(:,i),nxt_xloc(:,i))
        cycle
     endif
     if(it==1 .or. .not.nxt_have_mc) then
        do ig=1,nxt_ng
           jg(ig) = bplanck(nxt_nug(ig),trad(i))
        enddo
        jb = 0d0
        ljb = .false.
        nx_dis = 0d0
        call nx_lte(temp(i),max(neoth(i),1d-3*ntot),x)
     else
        jg = nxt_jg(:,l)
        jb = nxt_jb(:,l)
        ljb = nxt_ljb(:,l)
        x = nxt_xloc(:,i)
     endif
!-- (nxt_esrates: the absorption in the cores of the Sobolev lines from the
!-- transport, nxt_esfar; the static lines of a slow cell keep bes J_es)
     lfar = nxt_esrates .and. nxt_have_mc .and. it>1
     if(lfar) then
        efc = nxt_efsum(:,l)
        if(nxt_dopmc>0d0 .and. nxt_anystat) then
           if(nxt_slow(l)) then
              do k=1,nxa_nlin
                 if(nxt_lstat(k)) efc(k) = -1d0
              enddo
           endif
        endif
     endif
     xold = x
     call nxt_ntrates(nelv(:,i),x,edep(i),nti,ntx)
!-- (level dissolution: the cell's parameters with the populations from
!-- which the solution starts)
     nx_dis = 0d0
     if(nxd_on) call nxd_cellx(neoth(i)+free_el(nelv(:,i),x),nelv(:,i),x,temp(i),nx_dis)
     nx_vtb = -1d0
     if(present(vel)) call nxt_es_cell(l,vel(i))
     if(present(kap) .and. present(fcl) .and. present(dr) .and. it>1 .and. nxt_have_mc) &
          call nxt_ali_cell(nxt_ng,kap(:,i),fcl(i),dr(i),lali=nxt_alion)
     if(nxt_dtinv>0d0 .and. lfar) then
        call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
             nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,xprev=nxt_xprev(:,l),dtinv=nxt_dtinv, &
             efar=efc)
     elseif(nxt_dtinv>0d0) then
        call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
             nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,xprev=nxt_xprev(:,l),dtinv=nxt_dtinv)
     elseif(lfar) then
        call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
             nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,efar=efc)
     else
        call nx_solve_cell(temp(i),trad(i),nelv(:,i),neoth(i),qb(i),qv(i),nxt_ng, &
             nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr)
     endif
     call nxt_ali_off
     nx_vtb = -1d0
     nx_dis = 0d0
!-- (a solution that is not finite counts as failed; the guess is kept)
     if(any(x/=x) .or. any(abs(x)>huge(x))) then
        if(nnan<5) write(0,'(a,i6,a,es11.3,a,es11.3,a,i2,a,l2)') 'nxt_gas_solve: populations not finite in cell', &
             l,', T =',temp(i),' K, n_atom =',ntot,' cm^-3, ierr',ierr,'; guess finite:',all(xold==xold)
        nnan = nnan + 1
        ierr = 1
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
!-- (convergence: largest relative change of a state above 1e-4)
     nc = nc + 1
     dxc(nc) = 0d0
     do k=1,nxa_nst
        if(x(k) > 1d-4 .and. xold(k) > 0d0) dxc(nc) = max(dxc(nc),abs(x(k)/xold(k) - 1d0))
     enddo
  enddo
  if(nfail>0) write(0,*) 'nxt_gas_solve: no solution in cells:',nfail
!-- median and largest change over the cells (insertion sort)
  nxt_dxstat = 0d0
  if(nc>0) then
     do i=2,nc
        w = dxc(i)
        k = i - 1
        do while(k>=1)
           if(dxc(k) <= w) exit
           dxc(k+1) = dxc(k)
           k = k - 1
        enddo
        dxc(k+1) = w
     enddo
     nxt_dxstat = [dxc((nc+1)/2),dxc(nc)]
  endif
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
pure function nxt_vtb(v) result(vt)
!-- microturbulent velocity [cm/s] at flow velocity v for the electron
!-- scattering of trapped line photons (-1: off; nx_vtb)
  real*8,intent(in) :: v
  real*8 :: vt
  vt = -1d0
  if(.not.nxt_estrap) return
  vt = nxt_vturbv(v)
end function nxt_vtb
!
!
pure function nxt_vturbv(v) result(vt)
!-- microturbulent velocity [cm/s] at flow velocity v (in_nltex_vturb, 0
!-- if not given): linear in v between nxt_vturb(1) and nxt_vturb(2)
  real*8,intent(in) :: v
  real*8 :: vt
  vt = nxt_vturb(1)
  if(nxt_vmax > 0d0) vt = vt + (nxt_vturb(2)-nxt_vturb(1))*min(abs(v)/nxt_vmax,1d0)
  vt = max(vt,0d0)
end function nxt_vturbv
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
subroutine nxt_ldb(ncell,ng,vol,dr,tempinv,cap,capemit,planck,spk,nlump,glump)
!-----------------------------------------------------------------------
! Mean intensity of the groups for the statistical equilibrium, after
! nxt_store_jg.
! * DDMC packets in the lump of a cell's optically thick groups stand for
!   radiation in equilibrium with the source function there: their tally
!   (nxt_jlsum) is shared among the lump's groups as the DDMC source
!   spectrum (spk, proportional to S dnu), not given to the group that the
!   packet happens to carry.
! * Local detailed balance (nxt_ldbtau > 0): a group with absorption depth
!   tau = kappa dr >= nxt_ldbtau across the cell gets
!      J = Lambda* S + J_ext,
!   S = (kappa_th B + eta_nbf)/kappa the source function of the populations
!   with which the packets ran, Lambda* (nxt_lsldb) the share of its own
!   emission that a slab of depth tau keeps, and J_ext the radiation from
!   elsewhere (nxt_jxsum: packets on IMC steps outside the cell of their
!   last thermal or bound-free emission; the DDMC lump is local).  With
!   DDMC (nxt_ldbone) these groups are DDMC groups, in the diffusion limit:
!   Lambda* = 1, J = S.  The
!   Monte Carlo estimate of such a group is the cell's own emission as
!   sampled by the few packets emitted in it: zero where the group's share
!   of the emission is small (the Lyman continuum of a cool wind), which
!   turns the on-the-spot cycle into net recombination.  The statistical
!   equilibrium takes the cell's own cycle out with the same Lambda*
!   (nxt_ali_cell).
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell,ng
  real*8,intent(in) :: vol(ncell),dr(ncell),tempinv(ncell)
  real*4,intent(in) :: cap(ng,ncell),capemit(ng,ncell)
  real*8,intent(in) :: planck(ng,ncell),spk(ng,ncell)
  integer,intent(in) :: nlump(ncell)
  integer*2,intent(in) :: glump(ng,ncell)
  real*8,parameter :: sb = 5.670374419d-5
  integer :: i,ig,k,ib
  real*8 :: w,tau,s,u,ls,jx,wl
  logical :: ldiag
!
  ldiag = allocated(nxt_nug)
  if(ldiag) ldiag = size(nxt_nug)==ng
  if(allocated(nxt_ldbd)) then
     if(size(nxt_ldbd,3)/=ncell) deallocate(nxt_ldbd)
  endif
  if(.not.allocated(nxt_ldbd)) allocate(nxt_ldbd(3,2,ncell))
  nxt_ldbd = 0d0
  do i=1,ncell
     if(vol(i)<=0d0) cycle
!-- DDMC lump: shared as the source spectrum
     if(nxt_jlsum(i)>0d0 .and. nlump(i)>0) then
        w = 0d0
        do k=1,nlump(i)
           w = w + spk(glump(k,i),i)
        enddo
        do k=1,nlump(i)
           ig = glump(k,i)
           if(w>0d0) then
              nxt_jg(ig,i) = nxt_jg(ig,i) + nxt_jlsum(i)*spk(ig,i)/w/(4d0*pi*vol(i)*nxt_dnug(ig))
           else
              nxt_jg(ig,i) = nxt_jg(ig,i) + nxt_jlsum(i)/dble(nlump(i))/(4d0*pi*vol(i)*nxt_dnug(ig))
           endif
        enddo
     endif
!-- local detailed balance
     if(nxt_ldbtau<=0d0) cycle
     u = sb/(pi*tempinv(i)**4)
     do ig=1,ng
        tau = cap(ig,i)*dr(i)
        if(tau<nxt_ldbtau) cycle
        s = (planck(ig,i)*capemit(ig,i) + grd_emnbf(ig,i))/cap(ig,i)*u/nxt_dnug(ig)
        if(nxt_ldbone) then
           ls = 1d0
        else
           ls = nxt_lsldb(tau)
        endif
        jx = nxt_jxsum(ig,i)/(4d0*pi*vol(i)*nxt_dnug(ig))
        ib = 0
        if(ldiag) then
           wl = cc/nxt_nug(ig)
           if(wl>=228d-8 .and. wl<504d-8) ib = 1
           if(wl>=504d-8 .and. wl<912d-8) ib = 2
        endif
        if(ib>0) nxt_ldbd(:,ib,i) = nxt_ldbd(:,ib,i) + 4d0*pi*nxt_dnug(ig)*[ls*s,jx,nxt_jg(ig,i)]
        nxt_jg(ig,i) = ls*s + jx
     enddo
  enddo
end subroutine nxt_ldb
!
!
pure function nxt_lsldb(tau) result(ls)
!-- share of the isotropic emission of a slab of absorption depth tau
!-- absorbed in it (mean over the slab), 1 - (1 - 2 E_3(tau))/(2 tau)
  real*8,intent(in) :: tau
  real*8 :: ls
  if(tau < 1d-6) then
     ls = 0d0
  else
     ls = 1d0 - (1d0 - 2d0*nxt_e3(tau))/(2d0*tau)
  endif
end function nxt_lsldb
!
!
pure function nxt_e3(x) result(e3)
!-- exponential integral E_3(x), x > 0: E_1 from Abramowitz & Stegun 5.1.53
!-- (x <= 1) and 5.1.56, then the recurrence E_n+1 = (e^-x - x E_n)/n up to
!-- x = 5 and the asymptotic series above (where the recurrence cancels)
  real*8,intent(in) :: x
  real*8 :: e3,e1,e2,y
  if(x > 5d0) then
     y = 1d0/(x+3d0)
     e3 = exp(-x)*y*(1d0 + 3d0*y**2 + 3d0*(3d0-2d0*x)*y**4)
     return
  endif
  if(x <= 1d0) then
     e1 = -log(x) - 0.57721566d0 + x*(0.99999193d0 + x*(-0.24991055d0 + x*(0.05519968d0 + &
          x*(-0.00976004d0 + x*0.00107857d0))))
  else
     e1 = exp(-x)/x*(x**4 + 8.5733287401d0*x**3 + 18.0590169730d0*x**2 + 8.6347608925d0*x + 0.2677737343d0)/ &
          (x**4 + 9.5733223454d0*x**3 + 25.6329561486d0*x**2 + 21.0996530827d0*x + 3.9584969228d0)
  endif
  e2 = exp(-x) - x*e1
  e3 = 0.5d0*(exp(-x) - x*e2)
end function nxt_e3
!
!
subroutine nxt_tally_reset
  nxt_jbsum = 0d0
  nxt_jbcnt = 0d0
  nxt_jxsum = 0d0
  nxt_jlsum = 0d0
  nxt_jlosum = 0d0
  nxt_jwsum = 0d0
  nxt_jwsq = 0d0
  nxt_haem = 0d0
  nxt_haab = 0d0
  nxt_efsum = 0d0
end subroutine nxt_tally_reset
!
!
subroutine nxt_tally_finish(ncell,vol)
!-----------------------------------------------------------------------
! Mean intensity incident on each NLTE line (after nxt_store_jg and
! nxt_ldb): the resonance tally (sum of e c/(nu |Q|)) where the line has
! nxt_njbmin crossings, else the J of the line's group without the photons
! of its lines that are still to the red of their line in the comoving
! frame (nxt_jlosum; packet%ilem).  In an expanding flow those only
! redshift away from it: the Sobolev incident intensity is that of the continuum
! (CMFGEN's Sobolev option).  The whole group J, which the line's own
! escaping photons dominate where it is strong (Halpha in the outer wind of
! a IIn: a few crossings per cell; 4-5 times the continuum in a 1000 km/s
! group), pumped the upper level.  With nxt_jbwin, first the window
! estimator (nxt_window_tally) where it has nxt_njbmin effective samples.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: vol(ncell)
  integer :: i,k,g
  do i=1,ncell
     nxt_jb(:,i) = nxt_jbsum(:,i)/(4d0*pi*vol(i))
     nxt_ljb(:,i) = nxt_jbcnt(:,i) >= nxt_njbmin
!-- (nxt_jbwin: the window estimator where it has nxt_njbmin effective samples)
     if(nxt_jbwin>0d0) then
        do k=1,nxa_nlin
           if(nxt_jwsq(k,i) <= 0d0) cycle
           if(nxt_jwsum(k,i)**2 < nxt_njbmin*nxt_jwsq(k,i)) cycle
           nxt_jb(k,i) = nxt_jwsum(k,i)/(4d0*pi*vol(i))
           nxt_ljb(k,i) = .true.
        enddo
     endif
     do k=1,nxa_nlin
        if(nxt_ljb(k,i)) cycle
        g = nxt_lgrp(k)
        if(g<1) cycle
        nxt_jb(k,i) = max(nxt_jg(g,i) - nxt_jlosum(g,i)/(4d0*pi*vol(i)*nxt_dnug(g)),0d0)
        nxt_ljb(k,i) = .true.
     enddo
!-- (the far-side core absorption of the electron-scattered photons, nxt_esfar)
     if(nxt_esrates) then
        do k=1,nxa_nlin
           if(nxt_kq(k,i) > 0d0) then
              nxt_efsum(k,i) = nxt_efsum(k,i)*cc/(4d0*pi*vol(i)*nxa_lnu(k)*nxt_kq(k,i))
           else
              nxt_efsum(k,i) = 0d0
           endif
        enddo
     endif
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
  real*8 :: bnd(nxa_nlv),bnb(nxa_nlv) !bound-free weights and their bounds in a block of groups
  integer :: ixb(nxa_nlv)
  real*8 :: jg(nxt_ng),jb(nxa_nlin),jl(nxa_nlin),btab(0:nx_nbt),etab(0:nx_nbt),sig,vt,bu(nxa_nlin)
  logical :: ljb(nxa_nlin),lst
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
     if(nxt_dopmc > 0d0) then
        if(.not.allocated(nxt_slow)) allocate(nxt_slow(ncell),nxt_vdl(nxa_nel,ncell), &
             nxt_vdlmax(ncell),nxt_vdlmin(ncell))
        vt = nxt_vturbv(0.5d0*(vxarr(ic)+vxarr(ic+1)))
        do iel=1,nxa_nel
           nxt_vdl(iel,ic) = nx_vdop(iel,t,vt)
        enddo
        nxt_vdlmax(ic) = maxval(nxt_vdl(:,ic))
        nxt_vdlmin(ic) = minval(nxt_vdl(:,ic))
        nxt_slow(ic) = max(abs(vxarr(ic+1)-vxarr(ic)),abs(qv)*(xarr(ic+1)-xarr(ic))) < &
             nxt_dopmc*nxt_vdlmin(ic)
     endif
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
!-- level dissolution: the cell's parameters, the dissolved shares, the
!-- lowest frequency of each level's dissolved continuum
     nx_dis = 0d0
     if(nxd_on) then
        if(.not.allocated(nxt_disp)) allocate(nxt_disp(3,ncell),nxt_dnuc(nxa_nlv,ncell),nxt_fdis(nxa_nlv,ncell))
        call nxd_cellx(ne,nxt_nel(:,ic),nxt_x(:,ic),t,nx_dis)
        nxt_disp(:,ic) = nx_dis
        nxt_fdis(:,ic) = 0.
        do lv=1,nxa_nlv
           nxt_dnuc(lv,ic) = nxa_lvnu0(lv)
           if(nxd_ga(lv) > 0) nxt_dnuc(lv,ic) = nxd_numin(nxa_lvnu0(lv),nxd_z(lv),nx_dis)
        enddo
     endif
     call nxd_weights(t,nx_dis,nxt_w(:,ic))
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
!-- electron scattering of the trapped line photons (volume-mean Thomson
!-- opacity with the volume-mean K); in the lines of nxt_useb the macro-atom
!-- emits with beta + bes and absorbs (beta Jb + bes Jes)
     nxt_bes(:,ic) = 0d0
     if(nxt_estrap .and. nxt_tau0<=0d0) then
        call nx_es_table(qb,qv,etab)
        sig = 6.6524587d-25*nxt_ne(ic)
        if(allocated(nxt_fcl)) sig = sig*nxt_fcl(ic)
        vt = nxt_vtb(0.5d0*(vxarr(ic)+vxarr(ic+1)))
        do k=1,nxa_nlin
           nxt_bes(k,ic) = nx_bes(nxt_kq(k,ic),sig*nx_vdop(nxa_lel(k),t,vt), &
                merge(nx_wq(k,t),0d0,nxt_esdamp),nxt_beta(k,ic),qb,qv,etab)
           if(.not.nxt_useb(k,ic)) cycle
!-- (nxt_esrates: the core absorption of the transport, nxt_esfar, for the
!-- Sobolev lines; the static lines of a slow cell keep bes J_es)
           lst = .false.
           if(nxt_dopmc>0d0 .and. nxt_anystat) lst = nxt_slow(ic) .and. nxt_lstat(k)
           if(nxt_esrates .and. nxt_have_mc .and. .not.lst) then
              jl(k) = (nxt_beta(k,ic)*jl(k) + nxt_efsum(k,ic))/(nxt_beta(k,ic) + nxt_bes(k,ic))
              cycle
           endif
           if(nxt_bes(k,ic) <= 0d0) cycle
           jl(k) = (nxt_beta(k,ic)*jl(k) + nxt_bes(k,ic)*nx_jes(k,t,nxt_ng,nxt_nug,nxt_dnug,jg))/ &
                (nxt_beta(k,ic) + nxt_bes(k,ic))
        enddo
     endif
     if(nxt_tau0>0d0 .and. nxt_kha>0) then
!-- benchmark: only H-alpha, tau = tau0 for a gradient v/r
        nxt_kq(:,ic) = 0d0
        nxt_kq(nxt_kha,ic) = nxt_tau0*qv
     endif
!-- escape probability of the macro-atom
     do k=1,nxa_nlin
        bu(k) = nxt_beta(k,ic)
        if(nxt_useb(k,ic)) bu(k) = bu(k) + nxt_bes(k,ic)
     enddo
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
             jl(l0:l1),bu(l0:l1),q,tot,emt,kex,kin,chw,ch0)
!-- recombination cooling per level of this element (from level 1 of the next ion)
        do k=1,nxa_el(iel)%nlv
           lv = lo + k
           kw(1+lv) = nxt_nel(iel,ic)*nxt_x(o+nxa_lvc(lv),ic)*r%wc(k)*ne*nxa_rcool(lv,t)*nxd_lvw(lv,nx_dis)
!-- (the share of the recombinations into the dissolved continuum)
           if(nxd_on) then
              if(r%arl(k) > 0d0) nxt_fdis(lv,ic) = real(r%ard(k)/r%arl(k))
           endif
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
             nxa_la(k)*bu(k)*hh*nxa_lnu(k)
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
  nx_dis = 0d0
!-- bound-free absorption: per block of groups the levels of largest bound
!-- of their weight, their weight at each group of the block, and the bound
!-- of all others
  if(allocated(nxt_bfl)) then
     if(size(nxt_bfl,1)/=nxt_bfk .or. size(nxt_bfl,3)/=ncell) &
          deallocate(nxt_bfl,nxt_bfn,nxt_bfm,nxt_bfr,nxt_bfs,nxt_smaxb)
  endif
  if(.not.allocated(nxt_bfl)) then
     nxt_bsz = (nxt_ng + nxt_nbb - 1)/nxt_nbb
     nxt_nwm = (nxa_nlv + 63)/64
     allocate(nxt_smaxb(nxa_nlv,nxt_nbb),nxt_bfl(nxt_bfk,nxt_nbb,ncell),nxt_bfn(nxt_nbb,ncell), &
          nxt_bfm(nxt_nwm,nxt_nbb,ncell),nxt_bfr(nxt_nbb,ncell),nxt_bfs(nxt_ng,ncell))
     nxt_smaxb = 0.
     do lv=1,nxa_nlv
        do ig=1,min(nxa_g1(lv),nxt_ng)
           j = (ig-1)/nxt_bsz + 1
           nxt_smaxb(lv,j) = max(nxt_smaxb(lv,j),nxa_sgm(nxa_sgo(lv)+ig))
        enddo
!-- (level dissolution: the dissolved continuum below the edge, at most the
!-- threshold cross-section)
        if(nxd_on) then
           if(nxd_ga(lv) > 0) then
              do ig=nxd_ga(lv),min(nxd_gb(lv),nxt_ng)
                 j = (ig-1)/nxt_bsz + 1
                 nxt_smaxb(lv,j) = max(nxt_smaxb(lv,j),real(nxd_sth(lv)))
              enddo
           endif
        endif
     enddo
  endif
  nxt_bfstat = 0d0
  do ic=1,ncell
     do lv=1,nxa_nlv
        iel = nxa_lvel(lv)
        bnd(lv) = nxt_nel(iel,ic)*nxt_x(nxa_el(iel)%xoff+nxa_lvs(lv),ic)*nxt_w(lv,ic)
     enddo
     do j=1,nxt_nbb
        do lv=1,nxa_nlv
           bnb(lv) = bnd(lv)*nxt_smaxb(lv,j)
           ixb(lv) = lv
        enddo
        k = min(nxt_bfk,nxa_nlv)
        call top_k(bnb,k,ixb)
        nxt_bfr(j,ic) = sum(bnb(ixb(k+1:nxa_nlv)))
        call sort_desc(bnb,ixb(1:k))
!-- (levels without weight in the block stay out)
        do while(k>0)
           if(bnb(ixb(k))>0d0) exit
           k = k - 1
        enddo
        nxt_bfn(j,ic) = k
        nxt_bfl(1:k,j,ic) = ixb(1:k)
        nxt_bfm(:,j,ic) = 0_8
        do m=1,k
           lv = ixb(m)
           nxt_bfm((lv-1)/64+1,j,ic) = ibset(nxt_bfm((lv-1)/64+1,j,ic),mod(lv-1,64))
        enddo
        do ig=(j-1)*nxt_bsz+1,min(j*nxt_bsz,nxt_ng)
           nxt_bfs(ig,ic) = 0d0
           do m=1,k
              lv = ixb(m)
              nxt_bfs(ig,ic) = nxt_bfs(ig,ic) + bnd(lv)*bfsig(lv,ig,ic)
           enddo
!-- (share of the absorptions at the group that walk the other levels)
           if(nxt_bfr(j,ic)>0d0) then
              ww = nxt_bfr(j,ic)/(nxt_bfr(j,ic) + nxt_bfs(ig,ic))
              nxt_bfstat(1) = nxt_bfstat(1) + 1d0
              nxt_bfstat(2) = nxt_bfstat(2) + ww
              if(ww>0.1d0) nxt_bfstat(3) = nxt_bfstat(3) + 1d0
           elseif(nxt_bfs(ig,ic)>0d0) then
              nxt_bfstat(1) = nxt_bfstat(1) + 1d0
           endif
        enddo
     enddo
  enddo
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
        nxt_tmact(nxt_nact) = ww
        nxt_isrt(k) = nxt_nact
     endif
  enddo
  call nxt_act_lookup
end subroutine nxt_tables
!
!
subroutine nxt_act_lookup
!-- lookup of the active-line list by ln(wavelength) (nxt_neighbors), and
!-- the lines with static profiles in slow cells (nxt_dopz)
  integer :: k,j
  real*8 :: e
  if(.not.allocated(nxt_lstat)) allocate(nxt_lstat(nxa_nlin),nxt_sjlo(0:nxa_nlin+1),nxt_sjhi(0:nxa_nlin+1))
  do k=1,nxa_nlin
     nxt_lstat(k) = nxa_el(nxa_lel(k))%z >= nxt_dopz
  enddo
  nxt_sjlo(0) = 0
  do j=1,nxt_nact
     nxt_sjlo(j) = merge(nxt_sjlo(j-1),j,nxt_lstat(nxt_act(j)))
  enddo
  nxt_sjhi(nxt_nact+1) = 0
  do j=nxt_nact,1,-1
     nxt_sjhi(j) = merge(nxt_sjhi(j+1),j,nxt_lstat(nxt_act(j)))
  enddo
  nxt_anystat = any(nxt_lstat(nxt_act(1:nxt_nact)))
!-- windows of the incident-intensity estimator (nxt_window_tally)
  if(.not.allocated(nxt_wwb)) allocate(nxt_wwb(nxa_nlin),nxt_wwr(nxa_nlin))
  nxt_wwb = 0d0
  nxt_wwr = 0d0
  if(nxt_jbwin>0d0) then
!-- (up to the nearest line with a Sobolev depth of nxt_tauwin somewhere)
     e = 0d0
     do j=1,nxt_nact
        nxt_wwb(j) = nxt_jbwin
        if(e>0d0) nxt_wwb(j) = min(nxt_jbwin,(1d0 - e/nxt_wlact(j))*cc)
        if(nxt_tmact(j) >= nxt_tauwin) e = nxt_wlact(j)
     enddo
     e = 0d0
     do j=nxt_nact,1,-1
        nxt_wwr(j) = nxt_jbwin
        if(e>0d0) nxt_wwr(j) = min(nxt_jbwin,(e/nxt_wlact(j) - 1d0)*cc)
        if(nxt_tmact(j) >= nxt_tauwin) e = nxt_wlact(j)
     enddo
  endif
  if(.not.allocated(nxt_wcnt)) allocate(nxt_wcnt(0:nxt_nwt))
  nxt_wcnt = 0
  if(nxt_nact<1) return
  nxt_lw0 = log(nxt_wlact(1))
  nxt_dlw = max(log(nxt_wlact(nxt_nact)) - nxt_lw0,1d-6)/nxt_nwt
  j = 0
  do k=0,nxt_nwt
     e = exp(nxt_lw0 + k*nxt_dlw)
     do while(j < nxt_nact)
        if(nxt_wlact(j+1) > e) exit
        j = j + 1
     enddo
     nxt_wcnt(k) = j
  enddo
end subroutine nxt_act_lookup
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
pure subroutine nxt_static_line(ic,wlc,kown,r,kap,kpick)
!-----------------------------------------------------------------------
! Opacity [1/cm] of the active lines within nxt_dopx Doppler widths of the
! comoving wavelength wlc in slow cell ic, static Doppler profiles
! K/(sqrt(pi) v_D) exp(-x^2), of the lines with static profiles
! (nxt_lstat) without the line kown (global index; the packet's own line,
! 0 if none).  With r in [0,1) also the line that
! absorbs (kpick, global index; chosen by its share of the opacity).
!-----------------------------------------------------------------------
  integer,intent(in) :: ic,kown
  real*8,intent(in) :: wlc,r
  real*8,intent(out) :: kap
  integer,intent(out) :: kpick
  real*8,parameter :: c = 2.99792458d10, spi = 1.7724538509055159d0
  integer :: jlo,jhi,j,k,pass
  real*8 :: wmin,wmax,x,kl,cum,target
  kap = 0d0
  kpick = 0
  if(nxt_nact <= 0) return
  wmin = wlc/(1d0 + nxt_dopx*nxt_vdlmax(ic)/c)
  wmax = wlc/(1d0 - nxt_dopx*nxt_vdlmax(ic)/c)
  call nxt_neighbors(wlc,0,0,jlo,jhi)
!-- (jlo, jhi: the lines just below and above, 0 if none)
!-- the sum, then (r >= 0) the line where the cumulative sum passes r*kap
  target = -1d0
  cum = 0d0
  do pass=1,2
     if(pass == 2) then
        if(r < 0d0 .or. kap <= 0d0) return
        target = r*kap
     endif
     do j=jlo,1,-1
        if(nxt_wlact(j) < wmin) exit
        k = nxt_act(j)
        if(k == kown .or. .not.nxt_lstat(k)) cycle
        x = (wlc/nxt_wlact(j) - 1d0)*c/nxt_vdl(nxa_lel(k),ic)
        if(abs(x) > nxt_dopx) cycle
        kl = max(nxt_kq(k,ic),0d0)*exp(-x*x)/(spi*nxt_vdl(nxa_lel(k),ic))
        if(pass == 1) then
           kap = kap + kl
        else
           cum = cum + kl
           kpick = k
           if(cum >= target) return
        endif
     enddo
     do j=merge(jhi,nxt_nact+1,jhi>0),nxt_nact
        if(nxt_wlact(j) > wmax) exit
        k = nxt_act(j)
        if(k == kown .or. .not.nxt_lstat(k)) cycle
        x = (wlc/nxt_wlact(j) - 1d0)*c/nxt_vdl(nxa_lel(k),ic)
        if(abs(x) > nxt_dopx) cycle
        kl = max(nxt_kq(k,ic),0d0)*exp(-x*x)/(spi*nxt_vdl(nxa_lel(k),ic))
        if(pass == 1) then
           kap = kap + kl
        else
           cum = cum + kl
           kpick = k
           if(cum >= target) return
        endif
     enddo
  enddo
end subroutine nxt_static_line
!
!
pure subroutine nxt_static_tally(ic,wlc,ed,kown,q,dl,jbsum,jbcnt)
!-----------------------------------------------------------------------
! Incident intensity of the lines of slow cell ic, taken where the Sobolev
! theory takes it: outside the profile, on the side the comoving frequency
! of the flight drifts in from (blue for a positive projected velocity
! gradient q of the flight, red for a negative one).  A flight at comoving
! wavelength wlc with energy times length ed (comoving) between nxt_dopx
! and nxt_dopx+nxt_dopw Doppler widths from a line on that side adds
! ed c/(nu v_D nxt_dopw), the path-length estimator of J_nu over the
! window; a drift through it adds the Sobolev crossing e c/(nu |q|).
! (Inside the profile of a thick line the photons drifting in are absorbed
! in the wing, so the intensity there is about 1/tau of the incident one,
! and the rate equations apply the Sobolev escape to J_b themselves.)
! The packet's own line (kown) is left out.  The count of a line (jbcnt,
! for nxt_njbmin) grows by the share of the window the flight of comoving
! length dl drifts across, |q| dl/(v_D nxt_dopw) (at most 1), as a
! Sobolev crossing counts 1: the flights of a slow cell are short, and
! counted one each they let a few packets pass the threshold.  The
! tallies (jbsum, jbcnt: nxt_jbsum, nxt_jbcnt) are arguments, so that a
! threaded caller passes its reduction copies.
!-----------------------------------------------------------------------
  integer,intent(in) :: ic,kown
  real*8,intent(in) :: wlc,ed,q,dl
  real*8,intent(inout) :: jbsum(:,:),jbcnt(:,:)
  real*8,parameter :: c = 2.99792458d10
  integer :: jlo,jhi,j,k
  real*8 :: wlim,x,vd
  if(nxt_nact <= 0) return
  call nxt_neighbors(wlc,0,0,jlo,jhi)
  if(q >= 0d0) then
!-- the lines to the red of the flight
     wlim = wlc/(1d0 - (nxt_dopx+nxt_dopw)*nxt_vdlmax(ic)/c)
     do j=merge(jhi,nxt_nact+1,jhi>0),nxt_nact
        if(nxt_wlact(j) > wlim) exit
        k = nxt_act(j)
        if(k == kown .or. .not.nxt_lstat(k)) cycle
        vd = nxt_vdl(nxa_lel(k),ic)
        x = (wlc/nxt_wlact(j) - 1d0)*c/vd
        if(x >= -nxt_dopx .or. x < -nxt_dopx-nxt_dopw) cycle
        jbsum(k,ic) = jbsum(k,ic) + ed*c/(nxa_lnu(k)*vd*nxt_dopw)
        jbcnt(k,ic) = jbcnt(k,ic) + min(1d0,abs(q)*dl/(vd*nxt_dopw))
     enddo
  else
!-- the lines to the blue of the flight
     wlim = wlc/(1d0 + (nxt_dopx+nxt_dopw)*nxt_vdlmax(ic)/c)
     do j=jlo,1,-1
        if(nxt_wlact(j) < wlim) exit
        k = nxt_act(j)
        if(k == kown .or. .not.nxt_lstat(k)) cycle
        vd = nxt_vdl(nxa_lel(k),ic)
        x = (wlc/nxt_wlact(j) - 1d0)*c/vd
        if(x <= nxt_dopx .or. x > nxt_dopx+nxt_dopw) cycle
        jbsum(k,ic) = jbsum(k,ic) + ed*c/(nxa_lnu(k)*vd*nxt_dopw)
        jbcnt(k,ic) = jbcnt(k,ic) + min(1d0,abs(q)*dl/(vd*nxt_dopw))
     enddo
  endif
end subroutine nxt_static_tally
!
!
pure subroutine nxt_window_tally(ic,wla,wlb,ed,jwsum,jwsq)
!-----------------------------------------------------------------------
! Incident intensity of the Sobolev lines of cell ic (nxt_jbwin > 0): the
! path-length estimator of J_nu in a window next to each line, on the side
! the comoving frequency of the flight drifts in from (to the blue of the
! line for a flight whose comoving wavelength grows from wla to wlb, to the
! red for one whose wavelength shrinks), nxt_jbwin wide or up to the next
! active line (nxt_wwb, nxt_wwr).  The part f of the flight (comoving energy
! times length ed) in a window of width w adds ed f c/(nu w), and its square
! goes to jwsq for the effective number of samples, jwsum^2/jwsq: a drift
! across the window adds e c/(nu |Q|), as a resonance crossing, and a flight
! that stays in it its path.  Where the flights drift a small part of the
! window between electron scatterings (a slow, electron-scattering CSM) it
! has many samples for each crossing, and unlike the group J without the
! line's own photons (nxt_tally_finish) it counts the photons of the line
! that drift back into it: those scattered by electrons to the incident
! side, and in a decelerating or non-monotonic flow those emitted elsewhere
! whose projected velocity gradient changes sign along their path.  The
! line's own escaping photons drift away from it on the other side.
! (Lines with static profiles in a slow cell: nxt_static_tally.)
!-----------------------------------------------------------------------
  integer,intent(in) :: ic
  real*8,intent(in) :: wla,wlb,ed
  real*8,intent(inout) :: jwsum(:,:),jwsq(:,:)
  integer :: jlo,jhi,j,k
  real*8 :: wl0,w,a,b,v
  logical :: lst
  if(nxt_nact <= 0 .or. wla == wlb .or. ed <= 0d0) return
  lst = .false.
  if(nxt_dopmc>0d0 .and. nxt_anystat) lst = nxt_slow(ic)
  call nxt_neighbors(wla,0,0,jlo,jhi)
  if(wlb > wla) then
!-- the windows to the blue of the lines to the red of the flight
     if(jhi <= 0) return
     do j=jhi,nxt_nact
        wl0 = nxt_wlact(j)
        if(wl0*(1d0 - nxt_jbwin/cc) >= wlb) exit
        w = nxt_wwb(j)
        if(w <= 0d0) cycle
        k = nxt_act(j)
        if(lst) then
           if(nxt_lstat(k)) cycle
        endif
        a = max(wla,wl0*(1d0 - w/cc))
        b = min(wlb,wl0)
        if(b <= a) cycle
        v = ed*(b - a)/(wlb - wla)*cc/(nxa_lnu(k)*w)
        jwsum(k,ic) = jwsum(k,ic) + v
        jwsq(k,ic) = jwsq(k,ic) + v*v
     enddo
  else
!-- the windows to the red of the lines to the blue of the flight
     if(jlo <= 0) return
     do j=jlo,1,-1
        wl0 = nxt_wlact(j)
        if(wl0*(1d0 + nxt_jbwin/cc) <= wlb) exit
        w = nxt_wwr(j)
        if(w <= 0d0) cycle
        k = nxt_act(j)
        if(lst) then
           if(nxt_lstat(k)) cycle
        endif
        a = max(wlb,wl0)
        b = min(wla,wl0*(1d0 + w/cc))
        if(b <= a) cycle
        v = ed*(b - a)/(wla - wlb)*cc/(nxa_lnu(k)*w)
        jwsum(k,ic) = jwsum(k,ic) + v
        jwsq(k,ic) = jwsq(k,ic) + v*v
     enddo
  endif
end subroutine nxt_window_tally
!
!
pure subroutine nxt_esfar(ic,wlc,q,vflow,temp,sig,r1,r2,kres,ksk)
!-----------------------------------------------------------------------
! (nxt_esrates) A photon scattered by a thermal electron to comoving
! wavelength wlc in cell ic (Thomson opacity sig [1/cm]), with projected
! velocity gradient q along its new direction, within the Doppler profile
! of a line: the Sobolev point resonance would absorb it only where its
! comoving frequency drifts through the line centre, and not at all
! behind the line, but the profile absorbs it where it lands.  At x Doppler
! widths from the centre, positive on the side it drifts in from, it
! drifts through the rest of the profile while electrons scatter it again
! (a = sig v_D/|q| per Doppler width of drift); the depth of the line
! before its next scattering is
!   tau = (K/|q|) Int_0^inf phi(x - u) exp(-a u) du
!       = (K/2|q|) exp(-x^2) erfcx(a/2 - x).
! Without the scattering (a = 0) this is the Sobolev depth K/|q| ahead of
! the line and (K/|q|) erfc(|x|)/2 behind it, and with frequent scattering
! (a >> 1) K phi(x)/(sig v_D): the profile against the electrons.  A slow,
! electron-scattering flow (the deep wind, a dense CSM) scatters a photon
! many times before it drifts a Doppler width, so it never reaches the
! centre: there the profile absorbs the photons that electrons scatter
! into the cores.  It is absorbed with probability 1 - exp(-sum of the
! depths), by a line in proportion to its depth (kres; 0: not absorbed).
! The lines ahead of the photon (x > 0; ksk, at most four) have then been
! tested for its whole path: the transport passes their resonances without
! absorption until the next scattering.  The rates take this absorption
! from its tally (nxt_efsum) instead of bes J_es.  (Lines with static
! profiles in a slow cell absorb in their profiles.)
!-----------------------------------------------------------------------
  integer,intent(in) :: ic
  real*8,intent(in) :: wlc,q,vflow,temp,sig,r1,r2
  integer,intent(out) :: kres,ksk(4)
  real*8,parameter :: xmax = 5d0
  integer,parameter :: mk = 16
  integer :: jlo,jhi,j,k,n,iel,iside,nsk
  integer :: kk(mk)
  real*8 :: tt(mk),vt,vdm,vd,x,aq,tsum,s,dv,a
  logical :: lst
  kres = 0
  ksk = 0
  nsk = 0
  if(nxt_nact <= 0) return
  aq = max(abs(q),1d-30)
  vt = nxt_vtb(vflow)
  vdm = 0d0
  do iel=1,nxa_nel
     vdm = max(vdm,nx_vdop(iel,temp,vt))
  enddo
  lst = .false.
  if(nxt_dopmc>0d0 .and. nxt_anystat) lst = nxt_slow(ic)
  call nxt_neighbors(wlc,0,0,jlo,jhi)
  n = 0
  tsum = 0d0
!-- the lines to the blue of the photon (iside = 1), then to the red (2)
  do iside=1,2
     if(iside == 1) then
        j = jlo
     else
        j = jhi
        if(j <= 0) exit
     endif
     do while(j >= 1 .and. j <= nxt_nact)
        dv = (wlc/nxt_wlact(j) - 1d0)*cc
        if(abs(dv) > xmax*vdm) exit
        k = nxt_act(j)
        if(iside == 1) then
           j = j - 1
        else
           j = j + 1
        endif
        if(lst) then
           if(nxt_lstat(k)) cycle
        endif
        vd = nx_vdop(nxa_lel(k),temp,vt)
!-- offset in Doppler widths, positive on the side the photon drifts in
!-- from (the blue side for q > 0: the photon to the red of a line is behind it)
        x = -dv/vd
        if(q < 0d0) x = -x
        if(abs(x) > xmax) cycle
        a = sig*vd/aq
        s = 0.5d0*max(nxt_kq(k,ic),0d0)/aq*exp(-x*x)*erfc_scaled(0.5d0*a - x)
        if(s < 1d-6) cycle
        n = n + 1
        kk(n) = k
        tt(n) = s
        tsum = tsum + s
        if(x > 0d0 .and. nsk < 4) then
           nsk = nsk + 1
           ksk(nsk) = k
        endif
        if(n == mk) exit
     enddo
     if(n == mk) exit
  enddo
  if(n == 0) return
  if(r1 >= 1d0 - exp(-tsum)) return
  s = r2*tsum
  do j=1,n
     s = s - tt(j)
     if(s <= 0d0) exit
  enddo
  kres = kk(min(j,n))
  ksk = 0
end subroutine nxt_esfar
!
!
pure subroutine nxt_sobolev_neighbors(jlo,jhi)
!-- in a slow cell: from the neighbors jlo, jhi of nxt_neighbors, the
!-- nearest lines that keep their Sobolev resonances there (0: none)
  integer,intent(inout) :: jlo,jhi
  if(jlo>0) jlo = nxt_sjlo(jlo)
  if(jhi>0) jhi = nxt_sjhi(jhi)
end subroutine nxt_sobolev_neighbors
!
!
pure logical function nxt_useb(k,ic)
!-- the electron scattering of the photons trapped in line k is in the
!-- rates and the macro-atom of cell ic: all lines with nxt_esrates, the
!-- static lines of the slow cells with nxt_dopes
  integer,intent(in) :: k,ic
  nxt_useb = nxt_esrates
  if(nxt_useb .or. .not.nxt_dopes) return
  if(.not.allocated(nxt_slow) .or. .not.allocated(nxt_lstat)) return
  nxt_useb = nxt_slow(ic) .and. nxt_lstat(k)
end function nxt_useb
!
!
subroutine nxt_es_cell(l,vel)
!-- the electron scattering of trapped line photons in the statistical
!-- equilibrium of cell l (flow velocity vel): nltexmod's nx_vtb and the
!-- lines that have it, nx_esm (nxt_useb)
  integer,intent(in) :: l
  real*8,intent(in) :: vel
  integer :: k
  nx_vtb = -1d0
  if(nxt_esrates) then
     nx_vtb = nxt_vtb(vel)
     if(allocated(nx_esm)) nx_esm = .true.
  elseif(nxt_dopes) then
     if(.not.allocated(nxt_slow)) return
     if(.not.nxt_slow(l)) return
     nx_vtb = nxt_vtb(vel)
     if(.not.allocated(nx_esm)) allocate(nx_esm(nxa_nlin))
     do k=1,nxa_nlin
        nx_esm(k) = nxt_useb(k,l)
     enddo
  endif
end subroutine nxt_es_cell
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
  integer :: a,b,m,k
  real*8 :: u
  if(iline>0) then
     if(ilside>0) then
        jlo = iline
        jhi = iline + 1
     else
        jlo = iline - 1
        jhi = iline
     endif
  else
!-- bracket from the lookup, one bin wider on each side (rounding):
!-- nxt_wlact(a) <= wlc < nxt_wlact(b)
     a = 0
     b = nxt_nact + 1
     if(allocated(nxt_wcnt) .and. nxt_nact>0) then
        u = min(max((log(wlc) - nxt_lw0)/nxt_dlw,-1d0),dble(nxt_nwt+2))
        k = int(u)
        if(u >= 1d0) a = nxt_wcnt(min(k-1,nxt_nwt))
        if(k+2 <= nxt_nwt) b = nxt_wcnt(max(k+2,0)) + 1
     endif
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
pure subroutine nxt_bffreq(lv,t,state,nu,ic)
!-----------------------------------------------------------------------
! Frequency of a recombination photon to global level lv at temperature t:
! energy spectrum j_nu ~ sigma(nu) nu^3 e^{-h nu/kT}, by rejection from
! nu = nu0 + y kT/h, y ~ e^{-y}, with the running maximum of sigma u^3.
! With the cell ic (a recombination of the macro-atom; level dissolution):
! into the dissolved continuum below the edge with its share nxt_fdis.
!-----------------------------------------------------------------------
  integer,intent(in) :: lv
  real*8,intent(in) :: t
  type(rnd_t),intent(inout) :: state
  real*8,intent(out) :: nu
  integer,intent(in),optional :: ic
  real*8 :: r,y,u,f,fmax,nu0,r2
  integer :: i,k,it,j0,j1,jc,a,b,m,ic0
  logical :: ok
  if(present(ic) .and. nxd_on) then
     if(allocated(nxt_fdis)) then
        if(nxt_fdis(lv,ic) > 0.) then
           call rnd_r(r,state)
           if(r < nxt_fdis(lv,ic)) then
              do it=1,1000
                 call rnd_r(r,state)
                 call rnd_r(r2,state)
                 call nxd_sample(lv,nxt_disp(:,ic),t,nxt_dnuc(lv,ic),r,r2,nu,ok)
                 if(ok) return
              enddo
              nu = nxa_lvnu0(lv)
              return
           endif
        endif
     endif
  endif
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
     call nxt_bffreq(idx,1d0/grd_tempinv(ic),state,nu,ic)
     wlc = cc/nu
  case default
     call rnd_r(r,state)
     wlc = cc/(nlteh_twophoton_y(r)*nxt_2pnu(idx))
  end select
end subroutine emit_state
!
!
pure function nxt_emloc(em) result(l)
!-- emission (packet%em) that is part of the source function of the cell's
!-- groups (thermal or bound-free; nxt_ldb)
  integer,intent(in) :: em
  logical :: l
  l = mod(em,100)==2 .or. mod(em,100)==3
end function nxt_emloc
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
  case(em_line,em_lines)
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
pure function bfsig(lv,ig,ic) result(s)
!-- bound-free cross-section of global level lv, mean over group ig in
!-- cell ic: above its edge, and below it its dissolved continuum (level
!-- dissolution, nltedismod)
  integer,intent(in) :: lv,ig,ic
  real*8 :: s,dm,dh,dhs,dsp
  s = 0d0
  if(ig <= nxa_g1(lv)) s = nxa_sgm(nxa_sgo(lv)+ig)
  if(.not.nxd_on) return
  if(nxd_ga(lv) <= 0 .or. ig < nxd_ga(lv) .or. ig > nxd_gb(lv)) return
  if(.not.allocated(nxt_disp)) return
  call nxd_gint(lv,ig,nxt_disp(:,ic),1d4,nxt_dnuc(lv,ic),dm,dh,dhs,dsp)
  s = s + dm
end function bfsig
!
!
pure subroutine nxt_bfabs(ic,ig,nu,state,kind,idx,wlc)
!-----------------------------------------------------------------------
! NLTE bound-free absorption of a packet of comoving frequency nu (group
! ig) in cell ic: choose the level by n_l <sigma_l>_g (nxt_bfpick), then
! activate its continuum state with probability nu_l/nu (the ionization
! energy) or the thermal pool (the photo-electron's kinetic energy).
!-----------------------------------------------------------------------
  integer,intent(in) :: ic,ig
  real*8,intent(in) :: nu
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: kind,idx
  real*8,intent(out) :: wlc
  real*8 :: r
  integer :: lsel
  call nxt_bfpick(ic,ig,state,lsel)
  if(lsel==0) then
     call nxt_kpacket(ic,state,kind,idx,wlc)
     return
  endif
  call rnd_r(r,state)
  if(r < nxa_lvnu0(lsel)/nu) then
     call nxt_activate(ic,nxa_lvel(lsel),nxa_lvc(lsel),state,kind,idx,wlc)
  else
     call nxt_kpacket(ic,state,kind,idx,wlc)
  endif
end subroutine nxt_bfabs
!
!
pure subroutine nxt_bfpick(ic,ig,state,lsel)
!-----------------------------------------------------------------------
! Level of a bound-free absorption in group ig of cell ic, with
! probability n_l <sigma_l>_g over the sum of all levels (0: none absorbs).
! The choice is exact and mostly from the list of the group's block
! (nxt_tables), of weight s at the group: a uniform number on [0, s + R),
! with R the bound of the weight of all other levels, selects from the
! list below s.  Above s the other levels are walked with their weights,
! which fill [s, s + s') with s' <= R; a number beyond s + s' is drawn
! again on [0, s + s').  The other levels are walked with probability
! R/(s + R).
!-----------------------------------------------------------------------
  integer,intent(in) :: ic,ig
  type(rnd_t),intent(inout) :: state
  integer,intent(out) :: lsel
  real*8 :: s,so,r
  integer :: ib,nk
!
  lsel = 0
  if(allocated(nxt_bfl)) then
     ib = (ig-1)/nxt_bsz + 1
     nk = nxt_bfn(ib,ic)
     s = nxt_bfs(ig,ic)
     call rnd_r(r,state)
     r = r*(s + nxt_bfr(ib,ic))
     if(r < s) then
        lsel = from_list(r)
     else
        call others(r - s,lsel,so)
     endif
  else
!-- (no lists: all levels)
     ib = 1
     nk = 0
     s = 0d0
     call others(huge(1d0),lsel,so)
     lsel = 0
  endif
  if(lsel==0) then
     if(s + so <= 0d0) return
     call rnd_r(r,state)
     r = r*(s + so)
     if(r < s) then
        lsel = from_list(r)
     else
        call others(r - s,lsel,so)
        if(lsel==0) lsel = lastpos()
     endif
  endif
contains
  pure real*8 function wgt(lv_)
!-- n_l <sigma_l>_g of level lv_ (0 below its edge but for its dissolved
!-- continuum, bfsig)
    integer,intent(in) :: lv_
    integer :: iel_
    wgt = 0d0
    iel_ = nxa_lvel(lv_)
    wgt = nxt_nel(iel_,ic)*nxt_x(nxa_el(iel_)%xoff+nxa_lvs(lv_),ic)*nxt_w(lv_,ic)*bfsig(lv_,ig,ic)
  end function wgt
  pure integer function from_list(rr)
!-- listed level where the cumulative weight passes rr (< s)
    real*8,intent(in) :: rr
    integer :: k_,l_
    real*8 :: w_,c_
    from_list = 0
    c_ = 0d0
    do k_=1,nk
       l_ = nxt_bfl(k_,ib,ic)
       w_ = wgt(l_)
       if(w_ <= 0d0) cycle
       c_ = c_ + w_
       from_list = l_
       if(c_ > rr) return
    enddo
  end function from_list
  pure subroutine others(rr,lsel_,so_)
!-- other level where their cumulative weight passes rr (0: not reached;
!-- so_ is then their total weight)
    real*8,intent(in) :: rr
    integer,intent(out) :: lsel_
    real*8,intent(out) :: so_
    integer :: l_,iw_
    real*8 :: w_
    lsel_ = 0
    so_ = 0d0
    do l_=1,nxa_nlv
       if(ig > nxa_g1(l_)) cycle
       if(nk>0) then
          iw_ = (l_-1)/64 + 1
          if(btest(nxt_bfm(iw_,ib,ic),mod(l_-1,64))) cycle
       endif
       w_ = wgt(l_)
       if(w_ <= 0d0) cycle
       so_ = so_ + w_
       if(so_ > rr) then
          lsel_ = l_
          return
       endif
    enddo
  end subroutine others
  pure integer function lastpos()
!-- (rounding) the last level with weight
    integer :: l_
    lastpos = 0
    do l_=nxa_nlv,1,-1
       if(wgt(l_) > 0d0) then
          lastpos = l_
          return
       endif
    enddo
  end function lastpos
end subroutine nxt_bfpick
!
!
pure subroutine top_k(a,k,ix)
!-- partial ordering of the index ix (quickselect): the k largest a(ix)
!-- first, in no particular order
  real*8,intent(in) :: a(:)
  integer,intent(in) :: k
  integer,intent(inout) :: ix(:)
  integer :: lo,hi,i,j,t
  real*8 :: p
  lo = 1
  hi = size(ix)
  if(k<1 .or. k>=hi) return
  do while(hi>lo)
     p = a(ix((lo+hi)/2))
     i = lo
     j = hi
     do while(i<=j)
        do while(a(ix(i)) > p)
           i = i + 1
        enddo
        do while(a(ix(j)) < p)
           j = j - 1
        enddo
        if(i<=j) then
           t = ix(i)
           ix(i) = ix(j)
           ix(j) = t
           i = i + 1
           j = j - 1
        endif
     enddo
!-- a(ix(lo:j)) >= p >= a(ix(i:hi)), and = p between
     if(k<=j) then
        hi = j
     elseif(k>=i) then
        lo = i
     else
        exit
     endif
  enddo
end subroutine top_k
!
!
pure subroutine sort_desc(a,ix)
!-- heap sort of the index ix by decreasing a(ix)
  real*8,intent(in) :: a(:)
  integer,intent(inout) :: ix(:)
  integer :: n,i,j,l,ir,t
  n = size(ix)
  if(n<2) return
  l = n/2 + 1
  ir = n
  do
     if(l>1) then
        l = l - 1
        t = ix(l)
     else
        t = ix(ir)
        ix(ir) = ix(1)
        ir = ir - 1
        if(ir==1) then
           ix(1) = t
           exit
        endif
     endif
     i = l
     j = l + l
     do while(j<=ir)
        if(j<ir) then
           if(a(ix(j)) > a(ix(j+1))) j = j + 1
        endif
        if(a(t) > a(ix(j))) then
           ix(i) = ix(j)
           i = j
           j = j + j
        else
           j = ir + 1
        endif
     enddo
     ix(i) = t
  enddo
end subroutine sort_desc
!
!
subroutine nxt_write(ncell,rcell,vcell,tempinv,radtemp,vol)
!-----------------------------------------------------------------------
! Per-cell NLTE diagnostics: output.nltex (ion fractions of each element,
! departure coefficients of the first levels of each ion, line depths of
! selected lines) and, for hydrogen, output.nlteh in the format of the
! H-only treatment, with the estimators of the Halpha incident intensity
! of the last Monte Carlo step (nxt_tally_finish: resonance crossings and
! their number, the window estimate and its effective number of samples,
! the group J without the line's own photons), its emission by the
! macro-atom and absorption at the resonances in the Monte Carlo [erg/s],
! and the absorption of the electron-scattered photons in its core:
! bes J_es and (nxt_esrates) the far-side tally, both beta-weighted.
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: rcell(ncell),vcell(ncell),tempinv(ncell),radtemp(ncell),vol(ncell)
  integer :: ic,iel,j,i,o,k,kla,kha,khb,ih,g
  real*8 :: t,xl(nxa_nst),qv,b(3),fr(nxa_maxion+1),b8(8),trw(2),bp,jd(5),jes
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
!-- local detailed balance: the parts of J in its groups
  if(allocated(nxt_ldbd)) then
     open(4,file='output.nltex_ldb',status='replace')
     write(4,'(a)') '# groups in local detailed balance: 4 pi Int dnu [erg/s/cm^2] of Lambda* S, '// &
          'J_ext and the Monte Carlo estimate, at 228-504 and 504-912 A'
     write(4,'(a1,a11,6a12)') '#','r[cm]','LS_228','Jx_228','JMC_228','LS_504','Jx_504','JMC_504'
     do ic=1,ncell
        write(4,'(1p,7e12.4)') rcell(ic),nxt_ldbd(:,:,ic)
     enddo
     close(4)
  endif
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
  write(4,'(a1,a11,29a12)') '#','r[cm]','v[cm/s]','T[K]','Trad[K]','nH','ne','x_HII','x_1','x_2','x_3', &
       'b_1','b_2','b_3','tau_Lya','tau_Ha','tau_Hb','beta_Ha','Jb_Ha/B','ljb_Ha','bes_Ha', &
       'nres_Ha','Jres_Ha/B','nwin_Ha','Jwin_Ha/B','Jgrp_Ha/B','LHa_emMC','LHa_abMC','besJes/B','Jfar_Ha/B'
  g = nxt_lgrp(kha)
  do ic=1,ncell
     t = 1d0/tempinv(ic)
     call nx_lte(t,nxt_ne(ic),xl)
     do k=1,3
        b(k) = (nxt_x(o+k,ic)/max(nxt_x(o+nxa_el(iel)%ns,ic),1d-300))/ &
             (xl(o+k)/max(xl(o+nxa_el(iel)%ns),1d-300))
     enddo
     qv = max(vcell(ic)/rcell(ic),1d-30)
     bp = max(bplanck(nxa_lnu(kha),t),1d-300)
     jd(1) = nxt_jbcnt(kha,ic)
     jd(2) = nxt_jbsum(kha,ic)/(4d0*pi*vol(ic)*bp)
     jd(3) = 0d0
     if(nxt_jwsq(kha,ic)>0d0) jd(3) = nxt_jwsum(kha,ic)**2/nxt_jwsq(kha,ic)
     jd(4) = nxt_jwsum(kha,ic)/(4d0*pi*vol(ic)*bp)
     jd(5) = 0d0
     if(g>0) jd(5) = max(nxt_jg(g,ic) - nxt_jlosum(g,ic)/(4d0*pi*vol(ic)*nxt_dnug(g)),0d0)/bp
     jes = 0d0
     if(nxt_have_mc) jes = nxt_bes(kha,ic)*nx_jes(kha,t,nxt_ng,nxt_nug,nxt_dnug,nxt_jg(:,ic))/bp
     if(jes < 1d-99) jes = 0d0
     write(4,'(1p,29e12.4)') rcell(ic),vcell(ic),t,radtemp(ic),nxt_nel(iel,ic),nxt_ne(ic), &
          nxt_x(o+nxa_el(iel)%ns,ic),nxt_x(o+1,ic),nxt_x(o+2,ic),nxt_x(o+3,ic),b(1),b(2),b(3), &
          nxt_kq(kla,ic)/qv,nxt_kq(kha,ic)/qv,nxt_kq(khb,ic)/qv,nxt_beta(kha,ic), &
          nxt_jb(kha,ic)/bp,merge(1d0,0d0,nxt_jbcnt(kha,ic)>=nxt_njbmin), &
          merge(nxt_bes(kha,ic),0d0,nxt_bes(kha,ic)>1d-99),jd,nxt_haem(ic),nxt_haab(ic),jes, &
          merge(nxt_efsum(kha,ic)/bp,0d0,nxt_esrates .and. nxt_efsum(kha,ic)/bp > 1d-99)
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
subroutine nxt_ali_cell(ng,kap,fcl,dr,lali)
!-----------------------------------------------------------------------
! Accelerated lambda iteration in the next statistical-equilibrium solves
! (nltexmod): the absorption opacity kap of the cell's groups with which
! the Monte Carlo J was made (per clump volume), and the local operator
! Lambda* of the cell with tau = fcl kap dr (the volume-mean opacity
! across the cell width dr): in the groups in local detailed balance
! (tau >= nxt_ldbtau, nxt_ldb) the Lambda* of their J, in the others
! (with lali, default true) 1 - (1 - e^-tau)/tau, at most 0.99, else 0.
!-----------------------------------------------------------------------
  integer,intent(in) :: ng
  real*4,intent(in) :: kap(ng)
  real*8,intent(in) :: fcl,dr
  logical,intent(in),optional :: lali
  integer :: ig
  real*8 :: tau
  logical :: la
  la = .true.
  if(present(lali)) la = lali
  if(allocated(nx_ali_lst)) then
     if(size(nx_ali_lst)/=ng) deallocate(nx_ali_lst,nx_ali_kap)
  endif
  if(.not.allocated(nx_ali_lst)) allocate(nx_ali_lst(ng),nx_ali_kap(ng))
  do ig=1,ng
     nx_ali_kap(ig) = kap(ig)
     tau = fcl*kap(ig)*dr
     if(nxt_ldbtau>0d0 .and. tau>=nxt_ldbtau) then
        nx_ali_lst(ig) = nxt_lsldb(tau)
        if(nxt_ldbone) nx_ali_lst(ig) = 1d0
     elseif(.not.la) then
        nx_ali_lst(ig) = 0d0
     elseif(tau < 1d-4) then
        nx_ali_lst(ig) = 0.5d0*tau
     else
        nx_ali_lst(ig) = min(1d0 - (1d0 - exp(-tau))/tau,0.99d0)
     endif
  enddo
  nx_ali_on = .true.
end subroutine nxt_ali_cell
!
!
subroutine nxt_ali_off
  nx_ali_on = .false.
end subroutine nxt_ali_off
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
!-- the radiation field of the last iteration (warm start)
  call nxt_write_jmc(ncell,mext)
end subroutine nxt_write_state
!
!
subroutine nxt_write_jmc(ncell,mext)
!-----------------------------------------------------------------------
! output.nltex_jmc: the Monte Carlo estimators of the last iteration that
! the next iteration's tables take (nxt_tables), the J of the groups and
! the incident J of the lines (-1: none), per cell with the exterior mass
! of the cell centre.  A warm start from output.nltex_state reads them
! from this file next to it (nxt_read_init), so that its first iteration
! has the macro-atom of the saved state instead of J = B(T_R).
!-----------------------------------------------------------------------
  integer,intent(in) :: ncell
  real*8,intent(in) :: mext(ncell)
  integer :: ic,ig
  open(4,file='output.nltex_jmc',status='replace')
  write(4,'(a)') '# SuperLite Monte Carlo radiation field per cell (warm start: in_nltex_init)'
  write(4,'(a,3i9)') 'ng nlin ncell',nxt_ng,nxa_nlin,ncell
  write(4,'(a)') '# group boundaries [Hz] (decreasing)'
  write(4,'(1p,6es22.14e3)') nxt_nug(1)+0.5d0*nxt_dnug(1),(nxt_nug(ig)-0.5d0*nxt_dnug(ig),ig=1,nxt_ng)
  write(4,'(a)') '# per cell: exterior mass of the cell centre [g]; J of the groups '// &
       '[erg/s/cm^2/Hz/sr]; incident J of the lines (-1: none)'
  do ic=1,ncell
     write(4,'(1p,es16.7e3)') mext(ic)
     write(4,'(1p,10es13.4e3)') max(nxt_jg(:,ic),0d0)
     write(4,'(1p,10es13.4e3)') merge(max(nxt_jb(:,ic),0d0),-1d0,nxt_ljb(:,ic))
  enddo
  close(4)
end subroutine nxt_write_jmc
!
!
subroutine nxt_read_prev(fname,ncell,mass,dt,lprint)
!-----------------------------------------------------------------------
! Populations of the previous model (its output.nltex_state) at the
! exterior mass coordinate of each cell centre (nxt_read_state), for the
! time-dependent statistical equilibrium with dt [s] since that model.
! The exterior mass is Lagrangian when both models end at the same mass
! shell (the outer edge of the CSM); the model atoms must be the same.
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname
  integer,intent(in) :: ncell
  real*8,intent(in) :: mass(ncell),dt
  logical,intent(in) :: lprint
  real*8 :: t(ncell),ne(ncell)
!
  if(dt<=0d0) stop 'nxt_read_prev: in_nltex_dtprev must be > 0'
  if(allocated(nxt_xprev)) deallocate(nxt_xprev)
  allocate(nxt_xprev(nxa_nst,ncell))
  call nxt_read_state(fname,ncell,mass,nxt_xprev,t,ne,lprint,'time-dependent: previous model')
  nxt_dtinv = 1d0/dt
  if(lprint) write(6,'(1x,a,es10.3,a)') 'NLTE time-dependent: dt =',dt,' s'
end subroutine nxt_read_prev
!
!
subroutine nxt_read_init(fname,ncell,mass,lprint)
!-----------------------------------------------------------------------
! Warm start: the populations, T_e and n_e of a saved state (output.
! nltex_state of a run on the same structure and model atoms) at the
! exterior mass coordinate of each cell centre (nxt_read_state).  The
! first iteration keeps them (nxt_gas_solve) and the transport starts
! from that state, with its radiation field (nxt_read_jmc), so a rerun
! converges in a few iterations.
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname
  integer,intent(in) :: ncell
  real*8,intent(in) :: mass(ncell)
  logical,intent(in) :: lprint
  integer :: k
  if(allocated(nxt_xinit)) deallocate(nxt_xinit,nxt_tinit,nxt_neinit)
  allocate(nxt_xinit(nxa_nst,ncell),nxt_tinit(ncell),nxt_neinit(ncell))
  call nxt_read_state(fname,ncell,mass,nxt_xinit,nxt_tinit,nxt_neinit,lprint,'warm start from')
  nxt_warm = .true.
!-- the radiation field of that state (output.nltex_jmc next to it)
  k = index(fname,'nltex_state',back=.true.)
  if(k>0) call nxt_read_jmc(fname(:k-1)//'nltex_jmc'//fname(k+11:),ncell,mass,lprint)
end subroutine nxt_read_init
!
!
subroutine nxt_read_jmc(fname,ncell,mass,lprint)
!-----------------------------------------------------------------------
! Warm start: the Monte Carlo estimators of the saved state (output.
! nltex_jmc, nxt_write_jmc) at the exterior mass coordinate of each cell
! centre, linear in it between that model's cell centres as the
! populations (nxt_read_state), the J of this model's groups averaged over
! the overlap with that model's groups.  The first iteration then makes
! the tables of the saved state (nxt_have_mc); without the file it makes
! them with J = B(T_R) and jb = 0, and its Monte Carlo J lacks most of the
! re-emission in the lines.
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname
  integer,intent(in) :: ncell
  real*8,intent(in) :: mass(ncell)
  logical,intent(in) :: lprint
  integer :: iu,ios,ngp,nl,npc,ic,ig,k,kk,a,b,m
  character(256) :: line
  real*8,allocatable :: nub(:),mp(:),jgp(:,:),jbp(:,:),jgn(:,:)
  real*8 :: mext(ncell),w,nlo,nhi,s,o
!
  open(newunit=iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) then
     if(lprint) write(6,'(1x,a,1x,a,a)') 'NLTE warm start: no',trim(fname), &
          '; the first iteration makes its tables with J = B(T_R)'
     return
  endif
  call skipc
  read(iu,'(a)') line
  if(index(line,'ng nlin ncell')/=1) stop 'nxt_read_jmc: not an output.nltex_jmc file'
  read(line(14:),*) ngp,nl,npc
  if(nl/=nxa_nlin) stop 'nxt_read_jmc: other model atoms than this model'
  allocate(nub(ngp+1),mp(npc),jgp(ngp,npc),jbp(nl,npc),jgn(nxt_ng,npc))
  call skipc
  read(iu,*) nub
  do ic=1,npc
     call skipc
     read(iu,*) mp(ic)
     read(iu,*) jgp(:,ic)
     read(iu,*) jbp(:,ic)
  enddo
  close(iu)
!-- this model's groups (decreasing frequency, as that model's)
  k = 1
  do ig=1,nxt_ng
     nhi = nxt_nug(ig) + 0.5d0*nxt_dnug(ig)
     nlo = nxt_nug(ig) - 0.5d0*nxt_dnug(ig)
     do while(k<ngp)
        if(nub(k+1) < nhi) exit
        k = k + 1
     enddo
     jgn(ig,:) = 0d0
     s = 0d0
     do kk=k,ngp
        if(nub(kk) <= nlo) exit
        o = min(nhi,nub(kk)) - max(nlo,nub(kk+1))
        if(o <= 0d0) cycle
        jgn(ig,:) = jgn(ig,:) + o*jgp(kk,:)
        s = s + o
     enddo
     if(s > 0d0) then
        jgn(ig,:) = jgn(ig,:)/s
     elseif(nlo >= nub(1)) then
        jgn(ig,:) = jgp(1,:)
     else
        jgn(ig,:) = jgp(ngp,:)
     endif
  enddo
!-- exterior mass of the cell centres of this model (mp decreases outward)
  mext(ncell) = 0.5d0*mass(ncell)
  do ic=ncell-1,1,-1
     mext(ic) = mext(ic+1) + 0.5d0*(mass(ic+1) + mass(ic))
  enddo
  do ic=1,ncell
     if(mext(ic) >= mp(1)) then
        a = 1; b = 1; w = 0d0
     elseif(mext(ic) <= mp(npc)) then
        a = npc; b = npc; w = 0d0
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
     endif
     nxt_jg(:,ic) = (1d0-w)*jgn(:,a) + w*jgn(:,b)
     do k=1,nxa_nlin
        if(jbp(k,a)>=0d0 .and. jbp(k,b)>=0d0) then
           nxt_jb(k,ic) = (1d0-w)*jbp(k,a) + w*jbp(k,b)
        else
           nxt_jb(k,ic) = max(jbp(k,a),jbp(k,b))
        endif
        nxt_ljb(k,ic) = nxt_jb(k,ic) >= 0d0
        nxt_jb(k,ic) = max(nxt_jb(k,ic),0d0)
     enddo
  enddo
  nxt_have_mc = .true.
  if(lprint) write(6,'(1x,a,1x,a,a,i5,a,i6,a)') 'NLTE warm start: Monte Carlo radiation field from', &
       trim(fname),':',npc,' cells,',ngp,' groups'
  deallocate(nub,mp,jgp,jbp,jgn)
contains
  subroutine skipc
!-- to the next record that is not a comment
    do
       read(iu,'(a)',iostat=ios) line
       if(ios/=0) stop 'nxt_read_jmc: file truncated'
       if(line(1:1)/='#' .and. len_trim(line)>0) exit
    enddo
    backspace(iu)
  end subroutine skipc
end subroutine nxt_read_jmc
!
!
subroutine nxt_read_state(fname,ncell,mass,x,t,ne,lprint,what)
!-----------------------------------------------------------------------
! A model's output.nltex_state at the exterior mass coordinate of each
! cell centre of this model: the state fractions, T_e and n_e, linear in
! that coordinate between that model's cell centres (the end values
! beyond); each element's fractions renormalized.  The model atoms must
! be the same.
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname,what
  integer,intent(in) :: ncell
  real*8,intent(in) :: mass(ncell)
  real*8,intent(out) :: x(nxa_nst,ncell),t(ncell),ne(ncell)
  logical,intent(in) :: lprint
  integer :: iu,ios,nel,nst,npc,ic,iel,a,b,m,k,o
  integer :: zs(200)
  character(4096) :: line
  real*8,allocatable :: mp(:),xp(:,:),tp(:),np(:)
  real*8 :: mext(ncell),rd(5),w
!
  open(newunit=iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) stop 'nxt_read_state: cannot open the state file'
  call nextrec
  if(index(line,'nel nst ncell')/=1) stop 'nxt_read_state: not an output.nltex_state file'
  read(line(14:),*) nel,nst,npc
  if(nel/=nxa_nel .or. nst/=nxa_nst) stop 'nxt_read_state: other model atoms than this model'
  call nextrec
  if(nel>100) stop "nxt_read_state: too many elements"
  read(line(14:),*) zs(:2*nel)
  do iel=1,nel
     if(zs(2*iel-1)/=nxa_el(iel)%z .or. zs(2*iel)/=nxa_el(iel)%ns) &
          stop 'nxt_read_state: other model atoms than this model'
  enddo
  allocate(mp(npc),xp(nst,npc),tp(npc),np(npc))
  do ic=1,npc
     call nextrec
     read(line,*) rd
     read(iu,*) xp(:,ic)
     mp(ic) = rd(1)
     tp(ic) = rd(4)
     np(ic) = rd(5)
  enddo
  close(iu)
!-- exterior mass of the cell centres of this model
  mext(ncell) = 0.5d0*mass(ncell)
  do ic=ncell-1,1,-1
     mext(ic) = mext(ic+1) + 0.5d0*(mass(ic+1) + mass(ic))
  enddo
!-- mp decreases outward (index npc: the outermost cell)
  do ic=1,ncell
     if(mext(ic) >= mp(1)) then
        x(:,ic) = xp(:,1); t(ic) = tp(1); ne(ic) = np(1)
     elseif(mext(ic) <= mp(npc)) then
        x(:,ic) = xp(:,npc); t(ic) = tp(npc); ne(ic) = np(npc)
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
        x(:,ic) = (1d0-w)*xp(:,a) + w*xp(:,b)
        t(ic) = (1d0-w)*tp(a) + w*tp(b)
        ne(ic) = (1d0-w)*np(a) + w*np(b)
     endif
     do iel=1,nxa_nel
        o = nxa_el(iel)%xoff
        k = nxa_el(iel)%ns
        w = sum(x(o+1:o+k,ic))
        if(w>0d0) x(o+1:o+k,ic) = x(o+1:o+k,ic)/w
     enddo
  enddo
  if(lprint) write(6,'(1x,a,1x,a,1x,a,a,i5,a,es10.3,a,es10.3,a)') 'NLTE',what,trim(fname),':',npc, &
       ' cells, exterior mass',mp(1),' g (this model',mext(1),' g)'
  deallocate(mp,xp,tp,np)
contains
  subroutine nextrec
!-- next record that is not a comment
    do
       read(iu,'(a)',iostat=ios) line
       if(ios/=0) stop 'nxt_read_state: file truncated'
       if(line(1:1)/='#' .and. len_trim(line)>0) exit
    enddo
  end subroutine nextrec
end subroutine nxt_read_state
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
