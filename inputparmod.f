*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
*See LANL_COPYING and LANL_README for details of LANL copyright assertion.
      module inputparmod
c     ------------------
      implicit none
************************************************************************
* input parameters
* New parameters need to be added in this routine in four places:
*  1) variable declaration
*  2) namelist
*  3) pointer array
*  4) sanity conditions
************************************************************************
c-- model id and comments
      character(40) :: in_name = "spn"  !simulation name/title, for post-processing identification
      character(80) :: in_comment = ""  !comment if any
c-- parallelization
      integer :: in_nomp = 1  !number of openmp threads
c
c-- experimental iterative approach
      integer :: in_niter = 1 !number of iterations
c
c-- grid geometry and dimensions
      integer :: in_grd_igeom = 11 !only 11=1Dsph is accepted for now
      integer :: in_ndim(3) = [1, 1, 1]  ![nx,1,1] nx: number of x-direction cells
c
c
c
c-- group structure
      integer :: in_grp_ng = -1      !number of groups: 0 read wl-grid from file
      integer :: in_grp_ngs = 1      !number of subgroups per opacity group
      real*8 :: in_grp_wlmin =   100d-8  !lower wavelength boundary [cm]
      real*8 :: in_grp_wlmax = 32000d-8  !upper wavelength boundary [cm]
c
c
c-- outbound flux group and direction bins
      integer :: in_flx_ndim(3) = [0, 1, 1]
      real*8 :: in_flx_wlmin =  1000d-8  !lower wavelength flux boundary [cm]
      real*8 :: in_flx_wlmax = 32000d-8  !upper wavelength flux boundary [cm]
c
c-- interior source
      integer :: in_src_ns = 0     !number of source particles generated per time step (total over all ranks)
      integer :: in_src_n2s = -1   !2^n source particles generated per time step (total over all ranks)
c-- the first in_niter_early iterations with 2^in_src_n2s_early source
c-- particles (fast convergence of the gas state before the final ones)
      integer :: in_niter_early = 0
      integer :: in_src_n2s_early = -1
      real*8 :: in_srcepwr = 1d0  !source particle number-energy slope, 1 is linear, equal number of packets per erg.
      real*8 :: in_L_bol = 1d0 ![erg/s] !bolometric luminosity
      real*8 :: in_L_inner = 0d0 ![erg/s] !experimental inner boundary source
c-- boundary source
      character(4) :: in_surfsrcmu = 'isot' !isot|beam: surface source direction distribution
      real*8 :: in_srcmax = 0d0 !peak source strength
c
c-- transport
      logical :: in_trn_errorfatal = .true. !stop on transport error, disable for production runs
      real*8 :: in_trn_tauddmc = 5d0  !number of mean free paths per cell required for DDMC
      real*8 :: in_taulump = 10d0 !number of of mean free paths needed to lump DDMC groups
      logical :: in_puretran = .false. !use IMC only instead of IMC+DDMC hybrid
      logical :: in_trn_noamp = .true.  !disable amplification factor
      logical :: in_thm_aniso = .false. !experimental anisotropic scattering
c-- packets crossing the inner boundary inward: 'refl' re-enter outward at
c-- the same wavelength (legacy), 'thrm' are re-emitted with the boundary
c-- blackbody (T = in_srcmax or the inner cell's T_rad), 'absb' are absorbed
      character(4) :: in_innerbc = 'refl'
c-- the luminosity of the inner boundary ('lumg') injected with the
c-- spectrum of the diffusion approximation, as CMFGEN's inner boundary:
c-- dB/dT over the extinction of the innermost cell (Rosseland weighting),
c-- instead of the Planck spectrum, which has (x/4)e^x/(e^x - 1) times less
c-- flux for a grey opacity (x = h nu/kT: 2.1 times less at 1000 A and
c-- 17 000 K, 2 times more at 6000 A)
      logical :: in_srcdiff = .false.
      logical :: in_thm_thermal = .false. !Thomson scattering off thermal (Maxwellian) electrons:
!                                          line broadening and Comptonization by electron scattering
      logical :: in_trn_doplegacy = .false. !legacy Doppler shifts: divergence-based distance and
!                                            lab wavelength reset at each group crossing
c
c-- debugging
      logical :: in_noeos = .false.     !don't use the EOS
c
c-- physical opacities
      real*8 :: in_opacmixrossel = 0d0 !mix rosseland with planck average, 1=pure rosseland
c-- physical opacities debuging
      logical :: in_nobbopac = .false.  !turn off bound-bound opacity
      logical :: in_nobfopac = .false.  !turn off bound-bound opacity
      logical :: in_noffopac = .false.  !turn off bound-bound opacity
      logical :: in_nothmson = .false.  !turn off thomson scattering
c
c-- experimental nebular approximation
      logical :: in_nebular = .false. ! turn on nebular approximation
c-- photospheric radius for dilution factor
      real*8 :: in_R_phot = -1 ![cm] photospheric radius where Rosseland mean optical depth = 1
c
c-- NLTE treatment for hydrogen b-b opacity calculations
      logical :: in_nlte = .false. ! turn on NLTE
      integer :: in_nlte_nelem = 1 ! number of elements treated with nLTE ! for now keep =1 for h
c-- NLTE hydrogen with Sobolev line resonances and a macro-atom (requires in_nlte)
      logical :: in_nlteh = .false.
      integer :: in_nlteh_nlev = 30 !bound levels n = 1..in_nlteh_nlev
      real*8 :: in_nlteh_damp = 0.5d0 !weight of the new populations per iteration
c-- benchmark: in_nlteh_tau0 > 0 makes H-alpha a pure-scattering Sobolev
c-- line with tau = in_nlteh_tau0 (for v/r = dv/dr) and no other H lines
      real*8 :: in_nlteh_tau0 = 0d0
c-- NLTE ions from model-atom files (Tools/cmfgen2superlite_atom.py),
c-- names separated by blanks; with in_nlteh hydrogen is built in.
c-- All NLTE elements share the macro-atom and statistical equilibrium.
      character(1024) :: in_nltex_atoms = ' '
c-- directory of the model-atom files (prefixed to relative names)
      character(256) :: in_nltex_dir = ' '
c-- charge-exchange reactions among the NLTE elements
c-- (Tools/cmfgen2superlite_chgexch.py; relative to in_nltex_dir; 'none':
c-- no charge exchange)
      character(256) :: in_nltex_chgexch = 'chgexch.dat'
c-- NLTE treatment active (in_nlteh or model atoms; set, not read)
      logical :: in_nltex = .false.
c-- interior source: 'surf' blackbody surfaces at T_rad (legacy),
c-- 'lumg' energy injected where the luminosity L(r) of the input ('lum'
c-- column, at the cell's inner edge) grows outward, emitted isotropically
c-- with the local emissivity, plus L at the inner boundary
      character(4) :: in_srctype = 'surf'
c-- shocks (velocity dropping outward by at least in_shock_dvmin): the
c-- injected energy of shock cells is emitted as a radiatively cooling
c-- flow, T uniform between the gas temperature and the post-shock
c-- temperature (capped at in_shock_tmax), bremsstrahlung shape at each T
c-- (or the CIE spectrum, in_shock_cie)
      real*8 :: in_shock_dvmin = 3d7 ![cm/s]
      real*8 :: in_shock_tmax = 1d7 ![K]
c-- CIE emissivity table (Tools/chianti_cie_spectra.py): the spectrum at
c-- each T of the cooling flow is the CIE spectrum of the cell's
c-- composition instead of bremsstrahlung ('none': bremsstrahlung)
      character(256) :: in_shock_cie = 'none'
c-- time-dependent statistical equilibrium of the NLTE elements (lumg
c-- source): output.nltex_state of the previous model of the sequence
c-- ('none': steady state), matched by the exterior mass coordinate, and
c-- the time since that model [s]
      character(256) :: in_nltex_prev = 'none'
      real*8 :: in_nltex_dtprev = 0d0
c-- warm start of the NLTE elements: output.nltex_state of a run on the
c-- same structure and model atoms ('none': start from J = B(T_R)); its
c-- populations (and T_e with in_tesolve) are kept in the first iteration
      character(256) :: in_nltex_init = 'none'
c-- (diagnostic) the populations of in_nltex_init in every iteration, no
c-- statistical equilibrium: the transport and the formal integral on
c-- given populations, e.g. those of a CMFGEN model
c-- (Tools/cmfgen2superlite_state.py)
      logical :: in_nltex_hold = .false.
c-- lines of the line list of NLTE ions whose upper level lies above the
c-- model atom, and of the top ions, in the Sobolev expansion opacity
c-- (lower level from the highest model level below it, Boltzmann above);
c-- f: only the lines of the model atoms (as a CMFGEN model with them)
      logical :: in_nltex_bbext = .true.
c-- continuum destruction of the trapped photons of the NLTE lines
c-- (Hummer & Rybicki 1985): the absorption continuum of the line's group
c-- takes a photon before the line takes it back with the probability
c-- P(kappa_c sqrt(pi) v_D/K), in the rates, the macro-atom and the source
c-- function of the groups in local detailed balance (nltexmod nx_pcont)
      logical :: in_nltex_cdest = .true.
c-- incident intensity of the Sobolev NLTE lines in the rates: the terms of
c-- the resonance estimator weighted by the escape probability of their
c-- direction over the mean one, beta(mu)/beta_bar, so that beta_bar J_b
c-- is Sobolev's <beta(mu) I_inc>; f: the plain mean <I_inc> (as before),
c-- which takes the radiation of a core along the directions of small
c-- |Q(mu)| (a wind of nearly constant velocity) for isotropic
c-- (nltextrmod nxt_jbw).  Homologous flows: the same.  Off by default:
c-- exact for an isolated pure-scattering Sobolev line, but in the CMFGEN
c-- H/He benchmark wind it puts b1 1.4-4.6 times CMFGEN's comoving-frame
c-- value against 0.86-1.5 times without it (CSM_PHYSICS.md 4.4).
      logical :: in_nltex_jbbeta = .false.
c-- adaptive damping of the populations between iterations: in each cell
c-- the weight of the new solution 1/(1 - g), g the slope of the solve's
c-- output against its input along the last step (secant), in [0.1, 1]
c-- (nltextrmod dmp_weight); f: the fixed weight in_nlteh_damp
      logical :: in_nltex_dmpsec = .false.
c-- weight of the previous iterations in the Monte Carlo group J of the
c-- NLTE rates (exponential average against the noise of thin cells; 0:
c-- none)
      real*8 :: in_nltex_javg = 0d0
c-- the transport steps at the Sobolev resonances of the lines whose depth
c-- exceeds this in some cell (active lines); the others are left out
      real*8 :: in_nltex_taumin = 1d-4
c-- groups without Monte Carlo flights in a cell take the mean J over
c-- this many groups on each side
c-- (0: none, J = 0 where no packet flew)
      integer :: in_nltex_jfill = 0
c-- emission of the NLTE macro-atoms and k-packets biased into the band
c-- [in_nltex_fwl1, in_nltex_fwl2] (A) by this factor, the packets
c-- weighted (<= 1: off); for a far UV of a thin wind no packet reaches
      real*8 :: in_nltex_fbias = 0d0
c-- the groups in local detailed balance with Lambda* = 1 (DDMC) take J
c-- from the bound-free emission and opacity of the current populations of
c-- all NLTE elements in every iteration of the statistical equilibrium
c-- (instead of the per-level ALI around the J of the last transport)
      logical :: in_nltex_ldbse = .false.
c-- X-ray inner-shell photoionization of the NLTE ions followed by an Auger
c-- electron in the rates: an ion with at least 4 electrons up two stages,
c-- a Li-like one up one (the core shells of the X-ray opacity; off with
c-- in_nobfopac)
      logical :: in_nltex_auger = .true.
      real*8 :: in_nltex_fwl1 = 912d0
      real*8 :: in_nltex_fwl2 = 1064d0
c-- incident intensity of the Sobolev NLTE lines in the statistical
c-- equilibrium: path-length estimator of J_nu in a window this wide
c-- [km/s] next to each line, on the side the comoving frequency of a
c-- flight drifts in from (up to the next line), where it has 30 effective
c-- samples (0: off; the resonance crossings where a line has 30, else the
c-- J of its group without the line's own photons).  It counts the line
c-- photons that drift back into the line: scattered by electrons to the
c-- incident side, or in a decelerating CSM (nltextrmod, nxt_window_tally)
      real*8 :: in_nltex_jbwin = 0d0
c-- level dissolution by the plasma microfield (Hummer & Mihalas 1988
c-- occupation probabilities as CMFGEN's DO_DIS, nltedismod): the LTE
c-- populations carry w, and the bound-free continua go on below the edges
      logical :: in_nltex_dis = .false.
c-- accelerated lambda iteration of the NLTE bound-free transitions (local
c-- operator of each cell and group; nltexmod)
      logical :: in_nltex_ali = .false.
c-- local detailed balance of the NLTE elements in the groups whose
c-- absorption depth across a cell is at least this (0: off): J = Lambda* S
c-- + the radiation from elsewhere instead of the Monte Carlo estimate
c-- (nltextrmod, nxt_ldb).  Default (< 0): 30 when these groups are DDMC
c-- groups (Lambda* = 1: not in_puretran, in_trn_tauddmc <= 30), else off
      real*8 :: in_nltex_ldb = -1d0
c-- electron scattering of NLTE line photons trapped in their Sobolev
c-- resonance zones (nltexmod, nx_bes): microturbulent velocity [km/s] of
c-- the line profiles at v = 0 and at the outer edge (linear in v; < 0:
c-- the value at v = 0), added in quadrature to the thermal velocity.  < 0:
c-- off (point resonances).  Default: in the formal integral, from the
c-- populations of the iterations (as CMFGEN's observer spectrum with
c-- incoherent scattering after iterations with coherent scattering).
c-- in_nltex_esrates (experimental): also in the statistical equilibrium,
c-- the macro-atom and the transport (escape probability beta + bes, and
c-- the photons scattered into the line cores)
      real*8 :: in_nltex_vturb = -1d0, in_nltex_vturbmax = -1d0
      logical :: in_nltex_esrates = .false.
c-- (the escape of the scattered photons through the Lorentz wings of the
c-- lines, radiative damping; f: Doppler profiles only, as CMFGEN's
c-- iterations)
      logical :: in_nltex_esdamp = .true.
c-- Doppler-profile NLTE lines in the Monte Carlo transport of the cells
c-- whose velocity changes by less than in_nltex_dopmc Doppler widths
c-- (static profiles within 4 widths of a packet's comoving frequency:
c-- the line overlap of a slow inner wind); 0: Sobolev resonances.  In
c-- homologous ejecta no cell is that slow
      real*8 :: in_nltex_dopmc = 2d0
c-- (only the lines of the elements from Z = in_nltex_dopz: the H and He
c-- lines keep their Sobolev resonances, whose point absorption and
c-- macro-atom escape are consistent with each other in a dense wind)
      integer :: in_nltex_dopz = 3
c-- (the static lines escape their trapping also by the electron
c-- scattering of the trapped photons and absorb the photons scattered
c-- into their cores, in the rates and the macro-atom: the transport
c-- scatters photons into the static cores explicitly; needs
c-- in_nltex_vturb >= 0 and in_thm_thermal)
      logical :: in_nltex_dopes = .true.
c-- overlapping NLTE lines: a photon trapped in the resonance zone of a
c-- line is absorbed also by the lines whose Doppler profiles (thermal
c-- and in_nltex_vturb) overlap its own, within this separation [km/s]
c-- (0: off).  Sobolev zone of the overlapping lines in the rates
c-- (escape, absorption by the partners, the partners' photons) and in
c-- the macro-atom (net photons to the partners activate their macro-
c-- atoms): O I 1025.76 in Lyman beta (Bowen fluorescence; nltexmod)
      real*8 :: in_nltex_overlap = 0d0
c-- (in_nltex_overlap) partners only for the lines whose lower level is
c-- among the first this many levels of their ion (0: all lines)
      integer :: in_nltex_ovllmax = 0
c-- formal integral of the emergent spectrum after the last iteration
c-- (formalmod): range [cm], spacing [km/s], rays through the inner
c-- boundary
      logical :: in_fi = .false.
      real*8 :: in_fi_wl1 = 1000d-8, in_fi_wl2 = 10000d-8
      real*8 :: in_fi_dv = 20d0
      integer :: in_fi_npc = 20
c-- (the Monte Carlo J of the last in_fi_navg iterations, averaged: less
c-- noise in the electron-scattering emissivity of a converged model)
      integer :: in_fi_navg = 1
c-- (NLTE lines with Doppler profiles, thermal and the microturbulence
c-- in_nltex_vturb..in_nltex_vturbmax, as CMFGEN's cmf_flux; f: Sobolev
c-- resonances, no width)
      logical :: in_fi_doppler = .true.
c-- (only the formal integral, with the populations and the Monte Carlo J
c-- of the saved state of in_nltex_init and its output.nltex_jmc; no
c-- transport: as CMFGEN's cmf_flux after a model)
      logical :: in_fi_only = .false.
c-- (passes of the electron scattering with the J of the formal solution,
c-- from the Monte Carlo J, as CMFGEN's cmf_flux: line photons trapped in
c-- their resonance zones and their scattering through the flow, no Monte
c-- Carlo noise in the scattering emissivity; 0: the Monte Carlo J)
      integer :: in_fi_nes = 0
c-- (spacing [km/s] of the comoving J and the rays of the passes but the
c-- last two, which are on the grid of the spectrum: the scattering
c-- emissivity is smooth over the electron kernel; the wings of strong
c-- lines then within ~5 %, for exploration; 0: the spacing of the
c-- spectrum)
      real*8 :: in_fi_nesdv = 0d0
c-- the J of the groups between in_nltex_fijwl1 and
c-- in_nltex_fijwl2 [A] in the rates from the formal solution, every
c-- iteration after the transport: in_nltex_fijnes electron-scattering
c-- passes on a comoving grid of in_nltex_fijdv [km/s], with the source
c-- functions of the iteration (CMFGEN's J where the Monte Carlo is
c-- starved: the X-rays and the He II continuum of a thin wind).  The
c-- groups in local detailed balance keep Lambda* S + J_ext.
c-- 1: on, 0: off; default (< 0): on with an X-ray source
c-- (in_xray_table), where the Monte Carlo has almost no X-ray packets
      integer :: in_nltex_fij = -1
      real*8 :: in_nltex_fijwl1 = 5d0, in_nltex_fijwl2 = 912d0
      real*8 :: in_nltex_fijdv = 300d0
      integer :: in_nltex_fijnes = 3
c-- X-ray emission of wind shocks as CMFGEN's INC_XRAYS (xraysrcmod):
c-- plasma table in CMFGEN's RS_XRAY_FLUXES format ('none': off); two
c-- components: filling factor, shock temperature [K], onset velocity
c-- [km/s]; smooth-wind densities (CMFGEN's X_SM_WIND)
      character(256) :: in_xray_table = 'none'
      real*8 :: in_xray_ff1 = 0d0, in_xray_ff2 = 0d0
      real*8 :: in_xray_t1 = 0d0, in_xray_t2 = 0d0
      real*8 :: in_xray_v1 = 0d0, in_xray_v2 = 0d0
      logical :: in_xray_smooth = .true.
c-- (share of the packets of a cell with X-ray emission that start as
c-- X-rays, at least: weighted, the X-rays p/q and the others
c-- (1-p)/(1-q) with p their share of the energy, so that the X-ray J of a
c-- wind whose X-rays are 1e-6 of its luminosity is sampled; 0: q = p)
      real*8 :: in_xray_q = 0.05d0
c-- bound-bound opacity of the LTE elements: in_bbsobolev = t uses the
c-- Sobolev expansion opacity with the local velocity gradient (saturating
c-- thick lines) instead of line strengths smeared over the group;
c-- in_bbeps is the share of line interactions that thermalize the photon
c-- (1: pure absorption, the released code; <0: per line from collisional
c-- de-excitation (van Regemorter) against Sobolev escape), the rest is
c-- coherent scattering in the comoving frame
      logical :: in_bbsobolev = .false.
      real*8 :: in_bbeps = 1d0
c-- Rayleigh scattering by H I in the scattering part of the group opacity:
c-- 'lyman' with the resonances of the Lyman lines (rayleighmod, as CMFGEN's
c-- RAYLEIGH_SCAT; the Lorentz wings of Lyman alpha), 'kurucz' the
c-- far-from-resonance fit with the frequency capped at Lyman alpha
c-- (Kurucz's HRAYOP), 'none'; 'lymes' (benchmark) the Lyman resonances
c-- with the Rayleigh-scattered photons redistributed by the thermal
c-- electrons, as CMFGEN's iterations with incoherent electron scattering
c-- (COH_ES = F) treat its Rayleigh opacity
      character(6) :: in_rayleigh = 'kurucz'
c-- electron temperature from the thermal balance of the gas with the
c-- Monte Carlo radiation field (requires in_nlteh), mixed with the old one
c-- with weight in_te_damp per iteration; the LTE elements are ionized and
c-- excited at the input temperature (in_tion = 'hydr'), at the electron
c-- temperature ('elec'), or in the modified nebular approximation with
c-- the temperature and dilution of the Monte Carlo radiation field ('nebu')
      logical :: in_tesolve = .false.
      real*8 :: in_te_damp = 0.5d0
c-- NLTE elements re-solved at each trial temperature of the thermal balance:
c-- Z <= in_te_zsolve (99: all; 2: H and He, the others held)
      integer :: in_te_zsolve = 99
c-- T_e step per cell from the secant of the last two iterations (root
c-- against old T; tebalmod te_stepsec) instead of the sign-reversal
c-- damping (te_step)
      logical :: in_te_secant = .false.
      character(4) :: in_tion = 'hydr'
c-- photons above in_nt_emin [eV] absorbed by the gas make fast
c-- photoelectrons that ionize and excite H (Shull & van Steenberg 1985;
c-- requires in_nlteh); 0: all absorbed energy heats the gas
      real*8 :: in_nt_emin = 0d0
c-- the energy of fast electrons (in_nt_emin, in_nt_lumg) shared among
c-- heating, ionization and excitation of all NLTE species by the
c-- Spencer-Fano degradation of electrons injected at in_nt_sfemax [eV] on
c-- a grid of in_nt_sfde [eV] (nltentmod); f: Shull & van Steenberg's
c-- shares for H I and He I
      logical :: in_nt_sf = .false.
      real*8 :: in_nt_sfemax = 1000d0
      real*8 :: in_nt_sfde = 1d0
c-- share of the 'lumg' energy injected in a cell that is deposited by
c-- fast electrons (gamma-rays and positrons of radioactive decay; needs
c-- in_nt_sf), the rest is thermal
      real*8 :: in_nt_lumg = 0d0
c-- share of the 'lumg' energy of a cell that is stored radiation released
c-- by diffusion (a supernova's L_bol above its deposition), not deposited
c-- in the gas: emitted with a Planck spectrum at the cell's radiation
c-- temperature and left out of the heating of the thermal balance
      real*8 :: in_lumg_rad = 0d0
c-- free-free opacity with each element as a fully ionized atom of
c-- charge Z (released code) instead of its ions from the EOS
      logical :: in_fflegacy = .false.
c
c-- analytic opacities
      character(4) :: in_opacanaltype = 'none'    !none|grey|mono|line: group opacity structure type
c-- line specific group structure
      real*8 :: in_ldisp1 = 1d0  !loosely speaking, the analytic odd group line strength
      real*8 :: in_ldisp2 = 1d0  !loosely speaking, the analytic even group line strength
c-- scattering terms:
      real*8 :: in_gas_sigcoef = 0d0  !power law absorption opacity coefficient
      real*8 :: in_gas_sigtpwr = 0d0  !power law absorption opacity temperature exponent
      real*8 :: in_gas_sigrpwr = 0d0  !power law absorption opacity density exponent
c-- absorption terms:
      real*8 :: in_gas_capcoef = 0d0  !power law absorption opacity coefficient
      real*8 :: in_gas_captpwr = 0d0  !power law absorption opacity temperature exponent
      real*8 :: in_gas_caprpwr = 0d0  !power law absorption opacity density exponent
c
c
c-- output
      logical :: in_io_grabstdout = .false.  !write stdout to file
      logical :: in_io_dogrdtally = .false. !write transport tallies per grid cell
      logical :: in_io_nogriddump = .false.  !don't write grid cell variables
      character(4) :: in_io_opacdump = 'off '    !off|one|all: write opacity data to file
      character(4) :: in_io_pdensdump = 'off '   !off|on: write partial densities to file
      logical :: in_io_profdump = .false.  !write profile to file at the end
      logical :: in_io_debugging = .false. !sets nogriddump to false, opacdump,pdensdump,dogrdtally,profdump to true
c
c-- runtime parameter namelist
      namelist /inputpars/
     & in_name,in_comment,
     & in_nomp,
     & in_niter,
!grd
     & in_grd_igeom,in_ndim,
!grp
     & in_grp_ng,in_grp_ngs,in_grp_wlmin,in_grp_wlmax,
!flx
     & in_flx_ndim,in_flx_wlmin,in_flx_wlmax,
!src
     & in_src_ns,in_src_n2s,in_srcepwr,in_srcmax,in_surfsrcmu,
     & in_niter_early,in_src_n2s_early,
     & in_L_bol,in_L_inner,
!trn
     & in_trn_errorfatal,in_puretran,
     & in_trn_noamp,
     & in_trn_tauddmc,in_taulump,
     & in_thm_aniso,in_thm_thermal,in_trn_doplegacy,in_innerbc,
     & in_srcdiff,
!gas
     & in_noeos,
     & in_nebular,in_R_phot,
     & in_nlte,in_nlte_nelem,in_nlteh,in_nlteh_nlev,in_nlteh_damp,
     & in_nlteh_tau0,in_nltex_atoms,in_nltex_dir,in_nltex_chgexch,
     & in_srctype,in_shock_dvmin,in_shock_tmax,in_shock_cie,
     & in_nltex_prev,in_nltex_dtprev,in_nltex_init,in_nltex_hold,
     & in_nltex_bbext,in_nltex_cdest,in_nltex_jbbeta,in_nltex_jbwin,
     & in_nltex_dmpsec,in_nltex_javg,in_nltex_taumin,in_nltex_jfill,
     & in_nltex_fbias,in_nltex_fwl1,in_nltex_fwl2,in_nltex_ldbse,
     & in_nltex_auger,in_nltex_ovllmax,
     & in_nltex_dis,
     & in_nltex_ali,
     & in_nltex_ldb,in_nltex_vturb,in_nltex_vturbmax,in_nltex_esrates,
     & in_nltex_esdamp,in_nltex_dopmc,in_nltex_dopz,in_nltex_dopes,
     & in_nltex_overlap,
     & in_fi,in_fi_wl1,in_fi_wl2,in_fi_dv,in_fi_npc,in_fi_navg,
     & in_fi_doppler,in_fi_only,in_fi_nes,in_fi_nesdv,
     & in_nltex_fij,in_nltex_fijwl1,in_nltex_fijwl2,in_nltex_fijdv,
     & in_nltex_fijnes,
     & in_xray_table,in_xray_ff1,in_xray_ff2,in_xray_t1,in_xray_t2,
     & in_xray_v1,in_xray_v2,in_xray_smooth,in_xray_q,
     & in_bbsobolev,in_bbeps,in_rayleigh,in_tesolve,in_te_damp,
     & in_te_zsolve,in_te_secant,
     & in_tion,
     & in_nt_emin,in_nt_sf,in_nt_sfemax,in_nt_sfde,in_nt_lumg,
     & in_lumg_rad,
     & in_fflegacy,
!phys opac
     & in_opacmixrossel,
     & in_nobbopac,in_nobfopac,in_noffopac,in_nothmson,
!anal opac
     & in_opacanaltype,
     & in_ldisp1,in_ldisp2,
     & in_gas_sigcoef,in_gas_sigtpwr,in_gas_sigrpwr,
     & in_gas_capcoef,in_gas_captpwr,in_gas_caprpwr,
!io
     & in_io_grabstdout,
     & in_io_nogriddump,in_io_dogrdtally,
     & in_io_opacdump,in_io_pdensdump,
     & in_io_profdump,in_io_debugging
c
c-- pointers
c
      integer,parameter,private :: npointers = 100
c
      type lptr
       logical,pointer :: p
      endtype lptr
      type(lptr) :: in_l(npointers)
c
      type iptr
       integer,pointer :: p
      endtype iptr
      type(iptr) :: in_i(npointers)
c
      type rptr
       real*8,pointer :: p
      endtype rptr
      type(rptr) :: in_r(npointers)
c
      type cptr
       character(4),pointer :: p
      endtype cptr
      type(cptr) :: in_c(npointers)
c
      public
      private inputpars
      save
c
c
      contains
c
      subroutine inputpar_create_pointers(il,ii,ir,ic)
c     ------------------------------------------------
      implicit none
************************************************************************
* create pointer arrays (not to be confused with a array pointers) to
* all input parameters.  These arrays are used in mpimod to broadcast
* all input parameters at once.
************************************************************************
      integer,intent(out) :: il,ii,ir,ic
c
c-- init
      il=0
      ii=0
      ir=0
      ic=0
c
      call insertl(in_io_grabstdout,in_l,il)
      call inserti(in_nomp,in_i,ii)
      call inserti(in_niter,in_i,ii)
      call inserti(in_grd_igeom,in_i,ii)
      call inserti(in_ndim(1),in_i,ii)
      call inserti(in_ndim(2),in_i,ii)
      call inserti(in_ndim(3),in_i,ii)
      call inserti(in_flx_ndim(1),in_i,ii)
      call inserti(in_flx_ndim(2),in_i,ii)
      call inserti(in_flx_ndim(3),in_i,ii)
      call insertr(in_flx_wlmin,in_r,ir)
      call insertr(in_flx_wlmax,in_r,ir)
      call insertl(in_io_nogriddump,in_l,il)
      call insertl(in_io_dogrdtally,in_l,il)
      call inserti(in_src_ns,in_i,ii)
      call inserti(in_src_n2s,in_i,ii)
      call inserti(in_niter_early,in_i,ii)
      call inserti(in_src_n2s_early,in_i,ii)
      call insertr(in_L_bol,in_r,ir)
      call insertr(in_L_inner,in_r,ir)
      call insertl(in_puretran,in_l,il)
      call insertl(in_trn_errorfatal,in_l,il)
      call insertl(in_trn_noamp,in_l,il)
      call insertr(in_trn_tauddmc,in_r,ir)
      call insertr(in_taulump,in_r,ir)
      call insertl(in_thm_aniso,in_l,il)
      call insertl(in_trn_doplegacy,in_l,il)
      call insertl(in_thm_thermal,in_l,il)
      call insertc(in_innerbc,in_c,ic)
      call insertl(in_srcdiff,in_l,il)
      call inserti(in_grp_ng,in_i,ii)
      call inserti(in_grp_ngs,in_i,ii)
      call insertr(in_grp_wlmin,in_r,ir)
      call insertr(in_grp_wlmax,in_r,ir)
      call insertr(in_opacmixrossel,in_r,ir)
      call insertc(in_opacanaltype,in_c,ic)
      call insertr(in_ldisp1,in_r,ir)
      call insertr(in_ldisp2,in_r,ir)
      call insertr(in_gas_sigcoef,in_r,ir)
      call insertr(in_gas_sigtpwr,in_r,ir)
      call insertr(in_gas_sigrpwr,in_r,ir)
      call insertr(in_gas_capcoef,in_r,ir)
      call insertr(in_gas_captpwr,in_r,ir)
      call insertr(in_gas_caprpwr,in_r,ir)
      call insertc(in_surfsrcmu,in_c,ic)
      call insertr(in_srcmax,in_r,ir)
      call insertr(in_srcepwr,in_r,ir)
      call insertc(in_io_opacdump,in_c,ic)
      call insertc(in_io_pdensdump,in_c,ic)
      call insertl(in_noeos,in_l,il)
      call insertl(in_nebular,in_l,il)
      call insertr(in_R_phot,in_r,ir)
      call insertl(in_nlte,in_l,il)
      call inserti(in_nlte_nelem,in_i,ii)
      call insertl(in_nlteh,in_l,il)
      call inserti(in_nlteh_nlev,in_i,ii)
      call insertr(in_nlteh_damp,in_r,ir)
      call insertr(in_nlteh_tau0,in_r,ir)
      call insertc(in_srctype,in_c,ic)
      call insertr(in_shock_dvmin,in_r,ir)
      call insertr(in_shock_tmax,in_r,ir)
      call insertr(in_nltex_dtprev,in_r,ir)
      call insertr(in_xray_ff1,in_r,ir)
      call insertr(in_xray_ff2,in_r,ir)
      call insertr(in_xray_t1,in_r,ir)
      call insertr(in_xray_t2,in_r,ir)
      call insertr(in_xray_v1,in_r,ir)
      call insertr(in_xray_v2,in_r,ir)
      call insertl(in_xray_smooth,in_l,il)
      call insertr(in_xray_q,in_r,ir)
      call insertl(in_nltex_ali,in_l,il)
      call insertr(in_nltex_ldb,in_r,ir)
      call insertr(in_nltex_vturb,in_r,ir)
      call insertr(in_nltex_vturbmax,in_r,ir)
      call insertl(in_nltex_esrates,in_l,il)
      call insertl(in_nltex_hold,in_l,il)
      call insertl(in_nltex_bbext,in_l,il)
      call insertl(in_nltex_cdest,in_l,il)
      call insertl(in_nltex_jbbeta,in_l,il)
      call insertl(in_nltex_dmpsec,in_l,il)
      call insertr(in_nltex_javg,in_r,ir)
      call insertr(in_nltex_taumin,in_r,ir)
      call inserti(in_nltex_jfill,in_i,ii)
      call insertr(in_nltex_fbias,in_r,ir)
      call insertl(in_nltex_ldbse,in_l,il)
      call insertl(in_nltex_auger,in_l,il)
      call inserti(in_nltex_ovllmax,in_i,ii)
      call insertr(in_nltex_fwl1,in_r,ir)
      call insertr(in_nltex_fwl2,in_r,ir)
      call insertr(in_nltex_jbwin,in_r,ir)
      call insertl(in_nltex_dis,in_l,il)
      call insertl(in_nltex_esdamp,in_l,il)
      call insertr(in_nltex_dopmc,in_r,ir)
      call inserti(in_nltex_dopz,in_i,ii)
      call insertl(in_nltex_dopes,in_l,il)
      call insertr(in_nltex_overlap,in_r,ir)
      call insertl(in_fi,in_l,il)
      call insertr(in_fi_wl1,in_r,ir)
      call insertr(in_fi_wl2,in_r,ir)
      call insertr(in_fi_dv,in_r,ir)
      call inserti(in_fi_npc,in_i,ii)
      call inserti(in_fi_navg,in_i,ii)
      call insertl(in_fi_doppler,in_l,il)
      call insertl(in_fi_only,in_l,il)
      call inserti(in_fi_nes,in_i,ii)
      call insertr(in_fi_nesdv,in_r,ir)
      call inserti(in_nltex_fij,in_i,ii)
      call insertr(in_nltex_fijwl1,in_r,ir)
      call insertr(in_nltex_fijwl2,in_r,ir)
      call insertr(in_nltex_fijdv,in_r,ir)
      call inserti(in_nltex_fijnes,in_i,ii)
      call insertl(in_bbsobolev,in_l,il)
      call insertr(in_bbeps,in_r,ir)
      call insertl(in_tesolve,in_l,il)
      call insertr(in_te_damp,in_r,ir)
      call inserti(in_te_zsolve,in_i,ii)
      call insertl(in_te_secant,in_l,il)
      call insertc(in_tion,in_c,ic)
      call insertr(in_nt_emin,in_r,ir)
      call insertl(in_nt_sf,in_l,il)
      call insertr(in_nt_sfemax,in_r,ir)
      call insertr(in_nt_sfde,in_r,ir)
      call insertr(in_nt_lumg,in_r,ir)
      call insertr(in_lumg_rad,in_r,ir)
      call insertl(in_fflegacy,in_l,il)
      call insertl(in_nobbopac,in_l,il)
      call insertl(in_nobfopac,in_l,il)
      call insertl(in_noffopac,in_l,il)
      call insertl(in_nothmson,in_l,il)
      call insertl(in_io_profdump,in_l,il)
      call insertl(in_io_debugging,in_l,il)
c
      contains
c
      subroutine insertl(par,arr,i)
      implicit none
      logical,intent(in),target :: par
      type(lptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'insertl: i>npointers'
      arr(i)%p => par
      end subroutine insertl
c
      subroutine inserti(par,arr,i)
      implicit none
      integer,intent(in),target :: par
      type(iptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'inserti: i>npointers'
      arr(i)%p => par
      end subroutine inserti
c
      subroutine insertr(par,arr,i)
      implicit none
      real*8,intent(in),target :: par
      type(rptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'insertr: i>npointers'
      arr(i)%p => par
      end subroutine insertr
c
      subroutine insertc(par,arr,i)
      implicit none
      character(4),intent(in),target :: par
      type(cptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'inserti: i>npointers'
      arr(i)%p => par
      end subroutine insertc

      end subroutine inputpar_create_pointers
c
c
      subroutine read_inputpars
c     -------------------------
      implicit none
************************************************************************
* read the input parameter namelist
************************************************************************
      character(15),parameter :: fname='input.par'
c
c-- read namelist
      open(4,file=fname,status='old',err=66)
      read(4,nml=inputpars,end=67,err=68)
      close(4)
c
      return
66    stop 'read_inputpars: namelist input file missing: input.par'
67    stop 'read_inputpars: namelist missing or bad in input.par'
68    stop 'read_inputpars: ivalid parameters or values in namelist'

      end subroutine read_inputpars
c
c
c
      subroutine parse_inputpars(nmpi)
c     --------------------------------
      use miscmod, only:warn
c$    use omp_lib
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* parse the input parameter namelist
************************************************************************
      integer :: istat
c$    integer :: i
c
c-- dump namelist to stdout
      write(6,*) 'namelist read:'
      write(6,nml=inputpars)
      write(6,*)
c
c-- write simulation name to file
      open(4,file='output.name',iostat=istat)
      if(istat/=0) stop 'parse_inputpars: open output.name error'
      write(4,'(a)') trim(in_name)
      close(4)
c
c-- check input parameter validity
      if(in_nomp<0) stop 'in_nomp invalid'
      if(in_nomp==0 .and. nmpi>1) stop 'no in_nomp==0 in mpi mode'
c
      if(any(in_ndim<1)) stop 'in_ndim invalid'
c
      select case(in_grd_igeom)
      case(11)
       if(in_ndim(2)>1 .or. in_ndim(3)>1) stop 'in_ndim invalid'
       if(in_flx_ndim(2)/=1) stop 'in_flx_ndim(2) inval'
       if(in_flx_ndim(3)/=1) stop 'in_flx_ndim(3) inval'
      case default
       stop 'in_grd_igeom invalid'
      endselect
c
      if(in_io_debugging) then
        in_io_nogriddump = .false.
        in_io_profdump = .true.
        in_io_opacdump = 'one'
        in_io_pdensdump = 'on'
        in_io_dogrdtally = .true.
      endif
c
      if(in_io_nogriddump .and. in_io_dogrdtally) stop
     &   'dogridtally and !griddump'
      if(in_io_nogriddump .and. in_io_opacdump.ne.'off') stop
     &   'opacdump and !griddump'
c
      if(in_grp_ng<0) stop 'in_grp_ng invalid'
      if(in_grp_ngs<=0) stop 'in_grp_ngs invalid'
      if(in_grp_wlmin<0) stop 'in_grp_wlmin invalid'
      if(in_grp_wlmax<=in_grp_wlmin) stop 'in_grp_wlmax invalid'
c
      if(in_srcepwr<=0d0) stop 'in_srcpwr <= 0'
      if(in_src_ns<=0 .eqv. in_src_n2s<0) stop
     &  'use in_src_ns or in_src_n2s'
      if(in_niter_early>0 .and. (in_src_n2s_early<0 .or.
     &  in_src_n2s_early>in_src_n2s)) stop
     &  'in_niter_early: 0 <= in_src_n2s_early <= in_src_n2s'
c
      if(in_L_bol.le.0) stop 'in_L_bol<=0'
c
      if(in_taulump<=.05d0*in_trn_tauddmc) stop 'in_taulump too small' !don't let scattering dominate
c
      select case(in_opacanaltype)
      case('none')
        if(in_nobbopac.and.in_nobfopac.and.in_noffopac)
     &   call warn('inputparmod','no phys opac + in_opacanaltype==none')
      case('mono')
      case('grey')
      case('line')
      case default
       stop 'in_opacanaltype unknown'
      end select
c-- disallow physical opacities when analytic opacities are selected
      if(in_opacanaltype/='none' .and. .not.(in_nobbopac .and.
     &  in_nobfopac .and. in_noffopac)) then
       stop 'in_no??opac: no physical opac allowed with analytic opac'
      endif
c
      if(in_nobbopac) call warn('read_inputpars','bb opacity disabled!')
      if(in_nobfopac) call warn('read_inputpars','bf opacity disabled!')
      if(in_noffopac) call warn('read_inputpars','ff opacity disabled!')
      if(in_nothmson) call warn('read_inputpars','Thomson disabled')
c
      select case(in_io_opacdump)
      case('off ')
      case('one ')
      case('all ') !-- for experimental iterative method
       call warn('read_inputpars',
     &   "in_io_opacdump=='all' will generate a big data file!")
      case default
       stop 'in_io_opacdump invalid'
      end select
c
      if(in_opacmixrossel<0d0 .or. in_opacmixrossel>1d0) then
       stop 'in_opacmixrossel invalid'
      endif
c
      if(in_nebular.and.in_R_phot.le.0) stop 'in_R_phot<=0'
      in_nltex = in_nlteh .or. len_trim(in_nltex_atoms)>0
      if(in_nltex.and..not.in_nlte) stop 'NLTE ions require in_nlte'
      if(in_nltex.and.in_grp_ngs/=1) stop 'NLTE ions: in_grp_ngs=1'
      if(in_innerbc/='refl'.and.in_innerbc/='thrm'.and.
     &   in_innerbc/='absb') stop 'in_innerbc: refl|thrm|absb'
      if(in_nlteh_nlev<4 .or. in_nlteh_nlev>60) stop 'in_nlteh_nlev'
      if(in_srctype/='surf' .and. in_srctype/='lumg')
     &  stop 'in_srctype must be surf or lumg'
      if(in_bbeps>1d0) stop 'in_bbeps > 1'
      if(in_rayleigh/='lyman'.and.in_rayleigh/='kurucz'.and.
     &  in_rayleigh/='lymes'.and.in_rayleigh/='none')
     &  stop 'in_rayleigh: lyman|kurucz|lymes|none'
      if(in_rayleigh=='lymes'.and..not.in_thm_thermal)
     &  stop 'in_rayleigh = lymes: needs in_thm_thermal'
      if(in_tesolve.and..not.in_nltex) stop 'in_tesolve: NLTE ions'
      if(in_tion/='hydr'.and.in_tion/='elec'.and.in_tion/='nebu')
     &  stop 'in_tion: hydr|elec|nebu'
      if(in_tion=='nebu'.and..not.in_nltex) stop 'in_tion: NLTE ions'
      if(in_te_damp<=0d0.or.in_te_damp>1d0) stop 'in_te_damp: (0,1]'
      if(trim(in_nltex_prev)/='none') then
       if(in_srctype/='lumg') stop 'in_nltex_prev: needs lumg'
       if(len_trim(in_nltex_atoms)==0) stop 'in_nltex_prev: no atoms'
       if(in_nltex_dtprev<=0d0) stop 'in_nltex_dtprev: must be > 0'
      endif
      if(in_nt_emin>0d0.and..not.in_nltex) stop 'in_nt_emin: NLTE'
      if(in_nt_sf.and..not.in_nltex) stop 'in_nt_sf: NLTE'
      if(in_nt_lumg>0d0.and..not.in_nt_sf) stop 'in_nt_lumg: in_nt_sf'
      if(in_nt_lumg>0d0.and.in_srctype/='lumg')
     &  stop 'in_nt_lumg: in_srctype=lumg'
      if(in_nt_lumg<0d0.or.in_nt_lumg>1d0) stop 'in_nt_lumg: [0,1]'
      if(in_nltex_javg<0d0.or.in_nltex_javg>=1d0)
     &  stop 'in_nltex_javg: [0,1)'
      if(in_nobfopac) in_nltex_auger = .false.
      if(in_xray_q<0d0.or.in_xray_q>=1d0) stop 'in_xray_q: [0,1)'
      if(in_nltex_fij<0) then
       in_nltex_fij = 0
       if(trim(in_xray_table)/='none') in_nltex_fij = 1
      endif
      if(in_nltex_fij>0.and.(in_nltex_fijnes<1.or.
     &  in_nltex_fijwl1>=in_nltex_fijwl2.or.in_nltex_fijdv<=0d0))
     &  stop 'in_nltex_fij: fijnes>=1, fijwl1<fijwl2, fijdv>0'
      if(in_nltex_ldb<0d0) then
       in_nltex_ldb = 0d0
       if(.not.in_puretran.and.in_trn_tauddmc<=30d0) in_nltex_ldb = 30d0
      endif
      if(in_nltex_fbias>1d0.and.in_nltex_fwl1>=in_nltex_fwl2)
     &  stop 'in_nltex_fwl1 < in_nltex_fwl2'
      if(in_lumg_rad<0d0.or.in_nt_lumg+in_lumg_rad>1d0)
     &  stop 'in_lumg_rad: [0,1-in_nt_lumg]'
      if(in_lumg_rad>0d0.and.in_srctype/='lumg')
     &  stop 'in_lumg_rad: in_srctype=lumg'
      if(in_nltex_hold.and.trim(in_nltex_init)=='none')
     &  stop 'in_nltex_hold: needs in_nltex_init'
      if(in_nltex_jbwin<0d0) stop 'in_nltex_jbwin < 0'
      if(in_nltex_overlap<0d0) stop 'in_nltex_overlap < 0'
c
      if(in_nlte.and.in_nlte_nelem.ne.1) stop 'in_nlte_nelem!=1'
c
      select case(in_io_pdensdump)
      case('off ')
      case('on ')
      case default
       stop 'in_io_pdensdump invalid'
      end select
c
c-- set the number of threads
c$    if(.false.) then
c-- serial run
       in_nomp = 1
c$    else
c-- openmp run
c$     if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c$omp parallel shared(in_nomp) private(i)
c$     i = omp_get_num_threads()
c$     if(in_nomp/=0 .and. i/=in_nomp)
c$   &   stop 'read_inputpars: in_nomp error'
c$     in_nomp = i
c$omp end parallel
c$    endif
      write(6,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &  nmpi,in_nomp,nmpi*in_nomp
      if(in_io_grabstdout) then
       write(0,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &   nmpi,in_nomp,nmpi*in_nomp
      endif
c
      end subroutine parse_inputpars
c
c
c
      subroutine provide_inputpars(nmpi)
c     -------------------------------------
      use physconstmod
      use particlemod
      use sourcemod
      use transportmod
      use fluxmod
      use groupmod
      use gasmod
      use gridmod
      use nltemod
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* Distribute the input parameter values to the respective modules.
* This needs to be called AFTER the values are bcast.
************************************************************************
      integer :: mpart
      integer :: ns
c
      mpart = int(int(2,8)**in_src_n2s/nmpi) !max number of particles
      mpart = max(mpart,in_src_ns/nmpi)
      prt_npartmax = mpart
c
      ns = int(int(2,8)**in_src_n2s/nmpi)
      ns = max(ns,in_src_ns/nmpi)
      src_ns = ns
c
      trn_tauddmc = in_trn_tauddmc
      trn_taulump = in_taulump
      trn_thm_aniso = in_thm_aniso
      trn_dopexact = .not.in_trn_doplegacy
      trn_thm_thermal = in_thm_thermal
      trn_raythm = in_rayleigh=='lymes'
      trn_innerbc = in_innerbc
      trn_tinner = in_srcmax
      trn_errorfatal = in_trn_errorfatal
      trn_noampfact = in_trn_noamp
c
      flx_ndim  = in_flx_ndim
      flx_wlmin = in_flx_wlmin
      flx_wlmax = in_flx_wlmax
c
      grp_ng    = in_grp_ng
      grp_ngs   = in_grp_ngs
      grp_wlmin = in_grp_wlmin
      grp_wlmax = in_grp_wlmax
c
      grd_igeom = in_grd_igeom
      grd_nx    = in_ndim(1)
      grd_ny    = in_ndim(2)
      grd_nz    = in_ndim(3)
c
      nlte_nelem = in_nlte_nelem
      R_phot = in_R_phot
c
      end subroutine provide_inputpars
c
      end module inputparmod
c vim: fdm=marker
