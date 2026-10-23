!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
program superlite

  use randommod
  use sourcemod
  use mpimod
  use transportmod
  use inputparmod
  use groupmod
  use gridmod
  use gasmod
  use particlemod
  use physconstmod
  use miscmod
  use totalsmod

  use inputstrmod
  use fluxmod

  use ionsmod, only:ions_read_data,ions_alloc_grndlev
  use nltemod
  use bfxsmod, only:bfxs_read_data
  use ffxsmod, only:ffxs_read_data
  use timingmod
  use countersmod
  use nltehmod, only:nlteh_set_collfits,nlteh_init
  use nlteatommod, only:nxa_nel,nxa_nion,nxa_ion,nxa_nlin,nxa_nst,nxa_read_cx
  use nltextrmod
  use xraysrcmod
  use formalmod, only:fi_spectrum
  use nltexmod, only:nx_esdamp,nx_ovl_dv,nxo_npair,nx_pc_on
  use nltedismod, only:nxd_on,nxd_setup
  use bbsobmod
  use tebalmod, only:te_damp,te_ntemin,te_zsolve,te_secant
  use nltentmod, only:nxn_on,nxn_emax,nxn_de,nxn_init
  use shockciemod, only:shc_read,shc_setup,shc_nel,shc_elz

  implicit none
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer :: ierr, it
  integer :: icell1, ncell !number of cells per rank (gas_ncell)
  real*8 :: t0, t1 !timing
  character(25) :: msg
!-- NLTE model-atom files
  integer :: i, nfatom, j, l
  real*8,allocatable :: shcx(:,:)
  real*8,allocatable :: mcell(:),xh(:),xhe(:)
  real*8,allocatable :: jgfi(:,:) !formal integral: J of the last in_fi_navg iterations
  integer :: nfi
  logical :: lclump = .false.
  character(256) :: fatom(150)
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
!
!-- initialize timing module
  call timingmod_init
  call countersmod_init
!
!--
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
  if(lmpi0) then
     t0 = t_time()
!-- startup message
     call banner
!-- read runtime parameters
     call read_inputpars
!-- redirect stdout to file if selected
     call open_logfiles
!-- parse and verify runtime parameters
     call parse_inputpars(nmpi)
!
!-- read input structure
     call read_inputstr(in_grd_igeom,in_ndim,nmpi)
!-- compressed domain, serialize non-void cells
     call inputstr_compress

!-- READ DATA
!-- read ion and level data
     call ions_read_data(gas_nelem)  !ion and level data
!-- read bbxs data
     if(.not.in_nobbopac) then!bound-bound cross section data
       if(in_nlte) then !NLTE
         call read_nlte_data(gas_nelem,in_nlte_nelem)
       else !lte & nebular approximation
         call read_bbxs_data(gas_nelem)
       endif
     endif
!-- read bfxs data
     if(.not.in_nobfopac) call bfxs_read_data !bound-free cross section data
!-- read ffxs data
     if(.not.in_noffopac) call ffxs_read_data !free-free cross section data
!
     t1 = t_time()
     t_setup = t1-t0
  endif !impi
!-- broadcast init info from impi0 rank to all others
  call bcast_permanent !MPI
  call provide_inputpars(nmpi)
!-- domain-decompose input structure
  call scatter_inputstruct(in_ndim,icell1,ncell) !MPI

!--
!-- setup remaining modules
!==========================
!-- wlgrid (before grid setup)
  call groupmod_init
  call fluxgrid_setup

!-- setup spatial grid
  call gridmod_init(lmpi0,grp_ng,str_nc,icell1,ncell)
  call grid_setup
!-- setup gas
  call gasmod_init(lmpi0,icell1,ncell,grp_ng)
  call gas_setup
!-- microclumping: volume filling factors of all cells
  allocate(mcell(grd_ncell))
  mcell = 0d0
  mcell(icell1:icell1+ncell-1) = gas_fcl
  call allreduce_sum(grd_ncell,mcell)
  where(mcell>0d0) grd_fcl = mcell
  deallocate(mcell)
  lclump = any(grd_fcl<1d0)
  if(lclump .and. lmpi0) write(6,'(1x,a,f7.4,a,f7.4)') 'microclumping: volume filling factor', &
       minval(grd_fcl),' -',maxval(grd_fcl)
!-- inputstr no longer needed
  call inputstr_dealloc

!-- create procedure pointers for the selected geometry
  call transportmod_init(grd_igeom)
!-- source particles
  call sourcemod_init(nmpi)
  call shock_temperatures
  if(lmpi0 .and. in_srctype=='lumg') write(6,'(1x,a,i4,a,es9.2,a)') 'shock cells:', &
       count(src_tsh>0d0),', max post-shock T',maxval(src_tsh),' K'
!-- CIE spectra of the shocks: table and the composition of the shock cells
  if(in_srctype=='lumg' .and. trim(in_shock_cie)/='none' .and. any(src_tsh>0d0)) then
     call shc_read(in_shock_cie)
     allocate(shcx(shc_nel,grd_ncell))
     shcx = 0d0
     do i=1,gas_ncell
        l = icell1 + i - 1
        if(src_tsh(l) <= 0d0) cycle
        do j=1,shc_nel
           if(shc_elz(j) <= gas_nelem) shcx(j,l) = gas_natomfr(shc_elz(j),i)
        enddo
     enddo
     call allreduce_sum(size(shcx),shcx)
     call shc_setup(grd_ncell,src_tsh,shcx)
     deallocate(shcx)
     if(lmpi0) write(6,'(1x,2a,i3,a)') 'shocks: CIE spectra from ',trim(in_shock_cie),shc_nel,' elements'
  endif

!-- allocate arrays of sizes retreived in bcast_permanent
  call ions_alloc_grndlev(gas_nelem,gas_ncell)  !ground state occupation numbers
  if(in_nlte) call nlte_alloc_nlev(nlte_nelem,gas_ncell) !NLTE
!-- NLTE elements with line resonances and macro-atom: built-in hydrogen
!-- (in_nlteh) and the ions of the model-atom files
  if(in_nltex) then
     if(in_nlteh) then
        call nlteh_init(in_nlteh_nlev)
!-- H collision fits of data.nlte.h1 (tabulated by the built-in atom)
        if(allocated(nlte_data)) call nlteh_set_collfits(nlte_data(1)%i(1)%nlev,nlte_data(1)%i(1)%nplev, &
             nlte_data(1)%i(1)%glev,nlte_data(1)%i(1)%coll_nlin, &
             nlte_data(1)%i(1)%llw_coll,nlte_data(1)%i(1)%lup_coll, &
             nlte_data(1)%i(1)%C0,nlte_data(1)%i(1)%C1,nlte_data(1)%i(1)%C2, &
             nlte_data(1)%i(1)%C3,nlte_data(1)%i(1)%delE)
     endif
     call atomfiles(in_nltex_atoms,nfatom,fatom)
     if(len_trim(in_nltex_dir)>0) then
        do i=1,nfatom
           if(fatom(i)(1:1)/='/') fatom(i) = trim(in_nltex_dir)//'/'//trim(fatom(i))
        enddo
     endif
!-- (overlapping lines: the partners within in_nltex_overlap, nxt_init)
     nx_ovl_dv = in_nltex_overlap*1d5
     call nxt_init(grd_ncell,gas_ncell,grp_ng,grp_wl,nfatom,fatom(:nfatom),merge(in_nlteh_nlev,0,in_nlteh))
     if(lmpi0 .and. nxo_npair>0) write(6,'(1x,a,f7.1,a,i8,a)') 'NLTE lines: overlapping resonance zones '// &
          'within',in_nltex_overlap,' km/s (',nxo_npair/2,' pairs)'
!-- charge exchange among the NLTE elements
     if(trim(in_nltex_chgexch)/='none' .and. len_trim(in_nltex_chgexch)>0) then
        fatom(1) = in_nltex_chgexch
        if(len_trim(in_nltex_dir)>0 .and. fatom(1)(1:1)/='/') fatom(1) = trim(in_nltex_dir)//'/'//trim(fatom(1))
        call nxa_read_cx(fatom(1),lmpi0)
     endif
!-- time-dependent statistical equilibrium: populations of the previous model
     if(trim(in_nltex_prev)/='none') then
        allocate(mcell(grd_ncell))
        mcell = 0d0
        mcell(icell1:icell1+ncell-1) = gas_mass
        call allreduce_sum(grd_ncell,mcell)
        call nxt_read_prev(in_nltex_prev,grd_ncell,mcell,in_nltex_dtprev,lmpi0)
        deallocate(mcell)
     endif
!-- warm start: the populations (and T_e with in_tesolve) of a saved state
     if(trim(in_nltex_init)/='none') then
        allocate(mcell(grd_ncell))
        mcell = 0d0
        mcell(icell1:icell1+ncell-1) = gas_mass
        call allreduce_sum(grd_ncell,mcell)
        call nxt_read_init(in_nltex_init,grd_ncell,mcell,lmpi0)
        if(in_fi_only .and. .not.nxt_have_mc) stop 'in_fi_only: needs the output.nltex_jmc of in_nltex_init'
        nxt_hold = in_nltex_hold
!-- (in_niter_early is meant for runs from scratch: in a warm start it
!-- reduces the packets of the first iterations all the same)
        if(lmpi0 .and. in_niter_early>0) write(6,'(1x,a,i0,a,i0,a)') &
             'warning: warm start with in_niter_early = ',in_niter_early, &
             ': these first iterations have 2^',in_src_n2s_early,' source particles (in_src_n2s_early)'
        if(lmpi0 .and. nxt_hold) write(6,'(1x,a)') &
             'NLTE populations held at the warm start in every iteration (in_nltex_hold)'
        deallocate(mcell)
        do i=1,ncell
           nxt_xloc(:,i) = nxt_xinit(:,icell1+i-1)
           if(in_tesolve) gas_temp(i) = nxt_tinit(icell1+i-1)
        enddo
     endif
     nxt_damp = in_nlteh_damp
     nxt_ldbtau = in_nltex_ldb
     nxt_bbext = in_nltex_bbext
     if(lmpi0 .and. .not.nxt_bbext) write(6,'(1x,a)') &
          'NLTE elements: only the lines of the model atoms (in_nltex_bbext = f)'
     nx_pc_on = in_nltex_cdest
     if(lmpi0) write(6,'(1x,a,l2)') 'NLTE lines: continuum destruction of the trapped photons (in_nltex_cdest):',nx_pc_on
     nxt_jbbeta = in_nltex_jbbeta
     nxt_dmpsec = in_nltex_dmpsec
     nxt_javg = in_nltex_javg
     nxt_taumin = in_nltex_taumin
     nxt_jfill = in_nltex_jfill
     if(lmpi0 .and. nxt_javg>0d0) write(6,'(1x,a,f5.2,a)') &
          'NLTE group J: weight',nxt_javg,' of the previous iterations (in_nltex_javg)'
     if(lmpi0 .and. nxt_dmpsec) write(6,'(1x,a)') 'NLTE populations: adaptive damping between iterations (in_nltex_dmpsec)'
     if(lmpi0) write(6,'(1x,a,l2)') 'NLTE lines: incident intensity weighted by beta(mu)/beta_bar (in_nltex_jbbeta):', &
          nxt_jbbeta
!-- level dissolution (occupation probabilities, CMFGEN's DO_DIS)
     nxd_on = in_nltex_dis
     if(nxd_on) then
        call nxd_setup(grp_ng,grp_wl)
        if(lmpi0) write(6,'(1x,a)') 'NLTE levels: dissolution by the plasma microfield (in_nltex_dis)'
     endif
     nxt_jbwin = in_nltex_jbwin*1d5
     if(lmpi0 .and. nxt_jbwin>0d0) write(6,'(1x,a,f7.1,a)') &
          'NLTE lines: incident intensity from windows of',in_nltex_jbwin,' km/s on their incident side'
!-- electron scattering of line photons trapped in their resonance zones
!-- (the electrons' Doppler shifts: thermal Thomson scattering)
     nxt_estrap = in_nltex_vturb>=0d0 .and. in_thm_thermal
     nxt_esrates = nxt_estrap .and. in_nltex_esrates
     call nxt_alloc_opt
     nxt_esdamp = in_nltex_esdamp
     nx_esdamp = in_nltex_esdamp
!-- (the microturbulence also widens the line profiles of the formal integral)
     if(in_nltex_vturb>=0d0) then
        nxt_vturb = in_nltex_vturb*1d5
        if(in_nltex_vturbmax>=0d0) nxt_vturb(2) = in_nltex_vturbmax*1d5
        nxt_vmax = grd_vxarr(grd_nx+1)
     endif
     if(nxt_estrap) then
        if(lmpi0) write(6,'(1x,a,f7.2,a,f7.2,a,l2)') 'NLTE lines: electron scattering in the resonance zones, '// &
             'microturbulence',nxt_vturb(1)*1d-5,' to',nxt_vturb(2)*1d-5,' km/s; in the rates:',nxt_esrates
     elseif(in_nltex_vturb>=0d0 .and. lmpi0) then
        write(6,'(1x,a,f7.2,a,f7.2,a)') 'NLTE lines: microturbulence',nxt_vturb(1)*1d-5,' to', &
             nxt_vturb(2)*1d-5,' km/s in the formal integral (the trapped-photon scattering needs in_thm_thermal)'
     endif
     nxt_ldbone = .not.in_puretran .and. in_trn_tauddmc<=in_nltex_ldb
     nxt_alion = in_nltex_ali
!-- Doppler-profile lines in the transport of slow cells
     nxt_dopmc = max(in_nltex_dopmc,0d0)
     nxt_dopz = max(in_nltex_dopz,1)
     nxt_dopes = nxt_dopmc>0d0 .and. nxt_estrap .and. in_nltex_dopes
     if(nxt_dopmc>0d0 .and. lmpi0) write(6,'(1x,a,i3,a,f6.2,a,l2)') 'NLTE lines: static Doppler profiles of the '// &
          'lines from Z =',nxt_dopz,' in the transport of the cells whose velocity changes by less than', &
          nxt_dopmc,' Doppler widths; their trapped-photon electron scattering in the rates:',nxt_dopes
     if(lclump) then
        allocate(nxt_fcl(grd_ncell))
        nxt_fcl = grd_fcl
     endif
     nxt_tau0 = in_nlteh_tau0
     nxt_ntemin = in_nt_emin*pc_ev
!-- fast electrons: Spencer-Fano degradation with all NLTE species
     if(in_nt_sf) then
        nxn_on = .true.
        nxn_emax = in_nt_sfemax
        nxn_de = in_nt_sfde
        call nxn_init(in_nltex_dir,grd_ncell,ncell,lmpi0)
        if(lmpi0 .and. in_nt_lumg>0d0) write(6,'(1x,a,f6.3,a)') 'fast electrons deposit',in_nt_lumg, &
             ' of the lumg energy of each cell (in_nt_lumg)'
     endif
     if(lmpi0) then
        write(6,'(1x,a,i3,a,i3,a,i6,a,i6,a)') 'NLTE:',nxa_nel,' elements,',nxa_nion,' ions,', &
             nxa_nst,' states,',nxa_nlin,' lines (Sobolev resonances, macro-atom)'
        do i=1,nxa_nion
           if(nxa_ion(i)%nsl < nxa_ion(i)%nlev) then
              write(6,'(3x,a8,i5,a,i5,a,i6,a)') nxa_ion(i)%name,nxa_ion(i)%nlev,' levels in', &
                   nxa_ion(i)%nsl,' super-levels',nxa_ion(i)%nlin,' lines'
           else
              write(6,'(3x,a8,i5,a,i6,a)') nxa_ion(i)%name,nxa_ion(i)%nlev,' levels',nxa_ion(i)%nlin,' lines'
           endif
        enddo
     endif
  endif
!-- fast photoelectrons from X-ray absorption
  te_ntemin = in_nt_emin*pc_ev
!-- electron temperature from the thermal balance
  if(in_tesolve) then
     te_damp = in_te_damp
     te_zsolve = in_te_zsolve
     te_secant = in_te_secant
     if(lmpi0) write(6,'(1x,a,f5.2,2a)') 'T_e from the thermal balance, damping', &
          in_te_damp,', LTE elements at T_ion = ',in_tion
  endif
!-- lines of the LTE elements: Sobolev expansion opacity and line scattering
  bbs_sobolev = in_bbsobolev
  bbs_eps = in_bbeps
  call bbs_init(grp_ng,grd_ncell,gas_ncell)
  if(lmpi0 .and. bbs_on) then
     if(bbs_eps<0d0) then
        write(6,'(1x,a,l2,a)') 'bb lines: Sobolev opacity',bbs_sobolev, &
             ', thermalization from van Regemorter rates'
     else
        write(6,'(1x,a,l2,a,f6.3)') 'bb lines: Sobolev opacity',bbs_sobolev, &
             ', thermalization share',bbs_eps
     endif
  endif
  call particle_alloc(lmpi0)

!-- initialize random number generator, use different seeds for each rank
  call rnd_init(in_nomp,impi)
!-- use extra stream for non-threaded code
  rnd_state = rnd_states(1)

!-- volume
  if(str_lvol) then ! gather grd_vol from gas_vol if vol in input.str
    call allgather_gridvolume
  else ! calculate grd_vol using grd_xarr
    call grid_volume(grd_igeom)
  endif

!-- gather initial radiation temperature
  call allgather_radtemp
!-- X-ray emission of wind shocks (CMFGEN's INC_XRAYS)
  if(trim(in_xray_table)/='none') then
     call xs_read(in_xray_table,[in_xray_t1,in_xray_t2])
     allocate(mcell(grd_ncell),xh(grd_ncell),xhe(grd_ncell))
     mcell = 0d0; xh = 0d0; xhe = 0d0
     mcell(icell1:icell1+ncell-1) = gas_natom
     xh(icell1:icell1+ncell-1) = gas_natomfr(1,:)
     xhe(icell1:icell1+ncell-1) = gas_natomfr(2,:)
     call allreduce_sum(grd_ncell,mcell)
     call allreduce_sum(grd_ncell,xh)
     call allreduce_sum(grd_ncell,xhe)
     call xs_setup(grd_ncell,0.5d0*(grd_vxarr(:grd_ncell)+grd_vxarr(2:)),grd_vol,grd_fcl,mcell,xh,xhe, &
          [in_xray_ff1,in_xray_ff2],1d5*[in_xray_v1,in_xray_v2],in_xray_smooth)
     deallocate(mcell,xh,xhe)
     if(lmpi0) write(6,'(1x,a,es10.3,a,2(f6.3,es9.2,a))') 'X-ray emission of wind shocks:',sum(xs_pow), &
          ' erg/s (f, T:',in_xray_ff1,in_xray_t1,' K;',in_xray_ff2,in_xray_t2,' K)'
  endif


!-- memory usage
  if(lmpi0) then
     msg = 'post setup:'
     write(6,*) 'memusg: ',msg,memusg()
  endif

!-- Monte Carlo loop
!=================
  if(lmpi0) then
     write(6,*)
     write(6,*) "MC calculations:"
     write(6,*) "===================="
  endif
!
  do it=1,in_niter ! iterations
     t_timelin(1) = t_time() !timeline

!-- source particles of this iteration
     if(it<=in_niter_early) then
        src_ns = int(int(2,8)**in_src_n2s_early/nmpi)
     else
        src_ns = int(max(int(2,8)**in_src_n2s/nmpi,int(in_src_ns/nmpi,8)))
     endif
!-- write timestep
     if(lmpi0) then !NLTE
      write(6,*)
      write(6,'(1x,a4,1x,i2,a,i11)') 'it: ',it,'  source particles:',int(src_ns,8)*nmpi
     endif
!-- update all non-permanent variables
     call gas_update(it)

     call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
     t_timelin(2) = t_time() !timeline

!-- gather from gas workers and broadcast to world ranks
     call bcast_nonpermanent !MPI
!-- NLTE elements: Sobolev strengths, macro-atom and k-packet tables
     if(in_nltex) call nxt_tables(grd_ncell,grd_tempinv,grd_radtemp, &
          grd_capemitgrey,grd_xarr,grd_vxarr)
     if(in_nltex .and. nxt_dopmc>0d0 .and. lmpi0) then
        if(allocated(nxt_slow)) write(6,'(1x,a,i5,a,i5,a)') 'NLTE lines: static Doppler profiles in', &
             count(nxt_slow),' of',grd_ncell,' cells'
     endif
     if(in_nltex .and. lmpi0 .and. nxt_bfstat(1)>0d0) write(6,'(1x,a,f7.4,a,f7.4)') &
          'NLTE bound-free sampling: mean share of full walks',nxt_bfstat(2)/nxt_bfstat(1), &
          ', groups above 0.1:',nxt_bfstat(3)/nxt_bfstat(1)
!-- microclumping: the transport sees the volume-mean opacities
     if(lclump) call clump_scale
!-- source energy
     call sourceenergy
!-- the X-ray source and the shock cells in the source functions of the
!-- thick groups (local detailed balance, DDMC spectra) and in the formal
!-- integral
     if(in_nltex .and. xs_on) call xray_emnbf
     if(in_nltex .and. in_srctype=='lumg') then
        if(any(src_tsh>0d0)) call shock_emnbf
     endif
!-- only the formal integral: the saved state and its Monte Carlo J
     if(in_fi_only) then
        call fi_spectrum(grd_ncell,grd_xarr,grd_vxarr, &
             grd_tempinv,in_fi_wl1,in_fi_wl2,in_fi_dv,in_fi_npc,impi,nmpi,lmpi0,dop=in_fi_doppler,nes=in_fi_nes,nesdv=in_fi_nesdv)
        exit
     endif

!-- update photospheric radius for NLTE calculations (experimental iterative method)
     if(in_nlte) call get_rphot
     if(lmpi0.and.in_nlte) then
       write(6,'(1x,a8,1x,es8.2)') "R_phot: ", R_phot
       write(6,*)
     endif
!
     if(lmpi0) write(6,*) "--> Calculating leakage opacity and Emission probability"
!-- source energy
     call ddmc_spectra          !Planck and DDMC source spectra of all cells
     call leakage_opacity       !IMC-DDMC albedo coefficients and DDMC leakage opacities
     call emission_probability  !emission probabilities for ep-group in each cell
     call allgather_leakage !MPI
!-- DDMC group lumping of all cells
     if(.not.in_puretran) call ddmc_lumps

     if(lmpi0) write(6,*) "--> Assign source particles"
     call sourcenumbers         !number of source prt_particles

     if(src_nnew>0) then
        allocate(src_ivacant(src_nnew))
        call vacancies          !Storing vacant "prt_particles" indexes in ordered array "src_ivacant"
        call boundary_source    !properties of prt_particles on domain boundary
        call interior_source!(it) !properties of prt_particles emitted in domain interior
        call source_transformdirection
        deallocate(src_ivacant)
     endif
!-- tally source particles for output
     call sourcetally
!-- write output for source particle distribution
     call reduce_srctally
     if(lmpi0) call output_source
!
!
!-- advance particles
     if(lmpi0) write(6,*) "--> Particle MC tranport"
     t_timelin(3) = t_time() !timeline
     if(in_nltex) call nxt_tally_reset
     call particle_advance
     call reduce_gridtally !MPI  !collect transport results from all workers
     call allreduce_jrad !MPI !NLTE !collect radiation intensity from all workers
!-- NLTE elements: Monte Carlo estimators for the next iteration
     if(in_nltex) then
        call allreduce_nlteh !MPI
        call nxt_store_jg(grd_ncell,grp_ng,grd_jrad,grd_vol)
        call nxt_ldb(grd_ncell,grp_ng,grd_vol,grd_xarr(2:grd_ncell+1)-grd_xarr(1:grd_ncell), &
             grd_tempinv,grd_cap,grd_capemit,grd_planck,grd_spk,grd_nlump,grd_glump)
        call nxt_jfill_cells(grd_ncell,grp_ng)
        call nxt_tally_finish(grd_ncell,grd_vol)
!-- formal integral of the emergent spectrum with the radiation field of the
!-- last in_fi_navg iterations
        if(in_fi .and. it>in_niter-max(in_fi_navg,1)) then
           if(.not.allocated(jgfi)) then
              allocate(jgfi(grp_ng,grd_ncell))
              jgfi = 0d0
              nfi = 0
           endif
           jgfi = jgfi + nxt_jg
           nfi = nfi + 1
        endif
        if(in_fi .and. it==in_niter) then
           jgfi = jgfi/nfi
           call fi_spectrum(grd_ncell,grd_xarr,grd_vxarr, &
                grd_tempinv,in_fi_wl1,in_fi_wl2,in_fi_dv,in_fi_npc,impi,nmpi,lmpi0,jgfi,in_fi_doppler,in_fi_nes,in_fi_nesdv)
           if(lmpi0 .and. nfi>1) write(6,'(1x,a,i3,a)') 'formal integral with the J of the last',nfi,' iterations'
           deallocate(jgfi)
        endif
     endif
     call scatter_jrad !MPI

!-- estimated radiation energy
     if(in_nlte) call nlte_rad_est

!-- print packet advance load-balancing info
     if(lmpi0) then
        call timereg(t_pcktmin,t_pckt_stat(1))
        call timereg(t_pcktmea,t_pckt_stat(2))
        call timereg(t_pcktmax,t_pckt_stat(3))
     endif
     t_timelin(4) = t_time() !timeline

!-- tally flux packets
     call fluxtally
!-- output
     if(lmpi0) write(6,*) "--> Tally flux packets"
     call reduce_fluxtally !MPI  !collect flux results from all workers
     if(lmpi0) call output_flux
!
     if(lmpi0) write(6,*) "--> Print output files"

     call reduce_gasprop
     call reduce_totals

!-- output
     if(lmpi0) then

!-- write output
        call output_grid
        call output_grp
        if(it>1.and.in_io_opacdump=='one') in_io_opacdump = 'off'
        if(it==in_niter.and.in_io_profdump) call output_profile
!-- (NLTE populations and radiation field: every iteration, to follow the
!-- convergence)
        if(in_nltex) call nxt_write(grd_ncell,grd_rcell, &
             0.5d0*(grd_vxarr(:grd_ncell)+grd_vxarr(2:)),grd_tempinv,grd_radtemp,grd_vol)
!-- (the state and its radiation field every iteration as well: a run that
!-- stops can be continued from the last one with in_nltex_init)
        if(in_nltex) call nxt_write_state(grd_ncell,grd_mass,grd_vol,grd_tempinv)

!-- write stdout
        write(6,*)
        write(6,*) "MC Results:"
        write(6,*) "===================="
        write(6,'(2(1x,a7))') 'nsrc','nflux'
        if(ct_nnonvacant(2)<1000000) then
           write(6,'(2(1x,i7))') ct_npcreate(2),ct_npflux(2)
        else
           write(6,'(2(i7,"k"))') ct_npcreate(2)/1000,ct_npflux(2)/1000
        endif
        if(flx_lumout>0d0) write(6,'(1x,a,es10.3,a,f7.4,a)') 'escaped outside the flux bins:', &
             flx_lumout,' erg/s (',flx_lumout/max(flx_lumout+sum(flx_luminos),1d-300),' of the total)'
        if(it==in_niter) call write_emitchan
        write(6,*)
     endif !impi
!-- write timestep timing to file
     call timing_cycle(impi)
     call counters_cycle(impi)
     t_timelin(5) = t_time() !timeline
     t_timeline = t_timeline + (t_timelin(2:) - t_timelin(:4))

  enddo ! iterations
!
!
!--
!-- FINISH UP:
!=============
  call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output
  if(lmpi0) then
!
!-- print memory usage
     msg = 'post loop:'
     write(6,*)
     write(6,*) 'memusg: ',msg,memusg()

!-- print cpu timing usage
     t1 = t_time()
     t_all = t1 - t0
     call print_counters
     call print_timing
     write(6,*)
     write(6,*) 'SuperLite finished'
     if(in_io_grabstdout) write(0,'(a,f8.2,"s")')'SuperLite finished',t_all!repeat to stderr
  endif
!-- Clean up memory. (This helps to locate memory leaks)
  call dealloc_all
  call mpi_finalize(ierr) !MPI

contains
!
subroutine clump_scale
!-- microclumping (volume filling factor f): the opacities and emissivities
!-- of the gas are per unit clump volume (clump densities, Sobolev depths
!-- of the mean opacity); the transport takes their volume means, f times
!-- these (the NLTE Sobolev strengths are means already, nxt_tables)
  use bbsobmod, only:bbs_on,grd_capsc
  integer :: l
  real*4 :: f4
  do l=1,grd_ncell
     if(grd_fcl(l)>=1d0) cycle
     f4 = real(grd_fcl(l))
     grd_cap(:,l) = f4*grd_cap(:,l)
     grd_capemit(:,l) = f4*grd_capemit(:,l)
     grd_emiss(:,l) = f4*grd_emiss(:,l)
     grd_sig(l) = grd_fcl(l)*grd_sig(l)
     grd_capgrey(l) = grd_fcl(l)*grd_capgrey(l)
     grd_capemitgrey(l) = grd_fcl(l)*grd_capemitgrey(l)
     grd_capross(l) = grd_fcl(l)*grd_capross(l)
     if(in_nltex) grd_capnbf(:,l) = f4*grd_capnbf(:,l)
     if(in_nltex) grd_emnbf(:,l) = f4*grd_emnbf(:,l)
     if(bbs_on) grd_capsc(:,l) = f4*grd_capsc(:,l)
  enddo
end subroutine clump_scale
!
subroutine xray_emnbf
!-- add the emission of the X-ray source to grd_emnbf (the emissivity of the
!-- source functions, in the units of the Planck group integrals: Int eta
!-- dnu/(sigma T^4/pi)).  In a group so thick that the radiation is in
!-- local detailed balance (in_nltex_ldb, DDMC) the photoionizations then
!-- take the X-ray photons emitted there, as CMFGEN's source function
!-- holds its X-ray emissivity: without it J = S of the gas alone, and an
!-- X-ray-irradiated wind kept the He ionization of the wind without X-rays.
  integer :: l
  real*8,allocatable :: em(:,:)
  allocate(em(grp_ng,grd_ncell))
  call xs_group_emis(grd_ncell,grp_ng,grp_wl,grd_vol,em)
  do l=1,grd_ncell
     grd_emnbf(:,l) = grd_emnbf(:,l) + real(em(:,l)*pc_pi*grd_tempinv(l)**4/pc_sb)
  enddo
  deallocate(em)
end subroutine xray_emnbf
!
subroutine shock_emnbf
!-- add the emission of the shock cells to grd_emnbf as xray_emnbf: their
!-- injected energy (without the share of the wind X-rays), emitted as the
!-- cooling flow of interior_source, uniform in T between the gas and the
!-- post-shock temperatures (in_shock_tmax at most), CIE spectra
!-- (in_shock_cie) or the bremsstrahlung shape
  use shockciemod, only:shc_group_frac
  integer :: l
  real*8 :: tlo,thi,p
  real*8,allocatable :: f(:)
  allocate(f(grp_ng))
  do l=1,grd_ncell
     if(src_tsh(l)<=0d0 .or. grd_emit(l)<=0d0 .or. grd_vol(l)<=0d0) cycle
     p = grd_emit(l)
     if(xs_on) p = p*(1d0 - xs_px(l))
     tlo = 1d0/grd_tempinv(l)
     thi = max(tlo,min(src_tsh(l),in_shock_tmax))
     call shc_group_frac(l,tlo,thi,grp_ng,grp_wl,f)
     grd_emnbf(:,l) = grd_emnbf(:,l) + real(p*f/(4d0*pc_pi*grd_vol(l))*pc_pi*grd_tempinv(l)**4/pc_sb)
  enddo
  deallocate(f)
end subroutine shock_emnbf
!
subroutine write_emitchan
!-- escaped luminosity by last emission process and element, per band
  integer :: k,z,b
  character(12),parameter :: kn(0:6) = [character(12) :: 'unknown','inner_bound', &
       'thermal','bound-free','line','two-photon','shock']
  open(4,file='output.emitchan',status='replace')
  write(4,'(a)') '# escaped luminosity [erg/s] by last emission process and element (Z), per band [A]'
  write(4,'(a12,a4,9a11)') '# process','Z','<912','912-3e3','3e3-5e3','5e3-7e3','7e3-1e4', &
       '1e4-2.5e4','2.5e4-1e5','>1e5','total'
  do k=0,6
     do z=0,30
        if(sum(flx_emch(k,z,:))<=0d0) cycle
        write(4,'(a12,i4,1p,9e11.3)') kn(k),z,(flx_emch(k,z,b),b=1,flx_nemb),sum(flx_emch(k,z,:))
     enddo
  enddo
  close(4)
end subroutine write_emitchan


subroutine atomfiles(str,n,files)
!-- split a blank-separated list of model-atom files
  character(*),intent(in) :: str
  integer,intent(out) :: n
  character(*),intent(out) :: files(:)
  integer :: k,j
  n = 0
  k = 1
  do while(k <= len_trim(str))
     if(str(k:k)==' ' .or. str(k:k)==',') then
        k = k + 1
        cycle
     endif
     j = k
     do while(j <= len_trim(str))
        if(str(j:j)==' ' .or. str(j:j)==',') exit
        j = j + 1
     enddo
     n = n + 1
     if(n > size(files)) stop 'in_nltex_atoms: too many files'
     files(n) = str(k:j-1)
     k = j
  enddo
end subroutine atomfiles

end program superlite
! vim: fdm=marker
