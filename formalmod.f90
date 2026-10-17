!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module formalmod
!-----------------------------------------------------------------------
! Formal integral of the emergent spectrum, as CMFGEN's cmf_flux, after
! the last iteration (in_fi).  The observer-frame intensity is integrated
! along rays of impact parameter p through the spherical grid, with the
! source functions of the last iteration:
!  - continuum (per group of the cell): absorption and thermal emission of
!    the LTE opacities (emission opacity times B(T_e)), the recombination
!    emissivity of the NLTE levels, electron scattering of the Monte Carlo
!    J redistributed by the thermal electron motions (the angle-averaged
!    kernel of thermal Thomson scattering, isotropic or dipole as in the
!    transport: esk), the scattering part of the expansion opacity of the
!    LTE lines (J);
!  - NLTE lines with the source function of the level populations
!    (complete redistribution) and Doppler profiles, the thermal width of
!    the ion and the microturbulence (in_nltex_vturb, linear in v up to
!    in_nltex_vturbmax): K/(sqrt(pi) v_D) exp(-x^2) per unit length, which
!    is the Sobolev depth K/|Q(mu)| where the flow crosses the profile
!    quickly, and in a slow flow (a coasting wind: the narrow P-Cygni
!    absorption) the depth of the profile over the path (as CMFGEN's
!    cmf_flux); or (dop absent/false) Sobolev resonances where the comoving
!    frequency crosses the line (several in a non-monotonic flow);
!  - with nxt_estrap, the line photons scattered by electrons while trapped
!    in their resonance zones (nx_bes), n_u A bes h nu/(4 pi) per line
!    with the populations of the iterations (as CMFGEN's cmf_flux adds the
!    incoherent scattering to the populations of iterations with coherent
!    scattering), redistributed by the thermal electrons as the electron
!    scattering of J, on a grid of the comoving frequency (30 km/s or the
!    spacing of the spectrum).  These photons are not in the Monte Carlo J,
!    so their next scattering is added to first order: the mean intensity
!    of the optically thin shells of their emission, J(r) = L/(32 pi^2 r r')
!    ln|(r+r')/(r-r')| (no Doppler shifts between the shells), scattered
!    with the same kernel (rescat).
! The invariant I/nu^3 is integrated with the comoving source functions
! (first order in v/c).  Rays through the inner boundary start with the
! Monte Carlo J of the innermost cell.  L_nu = 8 pi^2 int I p dp on a
! grid uniform in ln(wavelength); the frequencies are shared out among
! the MPI ranks.
!-----------------------------------------------------------------------
  implicit none
  private
  public :: fi_spectrum
  real*8,parameter :: pi = 3.141592653589793d0, cc = 2.99792458d10, hh = 6.6260693d-27
  real*8,parameter :: kb = 1.3806505d-16, me = 9.1093826d-28, sb = 5.670400d-5
!-- thermal electron scattering of a monochromatic beam, averaged over the
!-- scattering angle, against u = |dnu|/(nu sqrt(2kT/m_e c^2)) on [0, ukmax]
  integer,parameter :: nkt = 350
  real*8,parameter :: ukmax = 7d0
  real*8 :: ktab(0:nkt)
!
contains
!
!
subroutine fi_spectrum(ncell,xarr,vxarr,tempinv,wl1,wl2,dvkms,npc,impi,nmpi,lmpi0,jgavg,dop,nes)
  use gridmod, only:grd_cap,grd_sig,grd_capemit
  use bbsobmod, only:grd_capsc
  use groupmod, only:grp_ng,grp_wl
  use nltextrmod, only:nxt_jg,nxt_nug,nxt_dnug,grd_emnbf,nxt_estrap,nxt_vturbv
  use nlteatommod, only:nxa_nel
  use nltexmod, only:nx_vdop
  use mpimod, only:allreduce_sum
  use transportmod, only:trn_thm_aniso
  integer,intent(in) :: ncell,npc,impi,nmpi
  real*8,intent(in) :: xarr(ncell+1),vxarr(ncell+1),tempinv(ncell),wl1,wl2,dvkms
  logical,intent(in) :: lmpi0
  real*8,intent(in),optional :: jgavg(:,:) !(ng,ncell) J of several iterations (else nxt_jg)
  logical,intent(in),optional :: dop       !Doppler line profiles (else Sobolev resonances)
  integer,intent(in),optional :: nes       !passes of the electron-scattering J (else the Monte Carlo J)
  integer :: nnu,i,j,ic,ig,ng,np,k,jg,nt,nesp,nj,mx,ipass,jj
  real*8 :: dln,t,x,bnu,se,sig2,wsum,lnu0,dl,lt0,dlt,lext,lj0,wdep,dmax,lprev,nga
  real*8,allocatable :: wl(:),lnu(:),chi(:,:),eta(:,:),p(:),wp(:),lng(:),etes(:,:),jgx(:,:)
  real*8,allocatable :: vdc(:,:),vdmax(:),jfi(:,:),jnw(:,:),etj(:,:),vol(:),jold(:,:),jmc0(:,:)
  real*4,allocatable :: etr(:,:)
  logical :: ldop,lj,ldep
!
  ng = grp_ng
!-- Doppler widths of the lines per element and cell: thermal and the
!-- microturbulence
  ldop = .false.
  if(present(dop)) ldop = dop
  nesp = 0
  if(present(nes)) nesp = max(nes,0)
  lj = nesp > 0
  ldep = .false.
  wdep = 0d0
  allocate(vdc(nxa_nel,ncell),vdmax(ncell))
  do ic=1,ncell
     x = nxt_vturbv(0.5d0*(vxarr(ic)+vxarr(ic+1)))
     do k=1,nxa_nel
        vdc(k,ic) = nx_vdop(k,1d0/tempinv(ic),x)
     enddo
     vdmax(ic) = maxval(vdc(:,ic))
  enddo
!-- observer wavelengths
  dln = dvkms*1d5/cc
  nnu = int(log(wl2/wl1)/dln) + 1
  allocate(wl(nnu),lnu(nnu))
  do i=1,nnu
     wl(i) = wl1*exp((i-1)*dln)
  enddo
  call es_kernel(trn_thm_aniso)
  allocate(jgx(ng,ncell))
  if(present(jgavg)) then
     jgx = jgavg
  else
     jgx = nxt_jg
  endif
!-- continuum opacity and emissivity per group and cell (comoving frame);
!-- the electron scattering of J (etes) is smooth in frequency and is
!-- interpolated between the group centres (etesi), the rest is constant
!-- within a group (edges at group boundaries stay sharp)
  allocate(chi(ng,ncell),eta(ng,ncell),lng(ng),etes(ng,ncell))
  lng = log(nxt_nug)
  do ic=1,ncell
     t = 1d0/tempinv(ic)
!-- electron scattering of J redistributed in frequency (thermal electrons)
     se = sqrt(2d0*kb*t/me)/cc
     do ig=1,ng
        x = hh*nxt_nug(ig)/(kb*t)
        bnu = 0d0
        if(x < 700d0) bnu = 2d0*hh*nxt_nug(ig)**3/cc**2/em1(x)
        chi(ig,ic) = grd_cap(ig,ic) + grd_sig(ic)
        eta(ig,ic) = grd_capemit(ig,ic)*bnu + grd_emnbf(ig,ic)*sb*t**4/pi/nxt_dnug(ig)
        if(allocated(grd_capsc)) then
           chi(ig,ic) = chi(ig,ic) + grd_capsc(ig,ic)
           eta(ig,ic) = eta(ig,ic) + grd_capsc(ig,ic)*jgx(ig,ic)
        endif
!-- J averaged over the redistribution kernel
        sig2 = 0d0
        wsum = 0d0
        do jg=ig,1,-1
           dl = (lng(jg) - lng(ig))/se
           if(dl > ukmax) exit
           x = esk(dl)*nxt_dnug(jg)/nxt_nug(jg)
           sig2 = sig2 + x*jgx(jg,ic)
           wsum = wsum + x
        enddo
        do jg=ig+1,ng
           dl = (lng(ig) - lng(jg))/se
           if(dl > ukmax) exit
           x = esk(dl)*nxt_dnug(jg)/nxt_nug(jg)
           sig2 = sig2 + x*jgx(jg,ic)
           wsum = wsum + x
        enddo
        etes(ig,ic) = 0d0
        if(wsum > 0d0) etes(ig,ic) = grd_sig(ic)*sig2/wsum
     enddo
  enddo
!-- electron scattering of the trapped line photons [erg/s/cm^3/Hz/sr] at
!-- comoving ln(nu) = lt0 + (j-1) dlt
  nt = 0
  dlt = max(dln,1d-4)
  lext = 0.02d0 + ukmax*maxval(sqrt(2d0*kb/(me*tempinv)))/cc + maxval(abs(vxarr))/cc
  lt0 = log(cc/wl2) - lext
  if(nxt_estrap .and. .not.lj) then
     nt = int((log(cc/wl1) + lext - lt0)/dlt) + 2
     allocate(etr(nt,ncell))
     call trapped_es(etr)
     call rescat(etr)
!-- (its luminosity in the table: 4 pi Int eta dnu dV)
     if(lmpi0) then
        x = 0d0
        do ic=1,ncell
           x = x + 4d0*pi*4d0*pi/3d0*(xarr(ic+1)**3-xarr(ic)**3)* &
                sum(dble(etr(:,ic))*exp(lt0 + (/(j-1,j=1,nt)/)*dlt))*dlt
        enddo
        write(6,'(1x,a,es11.4,a)') 'formal integral: electron scattering of trapped line photons',x,' erg/s'
     endif
  endif
!-- impact parameters: core rays uniform in p^2, then two per cell
  np = npc + 2*ncell
  allocate(p(np),wp(np))
  do j=1,npc
     p(j) = xarr(1)*sqrt((j-0.5d0)/npc)
  enddo
  do ic=1,ncell
     p(npc+2*ic-1) = xarr(ic) + 0.25d0*(xarr(ic+1)-xarr(ic))
     p(npc+2*ic) = xarr(ic) + 0.75d0*(xarr(ic+1)-xarr(ic))
  enddo
!-- trapezoid weights of int I p dp (I -> 0 at the outer edge)
  wp = 0d0
  do j=1,np
     if(j>1) wp(j) = wp(j) + 0.5d0*(p(j)-p(j-1))*p(j)
     if(j<np) then
        wp(j) = wp(j) + 0.5d0*(p(j+1)-p(j))*p(j)
     else
        wp(j) = wp(j) + 0.5d0*(xarr(ncell+1)-p(j))*p(j)
     endif
  enddo
  wp(1) = wp(1) + 0.5d0*p(1)*p(1)
!-- the frequencies of this rank
  lnu = 0d0
  if(.not.lj) then
     do i=impi+1,nnu,nmpi
        lnu0 = 0d0
        do j=1,np
           lnu0 = lnu0 + wp(j)*ray(p(j),cc/wl(i))
        enddo
        lnu(i) = 8d0*pi*pi*lnu0
     enddo
  else
!-- electron scattering with the J of the formal solution (as CMFGEN's
!-- cmf_flux): comoving J per cell on the grid of the spectrum, extended by
!-- the kernel and the flow, from the Monte Carlo J; each pass integrates
!-- the rays of all its frequencies with the emissivity of the previous J
!-- and tallies the new J along them (parallel rays sample every direction
!-- of a spherical shell: J = 1/V sum_p 2 pi p dp Int I dz)
     mx = int(lext/dln) + 1
     nj = nnu + 2*mx
     lj0 = log(cc/wl(nnu)) - mx*dln
     allocate(jfi(nj,ncell),jnw(nj,ncell),etj(nj,ncell),vol(ncell),jold(nj,ncell))
!-- the volume of a cell as the rays sample it (sum_p 2 pi p dp of their
!-- chords on both sides), so that a uniform I gives J = I: the shell volume
!-- differs by 0.2-7 % (70 % in a thin innermost cell), an excess that a
!-- scattering wind amplifies pass by pass
     do ic=1,ncell
        vol(ic) = 0d0
        do j=1,np
           vol(ic) = vol(ic) + 4d0*pi*wp(j)*(sqrt(max(xarr(ic+1)**2-p(j)**2,0d0)) - &
                sqrt(max(xarr(ic)**2-p(j)**2,0d0)))
        enddo
     enddo
     call jinit
     allocate(jmc0(nj,ncell))
     jmc0 = jfi
     lprev = 0d0
     do ipass=1,nesp
        call jemis
        jnw = 0d0
        lnu = 0d0
        ldep = .true.
        do jj=impi+1,nj,nmpi
           lnu0 = 0d0
           do j=1,np
              wdep = 2d0*pi*wp(j)
              lnu0 = lnu0 + wp(j)*ray(p(j),exp(lj0 + (jj-1)*dln))
           enddo
           i = nnu + mx + 1 - jj
           if(i >= 1 .and. i <= nnu) lnu(i) = 8d0*pi*pi*lnu0
        enddo
        ldep = .false.
        call allreduce_sum(nj*ncell,jnw)
        call allreduce_sum(nnu,lnu)
        dmax = 0d0
        do ic=1,ncell
           jnw(:,ic) = jnw(:,ic)/vol(ic)
           x = sum(jfi(mx+1:mx+nnu,ic))
           if(x > 0d0) dmax = max(dmax,abs(sum(jnw(mx+1:mx+nnu,ic))/x - 1d0))
        enddo
!-- every third pass Ng's acceleration from the last three J (the passes
!-- are a linear iteration; in a thick scattering wind it converges slowly)
        nga = 0d0
        if(mod(ipass,3) == 0 .and. ipass < nesp) call ngacc(nga)
        jold = jfi
        jfi = jnw
        x = sum(lnu(2:)*(cc/wl(1:nnu-1)-cc/wl(2:nnu)))
        if(lmpi0) write(6,'(1x,a,i3,a,es11.4,a,f8.4,a,es9.2,a,f8.3)') 'formal integral: e.s. pass',ipass, &
             ', L =',x,' erg/s, change',x/max(lprev,1d-300)-1d0,', largest change of a cell J',dmax, &
             ', Ng',nga
        lprev = x
        if(ipass < nesp) lnu = 0d0
     enddo
!-- J of the last pass against the Monte Carlo J per cell, Int J dnu over
!-- the comoving frequencies of the spectrum and its thirds
     if(lmpi0) then
        open(4,file='output.fi_jes',status='replace')
        write(4,'(a)') '# formal integral (in_fi_nes): Int J dnu of the last pass / Monte Carlo, per cell'
        write(4,'(a)') '# cell  r_centre[cm]  all  blue third  middle third  red third'
        do ic=1,ncell
           write(4,'(i5,es13.5,4f9.4)') ic,0.5d0*(xarr(ic)+xarr(ic+1)),jratio(ic,mx+1,mx+nnu), &
                jratio(ic,mx+2*nnu/3+1,mx+nnu),jratio(ic,mx+nnu/3+1,mx+2*nnu/3),jratio(ic,mx+1,mx+nnu/3)
        enddo
        close(4)
     endif
     deallocate(jfi,jnw,etj,vol,jold,jmc0)
  endif
  if(.not.lj) call allreduce_sum(nnu,lnu)
!-- L_lambda [erg/s/A]
  if(lmpi0) then
     open(4,file='output.fi_spectrum',status='replace')
     write(4,'(a)') '# formal integral (in_fi): wavelength [A], L_lambda [erg/s/A], L_nu [erg/s/Hz]'
     write(4,'(a,i8,a,f8.2,a,i5,a)') '#',nnu,' points,',dvkms,' km/s spacing,',np,' rays'
     do i=1,nnu
        write(4,'(f14.4,2es14.6)') wl(i)*1d8,lnu(i)*cc/wl(i)**2*1d-8,lnu(i)
     enddo
     close(4)
     write(6,'(1x,a,es12.4,a)') 'formal integral: L =',sum(lnu(2:)*(cc/wl(1:nnu-1)-cc/wl(2:nnu))), &
          ' erg/s in the range (output.fi_spectrum)'
  endif
  deallocate(wl,lnu,chi,eta,lng,p,wp,etes,jgx,vdc,vdmax)
  if(allocated(etr)) deallocate(etr)
!
contains
!
  real*8 function jratio(c,ja,jb)
!-- Int J dnu of the passes over Int of the Monte Carlo J, J grid points ja..jb
    integer,intent(in) :: c,ja,jb
    real*8 :: s1,s2,w
    integer :: jl
    s1 = 0d0
    s2 = 0d0
    do jl=ja,jb
       w = exp(lj0 + (jl-1)*dln)
       s1 = s1 + w*jfi(jl,c)
       s2 = s2 + w*jmc0(jl,c)
    enddo
    jratio = 0d0
    if(s2 > 0d0) jratio = s1/s2
  end function jratio
!
  subroutine trapped_es(et)
!-- emissivity of the line photons scattered by electrons in the resonance
!-- zones: the thermal-electron kernel (as the e.s. of J above) about the
!-- line, per cell (volume mean: microclumping)
    use nltextrmod, only:nxt_bes,nxt_x,nxt_nel,nxt_w,nxt_fcl
    use nlteatommod
    real*4,intent(out) :: et(nt,ncell)
    integer :: c,kl,iel,o,j0,j1,jj
    real*8 :: em,emax,sel,lk,ww,dk
    et = 0d0
    do c=1,ncell
       sel = sqrt(2d0*kb/(me*tempinv(c)))/cc
       emax = 0d0
       do kl=1,nxa_nlin
          if(nxt_bes(kl,c) <= 0d0) cycle
          iel = nxa_lel(kl)
          o = nxa_el(iel)%xoff
!-- n_u A bes h nu / (4 pi)
          em = nxt_nel(iel,c)*nxt_x(o+nxa_lsu(kl),c)*nxt_w(nxa_luv(kl),c)* &
               nxa_la(kl)*nxt_bes(kl,c)*hh*nxa_lnu(kl)/(4d0*pi)
          if(allocated(nxt_fcl)) em = em*nxt_fcl(c)
          if(em <= 1d-8*emax) cycle
          emax = max(emax,em)
          lk = log(nxa_lnu(kl))
          j0 = max(1,int((lk - ukmax*sel - lt0)/dlt) + 1)
          j1 = min(nt,int((lk + ukmax*sel - lt0)/dlt) + 2)
          do jj=j0,j1
             dk = abs(lt0 + (jj-1)*dlt - lk)/sel
             ww = esk(dk)/(sel*exp(lt0 + (jj-1)*dlt))
             et(jj,c) = et(jj,c) + real(em*ww)
          enddo
       enddo
    enddo
  end subroutine trapped_es
!
  real*8 function ray(pp,nu)
!-- observer-frame intensity I_nu of the ray with impact parameter pp
    real*8,intent(in) :: pp,nu
    integer :: ie,ne
    real*8 :: zs(2*ncell+2),ui
    ui = 0d0
!-- crossings with the shell edges, from the back to the front
    ne = 0
    if(pp < xarr(1)) then
!-- core ray: in a pass of the electron-scattering J also its part behind
!-- the inner boundary, toward it (only tallied: those directions of J)
       if(ldep) then
          do ie=ncell+1,1,-1
             ne = ne + 1
             zs(ne) = -sqrt(xarr(ie)**2 - pp*pp)
          enddo
          call segments(ui,pp,nu,zs,ne)
          ui = 0d0
          ne = 0
       endif
!-- from the inner boundary on the observer's side
       ne = 1
       zs(1) = sqrt(xarr(1)**2 - pp*pp)
       ui = core_i(zs(1)/xarr(1),nu)
       do ie=2,ncell+1
          ne = ne + 1
          zs(ne) = sqrt(xarr(ie)**2 - pp*pp)
       enddo
    else
       do ie=ncell+1,1,-1
          if(xarr(ie) <= pp) exit
          ne = ne + 1
          zs(ne) = -sqrt(xarr(ie)**2 - pp*pp)
       enddo
       ne = ne + 1
       zs(ne) = 0d0
       do ie=1,ncell+1
          if(xarr(ie) <= pp) cycle
          ne = ne + 1
          zs(ne) = sqrt(xarr(ie)**2 - pp*pp)
       enddo
    endif
    call segments(ui,pp,nu,zs,ne)
    ray = ui*nu**3
  end function ray
!
  subroutine segments(ui,pp,nu,zs,ne)
!-- the intervals between the crossings zs (increasing z), each inside one
!-- cell, in substeps by the change of the line-of-sight velocity (50 km/s)
    real*8,intent(inout) :: ui
    real*8,intent(in) :: pp,nu,zs(ne)
    integer,intent(in) :: ne
    integer :: iz,c,ns,is
    real*8 :: za,zb,zm,ra,rb,rm,dvdr,va,vb,nua,nub,vlos(2)
    do iz=1,ne-1
       za = zs(iz)
       zb = zs(iz+1)
       if(zb <= za) cycle
       zm = 0.5d0*(za+zb)
       rm = sqrt(pp*pp + zm*zm)
       c = cellof(rm)
       dvdr = (vxarr(c+1)-vxarr(c))/(xarr(c+1)-xarr(c))
       ra = sqrt(pp*pp + za*za)
       rb = sqrt(pp*pp + zb*zb)
       vlos(1) = za/max(ra,1d-300)*vel(c,ra)
       vlos(2) = zb/max(rb,1d-300)*vel(c,rb)
       ns = max(2,min(200,int(abs(vlos(2)-vlos(1))/5d6)+1))
       do is=1,ns
          va = za + (zb-za)*(is-1)/dble(ns)
          vb = za + (zb-za)*is/dble(ns)
          ra = sqrt(pp*pp + va*va)
          rb = sqrt(pp*pp + vb*vb)
          nua = nu*(1d0 - va/max(ra,1d-300)*vel(c,ra)/cc)
          nub = nu*(1d0 - vb/max(rb,1d-300)*vel(c,rb)/cc)
          call substep(ui,nu,pp,c,va,vb,nua,nub,dvdr)
       enddo
    enddo
  end subroutine segments
!
  subroutine rescat(et)
!-- first-order electron scattering of the trapped-photon emission et:
!-- J of each cell from all cells' emission as optically thin shells, each
!-- weighted by exp(-tau_es) to the outer edge (the emission of the optically
!-- thick layers stays there), then sigma times J redistributed with the
!-- thermal-electron kernel, added to et
    real*4,intent(inout) :: et(nt,ncell)
    real*8,allocatable :: jt(:,:),dd(:,:)
    real*8 :: rc(ncell),ra,rb,lj,sel,u,tauo(ncell+1)
    integer :: c,c2,jj,j2,m
    allocate(jt(nt,ncell),dd(ncell,ncell))
    tauo(ncell+1) = 0d0
    do c=ncell,1,-1
       rc(c) = 0.5d0*(xarr(c)+xarr(c+1))
       tauo(c) = tauo(c+1) + grd_sig(c)*(xarr(c+1)-xarr(c))
    enddo
!-- J per unit luminosity of a shell [1/cm^2/sr]; within a cell the mean
!-- over two sub-shells (the logarithm diverges at r = r')
    do c=1,ncell
       do c2=1,ncell
          dd(c,c2) = 0d0
          do m=1,2
             ra = rc(c)
             rb = xarr(c2) + (0.25d0 + 0.5d0*(m-1))*(xarr(c2+1)-xarr(c2))
             dd(c,c2) = dd(c,c2) + 0.5d0*log(abs((ra+rb)/(ra-rb)))/(32d0*pi*pi*ra*rb)
          enddo
       enddo
    enddo
!-- J of the trapped emission: shell luminosity 4 pi V eta per Hz
    jt = 0d0
    do c2=1,ncell
       ra = 4d0*pi*4d0*pi/3d0*(xarr(c2+1)**3 - xarr(c2)**3)*exp(-0.5d0*(tauo(c2)+tauo(c2+1)))
       if(ra <= 0d0) cycle
       do c=1,ncell
          jt(:,c) = jt(:,c) + dd(c,c2)*ra*dble(et(:,c2))
       enddo
    enddo
!-- its electron scattering with the thermal-electron kernel (bins of
!-- width dlt in ln nu; the kernel per unit u = dln(nu)/sel)
    do c=1,ncell
       sel = sqrt(2d0*kb/(me*tempinv(c)))/cc
       m = int(ukmax*sel/dlt) + 1
       do jj=1,nt
          lj = 0d0
          do j2=max(1,jj-m),min(nt,jj+m)
             if(jt(j2,c) <= 0d0) cycle
             u = abs(jj-j2)*dlt/sel
             lj = lj + esk(u)*jt(j2,c)
          enddo
          et(jj,c) = et(jj,c) + real(grd_sig(c)*lj*dlt/sel)
       enddo
    enddo
    deallocate(jt,dd)
  end subroutine rescat
!
  subroutine es_kernel(ldip)
!-- the kernel table: a scattering angle theta shifts by u ~ N(0, t),
!-- t = 1-cos theta (u in units of sqrt(2kT/m_e c^2)), t uniform on [0,2]
!-- (isotropic) or weighted by 1+cos^2 theta (dipole); midpoint rule in
!-- s = sqrt(t)
    logical,intent(in) :: ldip
    integer,parameter :: ns = 400
    integer :: j,k
    real*8 :: u,sq,p
    do j=0,nkt
       u = ukmax*j/dble(nkt)
       ktab(j) = 0d0
       do k=1,ns
          sq = sqrt(2d0)*(k-0.5d0)/ns
          if(ldip) then
             p = 3d0/8d0*(1d0 + (1d0-sq*sq)**2)
          else
             p = 0.5d0
          endif
          ktab(j) = ktab(j) + p*exp(-u*u/(2d0*sq*sq))
       enddo
       ktab(j) = ktab(j)*sqrt(2d0/pi)*sqrt(2d0)/ns
    enddo
  end subroutine es_kernel
!
  real*8 function esk(u)
!-- the kernel at u >= 0 (per unit u; normalized over -inf < u < inf)
    real*8,intent(in) :: u
    real*8 :: w
    integer :: j
    esk = 0d0
    if(u >= ukmax) return
    w = u/ukmax*nkt
    j = int(w)
    w = w - j
    esk = (1d0-w)*ktab(j) + w*ktab(j+1)
  end function esk
!
  real*8 function etesi(nuc,g,c)
!-- electron-scattering emissivity of J at comoving frequency nuc in group g
!-- of cell c, linear in ln(nu) between the group centres
    real*8,intent(in) :: nuc
    integer,intent(in) :: g,c
    real*8 :: l,w
    integer :: g2
    l = log(nuc)
    if(l < lng(g)) then
       g2 = min(g+1,ng)       !(groups in order of wavelength)
    else
       g2 = max(g-1,1)
    endif
    etesi = etes(g,c)
    if(g2 == g) return
    w = (l - lng(g))/(lng(g2) - lng(g))
    etesi = (1d0-w)*etes(g,c) + w*etes(g2,c)
  end function etesi
!
  real*8 function etrap(nuc,c)
!-- trapped-photon e.s. emissivity at comoving frequency nuc in cell c
    real*8,intent(in) :: nuc
    integer,intent(in) :: c
    real*8 :: w
    integer :: j
    etrap = 0d0
    if(nt <= 0) return
    w = (log(nuc) - lt0)/dlt
    j = int(w) + 1
    if(j < 1 .or. j >= nt) return
    w = w - (j-1)
    etrap = (1d0-w)*etr(j,c) + w*etr(j+1,c)
  end function etrap
!
  real*8 function core_i(mu,nu)
!-- invariant intensity I/nu^3 leaving the inner boundary: the Monte Carlo
!-- J of the innermost cell at the comoving frequency
    real*8,intent(in) :: mu,nu
    real*8 :: nuc
    integer :: g
    nuc = nu*(1d0 - mu*vxarr(1)/cc)
    g = groupof(nuc)
    core_i = 0d0
    if(g>0) core_i = jgx(g,1)/nuc**3
  end function core_i
!
  subroutine substep(ui,nu,pp,c,za,zb,nua,nub,dvdr)
!-- continuum and the NLTE line resonances between comoving frequencies
!-- nua (at za) and nub (at zb) in cell c
    use nltextrmod, only:nxt_nact,nxt_act,nxt_wlact,nxt_kq
    real*8,intent(inout) :: ui
    real*8,intent(in) :: nu,pp,za,zb,nua,nub,dvdr
    integer,intent(in) :: c
    integer :: g,ja,jb,j,jstep,k
    real*8 :: tprev,tk,ds,ch,sc,num,dtau,zk,rk,mu,q,tau,sl,nul
    num = 0.5d0*(nua+nub)
    g = groupof(num)
    ch = 0d0
    sc = 0d0
    if(g>0) then
       ch = chi(g,c)*num/nu
       if(chi(g,c) > 0d0) sc = (eta(g,c) + esem(num,g,c) + etrap(num,c))/chi(g,c)/num**3
    endif
    if(ldop) then
       call substep_dop(ui,nu,c,za,zb,nua,nub,ch,sc)
       return
    endif
!-- active lines with comoving frequency between nua and nub
    ja = 0
    jb = -1
    if(nxt_nact>0 .and. nua/=nub) then
       ja = lineidx(cc/max(nua,nub),.true.)
       jb = lineidx(cc/min(nua,nub),.false.)
    endif
    tprev = 0d0
    if(jb >= ja .and. ja >= 1) then
!-- frequency falling along the ray: increasing wavelength
       if(nub < nua) then
          j = ja
          jstep = 1
       else
          j = jb
          jstep = -1
       endif
       do
          k = nxt_act(j)
          nul = cc/nxt_wlact(j)
          tk = (nul - nua)/(nub - nua)
          tk = min(max(tk,tprev),1d0)
          ds = (tk - tprev)*(zb - za)
          call contd(ui,ch*ds,sc,c,nua+(nub-nua)*0.5d0*(tprev+tk),ds)
!-- Sobolev resonance
          zk = za + tk*(zb - za)
          rk = sqrt(pp*pp + zk*zk)
          mu = zk/max(rk,1d-300)
          q = abs(mu*mu*dvdr + (1d0-mu*mu)*vel(c,rk)/rk)
          tau = nxt_kq(k,c)/max(q,1d-30)
          if(tau > 0d0) then
             sl = srcline(k,c,nul)
             if(tau > 50d0) then
                ui = sl/nul**3
             else
                ui = ui*exp(-tau) + sl/nul**3*(1d0 - exp(-tau))
             endif
          endif
          tprev = tk
          if(j == merge(jb,ja,jstep==1)) exit
          j = j + jstep
       enddo
    endif
    ds = (1d0 - tprev)*(zb - za)
    call contd(ui,ch*ds,sc,c,nua+(nub-nua)*0.5d0*(tprev+1d0),ds)
  end subroutine substep
!
  subroutine substep_dop(ui,nu,c,za,zb,nua,nub,ch,sc)
!-- the substep with Doppler line profiles: opacity K/(sqrt(pi) v_D)
!-- exp(-x^2) per unit length, x = (nu_cmf/nu_l - 1) c/v_D (over the whole
!-- profile K/|Q|, the Sobolev depth), and the source function of the level
!-- populations (complete redistribution).  The comoving frequency is linear
!-- along the substep, so the depth of a line over a piece is exact (erf).
!-- The pieces march along the substep: within the profile of a line
!-- (|x| < xpmax) at most dxp in its x, and they stop where a profile
!-- begins; line and continuum share the source function within a piece.
!-- A line whose profile the ray crosses within a tenth of the cell width
!-- (path per Doppler width ds/|dx|: the Sobolev length) is a point
!-- resonance as in substep (fast flows: the ejecta of a supernova).
    use nltextrmod, only:nxt_nact,nxt_act,nxt_wlact,nxt_kq
    use nlteatommod, only:nxa_lel
    real*8,intent(inout) :: ui
    real*8,intent(in) :: nu,za,zb,nua,nub,ch,sc
    integer,intent(in) :: c
    integer,parameter :: mxl = 400
    real*8,parameter :: xpmax = 4.5d0, dxp = 0.5d0, fsob = 0.1d0
    integer :: ja,jb,j,k,nl,il,ns,is,m
    real*8 :: numin,numax,nul,vd,xa,xb,ds,t0,t1,dtc,dtl,st,x0,x1,w,dcell
    real*8 :: lxa(mxl),lxb(mxl),lsl(mxl),lkv(mxl),lta(mxl),ltb(mxl),ldt(mxl)
    real*8 :: stk(mxl),sttau(mxl),stsl(mxl)
    numin = min(nua,nub)
    numax = max(nua,nub)
    ds = zb - za
    dcell = xarr(c+1) - xarr(c)
    nl = 0
    ns = 0
    if(nxt_nact>0) then
       ja = lineidx(cc/(numax*(1d0 + xpmax*vdmax(c)/cc)),.true.)
       jb = lineidx(cc/(numin*(1d0 - xpmax*vdmax(c)/cc)),.false.)
       do j=max(ja,1),jb
          k = nxt_act(j)
          if(nxt_kq(k,c) <= 0d0) cycle
          nul = cc/nxt_wlact(j)
          vd = vdc(nxa_lel(k),c)
          xa = (nua/nul - 1d0)*cc/vd
          xb = (nub/nul - 1d0)*cc/vd
          if(min(xa,xb) >= xpmax .or. max(xa,xb) <= -xpmax) cycle
!-- (short Sobolev length: a point resonance where x = 0, depth K ds/(v_D |dx|))
          if(abs(xb - xa) > 1d-12) then
             if(ds/abs(xb - xa) < fsob*dcell) then
                w = -xa/(xb - xa)
                if(w < 0d0 .or. w >= 1d0 .or. ns == mxl) cycle
                ns = ns + 1
                stk(ns) = w
                sttau(ns) = nxt_kq(k,c)*ds/(vd*abs(xb - xa))
                stsl(ns) = srcline(k,c,nul)/nul**3
                cycle
             endif
          endif
          if(nl == mxl) cycle
          nl = nl + 1
          lxa(nl) = xa
          lxb(nl) = xb
          lsl(nl) = srcline(k,c,nul)/nul**3
          lkv(nl) = nxt_kq(k,c)/vd
!-- (the part of the substep within the profile, and the step in it)
          if(abs(xb - xa) > 1d-12) then
             x0 = (-xpmax - xa)/(xb - xa)
             x1 = (xpmax - xa)/(xb - xa)
             lta(nl) = max(min(x0,x1),0d0)
             ltb(nl) = min(max(x0,x1),1d0)
             ldt(nl) = dxp/abs(xb - xa)
          else
             lta(nl) = 0d0
             ltb(nl) = 1d0
             ldt(nl) = 1d0
          endif
       enddo
    endif
    if(nl == 0 .and. ns == 0) then
       call contd(ui,ch*ds,sc,c,0.5d0*(nua+nub),ds)
       return
    endif
!-- the point resonances in the order along the substep
    do is=2,ns
       w = stk(is); x0 = sttau(is); x1 = stsl(is)
       m = is - 1
       do while(m >= 1)
          if(stk(m) <= w) exit
          stk(m+1) = stk(m); sttau(m+1) = sttau(m); stsl(m+1) = stsl(m)
          m = m - 1
       enddo
       stk(m+1) = w; sttau(m+1) = x0; stsl(m+1) = x1
    enddo
    is = 1
    t0 = 0d0
    do while(t0 < 1d0)
!-- the end of this piece
       t1 = 1d0
       if(is <= ns) t1 = stk(is)
       do il=1,nl
          if(ltb(il) <= t0) cycle
          if(lta(il) > t0) then
             t1 = min(t1,lta(il))
          else
             t1 = min(t1,t0+ldt(il),ltb(il))
          endif
       enddo
       if(t1 > t0) then
          dtc = ch*(t1 - t0)*ds
          dtl = 0d0
          st = dtc*sc
          do il=1,nl
             if(ltb(il) <= t0 .or. lta(il) >= t1) cycle
             x0 = lxa(il) + (lxb(il) - lxa(il))*t0
             x1 = lxa(il) + (lxb(il) - lxa(il))*t1
             if(abs(x1 - x0) > 1d-6) then
                w = 0.5d0*(erf(x1) - erf(x0))/(x1 - x0)
             else
                w = exp(-(0.5d0*(x0+x1))**2)/sqrt(pi)
             endif
!-- (opacity in the observer frame: times nu_cmf/nu)
             w = lkv(il)*w*(t1 - t0)*ds*(nua + (nub - nua)*0.5d0*(t0 + t1))/nu
             dtl = dtl + w
             st = st + w*lsl(il)
          enddo
          if(dtc + dtl > 0d0) call contd(ui,dtc+dtl,st/(dtc+dtl),c,nua+(nub-nua)*0.5d0*(t0+t1),(t1-t0)*ds)
       endif
!-- point resonances at the end of the piece
       do while(is <= ns)
          if(stk(is) > t1) exit
          if(sttau(is) > 50d0) then
             ui = stsl(is)
          else
             ui = ui*exp(-sttau(is)) + stsl(is)*(1d0 - exp(-sttau(is)))
          endif
          is = is + 1
       enddo
       t0 = max(t1,min(t0+1d-9,1d0))
    enddo
  end subroutine substep_dop
!
  subroutine cont(ui,dtau,sc)
    real*8,intent(inout) :: ui
    real*8,intent(in) :: dtau,sc
    if(dtau <= 0d0) return
    if(dtau > 50d0) then
       ui = sc
    elseif(dtau < 1d-8) then
       ui = ui + (sc - ui)*dtau
    else
       ui = ui*exp(-dtau) - sc*em1(-dtau)
    endif
  end subroutine cont
!
  subroutine contd(ui,dtau,sc,c,nuc,ds)
!-- cont, and in a pass of the electron-scattering J the tally of the mean
!-- intensity over the piece (length ds, comoving frequency nuc, cell c)
!-- into the two nearest points of the J grid
    real*8,intent(inout) :: ui
    real*8,intent(in) :: dtau,sc,nuc,ds
    integer,intent(in) :: c
    real*8 :: um,w
    integer :: jl
    if(ldep .and. ds > 0d0) then
       if(dtau > 50d0) then
          um = sc + (ui - sc)/dtau
       elseif(dtau > 1d-8) then
          um = sc - (ui - sc)*em1(-dtau)/dtau
       else
          um = ui + 0.5d0*(sc - ui)*dtau
       endif
       w = (log(nuc) - lj0)/dln + 1d0
       jl = int(w)
       if(jl >= 1 .and. jl < nj) then
          w = w - jl
          um = wdep*ds*um*nuc**3
          jnw(jl,c) = jnw(jl,c) + (1d0-w)*um
          jnw(jl+1,c) = jnw(jl+1,c) + w*um
       endif
    endif
    call cont(ui,dtau,sc)
  end subroutine contd
!
  subroutine ngacc(a)
!-- Ng's acceleration (first order) of the J of the passes from the last
!-- three iterates jold, jfi, jnw: jnw <- (1-a) jnw + a jfi with a from
!-- the least squares of the relative changes (skipped when it would
!-- extrapolate by more than a factor 100; J stays >= 0)
    real*8,intent(out) :: a
    real*8 :: s1,s2,d0,d1,w
    integer :: c,jl
    a = 0d0
    s1 = 0d0
    s2 = 0d0
    do c=1,ncell
       do jl=1,nj
          if(jnw(jl,c) <= 0d0) cycle
          w = 1d0/jnw(jl,c)**2
          d0 = jnw(jl,c) - jfi(jl,c)
          d1 = jfi(jl,c) - jold(jl,c)
          s1 = s1 + w*d0*(d0 - d1)
          s2 = s2 + w*(d0 - d1)**2
       enddo
    enddo
    if(s2 <= 0d0) return
    a = s1/s2
    if(abs(a) > 100d0) then
       a = 0d0
       return
    endif
    jnw = max((1d0-a)*jnw + a*jfi, 0d0)
  end subroutine ngacc
!
  subroutine jinit
!-- the Monte Carlo J of the groups on the J grid (linear in ln nu between
!-- the group centres)
    integer :: c,jl,g,g2
    real*8 :: l,w
    do jl=1,nj
       l = lj0 + (jl-1)*dln
       g = groupof(exp(l))
       if(g <= 0) then
          jfi(jl,:) = 0d0
          cycle
       endif
       if(l < lng(g)) then
          g2 = min(g+1,ng)
       else
          g2 = max(g-1,1)
       endif
       w = 0d0
       if(g2 /= g) w = (l - lng(g))/(lng(g2) - lng(g))
       do c=1,ncell
          jfi(jl,c) = (1d0-w)*jgx(g,c) + w*jgx(g2,c)
       enddo
    enddo
  end subroutine jinit
!
  subroutine jemis
!-- electron-scattering emissivity of the J of the grid: sigma times J
!-- averaged over the thermal-electron kernel (weights renormalized at the
!-- ends of the grid)
    integer :: c,jl,j2,m
    real*8 :: sel,s1,s2,w
    real*8,allocatable :: kw(:)
    do c=1,ncell
       sel = sqrt(2d0*kb/(me*tempinv(c)))/cc
       m = int(ukmax*sel/dln) + 1
       allocate(kw(0:m))
       do j2=0,m
          kw(j2) = esk(j2*dln/sel)
       enddo
       do jl=1,nj
          s1 = 0d0
          s2 = 0d0
          do j2=max(1,jl-m),min(nj,jl+m)
             w = kw(abs(jl-j2))
             s1 = s1 + w*jfi(j2,c)
             s2 = s2 + w
          enddo
          etj(jl,c) = 0d0
          if(s2 > 0d0) etj(jl,c) = grd_sig(c)*s1/s2
       enddo
       deallocate(kw)
    enddo
  end subroutine jemis
!
  real*8 function esem(nuc,g,c)
!-- electron-scattering emissivity at comoving frequency nuc: of the J of
!-- the passes on its grid, else of the Monte Carlo J (etesi)
    real*8,intent(in) :: nuc
    integer,intent(in) :: g,c
    real*8 :: w
    integer :: jl
    if(lj) then
       w = (log(nuc) - lj0)/dln + 1d0
       jl = int(w)
       if(jl >= 1 .and. jl < nj) then
          w = w - jl
          esem = (1d0-w)*etj(jl,c) + w*etj(jl+1,c)
          return
       endif
    endif
    esem = etesi(nuc,g,c)
  end function esem
!
  real*8 function srcline(k,c,nul)
!-- source function of global line k in cell c from the level populations
    use nltextrmod, only:nxt_x,nxt_nel,nxt_w
    use nlteatommod
    integer,intent(in) :: k,c
    real*8,intent(in) :: nul
    integer :: i,iel,o,l,u
    real*8 :: nl,nup,rr
    srcline = 0d0
    i = nxa_lion(k)
    iel = nxa_lel(k)
    o = nxa_el(iel)%xoff
    l = nxa_ion(i)%ll(nxa_lloc(k))
    u = nxa_ion(i)%lu(nxa_lloc(k))
    nl = nxt_x(o+nxa_lsl(k),c)*nxt_w(nxa_llv(k),c)
    nup = nxt_x(o+nxa_lsu(k),c)*nxt_w(nxa_luv(k),c)
    if(nup <= 0d0) return
    rr = nxa_ion(i)%g(u)*nl/(nxa_ion(i)%g(l)*nup) - 1d0
    if(rr <= 0d0) return
    srcline = 2d0*hh*nul**3/cc**2/rr
  end function srcline
!
  integer function lineidx(w,lfirst)
!-- first active line at or above wavelength w (lfirst), or last at or
!-- below it
    use nltextrmod, only:nxt_nact,nxt_wlact
    real*8,intent(in) :: w
    logical,intent(in) :: lfirst
    integer :: a,b,m
    a = 0
    b = nxt_nact + 1
    do while(b-a>1)
       m = (a+b)/2
       if(nxt_wlact(m) < w) then
          a = m
       else
          b = m
       endif
    enddo
!-- nxt_wlact(a) < w <= nxt_wlact(b)
    if(lfirst) then
       lineidx = b
    else
       lineidx = a
       if(b <= nxt_nact) then
          if(nxt_wlact(b) == w) lineidx = b
       endif
    endif
  end function lineidx
!
  integer function groupof(nuc)
!-- group of comoving frequency nuc (0: outside the grid)
    real*8,intent(in) :: nuc
    real*8 :: w
    integer :: a,b,m
    w = cc/nuc
    groupof = 0
    if(w < grp_wl(1) .or. w >= grp_wl(ng+1)) return
    a = 1
    b = ng + 1
    do while(b-a>1)
       m = (a+b)/2
       if(grp_wl(m) <= w) then
          a = m
       else
          b = m
       endif
    enddo
    groupof = a
  end function groupof
!
  integer function cellof(r)
    real*8,intent(in) :: r
    integer :: a,b,m
    a = 1
    b = ncell + 1
    do while(b-a>1)
       m = (a+b)/2
       if(xarr(m) <= r) then
          a = m
       else
          b = m
       endif
    enddo
    cellof = a
  end function cellof
!
  real*8 function vel(c,r)
!-- velocity at radius r in cell c (linear in r)
    integer,intent(in) :: c
    real*8,intent(in) :: r
    vel = vxarr(c) + (vxarr(c+1)-vxarr(c))*(r - xarr(c))/(xarr(c+1)-xarr(c))
  end function vel
end subroutine fi_spectrum
!
!
pure real*8 function em1(x)
!-- exp(x) - 1 without cancellation
  real*8,intent(in) :: x
  if(abs(x) < 1d-5) then
     em1 = x*(1d0 + 0.5d0*x*(1d0 + x/3d0))
  else
     em1 = exp(x) - 1d0
  endif
end function em1
!
end module formalmod
! vim: fdm=marker
