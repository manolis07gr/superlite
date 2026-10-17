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
!    J redistributed by the thermal electron motions (Gaussian of
!    sqrt(2kT/m_e)/c in frequency), the scattering part of the expansion
!    opacity of the LTE lines (J);
!  - NLTE lines as Sobolev resonances where the comoving frequency crosses
!    the line (several in a non-monotonic flow): depth K/|Q(mu)| along the
!    ray and source function from the level populations.
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
!
contains
!
!
subroutine fi_spectrum(ncell,xarr,vxarr,tempinv,wl1,wl2,dvkms,npc,impi,nmpi,lmpi0)
  use gridmod, only:grd_cap,grd_sig,grd_capemit
  use bbsobmod, only:grd_capsc
  use groupmod, only:grp_ng,grp_wl
  use nltextrmod, only:nxt_jg,nxt_nug,nxt_dnug,grd_emnbf
  use mpimod, only:allreduce_sum
  integer,intent(in) :: ncell,npc,impi,nmpi
  real*8,intent(in) :: xarr(ncell+1),vxarr(ncell+1),tempinv(ncell),wl1,wl2,dvkms
  logical,intent(in) :: lmpi0
  integer :: nnu,i,j,ic,ig,ng,np,k,jg
  real*8 :: dln,t,x,bnu,se,sig2,wsum,lnu0,dl
  real*8,allocatable :: wl(:),lnu(:),chi(:,:),eta(:,:),p(:),wp(:),lng(:)
!
  ng = grp_ng
!-- observer wavelengths
  dln = dvkms*1d5/cc
  nnu = int(log(wl2/wl1)/dln) + 1
  allocate(wl(nnu),lnu(nnu))
  do i=1,nnu
     wl(i) = wl1*exp((i-1)*dln)
  enddo
!-- continuum opacity and emissivity per group and cell (comoving frame)
  allocate(chi(ng,ncell),eta(ng,ncell),lng(ng))
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
           eta(ig,ic) = eta(ig,ic) + grd_capsc(ig,ic)*nxt_jg(ig,ic)
        endif
!-- J averaged over the redistribution kernel
        sig2 = 0d0
        wsum = 0d0
        do jg=max(1,ig-60),min(ng,ig+60)
           dl = (lng(jg) - lng(ig))/se
           if(abs(dl) > 4d0) cycle
           x = exp(-0.5d0*dl*dl)*nxt_dnug(jg)/nxt_nug(jg)
           sig2 = sig2 + x*nxt_jg(jg,ic)
           wsum = wsum + x
        enddo
        if(wsum > 0d0) eta(ig,ic) = eta(ig,ic) + grd_sig(ic)*sig2/wsum
     enddo
  enddo
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
  do i=impi+1,nnu,nmpi
     lnu0 = 0d0
     do j=1,np
        lnu0 = lnu0 + wp(j)*ray(p(j),cc/wl(i))
     enddo
     lnu(i) = 8d0*pi*pi*lnu0
  enddo
  call allreduce_sum(nnu,lnu)
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
  deallocate(wl,lnu,chi,eta,lng,p,wp)
!
contains
!
  real*8 function ray(pp,nu)
!-- observer-frame intensity I_nu of the ray with impact parameter pp
    real*8,intent(in) :: pp,nu
    integer :: ie,ne,iz,c,ns,is
    real*8 :: zs(2*ncell+2),za,zb,zm,ra,rb,rm,dvdr,va,vb,nua,nub,vlos(2),ui
    ui = 0d0
!-- crossings with the shell edges, from the back to the front
    ne = 0
    if(pp < xarr(1)) then
!-- core ray: from the inner boundary on the observer's side
       ne = 1
       zs(1) = sqrt(xarr(1)**2 - pp*pp)
       c = 1
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
!-- the intervals: each inside one cell
    do iz=1,ne-1
       za = zs(iz)
       zb = zs(iz+1)
       if(zb <= za) cycle
       zm = 0.5d0*(za+zb)
       rm = sqrt(pp*pp + zm*zm)
       c = cellof(rm)
       dvdr = (vxarr(c+1)-vxarr(c))/(xarr(c+1)-xarr(c))
!-- substeps by the change of the line-of-sight velocity (50 km/s)
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
    ray = ui*nu**3
  end function ray
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
    if(g>0) core_i = nxt_jg(g,1)/nuc**3
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
       if(chi(g,c) > 0d0) sc = eta(g,c)/chi(g,c)/num**3
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
          call cont(ui,ch*ds,sc)
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
    call cont(ui,ch*ds,sc)
  end subroutine substep
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
