!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module tebalmod
!-----------------------------------------------------------------------
! Electron temperature from the thermal balance of the gas, with the
! radiation field of the previous Monte Carlo iteration:
!
!   f(T) = Lambda(T) - Gamma(T) = 0,
!
! per unit volume [erg/s/cm^3]:
!  * opacities of the LTE elements (free-free, bound-free, thermalized
!    share of the lines) from the last opacity calculation, absorption
!    and emission (bound-free emission from ground states in Saha
!    equilibrium with the next ion at T_e):
!      4 pi sum_g (kem_g B_g(T) - kabs_g J_g) dnu_g
!  * NLTE elements (nlteatommod), re-solved at each trial T (all elements by
!    default; with te_zsolve < 99 only Z <= te_zsolve, the others held: a
!    held population gives the wrong temperature dependence of the
!    collisional term of radiatively pumped levels):
!      photoionization heating of every level l (edge nu_l)
!        4 pi n_el sum_l sum_g (x_l - x_l* e^{-h(nu-nu_l)/kT}) J_g Int_g sigma_l (1 - nu_l/nu) dnu
!      recombination cooling  n_e sum_l n_c(l) rc_l(T)  (kinetic energy, Milne)
!      collisions  n_e sum_{l<u} (n_l q_lu - n_u q_ul) (E_u - E_l)
!                + n_e sum_l (n_l c_l - n_e n_c c3b_l) chi_l
!      (with super-levels: between states, nxa_coll_el; level populations
!      n_l = n_S w_l), and the lines within super-levels, thermal at T_e
!      (net emission beta (n_u A_ul - (n_l B_lu - n_u B_ul) J) h nu)
!  * Compton  sigma_T n_e/(m_e c^2) 4 pi sum_g J_g dnu_g (4kT - h nu_g (1 + n_g)),
!    n_g = c^2 J_g/(2 h nu_g^3) the photon occupation number (induced
!    scattering; the balance is exact for a Planck field at T)
!  * heating not carried by the radiation field: energy injected in the
!    cell outside shocks ('lumg' source)
!  * photons above te_ntemin make fast photoelectrons: only the heat share
!    fh (Shull & van Steenberg 1985, or the Spencer-Fano solution of the
!    cell, nltentmod) of their absorption heats the gas, the rest ionizes
!    and excites H and He, or all NLTE species (rates in the populations)
!
! f is bracketed in log T (steps of 5, 10, 20, 30 % ... from the old
! temperature) and the root found by secant steps (Illinois) to 0.1 %;
! each trial starts from the populations of the nearest temperature tried,
! solved to 1e-3; after the first trial only the elements with at least
! 1 % of the heating and cooling terms are re-solved.  The new temperature is mixed geometrically with the old one
! (te_damp) and its change per iteration limited to a factor te_maxfac;
! te_step damps per cell adaptively: the weight is halved when the proposed
! change reverses sign (the 2-cycles of a temperature update against an
! optically thick radiation field that responds to it) and relaxes back
! to te_damp while the changes keep their sign.
!-----------------------------------------------------------------------
  use nlteatommod
  use nltedismod
  use nltexmod, only:nx_solve_cell,nx_cx_heat,nx_jline,nx_ovl_cell,nx_ovl_keep,nxo_npair,nx_dis
  use nltextrmod, only:nxt_shares,nxt_elh,nxt_elhe
  use nltentmod, only:nxn_on,nxn_fhc
  implicit none
  private
  real*8,parameter :: pi = 3.141592653589793d0
  real*8,parameter :: hh = 6.6260693d-27, cc = 2.99792458d10, kb = 1.3806505d-16
  real*8,parameter :: sigt = 6.6524587d-25, mec2 = 9.1093826d-28*2.99792458d10**2
  real*8,public :: te_damp = 0.5d0              !weight of the new temperature
  real*8,public :: te_maxfac = 1.5d0            !maximum change per iteration
  real*8,public :: te_tmin = 1d3, te_tmax = 1d7
  real*8,public :: te_ntemin = 0d0              !fast photoelectrons above [erg], 0: off
!-- elements with Z <= te_zsolve are re-solved at each trial temperature;
!-- the populations of the heavier ones are held (their rates still follow T)
  integer,public :: te_zsolve = 99
  integer,parameter,public :: te_nterm = 6
  public :: te_balance,te_terms,te_step,te_heatshare
!
contains
!
!
pure function te_heatshare(nel,x) result(fh)
!-- heat share of the energy of fast photoelectrons (1 without NLTE H)
  real*8,intent(in) :: nel(nxa_nel),x(nxa_nst)
  real*8 :: fh,fih,feh,fihe,fehe,y
  fh = 1d0
  if(te_ntemin <= 0d0) return
!-- (Spencer-Fano with all NLTE species: the heat share of the cell)
  if(nxn_on) then
     fh = nxn_fhc
     return
  endif
  if(nxt_elh==0) return
  y = 0d0
  if(nxt_elhe>0 .and. nel(nxt_elh)>0d0) y = 12d0*nel(nxt_elhe)/nel(nxt_elh)
  call nxt_shares(x(nxa_el(nxt_elh)%xoff+nxa_el(nxt_elh)%ns),y,fh,fih,feh,fihe,fehe)
end function te_heatshare
!
!
subroutine te_pisums(ng,nug,dnug,jg,pa,pb)
!-- photoionization heating integrals of each level that do not depend on
!-- T: Sum_g J_g Int_g sigma (1 - nu_l/nu) dnu for groups below (pa) and
!-- above (pb) the fast-electron threshold te_ntemin
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng)
  real*8,intent(out) :: pa(nxa_nlv),pb(nxa_nlv)
  integer :: lv,ig
  real*8 :: w
  pa = 0d0
  pb = 0d0
  do lv=1,nxa_nlv
     do ig=1,nxa_g1(lv)
        w = jg(ig)*(nxa_sgm(nxa_sgo(lv)+ig)*dnug(ig) - hh*nxa_lvnu0(lv)*nxa_sgh(nxa_sgo(lv)+ig))
        if(te_ntemin > 0d0 .and. hh*nug(ig) > te_ntemin) then
           pb(lv) = pb(lv) + w
        else
           pa(lv) = pa(lv) + w
        endif
     enddo
  enddo
end subroutine te_pisums
!
!
subroutine te_terms(temp,ng,nug,dnug,jg,capabs,capem,nel,x,ne,gam_ext,terms,pa,pb,elc,beta,jl)
!-----------------------------------------------------------------------
! Heating and cooling terms at temperature temp for fixed populations.
! terms(1) LTE opacities, net cooling         terms(2) NLTE photoionization heating
! terms(3) NLTE recombination cooling         terms(4) NLTE collisional net cooling
!                                              (electrons, charge exchange, and
!                                              lines within super-levels)
! terms(5) Compton net cooling                terms(6) external heating
! f = t1 - t2 + t3 + t4 + t5 - t6
!-----------------------------------------------------------------------
  integer,intent(in) :: ng
  real*8,intent(in) :: temp,nug(ng),dnug(ng),jg(ng),nel(nxa_nel),x(nxa_nst),ne,gam_ext
  real*4,intent(in) :: capabs(ng),capem(ng)
  real*8,intent(out) :: terms(te_nterm)
  real*8,intent(in),optional :: pa(nxa_nlv),pb(nxa_nlv)  !te_pisums (computed if absent)
!-- (optional) per element: sum of the magnitudes of its photoionization
!-- heating, recombination cooling and collisional terms [erg/s/cm^3]
  real*8,intent(out),optional :: elc(nxa_nel)
!-- (optional) escape probabilities and incident intensities of the lines:
!-- the lines within super-levels (none without super-levels)
  real*8,intent(in),optional :: beta(nxa_nlin),jl(nxa_nlin)
  real*8 :: se(nxa_nel),sh0
  integer :: ig,lv,iel,o,i,j,k,ns,lo,s,u
  real*8 :: hkt,bg,sm,xs,fh,jh,nu,sh,sst,chi
  real*8 :: qa(nxa_nlv),qb(nxa_nlv),wg(nxa_nlv),wc
  real*8 :: pd(3),wo,nuc,dm,dh,dhs,dsp,ph,ps,pr
  logical :: lgo
  real*8,allocatable :: cq(:,:),eup(:,:),edn(:,:),cil(:),c3l(:)
!
  hkt = hh/(kb*temp)
  terms = 0d0
  fh = te_heatshare(nel,x)
!-- (level dissolution: the cell's occupation probabilities, nltedismod)
  pd = 0d0
  if(nxd_on) call nxd_cellx(ne,nel,x,temp,pd)
  call nxd_weights(temp,pd,wg)
!-- LTE opacities and Compton
  do ig=1,ng
     if(hkt*nug(ig) < 700d0) then
        bg = 2d0*hh*nug(ig)**3/cc**2/(exp(hkt*nug(ig)) - 1d0)
     else
        bg = 0d0
     endif
     jh = jg(ig)
     if(te_ntemin > 0d0 .and. hh*nug(ig) > te_ntemin) jh = fh*jg(ig)
     terms(1) = terms(1) + (capem(ig)*bg - capabs(ig)*jh)*dnug(ig)
     terms(5) = terms(5) + jg(ig)*dnug(ig)*(4d0*kb*temp - &
          hh*nug(ig)*(1d0 + cc**2*jg(ig)/(2d0*hh*nug(ig)**3)))
  enddo
  terms(1) = 4d0*pi*terms(1)
  terms(5) = 4d0*pi*sigt*ne/mec2*terms(5)
!-- NLTE photoionization heating (net of stimulated recombination) and
!-- recombination cooling
  if(present(pa) .and. present(pb)) then
     qa = pa
     qb = pb
  else
     call te_pisums(ng,nug,dnug,jg,qa,qb)
  endif
  sm = 0d0
  sh = 0d0
  se = 0d0
  do lv=1,nxa_nlv
     iel = nxa_lvel(lv)
     if(nel(iel) <= 0d0) cycle
     o = nxa_el(iel)%xoff
!-- the continuum: level 1 of the next ion
     wc = 1d0
     if(nxa_lvcl(lv)>0) wc = wg(nxa_lvcl(lv))
     wo = nxd_lvw(lv,pd)
     xs = x(o+nxa_lvc(lv))*wc*ne*nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)*wo
!-- Int_g sigma (1 - nu_l/nu) dnu = sgm dnu - h nu_l sgh; stimulated
!-- recombination from the edge up to h(nu-nu0) = 40 kT
     sst = 0d0
     do ig=nxa_g1(lv),1,-1
        nu = max(nug(ig),nxa_lvnu0(lv))
        if(hkt*(nu-nxa_lvnu0(lv)) > 40d0) exit
        jh = jg(ig)
        if(te_ntemin > 0d0 .and. hh*nug(ig) > te_ntemin) jh = fh*jg(ig)
        sst = sst + jh*exp(-hkt*(nu-nxa_lvnu0(lv)))* &
             (nxa_sgm(nxa_sgo(lv)+ig)*dnug(ig) - hh*nxa_lvnu0(lv)*nxa_sgh(nxa_sgo(lv)+ig))
     enddo
!-- (level dissolution: the dissolved continuum below the edge takes
!-- h(nu0 - nu) from the gas per photoionization and gives it per
!-- recombination, as CMFGEN's QUAD_MULTI)
     ph = 0d0
     ps = 0d0
     pr = 0d0
     lgo = pd(1) > 0d0
     if(lgo) lgo = nxd_ga(lv) > 0
     if(lgo) then
        nuc = nxd_numin(nxa_lvnu0(lv),nxd_z(lv),pd)
        do ig=nxd_ga(lv),min(nxd_gb(lv),ng)
           if(nug(ig) + 0.5d0*dnug(ig) <= nuc) exit
           call nxd_gint(lv,ig,pd,temp,nuc,dm,dh,dhs,dsp)
           if(dm <= 0d0) cycle
           ph = ph + jg(ig)*(dm*dnug(ig) - hh*nxa_lvnu0(lv)*dh)
           ps = ps + jg(ig)*hh*(nug(ig) - nxa_lvnu0(lv))*dhs
           pr = pr + hh*(nug(ig) - nxa_lvnu0(lv))*dsp
        enddo
     endif
     sh0 = nel(iel)*(x(o+nxa_lvs(lv))*wg(lv)*(qa(lv) + fh*qb(lv) + ph) - xs*(sst + ps))
     sm = sm + sh0
     se(iel) = se(iel) + 4d0*pi*abs(sh0)
     sh0 = nel(iel)*x(o+nxa_lvc(lv))*wc*(nxa_rcool(lv,temp)*wo + &
          4d0*pi*nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)*2d0/cc**2*pr*wo)
     sh = sh + sh0
     se(iel) = se(iel) + ne*abs(sh0)
  enddo
  terms(2) = 4d0*pi*sm
  terms(3) = ne*sh
!-- collisions: excitation/de-excitation between states and
!-- ionization/three-body recombination
  sm = 0d0
  do iel=1,nxa_nel
     if(nel(iel) <= 0d0) cycle
     o = nxa_el(iel)%xoff
     ns = nxa_el(iel)%ns
     lo = nxa_el(iel)%lvoff
     allocate(cq(ns,ns),eup(ns,ns),edn(ns,ns),cil(nxa_el(iel)%nlv),c3l(nxa_el(iel)%nlv))
     call nxa_coll_el(iel,temp,cq,eup,edn,cil,c3l)
     sh0 = sm
     do s=1,ns
        if(x(o+s) <= 0d0) cycle
        do u=1,ns
           if(u/=s) sm = sm + nel(iel)*x(o+s)*(eup(s,u) - edn(s,u))
        enddo
     enddo
     do k=1,nxa_el(iel)%nlv
        lv = lo + k
        i = nxa_lvion(lv)
        chi = hh*nxa_lvnu0(lv)
        wc = 1d0
        if(nxa_lvcl(lv)>0) wc = wg(nxa_lvcl(lv))
        sm = sm + nel(iel)*(x(o+nxa_lvs(lv))*wg(lv)*cil(k) - x(o+nxa_lvc(lv))*wc*ne*c3l(k))*chi
     enddo
     se(iel) = se(iel) + ne*abs(sm - sh0)
     deallocate(cq,eup,edn,cil,c3l)
  enddo
  terms(4) = ne*sm
  if(nxa_ncx>0) terms(4) = terms(4) - nx_cx_heat(temp,nel,x,wg)
!-- lines within super-levels: thermal at T_e
  if(nxa_anysl .and. present(beta) .and. present(jl)) then
     sm = 0d0
     do k=1,nxa_nlin
        if(nxa_lsl(k)/=nxa_lsu(k)) cycle
        iel = nxa_lel(k)
        if(nel(iel) <= 0d0) cycle
        o = nxa_el(iel)%xoff
        sh0 = nel(iel)*x(o+nxa_lsu(k))*beta(k)*hh*nxa_lnu(k)* &
             (wg(nxa_luv(k))*(nxa_la(k) + nxa_lbul(k)*jl(k)) - wg(nxa_llv(k))*nxa_lblu(k)*jl(k))
        sm = sm + sh0
        se(iel) = se(iel) + abs(sh0)
     enddo
     terms(4) = terms(4) + sm
  endif
  terms(6) = gam_ext
  if(present(elc)) elc = se
end subroutine te_terms
!
!
subroutine te_balance(told,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,capabs, &
     capem,gam_ext,x,tnew,ierr,tsolout,xprev,dtinv,scanlu)
!-----------------------------------------------------------------------
! Solve f(T) = 0 with the NLTE populations re-solved at each trial T; x is
! the initial guess for the populations (unchanged on return).  tnew is
! the damped and limited new temperature, tsolout (optional) the root;
! ierr = 1 if f has no sign change in [te_tmin, te_tmax] (the root is the
! bound).  xprev, dtinv (optional): time-dependent statistical equilibrium
! (nx_solve_cell).  scanlu (optional, diagnostics): after the solve, write
! f(T)/(sum of the magnitudes of the terms) at 41 temperatures within
! told*exp(+-0.3) to that unit.
!-----------------------------------------------------------------------
  real*8,intent(in) :: told,trad,nel(nxa_nel),neoth,qb,qv,gam_ext
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng),jb(nxa_nlin),nti(nxa_nst),ntx(nxa_nst)
  logical,intent(in) :: ljb(nxa_nlin)
  real*4,intent(in) :: capabs(ng),capem(ng)
  real*8,intent(in) :: x(nxa_nst)
  real*8,intent(out) :: tnew
  integer,intent(out) :: ierr
  real*8,intent(out),optional :: tsolout
  real*8,intent(in),optional :: xprev(nxa_nst),dtinv
  integer,intent(in),optional :: scanlu
  integer :: it,side,ielast
  real*8 :: nelast,xhlast
  real*8 :: a,b,c,fa,fb,fc,tsol,lo,hi,step,elc(nxa_nel),scale
  logical :: lsol(nxa_nel),lfirst
  real*8 :: pa(nxa_nlv),pb(nxa_nlv),jl(nxa_nlin),wg(nxa_nlv)
!
  ierr = 0
  call te_pisums(ng,nug,dnug,jg,pa,pb)
  do it=1,nxa_nlin
     jl(it) = nx_jline(it,ng,nug,jg,jb(it),ljb(it))
  enddo
!-- (overlapping lines: the zones of the cell at the old temperature, kept
!-- by the solves at the trial temperatures)
  if(nxo_npair > 0) then
     call nxd_weights(told,nx_dis,wg)
     call nx_ovl_cell(told,nel,x,wg,qb,qv)
     nx_ovl_keep = .true.
  endif
!-- first evaluation with every element (Z <= te_zsolve) re-solved; then
!-- only the elements with at least 1 % of the summed magnitudes of the
!-- NLTE heating and cooling terms (and of the other terms) are re-solved
  lsol = nxa_el(:)%z <= te_zsolve
  lfirst = .true.
!-- bracket in log T around the old temperature
  a = log(told)
  fa = fval(a)
  b = a
  fb = fa
  lo = log(te_tmin)
  hi = log(te_tmax)
  step = 0.05d0
  do it=1,40
     if(fa < 0d0) then
        b = min(a + step, hi)
     else
        b = max(a - step, lo)
     endif
     fb = fval(b)
     if(fa*fb <= 0d0) exit
     if(b <= lo .or. b >= hi) exit
     a = b
     fa = fb
     step = min(2d0*step, 0.3d0)
  enddo
  c = b
  if(fa*fb > 0d0) then
!-- no root: the bound
     ierr = 1
     tsol = exp(b)
  else
!-- Illinois (regula falsi) iterations
     side = 0
     do it=1,60
        c = (a*fb - b*fa)/(fb - fa)
        fc = fval(c)
        if(abs(b - a) < 1d-3 .or. fc == 0d0) exit
        if(fc*fb < 0d0) then
           a = b
           fa = fb
           b = c
           fb = fc
           side = 0
        else
           b = c
           fb = fc
           if(side == -1) fa = 0.5d0*fa
           side = -1
        endif
     enddo
     tsol = exp(c)
  endif
  if(present(tsolout)) tsolout = tsol
!-- damping and limiter
  tnew = told**(1d0-te_damp)*tsol**te_damp
  tnew = min(max(tnew, told/te_maxfac), told*te_maxfac)
  tnew = min(max(tnew, te_tmin), te_tmax)
!-- (diagnostics: the balance function around the old temperature)
  if(present(scanlu)) then
     if(scanlu > 0) then
        do it=0,40
           c = log(told) - 0.3d0 + 0.015d0*it
           fc = fval(c)
           write(scanlu,'(f10.1,es13.4,i3,2es13.4)') exp(c), fc/max(scale,1d-300), ielast, nelast, xhlast
        enddo
     endif
  endif
  nx_ovl_keep = .false.
contains
  real*8 function fval(lt)
    real*8,intent(in) :: lt
    real*8 :: t,xx(nxa_nst),ne,beta(nxa_nlin),terms(te_nterm)
    integer :: ie
    t = exp(lt)
!-- start from the populations at the old temperature: with the Lyman
!-- continuum in local detailed balance (J = S of the cell) the equations
!-- of a cell can have a neutral and an ionized solution besides the one
!-- that continues the old state, and a warm start from another trial
!-- temperature (formerly the nearest one) jumped between them
    xx = x
    call solve(t,xx,ne,beta,ie,100)
!-- (a solve that fails, singular or not converged, is retried with more
!-- iterations; if that fails too the terms are those of the old
!-- populations)
    if(ie /= 0) then
       xx = x
       call solve(t,xx,ne,beta,ie,400)
    endif
    ielast = ie
    nelast = ne
    xhlast = xx(nxa_el(1)%xoff + nxa_el(1)%ns)
    if(ie /= 0) then
       xx = x
       ne = neoth + free_el(xx)
    endif
    call te_terms(t,ng,nug,dnug,jg,capabs,capem,nel,xx,ne,gam_ext,terms,pa,pb,elc,beta,jl)
    fval = terms(1) - terms(2) + terms(3) + terms(4) + terms(5) - terms(6)
    if(lfirst) then
       lfirst = .false.
       scale = sum(elc) + abs(terms(1)) + abs(terms(5)) + abs(terms(6))
       lsol = lsol .and. elc >= 1d-2*scale
    endif
  end function fval
  subroutine solve(t,xx,ne,beta,ie,mit)
!-- statistical equilibrium at the trial temperature t from the guess xx
    real*8,intent(in) :: t
    real*8,intent(inout) :: xx(nxa_nst)
    real*8,intent(out) :: ne,beta(nxa_nlin)
    integer,intent(out) :: ie
    integer,intent(in) :: mit
    if(present(xprev) .and. present(dtinv)) then
       call nx_solve_cell(t,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,xx,ne,beta,ie, &
            lsolve=lsol,tol=1d-3,maxit=mit,xprev=xprev,dtinv=dtinv)
    else
       call nx_solve_cell(t,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,xx,ne,beta,ie, &
            lsolve=lsol,tol=1d-3,maxit=mit)
    endif
  end subroutine solve
  real*8 function free_el(xx)
    real*8,intent(in) :: xx(nxa_nst)
    integer :: je,jo
    free_el = 0d0
    do je=1,nxa_nel
       jo = nxa_el(je)%xoff
       free_el = free_el + nel(je)*sum(nxa_el(je)%chg*xx(jo+1:jo+nxa_el(je)%ns))
    enddo
  end function free_el
end subroutine te_balance
!
!
pure subroutine te_step(told,tsol,dprev,dloc,tnew)
!-- adaptive damping of one cell: dprev, the previous proposed change of
!-- ln T, and dloc, the cell's weight (start: 0 and te_damp), are updated
  real*8,intent(in) :: told,tsol
  real*8,intent(inout) :: dprev,dloc
  real*8,intent(out) :: tnew
  real*8 :: dl
  dl = log(tsol/told)
  if(dl*dprev < 0d0) then
     dloc = max(0.5d0*dloc, 0.05d0)
  else
     dloc = min(1.25d0*dloc, te_damp)
  endif
  dprev = dl
  tnew = told*exp(dloc*dl)
  tnew = min(max(tnew, told/te_maxfac), told*te_maxfac)
  tnew = min(max(tnew, te_tmin), te_tmax)
end subroutine te_step
!
end module tebalmod
! vim: fdm=marker
