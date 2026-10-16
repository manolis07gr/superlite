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
!  * NLTE elements (nlteatommod); the populations of H and He (Z <= te_zsolve)
!    are re-solved at each trial T, those of heavier elements held (their
!    rate coefficients follow T):
!      photoionization heating of every level l (edge nu_l)
!        4 pi n_el sum_l sum_g (x_l - x_l* e^{-h(nu-nu_l)/kT}) J_g Int_g sigma_l (1 - nu_l/nu) dnu
!      recombination cooling  n_e sum_l n_c(l) rc_l(T)  (kinetic energy, Milne)
!      collisions  n_e sum_{l<u} (n_l q_lu - n_u q_ul) (E_u - E_l)
!                + n_e sum_l (n_l c_l - n_e n_c c3b_l) chi_l
!  * Compton  sigma_T n_e/(m_e c^2) 4 pi sum_g J_g dnu_g (4kT - h nu_g (1 + n_g)),
!    n_g = c^2 J_g/(2 h nu_g^3) the photon occupation number (induced
!    scattering; the balance is exact for a Planck field at T)
!  * heating not carried by the radiation field: energy injected in the
!    cell outside shocks ('lumg' source)
!  * photons above te_ntemin make fast photoelectrons: only the heat share
!    fh (Shull & van Steenberg 1985) of their absorption heats the gas, the
!    rest ionizes and excites H and He (rates nti, ntx in the populations)
!
! f is bracketed in log T and the root found by bisection/secant steps
! (Illinois).  The new temperature is mixed geometrically with the old one
! (te_damp) and its change per iteration limited to a factor te_maxfac;
! te_step damps per cell adaptively: the weight is halved when the proposed
! change reverses sign (the 2-cycles of a temperature update against an
! optically thick radiation field that responds to it) and relaxes back
! to te_damp while the changes keep their sign.
!-----------------------------------------------------------------------
  use nlteatommod
  use nltexmod, only:nx_solve_cell
  use nltextrmod, only:nxt_shares,nxt_elh,nxt_elhe
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
  integer,public :: te_zsolve = 2
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
  if(te_ntemin <= 0d0 .or. nxt_elh==0) return
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
        w = jg(ig)*(nxa_sgm(ig,lv)*dnug(ig) - hh*nxa_lvnu0(lv)*nxa_sgh(ig,lv))
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
subroutine te_terms(temp,ng,nug,dnug,jg,capabs,capem,nel,x,ne,gam_ext,terms,pa,pb)
!-----------------------------------------------------------------------
! Heating and cooling terms at temperature temp for fixed populations.
! terms(1) LTE opacities, net cooling         terms(2) NLTE photoionization heating
! terms(3) NLTE recombination cooling         terms(4) NLTE collisional net cooling
! terms(5) Compton net cooling                terms(6) external heating
! f = t1 - t2 + t3 + t4 + t5 - t6
!-----------------------------------------------------------------------
  integer,intent(in) :: ng
  real*8,intent(in) :: temp,nug(ng),dnug(ng),jg(ng),nel(nxa_nel),x(nxa_nst),ne,gam_ext
  real*4,intent(in) :: capabs(ng),capem(ng)
  real*8,intent(out) :: terms(te_nterm)
  real*8,intent(in),optional :: pa(nxa_nlv),pb(nxa_nlv)  !te_pisums (computed if absent)
  integer :: ig,lv,iel,o,i,j,k,n,l,u
  real*8 :: hkt,bg,s,xs,fh,jh,nu,sh,eu,sst
  real*8 :: qa(nxa_nlv),qb(nxa_nlv)
  real*8,allocatable :: qup(:,:),qdn(:,:),cion(:),c3b(:)
!
  hkt = hh/(kb*temp)
  terms = 0d0
  fh = te_heatshare(nel,x)
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
  s = 0d0
  sh = 0d0
  do lv=1,nxa_nlv
     iel = nxa_lvel(lv)
     if(nel(iel) <= 0d0) cycle
     o = nxa_el(iel)%xoff
     xs = x(o+nxa_lvc(lv))*ne*nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)
!-- Int_g sigma (1 - nu_l/nu) dnu = sgm dnu - h nu_l sgh; stimulated
!-- recombination from the edge up to h(nu-nu0) = 40 kT
     sst = 0d0
     do ig=nxa_g1(lv),1,-1
        nu = max(nug(ig),nxa_lvnu0(lv))
        if(hkt*(nu-nxa_lvnu0(lv)) > 40d0) exit
        jh = jg(ig)
        if(te_ntemin > 0d0 .and. hh*nug(ig) > te_ntemin) jh = fh*jg(ig)
        sst = sst + jh*exp(-hkt*(nu-nxa_lvnu0(lv)))* &
             (nxa_sgm(ig,lv)*dnug(ig) - hh*nxa_lvnu0(lv)*nxa_sgh(ig,lv))
     enddo
     s = s + nel(iel)*(x(o+nxa_lvs(lv))*(qa(lv) + fh*qb(lv)) - xs*sst)
     sh = sh + nel(iel)*x(o+nxa_lvc(lv))*nxa_rcool(lv,temp)
  enddo
  terms(2) = 4d0*pi*s
  terms(3) = ne*sh
!-- collisions: excitation/de-excitation and ionization/three-body recombination
  s = 0d0
  do iel=1,nxa_nel
     if(nel(iel) <= 0d0) cycle
     o = nxa_el(iel)%xoff
     do j=1,nxa_el(iel)%nion
        i = nxa_el(iel)%ion(j)
        n = nxa_ion(i)%nlev
        allocate(qup(n,n),qdn(n,n),cion(n),c3b(n))
        call nxa_coll(nxa_ion(i),temp,qup,qdn,cion,c3b)
        do l=1,n-1
           do u=l+1,n
              eu = nxa_ion(i)%e(u) - nxa_ion(i)%e(l)
              s = s + nel(iel)*(x(o+nxa_ion(i)%soff+l)*qup(l,u) - x(o+nxa_ion(i)%soff+u)*qdn(u,l))*eu
           enddo
        enddo
        do k=1,n
           lv = nxa_ion(i)%lvoff + k
           s = s + nel(iel)*(x(o+nxa_ion(i)%soff+k)*cion(k) - x(o+nxa_lvc(lv))*ne*c3b(k))* &
                hh*nxa_ion(i)%nu0(k)
        enddo
        deallocate(qup,qdn,cion,c3b)
     enddo
  enddo
  terms(4) = ne*s
  terms(6) = gam_ext
end subroutine te_terms
!
!
subroutine te_balance(told,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,capabs, &
     capem,gam_ext,x,tnew,ierr,tsolout)
!-----------------------------------------------------------------------
! Solve f(T) = 0 with the NLTE populations re-solved at each trial T; x is
! the initial guess for the populations (unchanged on return).  tnew is
! the damped and limited new temperature, tsolout (optional) the root;
! ierr = 1 if f has no sign change in [te_tmin, te_tmax] (the root is the
! bound).
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
  integer :: it,side
  real*8 :: a,b,c,fa,fb,fc,tsol,lo,hi
  real*8 :: pa(nxa_nlv),pb(nxa_nlv)
!
  ierr = 0
  call te_pisums(ng,nug,dnug,jg,pa,pb)
!-- bracket in log T around the old temperature
  a = log(told)
  fa = fval(a)
  b = a
  fb = fa
  lo = log(te_tmin)
  hi = log(te_tmax)
  do it=1,40
     if(fa < 0d0) then
        b = min(a + 0.3d0, hi)
     else
        b = max(a - 0.3d0, lo)
     endif
     fb = fval(b)
     if(fa*fb <= 0d0) exit
     if(b <= lo .or. b >= hi) exit
     a = b
     fa = fb
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
        if(abs(b - a) < 1d-4 .or. fc == 0d0) exit
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
contains
  real*8 function fval(lt)
    real*8,intent(in) :: lt
    real*8 :: t,xx(nxa_nst),ne,beta(nxa_nlin),terms(te_nterm)
    integer :: ie
    t = exp(lt)
    xx = x
    call nx_solve_cell(t,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,xx,ne,beta,ie, &
         lsolve=nxa_el(:)%z <= te_zsolve)
    if(ie == 1) then
       xx = x
       ne = neoth + free_el(xx)
    endif
    call te_terms(t,ng,nug,dnug,jg,capabs,capem,nel,xx,ne,gam_ext,terms,pa,pb)
    fval = terms(1) - terms(2) + terms(3) + terms(4) + terms(5) - terms(6)
  end function fval
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
