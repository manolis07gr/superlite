!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nltexmod
!-----------------------------------------------------------------------
! Statistical equilibrium of the NLTE elements (model atoms: nlteatommod)
! in one cell.  All ions of an element are solved together: the levels of
! each NLTE ion and the top ion, with
!  - bound-bound radiative rates in the Sobolev approximation,
!       R_lu = beta B_lu Jb,   R_ul = beta (A_ul + B_ul Jb),
!    Jb the intensity incident on the resonance (Monte Carlo estimator) and
!    beta the escape probability for the direction-dependent velocity
!    gradient Q(mu) = mu^2 dv/dr + (1-mu^2) v/r (nlteh_beta)
!  - electron-impact excitation, de-excitation, ionization and three-body
!    recombination (nxa_coll)
!  - photoionization (to the ground state of the next ion) and stimulated
!    recombination by the Monte Carlo J of the groups, radiative
!    recombination from the Milne relation (nxa_alpha)
!  - two-photon decays
!  - non-thermal ionization and excitation of ground states (rates given)
!  - charge exchange with the other NLTE elements (nxa_cx; the partner's
!    populations of the current iteration)
! The electron density is iterated with the ionization of all elements;
! electrons of the other elements (neoth) come from the LTE/nebular EOS.
! Populations are fractions x of the element's number density, per state:
! a level, or a super-level whose levels share its departure coefficient
! (Boltzmann at T_e within it, shares nxa_wel); the rates between states
! are the sums over their levels with these shares, transitions within a
! super-level drop out.
!-----------------------------------------------------------------------
  use nlteatommod
  use nltehmod, only:nlteh_beta,bplanck
  implicit none
  private
  real*8,parameter :: pi = 3.141592653589793d0
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
!
  public :: nx_el_rates,nx_el_matrix,nx_el_step,nx_solve_cell,nx_lte,nx_kq,nx_jline
  public :: nx_beta_table,nx_beta,nx_cx_rates,nx_cx_heat
!-- rates of one element at a temperature (nx_el_rates); s: element-local
!-- state, k: element-local level
  type,public :: nxrt
     real*8,allocatable :: w(:)       !(nlv) share of level k in its state
     real*8,allocatable :: wc(:)      !(nlv) share of level 1 of the next ion (the continuum of k) in its state
     real*8,allocatable :: eps(:)     !(ns) state energies [erg] (mean of a super-level)
     real*8,allocatable :: cq(:,:)    !(ns,ns) electron-impact s -> s' per particle in s [cm^3/s]
     real*8,allocatable :: eup(:,:),edn(:,:) !(ns,ns) their energy taken from / given to the electrons [erg cm^3/s]
     real*8,allocatable :: cil(:),c3l(:) !(nlv) collisional ionization of k [cm^3/s]; three-body into k per continuum particle [cm^6/s]
     real*8,allocatable :: gaml(:)    !(nlv) photoionization of k [1/s]
     real*8,allocatable :: gstl(:),arl(:) !(nlv) stimulated, spontaneous recombination into k per continuum particle [cm^3/s]
!-- the same summed over the levels of each state (with the shares w; the
!-- recombination coefficients per particle in the continuum state)
     real*8,allocatable :: ci(:),c3(:),gam(:),gst(:),ar(:) !(ns)
  end type nxrt
!-- escape probability table of one cell: beta against ln(K/qs) on
!-- [nx_blt0, nx_blt0 + (nx_nbt-1) nx_bdlt], qs = max(|dv/dr|,|v/r|)
  integer,parameter,public :: nx_nbt = 401
  real*8,parameter :: nx_blt0 = -16.1180957d0          !ln 1e-7
  real*8,parameter :: nx_bdlt = (20.7232658d0 - nx_blt0)/(nx_nbt-1)   !to ln 1e9
!
contains
!
!
pure function nx_kq(k,nel,x,w) result(kq)
!-- Sobolev line strength K [1/s] of global line k (tau(mu) = K/|Q(mu)|);
!-- x: state populations of the line's element, w: the shares of its
!-- levels in their states (element-local), nel its number density
  integer,intent(in) :: k
  real*8,intent(in) :: nel,x(:),w(:)
  real*8 :: kq
  integer :: i,l,u,lo
  i = nxa_lion(k)
  l = nxa_ion(i)%ll(nxa_lloc(k))
  u = nxa_ion(i)%lu(nxa_lloc(k))
  lo = nxa_el(nxa_lel(k))%lvoff
  kq = nxa_fsob*nxa_lf(k)*(cc/nxa_lnu(k))*nel* &
       (x(nxa_lsl(k))*w(nxa_llv(k)-lo) - nxa_ion(i)%g(l)/nxa_ion(i)%g(u)*x(nxa_lsu(k))*w(nxa_luv(k)-lo))
end function nx_kq
!
!
subroutine nx_beta_table(qb,qv,tab)
!-----------------------------------------------------------------------
! Escape probabilities of a cell with gradients qb = dv/dr, qv = v/r
! (nlteh_beta) tabulated against ln(K/qs); tab(0) holds the large-depth
! limit <|Q(mu)|>/qs (beta -> <|Q|>/K).
!-----------------------------------------------------------------------
  real*8,intent(in) :: qb,qv
  real*8,intent(out) :: tab(0:nx_nbt)
  integer :: i
  real*8 :: qs,mu,q,qmin
  integer,parameter :: nmu = 64
  qs = max(abs(qb),abs(qv),1d-30)
  do i=1,nx_nbt
     tab(i) = nlteh_beta(exp(nx_blt0 + (i-1)*nx_bdlt)*qs,qb,qv)
  enddo
!-- <|Q|> with the floor of nlteh_beta
  qmin = 1d-3*qs
  tab(0) = 0d0
  do i=1,nmu
     mu = (i-0.5d0)/nmu
     q = max(abs(mu*mu*qb + (1d0-mu*mu)*qv),qmin)
     tab(0) = tab(0) + q
  enddo
  tab(0) = tab(0)/nmu/qs
end subroutine nx_beta_table
!
!
pure function nx_beta(kq,qb,qv,tab) result(beta)
!-- escape probability from the cell's table (nx_beta_table)
  real*8,intent(in) :: kq,qb,qv,tab(0:nx_nbt)
  real*8 :: beta,qs,lt,w
  integer :: i
  if(kq <= 0d0) then
     beta = 1d0
     return
  endif
  qs = max(abs(qb),abs(qv),1d-30)
  lt = log(kq/qs)
  if(lt <= nx_blt0) then
     beta = 1d0
  elseif(lt >= nx_blt0 + (nx_nbt-1)*nx_bdlt) then
     beta = min(1d0,tab(0)*qs/kq)
  else
     w = (lt - nx_blt0)/nx_bdlt
     i = min(int(w) + 1,nx_nbt-1)
     w = w - (i-1)
     beta = exp((1d0-w)*log(tab(i)) + w*log(tab(i+1)))
  endif
end function nx_beta
!
!
pure function nx_jline(k,ng,nug,jg,jb,ljb) result(jl)
!-- intensity incident on global line k: resonance estimator if available,
!-- else J of the group nearest to the line (nug: unused, the group of the
!-- line is set by nxa_groups for the same grid)
  integer,intent(in) :: k,ng
  real*8,intent(in) :: nug(ng),jg(ng),jb
  logical,intent(in) :: ljb
  real*8 :: jl
  if(ljb) then
     jl = jb
  else
     jl = jg(nxa_lgn(k))
  endif
end function nx_jline
!
!
subroutine nx_lte(temp,ne,x)
!-----------------------------------------------------------------------
! Saha-Boltzmann populations of all elements at temp and ne (fractions of
! each element; partition functions from the model levels; a super-level
! holds the sum of its levels).
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,ne
  real*8,intent(out) :: x(nxa_nst)
  integer :: iel,j,i,k,o,ns,lo,lv0
  real*8 :: lmax,lsaha
  real*8,allocatable :: lw(:)
!
  do iel=1,nxa_nel
     o = nxa_el(iel)%xoff
     ns = nxa_el(iel)%ns
     lv0 = nxa_el(iel)%lvoff
     allocate(lw(nxa_el(iel)%nlv))
!-- log of n_k/n_top: Saha-Boltzmann chain from the top ion down
     lsaha = 0d0
     do j=nxa_el(iel)%nion,1,-1
        i = nxa_el(iel)%ion(j)
        lo = nxa_ion(i)%lvoff - lv0
!-- n(level k of ion j)/n(ground of next) = ne phi_k T^-1.5 sahac e^{h nu0/kT}
        do k=1,nxa_ion(i)%nlev
           lw(lo+k) = lsaha + log(ne*nxa_ion(i)%g(k)/(2d0*nxa_ion(i)%gion)* &
                nxa_sahac*temp**(-1.5d0)) + hh*nxa_ion(i)%nu0(k)/(kb*temp)
        enddo
!-- relative to the top: next ion's ground state
        lsaha = lw(lo+1)
     enddo
     lmax = max(maxval(lw),0d0)
     x(o+1:o+ns) = 0d0
     do k=1,nxa_el(iel)%nlv
        x(o+nxa_lvs(lv0+k)) = x(o+nxa_lvs(lv0+k)) + exp(lw(k)-lmax)
     enddo
     x(o+ns) = exp(-lmax)
     x(o+1:o+ns) = x(o+1:o+ns)/sum(x(o+1:o+ns))
     deallocate(lw)
  enddo
end subroutine nx_lte
!
!
subroutine nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,r)
!-----------------------------------------------------------------------
! Rate coefficients of element iel that do not depend on the populations
! (type nxrt): the shares of the levels in their states, the state
! energies, electron-impact rates (nxa_coll_el), and per level
! photoionization and stimulated recombination by the Monte Carlo J of the
! groups (below the group grid: J with the shape of B(trad)) and
! spontaneous radiative recombination; then their sums per state.
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,ng
  real*8,intent(in) :: temp,trad,nug(ng),dnug(ng),jg(ng)
  type(nxrt),intent(inout) :: r
  integer :: j,i,k,lv,ig,igred,isub,ns,nl,lv0,kl,s
  integer,parameter :: nsub = 200
  real*8 :: sm,sst,nu,phit,numin,jred,dlt,xr,sig,hkt
!
  ns = nxa_el(iel)%ns
  nl = nxa_el(iel)%nlv
  lv0 = nxa_el(iel)%lvoff
  if(allocated(r%w)) then
     if(size(r%w)/=nl .or. size(r%eps)/=ns) deallocate(r%w,r%wc,r%eps,r%cq,r%eup,r%edn,r%cil,r%c3l, &
          r%gaml,r%gstl,r%arl,r%ci,r%c3,r%gam,r%gst,r%ar)
  endif
  if(.not.allocated(r%w)) allocate(r%w(nl),r%wc(nl),r%eps(ns),r%cq(ns,ns),r%eup(ns,ns),r%edn(ns,ns), &
       r%cil(nl),r%c3l(nl),r%gaml(nl),r%gstl(nl),r%arl(nl),r%ci(ns),r%c3(ns),r%gam(ns),r%gst(ns),r%ar(ns))
!-- shares and energies
  call nxa_wel(iel,temp,r%w)
  do kl=1,nl
     lv = nxa_lvcl(lv0+kl)
     if(lv>0) then
        r%wc(kl) = r%w(lv-lv0)
     else
        r%wc(kl) = 1d0
     endif
  enddo
  call nxa_eps(iel,r%w,r%eps)
!-- collisions
  call nxa_coll_el(iel,temp,r%cq,r%eup,r%edn,r%cil,r%c3l)
!-- radiative bound-free per level
  hkt = hh/(kb*temp)
  igred = minloc(nug,1)
  numin = nug(igred) - 0.5d0*dnug(igred)
  jred = jg(igred)
  do j=1,nxa_el(iel)%nion
     i = nxa_el(iel)%ion(j)
     do k=1,nxa_ion(i)%nlev
        lv = nxa_ion(i)%lvoff + k
        kl = lv - lv0
        r%arl(kl) = nxa_alpha(lv,temp)
!-- photoionization and stimulated recombination (group integrals)
        sm = 0d0
        do ig=1,nxa_g1(lv)
           sm = sm + jg(ig)*nxa_sgh(nxa_sgo(lv)+ig)
        enddo
!-- stimulated recombination: from the edge up to h(nu-nu0) = 40 kT
        sst = 0d0
        do ig=nxa_g1(lv),1,-1
           nu = max(nug(ig),nxa_lvnu0(lv))
           if(hkt*(nu-nxa_lvnu0(lv)) > 40d0) exit
           sst = sst + jg(ig)*nxa_sgh(nxa_sgo(lv)+ig)*exp(-hkt*(nu-nxa_lvnu0(lv)))
        enddo
!-- threshold below the group grid: J extrapolated with the shape of B(trad)
        if(nxa_lvnu0(lv) < numin .and. jred > 0d0) then
           do isub=1,nsub
              nu = nxa_lvnu0(lv)*(numin/nxa_lvnu0(lv))**((isub-0.5d0)/nsub)
              dlt = nxa_lvnu0(lv)*((numin/nxa_lvnu0(lv))**(dble(isub)/nsub) - &
                   (numin/nxa_lvnu0(lv))**(dble(isub-1)/nsub))
              xr = jred*bplanck(nu,trad)/max(bplanck(nug(igred),trad),1d-300)
              sig = nxa_sigma(nxa_ion(i),k,nu)
              sm = sm + xr*sig/(hh*nu)*dlt
              sst = sst + xr*sig/(hh*nu)*dlt*exp(-min(hkt*(nu-nxa_lvnu0(lv)),700d0))
           enddo
        endif
        phit = nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)
        r%gaml(kl) = 4d0*pi*sm
        r%gstl(kl) = 4d0*pi*phit*sst
     enddo
  enddo
!-- sums per state
  r%ci = 0d0; r%c3 = 0d0; r%gam = 0d0; r%gst = 0d0; r%ar = 0d0
  do kl=1,nl
     s = nxa_lvs(lv0+kl)
     r%ci(s) = r%ci(s) + r%w(kl)*r%cil(kl)
     r%gam(s) = r%gam(s) + r%w(kl)*r%gaml(kl)
     r%c3(s) = r%c3(s) + r%wc(kl)*r%c3l(kl)
     r%gst(s) = r%gst(s) + r%wc(kl)*r%gstl(kl)
     r%ar(s) = r%ar(s) + r%wc(kl)*r%arl(kl)
  enddo
end subroutine nx_el_rates
!
!
subroutine nx_el_matrix(iel,ne,r,jl,beta,nti,ntx,m)
!-----------------------------------------------------------------------
! Rate matrix of element iel: m(s,s') = rate s' -> s per particle in s',
! diagonal = -(total out).  r: nx_el_rates; jl, beta: incident intensity
! and escape probability of the element's lines (in the order of the
! global list); nti(s), ntx(s): non-thermal ionization of s and excitation
! from s (its lowest level) into the upper level of the strongest line
! from there [1/s].  Lines within a super-level drop out.
!-----------------------------------------------------------------------
  integer,intent(in) :: iel
  real*8,intent(in) :: ne
  type(nxrt),intent(in) :: r
  real*8,intent(in) :: jl(:),beta(:),nti(:),ntx(:)
  real*8,intent(out) :: m(:,:)
  integer :: j,i,o,k,kg,l,u,sc,ns,s,kl,lv0,k0,lo
!
  ns = nxa_el(iel)%ns
  lv0 = nxa_el(iel)%lvoff
  k0 = nxa_ion(nxa_el(iel)%ion(1))%loff
!-- collisions
  do l=1,ns
     do u=1,ns
        if(u/=l) m(u,l) = ne*r%cq(l,u)
     enddo
     m(l,l) = 0d0
  enddo
  do j=1,nxa_el(iel)%nion
     i = nxa_el(iel)%ion(j)
     o = nxa_ion(i)%soff
     lo = nxa_ion(i)%lvoff - lv0
!-- lines
     do k=1,nxa_ion(i)%nlin
        kg = nxa_ion(i)%loff + k
        kl = kg - k0
        l = nxa_lsl(kg)
        u = nxa_lsu(kg)
        if(l==u) cycle
        m(u,l) = m(u,l) + r%w(nxa_llv(kg)-lv0)*beta(kl)*nxa_lblu(kg)*jl(kl)
        m(l,u) = m(l,u) + r%w(nxa_luv(kg)-lv0)*beta(kl)*(nxa_la(kg) + nxa_lbul(kg)*jl(kl))
     enddo
!-- two-photon
     do k=1,nxa_ion(i)%ntwo
        l = o + nxa_ion(i)%sl(nxa_ion(i)%tl(k))
        u = o + nxa_ion(i)%sl(nxa_ion(i)%tu(k))
        if(l==u) cycle
        m(l,u) = m(l,u) + r%w(lo+nxa_ion(i)%tu(k))*nxa_ion(i)%ta(k)
     enddo
!-- bound-free with the state of level 1 of the next ion (top state for the last)
     sc = nxa_lvc(nxa_ion(i)%lvoff+1)
     do s=o+1,o+nxa_ion(i)%nsl
        m(sc,s) = m(sc,s) + r%gam(s) + ne*r%ci(s) + nti(s)
        m(s,sc) = m(s,sc) + ne*(r%ar(s) + r%gst(s)) + ne*ne*r%c3(s)
     enddo
!-- non-thermal excitation: from the lowest level of s into the upper level
!-- of its strongest line
     do s=o+1,o+nxa_ion(i)%nsl
        if(ntx(s) <= 0d0) cycle
        l = 0
        do k=1,nxa_ion(i)%nlev
           if(o+nxa_ion(i)%sl(k)==s) then
              l = k
              exit
           endif
        enddo
        u = 0
        do kl=1,nxa_ion(i)%nlin
           if(nxa_ion(i)%ll(kl)/=l) cycle
           if(u==0) then
              u = kl
           elseif(nxa_ion(i)%f(kl) > nxa_ion(i)%f(u)) then
              u = kl
           endif
        enddo
        if(u>0) then
           u = o + nxa_ion(i)%sl(nxa_ion(i)%lu(u))
           if(u/=s) m(u,s) = m(u,s) + ntx(s)
        endif
     enddo
  enddo
  do s=1,ns
     m(s,s) = 0d0
     m(s,s) = -sum(m(:,s))
  enddo
end subroutine nx_el_matrix
!
!
subroutine nx_el_step(iel,nel,ne,qb,qv,btab,jl,r,nti,ntx,it,x,beta,d,ierr, &
     ncx,cxf,cxt,cxr)
!-----------------------------------------------------------------------
! One iteration for the populations x of element iel (number density nel)
! at electron density ne: escape probabilities from x, then the linear
! statistical equilibrium; the first three iterations take the new
! solution, later ones its geometric mean with x.  d: largest relative
! change.  jl, beta: incident intensity and escape probability of the
! element's lines.  Optional: ncx further transitions cxf -> cxt at the
! rates cxr per particle [1/s] (charge exchange).
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,it
  real*8,intent(in) :: nel,ne,qb,qv,btab(0:nx_nbt),jl(:)
  type(nxrt),intent(in) :: r
  real*8,intent(in) :: nti(:),ntx(:)
  real*8,intent(inout) :: x(:)
  real*8,intent(out) :: beta(:),d
  integer,intent(out) :: ierr
  integer,intent(in),optional :: ncx,cxf(:),cxt(:)
  real*8,intent(in),optional :: cxr(:)
  integer :: ns,nl,k,kg,k0,info,ir
  real*8 :: m(nxa_el(iel)%ns,nxa_el(iel)%ns),xn(nxa_el(iel)%ns)
  integer :: ipiv(nxa_el(iel)%ns)
!
  ierr = 0
  d = 0d0
  ns = nxa_el(iel)%ns
  nl = size(jl)
  k0 = nxa_ion(nxa_el(iel)%ion(1))%loff
  do k=1,nl
     kg = k0 + k
     beta(k) = nx_beta(nx_kq(kg,nel,x,r%w),qb,qv,btab)
  enddo
  call nx_el_matrix(iel,ne,r,jl,beta,nti,ntx,m)
  if(present(ncx)) then
     do k=1,ncx
        m(cxt(k),cxf(k)) = m(cxt(k),cxf(k)) + cxr(k)
        m(cxf(k),cxf(k)) = m(cxf(k),cxf(k)) - cxr(k)
     enddo
  endif
!-- particle conservation replaces the equation of the most populated state
  ir = maxloc(x,1)
  m(ir,:) = 1d0
  xn = 0d0
  xn(ir) = 1d0
  call dgesv(ns,1,m,ns,ipiv,xn,ns,info)
  if(info/=0) then
     ierr = 1
     return
  endif
  xn = max(xn,1d-300)
  xn = xn/sum(xn)
  d = maxval(abs(xn-x)/(xn+1d-30),mask=xn>1d-20)
  if(it<=3) then
     x = xn
  else
     x = sqrt(x*xn)        !geometric damping keeps the fractions positive
     x = x/sum(x)
  endif
end subroutine nx_el_step
!
!
subroutine nx_solve_cell(temp,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,lsolve, &
     tol,maxit)
!-----------------------------------------------------------------------
! Statistical equilibrium of all NLTE elements in one cell.
!  temp, trad   electron and radiation temperature [K] (trad: shape of J
!               below the group grid)
!  nel(iel)     element number densities [cm^-3]
!  neoth        electrons from the other elements [cm^-3]
!  qb, qv       dv/dr and v/r [1/s]
!  jg           Monte Carlo J per group; jb, ljb: resonance J per global line
!  nti, ntx     non-thermal ionization/excitation rates per state (nxa_nst)
!  x            in: guess, out: populations (fractions of each element)
!  ne           out: electron density;  beta: escape probabilities
!  ierr         1: singular rate matrix, 2: not converged
!  lsolve       (optional) elements to solve; the others keep x (and
!               their electrons)
!  tol, maxit   (optional) convergence: largest relative change of the
!               populations (default 1e-6; n_e to 1e-2 tol) and iterations
!               (default 300)
! The escape probabilities, populations and electron density are iterated
! together (each iteration: one linear solve per element).
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,trad,nel(nxa_nel),neoth,qb,qv
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng),jb(nxa_nlin),nti(nxa_nst),ntx(nxa_nst)
  logical,intent(in) :: ljb(nxa_nlin)
  real*8,intent(inout) :: x(nxa_nst)
  real*8,intent(out) :: ne,beta(nxa_nlin)
  integer,intent(out) :: ierr
  logical,intent(in),optional :: lsolve(nxa_nel)
  real*8,intent(in),optional :: tol
  integer,intent(in),optional :: maxit
  integer :: iel,ns,o,l0,l1,it,ie,k,ncx,mit
  real*8 :: dtol
  real*8 :: neold,nenew,d,dmax,dne,dold,jl(nxa_nlin),btab(0:nx_nbt),wg(nxa_nlv)
  integer,parameter :: mcx = 4096
  integer :: cxf(mcx),cxt(mcx)
  real*8 :: cxr(mcx)
  logical :: lsol(nxa_nel)
  type(nxrt) :: r(nxa_nel)
!
  ierr = 0
  beta = 1d0
  lsol = .true.
  if(present(lsolve)) lsol = lsolve
  dtol = 1d-6
  if(present(tol)) dtol = tol
  mit = 300
  if(present(maxit)) mit = maxit
  call nx_beta_table(qb,qv,btab)
  call nxa_weights(temp,wg)
  do iel=1,nxa_nel
     if(.not.lsol(iel)) cycle
     call nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,r(iel))
  enddo
  do k=1,nxa_nlin
     jl(k) = nx_jline(k,ng,nug,jg,jb(k),ljb(k))
  enddo
  ne = neoth + chg(x)
  dold = 0d0
  do it=1,mit
     dmax = 0d0
     do iel=1,nxa_nel
        if(nel(iel) <= 0d0 .or. .not.lsol(iel)) cycle
        o = nxa_el(iel)%xoff
        ns = nxa_el(iel)%ns
        l0 = nxa_ion(nxa_el(iel)%ion(1))%loff + 1
        l1 = nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%loff + &
             nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%nlin
        ncx = 0
        if(nxa_ncx>0) call nx_cx_rates(iel,temp,nel,x,wg,ncx,cxf,cxt,cxr)
        call nx_el_step(iel,nel(iel),ne,qb,qv,btab,jl(l0:l1),r(iel), &
             nti(o+1:o+ns),ntx(o+1:o+ns),it,x(o+1:o+ns),beta(l0:l1),d,ie, &
             ncx,cxf,cxt,cxr)
        if(ie/=0) then
           ierr = ie
           return
        endif
        dmax = max(dmax,d)
     enddo
     neold = ne
     nenew = neoth + chg(x)
     dne = (nenew - neold)/max(neold,1d-300)
!-- geometric mean if the update overshoots (sign change)
     if(it>3 .and. dne*dold < 0d0) then
        ne = sqrt(neold*nenew)
     else
        ne = nenew
     endif
     dold = dne
     if(dmax < dtol .and. abs(dne) < 1d-2*dtol) exit
  enddo
  if(it>mit) ierr = 2
contains
  real*8 function chg(xx)
!-- free electrons from the NLTE elements [cm^-3]
    real*8,intent(in) :: xx(nxa_nst)
    integer :: je,jo
    chg = 0d0
    do je=1,nxa_nel
       jo = nxa_el(je)%xoff
       chg = chg + nel(je)*sum(nxa_el(je)%chg*xx(jo+1:jo+nxa_el(je)%ns))
    enddo
  end function chg
end subroutine nx_solve_cell
!
subroutine cx_reaction(j,temp,nel,x,wg,af,ar,np,w)
!-----------------------------------------------------------------------
! Charge-exchange reaction j at temperature temp: forward and reverse rate
! coefficients af, ar [cm^3/s], participant densities np(1:4) [cm^-3] and
! the shares w(l,k) of the levels in each participant (Boltzmann).  The
! reverse rate from LTE, n_A n_B/(n_C n_D) = G_A G_B/(G_C G_D)
! exp((chi_B - chi_C)/kT), G = sum g exp(-E/kT) over the term (E above the
! ion's ground level; the top state has g = gtop), chi the ionization
! energies of the ions of B and C.  wg: shares of the levels in their
! states (nxa_weights).
!-----------------------------------------------------------------------
  integer,intent(in) :: j
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst),wg(nxa_nlv)
  real*8,intent(out) :: af,ar,np(4),w(nxa_cxml,4)
  integer :: k,l,ie,io,lv,lo,kk
  real*8 :: kt,g(4),emin(4),ex
  kt = kb*temp
  do k=1,4
     ie = nxa_cx(j)%el(k)
     io = nxa_cx(j)%ion(k)
     lo = nxa_el(ie)%lvoff
     emin(k) = 0d0
     g(k) = 0d0
     np(k) = 0d0
     if(io>0) then
        kk = nxa_ion(io)%lvoff - lo
        emin(k) = minval(nxa_ion(io)%e(nxa_cx(j)%s(1:nxa_cx(j)%nl(k),k)-kk))
        do l=1,nxa_cx(j)%nl(k)
           lv = nxa_cx(j)%s(l,k)
           w(l,k) = nxa_ion(io)%g(lv-kk)*exp(-(nxa_ion(io)%e(lv-kk)-emin(k))/kt)
           g(k) = g(k) + w(l,k)
           np(k) = np(k) + nel(ie)*x(nxa_el(ie)%xoff+nxa_lvs(lo+lv))*wg(lo+lv)
        enddo
     else
        w(1,k) = nxa_el(ie)%gtop
        g(k) = w(1,k)
        np(k) = nel(ie)*x(nxa_el(ie)%xoff+nxa_el(ie)%ns)
     endif
     w(1:nxa_cx(j)%nl(k),k) = w(1:nxa_cx(j)%nl(k),k)/g(k)
  enddo
  af = nxa_cx_coef(j,temp)
  ex = (nxa_ion(nxa_cx(j)%ion(2))%chi - nxa_ion(nxa_cx(j)%ion(3))%chi &
       - emin(1) - emin(2) + emin(3) + emin(4))/kt
  ar = af*g(1)*g(2)/(g(3)*g(4))*exp(min(max(ex,-700d0),700d0))
end subroutine cx_reaction
!
!
subroutine nx_cx_rates(iel,temp,nel,x,wg,n,sf,st,rt)
!-----------------------------------------------------------------------
! Charge-exchange transitions of element iel with the populations x of
! the partner elements: n transitions sf -> st (element-local states) at
! the rates rt per particle in sf [1/s].  Products are distributed over
! the levels of a term with Boltzmann weights.  wg: shares of the levels
! in their states (nxa_weights).
!-----------------------------------------------------------------------
  integer,intent(in) :: iel
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst),wg(nxa_nlv)
  integer,intent(out) :: n,sf(:),st(:)
  real*8,intent(out) :: rt(:)
  integer :: j
  real*8 :: af,ar,np(4),w(nxa_cxml,4)
  n = 0
  do j=1,nxa_ncx
     if(nxa_cx(j)%el(1)/=iel .and. nxa_cx(j)%el(2)/=iel) cycle
     call cx_reaction(j,temp,nel,x,wg,af,ar,np,w)
     if(nxa_cx(j)%el(1)==iel) then
        call add(1,3,af*np(2))     !A -> C: alpha_f n_B
        call add(3,1,ar*np(4))     !C -> A: alpha_r n_D
     else
        call add(2,4,af*np(1))     !B -> D: alpha_f n_A
        call add(4,2,ar*np(3))     !D -> B: alpha_r n_C
     endif
  enddo
contains
  subroutine add(kf,kd,r)
    integer,intent(in) :: kf,kd
    real*8,intent(in) :: r
    integer :: a,b
    real*8 :: wa
    if(r<=0d0) return
    do a=1,nxa_cx(j)%nl(kf)
       wa = 1d0
       if(nxa_cx(j)%ion(kf)>0) wa = wg(nxa_el(iel)%lvoff+nxa_cx(j)%s(a,kf))
       do b=1,nxa_cx(j)%nl(kd)
          if(n>=size(sf)) stop 'nx_cx_rates: too many transitions'
          n = n + 1
          sf(n) = cxstate(kf,a)
          st(n) = cxstate(kd,b)
          rt(n) = r*w(b,kd)*wa
       enddo
    enddo
  end subroutine add
  integer function cxstate(kk,a)
!-- element-local state of level a of participant kk
    integer,intent(in) :: kk,a
    if(nxa_cx(j)%ion(kk)>0) then
       cxstate = nxa_lvs(nxa_el(iel)%lvoff+nxa_cx(j)%s(a,kk))
    else
       cxstate = nxa_el(iel)%ns
    endif
  end function cxstate
end subroutine nx_cx_rates
!
!
function nx_cx_heat(temp,nel,x,wg) result(h)
!-- net heating by charge exchange [erg/s/cm^3]: (forward - reverse) rate
!-- times the energy released by the forward reaction; wg: nxa_weights
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst),wg(nxa_nlv)
  real*8 :: h,af,ar,np(4),w(nxa_cxml,4)
  integer :: j
  h = 0d0
  do j=1,nxa_ncx
     call cx_reaction(j,temp,nel,x,wg,af,ar,np,w)
     h = h + (af*np(1)*np(2) - ar*np(3)*np(4))*nxa_cx(j)%q
  enddo
end function nx_cx_heat
!
!
end module nltexmod
! vim: fdm=marker
