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
! Populations are fractions x of the element's number density.
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
!-- escape probability table of one cell: beta against ln(K/qs) on
!-- [nx_blt0, nx_blt0 + (nx_nbt-1) nx_bdlt], qs = max(|dv/dr|,|v/r|)
  integer,parameter,public :: nx_nbt = 401
  real*8,parameter :: nx_blt0 = -16.1180957d0          !ln 1e-7
  real*8,parameter :: nx_bdlt = (20.7232658d0 - nx_blt0)/(nx_nbt-1)   !to ln 1e9
!
contains
!
!
pure function nx_kq(k,nel,x) result(kq)
!-- Sobolev line strength K [1/s] of global line k (tau(mu) = K/|Q(mu)|);
!-- x: populations of the line's element, nel its number density
  integer,intent(in) :: k
  real*8,intent(in) :: nel,x(:)
  real*8 :: kq
  integer :: i,l,u
  i = nxa_lion(k)
  l = nxa_ion(i)%ll(nxa_lloc(k))
  u = nxa_ion(i)%lu(nxa_lloc(k))
  kq = nxa_fsob*nxa_lf(k)*(cc/nxa_lnu(k))*nel* &
       (x(nxa_lsl(k)) - nxa_ion(i)%g(l)/nxa_ion(i)%g(u)*x(nxa_lsu(k)))
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
! each element; partition functions from the model levels).
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,ne
  real*8,intent(out) :: x(nxa_nst)
  integer :: iel,j,i,k,o,ns
  real*8 :: lw(nxa_nst),lmax,s,lsaha
!
  do iel=1,nxa_nel
     o = nxa_el(iel)%xoff
     ns = nxa_el(iel)%ns
!-- log of n_s/n_top: Saha-Boltzmann chain from the top ion down
     lw(o+ns) = 0d0
     lsaha = 0d0
     do j=nxa_el(iel)%nion,1,-1
        i = nxa_el(iel)%ion(j)
!-- n(level k of ion j)/n(ground of next) = ne phi_k T^-1.5 sahac e^{h nu0/kT}
        do k=1,nxa_ion(i)%nlev
           lw(o+nxa_ion(i)%soff+k) = lsaha + log(ne*nxa_ion(i)%g(k)/(2d0*nxa_ion(i)%gion)* &
                nxa_sahac*temp**(-1.5d0)) + hh*nxa_ion(i)%nu0(k)/(kb*temp)
        enddo
!-- relative to the top: next ion's ground state
        lsaha = lw(o+nxa_ion(i)%soff+1)
     enddo
     lmax = maxval(lw(o+1:o+ns))
     s = sum(exp(lw(o+1:o+ns)-lmax))
     x(o+1:o+ns) = exp(lw(o+1:o+ns)-lmax)/s
  enddo
end subroutine nx_lte
!
!
subroutine nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,cq,ci,c3,gam,gst,ar)
!-----------------------------------------------------------------------
! Rate coefficients of element iel that do not depend on the populations,
! indexed by the element-local state s (levels; the top state has none):
!  cq(s,s')  electron-impact transition s -> s' (same ion) [cm^3/s]
!  ci(s)     electron-impact ionization [cm^3/s]
!  c3(s)     three-body recombination into s per ion [cm^6/s]
!  gam(s)    photoionization rate [1/s]
!  gst(s)    stimulated recombination coefficient into s [cm^3/s]
!  ar(s)     spontaneous radiative recombination coefficient into s [cm^3/s]
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,ng
  real*8,intent(in) :: temp,trad,nug(ng),dnug(ng),jg(ng)
  real*8,intent(out) :: cq(:,:),ci(:),c3(:),gam(:),gst(:),ar(:)
  integer :: j,i,n,o,k,lv,ig,igred,isub
  integer,parameter :: nsub = 200
  real*8 :: s,sst,nu,phit,numin,jred,dlt,xr,sig,hkt
  real*8,allocatable :: qup(:,:),qdn(:,:),cion(:),c3b(:)
!
  cq = 0d0; ci = 0d0; c3 = 0d0; gam = 0d0; gst = 0d0; ar = 0d0
  hkt = hh/(kb*temp)
  igred = minloc(nug,1)
  numin = nug(igred) - 0.5d0*dnug(igred)
  jred = jg(igred)
  do j=1,nxa_el(iel)%nion
     i = nxa_el(iel)%ion(j)
     n = nxa_ion(i)%nlev
     o = nxa_ion(i)%soff
     allocate(qup(n,n),qdn(n,n),cion(n),c3b(n))
     call nxa_coll(nxa_ion(i),temp,qup,qdn,cion,c3b)
!-- qup(l,u) (l<u) and qdn(u,l) fill disjoint triangles: together s -> s'
     cq(o+1:o+n,o+1:o+n) = qup + qdn
     ci(o+1:o+n) = cion
     c3(o+1:o+n) = c3b
     deallocate(qup,qdn,cion,c3b)
     do k=1,n
        lv = nxa_ion(i)%lvoff + k
        ar(o+k) = nxa_alpha(lv,temp)
!-- photoionization and stimulated recombination (group integrals)
        s = 0d0
        do ig=1,nxa_g1(lv)
           s = s + jg(ig)*nxa_sgh(ig,lv)
        enddo
!-- stimulated recombination: from the edge up to h(nu-nu0) = 40 kT
        sst = 0d0
        do ig=nxa_g1(lv),1,-1
           nu = max(nug(ig),nxa_lvnu0(lv))
           if(hkt*(nu-nxa_lvnu0(lv)) > 40d0) exit
           sst = sst + jg(ig)*nxa_sgh(ig,lv)*exp(-hkt*(nu-nxa_lvnu0(lv)))
        enddo
!-- threshold below the group grid: J extrapolated with the shape of B(trad)
        if(nxa_lvnu0(lv) < numin .and. jred > 0d0) then
           do isub=1,nsub
              nu = nxa_lvnu0(lv)*(numin/nxa_lvnu0(lv))**((isub-0.5d0)/nsub)
              dlt = nxa_lvnu0(lv)*((numin/nxa_lvnu0(lv))**(dble(isub)/nsub) - &
                   (numin/nxa_lvnu0(lv))**(dble(isub-1)/nsub))
              xr = jred*bplanck(nu,trad)/max(bplanck(nug(igred),trad),1d-300)
              sig = nxa_sigma(nxa_ion(i),k,nu)
              s = s + xr*sig/(hh*nu)*dlt
              sst = sst + xr*sig/(hh*nu)*dlt*exp(-min(hkt*(nu-nxa_lvnu0(lv)),700d0))
           enddo
        endif
        phit = nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)
        gam(o+k) = 4d0*pi*s
        gst(o+k) = 4d0*pi*phit*sst
     enddo
  enddo
end subroutine nx_el_rates
!
!
subroutine nx_el_matrix(iel,ne,cq,ci,c3,gam,gst,ar,jl,beta,nti,ntx,m)
!-----------------------------------------------------------------------
! Rate matrix of element iel: m(s,s') = rate s' -> s per particle in s',
! diagonal = -(total out).  jl, beta: incident intensity and escape
! probability of the element's lines (in the order of the global list);
! nti(s), ntx(s): non-thermal ionization of s and excitation from s into
! the state's first allowed upper level (resonance line) [1/s].
!-----------------------------------------------------------------------
  integer,intent(in) :: iel
  real*8,intent(in) :: ne,cq(:,:),ci(:),c3(:),gam(:),gst(:),ar(:)
  real*8,intent(in) :: jl(:),beta(:),nti(:),ntx(:)
  real*8,intent(out) :: m(:,:)
  integer :: j,i,n,o,k,kg,l,u,sc,ns,s,kl
!
  ns = nxa_el(iel)%ns
  m = 0d0
  do j=1,nxa_el(iel)%nion
     i = nxa_el(iel)%ion(j)
     n = nxa_ion(i)%nlev
     o = nxa_ion(i)%soff
!-- collisions (all pairs of the ion)
     do l=1,n
        do u=1,n
           if(u/=l) m(o+u,o+l) = m(o+u,o+l) + ne*cq(o+l,o+u)
        enddo
     enddo
!-- lines
     do k=1,nxa_ion(i)%nlin
        kg = nxa_ion(i)%loff + k
        kl = kg - nxa_ion(nxa_el(iel)%ion(1))%loff
        l = o + nxa_ion(i)%ll(k)
        u = o + nxa_ion(i)%lu(k)
        m(u,l) = m(u,l) + beta(kl)*nxa_lblu(kg)*jl(kl)
        m(l,u) = m(l,u) + beta(kl)*(nxa_la(kg) + nxa_lbul(kg)*jl(kl))
     enddo
!-- two-photon
     do k=1,nxa_ion(i)%ntwo
        m(o+nxa_ion(i)%tl(k),o+nxa_ion(i)%tu(k)) = m(o+nxa_ion(i)%tl(k),o+nxa_ion(i)%tu(k)) + &
             nxa_ion(i)%ta(k)
     enddo
!-- bound-free with the ground state of the next ion (top state for the last)
     if(j<nxa_el(iel)%nion) then
        sc = nxa_ion(nxa_el(iel)%ion(j+1))%soff + 1
     else
        sc = ns
     endif
     do k=1,n
        s = o + k
        m(sc,s) = m(sc,s) + gam(s) + ne*ci(s) + nti(s)
        m(s,sc) = m(s,sc) + ne*(ar(s) + gst(s)) + ne*ne*c3(s)
     enddo
!-- non-thermal excitation: into the upper level of the strongest line from s
     do k=1,n
        s = o + k
        if(ntx(s) <= 0d0) cycle
        u = 0
        do kl=1,nxa_ion(i)%nlin
           if(nxa_ion(i)%ll(kl)/=k) cycle
           if(u==0) then
              u = kl
           elseif(nxa_ion(i)%f(kl) > nxa_ion(i)%f(u)) then
              u = kl
           endif
        enddo
        if(u>0) m(o+nxa_ion(i)%lu(u),s) = m(o+nxa_ion(i)%lu(u),s) + ntx(s)
     enddo
  enddo
  do s=1,ns
     m(s,s) = 0d0
     m(s,s) = -sum(m(:,s))
  enddo
end subroutine nx_el_matrix
!
!
subroutine nx_el_step(iel,nel,ne,qb,qv,btab,jl,cq,ci,c3,gam,gst,ar,nti,ntx,it,x,beta,d,ierr, &
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
  real*8,intent(in) :: cq(:,:),ci(:),c3(:),gam(:),gst(:),ar(:),nti(:),ntx(:)
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
     beta(k) = nx_beta(nx_kq(kg,nel,x),qb,qv,btab)
  enddo
  call nx_el_matrix(iel,ne,cq,ci,c3,gam,gst,ar,jl,beta,nti,ntx,m)
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
subroutine nx_solve_cell(temp,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,lsolve)
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
  integer :: iel,ns,o,l0,l1,it,ie,k,ncx
  real*8 :: neold,nenew,d,dmax,dne,dold,jl(nxa_nlin),btab(0:nx_nbt)
  integer,parameter :: mcx = 4096
  integer :: cxf(mcx),cxt(mcx)
  real*8 :: cxr(mcx)
  logical :: lsol(nxa_nel)
  type rt
     real*8,allocatable :: cq(:,:),ci(:),c3(:),gam(:),gst(:),ar(:)
  end type rt
  type(rt) :: r(nxa_nel)
!
  ierr = 0
  beta = 1d0
  lsol = .true.
  if(present(lsolve)) lsol = lsolve
  call nx_beta_table(qb,qv,btab)
  do iel=1,nxa_nel
     if(.not.lsol(iel)) cycle
     ns = nxa_el(iel)%ns
     allocate(r(iel)%cq(ns,ns),r(iel)%ci(ns),r(iel)%c3(ns),r(iel)%gam(ns),r(iel)%gst(ns),r(iel)%ar(ns))
     call nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,r(iel)%cq,r(iel)%ci,r(iel)%c3,r(iel)%gam, &
          r(iel)%gst,r(iel)%ar)
  enddo
  do k=1,nxa_nlin
     jl(k) = nx_jline(k,ng,nug,jg,jb(k),ljb(k))
  enddo
  ne = neoth + chg(x)
  dold = 0d0
  do it=1,300
     dmax = 0d0
     do iel=1,nxa_nel
        if(nel(iel) <= 0d0 .or. .not.lsol(iel)) cycle
        o = nxa_el(iel)%xoff
        ns = nxa_el(iel)%ns
        l0 = nxa_ion(nxa_el(iel)%ion(1))%loff + 1
        l1 = nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%loff + &
             nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%nlin
        ncx = 0
        if(nxa_ncx>0) call nx_cx_rates(iel,temp,nel,x,ncx,cxf,cxt,cxr)
        call nx_el_step(iel,nel(iel),ne,qb,qv,btab,jl(l0:l1), &
             r(iel)%cq,r(iel)%ci,r(iel)%c3,r(iel)%gam,r(iel)%gst,r(iel)%ar, &
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
     if(dmax < 1d-6 .and. abs(dne) < 1d-8) exit
  enddo
  if(it>300) ierr = 2
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
subroutine cx_reaction(j,temp,nel,x,af,ar,np,w)
!-----------------------------------------------------------------------
! Charge-exchange reaction j at temperature temp: forward and reverse rate
! coefficients af, ar [cm^3/s], participant densities np(1:4) [cm^-3] and
! the shares w(l,k) of the levels in each participant (Boltzmann).  The
! reverse rate from LTE, n_A n_B/(n_C n_D) = G_A G_B/(G_C G_D)
! exp((chi_B - chi_C)/kT), G = sum g exp(-E/kT) over the term (E above the
! ion's ground level; the top state has g = gtop), chi the ionization
! energies of the ions of B and C.
!-----------------------------------------------------------------------
  integer,intent(in) :: j
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst)
  real*8,intent(out) :: af,ar,np(4),w(nxa_cxml,4)
  integer :: k,l,m,ie,io,lv
  real*8 :: kt,g(4),emin(4),ex
  kt = kb*temp
  do k=1,4
     ie = nxa_cx(j)%el(k)
     io = nxa_cx(j)%ion(k)
     emin(k) = 0d0
     if(io>0) emin(k) = minval(nxa_ion(io)%e(nxa_cx(j)%s(1:nxa_cx(j)%nl(k),k)-nxa_ion(io)%soff))
     g(k) = 0d0
     np(k) = 0d0
     do l=1,nxa_cx(j)%nl(k)
        m = nxa_cx(j)%s(l,k)
        if(io>0) then
           lv = m - nxa_ion(io)%soff
           w(l,k) = nxa_ion(io)%g(lv)*exp(-(nxa_ion(io)%e(lv)-emin(k))/kt)
        else
           w(l,k) = nxa_el(ie)%gtop
        endif
        g(k) = g(k) + w(l,k)
        np(k) = np(k) + nel(ie)*x(nxa_el(ie)%xoff+m)
     enddo
     w(1:nxa_cx(j)%nl(k),k) = w(1:nxa_cx(j)%nl(k),k)/g(k)
  enddo
  af = nxa_cx_coef(j,temp)
  ex = (nxa_ion(nxa_cx(j)%ion(2))%chi - nxa_ion(nxa_cx(j)%ion(3))%chi &
       - emin(1) - emin(2) + emin(3) + emin(4))/kt
  ar = af*g(1)*g(2)/(g(3)*g(4))*exp(min(max(ex,-700d0),700d0))
end subroutine cx_reaction
!
!
subroutine nx_cx_rates(iel,temp,nel,x,n,sf,st,rt)
!-----------------------------------------------------------------------
! Charge-exchange transitions of element iel with the populations x of
! the partner elements: n transitions sf -> st at the rates rt per
! particle in sf [1/s].  Products are distributed over the levels of a
! term with Boltzmann weights.
!-----------------------------------------------------------------------
  integer,intent(in) :: iel
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst)
  integer,intent(out) :: n,sf(:),st(:)
  real*8,intent(out) :: rt(:)
  integer :: j
  real*8 :: af,ar,np(4),w(nxa_cxml,4)
  n = 0
  do j=1,nxa_ncx
     if(nxa_cx(j)%el(1)/=iel .and. nxa_cx(j)%el(2)/=iel) cycle
     call cx_reaction(j,temp,nel,x,af,ar,np,w)
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
    if(r<=0d0) return
    do a=1,nxa_cx(j)%nl(kf)
       do b=1,nxa_cx(j)%nl(kd)
          if(n>=size(sf)) stop 'nx_cx_rates: too many transitions'
          n = n + 1
          sf(n) = nxa_cx(j)%s(a,kf)
          st(n) = nxa_cx(j)%s(b,kd)
          rt(n) = r*w(b,kd)
       enddo
    enddo
  end subroutine add
end subroutine nx_cx_rates
!
!
function nx_cx_heat(temp,nel,x) result(h)
!-- net heating by charge exchange [erg/s/cm^3]: (forward - reverse) rate
!-- times the energy released by the forward reaction
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst)
  real*8 :: h,af,ar,np(4),w(nxa_cxml,4)
  integer :: j
  h = 0d0
  do j=1,nxa_ncx
     call cx_reaction(j,temp,nel,x,af,ar,np,w)
     h = h + (af*np(1)*np(2) - ar*np(3)*np(4))*nxa_cx(j)%q
  enddo
end function nx_cx_heat
!
!
end module nltexmod
! vim: fdm=marker
