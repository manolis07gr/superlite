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
!    gradient Q(mu) = mu^2 dv/dr + (1-mu^2) v/r (nlteh_beta), plus the
!    probability that a trapped photon is scattered out of the line by an
!    electron (nx_bes; with a microturbulent velocity nx_vtb >= 0)
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
  use nltedismod
  use nltehmod, only:nlteh_beta,bplanck
  use elemdatamod, only:elem_data
  use nltentmod, only:nxn_xn,nxn_xl,nxn_xu,nxn_xr
  implicit none
  private
  real*8,parameter :: pi = 3.141592653589793d0
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
!
  public :: nx_el_rates,nx_el_matrix,nx_el_step,nx_solve_cell,nx_lte,nx_kq,nx_jline
  public :: nx_beta_table,nx_beta,nx_betab_table,nx_bfac,nx_cx_rates,nx_cx_heat
  public :: nx_ytrap,nx_es_table,nx_bes,nx_vdop,nx_jes,nx_wq,nx_wwin,nx_bew,nx_pcont
  public :: nx_ovl_init,nx_ovl_zone,nx_ovl_cell,nx_ovl_terms,nx_ks
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
     real*8,allocatable :: ard(:)     !(nlv) the part of arl into the dissolved continuum below the edge (nx_dis)
!-- the same summed over the levels of each state (with the shares w; the
!-- recombination coefficients per particle in the continuum state)
     real*8,allocatable :: ci(:),c3(:),gam(:),gst(:),ar(:) !(ns)
  end type nxrt
!-- escape probability table of one cell: beta against ln(K/qs) on
!-- [nx_blt0, nx_blt0 + (nx_nbt-1) nx_bdlt], qs = max(|dv/dr|,|v/r|)
  integer,parameter,public :: nx_nbt = 401
!-- accelerated lambda iteration of the bound-free transitions (set by the
!-- caller per cell): the local operator of each group, the share of the
!-- group's intensity that the cell's own emission makes, and the
!-- absorption opacity of the group with which the Monte Carlo J was made
  logical,public :: nx_ali_on = .false.
  real*8,allocatable,public :: nx_ali_lst(:),nx_ali_kap(:) !(ng)
!-- (nx_ldbse) the groups with Lambda* = 1 (local detailed balance) take
!-- J = (kappa_0 J_0 + delta eta)/(kappa_0 + delta kappa) from the bound-free
!-- emission and opacity of the current populations of all elements in
!-- every iteration of the statistical equilibrium (nx_solve_cell), not
!-- the per-level ALI
  logical,public :: nx_ldbse = .false.
!-- (nx_auger) inner-shell photoionization by X-rays followed by an Auger
!-- electron: an ion with at least 4 electrons goes up two stages, a Li-like
!-- one up one (Verner & Yakovlev 1995 core shells, as the X-ray opacity of
!-- GAS/physical_opacity_nlte.f; the photon's energy is the fast electrons');
!-- nx_augr(s) the rate per particle of state s and nx_augt(s) the state of
!-- the element it goes to (0: none), per cell in nx_solve_cell
  logical,public :: nx_auger = .false.
  real*8,allocatable :: nx_augr(:),nx_augsg(:,:)  !(nst), (ng,nion) 4 pi sigma dnu/(h nu)
  integer,allocatable :: nx_augt(:)               !(nst)
  real*8,parameter :: nx_blt0 = -16.1180957d0          !ln 1e-7
  real*8,parameter :: nx_bdlt = (20.7232658d0 - nx_blt0)/(nx_nbt-1)   !to ln 1e9
!-- bins of |mu| (the cosine of the direction to the radius) for the
!-- angular distribution of the radiation of each group (nxt_jmsum), finer
!-- toward mu = 1, where a core is seen from far out
  integer,parameter,public :: nx_nmub = 8
  real*8,parameter,public :: nx_mub(0:nx_nmub) = [0d0,0.5d0,0.75d0,0.88d0,0.94d0,0.97d0,0.985d0,0.995d0,1d0]
!-- electron scattering of line photons trapped in their Sobolev resonance
!-- zones (nx_bes), set by the caller per cell: the microturbulent velocity
!-- [cm/s] added in quadrature to the thermal velocity for the Doppler width
!-- of the lines (< 0: off, point resonances), and whether the scattered
!-- photons must escape the Lorentz wings (nx_wq); nx_esm (when allocated)
!-- the lines that have it (the others keep beta alone)
  real*8,public :: nx_vtb = -1d0
  logical,public :: nx_esdamp = .true.
!-- (diagnostics, set by the caller per cell: nx_dump) the net flows into
!-- the ground state of H per H particle [1/s] from each state of H, the
!-- charge-exchange part, and of Lyman alpha the escape probability, the
!-- incident intensity of the rates and the source function, of the last
!-- step of nx_el_step
  logical,public :: nx_dump = .false.
  integer,public :: nx_dns = 0
  real*8,public :: nx_dflow(64) = 0d0,nx_dcx = 0d0,nx_dlya(9,7) = 0d0
!-- (and with nx_dump_el > 0: the flows into state nx_dump_st of that element
!-- from each of its states, charge exchange apart [1/s per particle of
!-- the element], nx_dflx(1:ns))
  integer,public :: nx_dump_el = 0, nx_dump_st = 0
  real*8,allocatable,public :: nx_dflx(:)
!-- and per element and ion stage j: the flows from ion j to ion j+1 (or the
!-- top state) and back without charge exchange, and the net charge exchange
!-- from j to j+1 [1/s per particle of the element]
  real*8,allocatable,public :: nx_dion(:,:,:)  !(3,nxa_maxion,nxa_nel)
!-- level dissolution (nltedismod), set by the caller per cell: the cell's
!-- parameters (nxd_cell; nx_dis(1) = 0: none)
  real*8,public :: nx_dis(3) = 0d0
  logical,allocatable,public :: nx_esm(:)    !(nlin)
!-- overlapping lines (nx_ovl_init, nx_ovl_cell): a photon trapped in the
!-- resonance zone of a line is absorbed also by the lines whose Doppler
!-- profiles overlap its own (O I 1025.76 in Lyman beta: Bowen
!-- fluorescence).  Partners of a line: the lines within nx_ovl_dv [cm/s]
!-- (0: off), listed per line (pairs nxo_ptr(k) .. nxo_ptr(k+1)-1: partner
!-- nxo_par, reverse pair nxo_rev, the partner's centre nxo_du [cm/s] to
!-- the blue of the line's).  Per cell (nx_ovl_cell; Doppler widths with
!-- the microturbulence nx_ovl_vtb [cm/s], set by the caller): the escape
!-- probability from the zone of the line and its partners relative to the
!-- line's own (nxo_rho; 1 without overlap), the probability per emission
!-- that the partner absorbs the photon (nxo_p), and the bluest line of the
!-- zone (nxo_jcl), whose incident intensity is the zone's.
  real*8,public :: nx_ovl_dv = 0d0, nx_ovl_vtb = 0d0
!-- (nx_ovl_lmax > 0) only the lines whose lower level is among the first
!-- nx_ovl_lmax levels of their ion have partners (the trapped resonance
!-- lines; the line list of an O-star wind has 3.5e7 pairs within 150 km/s)
  integer,public :: nx_ovl_lmax = 0
!-- (partners thinner than nx_ovl_tmin for the mean projected gradient are
!-- left out of a zone, nx_ovl_cell)
  real*8,public :: nx_ovl_tmin = 1d-3
  integer,public :: nxo_npair = 0
  integer,allocatable,public :: nxo_ptr(:)            !(nlin+1)
  integer,allocatable,public :: nxo_par(:),nxo_rev(:) !(npair)
  real*8,allocatable,public :: nxo_du(:)              !(npair)
  real*8,allocatable,public :: nxo_rho(:)             !(nlin)
  real*8,allocatable,public :: nxo_p(:)               !(npair)
  integer,allocatable,public :: nxo_jcl(:)            !(nlin)
!-- (nx_ovl_terms, per line with the current populations: the probability
!-- that a partner absorbs its photon, sum_j P(k,j), and sum_j K_j P(j,k),
!-- over all partners; over the partners of other elements the same
!-- probability and the partners' photons it absorbs, sum_j K_j S_j P(j,k)
!-- (the partners of the element are implicit in its rate matrix); the
!-- incident intensity of its zone.  nxo_use: in the rates of nx_el_step)
  real*8,allocatable,public :: nxo_out(:),nxo_ins(:),nxo_inw(:),nxo_jz(:),nxo_oux(:) !(nlin)
  logical,public :: nxo_use = .false.
!-- (nx_ovl_keep: nx_solve_cell keeps the zones the caller made for the
!-- cell, as the thermal balance's solves at trial temperatures)
  logical,public :: nx_ovl_keep = .false.
!-- (diagnostics: iterations of the last nx_solve_cell)
  integer,public :: nx_nits = 0
!-- continuum destruction of the trapped photons of the lines (nx_pcont):
!-- with nx_pc_on and the absorption opacity of the cell's groups nx_kcg
!-- [cm^-1] (set by the caller per cell, of the populations' volume as K),
!-- a photon emitted in a line is absorbed by the continuum before the line
!-- takes it back with the probability P(kappa_c sqrt(pi) v_D/K); the
!-- escape probability of the rates is beta + (1 - beta) P and their
!-- incident intensity takes P times the J of the line's group.  Per global
!-- line, while a cell is solved: kappa_c sqrt(pi) v_D [1/s] (0: none) and
!-- the J of its group
  logical,public :: nx_pc_on = .false.
  real*8,allocatable,public :: nx_kcg(:)     !(ng)
  real*8,public :: nx_pc_vtb = 0d0           !microturbulence [cm/s] of the Doppler widths
!-- (the J of the groups without the lines' own emission, which the lines
!-- see, set by the caller per cell; not allocated: jg)
  real*8,allocatable,public :: nx_jcg(:)     !(ng)
  logical :: nx_pc_use = .false.
  real*8,allocatable :: nx_pkc(:),nx_pjc(:)  !(nlin)
!-- y(tau) = tau G(tau) of nx_ytrap at log10(tau) = -2, -1.75, ..., 10
  integer,parameter :: nx_nyt = 49
  real*8,parameter :: nx_yt(nx_nyt) = [ &
       3.96959d-05, 1.25043d-04, 3.92700d-04, 1.22671d-03, 3.79593d-03, 1.15524d-02, 3.41525d-02, &
       9.60537d-02, 2.48580d-01, 5.64205d-01, 1.06469d+00, 1.63191d+00, 2.12760d+00, 2.54358d+00, &
       2.90899d+00, 3.23942d+00, 3.54347d+00, 3.82665d+00, 4.09272d+00, 4.34446d+00, 4.58392d+00, &
       4.81274d+00, 5.03220d+00, 5.24335d+00, 5.44706d+00, 5.64405d+00, 5.83494d+00, 6.02025d+00, &
       6.20045d+00, 6.37591d+00, 6.54701d+00, 6.71404d+00, 6.87727d+00, 7.03695d+00, 7.19330d+00, &
       7.34652d+00, 7.49678d+00, 7.64425d+00, 7.78908d+00, 7.93140d+00, 8.07133d+00, 8.20900d+00, &
       8.34450d+00, 8.47794d+00, 8.60940d+00, 8.73897d+00, 8.86673d+00, 8.99275d+00, 9.11710d+00]
!-- W(t) of nx_wwin at log10(t) = -2, -1.75, ..., 10
  real*8,parameter :: nx_wt(nx_nyt) = [ &
       4.97707d-02, 8.81906d-02, 1.55836d-01, 2.74024d-01, 4.77711d-01, 8.20344d-01, 1.37261d+00, &
       2.19866d+00, 3.28489d+00, 4.43746d+00, 5.32069d+00, 5.79788d+00, 6.06356d+00, 6.27177d+00, &
       6.45448d+00, 6.61969d+00, 6.77171d+00, 6.91330d+00, 7.04633d+00, 7.17220d+00, 7.29193d+00, &
       7.40634d+00, 7.51607d+00, 7.62164d+00, 7.72349d+00, 7.82199d+00, 7.91743d+00, 8.01009d+00, &
       8.10018d+00, 8.18791d+00, 8.27346d+00, 8.35697d+00, 8.43859d+00, 8.51843d+00, 8.59660d+00, &
       8.67321d+00, 8.74834d+00, 8.82207d+00, 8.89449d+00, 8.96564d+00, 9.03561d+00, 9.10443d+00, &
       9.17217d+00, 9.23887d+00, 9.30457d+00, 9.36930d+00, 9.43307d+00, 9.49590d+00, 9.55775d+00]
!-- ln(P/(1-P)) of nx_pcont at log10(r) = -12, -11.75, ..., 6
  integer,parameter :: nx_npt = 73
  real*8,parameter :: nx_pt(nx_npt) = [ &
       -25.8513112278d0, -25.2862151137d0, -24.7213469578d0, -24.1567168547d0, -23.5923355866d0, -23.0282146876d0, &
       -22.4643665150d0, -21.9008043306d0, -21.3375423907d0, -20.7745960490d0, -20.2119818723d0, -19.6497177726d0, &
       -19.0878231575d0, -18.5263191018d0, -17.9652285453d0, -17.4045765209d0, -16.8443904178d0, -16.2847002885d0, &
       -15.7255392073d0, -15.1669436901d0, -14.6089541900d0, -14.0516156836d0, -13.4949783686d0, -12.9390984953d0, &
       -12.3840393634d0, -11.8298725176d0, -11.2766791819d0, -10.7245519796d0, -10.1735969815d0, -9.6239361216d0, &
       -9.0757099931d0, -8.5290809834d0, -7.9842366029d0, -7.4413926674d0, -6.9007956607d0, -6.3627230630d0, &
       -5.8274795785d0, -5.2953859598d0, -4.7667554634d0, -4.2418510358d0, -3.7208146919d0, -3.2035606436d0, &
       -2.6896284485d0, -2.1780068829d0, -1.6669704561d0, -1.1540228318d0, -0.6360979226d0, -0.1101543363d0, &
       0.4259185815d0, 0.9725574047d0, 1.5285208350d0, 2.0916324587d0, 2.6596586923d0, 3.2308142121d0, &
       3.8038648502d0, 4.3780282288d0, 4.9528330587d0, 5.5280037303d0, 6.1033817803d0, 6.6788769750d0, &
       7.2544382127d0, 7.8300366422d0, 8.4056560032d0, 8.9812871401d0, 9.5569249008d0, 10.1325663868d0, &
       10.7082099679d0, 11.2838547272d0, 11.8595001491d0, 12.4351459436d0, 13.0107919477d0, 13.5864380695d0, &
       14.1620842575d0]
  real*8,parameter :: sigt = 6.6524587d-25, amu = 1.66053907d-24, me = 9.1093826d-28
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
subroutine nx_betab_table(qb,qv,tab)
!-----------------------------------------------------------------------
! Means of the escape probability over the bins of |mu| (nx_mub; 16-point
! midpoint rule in each, the floor of nlteh_beta), tabulated against
! ln(K/qs) as nx_beta_table; tab(0,b) holds the large-depth limit
! <|Q|>_b/qs.  With the share of each bin in the incident radiation they
! give Sobolev's <beta(mu) I(mu)> (nx_bfac).
!-----------------------------------------------------------------------
  real*8,intent(in) :: qb,qv
  real*8,intent(out) :: tab(0:nx_nbt,nx_nmub)
  integer :: i,j,b
  real*8 :: qs,mu,q(16),qmin,kq,t
  integer,parameter :: nmu = 16
  qs = max(abs(qb),abs(qv),1d-30)
  qmin = 1d-3*qs
  do b=1,nx_nmub
     do j=1,nmu
        mu = nx_mub(b-1) + (j-0.5d0)/nmu*(nx_mub(b) - nx_mub(b-1))
        q(j) = max(abs(mu*mu*qb + (1d0-mu*mu)*qv),qmin)
     enddo
     tab(0,b) = sum(q)/nmu/qs
     do i=1,nx_nbt
        kq = exp(nx_blt0 + (i-1)*nx_bdlt)*qs
        tab(i,b) = 0d0
        do j=1,nmu
           t = kq/q(j)
           if(t < 1d-6) then
              tab(i,b) = tab(i,b) + (1d0 - 0.5d0*t)
           else
              tab(i,b) = tab(i,b) + (1d0 - exp(-t))/t
           endif
        enddo
        tab(i,b) = tab(i,b)/nmu
     enddo
  enddo
end subroutine nx_betab_table
!
!
pure function nx_bfac(kq,qb,qv,tab,fb) result(f)
!-----------------------------------------------------------------------
! Sobolev's weighting of an incident radiation field whose shares in the
! bins of |mu| (nx_mub) are fb (sum 1): <beta(mu) I(mu)>/(beta_bar J) =
! sum_b fb beta_b / sum_b dmu_b beta_b, beta_b the bin means of the escape
! probability (nx_betab_table).  1 for an isotropic field (fb = dmu_b) or
! a thin line; for a thick line in a wind of constant velocity lit by a
! core within mu > 0.97, 0.045.
!-----------------------------------------------------------------------
  real*8,intent(in) :: kq,qb,qv,tab(0:nx_nbt,nx_nmub),fb(nx_nmub)
  real*8 :: f,qs,lt,w,bb(nx_nmub)
  integer :: i,b
  f = 1d0
  if(kq <= 0d0) return
  qs = max(abs(qb),abs(qv),1d-30)
  lt = log(kq/qs)
  if(lt <= nx_blt0) return
  if(lt >= nx_blt0 + (nx_nbt-1)*nx_bdlt) then
     bb = tab(0,:)
  else
     w = (lt - nx_blt0)/nx_bdlt
     i = min(int(w) + 1,nx_nbt-1)
     w = w - (i-1)
     do b=1,nx_nmub
        bb(b) = exp((1d0-w)*log(tab(i,b)) + w*log(tab(i+1,b)))
     enddo
  endif
  f = sum(fb*bb)/max(sum((nx_mub(1:nx_nmub) - nx_mub(0:nx_nmub-1))*bb),1d-300)
end function nx_bfac
!
!
pure function nx_ytrap(tau) result(y)
!-----------------------------------------------------------------------
! Electron scattering of a photon trapped in a line resonance of Sobolev
! depth tau: complete redistribution over a Doppler profile phi(x), x the
! comoving frequency in Doppler widths, which drifts by -s/L along a
! flight (L = v_D/|Q|, the Sobolev length).  Per emission, the probability
! of an electron scattering (opacity sigma, sigma L << 1) during a flight
! that would have ended in re-absorption by the line is sigma L G(tau),
!   G = int dx0 phi(x0) int_{-inf}^{x0} dx exp(-tau (Phi(x0)-Phi(x)))
!       (1 - exp(-tau Phi(x))),         Phi(x) = int_{-inf}^x phi,
! and y = tau G is tabulated (small tau: 0.397 tau^2; large tau:
! 2 sqrt(ln(tau/sqrt(pi))) - 0.35: before it escapes, the photon drifts
! across the core |x| < x_c, beyond which the line is thin).
!-----------------------------------------------------------------------
  real*8,intent(in) :: tau
  real*8 :: y,l,w
  integer :: i
  y = 0d0
  if(tau <= 0d0) return
  l = 4d0*(log10(tau) + 2d0)
  if(l <= 0d0) then
     y = nx_yt(1)*(1d2*tau)**2
  elseif(l >= nx_nyt-1) then
     y = 2d0*sqrt(log(tau/sqrt(pi))) - 0.35d0
  else
     i = int(l) + 1
     w = l - (i-1)
     y = exp((1d0-w)*log(nx_yt(i)) + w*log(nx_yt(i+1)))
  endif
end function nx_ytrap
!
!
pure function nx_wwin(t) result(w)
!-----------------------------------------------------------------------
! (nxt_esrates) Width [Doppler widths] of the window in which the transport
! (nxt_esfar) lets a Sobolev line of depth t = K/|Q| absorb the photons
! that electrons scatter to within five Doppler widths of its centre, where
! the flow drifts them across the profile before an electron scatters them
! again: at x Doppler widths (positive on the side they drift in from) the
! line ahead has the depth t erfc(-x)/2, so
!   W(t) = int_{-5}^{5} (1 - exp(-t erfc(-x)/2)) dx,
! tabulated (small t: 5 t; large t: 5 + sqrt(ln(t/(2 sqrt(pi)))) - 0.107,
! at most 10).  Wider than the core of nx_ytrap: the photons that land
! ahead of the core drift into it (W/y = 1.7 at t = 200, 1.1 at 1e8).
!-----------------------------------------------------------------------
  real*8,intent(in) :: t
  real*8 :: w,l,f
  integer :: i
  w = 0d0
  if(t <= 0d0) return
  l = 4d0*(log10(t) + 2d0)
  if(l <= 0d0) then
     w = nx_wt(1)*(1d2*t)
  elseif(l >= nx_nyt-1) then
     w = min(1d1,5d0 + sqrt(log(t/(2d0*sqrt(pi)))) - 0.107d0)
  else
     i = int(l) + 1
     f = l - (i-1)
     w = exp((1d0-f)*log(nx_wt(i)) + f*log(nx_wt(i+1)))
  endif
end function nx_wwin
!
!
pure function nx_pcont(r) result(p)
!-----------------------------------------------------------------------
! Continuum destruction (Hummer & Rybicki 1985): the probability that a
! photon emitted in a line with a Doppler profile is absorbed by the
! continuum before the line takes it back, for r = kappa_c/kappa_0 the
! ratio of the continuum's absorption opacity to the line's at its centre
! (kappa_0 = K/(sqrt(pi) v_D)):
!   P(r) = (1/sqrt(pi)) int r e^{-x^2}/(r + e^{-x^2}) dx = -r Li_1/2(-1/r),
! tabulated as ln(P/(1-P)) (small r: (2/sqrt(pi)) r sqrt(ln(1/r)); large
! r: 1 - 1/(sqrt(2) r)).  The lines of the metals inside the Lyman continuum of a
! dense wind (Fe III 858, r ~ 5; C II 858, r ~ 0.3) give their photons to
! H I this way within a few scatterings, where their escape from the
! resonance zone is ~1e-6.
!-----------------------------------------------------------------------
  real*8,intent(in) :: r
  real*8 :: p,l,f
  integer :: i
  p = 0d0
  if(r <= 0d0) return
  l = 4d0*(log10(r) + 12d0)
  if(l <= 0d0) then
     p = (r*1d12)*sqrt(log(1d0/r)/log(1d12))/(1d0 + exp(-nx_pt(1)))
  elseif(l >= nx_npt-1) then
     p = 1d0 - 1d0/(sqrt(2d0)*r)
  else
!-- (linear in the logit ln(P/(1-P)): to 7e-4)
     i = int(l) + 1
     f = l - (i-1)
     p = 1d0/(1d0 + exp(-((1d0-f)*nx_pt(i) + f*nx_pt(i+1))))
  endif
end function nx_pcont
!
!
subroutine nx_es_table(qb,qv,tab,wtab)
!-----------------------------------------------------------------------
! <y(K/|Q(mu)|)> over the directions (nx_ytrap) of a cell with gradients
! qb = dv/dr, qv = v/r against ln(K/qs) on the grid of nx_beta_table;
! tab(0) = <|Q|>/qs.  |Q| floored as in nlteh_beta.  wtab (optional): the
! same for the window of the transport's far-side absorption, <W(K/|Q|)>
! (nx_wwin; nx_bew), wtab(0) = tab(0).
!-----------------------------------------------------------------------
  real*8,intent(in) :: qb,qv
  real*8,intent(out) :: tab(0:nx_nbt)
  real*8,intent(out),optional :: wtab(0:nx_nbt)
  integer,parameter :: nmu = 16
  integer :: i,j
  real*8 :: qs,mu,kq,q(nmu)
  qs = max(abs(qb),abs(qv),1d-30)
  do j=1,nmu
     mu = (j-0.5d0)/nmu
     q(j) = max(abs(mu*mu*qb + (1d0-mu*mu)*qv),1d-3*qs)
  enddo
  tab(0) = sum(q)/nmu/qs
  do i=1,nx_nbt
     kq = exp(nx_blt0 + (i-1)*nx_bdlt)*qs
     tab(i) = 0d0
     do j=1,nmu
        tab(i) = tab(i) + nx_ytrap(kq/q(j))
     enddo
     tab(i) = tab(i)/nmu
  enddo
  if(.not.present(wtab)) return
  wtab(0) = tab(0)
  do i=1,nx_nbt
     kq = exp(nx_blt0 + (i-1)*nx_bdlt)*qs
     wtab(i) = 0d0
     do j=1,nmu
        wtab(i) = wtab(i) + nx_wwin(kq/q(j))
     enddo
     wtab(i) = wtab(i)/nmu
  enddo
end subroutine nx_es_table
!
!
pure function nx_bes(kq,sv,wq,beta,qb,qv,tab) result(bes)
!-----------------------------------------------------------------------
! Probability per emission that a photon of a line of Sobolev strength kq
! (and escape probability beta) is scattered out of the line by an
! electron while trapped in its resonance zone: sigma v_D <y(K/|Q(mu)|)>/K
! (nx_ytrap; tab from nx_es_table), sv = sigma v_D [1/s], times the
! probability that the photon, shifted by the electron's Doppler motion,
! escapes the Lorentz wings of the line, exp(-K wq/<|Q|>) with wq from
! nx_wq (Lyman alpha of a dense wind keeps its photons in its damping
! wings); bounded by the probability 1 - beta that the photon is absorbed
! in the line again.  With the escape probability, beta + bes is the
! effective escape probability of the line (as for a continuum in the
! resonance zone, Hummer & Rybicki 1985, ApJ 293, 258).  The saturation
! (1 - exp) keeps it right where electrons scatter a photon before it
! drifts a Doppler width (sigma v_D > |Q|): against a Monte Carlo of one
! resonance (test_estrap) within 2.5 % for sigma v_D/|Q| = 0.03 - 10.
!-----------------------------------------------------------------------
  real*8,intent(in) :: kq,sv,wq,beta,qb,qv,tab(0:nx_nbt)
  real*8 :: bes,qs,lt,w,y
  integer :: i
  bes = 0d0
  if(kq <= 0d0 .or. sv <= 0d0 .or. beta >= 1d0) return
  qs = max(abs(qb),abs(qv),1d-30)
  lt = log(kq/qs)
  if(lt <= nx_blt0) return
  if(lt >= nx_blt0 + (nx_nbt-1)*nx_bdlt) then
     y = nx_ytrap(kq/(tab(0)*qs))
  else
     w = (lt - nx_blt0)/nx_bdlt
     i = min(int(w) + 1,nx_nbt-1)
     w = w - (i-1)
     y = exp((1d0-w)*log(max(tab(i),1d-300)) + w*log(max(tab(i+1),1d-300)))
  endif
  w = sv*y/kq/(1d0-beta)*exp(-min(kq*wq/(tab(0)*qs),700d0))
  if(w < 1d-6) then
     bes = (1d0-beta)*w*(1d0-0.5d0*w)
  else
     bes = (1d0-beta)*(1d0-exp(-w))
  endif
end function nx_bes
!
!
pure function nx_bew(kq,sv,qb,qv,wtab) result(bew)
!-----------------------------------------------------------------------
! (nxt_esrates) The expected absorption by a Sobolev line of strength kq,
! per unit intensity J_es about it (nx_jes), of the photons that electrons
! scatter into its window (nx_wwin) as the transport absorbs them
! (nxt_esfar), as a beta-weighted intensity: sigma v_D <W(K/|Q(mu)|)>/K,
! sv = sigma v_D [1/s], wtab from nx_es_table.  It stands in for the
! transport's tally (nxt_efsum) where that has too few events: Lyman alpha
! in the outer wind of a IIn gets a few hundredths of a packet per cell
! and iteration.
!-----------------------------------------------------------------------
  real*8,intent(in) :: kq,sv,qb,qv,wtab(0:nx_nbt)
  real*8 :: bew,qs,lt,w,ww
  integer :: i
  bew = 0d0
  if(kq <= 0d0 .or. sv <= 0d0) return
  qs = max(abs(qb),abs(qv),1d-30)
  lt = log(kq/qs)
  if(lt <= nx_blt0) return
  if(lt >= nx_blt0 + (nx_nbt-1)*nx_bdlt) then
     ww = nx_wwin(kq/(wtab(0)*qs))
  else
     w = (lt - nx_blt0)/nx_bdlt
     i = min(int(w) + 1,nx_nbt-1)
     w = w - (i-1)
     ww = exp((1d0-w)*log(max(wtab(i),1d-300)) + w*log(max(wtab(i+1),1d-300)))
  endif
  bew = sv*ww/kq
end function nx_bew
!
!
pure function nx_wq(k,temp) result(wq)
!-----------------------------------------------------------------------
! Lorentz-wing parameter of line k at electron temperature temp: a photon
! scattered by an electron leaves the line centre by about Dx = v_e/v_D
! Doppler widths (v_e = sqrt(2kT/m_e)); on its way to the red the
! Lorentz wing of the line (damping a = Gamma lambda/(4 pi v_D)) has the
! depth tau a/(pi Dx) = (K/<|Q|>) wq, wq = Gamma lambda/(4 pi^2 v_e), with
! Gamma the radiative damping (nxa_lgam).
!-----------------------------------------------------------------------
  integer,intent(in) :: k
  real*8,intent(in) :: temp
  real*8 :: wq
  wq = nxa_lgam(k)*(cc/nxa_lnu(k))/(4d0*pi*pi*sqrt(2d0*kb*max(temp,1d0)/me))
end function nx_wq
!
!
pure function nx_jes(k,temp,ng,nug,dnug,jg) result(je)
!-----------------------------------------------------------------------
! Mean intensity from which electrons scatter photons into the core of line
! k: the group J averaged over the thermal-electron redistribution about
! the line (Gaussian of sqrt(2kT/m_e)/c in ln nu, as the electron
! scattering of formalmod).  The line's own escaped photons are part of it.
!-----------------------------------------------------------------------
  integer,intent(in) :: k,ng
  real*8,intent(in) :: temp,nug(ng),dnug(ng),jg(ng)
  real*8 :: je,se,lk,d,w,ws
  integer :: g,g0
  se = sqrt(2d0*kb*max(temp,1d0)/me)/cc
  lk = log(nxa_lnu(k))
  g0 = min(max(nxa_lgn(k),1),ng)
  je = 0d0
  ws = 0d0
!-- groups in order of wavelength: to the blue, then to the red
  do g=g0,1,-1
     d = (log(nug(g)) - lk)/se
     if(d > 3d0) exit
     w = exp(-0.5d0*d*d)*dnug(g)/nug(g)
     je = je + w*jg(g)
     ws = ws + w
  enddo
  do g=g0+1,ng
     d = (lk - log(nug(g)))/se
     if(d > 3d0) exit
     w = exp(-0.5d0*d*d)*dnug(g)/nug(g)
     je = je + w*jg(g)
     ws = ws + w
  enddo
  if(ws > 0d0) then
     je = je/ws
  else
     je = jg(g0)
  endif
end function nx_jes
!
!
subroutine nx_ovl_init(dv)
!-----------------------------------------------------------------------
! Partners of each NLTE line for the overlap of resonance zones: the
! lines whose centres lie within dv [cm/s] (nx_ovl_dv; 0: off), nearest
! first.
!-----------------------------------------------------------------------
  real*8,intent(in) :: dv
  integer :: i,j,k,l,m,n,mm
  integer,allocatable :: cnt(:)
  real*8 :: du
  if(allocated(nxo_ptr)) deallocate(nxo_ptr,nxo_par,nxo_rev,nxo_du)
  if(allocated(nxo_rho)) deallocate(nxo_rho,nxo_p,nxo_jcl,nxo_out,nxo_ins,nxo_inw,nxo_jz,nxo_oux)
  nx_ovl_dv = max(dv,0d0)
  nxo_npair = 0
  nxo_use = .false.
  if(nx_ovl_dv <= 0d0 .or. nxa_nlin < 2) return
  allocate(cnt(nxa_nlin),nxo_ptr(nxa_nlin+1))
!-- (nxa_lsrt: by wavelength, so the separation grows along the list)
  cnt = 0
  do i=1,nxa_nlin
     k = nxa_lsrt(i)
     do j=i+1,nxa_nlin
        l = nxa_lsrt(j)
        if(abs(1d0 - nxa_lnu(l)/nxa_lnu(k))*cc > nx_ovl_dv) exit
        if(.not.(lowl(k) .and. lowl(l))) cycle
        cnt(k) = cnt(k) + 1
        cnt(l) = cnt(l) + 1
     enddo
  enddo
  nxo_ptr(1) = 1
  do k=1,nxa_nlin
     nxo_ptr(k+1) = nxo_ptr(k) + cnt(k)
  enddo
  nxo_npair = nxo_ptr(nxa_nlin+1) - 1
  allocate(nxo_par(max(nxo_npair,1)),nxo_rev(max(nxo_npair,1)),nxo_du(max(nxo_npair,1)))
  cnt = 0
  do i=1,nxa_nlin
     k = nxa_lsrt(i)
     do j=i+1,nxa_nlin
        l = nxa_lsrt(j)
        if(abs(1d0 - nxa_lnu(l)/nxa_lnu(k))*cc > nx_ovl_dv) exit
        if(.not.(lowl(k) .and. lowl(l))) cycle
        m = nxo_ptr(k) + cnt(k)
        nxo_par(m) = l
        nxo_du(m) = (nxa_lnu(l)/nxa_lnu(k) - 1d0)*cc
        cnt(k) = cnt(k) + 1
        m = nxo_ptr(l) + cnt(l)
        nxo_par(m) = k
        nxo_du(m) = (nxa_lnu(k)/nxa_lnu(l) - 1d0)*cc
        cnt(l) = cnt(l) + 1
     enddo
  enddo
!-- nearest first (insertion sort of each list), then the reverse pairs
  do k=1,nxa_nlin
     do m=nxo_ptr(k)+1,nxo_ptr(k+1)-1
        n = nxo_par(m)
        du = nxo_du(m)
        mm = m - 1
        do while(mm >= nxo_ptr(k))
           if(abs(nxo_du(mm)) <= abs(du)) exit
           nxo_par(mm+1) = nxo_par(mm)
           nxo_du(mm+1) = nxo_du(mm)
           mm = mm - 1
        enddo
        nxo_par(mm+1) = n
        nxo_du(mm+1) = du
     enddo
  enddo
  do k=1,nxa_nlin
     do m=nxo_ptr(k),nxo_ptr(k+1)-1
        l = nxo_par(m)
        nxo_rev(m) = 0
        do mm=nxo_ptr(l),nxo_ptr(l+1)-1
           if(nxo_par(mm)==k) then
              nxo_rev(m) = mm
              exit
           endif
        enddo
        if(nxo_rev(m)==0) stop 'nx_ovl_init: pair without its reverse'
     enddo
  enddo
  allocate(nxo_rho(nxa_nlin),nxo_p(max(nxo_npair,1)),nxo_jcl(nxa_nlin))
  allocate(nxo_out(nxa_nlin),nxo_ins(nxa_nlin),nxo_inw(nxa_nlin),nxo_jz(nxa_nlin),nxo_oux(nxa_nlin))
  nxo_rho = 1d0
  nxo_p = 0d0
  nxo_jcl = [(k,k=1,nxa_nlin)]
  nxo_oux = 0d0
  nxo_out = 0d0
  nxo_ins = 0d0
  nxo_inw = 0d0
  nxo_jz = 0d0
contains
  pure logical function lowl(kk)
!-- the lower level of line kk among the first nx_ovl_lmax of its ion
    integer,intent(in) :: kk
    integer :: lv
    lowl = .true.
    if(nx_ovl_lmax <= 0) return
    lv = nxa_llv(kk)
    lowl = lv - nxa_ion(nxa_lvion(lv))%lvoff <= nx_ovl_lmax
  end function lowl
end subroutine nx_ovl_init
!
!
pure subroutine nx_ovl_zone(nl,u,vd,kq,k,qb,qv,esc,bet1,p,aext)
!-----------------------------------------------------------------------
! Resonance zone of nl overlapping lines in the Sobolev approximation with
! complete redistribution over Doppler profiles: centres u [cm/s] (to the
! blue), Doppler velocities vd [cm/s], Sobolev strengths kq [1/s] (>= 0).
! Along a ray the comoving frequency drifts to the red by |Q(mu)| per unit
! length.  On a grid of cells (a quarter of the narrowest Doppler width)
! each line has the depth tau phi du of its normalized profile.  A photon
! of line k emitted at a uniform place in a cell is absorbed in that cell
! or in the cells to its red, by each line in proportion to its depth
! there, or escapes beyond the reddest cell.  Per emission, averaged over
! the directions (|Q| floored at 1e-3 of the larger gradient, as
! nx_es_table): the escape probability esc, the escape probability of the
! line alone bet1 = <(1 - exp(-tau_k))/tau_k>, and the probability p(j)
! that line j absorbs the photon (p(k): the line itself).  esc + sum(p) =
! 1, and with A_k the absorption by line k of photons incident from the
! blue, tau_k (esc + sum_{j/=k} p(j)) = A_k + sum_{j/=k} tau_j p_j(k) per
! direction (detailed balance; test_overlap).  aext (optional): A_j
! averaged over the directions.
!-----------------------------------------------------------------------
  integer,intent(in) :: nl,k
  real*8,intent(in) :: u(nl),vd(nl),kq(nl),qb,qv
  real*8,intent(out) :: esc,bet1,p(nl)
  real*8,intent(out),optional :: aext(nl)
  integer,parameter :: nmu = 8, nmax = 4000
  integer :: n,i,j,im
  real*8 :: umin,umax,du,qs,mu,q,dd,t,s,e1,es,wk,em,tk
  real*8 :: tau(nl),pm(nl),a(nl),d(nl)
  real*8,allocatable :: phi(:,:)
  umin = minval(u - 5d0*vd)
  umax = maxval(u + 5d0*vd)
  du = 0.25d0*minval(vd)
  n = max(ceiling((umax - umin)/du),1)
  if(n > nmax) n = nmax
  du = (umax - umin)/n
  allocate(phi(nl,n))
  do i=1,n
     phi(:,i) = exp(-((umin + (i-0.5d0)*du - u)/vd)**2)
  enddo
  do j=1,nl
     phi(j,:) = phi(j,:)/max(sum(phi(j,:)),1d-300)
  enddo
  qs = max(abs(qb),abs(qv),1d-30)
  esc = 0d0
  bet1 = 0d0
  p = 0d0
  if(present(aext)) aext = 0d0
  do im=1,nmu
     mu = (im-0.5d0)/nmu
     q = max(abs(mu*mu*qb + (1d0-mu*mu)*qv),1d-3*qs)
     tau = max(kq,0d0)/q
!-- from the red end: a(j) the absorption by j of a photon entering the
!-- current cell from its blue side, es its escape
     a = 0d0
     es = 1d0
     em = 0d0
     pm = 0d0
     do i=1,n
        d = tau*phi(:,i)
        dd = sum(d)
        e1 = exp(-dd)
        if(dd > 1d-4) then
           t = (1d0 - e1)/dd
           s = (1d0 - t)/dd
        else
           t = 1d0 - dd*(0.5d0 - dd/6d0)
           s = 0.5d0 - dd/6d0
        endif
        wk = phi(k,i)
        if(wk > 0d0) then
           em = em + wk*t*es
           pm = pm + wk*(s*d + t*a)
        endif
        a = t*d + e1*a
        es = e1*es
     enddo
     esc = esc + em
     p = p + pm
     if(present(aext)) aext = aext + a
     tk = tau(k)
     if(tk > 1d-6) then
        bet1 = bet1 + (1d0 - exp(-tk))/tk
     else
        bet1 = bet1 + 1d0 - 0.5d0*tk
     endif
  enddo
  esc = esc/nmu
  bet1 = bet1/nmu
  p = p/nmu
  if(present(aext)) aext = aext/nmu
end subroutine nx_ovl_zone
!
!
subroutine nx_ovl_cell(temp,nel,x,wg,qb,qv)
!-----------------------------------------------------------------------
! The overlapping resonance zones of a cell (nx_ovl_zone): electron
! temperature temp, element densities nel, populations x (fractions of
! each element), level shares wg, gradients qb = dv/dr, qv = v/r.  For
! each line, its partners whose profiles reach its own (centres within 6
! times the sum of the Doppler widths, at most 64, nearest first) and that
! absorb noticeably: a depth of at least nx_ovl_tmin for the mean
! projected gradient (a partner of depth t takes about t of the photons a
! thick line would otherwise lose by escape, less than t of a thin
! line's): nxo_rho, nxo_p and nxo_jcl.
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,nel(nxa_nel),x(nxa_nst),wg(nxa_nlv),qb,qv
  integer,parameter :: mz = 64
  integer :: k,m,j,iel,nl,ii,o,ns,lo,nv,mm(mz)
  real*8 :: kq(nxa_nlin),vdl(nxa_nel),qtyp,esc,bet1
  real*8 :: uu(mz),vv(mz),kk(mz),pp(mz)
  if(nxo_npair <= 0) return
  qtyp = max((abs(qb) + 2d0*abs(qv))/3d0,1d-3*max(abs(qb),abs(qv)),1d-30)
  do iel=1,nxa_nel
     vdl(iel) = nx_vdop(iel,temp,nx_ovl_vtb)
  enddo
  do k=1,nxa_nlin
     iel = nxa_lel(k)
     o = nxa_el(iel)%xoff
     ns = nxa_el(iel)%ns
     lo = nxa_el(iel)%lvoff
     nv = nxa_el(iel)%nlv
     kq(k) = max(nx_kq(k,nel(iel),x(o+1:o+ns),wg(lo+1:lo+nv)),0d0)
  enddo
  do k=1,nxa_nlin
     nxo_rho(k) = 1d0
     nxo_jcl(k) = k
     if(nxo_ptr(k+1) <= nxo_ptr(k)) cycle
     nxo_p(nxo_ptr(k):nxo_ptr(k+1)-1) = 0d0
     nl = 1
     uu(1) = 0d0
     vv(1) = vdl(nxa_lel(k))
     kk(1) = kq(k)
     do m=nxo_ptr(k),nxo_ptr(k+1)-1
        j = nxo_par(m)
        if(kq(j) < nx_ovl_tmin*qtyp) cycle
        if(abs(nxo_du(m)) > 6d0*(vv(1) + vdl(nxa_lel(j)))) cycle
        if(nl >= mz) exit
        nl = nl + 1
        uu(nl) = nxo_du(m)
        vv(nl) = vdl(nxa_lel(j))
        kk(nl) = kq(j)
        mm(nl) = m
     enddo
     if(nl < 2) cycle
     call nx_ovl_zone(nl,uu(1:nl),vv(1:nl),kk(1:nl),1,qb,qv,esc,bet1,pp(1:nl))
     if(bet1 > 0d0) nxo_rho(k) = esc/bet1
     ii = 1
     do j=2,nl
        nxo_p(mm(j)) = pp(j)
        if(uu(j) > uu(ii)) ii = j
     enddo
     if(ii > 1) nxo_jcl(k) = nxo_par(mm(ii))
  enddo
end subroutine nx_ovl_cell
!
!
subroutine nx_ovl_terms(iel,nel,x,wg,jl)
!-----------------------------------------------------------------------
! The overlap terms of the lines of element iel for the rates (nx_el_step)
! with the current populations x of all elements (level shares wg,
! element densities nel) and the incident intensities jl of all lines:
! nxo_out, nxo_inw, nxo_oux, nxo_ins, nxo_jz (zones of nx_ovl_cell).
!-----------------------------------------------------------------------
  integer,intent(in) :: iel
  real*8,intent(in) :: nel(nxa_nel),x(nxa_nst),wg(nxa_nlv),jl(nxa_nlin)
  integer :: k,m,j,je,l0,l1,o,ns,lo,nv
  real*8 :: pjk
  l0 = nxa_ion(nxa_el(iel)%ion(1))%loff + 1
  l1 = nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%loff + &
       nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%nlin
  do k=l0,l1
     nxo_out(k) = 0d0
     nxo_oux(k) = 0d0
     nxo_ins(k) = 0d0
     nxo_inw(k) = 0d0
     nxo_jz(k) = jl(nxo_jcl(k))
     do m=nxo_ptr(k),nxo_ptr(k+1)-1
        j = nxo_par(m)
        je = nxa_lel(j)
        nxo_out(k) = nxo_out(k) + nxo_p(m)
        if(je /= iel) nxo_oux(k) = nxo_oux(k) + nxo_p(m)
        pjk = nxo_p(nxo_rev(m))
        if(pjk <= 0d0) cycle
        o = nxa_el(je)%xoff
        ns = nxa_el(je)%ns
        lo = nxa_el(je)%lvoff
        nv = nxa_el(je)%nlv
        nxo_inw(k) = nxo_inw(k) + pjk*max(nx_kq(j,nel(je),x(o+1:o+ns),wg(lo+1:lo+nv)),0d0)
        if(je /= iel) nxo_ins(k) = nxo_ins(k) + pjk*nx_ks(j,nel(je),x(o+1:o+ns),wg(lo+1:lo+nv))
     enddo
  enddo
end subroutine nx_ovl_terms
!
!
pure function nx_ks(k,nel,x,w) result(ks)
!-- Sobolev strength times source function of global line k (populations
!-- as nx_kq): fsob f lambda (g_l/g_u) n_u 2 h nu^3/c^2, finite for any
!-- populations; the photons the line emits per unit volume are 4 pi/(h c)
!-- times it
  integer,intent(in) :: k
  real*8,intent(in) :: nel,x(:),w(:)
  real*8 :: ks
  integer :: i,lo
  i = nxa_lion(k)
  lo = nxa_el(nxa_lel(k))%lvoff
  ks = nxa_fsob*nxa_lf(k)*(cc/nxa_lnu(k))*nel*nxa_ion(i)%g(nxa_ion(i)%ll(nxa_lloc(k)))/ &
       nxa_ion(i)%g(nxa_ion(i)%lu(nxa_lloc(k)))*x(nxa_lsu(k))*w(nxa_luv(k)-lo)* &
       2d0*hh*nxa_lnu(k)**3/cc**2
end function nx_ks
!
!
pure function nx_vdop(iel,temp,vtb) result(vd)
!-- Doppler velocity [cm/s] of the lines of element iel at temperature temp
!-- with microturbulence vtb [cm/s]
  integer,intent(in) :: iel
  real*8,intent(in) :: temp,vtb
  real*8 :: vd
  vd = sqrt(2d0*kb*temp/(elem_data(nxa_el(iel)%z)%m*amu) + max(vtb,0d0)**2)
end function nx_vdop
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
! holds the sum of its levels), with the occupation probabilities of the
! levels in the cell of nx_dis (level dissolution).
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
!-- (times the occupation probability, nx_dis)
        do k=1,nxa_ion(i)%nlev
           lw(lo+k) = lsaha + log(ne*nxa_ion(i)%g(k)/(2d0*nxa_ion(i)%gion)* &
                nxa_sahac*temp**(-1.5d0)*nxd_lvw(nxa_ion(i)%lvoff+k,nx_dis)) + hh*nxa_ion(i)%nu0(k)/(kb*temp)
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
subroutine nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,r,xo,nelo,neo)
!-----------------------------------------------------------------------
! Rate coefficients of element iel that do not depend on the populations
! (type nxrt): the shares of the levels in their states, the state
! energies, electron-impact rates (nxa_coll_el), and per level
! photoionization and stimulated recombination by the Monte Carlo J of the
! groups (below the group grid: J with the shape of B(trad)) and
! spontaneous radiative recombination; then their sums per state.
! With nx_ali_on and the populations xo of the element (density nelo) and
! electron density neo with which J was made: accelerated lambda iteration of the bound-free pairs
! (Rybicki & Hummer 1991 with a local operator).  The part of J that the
! cell's own recombination into level l makes, Lambda* eta_l/kappa, is
! taken out of the photoionization of l, and the same number of
! recombinations (those that the level itself re-absorbs, the share
! n_l sigma_l/kappa) out of its recombination: in a thick continuum the
! on-the-spot cycle then drops out of the equations, so that the
! populations follow the radiation from elsewhere at once instead of
! over many iterations.  The solution does not change (both corrections
! are equal at xo).
! Level dissolution (nx_dis, nltedismod): the LTE population of a level
! carries its occupation probability w, so its shares within a super-level,
! the recombination and three-body rates into it, and the electron-impact
! de-excitations (by the ratio of the states' mean w) have it as well; and
! its photoionization and recombination continue below its edge into the
! dissolved upper levels (nxd_gint, the Monte Carlo J of those groups).
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,ng
  real*8,intent(in) :: temp,trad,nug(ng),dnug(ng),jg(ng)
  type(nxrt),intent(inout) :: r
  real*8,intent(in),optional :: xo(:),nelo,neo
  integer :: j,i,k,lv,ig,igred,isub,ns,nl,lv0,kl,s
  integer,parameter :: nsub = 200
  real*8 :: sm,sst,nu,phit,numin,jred,dlt,xr,sig,hkt
  real*8 :: dal,gsf,nlo,arg,wo,nuc,dm,dh,dhs,dsp,sdis
  real*8 :: wocc(nxa_el(iel)%nlv),ws(nxa_el(iel)%ns)
  logical :: lali,ldis,lgo
!
  ns = nxa_el(iel)%ns
  nl = nxa_el(iel)%nlv
  lv0 = nxa_el(iel)%lvoff
  if(allocated(r%w)) then
     if(size(r%w)/=nl .or. size(r%eps)/=ns) deallocate(r%w,r%wc,r%eps,r%cq,r%eup,r%edn,r%cil,r%c3l, &
          r%gaml,r%gstl,r%arl,r%ard,r%ci,r%c3,r%gam,r%gst,r%ar)
  endif
  if(.not.allocated(r%w)) allocate(r%w(nl),r%wc(nl),r%eps(ns),r%cq(ns,ns),r%eup(ns,ns),r%edn(ns,ns), &
       r%cil(nl),r%c3l(nl),r%gaml(nl),r%gstl(nl),r%arl(nl),r%ard(nl),r%ci(ns),r%c3(ns),r%gam(ns),r%gst(ns),r%ar(ns))
!-- shares and energies (level dissolution: the shares of a super-level's
!-- levels in proportion to their LTE populations w g e^{-E/kT}; ws the
!-- states' mean occupation probability)
  call nxa_wel(iel,temp,r%w)
  ldis = nx_dis(1) > 0d0
  wocc = 1d0
  ws = 1d0
  if(ldis) then
     ws = 0d0
     do kl=1,nl
        wocc(kl) = nxd_lvw(lv0+kl,nx_dis)
        ws(nxa_lvs(lv0+kl)) = ws(nxa_lvs(lv0+kl)) + wocc(kl)*r%w(kl)
     enddo
     ws(ns) = 1d0
     do kl=1,nl
        s = nxa_lvs(lv0+kl)
        if(ws(s) > 0d0) r%w(kl) = r%w(kl)*wocc(kl)/ws(s)
     enddo
     where(ws <= 0d0) ws = 1d0
  endif
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
  if(ldis) then
!-- (detailed balance with the dissolved LTE populations: de-excitations
!-- by the ratio of the states' mean w, three-body recombination by w)
     do s=1,ns
        do kl=1,ns
           if(r%eps(s) > r%eps(kl)) then
              r%cq(s,kl) = r%cq(s,kl)*ws(kl)/ws(s)
              r%edn(s,kl) = r%edn(s,kl)*ws(kl)/ws(s)
           endif
        enddo
     enddo
     r%c3l = r%c3l*wocc
  endif
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
!-- dissolved continuum below the edge
        sdis = 0d0
        lgo = ldis
        if(lgo) lgo = nxd_ga(lv) > 0
        if(lgo) then
           nuc = nxd_numin(nxa_lvnu0(lv),nxd_z(lv),nx_dis)
           if(nuc < nxa_lvnu0(lv)) then
              do ig=nxd_ga(lv),min(nxd_gb(lv),ng)
                 if(nug(ig) + 0.5d0*dnug(ig) <= nuc) exit
                 call nxd_gint(lv,ig,nx_dis,temp,nuc,dm,dh,dhs,dsp)
                 sm = sm + jg(ig)*dh
                 sst = sst + jg(ig)*dhs
                 sdis = sdis + dsp
              enddo
           endif
        endif
        wo = wocc(kl)
        r%gaml(kl) = 4d0*pi*sm
        r%gstl(kl) = 4d0*pi*phit*sst*wo
        r%ard(kl) = 4d0*pi*phit*2d0/cc**2*sdis*wo
        r%arl(kl) = r%arl(kl)*wo + r%ard(kl)
     enddo
  enddo
!-- accelerated lambda iteration (per level): dal = recombinations into l per
!-- n_c n_e that l re-absorbs in the cell, 4 pi sum_g Lambda*_g (n_l sigma_g/
!-- kappa_g) eta_g/(n_c n_e h nu) dnu; the photoionization of l by them,
!-- dal n_c n_e/n_l, at the populations xo
  lali = .false.
  if(nx_ali_on .and. present(xo) .and. present(nelo) .and. present(neo)) lali = allocated(nx_ali_lst)
  if(lali) then
     do kl=1,nl
        lv = lv0 + kl
        s = nxa_lvs(lv)
        nlo = nelo*xo(s)*r%w(kl)
        if(nlo <= 0d0) cycle
        phit = nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)
        dal = 0d0
        do ig=min(nxa_g1(lv),ng),1,-1
           nu = max(nug(ig),nxa_lvnu0(lv))
           arg = hkt*(nu-nxa_lvnu0(lv))
           if(arg > 40d0) exit
           if(nx_ali_kap(ig) <= 0d0 .or. nx_ali_lst(ig) <= 0d0) cycle
           if(nx_ldbse .and. nx_ali_lst(ig) >= 0.99d0) cycle
           dal = dal + nx_ali_lst(ig)*min(1d0,nlo*nxa_sgm(nxa_sgo(lv)+ig)/nx_ali_kap(ig))* &
                2d0*hh*nu**3/cc**2*nxa_sgh(nxa_sgo(lv)+ig)*phit*exp(-arg)
        enddo
        dal = 4d0*pi*dal
!-- photoionization of l by the re-absorbed recombinations: dal n_c n_e/n_l
        gsf = dal*nelo*xo(nxa_lvc(lv))*r%wc(kl)*neo/nlo
        if(gsf <= 0d0) cycle
!-- (Monte Carlo noise: never more than the photoionization itself; and
!-- the group sum of dal never more than the recombination coefficient,
!-- which it can exceed by its quadrature error when Lambda* -> 1: both
!-- terms scaled together, so that they stay equal at xo)
        if(gsf > r%gaml(kl)) then
           dal = dal*r%gaml(kl)/gsf
           gsf = r%gaml(kl)
        endif
        if(dal > r%arl(kl)) then
           gsf = gsf*r%arl(kl)/dal
           dal = r%arl(kl)
        endif
        r%gaml(kl) = r%gaml(kl) - gsf
        r%arl(kl) = r%arl(kl) - dal
     enddo
  endif
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
! from there [1/s]; and the non-thermal excitations between states of the
! Spencer-Fano solution of the cell (nltentmod nxn_xl, nxn_xu, nxn_xr,
! global states, per particle of the lower one).  Lines within a
! super-level drop out.
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
!-- (nx_auger) inner-shell photoionization with the Auger electron
  if(nx_auger .and. allocated(nx_augr)) then
     o = nxa_el(iel)%xoff
     do s=1,ns
        if(nx_augt(o+s) > 0 .and. nx_augr(o+s) > 0d0) m(nx_augt(o+s),s) = m(nx_augt(o+s),s) + nx_augr(o+s)
     enddo
  endif
!-- non-thermal excitation of the Spencer-Fano solution
  o = nxa_el(iel)%xoff
  do k=1,nxn_xn
     l = nxn_xl(k) - o
     u = nxn_xu(k) - o
     if(l<1 .or. l>ns .or. u<1 .or. u>ns) cycle
     m(u,l) = m(u,l) + nxn_xr(k)
  enddo
  do s=1,ns
     m(s,s) = 0d0
     m(s,s) = -sum(m(:,s))
  enddo
end subroutine nx_el_matrix
!
!
subroutine nx_el_step(iel,nel,ne,qb,qv,btab,jl,r,nti,ntx,it,x,beta,d,ierr, &
     ncx,cxf,cxt,cxr,xp,dti,sv,etab,wqe,jes,efar,ewtab)
!-----------------------------------------------------------------------
! One iteration for the populations x of element iel (number density nel)
! at electron density ne: escape probabilities from x, then the linear
! statistical equilibrium; the first three iterations take the new
! solution, later ones its geometric mean with x.  d: largest relative
! change.  jl, beta: incident intensity and escape probability of the
! element's lines.  Optional: ncx further transitions cxf -> cxt at the
! rates cxr per particle [1/s] (charge exchange); xp, dti: the fractions
! of the previous model at the same mass coordinate and 1/dt [1/s] for the
! time-dependent equations (backward Euler in the Lagrangian frame,
! dx/dt = M x: (M - 1/dt) x = -xp/dt; particle conservation holds since
! sum xp = 1).  Optional: sv = sigma_es v_D [1/s], etab (nx_es_table), the
! Lorentz-wing parameters wqe (nx_wq) and the intensities jes (nx_jes) of
! the element's lines add the electron scattering of the trapped photons
! (nx_bes): escape probability beta + bes, incident intensity
! (beta Jb + bes Jes)/(beta + bes); only to the lines of nx_esm when it
! is allocated.  efar (optional, >= 0 for a line): the absorption of the
! photons scattered into its core in the transport (nxt_esfar), as a
! beta-weighted intensity in place of bes Jes; -2 (a tally with too few
! events): its expectation nx_bew Jes, with ewtab (nx_es_table).  With nxo_use the
! overlapping lines' terms of nx_ovl_terms (set by nx_solve_cell).
!-----------------------------------------------------------------------
  integer,intent(in) :: iel,it
  real*8,intent(in) :: nel,ne,qb,qv,btab(0:nx_nbt),jl(:)
  type(nxrt),intent(in) :: r
  real*8,intent(in) :: nti(:),ntx(:)
  real*8,intent(inout) :: x(:)
  real*8,intent(out) :: beta(:),d
  integer,intent(out) :: ierr
  integer,intent(in),optional :: ncx,cxf(:),cxt(:)
  real*8,intent(in),optional :: cxr(:),xp(:),dti,sv,etab(0:nx_nbt),wqe(:),jes(:),efar(:),ewtab(0:nx_nbt)
  integer :: ns,nl,k,kg,k0,info,ir,mo,j
  real*8 :: kq,bes,jle(size(jl)),b0,pin,pw,bjl,cx,dlb(4,size(jl)),wl,wu,gr,ees,pc
  logical :: les,lesk,lfar,lfark,lovl,lstv
  real*8 :: m(nxa_el(iel)%ns,nxa_el(iel)%ns),xn(nxa_el(iel)%ns)
  real*8 :: msv(nxa_el(iel)%ns,nxa_el(iel)%ns)
  integer :: ipiv(nxa_el(iel)%ns),i1
  logical :: ltd,ldump
!
  ierr = 0
  d = 0d0
  ns = nxa_el(iel)%ns
  nl = size(jl)
  k0 = nxa_ion(nxa_el(iel)%ion(1))%loff
  les = present(sv) .and. present(etab) .and. present(wqe) .and. present(jes)
  if(les) les = sv > 0d0
  lfar = present(efar)
  ldump = nx_dump
  do k=1,nl
     kg = k0 + k
     kq = nx_kq(kg,nel,x,r%w)
     b0 = nx_beta(kq,qb,qv,btab)
!-- overlapping lines (nxo_use; nx_ovl_terms): escape from the zone
!-- rho beta and absorption by the partners of other elements sum_j P(k,j),
!-- their photons absorbed, and the zone's incident intensity with the
!-- factor of detailed balance over all partners (LTE gives J = S); beta
!-- and beta*J of the rates.  The exchange with the partners of the element
!-- is in its rate matrix (below).
     lovl = nxo_use
     if(lovl) lovl = nxo_out(kg) > 0d0 .or. nxo_inw(kg) > 0d0
     if(lovl) then
        pin = 0d0
        pw = 0d0
        if(kq > 0d0) then
           pin = nxo_ins(kg)/kq
           pw = nxo_inw(kg)/kq
        endif
        beta(k) = nxo_rho(kg)*b0 + nxo_oux(kg)
        bjl = max(nxo_rho(kg)*b0 + nxo_out(kg) - pw,0d0)*nxo_jz(kg) + pin
     else
        beta(k) = b0
        bjl = b0*jl(k)
     endif
     if(ldump) dlb(:,k) = [beta(k),0d0,0d0,0d0]
     lesk = les
     if(lesk .and. allocated(nx_esm)) lesk = nx_esm(kg)
     if(lesk) then
        bes = nx_bes(kq,sv,wqe(k),b0,qb,qv,etab)
!-- (efar is read only when present: no short-circuit evaluation)
        lfark = lfar
        if(lfark) lfark = efar(k) >= 0d0
        lstv = lfar .and. present(ewtab)
        if(lstv) lstv = efar(k) < -1.5d0
        if(lfark) then
           ees = efar(k)
        elseif(lstv) then
!-- (a tally with too few events: its expectation)
           ees = nx_bew(kq,sv,qb,qv,ewtab)*jes(k)
        else
           ees = bes*jes(k)
        endif
        bjl = bjl + ees
        if(ldump) then
           dlb(3,k) = ees
           if(lfar .and. present(ewtab)) dlb(4,k) = nx_bew(kq,sv,qb,qv,ewtab)*jes(k)
        endif
        beta(k) = beta(k) + bes
        if(ldump) dlb(2,k) = bes
     endif
!-- (continuum destruction of the trapped photons, nx_pcont: what the
!-- continuum absorbs leaves the line, and the line takes as much of the
!-- continuum's J of its group)
     if(nx_pc_use .and. kq > 0d0) then
        if(nx_pkc(kg) > 0d0) then
           pc = (1d0 - min(beta(k),1d0))*nx_pcont(nx_pkc(kg)/kq)
           beta(k) = beta(k) + pc
           bjl = bjl + pc*nx_pjc(kg)
        endif
     endif
     jle(k) = bjl/max(beta(k),1d-300)
  enddo
  call nx_el_matrix(iel,ne,r,jle,beta,nti,ntx,m)
!-- overlapping lines of the element: the photons of line k absorbed by its
!-- partner j, n_u(k) A_k P(k,j), de-excite k and excite j (linear in the
!-- population of k's upper state, so the exchange is implicit)
  if(nxo_use) then
     do k=1,nl
        kg = k0 + k
        do mo=nxo_ptr(kg),nxo_ptr(kg+1)-1
           if(nxa_lel(nxo_par(mo)) /= iel .or. nxo_p(mo) <= 0d0) cycle
           cx = nxa_la(kg)*nxo_p(mo)*r%w(nxa_luv(kg)-nxa_el(iel)%lvoff)
           m(nxa_lsl(kg),nxa_lsu(kg)) = m(nxa_lsl(kg),nxa_lsu(kg)) + cx
           m(nxa_lsu(kg),nxa_lsu(kg)) = m(nxa_lsu(kg),nxa_lsu(kg)) - cx
           m(nxa_lsu(nxo_par(mo)),nxa_lsu(kg)) = m(nxa_lsu(nxo_par(mo)),nxa_lsu(kg)) + cx
           m(nxa_lsl(nxo_par(mo)),nxa_lsu(kg)) = m(nxa_lsl(nxo_par(mo)),nxa_lsu(kg)) - cx
        enddo
     enddo
  endif
  if(ldump) msv = m
  if(present(ncx)) then
     do k=1,ncx
        m(cxt(k),cxf(k)) = m(cxt(k),cxf(k)) + cxr(k)
        m(cxf(k),cxf(k)) = m(cxf(k),cxf(k)) - cxr(k)
     enddo
  endif
  ltd = .false.
  if(present(xp) .and. present(dti)) ltd = dti > 0d0
  if(ltd) then
     do k=1,ns
        m(k,k) = m(k,k) - dti
     enddo
  endif
!-- particle conservation replaces the equation of the most populated state
  ir = maxloc(x,1)
  m(ir,:) = 1d0
  xn = 0d0
  if(ltd) xn = -dti*xp
  xn(ir) = 1d0
  call dgesv(ns,1,m,ns,ipiv,xn,ns,info)
  if(info/=0) then
     ierr = 1
     return
  endif
  xn = max(xn,1d-300)
  xn = xn/sum(xn)
  d = maxval(abs(xn-x)/(xn+1d-30),mask=xn>1d-20)
!-- (diagnostics: flows with the new populations, charge exchange apart)
  if(ldump) call dump_flows
  if(ldump .and. iel == nx_dump_el .and. nx_dump_st >= 1 .and. nx_dump_st <= ns) then
     if(allocated(nx_dflx)) deallocate(nx_dflx)
     allocate(nx_dflx(ns))
     do k=1,ns
        nx_dflx(k) = msv(nx_dump_st,k)*xn(k) - msv(k,nx_dump_st)*xn(nx_dump_st)
     enddo
     nx_dflx(nx_dump_st) = 0d0
  endif
  if(ldump .and. nxa_el(iel)%z==1) then
     nx_dns = min(ns,size(nx_dflow))
     nx_dflow = 0d0
     do k=2,nx_dns
        nx_dflow(k) = msv(1,k)*xn(k) - msv(k,1)*xn(1)
     enddo
     nx_dcx = 0d0
     if(present(ncx)) then
        do k=1,ncx
           if(cxt(k)==1) nx_dcx = nx_dcx + cxr(k)*xn(cxf(k))
           if(cxf(k)==1) nx_dcx = nx_dcx - cxr(k)*xn(1)
        enddo
     endif
!-- the Lyman lines (the element's lines from the first state, at most six)
!-- and H alpha (states 2 - 3): wavelength [A], escape probability without
!-- and with the trapped-photon e.s. (bes), the incident intensities of the
!-- resonance estimator and of the e.s. term (per bes) and the one in the
!-- rates, S from the populations, the net radiative rate u -> l per
!-- particle, and the intensity of the e.s. term of the expectation (nx_bew,
!-- per bes; 0 without esrates), which stands in for a tally with too few events
     nx_dlya = 0d0
     j = 0
     do k=1,nl
        kg = k0 + k
        if(nxa_lsl(kg) == 1) then
           if(j >= size(nx_dlya,2) - 1) cycle
           j = j + 1
           i1 = j
        elseif(nxa_lsl(kg) == 2 .and. nxa_lsu(kg) == 3) then
           i1 = size(nx_dlya,2)
        else
           cycle
        endif
        wl = r%w(nxa_llv(kg)-nxa_el(iel)%lvoff)*xn(nxa_lsl(kg))
        wu = r%w(nxa_luv(kg)-nxa_el(iel)%lvoff)*xn(nxa_lsu(kg))
        gr = nxa_ion(nxa_lion(kg))%g(nxa_ion(nxa_lion(kg))%lu(nxa_lloc(kg)))/ &
             nxa_ion(nxa_lion(kg))%g(nxa_ion(nxa_lion(kg))%ll(nxa_lloc(kg)))
        nx_dlya(:,i1) = [1d8*cc/nxa_lnu(kg),dlb(1,k),dlb(1,k)+dlb(2,k),jl(k), &
             dlb(3,k)/max(dlb(2,k),1d-300),jle(k), &
             2d0*hh*nxa_lnu(kg)**3/cc**2/max(gr*wl/max(wu,1d-300) - 1d0,1d-300), &
             beta(k)*((nxa_la(kg) + nxa_lbul(kg)*jle(k))*wu - nxa_lblu(kg)*jle(k)*wl), &
             dlb(4,k)/max(dlb(2,k),1d-300)]
     enddo
  endif
  if(it<=3) then
     x = xn
  else
     x = sqrt(x*xn)        !geometric damping keeps the fractions positive
     x = x/sum(x)
  endif
contains
  subroutine dump_flows
!-- per ion stage: flows to the next stage and back (msv: no charge
!-- exchange), and the net charge exchange to the next stage
    integer :: j,ia,ib,s,t,a1,a2,b1,b2,kk
    if(.not.allocated(nx_dion)) allocate(nx_dion(3,nxa_maxion,nxa_nel))
    nx_dion(:,:,iel) = 0d0
    do j=1,nxa_el(iel)%nion
       ia = nxa_el(iel)%ion(j)
       a1 = nxa_ion(ia)%soff + 1
       a2 = nxa_ion(ia)%soff + nxa_ion(ia)%nsl
       if(j < nxa_el(iel)%nion) then
          ib = nxa_el(iel)%ion(j+1)
          b1 = nxa_ion(ib)%soff + 1
          b2 = nxa_ion(ib)%soff + nxa_ion(ib)%nsl
       else
          b1 = ns
          b2 = ns
       endif
       do s=a1,a2
          do t=b1,b2
             nx_dion(1,j,iel) = nx_dion(1,j,iel) + msv(t,s)*xn(s)
             nx_dion(2,j,iel) = nx_dion(2,j,iel) + msv(s,t)*xn(t)
          enddo
       enddo
       if(present(ncx)) then
          do kk=1,ncx
             if(cxf(kk)>=a1 .and. cxf(kk)<=a2 .and. cxt(kk)>=b1 .and. cxt(kk)<=b2) &
                  nx_dion(3,j,iel) = nx_dion(3,j,iel) + cxr(kk)*xn(cxf(kk))
             if(cxf(kk)>=b1 .and. cxf(kk)<=b2 .and. cxt(kk)>=a1 .and. cxt(kk)<=a2) &
                  nx_dion(3,j,iel) = nx_dion(3,j,iel) - cxr(kk)*xn(cxf(kk))
          enddo
       endif
    enddo
  end subroutine dump_flows
end subroutine nx_el_step
!
!
subroutine nx_solve_cell(temp,trad,nel,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr,lsolve, &
     tol,maxit,xprev,dtinv,efar)
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
!  xprev, dtinv (optional) time-dependent equations: the populations of
!               the previous model at the same mass coordinate and 1/dt
!               [1/s] (nx_el_step)
!  efar         (optional) per global line, the absorption of the photons
!               scattered into its core in the transport (>= 0; -2: a
!               tally with too few events, its expectation; nx_el_step)
! The escape probabilities, populations and electron density are iterated
! together (each iteration: one linear solve per element).  With nx_vtb >=
! 0 the escape probabilities include the electron scattering of trapped
! line photons (nx_bes, with the electron density of the iteration).
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
  real*8,intent(in),optional :: xprev(nxa_nst),dtinv
  real*8,intent(in),optional :: efar(nxa_nlin) !absorption in the line cores (nx_el_step)
  integer :: iel,ns,o,l0,l1,it,ie,k,ncx,mit
  logical :: lovl,lz
  real*8 :: xz(nxa_nst)
  real*8 :: dtol
  real*8 :: neold,nenew,d,dmax,dne,dold,jl(nxa_nlin),btab(0:nx_nbt),wg(nxa_nlv)
  real*8 :: etab(0:nx_nbt),ewtab(0:nx_nbt),sv,wq(nxa_nlin),jes(nxa_nlin)
  integer,parameter :: mcx = 4096
  integer :: cxf(mcx),cxt(mcx)
  real*8 :: cxr(mcx)
  logical :: lsol(nxa_nel)
  type(nxrt) :: r(nxa_nel)
!-- (nx_ldbse) groups with Lambda* = 1, the levels absorbing there, and the
!-- J, opacity, populations and rates with which the solve started
  logical :: lse
  integer :: nse,nsl
  integer,allocatable :: seg(:),sll(:),slk(:),sle(:)
  real*8,allocatable :: sej0(:),sek0(:),sej(:),ssg(:,:),ssh(:,:),sem(:,:),sga0(:),sgs0(:),sphi(:),xs0(:,:)
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
  etab = 0d0
  ewtab = 0d0
  if(nx_vtb >= 0d0) then
     if(present(efar)) then
        call nx_es_table(qb,qv,etab,ewtab)
     else
        call nx_es_table(qb,qv,etab)
     endif
  endif
  call nxd_weights(temp,nx_dis,wg)
  ne = neoth + chg(x)
  do iel=1,nxa_nel
     if(.not.lsol(iel)) cycle
     if(nx_ali_on) then
        o = nxa_el(iel)%xoff
        ns = nxa_el(iel)%ns
        call nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,r(iel),xo=x(o+1:o+ns),nelo=nel(iel),neo=ne)
     else
        call nx_el_rates(iel,temp,trad,ng,nug,dnug,jg,r(iel))
     endif
  enddo
  do k=1,nxa_nlin
     jl(k) = nx_jline(k,ng,nug,jg,jb(k),ljb(k))
  enddo
!-- continuum destruction: kappa_c sqrt(pi) v_D and the J of each line's group
  nx_pc_use = nx_pc_on .and. allocated(nx_kcg)
  if(nx_pc_use) nx_pc_use = size(nx_kcg)==ng
  if(nx_pc_use) then
     if(allocated(nx_pkc)) then
        if(size(nx_pkc)/=nxa_nlin) deallocate(nx_pkc,nx_pjc)
     endif
     if(.not.allocated(nx_pkc)) allocate(nx_pkc(nxa_nlin),nx_pjc(nxa_nlin))
     do k=1,nxa_nlin
        nx_pkc(k) = max(nx_kcg(nxa_lgn(k)),0d0)*sqrt(pi)*nx_vdop(nxa_lel(k),temp,nx_pc_vtb)
        nx_pjc(k) = jg(nxa_lgn(k))
        if(allocated(nx_jcg)) then
           if(size(nx_jcg)==ng) nx_pjc(k) = nx_jcg(nxa_lgn(k))
        endif
     enddo
  endif
  if(nx_vtb >= 0d0) then
     do k=1,nxa_nlin
        wq(k) = 0d0
        if(nx_esdamp) wq(k) = nx_wq(k,temp)
        jes(k) = nx_jes(k,temp,ng,nug,dnug,jg)
     enddo
  else
     wq = 0d0
     jes = 0d0
  endif
  dold = 0d0
  lovl = nxo_npair > 0
  if(nx_auger) call aug_rates
  lse = nx_ldbse .and. nx_ali_on
  if(lse) lse = allocated(nx_ali_lst) .and. allocated(nx_ali_kap)
  if(lse) lse = size(nx_ali_lst)==ng
  if(lse) call se_init
  do it=1,mit
     dmax = 0d0
     if(lse .and. it>1) call se_update
!-- overlapping lines: the zones at the first iteration, and every tenth
!-- while a population has moved by more than 10 % since they were made
!-- (unless the caller keeps its own: nx_ovl_keep); the terms of each
!-- element with the populations of the moment
     lz = lovl .and. .not.nx_ovl_keep .and. it==1
     if(lovl .and. .not.nx_ovl_keep .and. it>1 .and. mod(it,10)==0) &
          lz = any(abs(x - xz) > 0.1d0*xz .and. xz > 1d-10)
     if(lz) then
        call nx_ovl_cell(temp,nel,x,wg,qb,qv)
        xz = x
     endif
     nxo_use = lovl
     do iel=1,nxa_nel
        if(nel(iel) <= 0d0 .or. .not.lsol(iel)) cycle
        o = nxa_el(iel)%xoff
        ns = nxa_el(iel)%ns
        l0 = nxa_ion(nxa_el(iel)%ion(1))%loff + 1
        l1 = nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%loff + &
             nxa_ion(nxa_el(iel)%ion(nxa_el(iel)%nion))%nlin
        if(lovl) call nx_ovl_terms(iel,nel,x,wg,jl)
        ncx = 0
        if(nxa_ncx>0) call nx_cx_rates(iel,temp,nel,x,wg,ncx,cxf,cxt,cxr)
        sv = 0d0
        if(nx_vtb >= 0d0) sv = ne*sigt*nx_vdop(iel,temp,nx_vtb)
        if(present(xprev) .and. present(dtinv)) then
           if(present(efar)) then
              call nx_el_step(iel,nel(iel),ne,qb,qv,btab,jl(l0:l1),r(iel), &
                   nti(o+1:o+ns),ntx(o+1:o+ns),it,x(o+1:o+ns),beta(l0:l1),d,ie, &
                   ncx,cxf,cxt,cxr,xprev(o+1:o+ns),dtinv,sv=sv,etab=etab,wqe=wq(l0:l1),jes=jes(l0:l1), &
                   efar=efar(l0:l1),ewtab=ewtab)
           else
              call nx_el_step(iel,nel(iel),ne,qb,qv,btab,jl(l0:l1),r(iel), &
                   nti(o+1:o+ns),ntx(o+1:o+ns),it,x(o+1:o+ns),beta(l0:l1),d,ie, &
                   ncx,cxf,cxt,cxr,xprev(o+1:o+ns),dtinv,sv=sv,etab=etab,wqe=wq(l0:l1),jes=jes(l0:l1))
           endif
        elseif(present(efar)) then
           call nx_el_step(iel,nel(iel),ne,qb,qv,btab,jl(l0:l1),r(iel), &
                nti(o+1:o+ns),ntx(o+1:o+ns),it,x(o+1:o+ns),beta(l0:l1),d,ie, &
                ncx,cxf,cxt,cxr,sv=sv,etab=etab,wqe=wq(l0:l1),jes=jes(l0:l1),efar=efar(l0:l1),ewtab=ewtab)
        else
           call nx_el_step(iel,nel(iel),ne,qb,qv,btab,jl(l0:l1),r(iel), &
                nti(o+1:o+ns),ntx(o+1:o+ns),it,x(o+1:o+ns),beta(l0:l1),d,ie, &
                ncx,cxf,cxt,cxr,sv=sv,etab=etab,wqe=wq(l0:l1),jes=jes(l0:l1))
        endif
        if(ie/=0) then
           ierr = ie
           nxo_use = .false.
           nx_pc_use = .false.
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
  nxo_use = .false.
  nx_pc_use = .false.
  nx_nits = min(it,mit)
  if(it>mit) ierr = 2
contains
  subroutine aug_rates
!-- (nx_auger) the inner-shell photoionization rate of each ion by the J of
!-- the groups, and the state it goes to (two stages up, one for Li-like)
    use bfxsmod, only:bfxs_inner
    integer :: je,j,i,ne_,g,s,jt,st
    real*8 :: gam,en
    if(.not.allocated(nx_augr)) then
       allocate(nx_augr(nxa_nst),nx_augt(nxa_nst),nx_augsg(ng,nxa_nion))
       nx_augt = 0
       do je=1,nxa_nel
          do j=1,nxa_el(je)%nion
             i = nxa_el(je)%ion(j)
             ne_ = nxa_ion(i)%z - nxa_ion(i)%stage + 1
             do g=1,ng
                en = hh*nug(g)/1.602176634d-12
                nx_augsg(g,i) = 4d0*pi*1d-18*bfxs_inner(nxa_ion(i)%z,ne_,en)*dnug(g)/(hh*nug(g))
             enddo
             if(ne_ < 3) cycle
             st = 2
             if(ne_ == 3) st = 1
             jt = j + st
             if(jt <= nxa_el(je)%nion) then
                s = nxa_el(je)%xoff + nxa_ion(nxa_el(je)%ion(jt))%soff + 1
             elseif(jt == nxa_el(je)%nion + 1) then
                s = nxa_el(je)%xoff + nxa_el(je)%ns
             else
                cycle
             endif
             nx_augt(nxa_el(je)%xoff+nxa_ion(i)%soff+1:nxa_el(je)%xoff+nxa_ion(i)%soff+nxa_ion(i)%nsl) = &
                  s - nxa_el(je)%xoff
          enddo
       enddo
    endif
    if(size(nx_augsg,1) /= ng) return
    nx_augr = 0d0
    do je=1,nxa_nel
       do j=1,nxa_el(je)%nion
          i = nxa_el(je)%ion(j)
          gam = sum(jg*nx_augsg(:,i))
          s = nxa_el(je)%xoff + nxa_ion(i)%soff
          nx_augr(s+1:s+nxa_ion(i)%nsl) = gam
       enddo
    enddo
  end subroutine aug_rates
  subroutine se_init
!-- (nx_ldbse) the groups with Lambda* = 1 and an opacity, and the levels
!-- whose bound-free opacity there is at least 1e-8 of it at the start
    integer :: g,lv,je,kl,n,m
    real*8 :: nl0,nu,arg,hkt,phit
    nse = count(nx_ali_lst >= 0.99d0 .and. nx_ali_kap > 0d0)
    if(nse == 0) then
       lse = .false.
       return
    endif
    allocate(seg(nse),sej0(nse),sek0(nse),sej(nse))
    n = 0
    do g=1,ng
       if(nx_ali_lst(g) >= 0.99d0 .and. nx_ali_kap(g) > 0d0) then
          n = n + 1
          seg(n) = g
       endif
    enddo
    sej0 = jg(seg)
    sek0 = nx_ali_kap(seg)
    hkt = hh/(kb*temp)
!-- levels: count, then fill
    nsl = 0
    do m=1,2
       n = 0
       do je=1,nxa_nel
          if(nel(je) <= 0d0) cycle
          do kl=1,nxa_el(je)%nlv
             lv = nxa_el(je)%lvoff + kl
             if(nxa_g1(lv) < minval(seg)) cycle
             nl0 = nel(je)*x(nxa_el(je)%xoff+nxa_lvs(lv))*r(je)%w(kl)
             if(nl0 <= 0d0) cycle
             if(.not.any([(seg(g) <= nxa_g1(lv) .and. nl0*nxa_sgm(nxa_sgo(lv)+min(seg(g),nxa_g1(lv))) >= &
                  1d-8*sek0(g),g=1,nse)])) cycle
             n = n + 1
             if(m==2) then
                sll(n) = lv
                slk(n) = kl
                sle(n) = je
                phit = nxa_lvphi(lv)*nxa_sahac*temp**(-1.5d0)
                sphi(n) = phit
                sga0(n) = r(je)%gaml(kl)
                sgs0(n) = r(je)%gstl(kl)
                do g=1,nse
                   if(seg(g) > nxa_g1(lv)) cycle
                   ssg(g,n) = nxa_sgm(nxa_sgo(lv)+seg(g))
                   ssh(g,n) = nxa_sgh(nxa_sgo(lv)+seg(g))
                   nu = max(nug(seg(g)),nxa_lvnu0(lv))
                   arg = min(hkt*(nu-nxa_lvnu0(lv)),700d0)
!-- emission per n_c n_e [erg/s/cm^3/Hz/sr] / (n_c n_e): phi sigma 2 h nu^3/c^2 e^{-h(nu-nu0)/kT}
                   sem(g,n) = phit*ssg(g,n)*2d0*hh*nu**3/cc**2*exp(-arg)
                enddo
             endif
          enddo
       enddo
       if(m==1) then
          nsl = n
          allocate(sll(nsl),slk(nsl),sle(nsl),sphi(nsl),sga0(nsl),sgs0(nsl))
          allocate(ssg(nse,nsl),ssh(nse,nsl),sem(nse,nsl))
          ssg = 0d0; ssh = 0d0; sem = 0d0
       endif
    enddo
!-- the opacity and emission of these levels at the start (subtracted below)
    allocate(xs0(nsl,2))
    do n=1,nsl
       call se_pops(n,x,ne,xs0(n,1),xs0(n,2))
    enddo
  end subroutine se_init
  subroutine se_pops(n,xx,nee,nlv,ncne)
!-- the level population and its continuum's n_c n_e
    integer,intent(in) :: n
    real*8,intent(in) :: xx(nxa_nst),nee
    real*8,intent(out) :: nlv,ncne
    integer :: je,kl,lv,o
    je = sle(n)
    kl = slk(n)
    lv = sll(n)
    o = nxa_el(je)%xoff
    nlv = nel(je)*xx(o+nxa_lvs(lv))*r(je)%w(kl)
    ncne = nel(je)*xx(o+nxa_lvc(lv))*r(je)%wc(kl)*nee
  end subroutine se_pops
  subroutine se_update
!-- J of the groups from the current populations; the photoionization and
!-- stimulated recombination of the levels by its change since the start
    integer :: g,n,je,kl,lv,s
    real*8 :: dk(nse),de(nse),dj(nse),nlv,ncne,sm,sst,hkt,nu
    dk = 0d0
    de = 0d0
    do n=1,nsl
       call se_pops(n,x,ne,nlv,ncne)
       dk = dk + (nlv - xs0(n,1))*ssg(:,n)
       de = de + (ncne - xs0(n,2))*sem(:,n)
    enddo
    sej = max((sek0*sej0 + de)/max(sek0 + dk,0.1d0*sek0),0d0)
    dj = sej - sej0
    hkt = hh/(kb*temp)
    do n=1,nsl
       je = sle(n)
       kl = slk(n)
       lv = sll(n)
       sm = 0d0
       sst = 0d0
       do g=1,nse
          if(ssh(g,n) == 0d0) cycle
          nu = max(nug(seg(g)),nxa_lvnu0(lv))
          sm = sm + dj(g)*ssh(g,n)
          sst = sst + dj(g)*ssh(g,n)*exp(-min(hkt*(nu-nxa_lvnu0(lv)),700d0))
       enddo
       r(je)%gaml(kl) = max(sga0(n) + 4d0*pi*sm,0d0)
       r(je)%gstl(kl) = max(sgs0(n) + 4d0*pi*sphi(n)*sst,0d0)
    enddo
!-- the sums per state of the elements concerned
    do je=1,nxa_nel
       if(.not.any(sle==je)) cycle
       r(je)%gam = 0d0
       r(je)%gst = 0d0
       do kl=1,nxa_el(je)%nlv
          lv = nxa_el(je)%lvoff + kl
          s = nxa_lvs(lv)
          r(je)%gam(s) = r(je)%gam(s) + r(je)%w(kl)*r(je)%gaml(kl)
          r(je)%gst(s) = r(je)%gst(s) + r(je)%wc(kl)*r(je)%gstl(kl)
       enddo
    enddo
  end subroutine se_update
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
