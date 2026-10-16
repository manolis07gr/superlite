!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nltehmod
!-----------------------------------------------------------------------
! NLTE hydrogen: n-resolved (l-bundled) bound levels n = 1..nh_nlev and H II.
!
! Atomic data
!  - hydrogenic energies with the reduced-mass Rydberg constant
!  - oscillator strengths and bound-free Gaunt factors: Johnson (1972, ApJ
!    174, 227), accurate to a few per cent for all n
!  - photoionization: Kramers cross section times the Johnson Gaunt factor
!  - radiative recombination: Milne relation (consistent with the cross
!    sections, so photoionization/recombination obey detailed balance)
!  - electron-impact excitation: the l-resolved fits of Atoms/data.nlte.h1
!    summed over final and averaged over initial l (n <= 10); van Regemorter
!    with the effective Gaunt factor for pairs not covered
!  - electron-impact ionization: Seaton's formula (Mihalas 1978, eq. 5-79),
!    three-body recombination from detailed balance
!  - 2s-1s two-photon decay with l-mixed n=2 (A_eff = A_2s * 2/8)
!
! Statistical equilibrium (per cell)
!  - photoionization and stimulated recombination by the Monte Carlo mean
!    intensity J_nu in each group
!  - bound-bound radiative rates in the Sobolev approximation,
!       R_lu = beta*B_lu*Jb,   R_ul = beta*(A_ul + B_ul*Jb),
!    with Jb the intensity incident on the resonance (Monte Carlo
!    resonance estimator) and the angle-averaged escape probability for the
!    local, direction-dependent velocity gradient
!       Q(mu) = mu^2 dv/dr + (1-mu^2) v/r,   tau(mu) = K/|Q(mu)|,
!    which holds for non-homologous and non-monotonic flows (Q < 0 in
!    decelerating directions)
!  - electron density iterated with the H ionization; electrons from other
!    elements are taken from the LTE/nebular EOS
!
! Line wavelengths are in air above 2000 A (as the metal line list and
! observed spectra), vacuum below.
!-----------------------------------------------------------------------
  implicit none
  private
!
  real*8,parameter :: pi = 3.141592653589793d0
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
  real*8,parameter :: me = 9.1093826d-28, ee = 4.80320441d-10, ev = 1.60217653d-12
  real*8,parameter :: rydh = 109677.58d0                 !cm^-1, reduced-mass Rydberg for H
  real*8,parameter :: chih = hh*cc*rydh                  !erg, H ionization energy
  real*8,parameter :: fsob = pi*ee*ee/(me*cc)            !0.02654 cm^2/s
  real*8,parameter :: sahac = 4.1412d-16                 !(h^2/(2 pi me k))^1.5 [cm^3 K^1.5]
  real*8,parameter,public :: a2s = 8.2249d0              !2s-1s two-photon rate [1/s]
  real*8,parameter,public :: nh_chi = chih               !H ionization energy [erg]
!
  integer,public :: nh_nlev = 0        !bound levels
  integer,public :: nh_nlin = 0        !lines (all pairs n < n')
  real*8,allocatable,public :: nh_g(:)      !(nlev) statistical weight 2n^2
  real*8,allocatable,public :: nh_eexc(:)   !(nlev) excitation energy above n=1 [erg]
  real*8,allocatable,public :: nh_eion(:)   !(nlev) ionization energy from n [erg]
  real*8,allocatable,public :: nh_nu0(:)    !(nlev) threshold frequency [Hz]
  integer,allocatable,public :: nh_ll(:),nh_lu(:)  !(nlin) lower/upper level
  real*8,allocatable,public :: nh_nu(:),nh_wl(:)   !(nlin) frequency, wavelength [cm]
  real*8,allocatable,public :: nh_f(:),nh_a(:),nh_bul(:),nh_blu(:) !(nlin)
  integer,allocatable,public :: nh_lin(:,:) !(nlev,nlev) line index of (l,u), 0 if none
  integer,allocatable,public :: nh_srt(:)   !(nlin) line indices sorted by wavelength
  real*8,allocatable,public :: nh_wlsrt(:)  !(nlin) sorted wavelengths [cm]
!-- Gaunt factor coefficients
  real*8,allocatable :: gj0(:),gj1(:),gj2(:)
!-- collisional excitation from file: bundled fit list
  integer :: ncfit = 0
  integer,allocatable :: cf_l(:),cf_u(:)    !n of lower and upper level
  real*8,allocatable :: cf_w(:)             !g_nl/g_n weight
  real*8,allocatable :: cf_c(:,:)           !(4,ncfit) fit coefficients
  real*8,allocatable :: cf_de(:)            !transition energy [eV]
  logical,allocatable :: cf_have(:,:)       !(nlev,nlev) pair covered by the fits
!
  public :: nlteh_init,nlteh_sigma,nlteh_solve_cell,nlteh_lte,nlteh_air
  public :: nlteh_alpha,nlteh_qcoll,nlteh_cion,nlteh_beta,nlteh_gbf
  public :: nlteh_twophoton_y,nlteh_set_collfits,nlteh_cell_rates,nlteh_kq,bplanck
!
contains
!
!
pure function nlteh_air(wlvac) result(wl)
!-- vacuum to air wavelength (Morton 2000), cm; unchanged below 2000 A
  real*8,intent(in) :: wlvac
  real*8 :: wl,s2,n
  if(wlvac < 2d-5) then
     wl = wlvac
     return
  endif
  s2 = (1d-4/wlvac)**2              !(1/lambda[um])^2
  n = 1d0 + 8.34254d-5 + 2.406147d-2/(130d0 - s2) + 1.5998d-4/(38.9d0 - s2)
  wl = wlvac/n
end function nlteh_air
!
!
subroutine nlteh_init(nlev)
!-----------------------------------------------------------------------
! Set up the hydrogen model atom with nlev bound levels.
!-----------------------------------------------------------------------
  integer,intent(in) :: nlev
  integer :: n,l,u,k
  real*8 :: x,g
  integer,allocatable :: idx(:)
!
  nh_nlev = nlev
  nh_nlin = nlev*(nlev-1)/2
  allocate(nh_g(nlev),nh_eexc(nlev),nh_eion(nlev),nh_nu0(nlev))
  allocate(gj0(nlev),gj1(nlev),gj2(nlev))
  do n=1,nlev
     nh_g(n) = 2d0*n*n
     nh_eion(n) = chih/dble(n*n)
     nh_eexc(n) = chih - nh_eion(n)
     nh_nu0(n) = nh_eion(n)/hh
!-- Johnson (1972) Gaunt factor coefficients
     select case(n)
     case(1)
        gj0(n) = 1.1330d0; gj1(n) = -0.4059d0; gj2(n) = 0.07014d0
     case(2)
        gj0(n) = 1.0785d0; gj1(n) = -0.2319d0; gj2(n) = 0.02947d0
     case default
        gj0(n) = 0.9935d0 + 0.2328d0/n - 0.1296d0/n**2
        gj1(n) = -(0.6282d0 - 0.5598d0/n + 0.5299d0/n**2)/n
        gj2(n) = (0.3887d0 - 1.181d0/n + 1.470d0/n**2)/n**2
     end select
  enddo
!
  allocate(nh_ll(nh_nlin),nh_lu(nh_nlin),nh_nu(nh_nlin),nh_wl(nh_nlin))
  allocate(nh_f(nh_nlin),nh_a(nh_nlin),nh_bul(nh_nlin),nh_blu(nh_nlin))
  allocate(nh_lin(nlev,nlev))
  nh_lin = 0
  k = 0
  do l=1,nlev-1
     do u=l+1,nlev
        k = k + 1
        nh_ll(k) = l
        nh_lu(k) = u
        nh_lin(l,u) = k
        nh_lin(u,l) = k
        nh_nu(k) = (nh_eion(l) - nh_eion(u))/hh
        nh_wl(k) = nlteh_air(cc/nh_nu(k))
!-- Johnson (1972) oscillator strength
        x = 1d0 - (dble(l)/dble(u))**2
        g = gj0(l) + gj1(l)/x + gj2(l)/x**2
        nh_f(k) = 32d0/(3d0**1.5d0*pi)*dble(l)/dble(u)**3/x**3*g
!-- Einstein coefficients (J in erg/s/cm^2/Hz/sr)
        nh_a(k) = 8d0*pi**2*ee**2*nh_nu(k)**2/(me*cc**3)*nh_g(l)/nh_g(u)*nh_f(k)
        nh_bul(k) = nh_a(k)*cc**2/(2d0*hh*nh_nu(k)**3)
        nh_blu(k) = nh_bul(k)*nh_g(u)/nh_g(l)
     enddo
  enddo
!-- sorted by wavelength
  allocate(nh_srt(nh_nlin),nh_wlsrt(nh_nlin),idx(nh_nlin))
  idx = [(k,k=1,nh_nlin)]
  call sortidx(nh_wl,idx)
  nh_srt = idx
  nh_wlsrt = nh_wl(idx)
!
contains
  subroutine sortidx(a,ix)
!-- insertion sort of the index ix by a(ix)
    real*8,intent(in) :: a(:)
    integer,intent(inout) :: ix(:)
    integer :: i,j,t
    do i=2,size(ix)
       t = ix(i)
       j = i - 1
       do while(j>=1)
          if(a(ix(j)) <= a(t)) exit
          ix(j+1) = ix(j)
          j = j - 1
       enddo
       ix(j+1) = t
    enddo
  end subroutine sortidx
end subroutine nlteh_init
!
!
subroutine nlteh_set_collfits(nlev_nl,nplev,glev,ncoll,llw,lup,c0,c1,c2,c3,de)
!-----------------------------------------------------------------------
! Bundle the l-resolved electron-impact excitation fits of data.nlte.h1
! into n -> n' rate coefficients (sum over final, g-weighted average over
! initial l).  Called once after read_nlte_data.
!-----------------------------------------------------------------------
  integer,intent(in) :: nlev_nl,ncoll
  integer,intent(in) :: nplev(nlev_nl),llw(ncoll),lup(ncoll)
  real*8,intent(in) :: glev(nlev_nl),c0(ncoll),c1(ncoll),c2(ncoll),c3(ncoll),de(ncoll)
  integer :: i,n,np,k
!
  ncfit = 0
  do i=1,ncoll
     n = nplev(llw(i)); np = nplev(lup(i))
     if(n==np .or. max(n,np)>nh_nlev) cycle
     ncfit = ncfit + 1
  enddo
  allocate(cf_l(ncfit),cf_u(ncfit),cf_w(ncfit),cf_c(4,ncfit),cf_de(ncfit))
  allocate(cf_have(nh_nlev,nh_nlev))
  cf_have = .false.
  k = 0
  do i=1,ncoll
     n = nplev(llw(i)); np = nplev(lup(i))
     if(n==np .or. max(n,np)>nh_nlev) cycle
     k = k + 1
     cf_l(k) = min(n,np)
     cf_u(k) = max(n,np)
     if(n<np) then
        cf_w(k) = glev(llw(i))/nh_g(n)
     else
!-- fit given for the downward direction: weight of its lower (= our lower) level
        cf_w(k) = glev(lup(i))/nh_g(np)
     endif
     cf_c(:,k) = [c0(i),c1(i),c2(i),c3(i)]
     cf_de(k) = de(i)
     cf_have(cf_l(k),cf_u(k)) = .true.
  enddo
end subroutine nlteh_set_collfits
!
!
pure function nlteh_gbf(n,y) result(g)
!-- Johnson (1972) Gaunt factor at nu/nu_n = y (bound-free y>=1)
  integer,intent(in) :: n
  real*8,intent(in) :: y
  real*8 :: g
  g = gj0(n) + gj1(n)/y + gj2(n)/y**2
end function nlteh_gbf
!
!
pure function nlteh_sigma(n,nu) result(s)
!-- photoionization cross section of level n at frequency nu [cm^2]
  integer,intent(in) :: n
  real*8,intent(in) :: nu
  real*8 :: s,y
  y = nu/nh_nu0(n)
  if(y < 1d0) then
     s = 0d0
  else
     s = 7.9064d-18*dble(n)/y**3*nlteh_gbf(n,y)
  endif
end function nlteh_sigma
!
!
pure function nlteh_alpha(n,temp) result(alpha)
!-----------------------------------------------------------------------
! Spontaneous radiative recombination coefficient to level n [cm^3/s],
! Milne relation with the same cross section:
!   alpha_n = 4 pi Phi_n e^{-u} (kT/h) Int_0^inf sigma (2 nu^2/c^2) e^{-y} dy,
!   nu = nu_n + y kT/h,  Phi_n e^{-u} = n^2 sahac T^-1.5 (no overflow)
! Simpson's rule on y in [0,40] (400 intervals).
!-----------------------------------------------------------------------
  integer,intent(in) :: n
  real*8,intent(in) :: temp
  real*8 :: alpha
  integer,parameter :: ni = 400
  real*8,parameter :: ymax = 40d0
  real*8 :: s,nu,y,w,h
  integer :: i
  h = ymax/ni
  s = 0d0
  do i=0,ni
     y = i*h
     w = merge(1d0, merge(4d0, 2d0, mod(i,2)==1), i==0 .or. i==ni)
     nu = nh_nu0(n) + y*kb*temp/hh
     s = s + w*nlteh_sigma(n,nu)*2d0*nu*nu/cc**2*exp(-y)
  enddo
  s = s*h/3d0
  alpha = 4d0*pi*dble(n*n)*sahac*temp**(-1.5d0)*(kb*temp/hh)*s
end function nlteh_alpha
!
!
pure function e1x(x) result(e1)
!-- exponential integral E1(x), x>0 (Abramowitz & Stegun 5.1.53, 5.1.56)
  real*8,intent(in) :: x
  real*8 :: e1
  if(x <= 1d0) then
     e1 = -log(x) - 0.57721566d0 + x*(0.99999193d0 + x*(-0.24991055d0 + &
          x*(0.05519968d0 + x*(-0.00976004d0 + x*0.00107857d0))))
  else
     e1 = exp(-x)/x*(x**4 + 8.5733287401d0*x**3 + 18.0590169730d0*x**2 + &
          8.6347608925d0*x + 0.2677737343d0)/(x**4 + 9.5733223454d0*x**3 + &
          25.6329561486d0*x**2 + 21.0996530827d0*x + 3.9584969228d0)
  endif
end function e1x
!
!
pure subroutine nlteh_qcoll(temp,qup,qdn)
!-----------------------------------------------------------------------
! Electron-impact excitation qup(l,u) and de-excitation qdn(u,l) [cm^3/s]
! for all pairs l < u.  qup includes exp(-dE/kT); qdn is computed without
! it (no overflow at low T).
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp
  real*8,intent(out) :: qup(nh_nlev,nh_nlev),qdn(nh_nlev,nh_nlev)
  real*8 :: kt,x,f,lx,de,omega,gam,s(nh_nlev,nh_nlev)
  integer :: k,i,l,u
  kt = kb*temp/ev
  s = 0d0
!-- fits: q = 1.58e-5 (f/x) e^{-x} (kT[eV])^{-1.5}, f = exp(poly(ln x)), x = dE/kT
  do i=1,ncfit
     x = cf_de(i)/kt
     lx = log(x)
     f = exp(cf_c(1,i) + lx*(cf_c(2,i) + lx*(cf_c(3,i) + lx*cf_c(4,i))))
     s(cf_l(i),cf_u(i)) = s(cf_l(i),cf_u(i)) + cf_w(i)*1.58d-5*(f/x)*kt**(-1.5d0)
  enddo
  qup = 0d0
  qdn = 0d0
  do l=1,nh_nlev-1
     do u=l+1,nh_nlev
        de = (nh_eion(l) - nh_eion(u))/ev
        x = de/kt
        if(ncfit>0) then
           if(cf_have(l,u)) then
              qdn(u,l) = s(l,u)*nh_g(l)/nh_g(u)
              qup(l,u) = s(l,u)*exp(-x)
              cycle
           endif
        endif
!-- van Regemorter with the effective Gaunt factor for neutrals
        k = nh_lin(l,u)
        if(x > 50d0) then
           gam = 0.2d0
        else
           gam = max(0.2d0, 0.276d0*exp(x)*e1x(x))
        endif
        omega = 14.5d0*nh_f(k)*(chih/ev/de)**2*x*gam
        qdn(u,l) = 5.465d-11*sqrt(temp)*omega*nh_g(l)/nh_g(u)
        qup(l,u) = 5.465d-11*sqrt(temp)*omega*exp(-x)
     enddo
  enddo
end subroutine nlteh_qcoll
!
!
pure subroutine nlteh_cion(n,temp,cup,c3b)
!-----------------------------------------------------------------------
! Electron-impact ionization from n [cm^3/s] (Seaton; Mihalas eq. 5-79)
! and the three-body recombination coefficient per ion [cm^6/s] from
! detailed balance (rate per ion = ne^2 * c3b).
!-----------------------------------------------------------------------
  integer,intent(in) :: n
  real*8,intent(in) :: temp
  real*8,intent(out) :: cup,c3b
  real*8 :: u,base
  u = nh_eion(n)/(kb*temp)
  base = 1.55d13/sqrt(temp)*0.1d0*nlteh_sigma(n,nh_nu0(n))/u
  cup = base*exp(-u)
  c3b = base*dble(n*n)*sahac*temp**(-1.5d0)
end subroutine nlteh_cion
!
!
pure function nlteh_beta(kq,qb,qv) result(beta)
!-----------------------------------------------------------------------
! Angle-averaged Sobolev escape probability for tau(mu) = kq/|Q(mu)|,
! Q(mu) = mu^2*qb + (1-mu^2)*qv  (qb = dv/dr, qv = v/r), averaged over mu in
! [0,1] (the integrand is even in mu) with a 64-point midpoint rule, which
! copes with the kink where Q changes sign in decelerating flows.  |Q| is
! floored at 1e-3 of the largest gradient: where the projected gradient
! vanishes the line is optically thick, and the floor keeps the Sobolev
! length finite.
!-----------------------------------------------------------------------
  real*8,intent(in) :: kq,qb,qv
  real*8 :: beta
  real*8 :: q,t,qmin,mu
  integer :: i
  integer,parameter :: nmu = 64
  if(kq <= 0d0) then
     beta = 1d0          !no net absorption (or a maser): optically thin
     return
  endif
  qmin = 1d-3*max(abs(qb),abs(qv),1d-30)
  beta = 0d0
  do i=1,nmu
     mu = (i-0.5d0)/nmu
     q = max(abs(mu*mu*qb + (1d0-mu*mu)*qv), qmin)
     t = kq/q
     if(t < 1d-6) then
        beta = beta + (1d0 - 0.5d0*t)
     else
        beta = beta + (1d0 - exp(-t))/t
     endif
  enddo
  beta = beta/nmu
end function nlteh_beta
!
!
pure function nlteh_twophoton_y(r) result(y)
!-----------------------------------------------------------------------
! Fraction y = nu/nu(Ly-alpha) of a 2s-1s two-photon emission, sampled from
! the photon energy distribution y*psi(y) with psi the Nussbaumer & Schmutz
! (1984) fit, by inversion of a tabulated cumulative distribution.
! r: uniform deviate.
!-----------------------------------------------------------------------
  real*8,intent(in) :: r
  real*8 :: y
  integer,parameter :: nt = 200
  real*8 :: cdf(0:nt),yy,p,a
  integer :: i
  cdf(0) = 0d0
  do i=1,nt
     yy = (i-0.5d0)/nt
     a = yy*(1d0-yy)
     p = (a*(1d0-(4d0*a)**0.8d0) + 0.88d0*a**1.53d0*(4d0*a)**0.8d0)*yy
     cdf(i) = cdf(i-1) + p
  enddo
  cdf = cdf/cdf(nt)
  do i=1,nt
     if(cdf(i) >= r) exit
  enddo
  i = min(i,nt)
  if(cdf(i) > cdf(i-1)) then
     y = (i-1 + (r-cdf(i-1))/(cdf(i)-cdf(i-1)))/dble(nt)
  else
     y = (i-0.5d0)/nt
  endif
  y = min(max(y,1d-4),1d0-1d-4)
end function nlteh_twophoton_y
!
!
pure function bplanck(nu,t) result(b)
!-- Planck function per unit frequency [erg/s/cm^2/Hz/sr]
  real*8,intent(in) :: nu,t
  real*8 :: b,x
  x = hh*nu/(kb*t)
  if(x > 700d0) then
     b = 0d0
  elseif(x < 1d-6) then
     b = 2d0*nu*nu*kb*t/cc**2
  else
     b = 2d0*hh*nu**3/cc**2/(exp(x) - 1d0)
  endif
end function bplanck
!
!
pure subroutine nlteh_lte(temp,ne,x)
!-- Saha-Boltzmann fractions x(0:nlev) of n_H (0 = H II) at given ne
  real*8,intent(in) :: temp,ne
  real*8,intent(out) :: x(0:)
  real*8 :: lphi(nh_nlev),lmax,s
  integer :: n
  do n=1,nh_nlev
!-- log of n_n/(ne n_+)
     lphi(n) = log(dble(n*n)*sahac*temp**(-1.5d0)) + nh_eion(n)/(kb*temp) + log(ne)
  enddo
  lmax = max(0d0,maxval(lphi))
  s = exp(-lmax)
  do n=1,nh_nlev
     s = s + exp(lphi(n)-lmax)
  enddo
  x(0) = exp(-lmax)/s
  do n=1,nh_nlev
     x(n) = exp(lphi(n)-lmax)/s
  enddo
end subroutine nlteh_lte
!
!
subroutine nlteh_cell_rates(temp,trad,ng,nug,dnug,jg,jb,ljb,nonth, &
     gam,alp,alps,cup,c3b,qup,qdn,jl)
!-----------------------------------------------------------------------
! Rate coefficients of one cell that do not depend on the populations.
!  gam(n)   photoionization rate per atom in n [1/s] (incl. nonth for n=1)
!  alp(n)   spontaneous radiative recombination coefficient [cm^3/s]
!  alps(n)  stimulated recombination coefficient [cm^3/s]
!  cup(n)   collisional ionization coefficient [cm^3/s]
!  c3b(n)   three-body recombination coefficient per ion [cm^6/s]
!  qup(l,u), qdn(u,l)  collisional excitation/de-excitation [cm^3/s]
!  jl(k)    mean intensity incident on line k [erg/s/cm^2/Hz/sr]
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,trad
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng),jb(nh_nlin),nonth
  logical,intent(in) :: ljb(nh_nlin)
  real*8,intent(out) :: gam(nh_nlev),alp(nh_nlev),alps(nh_nlev),cup(nh_nlev),c3b(nh_nlev)
  real*8,intent(out) :: qup(nh_nlev,nh_nlev),qdn(nh_nlev,nh_nlev),jl(nh_nlin)
  integer :: n,k,ig,igred
  integer,parameter :: nsub = 200
  real*8 :: s,sst,sig,nu,dlt,phifac,numin,jred,xr
!
  do n=1,nh_nlev
     alp(n) = nlteh_alpha(n,temp)
     call nlteh_cion(n,temp,cup(n),c3b(n))
  enddo
  call nlteh_qcoll(temp,qup,qdn)
!-- (recombination to levels above the model atom is not included: it would
!-- need its reverse process to keep detailed balance; use enough levels,
!-- nh_nlev >= 30)
!-- photoionization and stimulated recombination from the group intensities;
!-- below the reddest group J is extrapolated with the shape of B_nu(trad)
  igred = minloc(nug,1)
  numin = nug(igred) - 0.5d0*dnug(igred)
  jred = jg(igred)
  do n=1,nh_nlev
     s = 0d0; sst = 0d0
     phifac = dble(n*n)*sahac*temp**(-1.5d0)   !Phi_n e^{-u_n}
     do ig=1,ng
        if(nug(ig)+0.5d0*dnug(ig) <= nh_nu0(n)) cycle
        nu = max(nug(ig), nh_nu0(n))
        sig = nlteh_sigma(n,nu)
        dlt = min(dnug(ig), nug(ig)+0.5d0*dnug(ig)-nh_nu0(n))
        s = s + jg(ig)*sig/(hh*nu)*dlt
        sst = sst + jg(ig)*sig/(hh*nu)*dlt*exp(-hh*(nu-nh_nu0(n))/(kb*temp))
     enddo
     if(nh_nu0(n) < numin .and. jred > 0d0) then
!-- threshold below the grid: integrate nu0..numin on a logarithmic sub-grid
        do ig=1,nsub
           nu = nh_nu0(n)*(numin/nh_nu0(n))**((ig-0.5d0)/nsub)
           dlt = nh_nu0(n)*((numin/nh_nu0(n))**(dble(ig)/nsub) - (numin/nh_nu0(n))**(dble(ig-1)/nsub))
           xr = jred*bplanck(nu,trad)/max(bplanck(nug(igred),trad),1d-300)
           sig = nlteh_sigma(n,nu)
           s = s + xr*sig/(hh*nu)*dlt
           sst = sst + xr*sig/(hh*nu)*dlt*exp(-hh*(nu-nh_nu0(n))/(kb*temp))
        enddo
     endif
     gam(n) = 4d0*pi*s
     alps(n) = 4d0*pi*phifac*sst
  enddo
  gam(1) = gam(1) + nonth
!-- line intensities
  do k=1,nh_nlin
     if(ljb(k)) then
        jl(k) = jb(k)
     else
        ig = minloc(abs(nug-nh_nu(k)),1)
        jl(k) = jg(ig)
     endif
  enddo
end subroutine nlteh_cell_rates
!
!
pure function nlteh_kq(k,nh,x) result(kq)
!-- Sobolev line strength K [1/s] of line k: tau(mu) = K/|Q(mu)|
  integer,intent(in) :: k
  real*8,intent(in) :: nh,x(0:nh_nlev)
  real*8 :: kq
  kq = fsob*nh_f(k)*(cc/nh_nu(k))*nh*(x(nh_ll(k)) - nh_g(nh_ll(k))/nh_g(nh_lu(k))*x(nh_lu(k)))
end function nlteh_kq
!
!
subroutine nlteh_solve_cell(temp,trad,nh,neoth,qb,qv,ng,nug,dnug,jg,jb,ljb,nonth,x,ne,beta,ierr)
!-----------------------------------------------------------------------
! Statistical equilibrium of H in one cell.
!  temp   electron temperature [K]
!  trad   radiation temperature [K]: shape of J beyond the red end of the group grid
!  nh     H number density [cm^-3]
!  neoth  electrons from other elements [cm^-3]
!  qb,qv  dv/dr and v/r at the cell centre [1/s]
!  ng,nug,dnug  group centre frequencies and widths [Hz]
!  jg     Monte Carlo mean intensity per group [erg/s/cm^2/Hz/sr]
!  jb     resonance mean intensity per line (used where ljb), else jg of the line's group
!  nonth  non-thermal ionization rate per H atom from the ground state [1/s]
!  x      in: initial guess (fractions of nh, 0 = H II); out: solution
!  ne     out: electron density
!  beta   out: Sobolev escape probability per line
!-----------------------------------------------------------------------
  real*8,intent(in) :: temp,trad,nh,neoth,qb,qv
  integer,intent(in) :: ng
  real*8,intent(in) :: nug(ng),dnug(ng),jg(ng),jb(nh_nlin),nonth
  logical,intent(in) :: ljb(nh_nlin)
  real*8,intent(inout) :: x(0:nh_nlev)
  real*8,intent(out) :: ne,beta(nh_nlin)
  integer,intent(out) :: ierr
!
  integer :: n,l,u,k,it,info,nn
  real*8 :: gam(nh_nlev),alp(nh_nlev),alps(nh_nlev),cup(nh_nlev),c3b(nh_nlev)
  real*8 :: qup(nh_nlev,nh_nlev),qdn(nh_nlev,nh_nlev)
  real*8 :: jl(nh_nlin)
  real*8 :: m(0:nh_nlev,0:nh_nlev),xn(0:nh_nlev)
  integer :: ipiv(0:nh_nlev)
  real*8 :: s,ne_old
!
  ierr = 0
  nn = nh_nlev
  call nlteh_cell_rates(temp,trad,ng,nug,dnug,jg,jb,ljb,nonth, &
       gam,alp,alps,cup,c3b,qup,qdn,jl)
!
!-- iterate escape probabilities and electron density
  ne = neoth + x(0)*nh
  do it=1,200
     do k=1,nh_nlin
        beta(k) = nlteh_beta(nlteh_kq(k,nh,x),qb,qv)
     enddo
!-- rate matrix: m(i,j) = rate j -> i, diagonal = -total out
     m = 0d0
     do k=1,nh_nlin
        l = nh_ll(k); u = nh_lu(k)
        m(u,l) = m(u,l) + beta(k)*nh_blu(k)*jl(k) + ne*qup(l,u)
        m(l,u) = m(l,u) + beta(k)*(nh_a(k) + nh_bul(k)*jl(k)) + ne*qdn(u,l)
     enddo
     m(1,2) = m(1,2) + 0.25d0*a2s
     do n=1,nn
        m(0,n) = m(0,n) + gam(n) + ne*cup(n)
        m(n,0) = m(n,0) + ne*(alp(n) + alps(n)) + ne*ne*c3b(n)
     enddo
     do n=0,nn
        m(n,n) = -(sum(m(:,n)) - m(n,n))
     enddo
!-- particle conservation replaces the H II equation
     m(0,:) = 1d0
     xn = 0d0
     xn(0) = 1d0
     call dgesv(nn+1,1,m,nn+1,ipiv,xn,nn+1,info)
     if(info/=0) then
        ierr = 1
        return
     endif
     xn = max(xn, 1d-300)
     xn = xn/sum(xn)
!-- converged?
     s = maxval(abs(xn-x)/(xn+1d-30), mask=xn>1d-20)
     ne_old = ne
     if(it<=3) then
        x = xn
     else
        x = sqrt(x*xn)       !geometric damping keeps fractions positive
        x = x/sum(x)
     endif
     ne = neoth + x(0)*nh
     if(s < 1d-6 .and. abs(ne-ne_old) <= 1d-8*ne) exit
  enddo
  if(it>200) ierr = 2
end subroutine nlteh_solve_cell
!
end module nltehmod
! vim: fdm=marker
