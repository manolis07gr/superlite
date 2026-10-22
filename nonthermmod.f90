!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nonthermmod
!-----------------------------------------------------------------------
! Degradation of fast electrons in a partly ionized gas (Spencer & Fano
! 1954; Kozma & Fransson 1992): the share of the energy of the electrons
! made by gamma-rays, X-rays or positrons that heats the thermal
! electrons, and the ionizations and excitations it makes per unit
! energy.  The degradation function y(E) (electrons per unit energy times
! path length: n_p sigma_p(E) y(E) dE collisions of process p in dE) on a
! uniform grid up to the injection energy is solved from the top down; a
! finished bin pushes its collisions to the bins the electrons land in
! (linear split between neighbours), so that the number of electrons is
! conserved exactly:
!  - excitation p (threshold eps): to E - eps;
!  - ionization p (threshold I): the primary to E - I - e and the
!    secondary to e, e in [0, (E-I)/2] with the distribution of Opal,
!    Peterson & Beaty (1971), 1/(1 + (e/J)^2), integrated over each bin;
!  - Coulomb losses to the thermal electrons, n_e (2 pi e^4/E) ln(2E/xi),
!    xi = hbar omega_p the plasma energy, and below 14 eV
!    ln(m v^3/(gamma e^2 omega_p)) (Kozma & Fransson 1992; as CMFGEN), as a
!    flux in energy (upwind).
! Electrons that fall below the grid, and the Coulomb losses, heat.
! Cross-sections: Lotz (1967) for the subshells of an ion, 4.5e-14 q
! ln(E/I)/(E I) cm^2 (10-15 % above the measured ones for H I), or the
! subshell fits of Arnaud & Rothenflug (1985); excitation (nt_xsec) by the
! effective Gaunt factor form of van Regemorter (1962) for allowed
! transitions, (8 pi^2/sqrt 3) a0^2 (R^2/(E eps)) f g(E/eps), with g of
! neutral atoms (0 at the threshold) or of ions (0.2 near it), or with a
! collision strength constant in energy, pi a0^2 (R/E) Omega/g_l, for the
! forbidden ones (as CMFGEN's non-thermal solver).
!-----------------------------------------------------------------------
  implicit none
  private
  public :: nt_grid, nt_lotz, nt_arfit, nt_xsec, nt_gbar, nt_solve
  public :: nt_ne, nt_de, nt_e
!-- kinds of excitation: allowed in a neutral atom, allowed in an ion,
!-- constant collision strength
  integer,parameter,public :: nt_kneu = 1, nt_kion = 2, nt_kom = 3
!
  integer :: nt_ne = 0
  real*8 :: nt_de = 1d0
  real*8,allocatable :: nt_e(:)
!
  real*8,parameter :: a0 = 5.29177210903d-9, ryd = 13.605693122994d0
  real*8,parameter :: pi = 3.1415926535897932d0
  real*8,parameter :: qe = 4.80320471d-10, evg = 1.602176634d-12, mel = 9.1093837d-28
  real*8,parameter :: hbar = 1.054571817d-27, eulg = 0.5772156649d0
!-- 2 pi e^4 [eV^2 cm^2]
  real*8,parameter :: e4 = 2d0*pi*qe**4/evg**2
!-- van Regemorter: (8 pi^2/sqrt 3) a0^2 R^2 [cm^2 eV^2]
  real*8,parameter :: cvr = 8d0*pi*pi/sqrt(3d0)*a0*a0*ryd*ryd
!
contains
!
!
subroutine nt_grid(emax,de)
!-- the energy grid: bin centres de, 2 de, ..., emax [eV]
  real*8,intent(in) :: emax,de
  integer :: i
  nt_de = de
  nt_ne = max(int(emax/de + 0.5d0),2)
  if(allocated(nt_e)) deallocate(nt_e)
  allocate(nt_e(nt_ne))
  do i=1,nt_ne
     nt_e(i) = i*de
  enddo
end subroutine nt_grid
!
!
pure subroutine nt_lotz(ni,thr,q,sig)
!-- ionization cross-section of an ion with subshells of thresholds thr [eV]
!-- and occupations q (Lotz 1967, a = 4.5e-14 cm^2 eV^2, b = c = 0)
  integer,intent(in) :: ni
  real*8,intent(in) :: thr(ni),q(ni)
  real*8,intent(out) :: sig(nt_ne)
  integer :: i,j
  sig = 0d0
  do j=1,ni
     if(thr(j) <= 0d0 .or. q(j) <= 0d0) cycle
     do i=1,nt_ne
        if(nt_e(i) <= thr(j)) cycle
        sig(i) = sig(i) + 4.5d-14*q(j)*log(nt_e(i)/thr(j))/(nt_e(i)*thr(j))
     enddo
  enddo
end subroutine nt_lotz
!
!
pure subroutine nt_arfit(thr,a,b,c,d,sig)
!-- direct ionization cross-section of Arnaud & Rothenflug (1985),
!-- 1e-14/(u I^2) [A(1-1/u) + B(1-1/u)^2 + C ln u + D ln(u)/u] cm^2, u = E/I
  real*8,intent(in) :: thr,a,b,c,d
  real*8,intent(out) :: sig(nt_ne)
  integer :: i
  real*8 :: u
  sig = 0d0
  do i=1,nt_ne
     u = nt_e(i)/thr
     if(u <= 1d0) cycle
     sig(i) = max(1d-14/(u*thr*thr)*(a*(1d0-1d0/u) + b*(1d0-1d0/u)**2 + c*log(u) + d*log(u)/u),0d0)
  enddo
end subroutine nt_arfit
!
!
pure function nt_gbar(kind,u) result(g)
!-- effective Gaunt factor of an allowed excitation at E = u eps (van
!-- Regemorter 1962), as a function of x = sqrt(u - 1): x (1 + x) 0.074
!-- near the threshold of a neutral atom, 0.2 near that of an ion (up to
!-- x = 1.2212, where the two meet), a quadratic in x up to x = 6, and the
!-- Bethe limit (sqrt 3/2 pi) ln u above
  integer,intent(in) :: kind
  real*8,intent(in) :: u
  real*8 :: g,x2,x
  g = 0d0
  x2 = u - 1d0
  if(x2 <= 0d0) return
  if(kind==nt_kion .and. x2 <= 1.4913889d0) then
     g = 0.2d0
  elseif(x2 <= 0.64d0) then
     x = sqrt(x2)
     g = 0.074d0*x*(1d0 + x)
  elseif(x2 <= 36d0) then
     x = sqrt(x2)
     g = -0.0745397d0 + x*(0.232715d0 - 0.00647558d0*x)
  else
     g = 0.105d0 + log(x2)/3.6276d0
  endif
  g = max(g,0d0)
end function nt_gbar
!
!
pure function nt_xsec(kind,eps,par,e) result(s)
!-- excitation cross-section [cm^2] at electron energy e [eV] of a
!-- transition of energy eps [eV]: allowed (kind nt_kneu, nt_kion; par =
!-- the absorption oscillator strength f) or with the constant collision
!-- strength (nt_kom; par = Omega/g_l)
  integer,intent(in) :: kind
  real*8,intent(in) :: eps,par,e
  real*8 :: s
  s = 0d0
  if(e <= eps .or. eps <= 0d0) return
  if(kind==nt_kom) then
     s = pi*a0*a0*ryd/e*par
  else
     s = cvr*par/(e*eps)*nt_gbar(kind,e/eps)
  endif
end function nt_xsec
!
!
subroutine nt_solve(nion,ithr,ijop,ion_n,ion_sig,nexc,xeps,xkind,xpar,exc_n,ne,e0,fheat,nion_ev,nexc_ev)
!-----------------------------------------------------------------------
! One electron injected at e0 (<= the top of the grid): the share of its
! energy that heats (fheat; the rest goes into ionization and excitation)
! and the ionizations and excitations of each process per eV of injected
! energy (nion_ev, nexc_ev).  Processes: nion ionizations (threshold ithr
! [eV], Opal width ijop [eV], density ion_n of the ion [cm^-3], cross-
! section ion_sig on the grid [cm^2]) and nexc excitations (energy xeps
! [eV], kind xkind and parameter xpar of nt_xsec, density exc_n of the
! lower level); the thermal electrons ne [cm^-3].  A process with zero
! density is left out.
!-----------------------------------------------------------------------
  integer,intent(in) :: nion,nexc
  real*8,intent(in) :: ithr(nion),ijop(nion),ion_n(nion),ion_sig(nt_ne,nion)
  real*8,intent(in) :: xeps(nexc),xpar(nexc),exc_n(nexc)
  integer,intent(in) :: xkind(nexc)
  real*8,intent(in) :: ne,e0
  real*8,intent(out) :: fheat,nion_ev(nion),nexc_ev(nexc)
  real*8 :: y(nt_ne),arr(nt_ne),stot(nt_ne),coul(nt_ne),cfl,num,emax,an,ea,eb,pa,de,st
  real*8 :: sx(nexc),da(nt_ne/2+1,nion)
  integer :: i,k,p,i0,kb
  real*8 :: fion,fexc
!
  de = nt_de
  i0 = min(max(int(e0/de + 0.5d0),1),nt_ne)
!-- ionization rate per unit path at each energy, the Coulomb loss, and
!-- the Opal integrals over the full bins of the secondary energy,
!-- atan((k+1/2) de/J) - atan((k-1/2) de/J) (from 0 for k = 1)
  stot = 0d0
  do p=1,nion
     if(ion_n(p) <= 0d0) cycle
     stot = stot + ion_n(p)*ion_sig(:,p)
     da(1,p) = atan(1.5d0*de/ijop(p))
     do k=2,size(da,1)
        da(k,p) = atan((k+0.5d0)*de/ijop(p)) - atan((k-0.5d0)*de/ijop(p))
     enddo
  enddo
  call nt_coulomb(ne,coul)
!
  y = 0d0
  arr = 0d0
  arr(i0) = 1d0
  nion_ev = 0d0
  nexc_ev = 0d0
  do i=i0,1,-1
     cfl = 0d0
     if(i < i0) cfl = coul(i+1)*y(i+1)
!-- the excitation rates at this energy
     st = stot(i)
     do p=1,nexc
        sx(p) = 0d0
        if(exc_n(p) <= 0d0 .or. nt_e(i) <= xeps(p)) cycle
        sx(p) = exc_n(p)*nt_xsec(xkind(p),xeps(p),xpar(p),nt_e(i))
        st = st + sx(p)
     enddo
!-- y de (stot + coul/de) = arrivals + the Coulomb flux from above
     y(i) = (arr(i) + cfl)/(de*st + coul(i))
     if(i == 1) exit
!-- excitations: to E - eps
     do p=1,nexc
        if(sx(p) <= 0d0) cycle
        num = sx(p)*y(i)*de
        nexc_ev(p) = nexc_ev(p) + num
        call push(nt_e(i)-xeps(p),num)
     enddo
!-- ionizations: the primary to E - I - e, the secondary to e
     do p=1,nion
        if(ion_n(p) <= 0d0 .or. nt_e(i) <= ithr(p)) cycle
        num = ion_n(p)*ion_sig(i,p)*y(i)*de
        if(num <= 0d0) cycle
        nion_ev(p) = nion_ev(p) + num
        emax = 0.5d0*(nt_e(i) - ithr(p))
        an = atan(emax/ijop(p))
!-- the secondary energy in the bins [k-1/2, k+1/2] de (from 0 in the first,
!-- clipped to emax in the last); each piece pushed at its mean energy
        kb = max(int(emax/de + 0.5d0),1)
        do k=1,kb
           ea = (k-0.5d0)*de
           if(k == 1) ea = 0d0
           eb = (k+0.5d0)*de
           if(eb < emax) then
              pa = da(k,p)/an
           else
              eb = emax
              if(eb <= ea) cycle
              pa = (an - atan(ea/ijop(p)))/an
           endif
           call push(0.5d0*(ea+eb),num*pa)
           call push(nt_e(i)-ithr(p)-0.5d0*(ea+eb),num*pa)
        enddo
     enddo
  enddo
!-- energy shares: ionization (the binding energy), excitation; the rest
!-- (Coulomb losses and electrons below the grid) heats
  fion = 0d0
  do p=1,nion
     fion = fion + ithr(p)*nion_ev(p)
  enddo
  fexc = 0d0
  do p=1,nexc
     fexc = fexc + xeps(p)*nexc_ev(p)
  enddo
  fheat = max(1d0 - (fion + fexc)/nt_e(i0),0d0)
  nion_ev = nion_ev/nt_e(i0)
  nexc_ev = nexc_ev/nt_e(i0)
contains
  subroutine push(et,nn)
!-- nn electrons at energy et: linear split onto the bin centres
    real*8,intent(in) :: et,nn
    real*8 :: xx,w
    integer :: kk
    if(et < nt_e(1)) then
       arr(1) = arr(1) + nn
       return
    endif
    xx = et/de
    kk = int(xx)
    w = xx - kk
    if(kk >= nt_ne) then
       arr(nt_ne) = arr(nt_ne) + nn
       return
    endif
    arr(kk) = arr(kk) + (1d0-w)*nn
    arr(kk+1) = arr(kk+1) + w*nn
  end subroutine push
end subroutine nt_solve
!
!
pure subroutine nt_coulomb(ne,coul)
!-- energy loss per unit path to the thermal electrons [eV/cm] on the grid:
!-- n_e (2 pi e^4/E) ln(2E/(hbar omega_p)) at E >= 14 eV, below that
!-- ln(m v^3/(gamma e^2 omega_p)), gamma = 0.5772 (Kozma & Fransson 1992)
  real*8,intent(in) :: ne
  real*8,intent(out) :: coul(nt_ne)
  real*8 :: wp,xi,v,e,lg
  integer :: i
  wp = sqrt(4d0*pi*max(ne,1d-30)*qe*qe/mel)
  xi = hbar*wp/evg
  do i=1,nt_ne
     e = nt_e(i)
     if(e >= 14d0) then
        lg = log(2d0*e/xi)
     else
        v = sqrt(2d0*e*evg/mel)
        lg = log(mel*v**3/(eulg*qe*qe*wp))
     endif
     coul(i) = ne*e4/e*max(lg,1d-3)
  enddo
end subroutine nt_coulomb
!
end module nonthermmod
! vim: fdm=marker
