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
!  - Coulomb losses to the thermal electrons, n_e (2 pi e^4/E)
!    ln(4E/zeta), zeta = 7.40e-11 (n_e/T_e)^1/2 eV (Xu & McCray 1991), as
!    a flux in energy (upwind).
! Electrons that fall below the grid, and the Coulomb losses, heat.
! Cross-sections: Lotz (1967) for the subshells of an ion, 4.5e-14 q
! ln(E/I)/(E I) cm^2 (10-15 % above the measured ones for H I), or the
! fits of Arnaud & Rothenflug (1985) where known; the effective Gaunt
! factor form of van Regemorter
! (1962) for allowed excitations, (8 pi/sqrt 3) pi a0^2 (R^2/(E eps)) f g,
! g = max(0.2, (sqrt 3/2 pi) ln(E/eps)).
!-----------------------------------------------------------------------
  implicit none
  private
  public :: nt_grid, nt_lotz, nt_arfit, nt_vreg, nt_solve
  public :: nt_ne, nt_de, nt_e
!
  integer :: nt_ne = 0
  real*8 :: nt_de = 1d0
  real*8,allocatable :: nt_e(:)
!
  real*8,parameter :: a0 = 5.29177210903d-9, ryd = 13.605693122994d0
  real*8,parameter :: pi = 3.1415926535897932d0
!-- 2 pi e^4 [eV^2 cm^2]
  real*8,parameter :: e4 = 2d0*pi*(4.80320471d-10)**4/(1.602176634d-12)**2
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
pure subroutine nt_vreg(eps,f,sig)
!-- excitation cross-section of an allowed transition of energy eps [eV] and
!-- absorption oscillator strength f (van Regemorter 1962)
  real*8,intent(in) :: eps,f
  real*8,intent(out) :: sig(nt_ne)
  integer :: i
  real*8 :: g
  sig = 0d0
  do i=1,nt_ne
     if(nt_e(i) <= eps) cycle
     g = max(0.2d0,sqrt(3d0)/(2d0*pi)*log(nt_e(i)/eps))
     sig(i) = 8d0*pi/sqrt(3d0)*pi*a0**2*ryd**2/(nt_e(i)*eps)*f*g
  enddo
end subroutine nt_vreg
!
!
subroutine nt_solve(nion,ithr,ijop,ion_n,ion_sig,nexc,xeps,exc_n,exc_sig,ne,te,e0,fheat,nion_ev,nexc_ev)
!-----------------------------------------------------------------------
! One electron injected at e0 (<= the top of the grid): the share of its
! energy that heats (fheat; the rest goes into ionization and excitation)
! and the ionizations and excitations of each process per eV of injected
! energy (nion_ev, nexc_ev).  Processes: nion ionizations (threshold ithr
! [eV], Opal width ijop [eV], density ion_n of the ion [cm^-3], cross-
! section ion_sig on the grid [cm^2]) and nexc excitations (energy xeps,
! density exc_n of the lower level, cross-section exc_sig); the thermal
! electrons ne [cm^-3] at te [K].  Only ratios of the densities matter.
!-----------------------------------------------------------------------
  integer,intent(in) :: nion,nexc
  real*8,intent(in) :: ithr(nion),ijop(nion),ion_n(nion),ion_sig(nt_ne,nion)
  real*8,intent(in) :: xeps(nexc),exc_n(nexc),exc_sig(nt_ne,nexc)
  real*8,intent(in) :: ne,te,e0
  real*8,intent(out) :: fheat,nion_ev(nion),nexc_ev(nexc)
  real*8 :: y(nt_ne),arr(nt_ne),stot(nt_ne),coul(nt_ne),zeta,cfl,num,emax,an,ea,eb,pa,pb,de
  integer :: i,k,p,i0,ka,kb
  real*8 :: fion,fexc
!
  de = nt_de
  i0 = min(max(int(e0/de + 0.5d0),1),nt_ne)
!-- collision rate per unit path at each energy, and the Coulomb loss
  stot = 0d0
  do p=1,nion
     stot = stot + ion_n(p)*ion_sig(:,p)
  enddo
  do p=1,nexc
     stot = stot + exc_n(p)*exc_sig(:,p)
  enddo
  zeta = 7.40d-11*sqrt(max(ne,1d-30)/max(te,1d0))
  do i=1,nt_ne
     coul(i) = ne*e4/nt_e(i)*log(max(4d0*nt_e(i)/zeta,1.0001d0))
  enddo
!
  y = 0d0
  arr = 0d0
  arr(i0) = 1d0
  nion_ev = 0d0
  nexc_ev = 0d0
  do i=i0,1,-1
     cfl = 0d0
     if(i < i0) cfl = coul(i+1)*y(i+1)
!-- y de (stot + coul/de) = arrivals + the Coulomb flux from above
     y(i) = (arr(i) + cfl)/(de*stot(i) + coul(i))
     if(i == 1) exit
!-- excitations: to E - eps
     do p=1,nexc
        num = exc_n(p)*exc_sig(i,p)*y(i)*de
        if(num <= 0d0) cycle
        nexc_ev(p) = nexc_ev(p) + num
        call push(nt_e(i)-xeps(p),num)
     enddo
!-- ionizations: the primary to E - I - e, the secondary to e
     do p=1,nion
        num = ion_n(p)*ion_sig(i,p)*y(i)*de
        if(num <= 0d0 .or. nt_e(i) <= ithr(p)) cycle
        nion_ev(p) = nion_ev(p) + num
        emax = 0.5d0*(nt_e(i) - ithr(p))
        an = atan(emax/ijop(p))
!-- the secondary energy in the bins [k-1/2, k+1/2] de (clipped to [0, emax]);
!-- each piece pushed at its mean energy (the bin centre, or the clipped
!-- interval's centre)
        ka = 1
        kb = max(int(emax/de + 0.5d0),1)
        do k=ka,kb
           ea = max((k-0.5d0)*de,0d0)
           if(k == ka) ea = 0d0
           eb = min((k+0.5d0)*de,emax)
           if(eb <= ea) cycle
           pa = (atan(eb/ijop(p)) - atan(ea/ijop(p)))/an
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
end module nonthermmod
! vim: fdm=marker
