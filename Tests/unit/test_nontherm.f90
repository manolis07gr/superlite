!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!-----------------------------------------------------------------------
! Degradation of fast electrons (nonthermmod):
!  - pure hydrogen (Arnaud & Rothenflug ionization, van Regemorter
!    excitations of the Lyman series), 3 keV electrons, against the shares of heating,
!    ionization and excitation of Shull & van Steenberg (1985, their fits
!    for H with He/H = 0.1) at fractional ionizations 1e-4 to 0.9, within
!    0.1 (they share the energy with He and use measured cross-sections:
!    the van Regemorter excitation of H here takes up to 0.08 more);
!  - the work per ion pair of neutral hydrogen, 36.3 eV (Dalgarno et al.
!    1999), and the heating share that goes to 1 as x_e -> 1;
!  - grid convergence: 1 eV and 0.5 eV bins;
!  - trace ions: the ionizations per particle of a route left out of the
!    degradation (zero density) equal those of the same route at 1e-6 of
!    the main target.
!-----------------------------------------------------------------------
program test_nontherm
  use nonthermmod
  implicit none
  integer,parameter :: nx = 5
  real*8,parameter :: lyf(5) = [0.4162d0, 0.0791d0, 0.0290d0, 0.0139d0, 0.0078d0]
  real*8 :: xev(7) = [1d-4, 1d-3, 1d-2, 3d-2, 0.1d0, 0.3d0, 0.9d0]
  real*8 :: ithr(1),ijop(1),ionn(1),exn(nx),eps(nx),fh,ni(1),nexc(nx),fi,fe,sh,si,se,dmax,w0
  real*8 :: fh1,fh2,lyfp(nx),whe
  integer :: kx(nx)
  real*8,allocatable :: sigi(:,:)
  integer :: i,k,nfail
!
  nfail = 0
  call setup(0.5d0)
  dmax = 0d0
  write(*,'(a)') '  x_e     heat   ion    exc   | Shull & van Steenberg 1985: heat   ion    exc'
  do i=1,7
     call run(xev(i),fh,fi,fe)
     sh = 0.9971d0*(1d0 - (1d0 - xev(i)**0.2663d0)**1.3163d0)
     si = 0.3908d0*(1d0 - xev(i)**0.4092d0)**1.7592d0
     se = 0.4766d0*(1d0 - xev(i)**0.2735d0)**1.5221d0
     write(*,'(es8.1,3f7.3,a,3f7.3)') xev(i),fh,fi,fe,'   |                             ',sh,si,se
     dmax = max(dmax,abs(fh-sh),abs(fi-si),abs(fe-se))
  enddo
!-- work per ion pair of neutral H
  call run(1d-6,fh,fi,fe)
  w0 = 13.6d0/fi
  write(*,'(a,f6.1,a)') '  work per ion pair of neutral hydrogen',w0,' eV (Dalgarno et al. 1999: 36.3 eV)'
!-- grid convergence at x_e = 0.01
  call run(1d-2,fh1,fi,fe)
  call setup(1d0)
  call run(1d-2,fh2,fi,fe)
  write(*,'(a,f7.4,a,f7.4)') '  heating share at x_e = 0.01 with 0.5 eV bins',fh1,', 1 eV bins',fh2
  write(*,'(a,f6.3)') '  largest difference from Shull & van Steenberg:',dmax
  if(dmax > 0.1d0 .or. abs(w0/36.3d0 - 1d0) > 0.1d0 .or. abs(fh1-fh2) > 0.01d0) nfail = nfail + 1
!-- neutral helium
  call he_run(whe)
  write(*,'(a,f6.1,a)') '  work per ion pair of neutral helium',whe,' eV (measured 41.3-42.3 eV; Dalgarno et al. 1999: 46 eV)'
  if(whe < 38d0 .or. whe > 50d0) nfail = nfail + 1
!-- trace ion
  call trace_run(fi,fe)
  write(*,'(a,es11.4,a,es11.4)') '  trace He in H: ionizations per particle and eV, left out',fi,', at 1e-6',fe
  if(abs(fi/fe - 1d0) > 1d-3) nfail = nfail + 1
  if(nfail > 0) then
     write(*,'(a)') 'test_nontherm: FAILED'
     stop 1
  endif
  write(*,'(a)') 'test_nontherm: PASS'
contains
  subroutine setup(de)
    real*8,intent(in) :: de
    call nt_grid(3000d0,de)
    if(allocated(sigi)) deallocate(sigi)
    allocate(sigi(nt_ne,1))
    call nt_arfit(13.6d0,22.8d0,-12.0d0,1.9d0,-22.6d0,sigi(:,1))
    do k=1,nx
       eps(k) = 13.6d0*(1d0 - 1d0/(k+1)**2)
    enddo
    lyfp = lyf
    kx = nt_kneu
  end subroutine setup
  subroutine run(xe,fheat,fion,fexc)
    real*8,intent(in) :: xe
    real*8,intent(out) :: fheat,fion,fexc
    ithr = 13.6d0
    ijop = 8d0
    ionn = 1d9*(1d0 - xe)
    exn = 1d9*(1d0 - xe)
    call nt_solve(1,ithr,ijop,ionn,sigi,nx,eps,kx,lyfp,exn,1d9*xe,3000d0,fheat,ni,nexc)
    fion = 13.6d0*ni(1)
    fexc = sum(eps*nexc)
  end subroutine run
  subroutine he_run(w)
!-- He I (x_e = 1e-6): ionization by the fit of Arnaud & Rothenflug, the
!-- allowed 1s2 1S - 1snp 1P lines (n = 2-5) and the forbidden 2 3S, 2 1S
!-- and 2 3P with constant collision strengths
    real*8,intent(out) :: w
    integer,parameter :: nh = 7
    real*8 :: he(nh),hp(nh),hn(nh),hx(nh),sg(1),fhh
    integer :: hk(nh)
    he = [21.218d0, 23.087d0, 23.742d0, 24.046d0, 19.820d0, 20.616d0, 20.964d0]
    hp = [0.2762d0, 0.0734d0, 0.0302d0, 0.0153d0, 0.068d0, 0.035d0, 0.022d0]
    hk = [nt_kneu, nt_kneu, nt_kneu, nt_kneu, nt_kom, nt_kom, nt_kom]
    hn = 1d9
    call nt_grid(3000d0,0.5d0)
    deallocate(sigi)
    allocate(sigi(nt_ne,1))
    call nt_arfit(24.587d0,17.8d0,-11.0d0,7.0d0,-23.2d0,sigi(:,1))
    call nt_solve(1,[24.587d0],[0.6d0*24.587d0],[1d9],sigi,nh,he,hk,hp,hn,1d3,3000d0,fhh,sg,hx)
    w = 24.587d0/(24.587d0*sg(1))
  end subroutine he_run
  subroutine trace_run(ppout,ppin)
!-- H at x_e = 0.01 with He I as a second route, left out or at 1e-6 of H
    real*8,intent(out) :: ppout,ppin
    real*8 :: n2(2),pp(2),nn(2)
    real*8,allocatable :: s2(:,:)
    call setup(1d0)
    allocate(s2(nt_ne,2))
    s2(:,1) = sigi(:,1)
    call nt_arfit(24.587d0,17.8d0,-11.0d0,7.0d0,-23.2d0,s2(:,2))
    exn = 0.99d9
    n2 = [0.99d9, 0d0]
    call nt_solve(2,[13.6d0,24.587d0],[8d0,0.6d0*24.587d0],n2,s2,nx,eps,kx,lyfp,exn,1d7,3000d0,fh,nn,nexc,pp)
    ppout = pp(2)
    n2(2) = 0.99d3
    call nt_solve(2,[13.6d0,24.587d0],[8d0,0.6d0*24.587d0],n2,s2,nx,eps,kx,lyfp,exn,1d7,3000d0,fh,nn,nexc,pp)
    ppin = nn(2)/n2(2)
  end subroutine trace_run
end program test_nontherm
