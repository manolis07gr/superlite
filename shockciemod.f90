!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module shockciemod
!-----------------------------------------------------------------------
! Emission of the radiative shocks ('lumg' source, src_tsh > 0) as a
! cooling flow in collisional ionization equilibrium (CIE).  Isobaric
! cooling releases the same energy per unit temperature, so T is uniform
! between the gas temperature and the post-shock temperature; at T the
! photon wavelength follows the CIE spectrum of the cell's composition,
!   eps_lambda(T) = sum_el n_el eps_el,lambda(T)   (per n_e),
! from per-element emissivity tables (lines, free-free, free-bound,
! two-photon) made from the CHIANTI database with ChiantiPy
! (Tools/chianti_cie_spectra.py; in_shock_cie).  Between the table
! temperatures one of the two neighbours is drawn with the linear weight in
! log T; outside the table the caller keeps the bremsstrahlung shape.
!-----------------------------------------------------------------------
  implicit none
  private
  logical,public :: shc_on = .false.
  integer,public :: shc_nel = 0
  integer,allocatable,public :: shc_elz(:)      !(nel) atomic numbers of the table
  integer :: nt = 0, nw = 0
  real*8,allocatable :: logt(:)                 !(nt) log10 T [K]
  real*8,allocatable :: ledg(:)                 !(nw+1) log10 of the bin edges [A]
  real*8,allocatable :: em(:,:,:)               !(nw,nt,nel) emissivity per bin per n_e n_el
!-- shock cells: cumulative normalized spectrum at each table temperature
  integer,allocatable :: shc_cell(:)            !(ncell) index of the shock cell, 0: none
  real*8,allocatable :: cum(:,:,:)              !(nw,nt,nshock)
  public :: shc_read,shc_setup,shc_sample,shc_tmin,shc_tmax
!
contains
!
!
subroutine shc_read(fname)
!-- read the CIE emissivity table (all ranks)
  character(*),intent(in) :: fname
  integer :: iu,ios,i,j,k,z,n
  character(64) :: key
  real*8 :: lt
  real*8,allocatable :: edg(:)
  open(newunit=iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) then
     write(0,*) 'shc_read: cannot open ',trim(fname)
     stop 'shc_read: CIE table missing'
  endif
  do
     read(iu,'(a)',iostat=ios) key
     if(ios/=0) stop 'shc_read: no header'
     if(key(1:1)/='#') exit
  enddo
  backspace(iu)
  read(iu,*) key,shc_nel
  allocate(shc_elz(shc_nel))
  read(iu,*) key,shc_elz
  read(iu,*) key,nt
  allocate(logt(nt))
  read(iu,*) key,logt
  read(iu,*) key,nw
  allocate(edg(nw+1),ledg(nw+1),em(nw,nt,shc_nel))
  read(iu,*) key,edg
  ledg = log10(edg)
  read(iu,*) key
  em = 0d0
  do n=1,shc_nel*nt
     read(iu,*,iostat=ios) z,lt
     if(ios/=0) stop 'shc_read: table truncated'
     k = findloc(shc_elz,z,1)
     j = minloc(abs(logt-lt),1)
     if(k<1 .or. abs(logt(j)-lt)>1d-3) stop 'shc_read: bad record'
     read(iu,*) (em(i,j,k),i=1,nw)
  enddo
  close(iu)
  em = max(em,0d0)
  shc_on = .true.
end subroutine shc_read
!
!
subroutine shc_setup(ncell,tsh,xfr)
!-- cumulative spectra of the shock cells (tsh > 0) for the number of atoms
!-- xfr(k,l) of table element k in cell l (any normalization)
  integer,intent(in) :: ncell
  real*8,intent(in) :: tsh(ncell),xfr(shc_nel,ncell)
  integer :: l,n,j,k
  real*8 :: s
  if(allocated(shc_cell)) deallocate(shc_cell)
  if(allocated(cum)) deallocate(cum)
  allocate(shc_cell(ncell))
  shc_cell = 0
  n = 0
  do l=1,ncell
     if(tsh(l)>0d0) then
        n = n + 1
        shc_cell(l) = n
     endif
  enddo
  allocate(cum(nw,nt,max(n,1)))
  cum = 0d0
  do l=1,ncell
     n = shc_cell(l)
     if(n==0) cycle
     do j=1,nt
        do k=1,shc_nel
           cum(:,j,n) = cum(:,j,n) + xfr(k,l)*em(:,j,k)
        enddo
        do k=2,nw
           cum(k,j,n) = cum(k,j,n) + cum(k-1,j,n)
        enddo
        s = cum(nw,j,n)
        if(s > 0d0) then
           cum(:,j,n) = cum(:,j,n)/s
        else
!-- no emitting element in the table: flat in log lambda
           cum(:,j,n) = [(dble(k)/nw,k=1,nw)]
        endif
     enddo
  enddo
end subroutine shc_setup
!
!
pure real*8 function shc_tmin()
!-- lowest temperature of the table [K]
  shc_tmin = 10d0**logt(1)
end function shc_tmin
!
!
pure real*8 function shc_tmax()
!-- highest temperature of the table [K]
  shc_tmax = 10d0**logt(nt)
end function shc_tmax
!
!
subroutine shc_sample(l,t,r1,r2,r3,wl)
!-- comoving wavelength [cm] of a photon of shock cell l emitted by gas at
!-- temperature t (table range), from three uniform random numbers
  integer,intent(in) :: l
  real*8,intent(in) :: t,r1,r2,r3
  real*8,intent(out) :: wl
  integer :: n,j,a,b,m
  real*8 :: x
  n = shc_cell(l)
!-- table temperature: neighbour in log T with the linear weight
  x = (log10(t) - logt(1))/(logt(2) - logt(1))
  j = min(max(int(x) + 1,1),nt-1)
  x = x - (j-1)
  if(r1 < x) j = j + 1
  j = min(max(j,1),nt)
!-- bin: first with cum >= r2
  a = 0
  b = nw
  do while(b-a>1)
     m = (a+b)/2
     if(cum(m,j,n) >= r2) then
        b = m
     else
        a = m
     endif
  enddo
!-- uniform in log lambda within the bin
  wl = 1d-8*10d0**(ledg(b) + r3*(ledg(b+1) - ledg(b)))
end subroutine shc_sample
!
end module shockciemod
! vim: fdm=marker
