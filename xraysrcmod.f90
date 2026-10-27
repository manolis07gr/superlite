!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module xraysrcmod
!-----------------------------------------------------------------------
! X-ray emission of shocks embedded in a wind, as CMFGEN's INC_XRAYS
! (comp_opac, rd_xray_spec): up to two components k, with filling factor
! f_k, shock temperature T_k and onset velocity V_k, emit per unit volume
!   eta_nu = sum_k exp(-V_k/v) f_k^2 Lambda_nu(T_k) n_e,X N_atom
! with Lambda_nu the cooling function per n_e n_H of a plasma table in
! CMFGEN's RS_XRAY_FLUXES format (bins equally spaced in energy, APEC
! spectra in atomic/misc/gal_xray_hr.dat; linear in log T), n_e,X =
! n_H + 2 n_He + 6 n_other (H and He fully ionized, 6 electrons from each
! other atom) and N_atom the density of nuclei.  The densities are those
! of the smooth wind (volume means; CMFGEN's X_SM_WIND), or with
! xs_smooth = .false. those of the clumps times the filling factor.  The
! X-ray energy is a source beyond the injected energy (mechanical).
!-----------------------------------------------------------------------
  implicit none
  private
  logical,public :: xs_on = .false.
  integer :: nbin = 0
  real*8 :: emin = 0d0, de = 0d0           !lower edge of bin 1, bin width [erg]
  real*8,allocatable :: cumk(:,:)           !(nbin,2) cumulative spectrum of each component
  real*8 :: ltot(2) = 0d0                   !4 pi Int eta dnu per n_e,X N_atom f_k^2 [erg cm^3/s]
  real*8,allocatable,public :: xs_pow(:)    !(ncell) X-ray power of each cell [erg/s]
  real*8,allocatable :: xs_w1(:)            !(ncell) share of component 1
  real*8,allocatable,public :: xs_px(:)     !(ncell) share of the X-rays in the cell's emission
  real*8,public :: xs_q = 0d0               !least share of the packets drawn as X-rays (weighted)
  public :: xs_read,xs_setup,xs_sample,xs_group_emis
!
contains
!
!
subroutine xs_read(fname,tk)
!-- the plasma table at the two shock temperatures tk [K] (0: component off)
  character(*),intent(in) :: fname
  real*8,intent(in) :: tk(2)
  integer :: iu,ios,nt,ned,it,k,j(2),i
  real*8 :: bmin,bsiz,ltmin,dlt,w(2),lt
  real*8,allocatable :: row(:),sp(:,:)
  character(256) :: line
  real*8,parameter :: ev = 1.60217653d-12
!
  open(newunit=iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) stop 'xs_read: cannot open the X-ray table (in_xray_table)'
  nbin = ikey('[N_BINS]')
  bmin = rkey('[BIN_MIN]')
  bsiz = rkey('[BIN_SIZE]')
  nt = ikey('[N_TEMP]')
  ltmin = rkey('[LOG_T_MIN]')
  dlt = rkey('[DEL_LOG_T]')
  ned = ikey('[N_ED]')
  if(ned<1) stop 'xs_read: N_ED'
!-- bracketing table temperatures of each component
  j = 0
  w = 0d0
  do k=1,2
     if(tk(k)<=0d0) cycle
     lt = log10(tk(k))
     j(k) = int((lt - ltmin)/dlt) + 1
     if(j(k)<1 .or. j(k)>nt-1) stop 'xs_read: shock temperature outside the table'
     w(k) = (lt - (ltmin + (j(k)-1)*dlt))/dlt
  enddo
  allocate(row(nbin),sp(nbin,2))
  sp = 0d0
!-- spectra of the first electron density (the table is assumed density
!-- independent, as in CMFGEN)
  do it=1,nt
     do
        read(iu,'(a)',iostat=ios) line
        if(ios/=0) stop 'xs_read: table truncated'
        if(index(line,'Temperature(K)')>0) exit
     enddo
     read(iu,*,iostat=ios) row
     if(ios/=0) stop 'xs_read: bad spectrum'
     do k=1,2
        if(j(k)==it) sp(:,k) = sp(:,k) + (1d0-w(k))*row
        if(j(k)+1==it) sp(:,k) = sp(:,k) + w(k)*row
     enddo
  enddo
  close(iu)
!-- 1e-23 erg cm^3/s per bin (per n_e n_H): cumulative, total
  emin = (bmin - 0.5d0*bsiz)*ev
  de = bsiz*ev
  allocate(cumk(nbin,2))
  do k=1,2
     sp(:,k) = 1d-23*max(sp(:,k),0d0)
     cumk(1,k) = sp(1,k)
     do i=2,nbin
        cumk(i,k) = cumk(i-1,k) + sp(i,k)
     enddo
     ltot(k) = cumk(nbin,k)
     if(ltot(k)>0d0) cumk(:,k) = cumk(:,k)/ltot(k)
  enddo
  xs_on = any(ltot>0d0)
contains
  integer function ikey(key)
    character(*),intent(in) :: key
    ikey = nint(rkey(key))
  end function ikey
  real*8 function rkey(key)
    character(*),intent(in) :: key
    rewind(iu)
    do
       read(iu,'(a)',iostat=ios) line
       if(ios/=0) stop 'xs_read: keyword missing'
       if(index(line,key)>0) exit
    enddo
    read(line,*) rkey
  end function rkey
end subroutine xs_read
!
!
subroutine xs_setup(ncell,vel,vol,fcl,nat,xh,xhe,ff,vsh,smooth)
!-- X-ray power of each cell: vel [cm/s] at the cell centre, vol [cm^3],
!-- fcl filling factor, nat nuclei per cell, xh, xhe number fractions of H
!-- and He; ff, vsh [cm/s] of the two components
  integer,intent(in) :: ncell
  real*8,intent(in) :: vel(ncell),vol(ncell),fcl(ncell),nat(ncell),xh(ncell),xhe(ncell)
  real*8,intent(in) :: ff(2),vsh(2)
  logical,intent(in) :: smooth
  integer :: l
  real*8 :: na,nex,p(2)
  if(allocated(xs_pow)) deallocate(xs_pow,xs_w1,xs_px)
  allocate(xs_pow(ncell),xs_w1(ncell),xs_px(ncell))
  xs_px = 0d0
  xs_pow = 0d0
  xs_w1 = 0d0
  do l=1,ncell
     if(nat(l)<=0d0 .or. vol(l)<=0d0 .or. vel(l)<=0d0) cycle
     na = nat(l)/vol(l)
     if(.not.smooth) na = na/fcl(l)
     nex = na*(xh(l) + 2d0*xhe(l) + 6d0*max(0d0,1d0-xh(l)-xhe(l)))
     p = ff**2*ltot*exp(-vsh/vel(l))*nex*na
     if(.not.smooth) p = p*fcl(l)
     xs_pow(l) = sum(p)*vol(l)
     if(xs_pow(l)>0d0) xs_w1(l) = p(1)/sum(p)
  enddo
end subroutine xs_setup
!
!
subroutine xs_sample(l,r1,r2,r3,wl)
!-- wavelength [cm] of an X-ray photon of cell l from three uniform
!-- random numbers
  integer,intent(in) :: l
  real*8,intent(in) :: r1,r2,r3
  real*8,intent(out) :: wl
  integer :: k,a,b,m
  real*8 :: e
  k = 2
  if(r1<xs_w1(l)) k = 1
  a = 0
  b = nbin
  do while(b-a>1)
     m = (a+b)/2
     if(cumk(m,k) >= r2) then
        b = m
     else
        a = m
     endif
  enddo
  e = emin + (b - 1 + r3)*de
  wl = 6.6260693d-27*2.99792458d10/e
end subroutine xs_sample
!
!
subroutine xs_group_emis(ncell,ng,wlg,vol,em)
!-- the X-ray emission of each cell in each group, Int eta_nu dnu over the
!-- group per steradian [erg/s/cm^3/sr]; wlg(ng+1) the group edges [cm],
!-- increasing.  (For the source functions of the optically thick groups,
!-- where the X-ray photons are absorbed near where they are emitted.)
  integer,intent(in) :: ncell,ng
  real*8,intent(in) :: wlg(ng+1),vol(ncell)
  real*8,intent(out) :: em(ng,ncell)
  real*8,parameter :: hc = 6.6260693d-27*2.99792458d10, pi = 3.1415926535897932d0
  integer :: l,ig,k
  real*8 :: f(ng,2)
  do k=1,2
     do ig=1,ng
        f(ig,k) = cumat(hc/wlg(ig),k) - cumat(hc/wlg(ig+1),k)
     enddo
  enddo
  em = 0d0
  do l=1,ncell
     if(xs_pow(l)<=0d0 .or. vol(l)<=0d0) cycle
     em(:,l) = xs_pow(l)/(4d0*pi*vol(l))*(xs_w1(l)*f(:,1) + (1d0-xs_w1(l))*f(:,2))
  enddo
contains
  real*8 function cumat(e,k)
!-- share of the emission of component k below the photon energy e
    real*8,intent(in) :: e
    integer,intent(in) :: k
    real*8 :: u,c0
    integer :: b
    u = (e - emin)/de
    if(u <= 0d0) then
       cumat = 0d0
    elseif(u >= dble(nbin)) then
       cumat = 1d0
    else
       b = int(u)
       c0 = 0d0
       if(b>0) c0 = cumk(b,k)
       cumat = c0 + (u - b)*(cumk(b+1,k) - c0)
    endif
  end function cumat
end subroutine xs_group_emis
!
end module xraysrcmod
! vim: fdm=marker
