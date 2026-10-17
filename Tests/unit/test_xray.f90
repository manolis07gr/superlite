!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!-----------------------------------------------------------------------
! X-ray source of wind shocks (xraysrcmod) on a synthetic plasma table in
! CMFGEN's RS_XRAY_FLUXES format (10 bins of 10 eV from 5 eV; a flat
! spectrum at 10^6 K, one rising with energy at 10^7 K):
!  - the power of each cell against sum_k f_k^2 Lambda_k exp(-V_k/v) n_e,X N_atom V
!    (the second component halfway in log T between the two spectra);
!  - the emission per group (xs_group_emis): the groups together hold the
!    whole power, 4 pi V Sum_g em = P; a group over half a table bin holds
!    half of it with each component's share; groups outside the table nothing.
!-----------------------------------------------------------------------
program test_xray
  use xraysrcmod
  implicit none
  real*8,parameter :: hc = 6.6260693d-27*2.99792458d10, ev = 1.60217653d-12, pi = 3.1415926535897932d0
  integer,parameter :: ng = 7
  real*8 :: vel(2),vol(2),fcl(2),nat(2),xh(2),xhe(2),ff(2),vsh(2)
  real*8 :: eg(ng+1),wlg(ng+1),em(ng,2),na,nex,p(2),ltot(2),w1,s1,s2
  integer :: iu,i,l,nfail
!
  nfail = 0
!-- the table
  open(newunit=iu,file='xray_test_table.dat',status='replace')
  write(iu,'(a)') '! synthetic X-ray table'
  write(iu,'(a)') '10            [N_BINS]'
  write(iu,'(a)') '10.0          [BIN_MIN]'
  write(iu,'(a)') '10.0          [BIN_SIZE]'
  write(iu,'(a)') '2             [N_TEMP]'
  write(iu,'(a)') '6.0           [LOG_T_MIN]'
  write(iu,'(a)') '1.0           [DEL_LOG_T]'
  write(iu,'(a)') '1             [N_ED]'
  write(iu,'(a)') '10.0          [LOG_ED_MIN]'
  write(iu,'(a)') '1.0           [DEL_LOG_ED]'
  write(iu,'(a)') 'Temperature(K) =  1.000E+06 Log(Ne)= 1.000E+01'
  write(iu,'(10es12.4)') (1d0,i=1,10)
  write(iu,'(a)') 'Temperature(K) =  1.000E+07 Log(Ne)= 1.000E+01'
  write(iu,'(10es12.4)') (dble(i),i=1,10)
  close(iu)
  call xs_read('xray_test_table.dat',[1d6,10d0**6.5d0])
!-- spectra per bin: flat (component 1), 0.5 + 0.5 b (component 2)
  ltot = 1d-23*[10d0, 0.5d0*10d0 + 0.5d0*55d0]
!
!-- two cells
  vel = [2d7,3d7]
  vol = [1d45,2d45]
  fcl = 1d0
  nat = [1d53,4d53]
  xh = 0.9d0
  xhe = 0.1d0
  ff = [0.1d0,0.2d0]
  vsh = [1d7,1d7]
  call xs_setup(2,vel,vol,fcl,nat,xh,xhe,ff,vsh,.true.)
  do l=1,2
     na = nat(l)/vol(l)
     nex = na*(xh(l) + 2d0*xhe(l))
     p = ff**2*ltot*exp(-vsh/vel(l))*nex*na
     call check('power of the cell',xs_pow(l),sum(p)*vol(l))
  enddo
!
!-- groups (energies decreasing, wavelengths increasing): above the table,
!-- 55-105, 50-55 (half of the bin 45-55 eV), 30-50, 5-30, below the table
  eg = [200d0,105d0,55d0,50d0,30d0,5d0,2d0,1d0]*ev
  wlg = hc/eg
  call xs_group_emis(2,ng,wlg,vol,em)
  w1 = ff(1)**2*ltot(1)/(ff(1)**2*ltot(1) + ff(2)**2*ltot(2))
  do l=1,2
     call check('all groups: the power',4d0*pi*vol(l)*sum(em(:,l)),xs_pow(l))
  enddo
!-- half of bin 5: 0.5/10 of the flat spectrum, 0.5*3/32.5 of the other
  s1 = 0.5d0/10d0
  s2 = 0.5d0*3d0/32.5d0
  call check('half a bin, both components',em(3,1)/sum(em(:,1)),w1*s1 + (1d0-w1)*s2)
!-- 55-105 eV: bins 6-10
  s1 = 5d0/10d0
  s2 = (0.5d0*5d0 + 0.5d0*40d0)/32.5d0
  call check('five bins',em(2,2)/sum(em(:,2)),w1*s1 + (1d0-w1)*s2)
  call check('above the table: nothing',em(1,1)+1d-300,1d-300)
  call check('below the table: nothing',em(6,1)+em(7,1)+1d-300,1d-300)
!
  open(newunit=iu,file='xray_test_table.dat',status='old')
  close(iu,status='delete')
  if(nfail == 0) then
     write(6,*) 'test_xray: PASS'
  else
     write(6,*) 'test_xray: FAIL',nfail
     stop 1
  endif
!
contains
!
  subroutine check(what,got,want)
    character(*),intent(in) :: what
    real*8,intent(in) :: got,want
    real*8 :: dev
    dev = abs(got/want - 1d0)
    write(6,'(1x,a,t40,2es13.5,a,es9.2)') what,got,want,'  rel. dev.',dev
    if(dev > 1d-10) nfail = nfail + 1
  end subroutine check
!
end program test_xray
