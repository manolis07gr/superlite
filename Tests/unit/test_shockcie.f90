!-- Test of shockciemod: a synthetic CIE table with two elements (Z = 1
!-- emits in bin 5 only, Z = 8 in bin 10 at log T = 5 and in bin 20 at
!-- log T = 6) and a shock cell with n(O)/n(H) = 3: the photons of each
!-- line follow the composition and the temperature weights, and fall in
!-- their bins.  The shares per group of the cooling flow (shc_group_frac,
!-- the source function of the thick groups): at one T on groups equal to
!-- the bins the same shares as the photons; half a bin in log lambda half
!-- of it; a T range of the table against the integral over T; above the
!-- table the bremsstrahlung shape against its integral; all summing to 1.
program test_shockcie
  use shockciemod
  implicit none
  integer,parameter :: nw = 30, nt = 2, ndraw = 400000
  real*8 :: edges(nw+1), r1, r2, r3, wl, t, f5, f10, f20, xfr(2,3), tsh(3)
  integer :: i, k, n5, n10, n20, nbad, nfail, iu
  real*8 :: wlg(nw+1), fg(nw), wlh(nw+2), fh(nw+1), tt, fx, tlo, thi, e1, e2
  real*8,parameter :: hck = 6.6260693d-27*2.99792458d10/1.3806505d-16
  nfail = 0
!-- table: bins of 0.1 dex from 1 A
  do i = 1, nw+1
     edges(i) = 10d0**(0.1d0*(i-1))
  enddo
  open(newunit=iu, file='test_shockcie.dat', status='replace')
  write(iu,'(a)') '# synthetic CIE table'
  write(iu,'(a)') 'nel 2'
  write(iu,'(a)') 'elements 1 8'
  write(iu,'(a)') 'nt 2'
  write(iu,'(a)') 'logt 5.0 6.0'
  write(iu,'(a,i0)') 'nw ', nw
  write(iu,'(a,99es14.6)') 'edges ', edges
  write(iu,'(a)') 'ne 1e10'
  write(iu,'(a)') '1 5.0'
  write(iu,'(99es12.4)') [(merge(2d0, 0d0, i==5), i=1,nw)]
  write(iu,'(a)') '1 6.0'
  write(iu,'(99es12.4)') [(merge(2d0, 0d0, i==5), i=1,nw)]
  write(iu,'(a)') '8 5.0'
  write(iu,'(99es12.4)') [(merge(1d0, 0d0, i==10), i=1,nw)]
  write(iu,'(a)') '8 6.0'
  write(iu,'(99es12.4)') [(merge(1d0, 0d0, i==20), i=1,nw)]
  close(iu)
  call shc_read('test_shockcie.dat')
!-- cell 2 is the shock: n_H = 1, n_O = 3
  tsh = [0d0, 1d7, 0d0]
  xfr = 0d0
  xfr(:,2) = [1d0, 3d0]
  call shc_setup(3, tsh, xfr)
!-- T at log T = 5.25: weight 0.75 on log T = 5, 0.25 on 6
  t = 10d0**5.25d0
  n5 = 0; n10 = 0; n20 = 0; nbad = 0
  do k = 1, ndraw
     call random_number(r1); call random_number(r2); call random_number(r3)
     call shc_sample(2, t, r1, r2, r3, wl)
     wl = wl*1d8
!-- (bins with a relative tolerance of 1e-6: the table holds the edges to 7
!-- digits)
     if(inbin(wl, 5)) then
        n5 = n5 + 1
     elseif(inbin(wl, 10)) then
        n10 = n10 + 1
     elseif(inbin(wl, 20)) then
        n20 = n20 + 1
     else
        nbad = nbad + 1
     endif
  enddo
!-- H share at each T: 2*1/(2*1 + 1*3) = 0.4
  f5 = dble(n5)/ndraw
  f10 = dble(n10)/ndraw
  f20 = dble(n20)/ndraw
  write(*,'(a,3f8.4,a,i0)') '  shares of bins 5, 10, 20:', f5, f10, f20, '  outside: ', nbad
  write(*,'(a,3f8.4)') '  expected:                 ', 0.4d0, 0.6d0*0.75d0, 0.6d0*0.25d0
  if(abs(f5-0.4d0) > 0.005d0 .or. abs(f10-0.45d0) > 0.005d0 .or. abs(f20-0.15d0) > 0.005d0 &
       .or. nbad > 0) nfail = nfail + 1
  write(*,'(a,2es10.3)') '  table range [K]:', shc_tmin(), shc_tmax()
  if(abs(shc_tmin()/1d5 - 1d0) > 1d-10 .or. abs(shc_tmax()/1d6 - 1d0) > 1d-10) nfail = nfail + 1
!
!-- shares per group: groups equal to the table bins, one T (to 1e-5: the
!-- table holds the edges to 7 digits)
  wlg = edges*1d-8
  call shc_group_frac(2, t, t, nw, wlg, fg)
  write(*,'(a,4f9.5)') '  group shares 5, 10, 20, sum:', fg(5), fg(10), fg(20), sum(fg)
  if(abs(fg(5)-0.4d0) > 1d-5 .or. abs(fg(10)-0.45d0) > 1d-5 .or. abs(fg(20)-0.15d0) > 1d-5 &
       .or. abs(sum(fg)-1d0) > 1d-12) nfail = nfail + 1
!-- bin 10 split at its middle in log lambda
  wlh(1:10) = wlg(1:10)
  wlh(11) = sqrt(wlg(10)*wlg(11))
  wlh(12:nw+2) = wlg(11:nw+1)
  call shc_group_frac(2, t, t, nw+1, wlh, fh)
  write(*,'(a,2f9.5)') '  halves of bin 10:', fh(10), fh(11)
  if(abs(fh(10)-0.225d0) > 1d-5 .or. abs(fh(11)-0.225d0) > 1d-5) nfail = nfail + 1
!-- T uniform over 10^5 - 10^5.5 K: bin 10 holds 0.6 (1 - log10(T) + 5)
!-- averaged over T (the 64 midpoints of the routine against 10^5 points)
  tlo = 1d5
  thi = 10d0**5.5d0
  call shc_group_frac(2, tlo, thi, nw, wlg, fg)
  fx = 0d0
  do k = 1, 100000
     tt = tlo + (k-0.5d0)/100000*(thi-tlo)
     fx = fx + 0.6d0*(1d0 - (log10(tt) - 5d0))/100000
  enddo
  write(*,'(a,3f9.5)') '  T range, bin 10 against the integral, sum:', fg(10), fx, sum(fg)
  if(abs(fg(10)/fx-1d0) > 1d-4 .or. abs(sum(fg)-1d0) > 1d-12) nfail = nfail + 1
!-- above the table: bremsstrahlung at 10^7 K, the share of 3-4 A, the
!-- tails beyond the grid in the first and last group
  tt = 1d7
  call shc_group_frac(2, tt, tt, nw, wlg, fg)
  e1 = exp(-hck/(wlg(6)*tt)) - exp(-hck/(wlg(5)*tt))
  e2 = exp(-hck/(wlg(1)*tt)) + exp(-hck/(wlg(2)*tt)) - exp(-hck/(wlg(1)*tt))
  write(*,'(a,4es12.4)') '  bremsstrahlung, group 5 and 1 against the integrals:', fg(5), e1, fg(1), e2
  if(abs(fg(5)/e1-1d0) > 1d-12 .or. abs(fg(1)/e2-1d0) > 1d-12 .or. abs(sum(fg)-1d0) > 1d-12) nfail = nfail + 1
  if(nfail > 0) then
     write(*,'(a,i3)') 'test_shockcie: FAILED checks:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_shockcie: PASS'
contains
  logical function inbin(w, i)
    real*8, intent(in) :: w
    integer, intent(in) :: i
    inbin = w >= edges(i)*(1d0 - 1d-6) .and. w < edges(i+1)*(1d0 + 1d-6)
  end function inbin
end program test_shockcie
