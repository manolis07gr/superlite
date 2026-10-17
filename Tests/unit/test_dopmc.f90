!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!-----------------------------------------------------------------------
! Static Doppler-profile lines of slow cells (in_nltex_dopmc):
!  - opacity of a line at its centre and one Doppler width off,
!    K/(sqrt(pi) v_D) exp(-x^2), none beyond nxt_dopx widths;
!  - the packet's own line left out, and the lines of H and He (they keep
!    their Sobolev resonances, nxt_dopz), which the resonance search of the
!    cell finds instead;
!  - the absorbing line chosen by its share of the opacity;
!  - the incident-intensity tally, over a window outside the profile on
!    the side the comoving frequency drifts in from (blue for Q > 0, red
!    for Q < 0): a packet drifting through it in short flights adds what
!    one Sobolev crossing adds, e c/(nu |Q|); nothing from inside the
!    profile (attenuated by the line) or for the packet's own line.
!-----------------------------------------------------------------------
program test_dopmc
  use nlteatommod
  use nltextrmod
  implicit none
  real*8,parameter :: c = 2.99792458d10, spi = 1.7724538509055159d0
  real*8 :: wl1,wl2,wl3,vd,kq,kap,expect,q,e,dx,x,s,cnt(3)
  integer :: k,i,n,nfail,kpick,jlo,jhi
!
  nfail = 0
!-- two lines of Fe in one cell, 30 km/s apart, and a line of H two
!-- Doppler widths to the blue of the first, which keeps its Sobolev
!-- resonance in slow cells (nxt_dopz 3): it must not change anything
!-- below but the resonance search
  wl1 = 5000d-8
  wl2 = wl1*(1d0 + 30d5/c)
  vd = 20d5
  wl3 = wl1*(1d0 - 2d0*vd/c)
  kq = 1d3
  nxa_nlin = 3
  nxa_nel = 2
  allocate(nxa_el(2))
  nxa_el(1)%z = 26
  nxa_el(2)%z = 1
  allocate(nxa_lel(3),nxa_lnu(3),nxa_lwl(3))
  nxa_lel = [1,1,2]
  nxa_lwl = [wl1,wl2,wl3]
  nxa_lnu = c/nxa_lwl
  nxt_nact = 3
  allocate(nxt_act(3),nxt_wlact(3),nxt_kq(3,1),nxt_jbsum(3,1),nxt_jbcnt(3,1))
  nxt_act = [3,1,2]
  nxt_wlact = [wl3,wl1,wl2]
  nxt_kq(:,1) = [kq,2d0*kq,5d0*kq]
  allocate(nxt_slow(1),nxt_vdl(2,1),nxt_vdlmax(1),nxt_vdlmin(1))
  nxt_slow = .true.
  nxt_vdl = vd
  nxt_vdlmax = vd
  nxt_vdlmin = vd
  nxt_dopz = 3
  call nxt_act_lookup
  if(nxt_lstat(1) .and. nxt_lstat(2) .and. .not.nxt_lstat(3) .and. nxt_anystat) then
     write(6,'(1x,a,t53,a)') 'static profiles for Z >= 3 only','ok'
  else
     write(6,'(1x,a,t53,a)') 'static profiles for Z >= 3 only','FAILED'
     nfail = nfail + 1
  endif
!-- the resonance search of a slow cell finds the H line only
  call nxt_neighbors(wl1*(1d0-vd/c),0,0,jlo,jhi)
  call nxt_sobolev_neighbors(jlo,jhi)
  call ints('between H and Fe: H below, none above',jlo,jhi,1,0)
  call nxt_neighbors(wl2*(1d0+vd/c),0,0,jlo,jhi)
  call nxt_sobolev_neighbors(jlo,jhi)
  call ints('above all lines: H below, none above',jlo,jhi,1,0)
  call nxt_neighbors(wl3*(1d0-vd/c),0,0,jlo,jhi)
  call nxt_sobolev_neighbors(jlo,jhi)
  call ints('below all lines: none below, H above',jlo,jhi,0,1)
!
!-- opacity at the centre of line 1, own line 2 left out
  call nxt_static_line(1,wl1,2,-1d0,kap,kpick)
  expect = kq/(spi*vd)
  call check('centre of line 1',kap,expect)
!-- one Doppler width to the red of line 1, line 2 left out
  call nxt_static_line(1,wl1*(1d0+vd/c),2,-1d0,kap,kpick)
  call check('one width off',kap,expect*exp(-1d0))
!-- one Doppler width above line 2, the highest (both lines count once)
  call nxt_static_line(1,wl2*(1d0+vd/c),0,-1d0,kap,kpick)
  x = (wl2*(1d0+vd/c)/wl1 - 1d0)*c/vd
  call check('above all lines',kap,(2d0*kq*exp(-1d0) + kq*exp(-x*x))/(spi*vd))
!-- beyond nxt_dopx widths of both lines
  call nxt_static_line(1,wl2*(1d0+5d0*vd/c),0,-1d0,kap,kpick)
  call check('beyond the profiles',kap+1d-300,1d-300)
!-- both lines at the midpoint (x = 0.75 for each)
  call nxt_static_line(1,wl1*(1d0+15d5/c),0,-1d0,kap,kpick)
  x = (wl1*(1d0+15d5/c)/wl1 - 1d0)*c/vd
  expect = kq*exp(-x*x)/(spi*vd)
  x = (wl1*(1d0+15d5/c)/wl2 - 1d0)*c/vd
  expect = expect + 2d0*kq*exp(-x*x)/(spi*vd)
  call check('midpoint, both lines',kap,expect)
!-- the absorbing line by its share: line 2 holds 2/3 of the opacity
  cnt = 0d0
  n = 30000
  do i=1,n
     call nxt_static_line(1,wl1*(1d0+15d5/c),0,(i-0.5d0)/n,kap,kpick)
     cnt(kpick) = cnt(kpick) + 1d0
  enddo
  call check('share of line 2',cnt(2)/n,2d0/3d0)
!
!-- tally: drifts in flights of 0.01 Doppler widths with gradient q (x in
!-- Doppler widths from line 1, positive to the red); sums against
!-- e c/(nu |q|)
  q = 1d-8
  e = 1d0
  dx = 0.01d0
  s = vd*dx/q
!-- expanding (q > 0): from 8 widths to the blue of line 1 to its centre,
!-- through the blue windows of both lines
  call drift(-8d0,0d0,0,q)
  call check('blue window of line 1 against Sobolev',nxt_jbsum(1,1),e*c/(nxa_lnu(1)*q))
  call check('blue window of line 2 against Sobolev',nxt_jbsum(2,1),e*c/(nxa_lnu(2)*q))
  call check('no static tally for the H line',nxt_jbsum(3,1)+1d-300,1d-300)
!-- the same with line 2 the packet's own
  call drift(-8d0,0d0,2,q)
  call check('line 1 with line 2 its own',nxt_jbsum(1,1),e*c/(nxa_lnu(1)*q))
  call check('own line not tallied',nxt_jbsum(2,1)+1d-300,1d-300)
!-- inside the profile of line 1 (attenuated there): nothing
  call drift(-3.9d0,3.9d0,2,q)
  call check('inside the profile not tallied',nxt_jbsum(1,1)+1d-300,1d-300)
!-- decelerating (q < 0): from line 1 to 8 widths to its red, through its
!-- red window
  call drift(0d0,8d0,2,-q)
  call check('red window for a negative gradient',nxt_jbsum(1,1),e*c/(nxa_lnu(1)*q))
!-- the blue window does not count for a negative gradient
  call drift(-8d0,0d0,2,-q)
  call check('blue window, negative gradient',nxt_jbsum(1,1)+1d-300,1d-300)
!
  if(nfail == 0) then
     write(6,*) 'test_dopmc: PASS'
  else
     write(6,*) 'test_dopmc: FAIL',nfail
     stop 1
  endif
!
contains
!
  subroutine drift(xa,xb,kown,qq)
!-- flights at the midpoints of dx steps from xa to xb, fresh tallies
    real*8,intent(in) :: xa,xb,qq
    integer,intent(in) :: kown
    real*8 :: xx
    nxt_jbsum = 0d0
    nxt_jbcnt = 0d0
    xx = xa + 0.5d0*dx
    do while(xx < xb)
       call nxt_static_tally(1,wl1*(1d0+xx*vd/c),e*s,kown,qq,nxt_jbsum,nxt_jbcnt)
       xx = xx + dx
    enddo
  end subroutine drift
!
  subroutine ints(what,ga,gb,wa,wb)
    character(*),intent(in) :: what
    integer,intent(in) :: ga,gb,wa,wb
    if(ga==wa .and. gb==wb) then
       write(6,'(1x,a,t53,a)') what,'ok'
    else
       write(6,'(1x,a,t53,a,4i4)') what,'FAILED',ga,gb,wa,wb
       nfail = nfail + 1
    endif
  end subroutine ints
!
  subroutine check(what,got,want)
    character(*),intent(in) :: what
    real*8,intent(in) :: got,want
    real*8 :: dev
    dev = abs(got/want - 1d0)
    write(6,'(1x,a,t44,2es13.5,a,es9.2)') what,got,want,'  rel. dev.',dev
    if(dev > 1d-2) nfail = nfail + 1
  end subroutine check
!
end program test_dopmc
