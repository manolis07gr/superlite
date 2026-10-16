!-- Unit test for dopplermod: exact comoving group crossings against a
!-- brute-force scan of the comoving wavelength along the path.
!-- Random cells with accelerating, decelerating and non-monotonic
!-- (a and b of opposite sign) linear velocity fields.
program test_doppler
  use dopplermod
  implicit none
  real*8,parameter :: c = 2.99792458d10
  integer,parameter :: ncase = 20000, nscan = 200000
  integer :: k, i, idir, idir_bf, nfail, nnone, nred, nblue, nnonmono
  real*8 :: xl, xr, vl, vr, r0, mu0, wl, wlo, whi, smax, s, s_bf, g, lam, a, b, ds, u(8), p2, r
  real*8 :: worst
  call random_seed()
  nfail = 0; nnone = 0; nred = 0; nblue = 0; nnonmono = 0; worst = 0d0
  do k = 1, ncase
     call random_number(u)
     xl = 1d14*(1d0 + 50d0*u(1))
     xr = xl*(1d0 + 0.3d0*u(2))
!-- edge velocities between 100 km/s and 20000 km/s, any ordering
     vl = 1d7 + 2d9*u(3)
     vr = 1d7 + 2d9*u(4)
     if(u(5) < 0.2d0) vr = vl*xr/xl                    ! homologous
     r0 = xl + (xr - xl)*u(6)
     mu0 = 2d0*u(7) - 1d0
     b = (vr - vl)/(xr - xl); a = vl - b*xl
     if(a*b < 0d0) nnonmono = nnonmono + 1
!-- path length to the cell boundary
     if(mu0 >= -sqrt(max(0d0,1d0-(xl/r0)**2))) then
        smax = abs(sqrt(xr**2 - (1d0-mu0**2)*r0**2) - mu0*r0)
     else
        smax = abs(sqrt(xl**2 - (1d0-mu0**2)*r0**2) + mu0*r0)
     endif
!-- comoving group of width 50-600 km/s around the current comoving wavelength
     g = dop_g(0d0, r0, mu0, a, b)
     call random_number(u)
     lam = 6.5d-5
     wl = lam*(1d0 - g/c)
     ds = (5d6 + 5.5d7*u(1))/c
     wlo = lam*(1d0 - ds*u(2)); whi = wlo + lam*ds      ! wlo <= lam < whi
     call doppler_crossing(r0, mu0, xl, xr, vl, vr, wl, wlo, whi, smax, c, s, idir)
!-- brute force
     idir_bf = 0; s_bf = huge(s_bf)
     do i = 1, nscan
        s_bf = smax*dble(i)/nscan
        g = dop_g(s_bf, r0, mu0, a, b)
        lam = wl/(1d0 - g/c)
        if(lam >= whi) then; idir_bf = 1; exit; endif
        if(lam < wlo) then; idir_bf = -1; exit; endif
     enddo
     if(idir_bf == 0) s_bf = huge(s_bf)
     if(idir == 0) nnone = nnone + 1
     if(idir == 1) nred = nred + 1
     if(idir == -1) nblue = nblue + 1
     if(idir /= idir_bf) then
        nfail = nfail + 1
        if(nfail <= 5) write(*,'(a,2i3,4es12.4)') 'direction mismatch', idir, idir_bf, s, s_bf, smax, mu0
     elseif(idir /= 0) then
        worst = max(worst, abs(s - s_bf)/smax)
        if(abs(s - s_bf) > 2d0*smax/nscan) then
           nfail = nfail + 1
           if(nfail <= 5) write(*,'(a,i3,3es14.6)') 'distance mismatch', idir, s, s_bf, smax
        endif
!-- at the crossing the comoving wavelength sits on the boundary
        g = dop_g(s, r0, mu0, a, b)
        lam = wl/(1d0 - g/c)
        if(idir == 1) r = lam/whi - 1d0
        if(idir == -1) r = lam/wlo - 1d0
        if(abs(r) > 1d-11) then
           nfail = nfail + 1
           if(nfail <= 5) write(*,'(a,es12.4)') 'boundary miss', r
        endif
     endif
  enddo
  write(*,'(a,i6,a,i6,a,i6,a,i6,a,i6)') 'cases', ncase, '  red', nred, '  blue', nblue, '  none', nnone, &
       '  non-monotonic cells', nnonmono
  write(*,'(a,es10.2)') 'worst |s - s_bruteforce|/smax', worst
  if(nfail > 0) then
     write(*,'(a,i6)') 'FAILED cases:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_doppler: PASS'
end program test_doppler
