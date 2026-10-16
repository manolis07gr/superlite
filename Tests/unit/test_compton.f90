!-- Unit test for comptonmod: one thermal Thomson scattering of a
!-- monochromatic photon.  Checks the analytic moments
!--   <dnu/nu>    = (4kT - h nu)/(m_e c^2)   (Kompaneets mean shift)
!--   <(dnu/nu)^2> = 2 kT/(m_e c^2)           (isotropic phase function)
!-- and that the scattered directions are isotropic in the comoving frame.
program test_compton
  use randommod
  use comptonmod
  implicit none
  real*8,parameter :: kb = 1.3806505d-16, me = 9.1093826d-28, c = 2.99792458d10, h = 6.6260693d-27
  integer,parameter :: n = 4000000
  type(rnd_t) :: st
  real*8 :: temps(3), wls(3), t, wl, mu_in, mu_out, ratio, d, s1, s2, sm1, sm2, th, xm, mean_exp, var_exp
  integer :: i, k, nfail
  temps = [1d4, 1d6, 1d4]
  wls   = [6.5628d-5, 6.5628d-5, 1d-8]          ! H-alpha, H-alpha, 1 A X-ray
  nfail = 0
  call rnd_init(1, 7)
  st = rnd_states(1)
  do k = 1, 3
     t = temps(k); wl = wls(k)
     s1 = 0; s2 = 0; sm1 = 0; sm2 = 0
     do i = 1, n
        mu_in = 2d0*(dble(mod(i,1000))+0.5d0)/1000d0 - 1d0
        call thermal_thomson(mu_in, wl, t, .false., st, mu_out, ratio)
        d = ratio - 1d0
        s1 = s1 + d; s2 = s2 + d*d
        sm1 = sm1 + mu_out; sm2 = sm2 + mu_out**2
     enddo
     th = kb*t/(me*c*c); xm = h*c/wl/(me*c*c)
!-- exact isotropic Compton recoil averages, y = x(1-cos) uniform on [0,2x]:
!--   <1/(1+y) - 1> = ln(1+2x)/(2x) - 1,  <(1/(1+y) - 1)^2> = (Y - 2 ln(1+Y) + Y/(1+Y))/Y, Y = 2x
!-- plus the thermal terms 4kT/mc^2 and 2kT/mc^2 (cross terms are O(kT x/mc^2))
     mean_exp = 4d0*th + log(1d0+2d0*xm)/(2d0*xm) - 1d0
     var_exp = 2d0*th + (2d0*xm - 2d0*log(1d0+2d0*xm) + 2d0*xm/(1d0+2d0*xm))/(2d0*xm)
     write(*,'(a,es9.2,a,es9.2,a)') 'T =', t, ' K, lambda =', wl*1d8, ' A'
     write(*,'(a,es12.4,a,es12.4,a,es9.2)') '  <dnu/nu>     ', s1/n, '  expected', mean_exp, &
          '  +-', sqrt(s2/n)/sqrt(dble(n))
     write(*,'(a,es12.4,a,es12.4)') '  <(dnu/nu)^2> ', s2/n, '  expected', var_exp
     write(*,'(a,f9.5,a,f9.5)') '  <mu_out>     ', sm1/n, '   <mu_out^2>', sm2/n
     if(abs(s1/n - mean_exp) > 5d0*sqrt(s2/n)/sqrt(dble(n)) + 0.02d0*abs(mean_exp)) nfail = nfail + 1
     if(abs(s2/n/var_exp - 1d0) > 0.01d0) nfail = nfail + 1
     if(abs(sm1/n) > 2d-3 .or. abs(sm2/n - 1d0/3d0) > 2d-3) nfail = nfail + 1
  enddo
  if(nfail > 0) then
     write(*,'(a,i3)') 'test_compton: FAILED checks:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_compton: PASS'
end program test_compton
