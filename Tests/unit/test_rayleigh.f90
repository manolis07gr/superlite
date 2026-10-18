!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
program test_rayleigh
!-----------------------------------------------------------------------
! Rayleigh scattering by H I (rayleighmod) against an independent
! evaluation of the Kramers-Heisenberg sum (Python/scipy: the same
! oscillator strengths, the continuum term by adaptive quadrature), the
! low-frequency limit sigma = 1.2654 sigma_T x^4 (the static polarizability
! of H, 4.5 a_0^3), the cap at Lyman alpha and the fit above the edge.
!-----------------------------------------------------------------------
  use rayleighmod
  implicit none
  real*8,parameter :: c = 2.99792458d10, nu1 = c*109678.764d0
  real*8,parameter :: lam(10) = [5000d0,1300d0,1260d0,1240d0,1230d0,1225d0,1205d0,1200d0,1150d0,1050d0]
  real*8,parameter :: ref(10) = [1.54808d-3,10.6694d0,36.3032d0,115.533d0,324.849d0,756.038d0, &
       542.26d0,246.838d0,10.6351d0,0.488864d0]
  real*8 :: x,s,dmax,sg(1),nb(2)
  integer :: i,nfail
  nfail = 0
!-- against the independent sum (wavelengths in A, vacuum: x = 0.75 lambda(Lya)/lambda)
  dmax = 0d0
  do i=1,size(lam)
     x = 0.75d0/(0.75d0*109678.764d-8)/lam(i)
     s = ray_sigma(x*nu1)
     dmax = max(dmax,abs(s/ref(i)-1d0))
     write(6,'(a,f8.1,a,es12.5,a,es12.5)') ' lambda',lam(i),' A: sigma/sigma_T',s,'  reference',ref(i)
  enddo
  write(6,'(a,es9.2)') ' largest relative deviation from the independent sum:',dmax
  if(dmax > 2d-3) nfail = nfail + 1
!-- Rayleigh's law with the static polarizability of H
  x = 0.02d0
  s = ray_sigma(x*nu1)/x**4
  write(6,'(a,f8.5,a)') ' sigma/(sigma_T x^4) at x = 0.02:',s,' (1.2654)'
  if(abs(s/1.2654d0-1d0) > 2d-3) nfail = nfail + 1
!-- the cap within about 1500 km/s of Lyman alpha: 1e4 f(Lya)^2
  s = ray_sigma(nu1*0.75d0*(1d0+3d-4))
  write(6,'(a,f9.2,a)') ' at Lyman alpha + 90 km/s:',s,' (1e4 f^2 = 1732.2)'
  if(abs(s/1732.2d0-1d0) > 1d-3) nfail = nfail + 1
!-- above the Lyman limit (Gavrila's fit)
  s = ray_sigma(1.5d0*nu1)
  if(abs(s-(1d0+1.66d0/sqrt(1.5d0))) > 1d-12) nfail = nfail + 1
!-- a group average equals the point value where the cross-section is flat
  nb = [nu1*0.75d0*(1d0-2d-4),nu1*0.75d0*(1d0+2d-4)]
  call ray_groups(2,nb,sg)
  write(6,'(a,f9.2)') ' group mean within +-60 km/s of Lyman alpha:',sg(1)
  if(abs(sg(1)/1732.2d0-1d0) > 1d-3) nfail = nfail + 1
  if(nfail > 0) then
     write(6,'(a,i2,a)') ' test_rayleigh:',nfail,' checks FAILED'
     stop 1
  endif
  write(6,'(a)') ' test_rayleigh: all checks passed'
end program test_rayleigh
