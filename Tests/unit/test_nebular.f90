!-- Unit test for ions_nebular_mc (modified nebular approximation with the
!-- Monte Carlo radiation temperature and dilution).  Run in a directory
!-- with data.ion and data.zeta (Data/).
!--  1. W = 1, T_e = T_R: the ion fractions must equal the LTE Saha solution
!--     (ions_solve_eos) for all elements, at several T and densities.
!--  2. Diluting the radiation field (W < 1) must lower the ionization of H
!--     and of the electron density monotonically.
!--  3. T_e >> T_R (delta factors of ~e^1000 for highly charged ions) must
!--     give finite fractions and electron density.
program test_nebular
  use ionsmod
  implicit none
  integer, parameter :: nel = 30
  real*8 :: natomfr(nel), nelec, nelec2, ndens, t, dmax, w, xh, xhold, neold
  real*8 :: frac(nel, nel+1)
  real*8, parameter :: ts(3) = [6d3, 1d4, 2d4], ns(2) = [1d9, 1d13], ws(4) = [1d0, 0.3d0, 0.1d0, 0.01d0]
  integer :: i, j, iz, ii, it, nfail
  nfail = 0
  call ions_read_data(nel)
!-- solar-like number fractions (H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe)
  natomfr = 0d0
  natomfr(1) = 0.92d0; natomfr(2) = 0.078d0; natomfr(6) = 3d-4; natomfr(7) = 8d-5; natomfr(8) = 5d-4
  natomfr(10) = 8d-5; natomfr(12) = 4d-5; natomfr(14) = 3d-5; natomfr(16) = 1.5d-5; natomfr(20) = 2d-6
  natomfr(26) = 3d-5
  natomfr = natomfr/sum(natomfr)
!-- 1. LTE limit
  do i = 1, 3
     do j = 1, 2
        t = ts(i); ndens = ns(j)
        nelec = 0.5d0
        call ions_solve_eos(natomfr, t, ndens, nelec, it)
        do iz = 1, nel
           do ii = 1, ion_el(iz)%ni
              frac(iz, ii) = ion_el(iz)%i(ii)%n
           enddo
        enddo
        nelec2 = 0.5d0
        call ions_nebular_mc(natomfr, t, t, 1d0, ndens, nelec2)
        dmax = 0d0
        do iz = 1, nel
           if(natomfr(iz) <= 0d0) cycle
           do ii = 1, ion_el(iz)%ni
              if(frac(iz, ii) > 1d-6) dmax = max(dmax, abs(ion_el(iz)%i(ii)%n/frac(iz, ii) - 1d0))
           enddo
        enddo
        write(6,'(a,f7.0,a,es8.1,a,2f9.5,a,es9.2)') 'T=', t, ' n=', ndens, '  n_e/n LTE, nebular W=1:', &
             nelec, nelec2, '   max rel. dev. of fractions >1e-6:', dmax
        if(abs(nelec2/nelec - 1d0) > 1d-6 .or. dmax > 1d-5) nfail = nfail + 1
     enddo
  enddo
!-- 2. dilution
  t = 1d4; ndens = 1d10
  xhold = 2d0; neold = 1d30
  do i = 1, 4
     w = ws(i)
     nelec = 0.5d0
     call ions_nebular_mc(natomfr, 0.9d0*t, t, w, ndens, nelec)
     xh = ion_el(1)%i(2)%n
     write(6,'(a,f6.3,a,es10.3,a,es10.3)') 'T_R=1e4, T_e=9e3, n=1e10, W=', w, '  x(H II)=', xh, '  n_e/n=', nelec
     if(xh > xhold .or. nelec > neold) nfail = nfail + 1
     xhold = xh; neold = nelec
  enddo
!-- 3. T_e >> T_R
  nelec = 0.5d0
  call ions_nebular_mc(natomfr, 3d4, 5d3, 1d-2, 1d12, nelec)
  dmax = 0d0
  do iz = 1, nel
     do ii = 1, ion_el(iz)%ni
        if(ion_el(iz)%i(ii)%n /= ion_el(iz)%i(ii)%n .or. ion_el(iz)%i(ii)%n > 1d0+1d-12) dmax = 1d0
     enddo
  enddo
  write(6,'(a,es10.3,a,es10.3)') 'T_e=3e4, T_R=5e3, W=0.01: n_e/n =', nelec, '  x(H II) =', ion_el(1)%i(2)%n
  if(dmax > 0d0 .or. nelec /= nelec .or. nelec <= 0d0) nfail = nfail + 1
  if(nfail == 0) then
     write(6,*) 'test_nebular: PASS'
  else
     write(6,*) 'test_nebular: FAIL', nfail
     stop 1
  endif
end program test_nebular
