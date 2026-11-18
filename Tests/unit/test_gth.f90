!-- Unit test for nx_gth_solve (nltexmod): stationary populations of a rate
!-- matrix by the GTH state reduction against LU (dgesv), on (1) a random
!-- dense matrix of moderate range, where both agree to round-off, and (2)
!-- an ionization ladder of 27 stages with populations spanning ~1e-250,
!-- where the exact ratios are known: detailed balance between neighbours,
!-- p(k+1)/p(k) = up(k)/down(k).  GTH must reproduce every one to 1e-10
!-- relative and stay positive; LU is only reported.
program test_gth
  use nltexmod, only:nx_gth_solve
  implicit none
  integer,parameter :: n1 = 40, n2 = 27
  real*8 :: m(n1,n1),p(n1),b(n1),u
  real*8 :: m2(n2,n2),p2(n2),b2(n2),up(n2),dn(n2),ex(n2),err
  integer :: i,j,ipiv(n1),info,nneg
  logical :: ok
  ok = .true.
  call random_seed(put=[(12345+i,i=1,64)])
  do j=1,n1
     do i=1,n1
        call random_number(u)
        m(i,j) = 10d0**(6d0*u - 3d0)
     enddo
     m(j,j) = 0d0
  enddo
  call nx_gth_solve(n1,m,7,p)
  do j=1,n1
     m(j,j) = -sum(m(:,j))
  enddo
  m(7,:) = 1d0
  b = 0d0
  b(7) = 1d0
  call dgesv(n1,1,m,n1,ipiv,b,n1,info)
  err = maxval(abs(p-b)/b)
  write(*,'(a,es10.2)') ' dense 40 states: largest relative difference GTH - LU',err
  if(err > 1d-10 .or. any(p <= 0d0)) ok = .false.
!-- ladder: up rates 1e-9 .. , down rates 1e2: each stage 1e-11 .. of the last
  m2 = 0d0
  do i=1,n2-1
     up(i) = 10d0**(-7d0 - 0.025d0*i)
     dn(i) = 10d0**(2d0 + 0.025d0*i)
     m2(i+1,i) = up(i)
     m2(i,i+1) = dn(i)
  enddo
  ex(1) = 1d0
  do i=1,n2-1
     ex(i+1) = ex(i)*up(i)/dn(i)
  enddo
  ex = ex/sum(ex)
  call nx_gth_solve(n2,m2,1,p2)
  err = maxval(abs(p2-ex)/ex)
  write(*,'(a,es10.2,a,es10.2)') ' ladder of 27 stages down to',minval(ex),': largest relative error GTH',err
  if(err > 1d-10 .or. any(p2 <= 0d0)) ok = .false.
  do j=1,n2
     m2(j,j) = -sum(m2(:,j))
  enddo
  m2(1,:) = 1d0
  b2 = 0d0
  b2(1) = 1d0
  call dgesv(n2,1,m2,n2,ipiv,b2,n2,info)
  nneg = count(b2 < 0d0)
  write(*,'(a,es10.2,a,i3)') '   (LU: largest relative error',maxval(abs(b2-ex)/ex),', negative',nneg,')'
  if(ok) then
     write(*,*) 'test_gth: PASS'
  else
     write(*,*) 'test_gth: FAIL'
     stop 1
  endif
end program test_gth
