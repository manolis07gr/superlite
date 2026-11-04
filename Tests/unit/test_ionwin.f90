!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
program test_ionwin
!-----------------------------------------------------------------------
! Ion window of the statistical equilibrium (nltexmod nx_win_tol): with H
! and an iron atom of many ions (fe1.atom ... in the model-atom directory,
! the first argument; skipped without), the populations of cells in
! photoionization and collisional equilibrium solved over the whole
! element and over the window agree (states above 1e-8 of the element
! within 1e-4), and the windowed solution takes less time.
!-----------------------------------------------------------------------
  use nlteatommod
  use nltexmod
  use nltextrmod
  use nltehmod, only:bplanck
  implicit none
  real*8,parameter :: km = 1d5, cl = 2.99792458d10
  integer,parameter :: ng = 4000, mf = 40
  real*8 :: wlg(ng+1),nug(ng),dnug(ng),jg(ng)
  real*8,allocatable :: x(:),xf(:),x0(:),beta(:),jb(:),nti(:),ntx(:),nel(:)
  logical,allocatable :: ljb(:)
  character(256) :: adir,files(mf)
  character(256) :: f
  integer :: i,k,nf,ierr,ic,iw,nfail,o,ns,nfe,ifull(2),nit(2)
  real*8 :: tolw
  character(32) :: ev
  logical :: lex
  real*8 :: t,trad,w,nh,ne,dev,tm(2),fion(30)
  real :: c0,c1
  real*8,parameter :: tc(3) = [1d4,2.5d4,3d5], tr(3) = [3d4,6d4,1d4], wc(3) = [1d-2,1d-3,1d-8]
  real*8,parameter :: nhc(3) = [1d10,1d9,1d9]
  nfail = 0
  tolw = 1d-12
  call get_environment_variable('SL_IWIN',ev,status=i)
  if(i==0 .and. len_trim(ev)>0) read(ev,*) tolw
  adir = ''
  if(command_argument_count() >= 1) call get_command_argument(1,adir)
  nf = 0
  if(len_trim(adir) > 0) then
     files(1) = trim(adir)//'/h1.atom'
     nf = 1
     do i=1,30
        write(f,'(a,a,i0,a)') trim(adir),'/fe',i,'.atom'
        inquire(file=f,exist=lex)
        if(.not.lex) exit
        nf = nf + 1
        files(nf) = f
     enddo
  endif
  nfe = nf - 1
  if(nfe < 6) then
     write(*,'(a)') ' (no iron atom of at least six ions in the model-atom directory: test skipped)'
     write(*,'(a)') ' test_ionwin: skipped'
     stop
  endif
  do i=1,ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i=1,ng
     nug(i) = 0.5d0*(cl/wlg(i) + cl/wlg(i+1))
     dnug(i) = cl/wlg(i) - cl/wlg(i+1)
  enddo
  call nxt_init(1,1,ng,wlg,nf,files(:nf),0)
  write(*,'(a,i3,a,i5,a)') ' H and Fe I -',nfe,' ions:',nxa_el(2)%ns,' states of Fe'
  allocate(x(nxa_nst),xf(nxa_nst),x0(nxa_nst),beta(nxa_nlin),jb(nxa_nlin),ljb(nxa_nlin))
  allocate(nti(nxa_nst),ntx(nxa_nst),nel(nxa_nel))
  nti = 0d0
  ntx = 0d0
  ljb = .true.
  o = nxa_el(2)%xoff
  ns = nxa_el(2)%ns
  do ic=1,3
     t = tc(ic)
     trad = tr(ic)
     w = wc(ic)
     nh = nhc(ic)
     do k=1,ng
        jg(k) = w*bplanck(nug(k),trad)
     enddo
     jb = [(w*bplanck(nxa_lnu(k),trad),k=1,nxa_nlin)]
     nel = [nh,3d-5*nh]
     call nx_lte(t,nh,x0)
     do iw=1,2
        nx_win_tol = merge(0d0,tolw,iw==1)
        nx_win_ns = 0
        x = x0
        call cpu_time(c0)
        call nx_solve_cell(t,trad,nel,0d0,1d0/1.7d6,1d0/1.7d6,ng,nug,dnug,jg,jb,ljb,nti,ntx,x,ne,beta, &
             ierr,maxit=1000)
        call cpu_time(c1)
        tm(iw) = c1 - c0
        ifull(iw) = ierr
        nit(iw) = nx_nits
        if(iw==1) xf = x
     enddo
     nx_win_tol = 0d0
     dev = 0d0
     do k=o+1,o+ns
        if(xf(k) > 1d-8) dev = max(dev,abs(x(k)/xf(k) - 1d0))
     enddo
!-- the ion fractions of the full solution
     fion = 0d0
     do i=1,nfe
        k = nxa_el(2)%ion(i)
        fion(i) = sum(xf(o+nxa_ion(k)%soff+1:o+nxa_ion(k)%soff+nxa_ion(k)%nsl))
     enddo
     write(*,'(a,f8.0,a,i3,a,es9.2,a,2f8.3,a,f6.3,a,2i2,a,2i5)') ' T =',t,' K: dominant Fe ion',maxloc(fion(:nfe),1), &
          '  max |x_win/x_full - 1| =',dev,'  seconds full, window',tm, &
          '  window states/all',real(nx_win_ns(1))/max(nx_win_ns(2),1),'  ierr',ifull,'  iterations',nit
     if(dev > 1d-4 .or. any(ifull/=0)) nfail = nfail + 1
  enddo
  if(nfail > 0) then
     write(*,'(a,i3)') ' test_ionwin: FAILED checks:',nfail
     stop 1
  endif
  write(*,'(a)') ' test_ionwin: all checks passed'
end program test_ionwin
