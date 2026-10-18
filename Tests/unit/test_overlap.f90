!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
program test_overlap
!-----------------------------------------------------------------------
! Overlapping resonance zones (nltexmod: nx_ovl_zone, nx_ovl_cell, the
! terms of nx_el_step):
!  1. one line: the escape probability (1 - exp(-tau))/tau
!  2. Lyman beta with O I 1025.76 (11.7 km/s to the red): the probability
!     that O I absorbs a Lyman-beta photon is the ratio of the depths, an
!     O I photon is absorbed by Lyman beta; esc + sum p = 1; detailed
!     balance tau_k (esc + p(j)) = A_k + tau_j p_j(k); against an
!     independent evaluation (Python, a five times finer grid)
!  3. a blend of two equal lines escapes as one line of the summed depth;
!     two thick lines three Doppler widths apart
!  4. with the H and O I model atoms (directory as the first argument;
!     skipped without): the LTE limit of the statistical equilibrium with
!     the overlap terms at high density in a Planck field, and in a dilute
!     field the Bowen pumping of O I 3d 3D by Lyman beta
!-----------------------------------------------------------------------
  use nlteatommod
  use nltexmod
  use nltextrmod
  use nltehmod, only:bplanck
  implicit none
  real*8,parameter :: km = 1d5, cl = 2.99792458d10
  integer,parameter :: ng = 6000
  real*8 :: esc,bet1,p(2),ax(2),e2,b2,p2(2),dev,lhs,rhs
  real*8 :: wlg(ng+1),nug(ng),dnug(ng),jg(ng)
  real*8,parameter :: taus(5) = [0.1d0,1d0,10d0,1d3,1d7]
  integer :: i,nfail,nargs
  character(256) :: adir
  nfail = 0
!-- 1. one line (Q the same in all directions: qb = qv)
  dev = 0d0
  do i=1,5
     call nx_ovl_zone(1,[0d0],[20d0*km],[taus(i)],1,1d0,1d0,esc,bet1,p(1:1))
     dev = max(dev,abs(esc/((1d0-exp(-taus(i)))/taus(i)) - 1d0),abs(esc+p(1)-1d0),abs(bet1/esc-1d0))
  enddo
  write(*,'(a,es10.3)') ' one line: max |escape/((1-exp(-tau))/tau) - 1|, sum rule:',dev
  if(dev > 1d-6) nfail = nfail + 1
!-- 2. Lyman beta (tau 5e6, 22.8 km/s) and O I 1025.76 (tau 650, 20.2 km/s)
  call nx_ovl_zone(2,[0d0,-11.7d0*km],[22.8d0*km,20.2d0*km],[5d6,650d0],1,1d0,1d0,esc,bet1,p,ax)
  call nx_ovl_zone(2,[0d0,-11.7d0*km],[22.8d0*km,20.2d0*km],[5d6,650d0],2,1d0,1d0,e2,b2,p2)
  write(*,'(a,es11.4,a,es11.4,a,f8.5,a,es10.3)') ' Lyb: escape',esc,'  absorbed by O I',p(2), &
       ' | O I: absorbed by Lyb',p2(1),'  escape',e2
  write(*,'(a)') '      (independent evaluation: 1.999e-07, 1.300e-04 | 0.99980, 4.93e-07)'
  if(abs(p(2)/1.300d-4 - 1d0) > 5d-3 .or. abs(p2(1) - 0.9998d0) > 2d-4 .or. abs(esc/1.999d-7 - 1d0) > 5d-3) &
       nfail = nfail + 1
  if(abs(esc + sum(p) - 1d0) > 1d-10 .or. abs(e2 + sum(p2) - 1d0) > 1d-10) nfail = nfail + 1
  lhs = 5d6*(esc + p(2))
  rhs = ax(1) + 650d0*p2(1)
  write(*,'(a,2es14.6)') ' detailed balance, Lyb: tau (esc + p_O) and A + tau_O p_O(H):',lhs,rhs
  if(abs(lhs/rhs - 1d0) > 1d-8) nfail = nfail + 1
!-- 3. a blend (escape 1/(tau1+tau2)), and lines 60 km/s apart (independent: 2.560e-2)
  call nx_ovl_zone(2,[0d0,0d0],[20d0*km,20d0*km],[1d3,1d3],1,1d0,1d0,esc,bet1,p)
  write(*,'(a,es11.4,a,f8.5)') ' blend of two lines of depth 1e3: escape',esc,' (5.000e-04), absorbed by the other',p(2)
  if(abs(esc/5d-4 - 1d0) > 1d-3 .or. abs(p(2) - 0.49975d0) > 1d-3) nfail = nfail + 1
  call nx_ovl_zone(2,[0d0,60d0*km],[20d0*km,20d0*km],[1d4,1d4],1,1d0,1d0,esc,bet1,p)
  write(*,'(a,es11.4,a)') ' lines 3 Doppler widths apart, depth 1e4: absorbed by the other',p(2),' (2.560e-02)'
  if(abs(p(2)/2.560d-2 - 1d0) > 1d-2) nfail = nfail + 1
!-- 4. with the model atoms
  nargs = command_argument_count()
  if(nargs >= 1) then
     call get_command_argument(1,adir)
     call atom_tests(trim(adir))
  else
     write(*,'(a)') ' (no model-atom directory: the statistical-equilibrium tests skipped)'
  endif
  if(nfail > 0) then
     write(*,'(a,i2,a)') ' test_overlap:',nfail,' checks FAILED'
     stop 1
  endif
  write(*,'(a)') ' test_overlap: all checks passed'
contains
  subroutine atom_tests(dir)
    character(*),intent(in) :: dir
    character(256) :: files(2)
    real*8,allocatable :: x(:),xl(:),xo(:),beta(:),jb(:),nti(:),ntx(:),nel(:)
    logical,allocatable :: ljb(:)
    real*8 :: t,ne,nh,w,b3d(2)
    integer :: k,ierr,o,s,kb,ko,ir
    logical :: lex
    files(1) = dir//'/h1.atom'
    files(2) = dir//'/o1.atom'
    inquire(file=files(2),exist=lex)
    if(.not.lex) then
       write(*,'(a)') ' (no o1.atom in the model-atom directory: the statistical-equilibrium tests skipped)'
       return
    endif
    do i=1,ng+1
       wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
    enddo
    do i=1,ng
       nug(i) = 0.5d0*(cl/wlg(i) + cl/wlg(i+1))
       dnug(i) = cl/wlg(i) - cl/wlg(i+1)
    enddo
    nx_ovl_dv = 150d0*km
    call nxt_init(1,1,ng,wlg,2,files,0)
    nx_ovl_vtb = 20d0*km
    allocate(x(nxa_nst),xl(nxa_nst),xo(nxa_nst),beta(nxa_nlin),jb(nxa_nlin),ljb(nxa_nlin))
    allocate(nti(nxa_nst),ntx(nxa_nst),nel(nxa_nel))
    nti = 0d0
    ntx = 0d0
!-- Lyman beta and O I 1025.76 (2p4 3P2 - 3d 3D3)
    kb = 0
    ko = 0
    do k=1,nxa_nlin
       if(nxa_el(nxa_lel(k))%z==1 .and. abs(cl/nxa_lnu(k)*1d8 - 1025.72d0) < 0.05d0) kb = k
       if(nxa_el(nxa_lel(k))%z==8 .and. abs(cl/nxa_lnu(k)*1d8 - 1025.762d0) < 0.002d0 .and. ko==0) ko = k
    enddo
    write(*,'(a,i6,a,i4,a,i6)') ' model atoms: pairs of overlapping lines within 150 km/s:',nxo_npair/2, &
         '   Lyb',kb,'  O I 1025.76',ko
    if(kb==0 .or. ko==0 .or. nxo_npair<=0) then
       nfail = nfail + 1
       return
    endif
!-- LTE limit: high density, Planck field, slow flow (thick lines)
    t = 8d3
    nh = 1d17
    do k=1,ng
       jg(k) = bplanck(nug(k),t)
    enddo
    jb = [(bplanck(nxa_lnu(k),t),k=1,nxa_nlin)]
    ljb = .true.
    nel = nh*1d-3
    nel(1) = nh
    call nx_lte(t,nh,x)
    call nx_solve_cell(t,t,nel,0d0,1d0/1.7d6,1d0/1.7d6,ng,nug,dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr)
    call nx_lte(t,ne,xl)
    dev = 0d0
    do s=1,nxa_nst
       if(xl(s) > 1d-10) dev = max(dev,abs(x(s)/xl(s) - 1d0))
    enddo
    write(*,'(a,es10.3,a,i2)') ' LTE limit with the overlap terms: max |n/n_LTE - 1| =',dev,'  ierr',ierr
    if(dev > 2d-2 .or. ierr/=0) nfail = nfail + 1
!-- dilute field (W = 1e-2 of 6000 K), n_H = 1e9, a coasting wind: O I 3d 3D
!-- relative to the O I ground level with and without the overlap
    t = 7.3d3
    nh = 1d9
    w = 1d-2
    do k=1,ng
       jg(k) = w*bplanck(nug(k),6d3)
    enddo
    jb = [(w*bplanck(nxa_lnu(k),6d3),k=1,nxa_nlin)]
    nel = nh*8.5d-4
    nel(1) = nh
    do ir=1,2
       if(ir==2) nxo_npair = 0
       call nx_lte(t,nh,x)
       call nx_solve_cell(t,6d3,nel,0d0,1d-3*2.3d7/3d15,2.3d7/3d15,ng,nug,dnug,jg,jb,ljb,nti,ntx,x,ne,beta, &
            ierr,maxit=1000)
       o = nxa_el(2)%xoff
       b3d(ir) = x(o+nxa_lsu(ko))/x(o+nxa_lsl(ko))
       if(ir==1) xo = x
    enddo
    nxo_npair = nxo_ptr(nxa_nlin+1) - 1
    write(*,'(a,es10.3,a,es10.3,a,f6.3,a,f6.3)') ' dilute field: n(O I 3d 3D)/n(O I 3P) with overlap',b3d(1), &
         '  without',b3d(2),'   x(H I) with',xo(1),' without',x(1)
    if(.not.(b3d(1) > 3d0*b3d(2))) nfail = nfail + 1
  end subroutine atom_tests
end program test_overlap
