!-- Super-levels (nlteatommod): the statistical equilibrium of the same Fe II-IV
!-- atoms (and He) with super-levels and with every level a state, in four
!-- regimes (photospheric, diluted nebular, dense, cool).  The Fe ion fractions
!-- and Fe II level populations agree to 10 %.
!--   usage: test_superlevel <atoms with super-levels> <the same without>
!-- (both directories with he1, he2, fe2, fe3, fe4.atom, made with
!-- Tools/make_nlte_atoms.sh: e.g. fe2 400:f_to_s_SN_A, fe3 300:f_to_s_105,
!-- fe4 300:f_to_s_100 and GFMIN=1e-4, and the same table without :f_to_s)
program test_superlevel
  use nltehmod
  use nlteatommod
  use nltexmod
  use nltextrmod
  implicit none
  real*8,parameter :: cc = 2.99792458d10
  integer,parameter :: ng = 6000, nc = 4
  real*8 :: wlg(ng+1), nug(ng), dnug(ng), jg(ng)
  character(256) :: files(5), dirs(2)
  real*8 :: res(2,nc,12), tt(nc), nn(nc), ww(nc), tr(nc)
  integer :: id, ic, i, k, ierr, ife, j, lv, nbad
  real*8,allocatable :: x(:), beta(:), jb(:), nti(:), ntx(:), nel(:), wv(:)
  logical,allocatable :: ljb(:)
  real*8 :: ne, q
  integer,parameter :: lsel(7) = [1, 6, 10, 30, 63, 100, 200]
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i = 1, ng
     nug(i) = 0.5d0*(cc/wlg(i) + cc/wlg(i+1))
     dnug(i) = cc/wlg(i) - cc/wlg(i+1)
  enddo
  if(command_argument_count() < 2) then
     write(*,'(a)') '  (no atom directories given: super-level test skipped)'
     stop
  endif
  call get_command_argument(1, dirs(1))
  call get_command_argument(2, dirs(2))
  tt = [1d4, 8d3, 1.2d4, 6d3]; nn = [1d10, 1d8, 1d13, 1d9]; ww = [0.5d0, 0.1d0, 1d0, 0.3d0]
  tr = [1.2d4, 1.5d4, 1.2d4, 9d3]
  do id = 1, 2
     files(1) = trim(dirs(id))//'/he1.atom'; files(2) = trim(dirs(id))//'/he2.atom'
     files(3) = trim(dirs(id))//'/fe2.atom'; files(4) = trim(dirs(id))//'/fe3.atom'
     files(5) = trim(dirs(id))//'/fe4.atom'
     call nxt_init(1, 1, ng, wlg, 5, files, 30)
     if(allocated(x)) deallocate(x, beta, jb, ljb, nti, ntx, nel, wv)
     allocate(x(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst), &
          nel(nxa_nel), wv(nxa_nlv))
     ife = 0
     do i = 1, nxa_nel
        if(nxa_el(i)%z == 26) ife = i
     enddo
     write(*,*) 'set', id, ' states', nxa_nst, ' Fe states', nxa_el(ife)%ns, ' lines', nxa_nlin
     do ic = 1, nc
        do k = 1, ng
           jg(k) = ww(ic)*bplanck(nug(k), tr(ic))
        enddo
        jb = 0d0; ljb = .false.
        nel = 0d0; nel(1) = nn(ic); nel(2) = 0.1d0*nn(ic); nel(ife) = 3d-5*nn(ic)
        call nx_lte(tt(ic), nn(ic), x)
        nti = 0d0; ntx = 0d0
        q = 1d0/1d6
        call nx_solve_cell(tt(ic), tr(ic), nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
        call nxa_weights(tt(ic), wv)
        res(id,ic,:) = 0d0
        do j = 1, nxa_el(ife)%nion
           i = nxa_el(ife)%ion(j)
           res(id,ic,j) = sum(x(nxa_el(ife)%xoff+nxa_ion(i)%soff+1:nxa_el(ife)%xoff+nxa_ion(i)%soff+nxa_ion(i)%nsl))
        enddo
        res(id,ic,4) = x(nxa_el(ife)%xoff+nxa_el(ife)%ns)
        res(id,ic,5) = ne
!-- Fe II levels 1, 6, 10, 30, 63, 100, 200 (populations, fraction of Fe)
        i = nxa_el(ife)%ion(1)
        do lv = 1, 7
           j = lsel(lv)
           res(id,ic,5+lv) = x(nxa_el(ife)%xoff+nxa_lvs(nxa_ion(i)%lvoff+j))*wv(nxa_ion(i)%lvoff+j)
        enddo
        write(*,'(a,i2,a,i2,a,4es10.3,a,es10.3)') ' set', id, ' case', ic, ' Fe II-V:', res(id,ic,1:4), ' ne', ne
     enddo
  enddo
  write(*,'(a)') ' ratio SL/full: Fe II, III, IV, top, ne, FeII lev 1,6,10,30,63,100,200'
  nbad = 0
  do ic = 1, nc
     write(*,'(a,i2,12f8.3)') ' case', ic, res(1,ic,:)/res(2,ic,:)
!-- ion fractions above 1e-3 and the level populations
     do k = 1, 12
        if(k <= 4 .and. res(2,ic,k) < 1d-3) cycle
        if(abs(res(1,ic,k)/res(2,ic,k) - 1d0) > 0.1d0) nbad = nbad + 1
     enddo
  enddo
  if(nbad > 0) then
     write(*,'(a,i3)') 'test_superlevel: FAILED checks:', nbad
     stop 1
  endif
  write(*,'(a)') 'test_superlevel: PASS'
end program test_superlevel
