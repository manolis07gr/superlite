!-- Tests of the multi-species NLTE treatment (nlteatommod, nltexmod,
!-- nltextrmod):
!--  1. built-in hydrogen reproduces the H-only treatment (nltehmod):
!--     recombination coefficients, photoionization rates, populations in the
!--     LTE limit and in a Case B nebula
!--  2. with He model atoms (directory given as the first argument; skipped
!--     without): LTE limit at high density, He I 5876/4471 in a Case B
!--     nebula, and the macro-atom: cumulative tables, and the 5876/4471
!--     energy ratio of recombination activations against the populations
program test_nltex
  use nltehmod
  use nlteatommod
  use nltexmod
  use nltextrmod
  use randommod
  use gridmod, only:grd_tempinv
  implicit none
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
  real*8,parameter :: pi_ = 3.141592653589793d0
  integer,parameter :: nlev = 30, ng = 6000
  real*8 :: wlg(ng+1), nug(ng), dnug(ng), jg(ng)
  real*8 :: xo(0:nlev), betao(nlev*(nlev-1)/2), jbo(nlev*(nlev-1)/2), neo
  logical :: ljbo(nlev*(nlev-1)/2)
  real*8,allocatable :: x(:), xl(:), beta(:), jb(:), nti(:), ntx(:), nel(:)
  real*8,allocatable :: cq(:,:), ci(:), c3(:), gam(:), gst(:), ar(:)
  logical,allocatable :: ljb(:)
  real*8 :: gamo(nlev), alpo(nlev), alpso(nlev), cupo(nlev), c3bo(nlev), qupo(nlev,nlev), qdno(nlev,nlev)
  real*8 :: jlo(nlev*(nlev-1)/2)
  type(nxion) :: ions(3)
  real*8 :: t, nh, w, t2, q, ne, dev, d, tl(3), nl(3)
  integer :: i, k, n, ierr, nfail, nargs
  character(256) :: adir
  nfail = 0
!-- group grid as in the A4 test (1 A - 30000 A, logarithmic)
  do i = 1, ng+1
     wlg(i) = 1d-8*(30d-5/1d-8)**((i-1d0)/ng)
  enddo
  do i = 1, ng
     nug(i) = 0.5d0*(cc/wlg(i) + cc/wlg(i+1))
     dnug(i) = cc/wlg(i) - cc/wlg(i+1)
  enddo
!
!-- 1. built-in hydrogen
  call nlteh_init(nlev)
  call nxa_builtin_h(nlev, ions(1))
  call nxa_setup(1, ions(1:1))
  call nxa_groups(ng, wlg)
  allocate(x(nxa_nst), xl(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst))
  allocate(nel(nxa_nel))
!-- recombination coefficients against a fine integration in ln u (the
!-- integrand of high levels at high T falls steeply above the edge)
  dev = 0d0
  d = 0d0
  tl = [5d3, 1d4, 3d4]
  do i = 1, 3
     t = tl(i)
     do n = 1, nlev, 4
        w = alpha_ref(n, t)
        dev = max(dev, abs(nxa_alpha(n, t)/w - 1d0))
        d = max(d, abs(nlteh_alpha(n, t)/w - 1d0))
     enddo
  enddo
  write(*,'(a,es10.3,a,es10.3,a)') '  max |alpha/alpha_ref - 1| =', dev, '   (nltehmod:', d, ')'
  if(dev > 3d-3) nfail = nfail + 1
!-- photoionization rates in a diluted 5e4 K field
  w = 1d-14; t2 = 5d4
  do k = 1, ng
     jg(k) = w*bplanck(nug(k), t2)
  enddo
  jbo = 0d0; ljbo = .false.
  call nlteh_cell_rates(1d4, t2, ng, nug, dnug, jg, jbo, ljbo, 0d0, gamo, alpo, alpso, cupo, c3bo, qupo, qdno, jlo)
  n = nxa_el(1)%ns
  allocate(cq(n,n), ci(n), c3(n), gam(n), gst(n), ar(n))
  call nx_el_rates(1, 1d4, t2, ng, nug, dnug, jg, cq, ci, c3, gam, gst, ar)
  dev = maxval(abs(gam(1:nlev)/gamo - 1d0))
  d = maxval(abs(cq(1,2:nlev)/qupo(1,2:nlev) - 1d0))
  write(*,'(a,es10.3,a,es10.3)') '  max |Gamma/Gamma(nltehmod) - 1| =', dev, '   collisions:', d
  if(dev > 1d-2 .or. d > 1d-2) nfail = nfail + 1
!-- LTE limit and Case B: populations against nlteh_solve_cell
  tl = [5d3, 1d4, 2d4]; nl = [1d17, 1d17, 1d18]
  do i = 1, 4
     if(i <= 3) then
        t = tl(i); nh = nl(i); t2 = t; q = 1d0/1.7d6
        do k = 1, ng
           jg(k) = bplanck(nug(k), t)
        enddo
        jbo = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]
        ljbo = .true.
     else
!-- Case B: nH = 1e4, T = 1e4 K, ionizing diluted 5e4 K blackbody only
        t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-9
        jg = 0d0
        do k = 1, ng
           if(nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nug(k), t2)
        enddo
        jbo = 0d0; ljbo = .true.
     endif
     xo = 0d0; xo(0) = 0.5d0; xo(1) = 0.5d0
     call nlteh_solve_cell(t, t2, nh, 0d0, q, q, ng, nug, dnug, jg, jbo, ljbo, 0d0, xo, neo, betao, ierr)
     x = 0d0; x(1) = 0.5d0; x(nxa_nst) = 0.5d0
     jb = jbo; ljb = ljbo; nti = 0d0; ntx = 0d0; nel = nh
     call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
     dev = max(abs(x(nxa_nst)/xo(0) - 1d0), maxval(abs(x(1:nlev)/xo(1:nlev) - 1d0)))
     write(*,'(a,i2,a,es9.2,a,f7.0,a,i2,a,es10.3,a,es10.3)') '  H case', i, ' nH=', nh, ' T=', t, '  ierr=', &
          ierr, '  x(HII)=', x(nxa_nst), '  max |x/x(nltehmod) - 1| =', dev
     if(i <= 3) then
        if(ierr /= 0 .or. dev > 1d-2) nfail = nfail + 1
     else
!-- Case B Balmer decrement (Storey & Hummer 1995: 2.86)
        k = ions(1)%lin(2,3); n = ions(1)%lin(2,4)
        w = (x(3)*nxa_la(k)*beta(k)*nxa_lnu(k))/(x(4)*nxa_la(n)*beta(n)*nxa_lnu(n))
        write(*,'(a,f7.3,a)') '  H Case B H-alpha/H-beta =', w, '   (Storey & Hummer 1995: 2.86)'
        if(ierr /= 0 .or. abs(w/2.86d0 - 1d0) > 0.05d0) nfail = nfail + 1
     endif
  enddo
  deallocate(x, xl, beta, jb, ljb, nti, ntx, nel, cq, ci, c3, gam, gst, ar)
!
!-- 2. He model atoms
  nargs = command_argument_count()
  if(nargs < 1) then
     write(*,'(a)') '  (no model-atom directory given: He tests skipped)'
  else
     call get_command_argument(1, adir)
     call he_tests(trim(adir))
     call all_lte(trim(adir))
  endif
  if(nfail > 0) then
     write(*,'(a,i3)') 'test_nltex: FAILED checks:', nfail
     stop 1
  endif
  write(*,'(a)') 'test_nltex: PASS'
contains
!
  subroutine he_tests(dir)
    character(*), intent(in) :: dir
    character(256) :: files(2)
    integer :: ihe, i2, o, k5876, k4471, ndraw, n58, n44, kind, idx, j, s, iel, ns
    real*8 :: yhe, r, rmc, rse, err, wlc, j58, j44, wls(6)
    type(rnd_t) :: st
    files(1) = dir//'/he1.atom'
    files(2) = dir//'/he2.atom'
    call nxt_init(1, 1, ng, wlg, 2, files, nlev)
    allocate(x(nxa_nst), xl(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst))
    allocate(nel(nxa_nel))
    ihe = 0
    do iel = 1, nxa_nel
       if(nxa_el(iel)%z == 2) ihe = iel
    enddo
    write(*,'(a,i3,a,i5,a,i6,a,i5)') '  He model atoms: elements', nxa_nel, '  states', nxa_nst, &
         '  lines', nxa_nlin, '  levels', nxa_nlv
!-- LTE limit at high density
    yhe = 0.1d0
    do j = 1, 2
       t = 1.2d4; if(j == 2) t = 3d4
       nh = 1d17
       do k = 1, ng
          jg(k) = bplanck(nug(k), t)
       enddo
       jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.
       nel = 0d0; nel(1) = nh; nel(ihe) = yhe*nh
       call nx_lte(t, nh, x)
       nti = 0d0; ntx = 0d0
       call nx_solve_cell(t, t, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
       call nx_lte(t, ne, xl)
       o = nxa_el(ihe)%xoff
       dev = 0d0
       do s = 1, nxa_el(ihe)%ns
          if(xl(o+s) > 1d-12) dev = max(dev, abs(x(o+s)/xl(o+s) - 1d0))
       enddo
       write(*,'(a,f7.0,a,i2,a,3es10.3,a,es10.3)') '  He LTE test T=', t, '  ierr=', ierr, '  He I,II,III:', &
            ionfr(o, 1), ionfr(o, 2), x(o+nxa_el(ihe)%ns), '  max |n/n_LTE - 1| =', dev
       if(ierr /= 0 .or. dev > 1d-2) nfail = nfail + 1
    enddo
!-- Case B nebula: He-ionizing diluted 5e4 K blackbody, resonance lines thick
    t = 1d4; nh = 1d4; w = 1d-15; t2 = 5d4; q = 1d-12
    jg = 0d0
    do k = 1, ng
       if(nug(k) > nh_nu0(1)) jg(k) = w*bplanck(nug(k), t2)
    enddo
    jb = 0d0; ljb = .true.
    nel = 0d0; nel(1) = nh; nel(ihe) = yhe*nh
    x = 0d0
    do iel = 1, nxa_nel
       x(nxa_el(iel)%xoff+nxa_el(iel)%ns) = 0.9d0
       x(nxa_el(iel)%xoff+1) = 0.1d0
    enddo
    nti = 0d0; ntx = 0d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
    o = nxa_el(ihe)%xoff
    k5876 = lnear(5875.6d-8); k4471 = lnear(4471.5d-8)
    j58 = popu(k5876)*nxa_la(k5876)*beta(k5876)*nxa_lnu(k5876)
    j44 = popu(k4471)*nxa_la(k4471)*beta(k4471)*nxa_lnu(k4471)
    r = j58/j44
    write(*,'(a,i2,a,3es10.3,a,es10.3)') '  He Case B: ierr=', ierr, '  He I,II,III:', ionfr(o,1), ionfr(o,2), &
         x(o+nxa_el(ihe)%ns), '  ne=', ne
    write(*,'(a,f7.3,a)') '  He I 5876/4471 =', r, '   (Case B at 1e4 K, low density: ~2.7-2.8)'
    wls = [584.33d-8, 10830.3d-8, 3888.6d-8, 5875.6d-8, 4471.5d-8, 6678.2d-8]
    do j = 1, 6
       k = lnear(wls(j))
       write(*,'(a,f9.1,a,es10.3,a,es10.3,a,es10.3)') '    He I', nxa_lwl(k)*1d8, ' tau=', nx_kq(k, nel(ihe), &
            x(o+1:o+nxa_el(ihe)%ns))/q, ' beta=', beta(k), ' x_u=', popu(k)
    enddo
    if(ierr /= 0 .or. r < 2.4d0 .or. r > 3.2d0) nfail = nfail + 1
!-- macro-atom against the populations (Lucy 2002): exact when the
!-- activations carry all energy input; recombination activations do so only
!-- without collisional excitation (which the thermal pool drives), so the
!-- comparison is made with the electron-impact rates off
    nxa_colfac = 0d0
    call nx_solve_cell(t, t2, nel, 0d0, q, q, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
    j58 = popu(k5876)*nxa_la(k5876)*beta(k5876)*nxa_lnu(k5876)
    j44 = popu(k4471)*nxa_la(k4471)*beta(k4471)*nxa_lnu(k4471)
!-- macro-atom tables of this cell
    nxt_x(:,1) = x; nxt_nel(:,1) = nel; nxt_ne(1) = ne
    nxt_jg(:,1) = jg; nxt_jb(:,1) = jb; nxt_ljb(:,1) = ljb; nxt_have_mc = .true.
    allocate(grd_tempinv(1)); grd_tempinv = 1d0/t
    call nxt_tables(1, [1d0/t], [t2], [0d0], [1d15, 1.001d15], [q*1d15, q*1.001d15])
    do iel = 1, nxa_nel
       ns = nxa_el(iel)%ns
       do s = 1, ns
          j = nxt_ecoff(iel) + (s-1)*(ns+1)
          if(abs(nxt_ecum(j+ns+1,1) - 1d0) > 1d-12 .or. &
               any(nxt_ecum(j+2:j+ns+1,1) < nxt_ecum(j+1:j+ns,1) - 1d-14)) then
             write(*,*) 'bad exit table: element', iel, ' state', s
             nfail = nfail + 1
          endif
       enddo
    enddo
    if(abs(nxt_kcum(nxt_nko,1) - 1d0) > 1d-12) nfail = nfail + 1
!-- recombination activations of He I: macro-atom in the He II ground state
    i2 = nxa_ion(nxa_el(ihe)%ion(2))%soff + 1
    call rnd_init(1, 7); st = rnd_states(1)
    ndraw = 2000000
    n58 = 0; n44 = 0
    do j = 1, ndraw
       call nxt_activate(1, ihe, i2, st, kind, idx, wlc)
       if(kind == em_line .and. idx == k5876) n58 = n58 + 1
       if(kind == em_line .and. idx == k4471) n44 = n44 + 1
    enddo
!-- packets carry equal energy: the count ratio is the energy ratio
    rmc = dble(n58)/dble(max(n44,1))
    err = rmc*sqrt(1d0/max(n58,1) + 1d0/max(n44,1))
    rse = j58/j44
    write(*,'(a,f8.4,a,f7.4,a,f8.4)') '  He I 5876/4471 energy from the macro-atom:', rmc, ' +-', err, &
         '   from the populations:', rse
    if(abs(rmc - rse) > 4d0*err + 0.02d0*rse) nfail = nfail + 1
!-- the same for H-alpha/H-beta (H recombination)
    k5876 = ions_lin(1, 2, 3); k4471 = ions_lin(1, 2, 4)
    n58 = 0; n44 = 0
    do j = 1, ndraw
       call nxt_activate(1, 1, nxa_el(1)%ns, st, kind, idx, wlc)
       if(kind == em_line .and. idx == k5876) n58 = n58 + 1
       if(kind == em_line .and. idx == k4471) n44 = n44 + 1
    enddo
    rmc = dble(n58)/dble(max(n44,1))
    err = rmc*sqrt(1d0/max(n58,1) + 1d0/max(n44,1))
    rse = (x(3)*nxa_la(k5876)*nxt_beta(k5876,1)*nxa_lnu(k5876))/(x(4)*nxa_la(k4471)*nxt_beta(k4471,1)*nxa_lnu(k4471))
    write(*,'(a,f8.4,a,f7.4,a,f8.4)') '  H-alpha/H-beta energy from the macro-atom:', rmc, ' +-', err, &
         '   from the populations:', rse
    if(abs(rmc - rse) > 4d0*err + 0.02d0*rse) nfail = nfail + 1
    nxa_colfac = 1d0
  end subroutine he_tests
!
  integer function ions_lin(iel, l, u)
!-- global index of the line (l,u) of the first ion of element iel
    integer, intent(in) :: iel, l, u
    integer :: ii
    ii = nxa_el(iel)%ion(1)
    ions_lin = nxa_ion(ii)%loff + nxa_ion(ii)%lin(l,u)
  end function ions_lin
!
  subroutine all_lte(dir)
!-- every model atom of the directory (list in atoms.txt, if present): the
!-- populations of all elements tend to Saha-Boltzmann at high density in a
!-- Planck field (detailed balance of the converted atomic data)
    character(*), intent(in) :: dir
    character(256) :: files(80), line
    integer :: nf, ios, iel, s, o, j, nbad, ncx
    real*8 :: devel, dmax
    integer :: cxf(4096), cxt(4096)
    real*8 :: cxr(4096), fnet(2000), fabs(2000)
    nf = 0
    open(9, file=dir//'/atoms.txt', status='old', iostat=ios)
    if(ios /= 0) then
       write(*,'(a)') '  (no atoms.txt in the model-atom directory: all-element LTE test skipped)'
       return
    endif
    do
       read(9,'(a)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       nf = nf + 1
       files(nf) = dir//'/'//trim(line)
    enddo
    close(9)
    call nxt_init(1, 1, ng, wlg, nf, files(:nf), nlev)
    call nxa_read_cx(dir//'/chgexch.dat', .true.)
    if(allocated(x)) deallocate(x, xl, beta, jb, ljb, nti, ntx, nel)
    allocate(x(nxa_nst), xl(nxa_nst), beta(nxa_nlin), jb(nxa_nlin), ljb(nxa_nlin), nti(nxa_nst), ntx(nxa_nst))
    allocate(nel(nxa_nel))
    write(*,'(a,i3,a,i3,a,i5,a,i7)') '  all atoms:', nxa_nel, ' elements', nxa_nion, ' ions', nxa_nst, &
         ' states', nxa_nlin
    nbad = 0
    do j = 1, 2
       t = 8d3; if(j == 2) t = 2.5d4
       nh = 1d17
       do k = 1, ng
          jg(k) = bplanck(nug(k), t)
       enddo
       jb = [(bplanck(nxa_lnu(k), t), k=1,nxa_nlin)]; ljb = .true.
       do iel = 1, nxa_nel
          nel(iel) = nh*1d-4
       enddo
       nel(1) = nh
       call nx_lte(t, nh, x)
       nti = 0d0; ntx = 0d0
       call nx_solve_cell(t, t, nel, 0d0, 1d0/1.7d6, 1d0/1.7d6, ng, nug, dnug, jg, jb, ljb, nti, ntx, x, ne, beta, ierr)
       call nx_lte(t, ne, xl)
!-- charge exchange: no net flow into any state at LTE (detailed balance)
       if(nxa_ncx > 0) then
          dmax = 0d0
          do iel = 1, nxa_nel
             o = nxa_el(iel)%xoff
             call nx_cx_rates(iel, t, nel, xl, ncx, cxf, cxt, cxr)
             fnet = 0d0; fabs = 0d0
             do k = 1, ncx
                fnet(cxt(k)) = fnet(cxt(k)) + xl(o+cxf(k))*cxr(k)
                fnet(cxf(k)) = fnet(cxf(k)) - xl(o+cxf(k))*cxr(k)
                fabs(cxf(k)) = fabs(cxf(k)) + xl(o+cxf(k))*cxr(k)
             enddo
             do k = 1, nxa_el(iel)%ns
                if(fabs(k) > 0d0) dmax = max(dmax, abs(fnet(k))/fabs(k))
             enddo
          enddo
          write(*,'(a,f7.0,a,i3,a,es9.2)') '  charge exchange at LTE, T=', t, ':', nxa_ncx, &
               ' reactions, max |net/out| per state =', dmax
          if(dmax > 1d-8) nbad = nbad + 1
       endif
       do iel = 1, nxa_nel
          o = nxa_el(iel)%xoff
          devel = 0d0
          do s = 1, nxa_el(iel)%ns
             if(xl(o+s) > 1d-10) devel = max(devel, abs(x(o+s)/xl(o+s) - 1d0))
          enddo
          if(devel > 2d-2 .or. ierr /= 0) then
             nbad = nbad + 1
             write(*,'(a,f7.0,a,i3,a,i2,a,es10.3)') '    T=', t, '  Z=', nxa_el(iel)%z, '  ierr=', ierr, &
                  '  max |n/n_LTE - 1| =', devel
          endif
       enddo
    enddo
    write(*,'(a,i3)') '  all-element LTE test, elements off LTE:', nbad
    if(nbad > 0) nfail = nfail + 1
  end subroutine all_lte
!
  real*8 function ionfr(o, j)
    integer, intent(in) :: o, j
    integer :: ii
    ionfr = 0d0
    do ii = 1, nxa_nion
       if(nxa_el(nxa_ion(ii)%el)%xoff /= o) cycle
       if(nxa_ion(ii)%stage /= j) cycle
       ionfr = sum(x(o+nxa_ion(ii)%soff+1:o+nxa_ion(ii)%soff+nxa_ion(ii)%nlev))
    enddo
  end function ionfr
!
  integer function lnear(wl)
!-- He line closest to wl [cm]
    real*8, intent(in) :: wl
    integer :: kk
    real*8 :: dmin
    dmin = 1d99
    lnear = 0
    do kk = 1, nxa_nlin
       if(nxa_el(nxa_lel(kk))%z /= 2) cycle
       if(abs(nxa_lwl(kk) - wl) < dmin) then
          dmin = abs(nxa_lwl(kk) - wl)
          lnear = kk
       endif
    enddo
  end function lnear
!
  real*8 function alpha_ref(nn, tt)
!-- Milne integral of nlteh_sigma in ln u (trapezoid, 200000 points to y = 60)
    integer, intent(in) :: nn
    real*8, intent(in) :: tt
    integer :: jj
    integer, parameter :: np = 200000
    real*8 :: nu0, ykt, lu, umax, u, y, f, fo, dl, ss
    nu0 = nh_nu0(nn)
    ykt = hh*nu0/(kb*tt)
    umax = 1d0 + 60d0/ykt
    dl = log(umax)/np
    ss = 0d0
    fo = 0d0
    do jj = 0, np
       lu = jj*dl
       u = exp(lu)
       y = (u - 1d0)*ykt
       f = nlteh_sigma(nn, u*nu0)*(u*nu0)**2*exp(-y)*u*nu0     !d nu = nu d ln u
       if(jj > 0) ss = ss + 0.5d0*(f + fo)*dl
       fo = f
    enddo
    alpha_ref = 4d0*pi_*dble(nn*nn)*4.1412d-16*tt**(-1.5d0)*2d0/cc**2*ss
  end function alpha_ref
!
  real*8 function popu(kk)
!-- population (fraction of the element) of the upper level of line kk
    integer, intent(in) :: kk
    popu = x(nxa_el(nxa_lel(kk))%xoff + nxa_lsu(kk))
  end function popu
end program test_nltex
