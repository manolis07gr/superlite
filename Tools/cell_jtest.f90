!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
program cell_jtest
!-----------------------------------------------------------------------
! SuperLite's statistical equilibrium of cells of a run in an external
! radiation field, against the same cells with the run's own Monte Carlo
! J: are the rate equations or the radiation field behind a difference
! from another code?  The external J(nu) per cell (cmfj_<cell>.dat: nu
! [Hz], J; Tools/cmfgen_j_cells.py writes CMFGEN's) gives the group J
! (mean over each group) and the incident intensity of each Sobolev line
! (J five Doppler widths, 20 km/s, to its blue).  The Lyman continuum and
! shortward stay in detailed balance with the current populations: the
! run's own J there (local detailed balance held it at S of the saved
! state) with Lambda* = 1, so the external J does not set the ground
! states of H, He directly.  Starts from the run's saved state
! (output.nltex_state, output.nltex_jmc; the run's model atoms, n_e as
! saved, no non-thermal or charge-exchange-free changes).
!
! usage: cell_jtest <run dir> <cell,cell,...> <J dir> [mode [wl_max [dis]]]
!   mode: all (default); fuv: the external J below wl_max [A] only (groups
!   and lines; default 1250); cont: the groups only; lines: the line
!   intensities only; hlines: those of the H lines only.  dis: level
!   dissolution (in_nltex_dis).  dump (7th argument): the net flows into
!   the ground state of H and the ionization flows of each element of
!   both solutions (as SL_RATEDUMP) to cell_jtest_dump_<cell>.
! build (as test_nebular in Tests/unit/run_tests.sh), with a CMake build B:
!   mpif90 -O2 -fallow-argument-mismatch -I$B Tools/cell_jtest.f90 \
!     $(find $B/CMakeFiles/superlite.dir -name '*.o' ! -name superlite.f90.o) -llapack -lblas
!-----------------------------------------------------------------------
  use nlteatommod
  use nltexmod
  use nltextrmod
  use nltedismod, only:nxd_on,nxd_setup,nxd_cellx
  implicit none
  integer,parameter :: mf = 64
  character(256) :: run,line,arg,jdir,mode,adir,cxfile,fname
  character(64) :: an(mf)
  character(256) :: files(mf)
  character(2048) :: atoms
  character(8) :: sym(40)
  character(2),parameter :: zs(30) = [character(2) :: 'h','he','li','be','b','c','n','o','f','ne','na','mg', &
       'al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn']
  real*8,parameter :: am(30) = [1.008d0,4.0026d0,6.94d0,9.012d0,10.81d0,12.011d0,14.007d0,15.999d0,18.998d0, &
       20.180d0,22.990d0,24.305d0,26.982d0,28.085d0,30.974d0,32.06d0,35.45d0,39.948d0,39.098d0,40.078d0, &
       44.956d0,47.867d0,50.942d0,51.996d0,54.938d0,55.845d0,58.933d0,58.693d0,63.546d0,65.38d0]
  real*8,parameter :: mu = 1.66053907d-24, cl = 2.99792458d10, vd = 20d5
  integer :: i,k,ic,ncell,nst,nel_,ng,ierr,nc,cells(200),ig,n,ios,j,iv,iel,ii,lv,nf,nsym,ncol,z
  integer,allocatable :: col(:)
  real*8,allocatable :: wlg(:),hdr(:,:),str(:,:),xs(:,:),x(:),beta(:),nti(:),ntx(:),nuc(:),jc(:),jg(:),jb(:),nel(:)
  real*4,allocatable :: kap(:)
  logical,allocatable :: ljb(:)
  real*8 :: qb,qv,ne,nu1,nu2,s,w,nu0,dr,xo(2),wlmax,ovl,vturb
  logical :: lg,ll,ldmp,cdest
  integer :: mxit,ierr1,nit1
!
!-- (SL_MAXIT: iteration limit of the statistical equilibrium, default 300)
  mxit = 300
  call get_environment_variable('SL_MAXIT',arg,status=ios)
  if(ios==0 .and. len_trim(arg)>0) read(arg,*) mxit
  call get_command_argument(1,run)
  call get_command_argument(2,arg)
  nc = count([(arg(i:i)==',',i=1,len_trim(arg))]) + 1
  read(arg,*) cells(1:nc)
  call get_command_argument(3,jdir)
  mode = 'all'
  if(command_argument_count() >= 4) call get_command_argument(4,mode)
  wlmax = 1250d0
  if(command_argument_count() >= 5) then
     call get_command_argument(5,arg)
     read(arg,*) wlmax
  endif
  nxd_on = .false.
  if(command_argument_count() >= 6) then
     call get_command_argument(6,arg)
     nxd_on = trim(arg)=='dis'
  endif
  ldmp = .false.
  if(command_argument_count() >= 7) then
     call get_command_argument(7,arg)
     ldmp = trim(arg)=='dump'
  endif
!-- (dump: the flows into one state, element and state index of the
!-- element as the 8th and 9th arguments)
  if(command_argument_count() >= 9) then
     call get_command_argument(8,arg)
     read(arg,*) nx_dump_el
     call get_command_argument(9,arg)
     read(arg,*) nx_dump_st
  endif
!-- the run's model atoms, charge exchange, overlap and microturbulence
  adir = ''
  atoms = ''
  cxfile = 'chgexch.dat'
  ovl = 0d0
  vturb = 0d0
  cdest = .true.
  open(9,file=trim(run)//'/input.par',status='old')
  do
     read(9,'(a)',iostat=ios) line
     if(ios /= 0) exit
     if(index(line,'in_nltex_dir') > 0) adir = quoted(line)
     if(index(line,'in_nltex_atoms') > 0) then
        backspace(9)
        read(9,'(a)') atoms
        atoms = quoted(atoms)
     endif
     if(index(line,'in_nltex_chgexch') > 0) cxfile = quoted(line)
     if(index(line,'in_nltex_overlap') > 0) ovl = value(line)
     if(index(line,'in_nltex_vturb ') > 0 .or. index(line,'in_nltex_vturb=') > 0) vturb = value(line)
     if(index(line,'in_nltex_cdest') > 0) cdest = index(line,'= t') > 0 .or. index(line,'=t') > 0 .or. &
          index(line,'.true.') > 0
  enddo
  close(9)
  nf = 0
  an = ''
  read(atoms,*,iostat=ios) an
  do i=1,mf
     if(len_trim(an(i)) == 0 .or. index(an(i),'.atom') == 0) exit
     nf = i
  enddo
  if(nf == 0) stop 'cell_jtest: no in_nltex_atoms in input.par'
  write(*,'(a,a,a,i3,a,f6.1,a,f6.1,a)') ' atoms ',trim(adir),':',nf,' files, overlap',ovl,' km/s, vturb',vturb,' km/s'
  write(*,'(a,a,a,f8.1,a,l2)') ' mode ',trim(mode),'  wl_max',wlmax,'  dissolution',nxd_on
  open(9,file=trim(run)//'/input.wlgrid',status='old')
  read(9,*) ng
  allocate(wlg(ng+1))
  read(9,*) wlg
  close(9)
  nx_ovl_dv = ovl*1d5
  open(9,file=trim(run)//'/output.nltex_state',status='old')
  read(9,'(a)') line
  read(9,*) line,line,line,nel_,nst,ncell
  read(9,'(a)') line
  read(9,'(a)') line
  allocate(hdr(5,ncell),xs(nst,ncell))
  do ic=1,ncell
     read(9,*) hdr(:,ic)
     read(9,*) xs(:,ic)
  enddo
  close(9)
  do i=1,nf
     files(i) = trim(adir)//'/'//trim(an(i))
  enddo
  call nxt_init(ncell,ncell,ng,wlg,nf,files(1:nf),0)
  if(trim(cxfile) /= 'none') call nxa_read_cx(trim(adir)//'/'//trim(cxfile),.false.)
  call nxt_read_init(trim(run)//'/output.nltex_state',ncell,hdr(2,:),.true.)
  if(nxd_on) call nxd_setup(ng,wlg)
!-- the structure and the mass fractions (input.str: columns after lum)
  open(9,file=trim(run)//'/input.str',status='old')
  read(9,'(a)') line
  read(9,'(a)') line
  read(9,'(a)') line
  sym = ''
  read(line(2:),*,iostat=ios) sym
  nsym = count(sym /= '')
  ncol = nsym
  allocate(str(ncol,ncell))
  do ic=1,ncell
     read(9,*) str(:,ic)
  enddo
  close(9)
  allocate(col(nxa_nel),nel(nxa_nel))
  col = 0
  do iel=1,nxa_nel
     z = nxa_el(iel)%z
     do k=9,nsym
        if(trim(sym(k)) == trim(zs(z))) col(iel) = k
     enddo
     if(col(iel) == 0) stop 'cell_jtest: an NLTE element is not in input.str'
  enddo
  allocate(x(nxa_nst),beta(nxa_nlin),nti(nxa_nst),ntx(nxa_nst),jg(ng),jb(nxa_nlin),ljb(nxa_nlin),kap(ng))
  nti = 0d0
  ntx = 0d0
  nx_ovl_vtb = vturb*1d5
  nxt_ldbtau = 30d0
  nxt_ldbone = .true.
!-- continuum destruction of the trapped line photons (in_nltex_cdest, on by
!-- default) with the bound-free opacity of the ground states as the
!-- continuum (the Lyman continuum inside it; the excited levels' continua
!-- are left out)
  nx_pc_on = cdest
  nx_pc_vtb = vturb*1d5
  write(*,'(a)') ' cell        r       x(H 1) own J   x(H 1) ext. J    x(H 2) own J  x(H 2) ext. J'
  do i=1,nc
     ic = cells(i)
     do iel=1,nxa_nel
        nel(iel) = hdr(3,ic)*str(col(iel),ic)/(am(nxa_el(iel)%z)*mu)
     enddo
     qb = (str(4,ic)-str(3,ic))/(str(2,ic)-str(1,ic))
     qv = (str(4,ic)+str(3,ic))/(str(2,ic)+str(1,ic))
     dr = str(2,ic) - str(1,ic)
!-- bound-free opacity of the ground states of all ions (the ALI's shares)
     kap = 0.
     do iel=1,nxa_nel
        do j=1,nxa_el(iel)%nion
           ii = nxa_el(iel)%ion(j)
           lv = nxa_ion(ii)%lvoff + 1
           do ig=1,min(nxa_g1(lv),ng)
              kap(ig) = kap(ig) + real(nel(iel)*xs(nxa_el(iel)%xoff+nxa_lvs(lv),ic)*nxa_sgm(nxa_sgo(lv)+ig))
           enddo
        enddo
     enddo
     if(nx_pc_on) then
        if(.not.allocated(nx_kcg)) allocate(nx_kcg(ng))
        nx_kcg = kap
     endif
!-- the external J at the cell
     write(fname,'(a,a,i0,a)') trim(jdir),'/cmfj_',ic,'.dat'
     open(9,file=trim(fname),status='old')
     n = 0
     do
        read(9,*,iostat=ios) w
        if(ios /= 0) exit
        n = n + 1
     enddo
     rewind(9)
     if(allocated(nuc)) deallocate(nuc,jc)
     allocate(nuc(n),jc(n))
     do k=1,n
        read(9,*) nuc(k),jc(k)
     enddo
     close(9)
     do ig=1,ng
        nu2 = cl/wlg(ig)
        nu1 = cl/wlg(ig+1)
!-- (every group with Lyman-continuum opacity: one that straddles the edge
!-- holds n_1 at the external J's S_LyC as firmly as the continuum itself)
        if(wlg(ig) < 911.76d-8) then
           jg(ig) = nxt_jg(ig,ic)
           cycle
        endif
        kap(ig) = 0.
        lg = trim(mode)=='all' .or. trim(mode)=='cont' .or. (trim(mode)=='fuv' .and. wlg(ig) < wlmax*1d-8)
        if(.not.lg) then
           jg(ig) = nxt_jg(ig,ic)
           cycle
        endif
        s = 0d0
        w = 0d0
        call seg(nu1,nu2,s,w)
        jg(ig) = s/max(w,1d-300)
     enddo
     do k=1,nxa_nlin
        ll = trim(mode)=='all' .or. trim(mode)=='lines' .or. (trim(mode)=='fuv' .and. cl/nxa_lnu(k) < wlmax*1d-8) &
             .or. (trim(mode)=='hlines' .and. nxa_el(nxa_lel(k))%z==1)
        if(ll) then
           nu0 = nxa_lnu(k)*(1d0 + 5d0*vd/cl)
           jb(k) = interp(nu0)
           ljb(k) = .true.
        else
           jb(k) = nxt_jb(k,ic)
           ljb(k) = nxt_ljb(k,ic)
        endif
     enddo
     do iv=1,2
        x = xs(:,ic)
        nx_dis = 0d0
        if(nxd_on) call nxd_cellx(hdr(5,ic),nel,x,hdr(4,ic),nx_dis)
        nx_dump = ldmp
        if(iv==1) then
           call nxt_ali_off
           call nx_solve_cell(hdr(4,ic),hdr(4,ic),nel,0d0,qb,qv,ng,nxt_nug,nxt_dnug,nxt_jg(:,ic),nxt_jb(:,ic), &
                nxt_ljb(:,ic),nti,ntx,x,ne,beta,ierr,maxit=mxit)
           xo = [x(1),x(2)]
           ierr1 = ierr
           nit1 = nx_nits
        else
           call nxt_ali_cell(ng,kap,1d0,dr,lali=.false.)
           call nx_solve_cell(hdr(4,ic),hdr(4,ic),nel,0d0,qb,qv,ng,nxt_nug,nxt_dnug,jg,jb,ljb,nti,ntx,x,ne,beta,ierr, &
                maxit=mxit)
           call nxt_ali_off
        endif
        if(ldmp) call dump(iv)
     enddo
     write(*,'(i5,es12.4,4es15.4,a,2i2,a,2i6)') ic,0.5d0*(str(1,ic)+str(2,ic)),xo(1),x(1),xo(2),x(2),'  ierr',ierr1,ierr, &
          '  its',nit1,nx_nits
  enddo
contains
  subroutine dump(jv)
!-- the flows of the last solution (nx_dump, as SL_RATEDUMP)
    integer,intent(in) :: jv
    integer :: lu,kk,jj
    character(64) :: fn
    real*8 :: wlv(nxa_nlv),netl(nxa_nlin),pl,pu,jlk
    nx_dump = .false.
    write(fn,'(a,i0)') 'cell_jtest_dump_',ic
    open(newunit=lu,file=trim(fn),position='append')
    if(jv==1) write(lu,'(a,i5,a,a,a,f8.1)') '##### cell',ic,'  mode ',trim(mode),'  wl_max',wlmax
    write(lu,'(a,a,a,es11.4,a,es11.4,a,2es11.4)') '# solution with the ',trim(merge('own J     ','external J',jv==1)), &
         '  T',hdr(4,ic),'  n_e',ne,'  x_1, x_2',x(1),x(2)
    write(lu,'(a)') '# Lyman lines and H alpha: lambda [A], beta, beta + bes, J resonance, J e.s., J rates, S, '// &
         'net u->l per H particle, J e.s. expectation'
    do kk=1,size(nx_dlya,2)
       if(nx_dlya(1,kk) > 0d0) write(lu,'(a,f9.2,8es12.4)') '#  ',nx_dlya(:,kk)
    enddo
    write(lu,'(a,es12.4)') '# charge exchange into the ground state:',nx_dcx
    write(lu,'(a)') '# state  net flow into the ground state [1/s per H particle], charge exchange apart'
    do kk=2,nx_dns
       write(lu,'(i6,es14.5)') kk,nx_dflow(kk)
    enddo
!-- lines inside the H Lyman continuum: their net radiative rate
!-- n_u A beta (1 - J_l/S) per volume, which H I absorbs (He I 584, the
!-- Fe III 3d6-3d5 4p lines at 800-870 A in a dense wind)
    write(lu,'(a)') '# lines inside the Lyman continuum: net radiative u->l [cm^-3 s^-1] '// &
         'per ion, and the largest lines (lambda [A], net, beta, J_l, S)'
    call nxa_weights(hdr(4,ic),wlv)
    netl = 0d0
    do kk=1,nxa_nlin
       if(cl/nxa_lnu(kk) >= 911.76d-8) cycle
       jj = nxa_lel(kk)
       pl = nel(jj)*x(nxa_el(jj)%xoff+nxa_lsl(kk))*wlv(nxa_llv(kk))
       pu = nel(jj)*x(nxa_el(jj)%xoff+nxa_lsu(kk))*wlv(nxa_luv(kk))
       jlk = nx_jline(kk,ng,nxt_nug,merge(nxt_jg(:,ic),jg,jv==1),merge(nxt_jb(kk,ic),jb(kk),jv==1), &
            merge(nxt_ljb(kk,ic),ljb(kk),jv==1))
       netl(kk) = beta(kk)*(pu*nxa_la(kk) - (pl*nxa_lblu(kk) - pu*nxa_lbul(kk))*jlk)
    enddo
    do kk=1,nxa_nion
       pl = 0d0
       do jj=1,nxa_nlin
          if(nxa_lion(jj) == kk) pl = pl + netl(jj)
       enddo
       if(pl /= 0d0) write(lu,'(a,a8,es12.4)') '#   ',nxa_ion(kk)%name,pl
    enddo
    do jj=1,10
       kk = maxloc(abs(netl),1)
       if(netl(kk) == 0d0) exit
       pl = nel(nxa_lel(kk))*x(nxa_el(nxa_lel(kk))%xoff+nxa_lsl(kk))*wlv(nxa_llv(kk))
       pu = nel(nxa_lel(kk))*x(nxa_el(nxa_lel(kk))%xoff+nxa_lsu(kk))*wlv(nxa_luv(kk))
       jlk = nx_jline(kk,ng,nxt_nug,merge(nxt_jg(:,ic),jg,jv==1),merge(nxt_jb(kk,ic),jb(kk),jv==1), &
            merge(nxt_ljb(kk,ic),ljb(kk),jv==1))
       write(lu,'(a,a8,f10.3,4es12.4)') '#     ',nxa_ion(nxa_lion(kk))%name,cl/nxa_lnu(kk)*1d8,netl(kk),beta(kk),jlk, &
            pu*nxa_la(kk)/max(pl*nxa_lblu(kk) - pu*nxa_lbul(kk),1d-300)
       netl(kk) = 0d0
    enddo
    if(nx_dump_el > 0 .and. allocated(nx_dflx)) then
       write(lu,'(a,i3,a,i4,a,es12.4,a)') '# element',nx_dump_el,' state',nx_dump_st,' (fraction', &
            x(nxa_el(nx_dump_el)%xoff+nx_dump_st),'): net flows into it from each state [1/s per particle], largest'
       do jj=1,12
          kk = maxloc(abs(nx_dflx),1)
          if(nx_dflx(kk) == 0d0) exit
          write(lu,'(a,i5,es14.5)') '#   ',kk,nx_dflx(kk)
          nx_dflx(kk) = 0d0
       enddo
    endif
    write(lu,'(a)') '# H states: fraction of H'
    write(lu,'(5es13.5)') x(nxa_el(1)%xoff+1:nxa_el(1)%xoff+nxa_el(1)%ns)
    write(lu,'(a)') '# per element and ion stage [1/s per particle of the element]: up, down, net charge exchange up'
    do kk=1,nxa_nel
       do jj=1,nxa_el(kk)%nion
          write(lu,'(2i4,a,3es13.4)') nxa_el(kk)%z,jj,'  ',nx_dion(:,jj,kk)
       enddo
    enddo
    close(lu)
  end subroutine dump
  function quoted(l) result(q)
    character(*),intent(in) :: l
    character(2048) :: q
    integer :: a,b
    q = ''
    a = index(l,"'")
    if(a == 0) return
    b = index(l(a+1:),"'")
    if(b == 0) return
    q = l(a+1:a+b-1)
  end function quoted
  real*8 function value(l)
    character(*),intent(in) :: l
    integer :: a,io
    character(64) :: t
    value = 0d0
    a = index(l,'=')
    if(a == 0) return
    t = adjustl(l(a+1:))
    a = index(t,'d')
    if(a > 0) t(a:a) = 'e'
    read(t,*,iostat=io) value
  end function value
  real*8 function interp(nu)
    real*8,intent(in) :: nu
    integer :: lo,hi,mid
    if(nu <= nuc(1)) then
       interp = jc(1)
       return
    endif
    if(nu >= nuc(n)) then
       interp = jc(n)
       return
    endif
    lo = 1
    hi = n
    do while(hi - lo > 1)
       mid = (lo + hi)/2
       if(nuc(mid) > nu) then
          hi = mid
       else
          lo = mid
       endif
    enddo
    interp = jc(lo) + (jc(hi) - jc(lo))*(nu - nuc(lo))/(nuc(hi) - nuc(lo))
  end function interp
  subroutine seg(a,b,s,w)
!-- integral of J over [a,b] (trapezoid on the external points) and the width
    real*8,intent(in) :: a,b
    real*8,intent(inout) :: s,w
    real*8 :: xa,ya,xb,yb
    integer :: q
    xa = a
    ya = interp(a)
    do q=1,n
       if(nuc(q) <= a) cycle
       if(nuc(q) >= b) exit
       xb = nuc(q)
       yb = jc(q)
       s = s + 0.5d0*(ya + yb)*(xb - xa)
       xa = xb
       ya = yb
    enddo
    yb = interp(b)
    s = s + 0.5d0*(ya + yb)*(b - xa)
    w = b - a
  end subroutine seg
end program cell_jtest
