!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nltentmod
!-----------------------------------------------------------------------
! Non-thermal ionization and excitation of the NLTE species by fast
! electrons (gamma-ray Compton electrons, positrons, X-ray photo- and Auger
! electrons), from the Spencer-Fano degradation of nonthermmod with all
! NLTE ions as targets, as in CMFGEN (Li, Dessart & Hillier 2012):
!  - ionization: the subshells of the ground configuration of each NLTE
!    ion, with the fits of Arnaud & Rothenflug (1985) from the file
!    nt_ar.dat of the model-atom directory (Tools/cmfgen2superlite_nt.py,
!    from CMFGEN's arnaud_rothenflug data) or else Lotz (1967) with the
!    subshell thresholds of Verner et al. (1996); the Opal width of the
!    secondaries 0.6 I.  A route acts on the levels of the ground term of
!    its ion and leads to level 1 of the next ion;
!  - excitation: from every level that holds more than nxn_omlev of its
!    element, in ions that hold more than nxn_omion of it, to every level
!    above with a line or a tabulated collision strength: by van
!    Regemorter where the line is allowed (A > 1e3 /s) or has no
!    collision strength, else with the collision strength of the table's
!    highest temperature, constant in energy;
!  - ionization routes of ions below nxn_omion of the most populated
!    target are left out of the degradation (at least the ten most
!    populated ions are kept); these trace ions are ionized by the
!    degradation spectrum of the others (their rate per particle; in the
!    IME layers of a SN Ia all of Fe, Co and Ni);
!  - one electron injected at nxn_emax (1000 eV) on a grid of nxn_de.
! nxn_rates solves a cell from its populations and gives the
! ionization rate per particle of each state, the excitation rates between
! states (nxn_xl, nxn_xu, nxn_xr; nltexmod nx_el_matrix) and the heat
! share (nxn_fhc; tebalmod te_heatshare), per unit deposited power, and
! stores the cumulative energy shares of the channels into which the
! Monte Carlo deposits a fast electron (nxn_pick): heat, the continuum
! state of each ionized ion, the upper state of each excitation.
!-----------------------------------------------------------------------
  use nlteatommod
  use nonthermmod
  implicit none
  private
  real*8,parameter :: evg = 1.602176634d-12, kbol = 1.380649d-16
!-- switches (inputparmod in_nt_sf, in_nt_sfemax, in_nt_sfde)
  logical,public :: nxn_on = .false.
  real*8,public :: nxn_emax = 1000d0           !injection energy, top of the grid [eV]
  real*8,public :: nxn_de = 1d0                !grid step [eV]
  real*8,public :: nxn_omion = 1d-3            !smallest share of an ion (see above)
  real*8,public :: nxn_omlev = 1d-4            !smallest share of a lower level of an excitation
!-- ionization routes
  integer,public :: nxn_nrt = 0
  integer,allocatable :: rt_ion(:)             !(nrt) NLTE ion
  real*8,allocatable :: rt_thr(:),rt_j(:)      !(nrt) threshold, Opal width [eV]
  real*8,allocatable :: rt_sig(:,:)            !(nt_ne,nrt) cross-section [cm^2]
  logical,allocatable :: rt_ar(:)              !(nrt) Arnaud & Rothenflug (else Lotz)
  logical,allocatable :: lv_gt(:)              !(nlv) level of the ground term of its ion
!-- excitation candidates in the order of the global lower level
  integer,public :: nxn_nxc = 0
  integer,allocatable :: xc_l(:),xc_u(:),xc_k(:) !(nxc) global lower, upper level; kind
  real*8,allocatable :: xc_e(:),xc_p(:)        !(nxc) energy [eV]; f or Omega/g_l
  integer,allocatable :: lv_x0(:),lv_x1(:)     !(nlv) candidates of a lower level
!-- global state -> element, element-local state; continuum state of an ion
  integer,allocatable :: st_el(:),st_s(:)      !(nst)
!-- the current cell: excitation rates per particle of the lower state
!-- [1/s] between global states, and the heat share
  integer,public :: nxn_xn = 0
  integer,allocatable,public :: nxn_xl(:),nxn_xu(:)
  real*8,allocatable,public :: nxn_xr(:)
  real*8,public :: nxn_fhc = 1d0
!-- cumulative energy shares of the deposition channels of each cell: heat
!-- (0) and the activation of each state (global index)
  real*8,allocatable,public :: nxn_qloc(:,:)   !(0:nst,gas_ncell)
  real*8,allocatable,public :: nxn_q(:,:)      !(0:nst,ncell)
!-- (diagnostics: per cell the heat share and the shares of ionization and
!-- excitation)
  real*8,allocatable,public :: nxn_sloc(:,:)   !(3,gas_ncell)
!
  public :: nxn_init,nxn_rates,nxn_done,nxn_pick
!
contains
!
!
subroutine nxn_init(dir,ncell,gncell,lmpi0)
!-----------------------------------------------------------------------
! The routes and the excitation candidates of the loaded model atoms.
!-----------------------------------------------------------------------
  use bfxsmod, only:bf_ph1,bfxs_read_data
  character(*),intent(in) :: dir
  integer,intent(in) :: ncell,gncell
  logical,intent(in) :: lmpi0
  integer,parameter :: mr = 4000
  integer :: i,k,u,kl,c,iu,ios,z,ne,n,l,ns,nar,nlz,nl,lv,iel,s,j
  real*8 :: thr,a,b,cc,d,om,e1,qocc(7),qa
  integer :: arz(mr),arne(mr)
  real*8 :: arp(5,mr)
  character(200) :: line
  integer,allocatable :: cmap(:,:)
  integer,allocatable :: tl(:),tu(:),tk(:)
  real*8,allocatable :: te(:),tp(:)
  logical :: lany
!
  call nt_grid(nxn_emax,nxn_de)
!-- levels of the ground term: the name of level 1 up to its [J]
  allocate(lv_gt(nxa_nlv))
  lv_gt = .false.
  do i=1,nxa_nion
     do k=1,nxa_ion(i)%nlev
        lv_gt(nxa_ion(i)%lvoff+k) = term(nxa_ion(i)%lname(k)) == term(nxa_ion(i)%lname(1))
     enddo
  enddo
!-- Arnaud & Rothenflug subshells: Z, electrons, n, l, I [eV], A, B, C, D
  nar = 0
  open(newunit=iu,file=trim(dir)//'/nt_ar.dat',status='old',action='read',iostat=ios)
  if(ios==0) then
     do
        read(iu,'(a)',iostat=ios) line
        if(ios/=0) exit
        if(len_trim(line)==0 .or. line(1:1)=='#') cycle
        read(line,*,iostat=ios) z,ne,n,l,thr,a,b,cc,d
        if(ios/=0) cycle
        if(nar==mr) stop 'nxn_init: nt_ar.dat too long'
        nar = nar + 1
        arz(nar) = z
        arne(nar) = ne
        arp(:,nar) = [thr,a,b,cc,d]
     enddo
     close(iu)
  endif
!-- (the subshells of Verner et al. for the Lotz cross-sections of the ions
!-- without fits)
  lany = .false.
  do i=1,nxa_nion
     ne = nxa_ion(i)%z - nxa_ion(i)%stage + 1
     if(count(arz(:nar)==nxa_ion(i)%z .and. arne(:nar)==ne) == 0) lany = .true.
  enddo
  if(lany .and. all(bf_ph1(1,:,8,8) == 0.)) call bfxs_read_data
!-- routes
  nxn_nrt = 0
  do i=1,nxa_nion
     ne = nxa_ion(i)%z - nxa_ion(i)%stage + 1
     n = count(arz(:nar)==nxa_ion(i)%z .and. arne(:nar)==ne)
     if(n==0) n = 7
     nxn_nrt = nxn_nrt + max(n,1)
  enddo
  allocate(rt_ion(nxn_nrt),rt_thr(nxn_nrt),rt_j(nxn_nrt),rt_ar(nxn_nrt),rt_sig(nt_ne,nxn_nrt))
  nxn_nrt = 0
  nlz = 0
  do i=1,nxa_nion
     z = nxa_ion(i)%z
     ne = z - nxa_ion(i)%stage + 1
     lany = .false.
     do k=1,nar
        if(arz(k)/=z .or. arne(k)/=ne) cycle
        nxn_nrt = nxn_nrt + 1
        rt_ion(nxn_nrt) = i
        rt_thr(nxn_nrt) = arp(1,k)
        rt_ar(nxn_nrt) = .true.
        call nt_arfit(arp(1,k),arp(2,k),arp(3,k),arp(4,k),arp(5,k),rt_sig(:,nxn_nrt))
        lany = .true.
     enddo
     if(lany .or. z>30) cycle
!-- Lotz with the subshells of Verner et al.: the outermost at the
!-- ionization energy of the model atom
     nlz = nlz + 1
     call occup(z,ne,qocc)
     e1 = nxa_ion(i)%chi/evg
     qa = 1d99
     do ns=1,7
        if(bf_ph1(1,ns,z,ne) > 0. .and. qocc(ns) > 0d0) qa = min(qa,dble(bf_ph1(1,ns,z,ne)))
     enddo
     do ns=1,7
        if(qocc(ns) <= 0d0) cycle
        thr = bf_ph1(1,ns,z,ne)
        if(thr <= 0d0) then
           if(ns /= 7 .and. ns /= 6) cycle
           thr = e1
        endif
        if(thr <= qa*(1d0+1d-6)) thr = e1
        nxn_nrt = nxn_nrt + 1
        rt_ion(nxn_nrt) = i
        rt_thr(nxn_nrt) = thr
        rt_ar(nxn_nrt) = .false.
        call nt_lotz(1,[thr],[qocc(ns)],rt_sig(:,nxn_nrt))
     enddo
  enddo
  rt_j = 0.6d0*rt_thr
!-- excitation candidates
  n = 0
  do i=1,nxa_nion
     n = n + nxa_ion(i)%nlin + nxa_ion(i)%ncol
  enddo
  allocate(tl(n),tu(n),tk(n),te(n),tp(n))
  allocate(lv_x0(nxa_nlv),lv_x1(nxa_nlv))
  lv_x0 = 1
  lv_x1 = 0
  nl = 0
  do i=1,nxa_nion
     allocate(cmap(nxa_ion(i)%nlev,nxa_ion(i)%nlev))
     cmap = 0
     do c=1,nxa_ion(i)%ncol
        k = nxa_ion(i)%cl(c)
        u = nxa_ion(i)%cu(c)
        if(k<1 .or. u<1) cycle
        cmap(min(k,u),max(k,u)) = c
     enddo
     do k=1,nxa_ion(i)%nlev-1
        lv = nxa_ion(i)%lvoff + k
        lv_x0(lv) = nl + 1
        do u=k+1,nxa_ion(i)%nlev
           kl = nxa_ion(i)%lin(k,u)
           c = cmap(k,u)
           om = 0d0
           if(c>0) om = nxa_ion(i)%com(nxa_ion(i)%ncolt,c)
           e1 = (nxa_ion(i)%e(u) - nxa_ion(i)%e(k))/evg
           if(e1 <= 0d0) cycle
           if(kl>0 .and. (om<=0d0 .or. nxa_ion(i)%a(max(kl,1))>1d3)) then
              if(nxa_ion(i)%f(kl) <= 0d0) cycle
              nl = nl + 1
              tk(nl) = nt_kion
              if(nxa_ion(i)%stage==1) tk(nl) = nt_kneu
              tp(nl) = nxa_ion(i)%f(kl)
           elseif(om > 0d0) then
              nl = nl + 1
              tk(nl) = nt_kom
              tp(nl) = om/nxa_ion(i)%g(k)
           else
              cycle
           endif
           tl(nl) = lv
           tu(nl) = nxa_ion(i)%lvoff + u
           te(nl) = e1
        enddo
        lv_x1(lv) = nl
     enddo
     deallocate(cmap)
  enddo
  nxn_nxc = nl
  allocate(xc_l(nl),xc_u(nl),xc_k(nl),xc_e(nl),xc_p(nl))
  xc_l = tl(:nl)
  xc_u = tu(:nl)
  xc_k = tk(:nl)
  xc_e = te(:nl)
  xc_p = tp(:nl)
  deallocate(tl,tu,tk,te,tp)
!-- states
  allocate(st_el(nxa_nst),st_s(nxa_nst))
  do iel=1,nxa_nel
     do s=1,nxa_el(iel)%ns
        st_el(nxa_el(iel)%xoff+s) = iel
        st_s(nxa_el(iel)%xoff+s) = s
     enddo
  enddo
  allocate(nxn_xl(nl),nxn_xu(nl),nxn_xr(nl))
  allocate(nxn_qloc(0:nxa_nst,gncell),nxn_q(0:nxa_nst,ncell),nxn_sloc(3,gncell))
  nxn_qloc = 1d0
  nxn_q = 1d0
  nxn_sloc = 0d0
  nxn_xn = 0
  if(lmpi0) then
     write(6,'(1x,a,f7.0,a,f5.2,a)') 'non-thermal electrons (Spencer-Fano): injected at',nxn_emax, &
          ' eV, grid',nxn_de,' eV'
     write(6,'(1x,a,i4,a,i4,a,i3,a,i7)') '  ionization routes',nxn_nrt,' (Arnaud & Rothenflug',count(rt_ar), &
          '; Lotz for',nlz,' ions), excitation candidates',nxn_nxc
     do j=1,nxa_nion
        if(count(rt_ion==j)==0) cycle
        write(6,'(3x,a8,a,10f8.2)') nxa_ion(j)%name,merge(' A&R ','Lotz ',any(rt_ar .and. rt_ion==j)), &
             pack(rt_thr,rt_ion==j)
     enddo
  endif
contains
  pure function term(nm) result(t)
    character(*),intent(in) :: nm
    character(len(nm)) :: t
    integer :: kk
    kk = index(nm,'[')
    t = nm
    if(kk>1) t = nm(:kk-1)
  end function term
  pure subroutine occup(zz,nn,q)
!-- subshell occupations 1s 2s 2p 3s 3p 3d 4s of the ground configuration
!-- (4s before 3d: 2 in neutral atoms, 1 in their ions, as Verner's shells)
    integer,intent(in) :: zz,nn
    real*8,intent(out) :: q(7)
    real*8,parameter :: cap(6) = [2d0,2d0,6d0,2d0,6d0,10d0]
    integer :: m
    real*8 :: r
    q = 0d0
    r = nn
    if(nn > 18 .and. bf_ph1(1,7,zz,nn) > 0.) then
       q(7) = merge(min(2d0,nn-18d0),1d0,zz==nn)
       r = r - q(7)
    endif
    do m=1,6
       q(m) = min(cap(m),r)
       r = r - q(m)
       if(r <= 0d0) exit
    enddo
  end subroutine occup
end subroutine nxn_init
!
!
subroutine nxn_rates(i,nel,x,ne,temp,edep,nti)
!-----------------------------------------------------------------------
! Fast electrons in local cell i, element densities nel [cm^-3], state
! populations x (fractions of the element), electron density ne and
! temperature temp, deposited power edep [erg/s/cm^3]: nti(s) the
! ionization rate per particle of state s [1/s]; the excitation rates
! between states (nxn_xn, nxn_xl, nxn_xu, nxn_xr), the heat share
! (nxn_fhc), the channel shares of the cell (nxn_qloc(:,i)) and its
! diagnostics (nxn_sloc(:,i)).
!-----------------------------------------------------------------------
  integer,intent(in) :: i
  real*8,intent(in) :: nel(nxa_nel),x(nxa_nst),ne,temp,edep
  real*8,intent(out) :: nti(nxa_nst)
  real*8 :: npop(nxa_nlv),nrt(nxn_nrt),nion(nxn_nrt),nst(nxa_nst),q(0:nxa_nst),ngt(nxn_nrt),npp(nxn_nrt)
  real*8 :: wsum(nxa_nst),fh,etot,nmax,nlim,w,c,ntot,nlow,fion,fexc
  real*8,allocatable :: ae(:),ap(:),an(:),aev(:)
  integer,allocatable :: ak(:),at(:)
  integer :: j,k,lv,s,iel,p,t,na,sl,su,nkeep,jj,ii
  real*8 :: ionpop(nxa_nion),elsum(nxa_nel)
  logical :: lion(nxa_nion),lexc(nxa_nion)
!
!-- level populations: Boltzmann at temp within the super-levels
  wsum = 0d0
  do j=1,nxa_nion
     iel = nxa_ion(j)%el
     do k=1,nxa_ion(j)%nlev
        s = nxa_el(iel)%xoff + nxa_ion(j)%soff + nxa_ion(j)%sl(k)
        w = nxa_ion(j)%g(k)*exp(-max(nxa_ion(j)%e(k)-nxa_ion(j)%eref(nxa_ion(j)%sl(k)),0d0)/(kbol*temp))
        npop(nxa_ion(j)%lvoff+k) = w
        wsum(s) = wsum(s) + w
     enddo
  enddo
  do iel=1,nxa_nel
     nst(nxa_el(iel)%xoff+1:nxa_el(iel)%xoff+nxa_el(iel)%ns) = &
          nel(iel)*x(nxa_el(iel)%xoff+1:nxa_el(iel)%xoff+nxa_el(iel)%ns)
  enddo
  ionpop = 0d0
  elsum = 0d0
  do j=1,nxa_nion
     iel = nxa_ion(j)%el
     do k=1,nxa_ion(j)%nlev
        lv = nxa_ion(j)%lvoff + k
        s = nxa_el(iel)%xoff + nxa_ion(j)%soff + nxa_ion(j)%sl(k)
        npop(lv) = nst(s)*npop(lv)/max(wsum(s),1d-300)
        ionpop(j) = ionpop(j) + npop(lv)
     enddo
     elsum(iel) = elsum(iel) + ionpop(j)
  enddo
!-- ionization routes: the population of the ground term; ions below
!-- nxn_omion of the largest left out, except the ten most populated
  nrt = 0d0
  do p=1,nxn_nrt
     j = rt_ion(p)
     do k=1,nxa_ion(j)%nlev
        lv = nxa_ion(j)%lvoff + k
        if(lv_gt(lv)) nrt(p) = nrt(p) + npop(lv)
     enddo
  enddo
  nmax = maxval(nrt)
  nlim = nxn_omion*nmax
  lion = .false.
  do p=1,nxn_nrt
     if(nrt(p) >= nlim .and. nrt(p) > 0d0) lion(rt_ion(p)) = .true.
  enddo
  nkeep = min(10,nxa_nion)
  do jj=1,nkeep
     c = -1d0
     ii = 0
     do p=1,nxn_nrt
        if(lion(rt_ion(p)) .or. nrt(p) <= c) cycle
        c = nrt(p)
        ii = rt_ion(p)
     enddo
     if(ii==0 .or. count(lion) >= nkeep) exit
     lion(ii) = .true.
  enddo
  nion = 0d0
  do p=1,nxn_nrt
     if(lion(rt_ion(p))) nion(p) = nrt(p)
  enddo
!-- excitations: from the levels above nxn_omlev of their element, in ions
!-- above nxn_omion of it
  do j=1,nxa_nion
     lexc(j) = ionpop(j) > nxn_omion*elsum(nxa_ion(j)%el)
  enddo
  na = 0
  do j=1,nxa_nion
     if(.not.lexc(j)) cycle
     ntot = elsum(nxa_ion(j)%el)
     do k=1,nxa_ion(j)%nlev-1
        lv = nxa_ion(j)%lvoff + k
        if(npop(lv) > nxn_omlev*ntot) na = na + max(lv_x1(lv)-lv_x0(lv)+1,0)
     enddo
  enddo
  allocate(ae(na),ap(na),an(na),aev(na),ak(na),at(na))
  na = 0
  do j=1,nxa_nion
     if(.not.lexc(j)) cycle
     ntot = elsum(nxa_ion(j)%el)
     do k=1,nxa_ion(j)%nlev-1
        lv = nxa_ion(j)%lvoff + k
        if(npop(lv) <= nxn_omlev*ntot) cycle
        do t=lv_x0(lv),lv_x1(lv)
           na = na + 1
           at(na) = t
           ae(na) = xc_e(t)
           ap(na) = xc_p(t)
           ak(na) = xc_k(t)
           an(na) = npop(lv)
        enddo
     enddo
  enddo
!
  ngt = nrt
  call nt_solve(nxn_nrt,rt_thr,rt_j,nion,rt_sig,na,ae,ak,ap,an,max(ne,1d-30),nxn_emax,fh,nrt,aev,npp)
!-- (nrt now: ionizations of each route per eV deposited; the trace routes
!-- from their rate per particle)
  do p=1,nxn_nrt
     if(nion(p) > 0d0 .or. ngt(p) <= 0d0) cycle
     nion(p) = ngt(p)
     nrt(p) = npp(p)*ngt(p)
  enddo
!
  etot = edep/evg
  nti = 0d0
  q = 0d0
  q(0) = fh
  fion = 0d0
  do p=1,nxn_nrt
     if(nrt(p) <= 0d0 .or. nion(p) <= 0d0) cycle
     j = rt_ion(p)
     iel = nxa_ion(j)%el
     do k=1,nxa_ion(j)%nlev
        lv = nxa_ion(j)%lvoff + k
        if(.not.lv_gt(lv)) cycle
        s = nxa_el(iel)%xoff + nxa_ion(j)%soff + nxa_ion(j)%sl(k)
        nti(s) = nti(s) + nrt(p)*npop(lv)/nion(p)
     enddo
!-- energy: the continuum state of the ion (level 1 of the next one)
     s = nxa_el(iel)%xoff + nxa_lvc(nxa_ion(j)%lvoff+1)
     q(s) = q(s) + rt_thr(p)*nrt(p)
     fion = fion + rt_thr(p)*nrt(p)
  enddo
  do s=1,nxa_nst
     if(nti(s) > 0d0) nti(s) = etot*nti(s)/max(nst(s),1d-300)
  enddo
  nxn_xn = 0
  fexc = 0d0
  do k=1,na
     if(aev(k) <= 0d0) cycle
     t = at(k)
     sl = nxa_el(nxa_lvel(xc_l(t)))%xoff + nxa_lvs(xc_l(t))
     su = nxa_el(nxa_lvel(xc_u(t)))%xoff + nxa_lvs(xc_u(t))
     q(su) = q(su) + xc_e(t)*aev(k)
     fexc = fexc + xc_e(t)*aev(k)
     if(sl==su) cycle
     nlow = nst(sl)
     if(nlow <= 0d0) cycle
     nxn_xn = nxn_xn + 1
     nxn_xl(nxn_xn) = sl
     nxn_xu(nxn_xn) = su
     nxn_xr(nxn_xn) = etot*aev(k)/nlow
  enddo
  nxn_fhc = fh
!-- cumulative channel shares
  do s=1,nxa_nst
     q(s) = q(s) + q(s-1)
  enddo
  nxn_qloc(:,i) = q/max(q(nxa_nst),1d-300)
  nxn_sloc(:,i) = [fh,fion,fexc]
  deallocate(ae,ap,an,aev,ak,at)
end subroutine nxn_rates
!
!
subroutine nxn_done
!-- the current cell's excitation list is no longer valid
  nxn_xn = 0
  nxn_fhc = 1d0
end subroutine nxn_done
!
!
pure subroutine nxn_pick(ic,r,iel,s)
!-- channel of a fast electron in cell ic for the uniform deviate r: heat
!-- (iel = 0) or the activation of element-local state s of element iel
  integer,intent(in) :: ic
  real*8,intent(in) :: r
  integer,intent(out) :: iel,s
  integer :: lo,hi,mid
  iel = 0
  s = 0
  if(r < nxn_q(0,ic)) return
  lo = 0
  hi = nxa_nst
!-- smallest index with q >= r
  do while(hi - lo > 1)
     mid = (lo + hi)/2
     if(nxn_q(mid,ic) >= r) then
        hi = mid
     else
        lo = mid
     endif
  enddo
  iel = st_el(hi)
  s = st_s(hi)
end subroutine nxn_pick
!
end module nltentmod
! vim: fdm=marker
