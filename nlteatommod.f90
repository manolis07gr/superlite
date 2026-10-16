!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
module nlteatommod
!-----------------------------------------------------------------------
! Model atoms of the NLTE species (statistical equilibrium: nltexmod;
! transport and macro-atom: nltextrmod).
!
! Each NLTE ion has levels (g, E), lines (f; A from f), collision
! strengths Omega(T), photoionization cross-sections to the ground state of
! the next ion (tabulated against u = nu/nu_edge, nu^-3 beyond the table)
! and two-photon decays.  Ions come from SuperLite model-atom files
! (Tools/cmfgen2superlite_atom.py converts CMFGEN atoms) or are built in
! (hydrogen: the n-resolved atom of nltehmod).
!
! Collision strengths missing from the tables follow CMFGEN (OMEGA_GEN_V3):
!  - f > 1e-5: van Regemorter, Omega = 47.972 s gbar f g_l/nu[1e15 Hz], with
!    gbar = 0.276 e^x E1(x) (x <= 14), 0.066 (1+1.5/x)/sqrt(x) (x > 14) for
!    neutrals (Auer & Mihalas 1973) and max(0.2 [0.7 for Delta n = 0],
!    0.276 e^x E1(x)) for ions; s = omegascale of the file
!  - otherwise Omega = omegaset of the file
!  - ionization: Seaton, C = 1.55e13 gbar sigma_0 T^-1/2 e^-u/u with
!    gbar = 0.1 min(zion,3)
! Collision rates: q_lu = 8.63e-6 Omega/(g_l sqrt T) e^{-dE/kT}.
!
! The ions of an element form one chain of states (element-local index):
! the levels of each NLTE ion in increasing stage, then the next ion
! ("top", one state).  State energies eps are counted from the ground state
! of the lowest NLTE ion.
!-----------------------------------------------------------------------
  use nltehmod, only:nlteh_init,nh_nlev,nh_nlin,nh_g,nh_eexc,nh_eion,nh_ll,nh_lu, &
       nh_f,nh_chi,nlteh_sigma,nlteh_qcoll,nlteh_air,a2s
  implicit none
  private
!
  real*8,parameter :: pi = 3.141592653589793d0
  real*8,parameter :: cc = 2.99792458d10, hh = 6.6260693d-27, kb = 1.3806505d-16
  real*8,parameter :: me = 9.1093826d-28, ee = 4.80320441d-10, ev = 1.60217653d-12
  real*8,parameter,public :: nxa_sahac = 4.1412d-16     !(h^2/(2 pi m_e k))^1.5 [cm^3 K^1.5]
  real*8,parameter,public :: nxa_fsob = pi*ee*ee/(me*cc) !0.02654 cm^2/s
  integer,parameter,public :: nxa_maxion = 8             !NLTE ions per element
!
  type,public :: nxion
     integer :: z = 0, stage = 0                 !atomic number, 1 = neutral
     character(8) :: name = ''
     real*8 :: zion = 1d0                        !charge seen by the outer electron
     real*8 :: chi = 0d0                         !ionization energy of level 1 [erg]
     real*8 :: gion = 1d0                        !g of the ground state of the next ion
     integer :: nlev = 0, nlin = 0, ntwo = 0, ncol = 0, ncolt = 0
     real*8,allocatable :: g(:),e(:),nu0(:)      !(nlev) weight, energy above level 1 [erg], edge [Hz]
     character(24),allocatable :: lname(:)
     integer,allocatable :: ll(:),lu(:)          !(nlin) lower, upper level
     real*8,allocatable :: f(:),a(:),nu(:)       !(nlin)
     logical,allocatable :: samen(:)             !(nlin) same n (CMFGEN SAME_N)
     integer,allocatable :: lin(:,:)             !(nlev,nlev) line of the pair, 0 if none
     integer,allocatable :: tl(:),tu(:)          !(ntwo) two-photon lower, upper level
     real*8,allocatable :: ta(:)                 !(ntwo) rate per atom in the upper level [1/s]
     real*8,allocatable :: colt(:)               !(ncolt) temperatures of the Omega table [K]
     integer,allocatable :: cl(:),cu(:)          !(ncol) pair; cu = 0: ionization of cl
     real*8,allocatable :: com(:,:)              !(ncolt,ncol) Omega
     real*8 :: omset = 0.1d0, omscale = 1d0
     integer,allocatable :: pn(:),pst(:)         !(nlev) cross-section table: points, start
     real*8,allocatable :: pu(:),ps(:)           !u = nu/nu0, sigma [cm^2]
     integer :: el = 0                           !element index
     integer :: loff = 0                         !global line index offset
     integer :: lvoff = 0                        !global level index offset
     integer :: soff = 0                         !element-local state offset of level 1
  end type nxion
!
  type,public :: nxelem
     integer :: z = 0, nion = 0
     integer :: ion(nxa_maxion) = 0              !NLTE ions, increasing stage
     integer :: ns = 0                           !states: all levels and the top ion
     integer :: xoff = 0                         !offset in the per-cell state arrays
     real*8 :: gtop = 1d0                        !weight of the top ion
     real*8,allocatable :: eps(:)                !(ns) energy of the states [erg]
     real*8,allocatable :: chg(:)                !(ns) ion charge of the states
  end type nxelem
!
  real*8,public :: nxa_colfac = 1d0              !scale of all electron-impact rates (tests)
  integer,public :: nxa_nion = 0, nxa_nel = 0
  type(nxion),allocatable,public :: nxa_ion(:)
  type(nxelem),allocatable,public :: nxa_el(:)
  integer,public :: nxa_nst = 0                  !states of all elements
!-- global lines (all ions)
  integer,public :: nxa_nlin = 0
  integer,allocatable,public :: nxa_lion(:)      !(nlin) ion
  integer,allocatable,public :: nxa_lel(:)       !(nlin) element
  integer,allocatable,public :: nxa_lsl(:),nxa_lsu(:) !(nlin) element-local lower, upper state
  integer,allocatable,public :: nxa_lloc(:)      !(nlin) line index within the ion
  real*8,allocatable,public :: nxa_lnu(:),nxa_lwl(:)  !(nlin) frequency [Hz], wavelength [cm]
  real*8,allocatable,public :: nxa_lf(:),nxa_la(:),nxa_lblu(:),nxa_lbul(:) !(nlin)
  integer,allocatable,public :: nxa_lsrt(:)      !(nlin) lines sorted by wavelength
!-- global levels (all ions)
  integer,public :: nxa_nlv = 0
  integer,allocatable,public :: nxa_lvion(:),nxa_lvel(:) !(nlv) ion, element
  integer,allocatable,public :: nxa_lvs(:),nxa_lvc(:)    !(nlv) element-local state, its continuum state
  real*8,allocatable,public :: nxa_lvnu0(:)               !(nlv) threshold [Hz]
  real*8,allocatable,public :: nxa_lvphi(:)               !(nlv) g/(2 g_+)
!-- group integrals of the cross-sections (set by nxa_groups)
  integer,public :: nxa_ng = 0
  real*8,allocatable,public :: nxa_sgm(:,:)   !(ng,nlv) group mean of sigma [cm^2]
  real*8,allocatable,public :: nxa_sgh(:,:)   !(ng,nlv) integral of sigma/(h nu) dnu over the group
  integer,allocatable,public :: nxa_g1(:)     !(nlv) reddest group above the edge (groups 1..g1)
  integer,allocatable,public :: nxa_lgn(:)    !(nlin) group with the centre nearest the line
!-- recombination tables on a logarithmic temperature grid
  integer,parameter,public :: nxa_ntt = 81
  real*8,public :: nxa_ltt(nxa_ntt)           !ln T
  real*8,allocatable,public :: nxa_lalp(:,:)  !(ntt,nlv) ln alpha [cm^3/s]
  real*8,allocatable,public :: nxa_lrc(:,:)   !(ntt,nlv) ln recombination cooling [erg cm^3/s]
!
!-- charge exchange  A + B -> C + D:  A = X^(n+1) and C = X^n of element X,
!-- B = Y^m and D = Y^(m+1) of element Y (CMFGEN data, Kingdon & Ferland
!-- 1996 fits); participant k = 1..4 (A, B, C, D) is a term of an NLTE ion
!-- (its levels, element-local states s) or the top state.  The reverse
!-- rate follows from detailed balance with the model levels.
  integer,parameter,public :: nxa_cxml = 16
  type,public :: nxcx
     integer :: el(4) = 0, ion(4) = 0          !element; NLTE ion (0: top state)
     integer :: nl(4) = 0, s(nxa_cxml,4) = 0   !states of the term
     integer :: typ = 2
     real*8 :: c(4) = 0d0, tlo = 0d0, thi = 1d2 !fit coefficients; T range [1e4 K]
     real*8 :: q = 0d0                          !energy released (forward), lowest levels [erg]
  end type nxcx
  integer,public :: nxa_ncx = 0
  type(nxcx),allocatable,public :: nxa_cx(:)
!
  public :: nxa_read,nxa_builtin_h,nxa_setup,nxa_groups,nxa_sigma,nxa_alpha,nxa_rcool
  public :: nxa_coll,nxa_air,nxa_read_cx,nxa_cx_coef
!
contains
!
!
pure function nxa_air(wlvac) result(wl)
!-- vacuum to air wavelength [cm] above 2000 A
  real*8,intent(in) :: wlvac
  real*8 :: wl
  wl = nlteh_air(wlvac)
end function nxa_air
!
!
subroutine getrec(iu,s,ios)
!-- next record that is not blank or a comment; text after '!' removed
  integer,intent(in) :: iu
  character(*),intent(out) :: s
  integer,intent(out) :: ios
  integer :: k
  do
     read(iu,'(a)',iostat=ios) s
     if(ios/=0) return
     k = index(s,'!')
     if(k>0) s(k:) = ' '
     s = adjustl(s)
     if(len_trim(s)==0) cycle
     if(s(1:1)=='#') cycle
     return
  enddo
end subroutine getrec
!
!
subroutine keyval(iu,key,s)
!-- next record must start with key; s returns the rest
  integer,intent(in) :: iu
  character(*),intent(in) :: key
  character(*),intent(out) :: s
  character(512) :: r
  integer :: ios
  call getrec(iu,r,ios)
  if(ios/=0) then
     write(0,*) 'nxa_read: end of file, expected ',key
     stop 'nxa_read: bad model atom'
  endif
  if(r(1:len(key))/=key) then
     write(0,*) 'nxa_read: expected ',key,', found: ',trim(r)
     stop 'nxa_read: bad model atom'
  endif
  s = r(len(key)+1:)
end subroutine keyval
!
!
subroutine nxa_read(fname,ion)
!-----------------------------------------------------------------------
! Read a SuperLite model-atom file (format 1).
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname
  type(nxion),intent(out) :: ion
  integer,parameter :: iu = 77
  character(512) :: s
  character(8) :: sym
  integer :: ios,i,k,l,u,isn,n,np,ist
  real*8 :: chi,lam,e,g,fval
!
  open(iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) then
     write(0,*) 'nxa_read: cannot open ',trim(fname)
     stop 'nxa_read: model atom missing'
  endif
  call keyval(iu,'element',s)
  read(s,*) ion%z,sym
  call keyval(iu,'ion',s)
  read(s,*) ion%stage
  write(ion%name,'(a,a)') trim(sym),trim(roman(ion%stage))
  call keyval(iu,'zion',s)
  read(s,*) ion%zion
  call keyval(iu,'chi',s)
  read(s,*) chi
  ion%chi = hh*cc*chi
  call keyval(iu,'gion',s)
  read(s,*) ion%gion
  call keyval(iu,'nlev',s)
  read(s,*) n
  ion%nlev = n
  allocate(ion%g(n),ion%e(n),ion%nu0(n),ion%lname(n))
  do i=1,n
     call getrec(iu,s,ios)
!-- index, name (may contain '/', which ends list-directed input), g, E
     s = adjustl(s)
     k = index(s,' ')
     s = adjustl(s(k:))
     k = index(s,' ')
     ion%lname(i) = s(:k-1)
     read(s(k:),*) g,e
     ion%g(i) = g
     ion%e(i) = hh*cc*e
  enddo
!-- energies above level 1
  ion%e = ion%e - ion%e(1)
  ion%nu0 = (ion%chi - ion%e)/hh
  if(any(ion%nu0<=0d0)) stop 'nxa_read: level above the ionization limit'
!-- lines
  call keyval(iu,'nline',s)
  read(s,*) n
  ion%nlin = n
  allocate(ion%ll(n),ion%lu(n),ion%f(n),ion%a(n),ion%nu(n),ion%samen(n))
  allocate(ion%lin(ion%nlev,ion%nlev))
  ion%lin = 0
  do k=1,n
     call getrec(iu,s,ios)
     read(s,*) l,u,fval,e,lam,isn
     ion%ll(k) = l
     ion%lu(k) = u
     ion%f(k) = fval
     ion%samen(k) = isn/=0
     ion%nu(k) = (ion%e(u) - ion%e(l))/hh
     if(ion%nu(k)<=0d0) then
        write(0,*) 'nxa_read: ',trim(fname),' line',k,l,u
        stop 'nxa_read: line with non-positive frequency'
     endif
     ion%lin(l,u) = k
     ion%lin(u,l) = k
  enddo
!-- A from f (consistent absorption and emission)
  do k=1,n
     ion%a(k) = 8d0*pi**2*ee**2*ion%nu(k)**2/(me*cc**3)*ion%g(ion%ll(k))/ion%g(ion%lu(k))*ion%f(k)
  enddo
!-- two-photon
  call keyval(iu,'ntwo',s)
  read(s,*) n
  ion%ntwo = n
  allocate(ion%tl(n),ion%tu(n),ion%ta(n))
  do k=1,n
     call getrec(iu,s,ios)
     read(s,*) ion%tl(k),ion%tu(k),ion%ta(k)
  enddo
!-- collisions
  call keyval(iu,'ncolt',s)
  read(s,*) n
  ion%ncolt = n
  allocate(ion%colt(n))
  call keyval(iu,'colt',s)
  read(s,*) ion%colt
  ion%colt = 1d4*ion%colt
  call keyval(iu,'omegaset',s)
  read(s,*) ion%omset
  call keyval(iu,'omegascale',s)
  read(s,*) ion%omscale
  call keyval(iu,'ncol',s)
  read(s,*) n
  ion%ncol = n
  allocate(ion%cl(n),ion%cu(n),ion%com(ion%ncolt,n))
  do k=1,n
     call getrec(iu,s,ios)
     read(s,*) ion%cl(k),ion%cu(k),ion%com(:,k)
  enddo
!-- photoionization
  call keyval(iu,'nphot',s)
  read(s,*) n
  if(n/=ion%nlev) stop 'nxa_read: nphot /= nlev'
  allocate(ion%pn(n),ion%pst(n))
  np = 0
  do i=1,n
     call getrec(iu,s,ios)
     read(s,*) k,ion%pn(i)
     ion%pst(i) = np + 1
     np = np + ion%pn(i)
     do k=1,ion%pn(i)
        call getrec(iu,s,ios)
     enddo
  enddo
  close(iu)
!-- second pass for the table (sizes known now)
  allocate(ion%pu(np),ion%ps(np))
  open(iu,file=fname,status='old',action='read')
  do
     call getrec(iu,s,ios)
     if(s(1:5)=='nphot') exit
  enddo
  ist = 0
  do i=1,n
     call getrec(iu,s,ios)
     do k=1,ion%pn(i)
        call getrec(iu,s,ios)
        ist = ist + 1
        read(s,*) ion%pu(ist),ion%ps(ist)
     enddo
  enddo
  close(iu)
contains
  function roman(k) result(r)
    integer,intent(in) :: k
    character(5) :: r
    character(5),parameter :: rn(12) = [character(5) :: 'I','II','III','IV','V','VI', &
         'VII','VIII','IX','X','XI','XII']
    r = rn(min(max(k,1),12))
  end function roman
end subroutine nxa_read
!
!
subroutine nxa_builtin_h(nlev,ion)
!-----------------------------------------------------------------------
! Hydrogen with n = 1..nlev from nltehmod (Johnson f-values and Gaunt
! factors, l-bundled levels).  The collision rates of nltehmod (fits of
! data.nlte.h1 if loaded, else van Regemorter) are tabulated as Omega(T) for
! all pairs; ionization follows Seaton (the nltehmod rule).
!-----------------------------------------------------------------------
  integer,intent(in) :: nlev
  type(nxion),intent(out) :: ion
  integer,parameter :: nt = 41, npu = 161
  integer :: n,k,l,u,it,ic,ip
  real*8 :: t,qup(nlev,nlev),qdn(nlev,nlev)
!
  if(nh_nlev/=nlev) call nlteh_init(nlev)
  ion%z = 1
  ion%stage = 1
  ion%name = 'HI'
  ion%zion = 1d0
  ion%chi = nh_chi
  ion%gion = 1d0
  ion%nlev = nlev
  allocate(ion%g(nlev),ion%e(nlev),ion%nu0(nlev),ion%lname(nlev))
  do n=1,nlev
     ion%g(n) = nh_g(n)
     ion%e(n) = nh_eexc(n)
     ion%nu0(n) = nh_eion(n)/hh
     write(ion%lname(n),'(i0,a)') n,'___'
  enddo
  ion%nlin = nh_nlin
  allocate(ion%ll(nh_nlin),ion%lu(nh_nlin),ion%f(nh_nlin),ion%a(nh_nlin),ion%nu(nh_nlin))
  allocate(ion%samen(nh_nlin),ion%lin(nlev,nlev))
  ion%lin = 0
  do k=1,nh_nlin
     ion%ll(k) = nh_ll(k)
     ion%lu(k) = nh_lu(k)
     ion%f(k) = nh_f(k)
     ion%samen(k) = .false.
     ion%nu(k) = (nh_eion(nh_ll(k)) - nh_eion(nh_lu(k)))/hh
     ion%a(k) = 8d0*pi**2*ee**2*ion%nu(k)**2/(me*cc**3)*nh_g(nh_ll(k))/nh_g(nh_lu(k))*nh_f(k)
     ion%lin(nh_ll(k),nh_lu(k)) = k
     ion%lin(nh_lu(k),nh_ll(k)) = k
  enddo
!-- two-photon 2s-1s, l-mixed n = 2
  ion%ntwo = 1
  allocate(ion%tl(1),ion%tu(1),ion%ta(1))
  ion%tl(1) = 1
  ion%tu(1) = 2
  ion%ta(1) = 0.25d0*a2s
!-- collision strengths of all pairs on 1e3..1e6 K
  ion%ncolt = nt
  allocate(ion%colt(nt))
  do it=1,nt
     ion%colt(it) = 10d0**(3d0 + 3d0*(it-1)/(nt-1d0))
  enddo
  ion%ncol = nlev*(nlev-1)/2
  allocate(ion%cl(ion%ncol),ion%cu(ion%ncol),ion%com(nt,ion%ncol))
  ic = 0
  do l=1,nlev-1
     do u=l+1,nlev
        ic = ic + 1
        ion%cl(ic) = l
        ion%cu(ic) = u
     enddo
  enddo
  do it=1,nt
     t = ion%colt(it)
     call nlteh_qcoll(t,qup,qdn)
     do ic=1,ion%ncol
        ion%com(it,ic) = qdn(ion%cu(ic),ion%cl(ic))*ion%g(ion%cu(ic))*sqrt(t)/8.63d-6
     enddo
  enddo
  ion%omset = 0d0
  ion%omscale = 1d0
!-- cross-sections on u = 1..1e3
  allocate(ion%pn(nlev),ion%pst(nlev),ion%pu(nlev*npu),ion%ps(nlev*npu))
  ip = 0
  do n=1,nlev
     ion%pn(n) = npu
     ion%pst(n) = ip + 1
     do k=1,npu
        ip = ip + 1
        ion%pu(ip) = 10d0**(3d0*(k-1)/(npu-1d0))
        ion%ps(ip) = nlteh_sigma(n,ion%pu(ip)*ion%nu0(n))
     enddo
  enddo
end subroutine nxa_builtin_h
!
!
subroutine nxa_setup(nion,ions)
!-----------------------------------------------------------------------
! Group the ions into elements, set the state indices, energies, global
! line and level lists, and the recombination tables.
!-----------------------------------------------------------------------
  integer,intent(in) :: nion
  type(nxion),intent(inout) :: ions(nion)
  integer :: i,j,k,iel,n,m,lv,it,ks
  integer :: order(nion)
  real*8 :: e0
  integer,allocatable :: idx(:)
!
  nxa_nion = nion
  if(allocated(nxa_ion)) deallocate(nxa_ion)
  allocate(nxa_ion(nion))
!-- set up again: drop the previous lists
  if(allocated(nxa_lion)) deallocate(nxa_lion,nxa_lel,nxa_lsl,nxa_lsu,nxa_lloc,nxa_lnu, &
       nxa_lwl,nxa_lf,nxa_la,nxa_lblu,nxa_lbul,nxa_lsrt)
  if(allocated(nxa_lvion)) deallocate(nxa_lvion,nxa_lvel,nxa_lvs,nxa_lvc,nxa_lvnu0,nxa_lvphi)
  if(allocated(nxa_lalp)) deallocate(nxa_lalp,nxa_lrc)
!-- sort by (z, stage)
  order = [(i,i=1,nion)]
  do i=2,nion
     k = order(i)
     j = i - 1
     do while(j>=1)
        if(ions(order(j))%z*100+ions(order(j))%stage <= ions(k)%z*100+ions(k)%stage) exit
        order(j+1) = order(j)
        j = j - 1
     enddo
     order(j+1) = k
  enddo
  do i=1,nion
     nxa_ion(i) = ions(order(i))
  enddo
!-- elements
  nxa_nel = 0
  do i=1,nion
     if(i==1) then
        nxa_nel = 1
     elseif(nxa_ion(i)%z/=nxa_ion(i-1)%z) then
        nxa_nel = nxa_nel + 1
     elseif(nxa_ion(i)%stage/=nxa_ion(i-1)%stage+1) then
        write(0,*) 'nxa_setup: ions of Z =',nxa_ion(i)%z,' not consecutive'
        stop 'nxa_setup: missing ion stage'
     endif
  enddo
  if(allocated(nxa_el)) deallocate(nxa_el)
  allocate(nxa_el(nxa_nel))
  iel = 0
  do i=1,nion
     if(i==1) then
        iel = 1
     elseif(nxa_ion(i)%z/=nxa_ion(i-1)%z) then
        iel = iel + 1
     endif
     nxa_el(iel)%z = nxa_ion(i)%z
     nxa_el(iel)%nion = nxa_el(iel)%nion + 1
     if(nxa_el(iel)%nion>nxa_maxion) stop 'nxa_setup: too many ions per element'
     nxa_el(iel)%ion(nxa_el(iel)%nion) = i
     nxa_ion(i)%el = iel
  enddo
!-- states and energies
  nxa_nst = 0
  do iel=1,nxa_nel
     n = 0
     do j=1,nxa_el(iel)%nion
        i = nxa_el(iel)%ion(j)
        nxa_ion(i)%soff = n
        n = n + nxa_ion(i)%nlev
     enddo
     nxa_el(iel)%ns = n + 1
     nxa_el(iel)%xoff = nxa_nst
     nxa_nst = nxa_nst + nxa_el(iel)%ns
     i = nxa_el(iel)%ion(nxa_el(iel)%nion)
     nxa_el(iel)%gtop = nxa_ion(i)%gion
     allocate(nxa_el(iel)%eps(n+1),nxa_el(iel)%chg(n+1))
     e0 = 0d0
     do j=1,nxa_el(iel)%nion
        i = nxa_el(iel)%ion(j)
        do k=1,nxa_ion(i)%nlev
           nxa_el(iel)%eps(nxa_ion(i)%soff+k) = e0 + nxa_ion(i)%e(k)
           nxa_el(iel)%chg(nxa_ion(i)%soff+k) = dble(nxa_ion(i)%stage - 1)
        enddo
        e0 = e0 + nxa_ion(i)%chi
     enddo
     nxa_el(iel)%eps(n+1) = e0
     nxa_el(iel)%chg(n+1) = dble(nxa_ion(i)%stage)
!-- recombination from the next NLTE ion comes from its level 1.  The
!-- cross-sections of CMFGEN lead to the ground term of the next ion (gion),
!-- which may be split in J: with its fine-structure levels thermalized
!-- relative to each other (dense gas), n_1/g_1 = n_term/g_term, so the
!-- rates are those of the term with g_1 for gion.
     do j=1,nxa_el(iel)%nion-1
        i = nxa_el(iel)%ion(j)
        k = nxa_el(iel)%ion(j+1)
        nxa_ion(i)%gion = nxa_ion(k)%g(1)
     enddo
  enddo
!-- global lines
  nxa_nlin = sum(nxa_ion(:)%nlin)
  allocate(nxa_lion(nxa_nlin),nxa_lel(nxa_nlin),nxa_lsl(nxa_nlin),nxa_lsu(nxa_nlin))
  allocate(nxa_lloc(nxa_nlin),nxa_lnu(nxa_nlin),nxa_lwl(nxa_nlin),nxa_lf(nxa_nlin))
  allocate(nxa_la(nxa_nlin),nxa_lblu(nxa_nlin),nxa_lbul(nxa_nlin),nxa_lsrt(nxa_nlin))
  m = 0
  do i=1,nion
     nxa_ion(i)%loff = m
     do k=1,nxa_ion(i)%nlin
        m = m + 1
        nxa_lion(m) = i
        nxa_lel(m) = nxa_ion(i)%el
        nxa_lloc(m) = k
        nxa_lsl(m) = nxa_ion(i)%soff + nxa_ion(i)%ll(k)
        nxa_lsu(m) = nxa_ion(i)%soff + nxa_ion(i)%lu(k)
        nxa_lnu(m) = nxa_ion(i)%nu(k)
        nxa_lwl(m) = nxa_air(cc/nxa_lnu(m))
        nxa_lf(m) = nxa_ion(i)%f(k)
        nxa_la(m) = nxa_ion(i)%a(k)
        nxa_lbul(m) = nxa_la(m)*cc**2/(2d0*hh*nxa_lnu(m)**3)
        nxa_lblu(m) = nxa_lbul(m)*nxa_ion(i)%g(nxa_ion(i)%lu(k))/nxa_ion(i)%g(nxa_ion(i)%ll(k))
     enddo
  enddo
  allocate(idx(nxa_nlin))
  idx = [(k,k=1,nxa_nlin)]
  call sortidx(nxa_lwl,idx)
  nxa_lsrt = idx
!-- global levels
  nxa_nlv = sum(nxa_ion(:)%nlev)
  allocate(nxa_lvion(nxa_nlv),nxa_lvel(nxa_nlv),nxa_lvs(nxa_nlv),nxa_lvc(nxa_nlv))
  allocate(nxa_lvnu0(nxa_nlv),nxa_lvphi(nxa_nlv))
  lv = 0
  do iel=1,nxa_nel
     do j=1,nxa_el(iel)%nion
        i = nxa_el(iel)%ion(j)
        nxa_ion(i)%lvoff = lv
        if(j<nxa_el(iel)%nion) then
           ks = nxa_ion(nxa_el(iel)%ion(j+1))%soff + 1
        else
           ks = nxa_el(iel)%ns
        endif
        do k=1,nxa_ion(i)%nlev
           lv = lv + 1
           nxa_lvion(lv) = i
           nxa_lvel(lv) = iel
           nxa_lvs(lv) = nxa_ion(i)%soff + k
           nxa_lvc(lv) = ks
           nxa_lvnu0(lv) = nxa_ion(i)%nu0(k)
           nxa_lvphi(lv) = nxa_ion(i)%g(k)/(2d0*nxa_ion(i)%gion)
        enddo
     enddo
  enddo
!-- recombination coefficients and cooling on 1e2..1e8 K
  allocate(nxa_lalp(nxa_ntt,nxa_nlv),nxa_lrc(nxa_ntt,nxa_nlv))
  do it=1,nxa_ntt
     nxa_ltt(it) = log(1d2) + (it-1)*log(1d6)/(nxa_ntt-1)
  enddo
  do lv=1,nxa_nlv
     i = nxa_lvion(lv)
     k = lv - nxa_ion(i)%lvoff
     do it=1,nxa_ntt
        call milne(nxa_ion(i),k,exp(nxa_ltt(it)),nxa_lalp(it,lv),nxa_lrc(it,lv))
     enddo
  enddo
  nxa_lalp = log(max(nxa_lalp,1d-300))
  nxa_lrc = log(max(nxa_lrc,1d-300))
contains
  subroutine sortidx(a,ix)
!-- insertion sort of the index ix by a(ix)
    real*8,intent(in) :: a(:)
    integer,intent(inout) :: ix(:)
    integer :: ii,jj,t
    do ii=2,size(ix)
       t = ix(ii)
       jj = ii - 1
       do while(jj>=1)
          if(a(ix(jj)) <= a(t)) exit
          ix(jj+1) = ix(jj)
          jj = jj - 1
       enddo
       ix(jj+1) = t
    enddo
  end subroutine sortidx
end subroutine nxa_setup
!
!
pure function nxa_sigma(ion,k,nu) result(s)
!-- photoionization cross-section of level k at frequency nu [cm^2]
  type(nxion),intent(in) :: ion
  integer,intent(in) :: k
  real*8,intent(in) :: nu
  real*8 :: s,u
  integer :: a,b,m,i0,i1
  s = 0d0
  u = nu/ion%nu0(k)
  if(u < 1d0) return
  i0 = ion%pst(k)
  i1 = i0 + ion%pn(k) - 1
  if(u >= ion%pu(i1)) then
     s = ion%ps(i1)*(ion%pu(i1)/u)**3
     return
  endif
  if(u <= ion%pu(i0)) then
     s = ion%ps(i0)
     return
  endif
!-- bracket
  a = i0
  b = i1
  do while(b-a>1)
     m = (a+b)/2
     if(ion%pu(m) <= u) then
        a = m
     else
        b = m
     endif
  enddo
  s = ion%ps(a) + (ion%ps(b) - ion%ps(a))*(u - ion%pu(a))/(ion%pu(b) - ion%pu(a))
end function nxa_sigma
!
!
subroutine milne(ion,k,temp,alpha,rc)
!-----------------------------------------------------------------------
! Spontaneous radiative recombination to level k (Milne relation with the
! tabulated cross-section) and its cooling coefficient (kinetic energy
! h(nu-nu0) per recombination):
!   alpha = 4 pi phi T^-1.5 sahac (2/c^2) Int sigma nu^2 e^{-h(nu-nu0)/kT} dnu
! phi = g/(2 g_+).  Integrated in y = h(nu-nu0)/kT to y = 50 with Simpson's
! rule on each interval of the cross-section table (subdivided to dy <= 0.05).
!-----------------------------------------------------------------------
  type(nxion),intent(in) :: ion
  integer,intent(in) :: k
  real*8,intent(in) :: temp
  real*8,intent(out) :: alpha,rc
  real*8 :: nu0,ykt,ya,yb,s0,s1,y,h,w,nu,sg,f0
  integer :: i0,i1,i,j,ns
  real*8,parameter :: ymax = 50d0, dymax = 0.05d0
!
  nu0 = ion%nu0(k)
  ykt = hh*nu0/(kb*temp)              !y per unit (u-1)
  i0 = ion%pst(k)
  i1 = i0 + ion%pn(k) - 1
  s0 = 0d0
  s1 = 0d0
  ya = 0d0
  i = i0
  do
!-- next node in y: table nodes, then geometric steps beyond the table
     if(i < i1) then
        yb = (ion%pu(i+1) - 1d0)*ykt
     else
        yb = max(ya*1.2d0, ya + 1d0)
     endif
     yb = min(yb,ymax)
     if(yb > ya) then
        ns = 2*max(1,ceiling((yb-ya)/dymax/2d0))
        h = (yb - ya)/ns
        do j=0,ns
           y = ya + j*h
           w = merge(1d0, merge(4d0, 2d0, mod(j,2)==1), j==0 .or. j==ns)
           nu = nu0*(1d0 + y/ykt)
           sg = nxa_sigma(ion,k,nu)
           f0 = w*sg*nu*nu*exp(-y)
           s0 = s0 + f0*h/3d0
           s1 = s1 + f0*y*h/3d0
        enddo
     endif
     ya = yb
     i = i + 1
     if(ya >= ymax) exit
  enddo
!-- dnu = (kT/h) dy
  alpha = 4d0*pi*ion%g(k)/(2d0*ion%gion)*nxa_sahac*temp**(-1.5d0)*2d0/cc**2*(kb*temp/hh)*s0
  rc = alpha/max(s0,1d-300)*s1*kb*temp
end subroutine milne
!
!
pure function nxa_alpha(lv,temp) result(a)
!-- radiative recombination coefficient of global level lv [cm^3/s]
  integer,intent(in) :: lv
  real*8,intent(in) :: temp
  real*8 :: a,lt,w
  integer :: it
  lt = min(max(log(temp),nxa_ltt(1)),nxa_ltt(nxa_ntt))
  w = (lt - nxa_ltt(1))/(nxa_ltt(2) - nxa_ltt(1))
  it = min(int(w) + 1,nxa_ntt-1)
  w = w - (it-1)
  a = exp((1d0-w)*nxa_lalp(it,lv) + w*nxa_lalp(it+1,lv))
end function nxa_alpha
!
!
pure function nxa_rcool(lv,temp) result(a)
!-- recombination cooling coefficient of global level lv [erg cm^3/s]
  integer,intent(in) :: lv
  real*8,intent(in) :: temp
  real*8 :: a,lt,w
  integer :: it
  lt = min(max(log(temp),nxa_ltt(1)),nxa_ltt(nxa_ntt))
  w = (lt - nxa_ltt(1))/(nxa_ltt(2) - nxa_ltt(1))
  it = min(int(w) + 1,nxa_ntt-1)
  w = w - (it-1)
  a = exp((1d0-w)*nxa_lrc(it,lv) + w*nxa_lrc(it+1,lv))
end function nxa_rcool
!
!
subroutine nxa_groups(ng,wlg)
!-----------------------------------------------------------------------
! Group integrals of the cross-sections, exact for the tabulated
! (piecewise linear in u) cross-section and its nu^-3 extension:
!   sgm(g) = (1/dnu_g) Int_g sigma dnu,   sgh(g) = Int_g sigma/(h nu) dnu
!-----------------------------------------------------------------------
  integer,intent(in) :: ng
  real*8,intent(in) :: wlg(ng+1)       !group boundaries [cm], increasing
  integer :: lv,i,k,ig,j,i0,i1
  real*8 :: nulo,nuhi,ua,ub,u1,u2,s1,s2,sl,nu0,dnu,ul,sl3
!
  nxa_ng = ng
  if(allocated(nxa_sgm)) deallocate(nxa_sgm,nxa_sgh,nxa_g1,nxa_lgn)
  allocate(nxa_sgm(ng,nxa_nlv),nxa_sgh(ng,nxa_nlv),nxa_g1(nxa_nlv),nxa_lgn(nxa_nlin))
!-- group of each line: nearest group centre in frequency
  do k=1,nxa_nlin
     nxa_lgn(k) = 1
     dnu = 1d99
     do ig=1,ng
        nulo = 0.5d0*(cc/wlg(ig) + cc/wlg(ig+1))
        if(abs(nulo - nxa_lnu(k)) < dnu) then
           dnu = abs(nulo - nxa_lnu(k))
           nxa_lgn(k) = ig
        endif
     enddo
  enddo
  nxa_sgm = 0d0
  nxa_sgh = 0d0
  do lv=1,nxa_nlv
     i = nxa_lvion(lv)
     k = lv - nxa_ion(i)%lvoff
     nu0 = nxa_ion(i)%nu0(k)
     i0 = nxa_ion(i)%pst(k)
     i1 = i0 + nxa_ion(i)%pn(k) - 1
     nxa_g1(lv) = 0
     do ig=ng,1,-1
        nulo = cc/wlg(ig+1)
        nuhi = cc/wlg(ig)
        if(nuhi <= nu0) cycle
        nxa_g1(lv) = max(nxa_g1(lv),ig)
        ua = max(nulo,nu0)/nu0
        ub = nuhi/nu0
        dnu = nuhi - nulo
!-- table segments overlapping [ua,ub]
        do j=i0,i1-1
           u1 = max(nxa_ion(i)%pu(j),ua)
           u2 = min(nxa_ion(i)%pu(j+1),ub)
           if(u2 <= u1) cycle
           sl = (nxa_ion(i)%ps(j+1) - nxa_ion(i)%ps(j))/(nxa_ion(i)%pu(j+1) - nxa_ion(i)%pu(j))
           s1 = nxa_ion(i)%ps(j) + sl*(u1 - nxa_ion(i)%pu(j))
           s2 = nxa_ion(i)%ps(j) + sl*(u2 - nxa_ion(i)%pu(j))
!-- Int sigma du (trapezoid exact for linear sigma); Int sigma/u du
           nxa_sgm(ig,lv) = nxa_sgm(ig,lv) + 0.5d0*(s1 + s2)*(u2 - u1)*nu0
           nxa_sgh(ig,lv) = nxa_sgh(ig,lv) + ((s1 - sl*u1)*log(u2/u1) + sl*(u2 - u1))/hh
        enddo
!-- below the first table point: its value
        ul = nxa_ion(i)%pu(i0)
        if(ua < ul) then
           u1 = ua
           u2 = min(ul,ub)
           s1 = nxa_ion(i)%ps(i0)
           nxa_sgm(ig,lv) = nxa_sgm(ig,lv) + s1*(u2 - u1)*nu0
           nxa_sgh(ig,lv) = nxa_sgh(ig,lv) + s1*log(u2/u1)/hh
        endif
!-- above the last point: sigma_l (u_l/u)^3
        ul = nxa_ion(i)%pu(i1)
        if(ub > ul) then
           u1 = max(ua,ul)
           u2 = ub
           sl3 = nxa_ion(i)%ps(i1)*ul**3
           nxa_sgm(ig,lv) = nxa_sgm(ig,lv) + sl3*0.5d0*(1d0/u1**2 - 1d0/u2**2)*nu0
           nxa_sgh(ig,lv) = nxa_sgh(ig,lv) + sl3*(1d0/u1**3 - 1d0/u2**3)/3d0/hh
        endif
        nxa_sgm(ig,lv) = nxa_sgm(ig,lv)/dnu
     enddo
  enddo
end subroutine nxa_groups
!
!
pure function ex_e1x(x) result(r)
!-- e^x E1(x), x > 0 (Abramowitz & Stegun 5.1.53, 5.1.56)
  real*8,intent(in) :: x
  real*8 :: r
  if(x <= 1d0) then
     r = exp(x)*(-log(x) - 0.57721566d0 + x*(0.99999193d0 + x*(-0.24991055d0 + &
          x*(0.05519968d0 + x*(-0.00976004d0 + x*0.00107857d0)))))
  else
     r = (x**4 + 8.5733287401d0*x**3 + 18.0590169730d0*x**2 + &
          8.6347608925d0*x + 0.2677737343d0)/(x**4 + 9.5733223454d0*x**3 + &
          25.6329561486d0*x**2 + 21.0996530827d0*x + 3.9584969228d0)/x
  endif
end function ex_e1x
!
!
pure subroutine nxa_coll(ion,temp,qup,qdn,cion,c3b)
!-----------------------------------------------------------------------
! Electron-impact rates of one ion at temperature temp:
!  qup(l,u), qdn(u,l)  excitation, de-excitation [cm^3/s] (all pairs l < u)
!  cion(k)  ionization from level k [cm^3/s]
!  c3b(k)   three-body recombination to k per ion [cm^6/s] (detailed balance)
!-----------------------------------------------------------------------
  type(nxion),intent(in) :: ion
  real*8,intent(in) :: temp
  real*8,intent(out) :: qup(ion%nlev,ion%nlev),qdn(ion%nlev,ion%nlev)
  real*8,intent(out) :: cion(ion%nlev),c3b(ion%nlev)
  real*8 :: om(ion%nlev,ion%nlev),omi(ion%nlev)
  real*8 :: x,gbar,nu15,f,w,sq,u,sig0,g1
  integer :: l,u_,k,ic,it,n
!
  n = ion%nlev
  om = 0d0
  omi = 0d0
!-- tabulated: log-log interpolation in T, constant beyond the table
  if(ion%ncol>0) then
     if(temp <= ion%colt(1)) then
        it = 1; w = 0d0
     elseif(temp >= ion%colt(ion%ncolt)) then
        it = ion%ncolt - 1; w = 1d0
     else
        it = 1
        do while(temp > ion%colt(it+1))
           it = it + 1
        enddo
        w = log(temp/ion%colt(it))/log(ion%colt(it+1)/ion%colt(it))
     endif
     do ic=1,ion%ncol
        if(ion%com(it,ic) > 0d0 .and. ion%com(it+1,ic) > 0d0) then
           x = exp((1d0-w)*log(ion%com(it,ic)) + w*log(ion%com(it+1,ic)))
        else
           x = (1d0-w)*ion%com(it,ic) + w*ion%com(it+1,ic)
        endif
        if(ion%cu(ic)==0) then
           omi(ion%cl(ic)) = x
        else
           om(ion%cl(ic),ion%cu(ic)) = x
        endif
     enddo
  endif
!-- missing pairs (CMFGEN OMEGA_GEN_V3)
  do l=1,n-1
     do u_=l+1,n
        if(om(l,u_) > 0d0) cycle
        k = ion%lin(l,u_)
        f = 0d0
        if(k>0) f = ion%f(k)
        if(f <= 1d-5) then
           om(l,u_) = ion%omset
        else
           nu15 = (ion%e(u_) - ion%e(l))/hh*1d-15
           x = (ion%e(u_) - ion%e(l))/(kb*temp)
           if(ion%zion <= 1.0001d0) then
              if(x <= 14d0) then
                 gbar = 0.276d0*ex_e1x(x)
              else
                 gbar = 0.066d0*(1d0 + 1.5d0/x)/sqrt(x)
              endif
           else
              g1 = 0.2d0
              if(ion%samen(k)) g1 = 0.7d0
              gbar = max(g1,0.276d0*ex_e1x(x))
           endif
           om(l,u_) = 47.972d0*ion%omscale*gbar*f*ion%g(l)/nu15
        endif
     enddo
  enddo
!-- rates
  sq = nxa_colfac*8.63d-6/sqrt(temp)
  qup = 0d0
  qdn = 0d0
  do l=1,n-1
     do u_=l+1,n
        x = (ion%e(u_) - ion%e(l))/(kb*temp)
        qdn(u_,l) = sq*om(l,u_)/ion%g(u_)
        qup(l,u_) = sq*om(l,u_)/ion%g(l)*exp(-min(x,700d0))
     enddo
  enddo
  do k=1,n
     u = hh*ion%nu0(k)/(kb*temp)
     if(omi(k) > 0d0) then
        cion(k) = sq*omi(k)/ion%g(k)*exp(-min(u,700d0))
        c3b(k) = sq*omi(k)/ion%g(k)*ion%g(k)/(2d0*ion%gion)*nxa_sahac*temp**(-1.5d0)
     else
!-- Seaton with gbar = 0.1 min(zion,3)
        sig0 = nxa_sigma(ion,k,ion%nu0(k))
        gbar = 0.1d0*min(ion%zion,3d0)
        cion(k) = nxa_colfac*1.55d13*gbar*sig0/sqrt(temp)/u*exp(-min(u,700d0))
        c3b(k) = nxa_colfac*1.55d13*gbar*sig0/sqrt(temp)/u*ion%g(k)/(2d0*ion%gion)*nxa_sahac*temp**(-1.5d0)
     endif
  enddo
end subroutine nxa_coll
!
subroutine nxa_read_cx(fname,lprint)
!-----------------------------------------------------------------------
! Charge-exchange reactions (Tools/cmfgen2superlite_chgexch.py) among the
! NLTE elements; call after nxa_setup.  A participant's term is matched by
! name (no '_', no J; alternatives separated by '/') with the levels of its
! NLTE ion, else it is the ground term; 'ION' or a stage above the NLTE
! ions is the top state.  Reactions with a participant outside the model
! are dropped.
!-----------------------------------------------------------------------
  character(*),intent(in) :: fname
  logical,intent(in) :: lprint
  integer,parameter :: mx = 200
  type(nxcx) :: t(mx)
  integer :: iu,ios,n,k,iel,i,j,z(4),st(4),nt,lv
  character(512) :: s
  character(64) :: w(20),nm(4)
  logical :: ok,lgr
  real*8 :: e(4)
  if(allocated(nxa_cx)) deallocate(nxa_cx)
  nxa_ncx = 0
  open(newunit=iu,file=fname,status='old',action='read',iostat=ios)
  if(ios/=0) then
     if(lprint) write(6,'(1x,2a)') 'NLTE: no charge-exchange data: ',trim(fname)
     allocate(nxa_cx(0))
     return
  endif
  n = 0
  do
     call getrec(iu,s,ios)
     if(ios/=0) exit
     call words(s,w,nt)
     if(nt<19) stop 'nxa_read_cx: short record'
     do k=1,4
        read(w(3*k-2),*) z(k)
        read(w(3*k-1),*) st(k)
        nm(k) = w(3*k)
     enddo
     if(n>=mx) stop 'nxa_read_cx: too many reactions'
     n = n + 1
     read(w(13),*) t(n)%typ
     read(w(14:17),*) t(n)%c
     read(w(18),*) t(n)%tlo
     read(w(19),*) t(n)%thi
     ok = z(1)==z(3) .and. z(2)==z(4) .and. z(1)/=z(2)
     do k=1,4
        if(.not.ok) exit
        call participant(z(k),st(k),nm(k),t(n),k,ok,lgr)
        if(ok.and.lgr.and.lprint) write(6,'(1x,a,i3,i2,1x,2a)') 'NLTE charge exchange: no level ', &
             z(k),st(k),trim(nm(k)),' (ground term taken)'
     enddo
     if(.not.ok) then
        n = n - 1
        cycle
     endif
!-- energy released by the forward reaction, lowest level of each term
     do k=1,4
        iel = t(n)%el(k)
        e(k) = minval(nxa_el(iel)%eps(t(n)%s(1:t(n)%nl(k),k)))
     enddo
     t(n)%q = e(1) - e(3) + e(2) - e(4)
  enddo
  close(iu)
  nxa_ncx = n
  allocate(nxa_cx(n))
  nxa_cx = t(1:n)
  if(lprint) then
     do j=1,n
        write(6,'(1x,a)',advance='no') 'NLTE charge exchange:'
        do k=1,4
           i = nxa_cx(j)%ion(k)
           if(k==3) write(6,'(a)',advance='no') ' ->'
           if(k==2.or.k==4) write(6,'(a)',advance='no') ' +'
           if(i>0) then
              lv = nxa_cx(j)%s(1,k) - nxa_ion(i)%soff
              write(6,'(1x,a,1x,a,a,i0,a)',advance='no') trim(nxa_ion(i)%name), &
                   trim(nxa_ion(i)%lname(lv)),'(',nxa_cx(j)%nl(k),')'
           else
              write(6,'(1x,a,i0,a)',advance='no') 'top(Z=',nxa_el(nxa_cx(j)%el(k))%z,')'
           endif
        enddo
        write(6,'(a,f7.3,a)') '  Q =',nxa_cx(j)%q/ev,' eV'
     enddo
  endif
contains
  subroutine words(str,wd,nw)
    character(*),intent(in) :: str
    character(64),intent(out) :: wd(:)
    integer,intent(out) :: nw
    integer :: p,q
    nw = 0
    p = 1
    do while(p<=len_trim(str))
       if(str(p:p)==' ') then
          p = p + 1
          cycle
       endif
       q = index(str(p:),' ')
       if(q==0) q = len_trim(str(p:)) + 1
       if(nw<size(wd)) then
          nw = nw + 1
          wd(nw) = str(p:p+q-2)
       endif
       p = p + q
    enddo
  end subroutine words
  subroutine participant(zz,stg,name,r,kk,found,lground)
    integer,intent(in) :: zz,stg,kk
    character(*),intent(in) :: name
    type(nxcx),intent(inout) :: r
    logical,intent(out) :: found,lground
    integer :: je,jj,ji,l,m,a0,a1
    character(64) :: alt,lnm,gnm
    found = .false.
    lground = .false.
    je = 0
    do jj=1,nxa_nel
       if(nxa_el(jj)%z==zz) je = jj
    enddo
    if(je==0) return
    r%el(kk) = je
    ji = 0
    do jj=1,nxa_el(je)%nion
       if(nxa_ion(nxa_el(je)%ion(jj))%stage==stg) ji = nxa_el(je)%ion(jj)
    enddo
    if(ji==0) then
!-- the top state
       if(stg/=nxa_ion(nxa_el(je)%ion(nxa_el(je)%nion))%stage+1) return
       r%ion(kk) = 0
       r%nl(kk) = 1
       r%s(1,kk) = nxa_el(je)%ns
       found = .true.
       return
    endif
    r%ion(kk) = ji
    m = 0
    if(name/='ION') then
       a0 = 1
       do
          a1 = index(name(a0:),'/')
          if(a1==0) then
             alt = name(a0:)
          else
             alt = name(a0:a0+a1-2)
          endif
          do l=1,nxa_ion(ji)%nlev
             lnm = normname(nxa_ion(ji)%lname(l))
             if(lnm==alt .and. m<nxa_cxml .and. .not.any(r%s(1:m,kk)==nxa_ion(ji)%soff+l)) then
                m = m + 1
                r%s(m,kk) = nxa_ion(ji)%soff + l
             endif
          enddo
          if(a1==0) exit
          a0 = a0 + a1
       enddo
    endif
    if(m==0) then
!-- ground term
       lground = name/='ION'
       gnm = normname(nxa_ion(ji)%lname(1))
       do l=1,nxa_ion(ji)%nlev
          if(normname(nxa_ion(ji)%lname(l))==gnm .and. m<nxa_cxml) then
             m = m + 1
             r%s(m,kk) = nxa_ion(ji)%soff + l
          endif
       enddo
    endif
    r%nl(kk) = m
    found = m>0
  end subroutine participant
  function normname(a) result(b)
!-- level name without '_' and the J value in brackets
    character(*),intent(in) :: a
    character(64) :: b
    integer :: p,q
    b = ' '
    q = 0
    do p=1,len_trim(a)
       if(a(p:p)=='[') exit
       if(a(p:p)=='_') cycle
       q = q + 1
       b(q:q) = a(p:p)
    enddo
  end function normname
end subroutine nxa_read_cx
!
!
pure function nxa_cx_coef(k,temp) result(alpha)
!-- forward rate coefficient of charge-exchange reaction k [cm^3/s]
  integer,intent(in) :: k
  real*8,intent(in) :: temp
  real*8 :: alpha,t4
  t4 = min(max(temp*1d-4,nxa_cx(k)%tlo),nxa_cx(k)%thi)
  t4 = max(t4,1d-4)
  select case(nxa_cx(k)%typ)
  case(1)
     alpha = nxa_cx(k)%c(1)*t4**nxa_cx(k)%c(2)
  case(2)
     alpha = nxa_cx(k)%c(1)*t4**nxa_cx(k)%c(2)*(1d0 + nxa_cx(k)%c(3)*exp(nxa_cx(k)%c(4)*t4))
  case(3)
     alpha = nxa_cx(k)%c(1)*t4**nxa_cx(k)%c(2)*exp(nxa_cx(k)%c(3)*t4)
  case default
     alpha = 0d0
  end select
  alpha = max(alpha,0d0)
end function nxa_cx_coef
!
!
end module nlteatommod
! vim: fdm=marker
