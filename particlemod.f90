!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
module particlemod

  implicit none

  type packet
     real*8 :: x, y, z
     real*8 :: mu, om
     real*8 :: e, e0, wl
!-- last emission (diagnostics, output.emitchan): 1 inner boundary,
!-- 2 thermal emissivity, 3 bound-free, 4 line, 5 two-photon, 6 shock;
!-- + 100 Z of the emitting element for 3-5
     integer :: em = 0
!-- cell of the last thermal or bound-free emission (part of the cell's
!-- source function), 0 after any other emission (nxt_ldb)
     integer :: icem = 0
!-- NLTE line of the last emission while the photon stays to the red of
!-- it in the comoving frame (thermal electron scattering may move it to
!-- the blue), 0 otherwise: in an expanding flow it is then not incident on
!-- that line (nxt_tally_finish)
     integer :: ilem = 0
!-- (nxt_esrates) 1: a photon just scattered out of a line core by an
!-- electron (em_lines at the emission of the packet): where it lands in
!-- the line profiles is tested on its first transport step (nxt_esfar)
     integer :: est = 0
  end type packet
!
!-- secondary particle properties
  type packet2
     real*8 :: dist          !particle travel distance
     character(4) :: stat    !particle status: live, cens, flux, dead
     integer :: ix, iy, iz   !positional cell indices
     integer :: ic, ig       !index into compressed domain arrays, group index
     integer :: itype        !IMC or DDMC type
     integer :: ipart, istep !particle number and transport step number
     integer :: idist        !transport distance identifier
     integer :: iline = 0    !sorted H line the packet sits on (comoving wavelength = line), 0 if none
     integer :: ilside = 0   !+1 red side, -1 blue side of that line
!-- (diagnostics: a line emitted by the macro-atom on this step, its cell and
!-- energy, and a line absorbed at its resonance; nxt_lemab)
     integer :: kem = 0, cem = 0, kab = 0
     real*8 :: eem = 0d0, eab = 0d0
!-- a line that absorbed the packet in its core after an electron
!-- scattering on this step, and the comoving energy (nxt_esfar, nxt_efsum)
     integer :: kfa = 0
     real*8 :: efa = 0d0
!-- (nxt_esrates) lines ahead of the packet whose profiles were tested for
!-- its whole path at its last electron scattering: their resonances are
!-- passed without absorption until the next scattering or emission
!-- (nxt_esfar); estest: a photon re-emitted as em_lines, to be tested on
!-- its next transport step
     integer :: kskip(4) = 0
     logical :: estest = .false.
  end type packet2
!
  type(packet),allocatable,target :: prt_particles(:)  !(prt_npartmax)
  logical,allocatable :: prt_isvacant(:)  !(prt_npartmax)
!
  integer :: prt_npartmax

  save

  contains

  subroutine particle_alloc(ltalk)
!--------------------------------------------------
    implicit none
    logical,intent(in) :: ltalk

    integer :: n

!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_particles(prt_npartmax),prt_isvacant(prt_npartmax))
    prt_isvacant = .true.
!
!-- print size only on master rank
    if(ltalk) then
      n = int(sizeof(prt_particles)/1024) !kB
      write(6,*) 'ALLOC particles:',n,"kB",n/1024,"MB",n/1024**2,"GB"
    endif
!
!-- output
    if(ltalk) then
       write(6,*)
       write(6,*) 'particle array:'
       write(6,*) '===================='
       write(6,*) 'npart :',prt_npartmax,nint(prt_npartmax/1000d0),'k'
       write(6,*)
    endif

  end subroutine particle_alloc


end module particlemod
! vim: fdm=marker
