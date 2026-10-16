*This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
      subroutine fluxtally
c     ------------------------
      use fluxmod
      use particlemod
      use miscmod
      use countersmod
      use timingmod
      implicit none
************************************************************************
* Go through all particles and tally the ones that have left the domain
* and end up in this flux time step.
************************************************************************
      integer :: ipart
      integer :: ig,imu,iom,ib,iek,iez
      real*8 :: help
      real*8 :: t0
      type(packet),pointer :: ptcl
c
c-- timer
      t0 = t_time()
c
c-- init
      flx_luminos = 0d0
      flx_lumdev = 0d0
      flx_lumnum = 0
      flx_lumout = 0d0
      flx_emch = 0d0
c
!!c$omp do schedule(static,1) !round-robin
      do ipart=1,prt_npartmax

c-- check vacancy
      if(prt_isvacant(ipart)) cycle
c
c-- active particle
      ptcl => prt_particles(ipart)
c
c-- check flux status
       if(ptcl%x/=huge(help)) cycle
c
c-- retrieving lab frame flux group, polar, azimuthal bin
       ig = binsrch(ptcl%wl,flx_wl,flx_ng+1,.true.)
       imu = binsrch(ptcl%mu,flx_mu,flx_nmu+1,.false.)
       iom = binsrch(ptcl%om,flx_om,flx_nom+1,.false.)
c
c-- by emission process
       ib = binsrch(ptcl%wl,flx_emband,flx_nemb+1,.false.)
       iek = min(max(mod(ptcl%em,100),0),6)
       iez = min(max(ptcl%em/100,0),30)
       flx_emch(iek,iez,ib) = flx_emch(iek,iez,ib) + ptcl%e
c-- tallying outbound luminosity (outside the bins: separate total)
       if(ig<1 .or. ig>flx_ng) then
        flx_lumout = flx_lumout+ptcl%e
       else
        flx_luminos(ig,imu,iom) = flx_luminos(ig,imu,iom)+ptcl%e
        flx_lumdev(ig,imu,iom) = flx_lumdev(ig,imu,iom)+ptcl%e**2
        flx_lumnum(ig,imu,iom) = flx_lumnum(ig,imu,iom)+1
       endif

c-- mark particle slot occupied or vacant
       prt_isvacant(ipart) = .true.

      enddo !ipart
!!c$omp end do nowait
c
c-- timing
      t0 = t_time() - t0
      call timereg(t_fluxtally,t0)
c
      end subroutine fluxtally
c vim: fdm=marker
