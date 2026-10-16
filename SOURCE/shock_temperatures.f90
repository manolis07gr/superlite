!This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
subroutine shock_temperatures
!-----------------------------------------------------------------------
! Shock cells: the velocity drops outward across the cell by more than
! 50 km/s and 2% (a gently decelerating wind is not a shock).  Contiguous
! shock cells form one shock (hydro codes spread a shock over a few zones);
! its total velocity drop dv, if at least in_shock_dvmin, sets the
! post-shock temperature of its cells, T = 3/16 mu m_p dv^2/k (strong
! shock, mu = 0.6).
!-----------------------------------------------------------------------
  use gridmod
  use sourcemod
  use inputparmod, only:in_shock_dvmin
  use physconstmod
  implicit none
  integer :: i,j
  real*8 :: dv
  real*8,parameter :: mu = 0.6d0
  src_tsh = 0d0
  i = 1
  do while(i<=grd_ncell)
     if(lshock(i)) then
        j = i
        do while(j<grd_ncell)
           if(.not.lshock(j+1)) exit
           j = j + 1
        enddo
        dv = grd_vxarr(i) - grd_vxarr(j+1)
        if(dv >= in_shock_dvmin) src_tsh(i:j) = 3d0/16d0*mu*pc_amu*dv**2/pc_kb
        i = j + 1
     else
        i = i + 1
     endif
  enddo
contains
  pure logical function lshock(k)
    integer,intent(in) :: k
    real*8 :: d
    d = grd_vxarr(k) - grd_vxarr(k+1)
    lshock = d > 5d6 .and. d > 0.02d0*grd_vxarr(k)
  end function lshock
end subroutine shock_temperatures
! vim: fdm=marker
