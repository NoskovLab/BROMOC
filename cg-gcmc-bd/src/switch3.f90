!    BROMOC  -  CG-GCMC-BD
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program.
!    Copyright (C) 2014 Pablo M. De Biase (pablodebiase@gmail.com)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

subroutine switch3(sw,dsw,z,plength2,p3,pcenter)
!INPUTS      
!plength2 -> half the pore length
!pcenter -> the pore center
!p3(itype) -> switchlength
!z -> Z-axis direction (pore axial)
implicit none
real sw,dsw,z,plength2,p3,pcenter
real delz,pcenter1,pcenter2,ip3,delzp3

pcenter1 = -plength2 + pcenter
pcenter2 = plength2 + pcenter

if (z.gt.pcenter1-p3.and.z.lt.pcenter2+p3) then
  if (z.ge.pcenter1 .and. z.le.pcenter2) then
    sw  = 0.0
    dsw = 0.0
  else
    if (z.gt.pcenter2) then
      delz = (z-pcenter2)
      ip3=1.0/p3
      delzp3=delz*ip3
      sw = (2.0*delzp3-3.0)*delzp3*delzp3
      dsw = 6.0*(delzp3-1.0)*delzp3*ip3
    else if (z.lt.pcenter1) then
      delz = (z-pcenter1)
      ip3=1.0/p3
      delzp3=delz*ip3
      sw = (-2.0*delzp3-3.0)*delzp3*delzp3
      dsw = -6.0*(delzp3+1.0)*ip3*delzp3
    endif
    sw = -sw
    dsw = -dsw
  endif
else
  sw = 1.0
  dsw = 0.0
endif

return
end
