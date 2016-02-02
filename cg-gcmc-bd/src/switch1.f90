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

subroutine switch1(sw,dsw,z,p1,ip1,zmemb1,zmemb2)
implicit none
real sw, dsw, z, p1
real zmemb1, zmemb2
real delz,delzp1,ip1

! zmemb1  ! lower limit of the membrane
! zmemb2  ! upper limit of the membrane

if (z.ge.(zmemb2+p1).or.z.le.(zmemb1-p1)) then ! outside membrane
  sw = 0.0
  dsw = 0.0
else
  if (z.gt.zmemb2) then ! z-zmemb2 > 0 ! in upper wall
    delz = z-zmemb2
    delzp1=delz*ip1
    sw = 1.0 + (2.0*delzp1-3.0)*delzp1**2
    dsw = -6.0*(delzp1-1.0)*ip1*delzp1
  elseif (z.lt.zmemb1) then ! z-zmemb1 < 0 ! in lower wall
    delz = zmemb1-z
    delzp1=delz*ip1
    sw = 1.0 + (2.0*delzp1-3.0)*delzp1**2
    dsw = 6.0*(delzp1-1.0)*ip1*delzp1
  else ! inside membrane
    sw = 1.0
    dsw = 0.0
  endif
endif
return
end subroutine
