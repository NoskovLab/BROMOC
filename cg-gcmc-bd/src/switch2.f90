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

subroutine switch2(sw,dsw,r2,radius,p2,ip2,r)
implicit none
real sw, dsw, r, radius, p2, r2
real delr,delrp2,ip2

if (r2.le.radius**2) then  ! inside pore
  sw  = 0.0
  dsw = 0.0
elseif (r2.lt.(radius+p2)**2) then ! inside pore wall
  r=sqrt(r2)
  delr = (radius+p2-r)
  delrp2=delr*ip2
  sw = 1.0 + (2.0*delrp2-3.0)*delrp2**2
  dsw = 6.0*(delrp2-1.0)*delrp2*ip2
else ! inside membrane
 sw = 1.0
 dsw = 0.0
endif

return
end subroutine
