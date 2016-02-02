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

real function angdih(vecR)
implicit none
!Input
real  vecR(3,3)
!Local variables
real  vecN(3,2), R2, sinphi, cosphi
integer i

vecN(1,1) = vecR(2,1)*vecR(3,2) - vecR(3,1)*vecR(2,2)
vecN(2,1) = vecR(3,1)*vecR(1,2) - vecR(1,1)*vecR(3,2)
vecN(3,1) = vecR(1,1)*vecR(2,2) - vecR(2,1)*vecR(1,2)
vecN(1,2) = vecR(2,2)*vecR(3,3) - vecR(3,2)*vecR(2,3)
vecN(2,2) = vecR(3,2)*vecR(1,3) - vecR(1,2)*vecR(3,3)
vecN(3,2) = vecR(1,2)*vecR(2,3) - vecR(2,2)*vecR(1,3)
R2 = sqrt(vecR(1,2)**2+vecR(2,2)**2+vecR(3,2)**2)
sinphi = 0.0
cosphi = 0.0
do i = 1, 3
  sinphi = sinphi + vecR(i,1)*vecN(i,2)
  cosphi = cosphi + vecN(i,1)*vecN(i,2)
enddo
sinphi = R2*sinphi
angdih = atan2(sinphi,cosphi)

return
end function
