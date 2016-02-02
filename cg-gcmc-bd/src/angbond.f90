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

real function angbond(vecR)
implicit none
!Input
real  vecR(3,3)
!Local variables
real  R(2), prod, cosphi
integer i

do i = 1, 2
  R(i) = sqrt(vecR(1,i)**2+vecR(2,i)**2+vecR(3,i)**2)
enddo
prod = 0.0
do i = 1, 3
  prod = prod + vecR(i,1)*vecR(i,2)
enddo
cosphi = prod/(R(1)*R(2))
angbond = acos(cosphi)

return
end function

