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

real Function M3(x)
!------------------------------------------------------------------------
implicit none
real   M2m,M2,X
if((x.ge.0.0).and.(x.le.2.0)) then
   M2m = 1.0 - abs(x - 1.0)
else
   M2m = 0.0
endif
if((x.ge.1.0).and.(x.le.3.0)) then
   M2  = 1.0 - abs(x - 2.0)
else
   M2  = 0.0
endif

if((x.ge.0.0).and.(x.le.3.0d0)) then
   M3 = x*0.5*M2m + (3.0-x)*0.5*M2
else
   M3 = 0.0
endif

return
end function
