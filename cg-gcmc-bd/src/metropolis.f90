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

subroutine metropolis(nmcm)
! Perform nmcm of metropolis monte carlo with constant number of particle

use grandmod
use constamod
implicit none
integer nmcm
integer iat, imove
!real*16 bltz,eold, enew
real bltz,eold, enew
real xnew, ynew, znew

do imove = 1, nmcm

!pick one atom and move it 
   call move(iat,xnew,ynew,znew)

!calculate new energy
   call interact(eold,x(iat),y(iat),z(iat),abs(typei(iat)),iat,.true.)
   call interact(enew,xnew,ynew,znew,abs(typei(iat)),iat,.false.)

   if (enew.le.eold) then    
      x(iat) = xnew       !accept the move
      y(iat) = ynew
      z(iat) = znew
      ener = ener + (enew-eold)
   else                            
      bltz = exp(-(enew-eold)*ikbt)
      if (rndm().lt.bltz) then
         x(iat) = xnew    !accept the move
         y(iat) = ynew
         z(iat) = znew
         ener = ener + (enew-eold)
      endif
   endif 

enddo

return
end
