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

SUBROUTINE MOVE(iat,xnew,ynew,znew)
use grandmod
use nucleotmod

implicit none
integer iat
real xnew, ynew, znew
real radius,natom2
logical*1 endok

!Pick one atom randomly and move it (if it moves outside the 
!limits, pick a different atom)

endok = .true.  
natom2 = float(natom - nsites)
do while (endok)
  iat = nsites + nfix + int(natom2*rndm()) + 1 ! [nsites+nfix+1,ntot]
  xnew = x(iat) + mcmax*(rndm()-0.5)
  ynew = y(iat) + mcmax*(rndm()-0.5)
  znew = z(iat) + mcmax*(rndm()-0.5)
  if (Qsphere) then
    radius = (xnew-cx)**2+(ynew-cy)**2+(znew-cz)**2
    endok = radius.gt.Rsphe2
  elseif (Qecyl) then
    endok = (((xnew-cx)*iecx)**2+((ynew-cy)*iecy)**2).gt.1.0.or.znew.lt.lz2m.or.znew.gt.lz2p
  else
    endok = xnew.lt.lx2m.or.xnew.gt.lx2p.or.ynew.lt.ly2m.or.ynew.gt.ly2p.or.znew.lt.lz2m.or.znew.gt.lz2p
  endif
enddo 
  
return
end
