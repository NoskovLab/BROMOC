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

subroutine insert(ib,xnew,ynew,znew)
!Create cartesian coordinates for a particle from buffer ib
use constamod
use grandmod
implicit none
integer ib
real xnew,ynew,znew
real radius2
logical*1 ok

ok = .false.
do while (.not.ok)
  xnew = lx*(rndm()-0.5)+cx ! [-lx/2,lx/2)
  ynew = ly*(rndm()-0.5)+cy ! [-ly/2,ly/2)
  znew = lzmin(ib)+(lzmax(ib)-lzmin(ib))*rndm() ! [lzmin,lzmax)
  ok = .true.
  if (Qsphere) then
    radius2 = (xnew-cx)**2+(ynew-cy)**2+(znew-cz)**2
    ok = radius2.gt.Rmin(ib)**2.and.radius2.lt.Rmax(ib)**2 ! (Rmin,Rmax)
  elseif (Qecyl) then
    ok = (((xnew-cx)*iecx)**2+((ynew-cy)*iecy)**2).lt.1.0
  endif
enddo

return
end
