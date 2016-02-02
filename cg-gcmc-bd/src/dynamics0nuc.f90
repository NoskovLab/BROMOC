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

SUBROUTINE DYNAMICS0NUC(ninit,nfinal)
use grandmod
use constamod
use stdiomod
use nucleotmod
use errormod

implicit none
integer ninit, nfinal, i
real,external :: rgauss
real delx, dely, delz

do i = ninit, nfinal
  if (stfree(i)) then
    delx = fx(i)*fact0n1
    dely = fy(i)*fact0n1
    delz = fz(i)*fact0n1
    if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
    if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
    if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)
    x(i) = x(i) + delx + fact0n2*rgauss()
    y(i) = y(i) + dely + fact0n2*rgauss()
    z(i) = z(i) + delz + fact0n2*rgauss()
    if (.not.Qdnafree) call fixcoor(x(i),y(i),z(i))
  endif 
enddo
return
end
