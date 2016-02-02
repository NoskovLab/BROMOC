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

SUBROUTINE DYNAMICS0(ninit,nfinal)
use grandmod
use constamod
use stdiomod
use errormod

implicit none
integer ninit, nfinal, itype, i
integer nfloc(1:ntype-nold,cntpts),nbloc(1:ntype-nold,cntpts)
real,external :: rgauss
real delx, dely, delz
real zold

!$omp parallel private(i,itype,zold,delx,dely,delz,nfloc,nbloc)
nfloc=0
nbloc=0
!$omp do
do i = ninit, nfinal
  itype = abs(typei(i))
  delx = fx(i)*fact1a(itype)
  dely = fy(i)*fact1a(itype)
  delz = fz(i)*fact1a(itype)
  if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
  if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
  if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)
  x(i) = x(i) + delx + fact2a(itype)*rgauss()
  y(i) = y(i) + dely + fact2a(itype)*rgauss()
  zold = z(i)
  z(i) = z(i) + delz + fact2a(itype)*rgauss()
! Keep track of the net flux of particle
  if (Qcountion) call countions(zold,z(i),itype-nold,nfloc,nbloc)
  call fixcoor(x(i),y(i),z(i))
enddo
!$omp end do
!$omp critical
nforward = nforward + nfloc
nbackward = nbackward + nbloc
!$omp end critical
!$omp end parallel
return
end subroutine

