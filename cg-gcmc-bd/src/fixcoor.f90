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

subroutine FIXCOOR(xi,yi,zi)
use grandmod
use constamod
implicit none
real  xi,yi,zi
real  r,r2,lzrad

! spherical system
if (Qsphere) then
  r2 = (xi-cx)**2+(yi-cy)**2+(zi-cz)**2
  if (r2.gt.Rsphe2) then
    lzrad=2.0*rsphe/sqrt(r2)-1.0
    xi = lzrad*(xi-cx)+cx
    yi = lzrad*(yi-cy)+cy
    zi = lzrad*(zi-cz)+cz
  endif
! cyllindrical system
elseif (Qecyl) then
  r2=((xi-cx)*iecx)**2+((yi-cy)*iecy)**2
  if (r2.gt.1.0) then
     r=0.999999999999999/sqrt(r2)
     xi = r*(xi-cx)+cx
     yi = r*(yi-cy)+cy
  endif
  if (zi.lt.lz2m) then
    zi = lz2m
  else if (zi.gt.lz2p) then
    zi = lz2p
  endif
! rectangular system
else
  if (xi.lt.lx2m) then
    xi = lx2m
  else if (xi.gt.lx2p) then !ge
    xi = lx2p
  endif
  if (yi.lt.ly2m) then
    yi = ly2m
  else if (yi.gt.ly2p) then !ge
    yi = ly2p
  endif
  if (zi.lt.lz2m) then
    zi = lz2m
  else if (zi.gt.lz2p) then
    zi = lz2p
  endif
endif

return
end
