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

SUBROUTINE DYNAMICS2(ninit,nfinal)
use grandmod
use constamod
use stdiomod
use nucleotmod
use errormod
use splinemod
use sevalmod

implicit none
integer itype,i,ninit,nfinal
real fact1, fact2
real rgauss
external rgauss
real delx, dely, delz, delDz
real zold 
real sw,dsw, idiffusion, didiffusion
real zz
integer nfloc(1:ntype-nold,cntpts),nbloc(1:ntype-nold,cntpts)

!$omp parallel private(i,itype,zz,sw,dsw,idiffusion,didiffusion,fact1,fact2,delDz,zold,delx,dely,delz,nfloc,nbloc)
nfloc=0
nbloc=0
!$omp do
do i=ninit,nfinal
  itype = abs(typei(i))
! space-dependent diffusion constant
  zz = z(i)
  if(xs(1).lt.xs(nspline)) then 
    if(zz.lt.xs(1)) then
      sw=ys(1)
      dsw=0.0
    elseif(zz.gt.xs(nspline)) then
      sw=ys(nspline)
      dsw=0.0
    else
      sw=seval(nspline,zz,xs,ys,b,c,d)
      dsw=sevald(nspline,zz,xs,b,c,d)
    endif
  elseif(xs(1).gt.xs(nspline)) then
    if(zz.gt.xs(1)) then
      sw=ys(1)
      dsw=0.0
    elseif(zz.lt.xs(nspline)) then
      sw=ys(nspline)
      dsw=0.0
    else
      sw=seval(nspline,zz,xs,ys,b,c,d)
      dsw=sevald(nspline,zz,xs,b,c,d)
    endif
  endif
!   sw  = seval (nspline,zz,xs,ys,b,c,d)
!   dsw = sevald(nspline,zz,xs,b,c,d)
  idiffusion = diffusion(itype)*sw
  didiffusion = diffusion(itype)*dsw
  fact1 = idiffusion*kBTdt
  delx  = fx(i)*fact1
  dely  = fy(i)*fact1
  delz  = fz(i)*fact1
  fact2 = sqrt(2.0*dt*idiffusion)
  delDz = didiffusion*dt
  if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
  if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
  if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)
  x(i)=x(i)+delx+fact2*rgauss()
  y(i)=y(i)+dely+fact2*rgauss()
  zold=z(i)
  z(i)=z(i)+delz+fact2*rgauss()+delDz

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
