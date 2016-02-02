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

subroutine dynamics1(ninit,nfinal)
use grandmod
use constamod
use stdiomod
use nucleotmod
use errormod

implicit none
integer itype,i,ninit,nfinal
real fact1, fact2
real rgauss
external rgauss
real delx, dely, delz, delDz
real zold
real sw, dsw, idiffusion, didiffusion
real zz, pp3, pore1, pore2,ipp3
integer nfloc(1:ntype-nold,cntpts),nbloc(1:ntype-nold,cntpts)

!$omp parallel private(i,itype,zz,pp3,pore1,pore2,sw,dsw,ipp3,idiffusion,didiffusion,fact1,fact2,delDz,zold,delx,dely,delz,nfloc,nbloc)
nfloc=0
nbloc=0
!$omp do
do i = ninit,nfinal
  itype=abs(typei(i))
! space-dependent diffusion constant
  zz = z(i)
  pp3 = p3(itype)
  pore1 = -plength2 + pcenter
  pore2 =  plength2 + pcenter
  if (zz.gt.(pore1-pp3) .and. zz.lt.(pore2+pp3)) then
    if (zz.ge.pore1 .and. zz.le.pore2) then
      sw  = 0.0
      dsw = 0.0
    else
      ipp3=1.0/pp3
      if (zz.gt.pore2) then
        delz = (zz-pore2)
        sw = 2.0*(delz*ipp3)**3-3.0*(delz*ipp3)**2
        dsw = 6.0*((delz*ipp3)**2-delz*ipp3)*ipp3
      else if (zz.lt.pore1) then
        delz = (zz-pore1)
        sw = -2.0*(delz*ipp3)**3-3.0*(delz*ipp3)**2
        dsw = -6.0*((delz*ipp3)**2+delz*ipp3)*ipp3
      endif
      sw = -sw
      dsw = -dsw
    endif
  else
    sw = 1.0
    dsw = 0.0
  endif
  idiffusion = diffusion(itype)*(ampl3(itype)+(1.0-ampl3(itype))*sw)
  didiffusion = diffusion(itype)*(1.0-ampl3(itype))*dsw
  fact1 = idiffusion*kBTdt
  delx  = fx(i)*fact1
  dely  = fy(i)*fact1
  delz  = fz(i)*fact1
  fact2 = sqrt(2.0*dt*idiffusion)
  delDz = didiffusion*dt

  if (abs(delx).gt.bdmax) delx = sign(bdmax,delx)
  if (abs(dely).gt.bdmax) dely = sign(bdmax,dely)
  if (abs(delz).gt.bdmax) delz = sign(bdmax,delz)

  x(i) = x(i) + delx + fact2*rgauss()
  y(i) = y(i) + dely + fact2*rgauss()
  zold = z(i)
  z(i) = z(i) + delz + fact2*rgauss() + delDz

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
