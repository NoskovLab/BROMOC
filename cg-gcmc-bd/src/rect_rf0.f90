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

subroutine rect_rf0(xj,yj,zj,j,jtype)
!-----------------------------------------------------------------------
!     This subroutine computes only the reaction field energy difference
!     due to jth ion in subroutine INTERACT in simul.f
!
use ioxmod
use constamod
use grandmod
use nucleotmod
use gsbpmod
!Inputs
implicit none
integer j,jtype
real  xj,yj,zj
!local variables
real  norm,coefi,coefj,coef2(ntpol)
real  xg,yg,zg,xxs,yys,zs
real  lpolx(ntot,xnpol),lpoly(ntot,ynpol),lpolz(ntot,znpol)
real  charge
integer i,ii,jj,ij,n
integer xpol,ypol,zpol,itype

!calculate Q_{m} coefficients
do ii = 1, ntpol
   coef(ii) = 0.0
enddo
do i = nsites+nfix+1, ntot
  if (i.eq.j) then
    itype = jtype
    xg = xj
    yg = yj
    zg = zj
  else
    itype = abs(typei(i))
    xg = x(i)             
    yg = y(i)             
    zg = z(i)             
  endif
  charge = cg(itype)
  if (zg.ge.zmin.and.zg.le.zmax .and.xg.ge.xmin.and.xg.le.xmax .and.yg.ge.ymin.and.yg.le.ymax) then
    xxs = xscal*xg
    yys = yscal*yg
    zs = zscal*zg
    lpolx(i,1) = 1.0
    lpoly(i,1) = 1.0
    lpolz(i,1) = 1.0
    lpolx(i,2) = xxs
    lpoly(i,2) = yys
    lpolz(i,2) = zs
    lpolx(i,3) = 0.5*(3.0*xxs*xxs-1.0)
    lpoly(i,3) = 0.5*(3.0*yys*yys-1.0)
    lpolz(i,3) = 0.5*(3.0*zs*zs-1.0)
    do n = 3, xnpol-1
      lpolx(i,n+1) = ((2.0*n-1)*xxs*lpolx(i,n)-(n-1.0)*lpolx(i,n-1))/n
    enddo
    do n = 3, ynpol-1
      lpoly(i,n+1) = ((2.0*n-1)*yys*lpoly(i,n)-(n-1.0)*lpoly(i,n-1))/n
    enddo
    do n = 3, znpol-1
      lpolz(i,n+1) = ((2.0*n-1)*zs*lpolz(i,n)-(n-1.0)*lpolz(i,n-1))/n
    enddo
    do ii = 1, ntpol
      xpol = lstpx(ii)
      ypol = lstpy(ii)
      zpol = lstpz(ii)
      norm = bnorm(ii)
      coef(ii) = coef(ii)+charge*norm*lpolx(i,xpol+1)*lpoly(i,ypol+1)*lpolz(i,zpol+1)
    enddo
  endif
enddo

!calculate Q_{m} coefficients for ntot- (jth ion)
charge = cg(jtype)
do ii = 1, ntpol
  xpol = lstpx(ii)
  ypol = lstpy(ii)
  zpol = lstpz(ii)
  norm = bnorm(ii)
  coef2(ii) = coef(ii)-charge*norm*lpolx(j,xpol+1)*lpoly(j,ypol+1)*lpolz(j,zpol+1)
enddo

!reaction field energy calculation g(ntot)-g(ntot-1)
egsbpb=0.0
do ii = 1, ntpol
  ij = (ii-1)*ntpol+ii
  coefi = coef(ii)
  coefj = coef2(ii)
  egsbpb = egsbpb+0.5*(coefi*coefi-coefj*coefj)*mmij(ij)
  do jj = ii+1, ntpol
    ij = (ii-1)*ntpol+jj
    egsbpb = egsbpb + (coefi*coef(jj)-coefj*coef2(jj))*mmij(ij)
  enddo
enddo
egsbpb = egsbpb*celec*rfscal

return
end subroutine
