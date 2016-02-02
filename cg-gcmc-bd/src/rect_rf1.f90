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

subroutine rect_rf1
!-----------------------------------------------------------------------
!     calculate the reaction field energy and forces on each ions
!
use ioxmod
use constamod
use grandmod
use nucleotmod
use gsbpmod     
!local variables
implicit none
real  rxnbfx,rxnbfy,rxnbfz
real  xg,yg,zg,dx,dy,dz,norm
real  ccc
real  xlpol,ylpol,zlpol,xxs,yys,zs
real  lpolx(ntot,xnpol),lpoly(ntot,ynpol),lpolz(ntot,znpol)
real  dlpolx(xnpol),dlpoly(ynpol),dlpolz(znpol)
real  mq(ntpol**2), charge
integer i,ii,jj,ij,mm,n
integer xpol,ypol,zpol,itype

!calculate q_{lm} coefficients
do ii = 1, ntpol
   coef(ii) = 0.0
enddo
do i = nsites+nfix+1, ntot
  itype = abs(typei(i))
  charge = cg(itype)
  xg = x(i)   
  yg = y(i)   
  zg = z(i)   
  if (xg.ge.xmin.and.xg.le.xmax .and.yg.ge.ymin.and.yg.le.ymax .and.zg.ge.zmin.and.zg.le.zmax) then
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
      lpolx(i,n+1)=((2.0*n-1)*xxs*lpolx(i,n)-(n-1.0)*lpolx(i,n-1))/n
    enddo
    do n = 3, ynpol-1
      lpoly(i,n+1)=((2.0*n-1)*yys*lpoly(i,n)-(n-1.0)*lpoly(i,n-1))/n
    enddo
    do n = 3, znpol-1
      lpolz(i,n+1)=((2.0*n-1)*zs*lpolz(i,n)-(n-1.0)*lpolz(i,n-1))/n
    enddo
    do ii = 1, ntpol
      xpol = lstpx(ii)
      ypol = lstpy(ii)
      zpol = lstpz(ii)
      norm = bnorm(ii)
      coef(ii) = coef(ii) + charge*norm*lpolx(i,xpol+1)*lpoly(i,ypol+1)*lpolz(i,zpol+1)
    enddo
  endif 
enddo

!construct MQ array to speed up the calculations
do ii = 1, ntpol
   mq(ii) = 0.0
   do jj = 1, ntpol
      ij = (ii-1)*ntpol+jj
      mq(ii) = mq(ii) + mmij(ij)*coef(jj)*rfscal
   enddo
enddo

!reaction field energy calculation    
egsbpb=0.0
do ii = 1, ntpol
   egsbpb = egsbpb + 0.5*coef(ii)*mq(ii)
enddo
egsbpb = egsbpb*celec

!reaction field force calculations     
if (Qforces) then 
  do mm = nsites+nfix+1, ntot
    itype = abs(typei(mm))
    charge = cg(itype)
    xg = x(mm)  
    yg = y(mm)  
    zg = z(mm)  
    if (xg.ge.xmin.and.xg.le.xmax .and.yg.ge.ymin.and.yg.le.ymax .and.zg.ge.zmin.and.zg.le.zmax) then          
      rxnbfx = 0.0
      rxnbfy = 0.0
      rxnbfz = 0.0
      ccc = charge*celec
      xxs = xscal*xg
      yys = yscal*yg
      zs = zscal*zg
      dlpolx(1) = 0.0
      dlpoly(1) = 0.0
      dlpolz(1) = 0.0
      dlpolx(2) = 1.0
      dlpoly(2) = 1.0
      dlpolz(2) = 1.0
      dlpolx(3) = 3.0*xxs
      dlpoly(3) = 3.0*yys
      dlpolz(3) = 3.0*zs
      do n = 3, xnpol-1
        dlpolx(n+1) = xxs*dlpolx(n) + n*lpolx(mm,n)
      enddo
      do n = 3, ynpol-1
        dlpoly(n+1) = yys*dlpoly(n) + n*lpoly(mm,n)
      enddo
      do n = 3, znpol-1
        dlpolz(n+1) = zs*dlpolz(n) + n*lpolz(mm,n)
      enddo
      do n = 1, xnpol
        dlpolx(n) = dlpolx(n)*xscal
      enddo
      do n = 1, ynpol
        dlpoly(n) = dlpoly(n)*yscal
      enddo
      do n = 1, znpol
        dlpolz(n) = dlpolz(n)*zscal
      enddo
      do ii = 1, ntpol
        xpol = lstpx(ii)
        ypol = lstpy(ii)
        zpol = lstpz(ii)
        norm = bnorm(ii)
        xlpol = lpolx(mm,xpol+1)
        ylpol = lpoly(mm,ypol+1)
        zlpol = lpolz(mm,zpol+1)
        dx = norm*dlpolx(xpol+1)*ylpol*zlpol
        dy = norm*xlpol*dlpoly(ypol+1)*zlpol
        dz = norm*xlpol*ylpol*dlpolz(zpol+1)
        rxnbfx = rxnbfx - dx*mq(ii)
        rxnbfy = rxnbfy - dy*mq(ii)
        rxnbfz = rxnbfz - dz*mq(ii)
      enddo
      fx(mm) = fx(mm) + rxnbfx*ccc
      fy(mm) = fy(mm) + rxnbfy*ccc
      fz(mm) = fz(mm) + rxnbfz*ccc
    endif 
  enddo
endif
return
end subroutine
