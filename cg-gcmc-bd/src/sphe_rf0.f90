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

subroutine sphe_rf0(xj,yj,zj,j,jtype)
! -----------------------------------------------------------------------
! this subroutine computes only the reaction field energy difference
! due to jth ion in subroutine interact in simul.f

use ioxmod
use constamod
use grandmod
use nucleotmod
use gsbpmod

implicit none
!input variables
integer j,jtype
real  xj,yj,zj
!local
integer i,ii,jj,ij,l,m,itype,lmax,mmax
real  norm,coefi,coefj,coef2(ntpol), charge, ir
real  sp,cp,st,ct,r,r2,xdiff,ydiff,zdiff,srdist2,ist
real  ar(ntot,0:24),ac(ntot,0:24),as(ntot,0:24),ap(ntot,0:24,0:24)

lmax = lstpl(ntpol)
mmax = abs(lstpm(ntpol))
srdist2 = srdist*srdist

!calculate q_{lm} coefficients for ntot
do ii = 1, ntpol
  coef(ii) = 0.0
enddo
do i = nsites+nfix+1, ntot
    if (i.eq.j) then
      itype = jtype
      xdiff = xj            
      ydiff = yj            
      zdiff = zj            
    else
      itype = abs(typei(i))
      xdiff = x(i)          
      ydiff = y(i)          
      zdiff = z(i)          
    endif
    charge = cg(itype)
    r2 = xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
    if (r2.le.srdist2) then
      r = sqrt(r2)
      ir = 1.0/r
      ct = zdiff*ir
      st = sqrt(1.0-ct*ct)
      ist=1.0/st
      cp = xdiff*ir*ist
      sp = ydiff*ir*ist
      if (r2.lt.rsmall) then                               ! in the origin
        ct = 0.0
        st = 0.0
        cp = 0.0
        sp = 0.0
      elseif (xdiff.gt.-rsmall.and.xdiff.lt.rsmall.and.ydiff.gt.-rsmall.and.ydiff.lt.rsmall) then   ! in the z-axis
        ct = 1.0
        if (zdiff.lt.0.0) ct = -1.0
        cp = 0.0
        sp = 0.0
      elseif (zdiff.gt.-rsmall.and.zdiff.lt.rsmall) then   ! in the xy plane
        ct = 0.0
        cp = xdiff*ir
        sp = ydiff*ir
      endif

      call rpowerl2(i,lmax,r,ar)           !  fill ar  (r^l   ) array
      call cosmphi2(i,mmax,cp,ac)          !  fill ac  (cos.. ) array
      call sinmphi2(i,mmax,cp,sp,as)       !  fill as  (sin.. ) array
      call alpol2(i,lmax,mmax,ct,ap)       !  fill ap  (p(lm) ) array

      do ii = 1, ntpol
        l = lstpl(ii)
        m = lstpm(ii)
        norm = bnorm(ii)
        if (l.ge.0.and.m.eq.0) then
          coef(ii) = coef(ii) + charge*norm*ar(i,l)*ap(i,l,m)
        elseif (l.gt.0.and.m.gt.0) then
          coef(ii) = coef(ii) + charge*norm*ar(i,l)*ac(i,m)*ap(i,l,m)
        elseif (l.gt.0.and.m.lt.0) then
          m = -m
          coef(ii) = coef(ii) + charge*norm*ar(i,l)*as(i,m)*ap(i,l,m)
        endif
      enddo
    endif
enddo

!calculate Q_{m} coefficients for ntot- (jth ion)
charge = cg(jtype)
do ii = 1, ntpol
  l = lstpl(ii)
  m = lstpm(ii)
  norm = bnorm(ii)
  if (l.ge.0.and.m.eq.0) then
    coef2(ii) = coef(ii) - charge*norm*ar(j,l)*ap(j,l,m)
  elseif (l.gt.0.and.m.gt.0) then
     coef2(ii) = coef(ii) - charge*norm*ar(j,l)*ac(j,m)*ap(j,l,m)
  elseif (l.gt.0.and.m.lt.0) then
      m = -m
      coef2(ii) = coef(ii) - charge*norm*ar(j,l)*as(j,m)*ap(j,l,m)
  endif
enddo

!reaction field energy calculation G(ntot)-G(ntot-1)
egsbpb = 0.0
do ii = 1, ntpol
  ij = (ii-1)*ntpol + ii
  coefi = coef(ii)
  coefj = coef2(ii)
  egsbpb = egsbpb + 0.5*(coefi*coefi-coefj*coefj)*mmij(ij)
  do jj = ii+1, ntpol
    ij = (ii-1)*ntpol + jj
    egsbpb = egsbpb + (coefi*coef(jj)-coefj*coef2(jj))*mmij(ij)
  enddo
enddo
egsbpb = egsbpb*celec

return
end subroutine
