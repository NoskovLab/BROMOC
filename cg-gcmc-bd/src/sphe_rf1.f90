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

subroutine sphe_rf1
!-----------------------------------------------------------------------
!calculate the reaction field energy and forces on each ions

use ioxmod
use constamod
use grandmod
use nucleotmod
use gsbpmod

implicit none
!local variables
integer i,ii,jj,ij,l,m,mm,itype,lmax,mmax
real  rxnbfx,rxnbfy,rxnbfz
real  norm, charge
real  ccc,rpl,cmp,smp,apl,ir,ist
real  sp,cp,st,ct,r,r2,xdiff,ydiff,zdiff,srdist2
real  ar(ntot,0:24),ac(ntot,0:24),as(ntot,0:24)
real  ap(ntot,0:24,0:24),adp(ntot,0:24,0:24)
real  mq(ntot)
real  spj(ntot),cpj(ntot),stj(ntot),ctj(ntot)
real  dx,dy,dz,dr,ddt,dp

lmax = lstpl(ntpol)
mmax = abs(lstpm(ntpol))
srdist2 = srdist*srdist

!calculate Q_{lm} coefficients
do ii =1,ntpol
  coef(ii) = 0.0
enddo
do i = nsites+nfix+1, ntot
  itype = abs(typei(i))
  charge = cg(itype)
  xdiff = x(i) 
  ydiff = y(i) 
  zdiff = z(i) 
  r2 = xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
  if (r2.le.srdist2) then
    r = sqrt(r2)
    ir = 1.0/r
    ct = zdiff*ir
    st = sqrt(1.0-ct*ct)
    ist=10/st
    cp = xdiff*ir*ist
    sp = ydiff*ir*ist
    if (r2.lt.rsmall) then  ! in the origin                             
      ct = 0.0
      st = 0.0
      cp = 0.0
      sp = 0.0
    elseif (xdiff.gt.-rsmall.and.xdiff.lt.rsmall.and.ydiff.gt.-rsmall.and.ydiff.lt.rsmall) then  ! in the z-axis
      ct = 1.0
      if( zdiff.lt.0.0) ct=-1.0
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
    call dalpol2(i,lmax,mmax,ct,ap,adp)  !  fill adp (dp(lm)) array
    spj(i) = sp
    cpj(i) = cp
    stj(i) = st
    ctj(i) = ct

    do ii=1,ntpol
      l = lstpl(ii)
      m = lstpm(ii)
      norm = bnorm(ii)
      if (l.ge.0.and.m.eq.0) then
        coef(ii) = coef(ii)+charge*norm*ar(i,l)*ap(i,l,m)
      elseif (l.gt.0.and.m.gt.0) then
        coef(ii) = coef(ii)+charge*norm*ar(i,l)*ac(i,m)*ap(i,l,m)
      elseif (l.gt.0.and.m.lt.0) then
        m = -m
        coef(ii) = coef(ii)+charge*norm*ar(i,l)*as(i,m)*ap(i,l,m)
      endif
    enddo
  endif 
enddo

!construct mq array to speed up the calcaulations
do ii = 1, ntpol
  mq(ii)=0.0
  do jj = 1, ntpol
    ij = (ii-1)*ntpol+jj
    mq(ii) = mq(ii)+mmij(ij)*coef(jj)
  enddo
enddo

!reaction field energy calculation    
egsbpb = 0.0
do ii = 1, ntpol
  egsbpb = egsbpb+0.5*coef(ii)*mq(ii)
enddo
egsbpb = egsbpb*celec

!reaction field force calculations     
if (Qforces) then
  do mm = nsites+nfix+1, ntot
    itype = abs(typei(mm))
    charge = cg(itype)
    xdiff = x(mm) 
    ydiff = y(mm) 
    zdiff = z(mm) 
    r2 = xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
    if (r2.le.srdist2) then
      rxnbfx = 0.0
      rxnbfy = 0.0
      rxnbfz = 0.0
      ccc=charge*celec
      sp = spj(mm)
      cp = cpj(mm)
      st = stj(mm)
      ct = ctj(mm)
      do ii = 1, ntpol
        l = lstpl(ii)
        m = lstpm(ii)
        norm = bnorm(ii)
        if (m.eq.0) then
          if (l.eq.0) then
            dr = 0.0
            ddt = 0.0
            dp = 0.0
          else
            rpl = ar(mm,l-1)
            dr = l*rpl*ap(mm,l,m)
            ddt = -rpl*adp(mm,l,m)*st
            dp = 0.0
          endif
        elseif( m.gt.0) then
          rpl = ar(mm,l-1)
          cmp = ac(mm,m)
          apl = ap(mm,l,m)
          dr = l*rpl*cmp*apl
          ddt = -rpl*cmp*adp(mm,l,m)*st
          dp = -rpl*m*as(mm,m)*apl*ist
          if (st.eq.0.0) dp = 0.0
        elseif (m.lt.0) then
          m = -m
          rpl = ar(mm,l-1)
          smp = as(mm,m)
          apl = ap(mm,l,m)
          dr = l*rpl*smp*apl
          ddt = -rpl*smp*adp(mm,l,m)*st
          dp = rpl*m*ac(mm,m)*apl*ist
          if (st.eq.0.0) dp = 0.0
        endif
        dx = norm*(dr*st*cp+ddt*ct*cp-dp*sp)
        dy = norm*(dr*st*sp+ddt*ct*sp+dp*cp)
        dz = norm*(dr*ct   -ddt*st         )
        rxnbfx = rxnbfx-dx*mq(ii)
        rxnbfy = rxnbfy-dy*mq(ii)
        rxnbfz = rxnbfz-dz*mq(ii)
      enddo
      fx(mm) = fx(mm)+rxnbfx*ccc
      fy(mm) = fy(mm)+rxnbfy*ccc
      fz(mm) = fz(mm)+rxnbfz*ccc
    endif
  enddo
endif
return
end subroutine
