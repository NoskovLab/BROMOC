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

subroutine vdwgd0spln(xj,yj,zj,j,jtype,jtype2,Qalert)
!-----------------------------------------------------------------------
!This subroutine computes only the repulsive potental energy
!for one particle used in subroutine INTERACT in simul.f

use ioxmod
use constamod
use stdiomod
use grandmod
use nucleotmod
use gsbpmod
implicit none
integer j, jtype, jtype2
real  xj, yj, zj
logical*1 Qalert,ok
!local 
integer ncyz,ncel3,ix,iy,iz,ifir,numb
integer k,l,m,ipx,ipy,ipz
real  xi,yi,zi,ai,bi,ci,fi,esvdw,phisum,phis
real  xc,yc,zc,m3

ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz
evdwgd = 0.0
ifir = 0

ok=xj.le.xbcen2+tranx2.and.xj.ge.xbcen2-tranx2.and. &
   yj.le.ybcen2+trany2.and.yj.ge.ybcen2-trany2.and. &
   zj.le.vzmax.and.zj.ge.vzmin

if (ok) then
!  ion cartesian coordinates in the local grid system
  xi = xj + tranx2-xbcen2
  yi = yj + trany2-ybcen2
  zi = zj + tranz2-zbcen2
!        if(xi.le.0.0.or.xi.gt.2.0*tranx2) goto 101
!        if(yi.le.0.0.or.yi.gt.2.0*trany2) goto 101
!        if(zi.le.0.0.or.zi.gt.2.0*tranz2) goto 101
  esvdw = svdw
  if (Qnmcden) then
    if (Qnucl .and. j.le.nsites) then
      if (namsite(j).eq.'S ') then
        ifir = 0
      else if (namsite(j).eq.'P ') then
        ifir = ncel3
      else if (namsite(j).eq.'Ab') then
        ifir = 2*ncel3
      else if (namsite(j).eq.'Tb') then
        ifir = 3*ncel3
      else if (namsite(j).eq.'Cb') then
        ifir = 4*ncel3
      else ! namsite(j).eq.'Gb'
        ifir = 5*ncel3
      endif
    else if (Qnucl .and. Qpar .and. j.gt.(nsites+nfix)) then
      if (istrs.eq.1) numb = jtype - (inuc+1)
      if (istrs.eq.2) numb = jtype - (2*inuc+1)
      ifir = (6 + numb)*ncel3
    else if (.not.Qnucl .and. Qpar .and. j.gt.nfix) then
      ifir = (jtype-1)*ncel3
    endif
  else
    if (Qsvdw) esvdw = esvdw * scal(jtype2)
  endif
  ix=nint(xi*idcel2)
  iy=nint(yi*idcel2)
  iz=nint(zi*idcel2)
  if(ix.eq.0) ix=1
  if(iy.eq.0) iy=1
  if(iz.eq.0) iz=1
  if(ix.eq.nclx2-1)ix=nclx2-2
  if(iy.eq.ncly2-1)iy=ncly2-2
  if(iz.eq.nclz2-1)iz=nclz2-2

  phisum=0.0
  do k = ix-1, ix+1
    ipx = k*ncyz
    xc = k*dcel2
    ai = 1.5 - (xi-xc)*idcel2
    ai = m3(ai)
    if (ai.ne.0.0) then 
      do l = iy-1, iy+1
        ipy = l*nclz2 + ipx
        yc = l*dcel2
        bi = 1.5 - (yi-yc)*idcel2
        bi = m3(bi)
        if (bi.ne.0.0) then
          do m = iz-1, iz+1
            ipz = m + ipy + 1
            zc = m*dcel2
            ci = 1.5 - (zi-zc)*idcel2
            fi = ai*bi*m3(ci)
            if (fi.ne.0.0) then
              phis=phiv(ipz+ifir)
              phisum = phisum + phis
            ! Repulsive Energy
              evdwgd = evdwgd + fi*esvdw*phis
            endif
          enddo
        endif
      enddo
    endif
  enddo
  if (phisum.ge.thold27) then
    evdwgd = 1.0e10
    if (Qalert) then
      warn(jtype2)=warn(jtype2)+1
      if (Qwarn) write(outu,'(a,i5,a,5f10.5)')'Warning in routine vdwgd0spln :: particle inside membrane or protein - ',j,'  '//atnam2(jtype2),xj,yj,zj,phisum,thold27
    endif
  endif
endif
return
end
