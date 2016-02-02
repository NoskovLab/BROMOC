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

subroutine vdwgd0trln(xj,yj,zj,j,jtype,jtype2,Qalert)
!-----------------------------------------------------------------------
!This subroutine computes only the repulsive potential energy
!for one particle used in subroutine INTERACT in simul.f

!Repulsive forces - VDWF

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
integer ncyz,ncel3,ix,iy,iz,n1,n2,n3,in3,ifir,numb
REAL  xi,yi,zi,ai,bi,ci,fi,esvdw
real  phisum,phis

ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz
evdwgd = 0.0
ifir = 0

ok=xj.le.xbcen2+tranx2.and.xj.ge.xbcen2-tranx2.and. &
   yj.le.ybcen2+trany2.and.yj.ge.ybcen2-trany2.and. &
   zj.le.vzmax.and.zj.ge.vzmin

if (ok) then
  !ion cartesian coordinates in the local grid system              
  xi = xj + tranx2-xbcen2
  yi = yj + trany2-ybcen2
  zi = zj + tranz2-zbcen2
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
  !integer counter for ion cartesian coordinates        
  ix = int(xi*idcel2)
  iy = int(yi*idcel2)
  iz = int(zi*idcel2)
  if (ix.eq.nclx2-1) ix=nclx2-2
  if (iy.eq.ncly2-1) iy=ncly2-2
  if (iz.eq.nclz2-1) iz=nclz2-2

 !Atom charge distribution by 8 adjacent grid points

  phisum = 0.0
  do n1 = ix, ix+1
    ai = xi - n1*dcel2
    ai = 1.0 - abs(ai)*idcel2
    do n2 = iy, iy+1
      bi = yi - n2*dcel2
      bi = 1.0 - abs(bi)*idcel2
      do n3=iz,iz+1
        ci = zi - n3*dcel2
        ci = 1.0 - abs(ci)*idcel2
        fi = ai*bi*ci
        in3 = n1*ncyz + n2*nclz2 + n3 + 1
        phis=phiv(in3+ifir)
        phisum = phisum + phis
        evdwgd = evdwgd + fi*esvdw*phis
      enddo ! n3
    enddo ! n2
  enddo ! n1
  if (phisum.ge.thold8) then
    evdwgd = 1.0e10
    if (Qalert) then
      warn(jtype2)=warn(jtype2)+1
      if (Qwarn) write(outu,'(a,i5,a,5f10.5)')'Warning in routine vdwgd0trln :: particle inside membrane or protein - ',j,'  '//atnam2(jtype2),xj,yj,zj,phisum,thold8
    endif
  endif
endif   
return
end
