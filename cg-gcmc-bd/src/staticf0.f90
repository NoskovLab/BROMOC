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

subroutine staticf0(xj,yj,zj,j,jtype)
!-----------------------------------------------------------------------
!This subroutine computes only the external static field energy
!for one particle used in subroutine INTERACT in simul.f   

use constamod
use grandmod
use nucleotmod
use gsbpmod
implicit none
!input variables
integer j, jtype
real  xj, yj, zj
!local variables
integer ncyz,ix,iy,iz,n1,n2,n3,in3
real  chi,xi,yi,zi,ai,bi,ci,fi
logical*1 Beion, Besite, ok

egsbpa = 0.0
Beion = Qpar .and. j.gt.nsites
Besite = .false.
if (Qnucl) then 
  if (j.le.nsites) Besite = namsite(j).eq.'P '
endif
ok = Beion .or. Besite

if (Beion) chi = cg(jtype)
if (Besite) chi = cgnuc
if (ok) ok=xj.le.xbcen1+tranx1.and.xj.ge.xbcen1-tranx1.and. &
           yj.le.ybcen1+trany1.and.yj.ge.ybcen1-trany1.and. &
           zj.le.zbcen1+tranz1.and.zj.ge.zbcen1-tranz1
if (ok) then
  ncyz = ncly1*nclz1
!  ion cartesian coordinates in the local grid system              
  xi = xj + tranx1-xbcen1
  yi = yj + trany1-ybcen1
  zi = zj + tranz1-zbcen1
!  integer*4 counter for ion cartesian coordinates        
  ix = int(xi*idcel1)
  iy = int(yi*idcel1)
  iz = int(zi*idcel1)
  if (ix.eq.nclx1-1) ix=nclx1-2
  if (iy.eq.ncly1-1) iy=ncly1-2
  if (iz.eq.nclz1-1) iz=nclz1-2

!Atom charge distribution by 8 adjacent grid points

  do n1 = ix, ix+1
    ai = xi - n1*dcel1
    ai = 1.0 - abs(ai)*idcel1
    do n2 = iy, iy+1
      bi = yi - n2*dcel1
      bi = 1.0 - abs(bi)*idcel1
      do n3 = iz, iz+1
        ci = zi - n3*dcel1
        ci = 1.0 - abs(ci)*idcel1
        fi = ai*bi*ci
        in3 = n1*ncyz + n2*nclz1 + n3 + 1
!Electrostatic Energy 
        egsbpa = egsbpa + (fi*chi*phix(in3)*celec)
      enddo ! n3
    enddo ! n2
  enddo ! n1
endif ! ok

return
end
