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

subroutine readmmij(unit,outu)
!-----------------------------------------------------------------------
! read INPUT matrix M* MMIJ for a spherical inner region
!
use ioxmod
use gsbpmod

implicit none
integer unit,outu
integer i,j,ij,ji,ntpol2


read(unit) shapes

if (shapes.eq.'RECTBOX ') then
  read(unit) xnpol,ynpol,znpol
  ntpol = xnpol*ynpol*znpol
  ntpol2 = ntpol*ntpol
  allocate(lstpl(ntpol),lstpm(ntpol),lstpx(ntpol),lstpy(ntpol),lstpz(ntpol),bnorm(ntpol),coef(ntpol),mmij(ntpol2))
  read(unit) (lstpx(i),i=1,ntpol)
  read(unit) (lstpy(i),i=1,ntpol)
  read(unit) (lstpz(i),i=1,ntpol)
  read(unit) (mmij(i), i=1,ntpol2)
  do i = 1, ntpol
    do j = 1, ntpol
      if (i.gt.j) then
        ij = (i-1)*ntpol + j
        ji = (j-1)*ntpol + i
        mmij(ij) = mmij(ji)
      endif
     enddo
  enddo        
 !Writting in output file          
  write(outu,101) 'Number of Legendre Polynomials in X  (XNPOL) = ',xnpol
  write(outu,101) 'Number of Legendre Polynomials in Y  (YNPOL) = ',ynpol
  write(outu,101) 'Number of Legendre Polynomials in Z  (ZNPOL) = ',znpol

elseif(shapes.eq.'SPHERE  '.or.shapes.eq.'NSPHERE ') then
  shapes = 'SPHERE  '
  read(unit) nmpol
  ntpol = nmpol*nmpol
  ntpol2 = ntpol*ntpol
  read(unit) (mmij(i),i=1,ntpol2)
  call sphe_svpol(nmpol,lstpl,lstpm)
  do i = 1, ntpol
    do j = 1, ntpol
      if (i.gt.j) then
        ij = (i-1)*ntpol + j
        ji = (j-1)*ntpol + i
        mmij(ij) = mmij(ji)
      endif  
    enddo 
  enddo  
!Writting in output file          
  write(outu,101) 'Number of Multipoles                 (NMPOL) = ',nmpol

endif
!formats
101  format(6x,a,i6)

return
end
