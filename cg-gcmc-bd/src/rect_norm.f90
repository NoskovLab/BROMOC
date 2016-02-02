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

subroutine rect_norm(ntpol,xscal,yscal,zscal,lstpx,lstpy,lstpz,bnorm)
!-----------------------------------------------------------------------
!calculate the normalization constants of the basis functions

implicit none
integer ntpol,lstpx(*),lstpy(*),lstpz(*)
real  bnorm(*),xscal,yscal,zscal
!local
real  ilxyz
integer n

ilxyz=xscal*yscal*zscal*0.125      ! inverse lxyz
do n=1,ntpol
   bnorm(n)=sqrt((2.0*lstpx(n)+1.0)*(2.0*lstpy(n)+1.0)*(2.0*lstpz(n)+1.0)*ilxyz)
enddo

return
end
