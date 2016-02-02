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

SUBROUTINE TESTFIRST(tol,delta)
use grandmod
use constamod
use stdiomod
implicit none
integer i
!real*16 fx2(ntot), fy2(ntot), fz2(ntot)
real fx2(ntot), fy2(ntot), fz2(ntot)
real x1,y1,z1
real delta , tol
logical*1 endlog

endlog = .false.

Qforces = .false.

do i=1,ntot
  x1=x(i)
  x(i)=x1-delta
  call energy
  fx2(i)=ener
  x(i)=x1+delta
  call energy
  fx2(i)=-(ener-fx2(i))/(2*delta)
  x(i)=x1

  y1=y(i)
  y(i)=y1-delta
  call energy
  fy2(i)=ener
  y(i)=y1+delta
  call energy
  fy2(i)=-(ener-fy2(i))/(2*delta)
  y(i)=y1

  z1=z(i)
  z(i)=z1-delta
  call energy
  fz2(i)=ener
  z(i)=z1+delta
  call energy
  fz2(i)=-(ener-fz2(i))/(2*delta)
  z(i)=z1
enddo   

Qforces = .true.
call energy

do i=1,ntot
  if(abs(fx(i)-fx2(i)).gt.tol)then
     write(outu,*) '**bad x-force**    atom',i
     write(outu,*) fx(i),fx2(i)
     endlog = .true.
  endif

  if(abs(fy(i)-fy2(i)).gt.tol)then
     write(outu,*) '**bad y-force**    atom',i
     write(outu,*) fy(i),fy2(i)
     endlog = .true.
  endif

  if(abs(fz(i)-fz2(i)).gt.tol)then
     write(outu,*) '**bad z-force**    atom',i
     write(outu,*) fz(i),fz2(i)
     endlog = .true.
  endif
enddo

Qforces = .false.
if (.not.endlog) then
  write(outu,'(6x,a)')'Forces test: first derivatives ok within tolerance'
endif  
write(outu,*)

return
end
