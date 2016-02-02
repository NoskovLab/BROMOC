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

!
!-----------------------------------------------------------------------
!
subroutine lualloc(u)
!
!.....lualloc - allocate a unique logical*1 unit number for i/o.
!
!     The unit number is returned through lualloc and u simultaneously.
!
use  stdiomod
use  strtoolsmod
use  ioxmod
use  errormod
implicit none
integer          u
!
!.....skip allocated units
!
u = 1
do while (u.le.maxopen.and.alloc(u)) 
  u = u + 1
enddo

if (u.gt.maxopen) then
   call error ('lualloc', 'exceded open file limit'//eol,faterr)
else
   alloc(u) = .true.
endif
return
end
