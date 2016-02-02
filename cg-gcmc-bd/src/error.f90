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

subroutine error (routine, message, errortype)

!.....error - proper handling of error messages.

use             stdiomod
use             errormod
implicit none

character*(*)       routine, message
integer             errortype

if (errortype.eq.faterr) then
   write (luwrite,100) 'fatal error', routine, message
   stop 1
else if (errortype.eq.warning) then
   write (luwrite,100) 'warning', routine, message
else if (errortype.eq.noerr) then
   write (luwrite,100) 'no error', routine, message
else
   write (luwrite,100) 'unknown type of error', routine, message
endif

100 format (/1x,'*** '/1x,'*** ',a,' at ',a,' :'/1x,'*** ',a/1x,'***')

end
