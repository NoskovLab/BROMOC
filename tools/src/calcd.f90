!    CALCD - Double Precision Calculator for BASH interface
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

program calcd
real*8 a,b,i
character c*1
read(*,*) a,c,b
if (c.eq.'+') i=a+b
if (c.eq.'-') i=a-b
if (c.eq.'*') i=a*b
if (c.eq.'%') i=a/b
if (c.eq.'^') i=a**b
write (*,*) i
end program
