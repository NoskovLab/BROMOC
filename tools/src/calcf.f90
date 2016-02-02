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

program calcreal
integer narg
real*8 a,b,i
character s*1,line*256
narg=COMMAND_ARGUMENT_COUNT()
call fullargument(line,narg)
if (narg.ge.1) then
  read(line,*) a,s,b
else
  read(*,*) a,s,b
endif

if (s.eq.'+') i=a+b
if (s.eq.'-') i=a-b
if (s.eq.'x') i=a*b
if (s.eq.'%') i=a/b
if (s.eq.'^') i=a**b
write (line,'(f30.14)') i
write (*,'(a)') trim(adjustl(line))

contains
! Cross Product: w = u x v
  function cross_product(u,v)
  implicit none
  real*8 u(3),v(3),w(3),cross_product(3)
  w(1)=u(2)*v(3)-u(3)*v(2)
  w(2)=u(3)*v(1)-u(1)*v(3)
  w(3)=u(1)*v(2)-u(2)*v(1)
  cross_product=w
  end function

end program

subroutine fullargument(line,narg)
implicit none
integer narg,i
character*256 line,argm*256
line=''
do i=1,narg
  call GET_COMMAND_ARGUMENT(i,argm)
  line=trim(line)//' '//trim(argm)
enddo
end subroutine

