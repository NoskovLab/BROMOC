!    DIGIT - to select numbers with particular digit precision
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

program digit
integer narg,j
real*8 a,b,i
character line*256
narg=COMMAND_ARGUMENT_COUNT()
call fullargument(line,narg)
if (narg.ge.1) then
  read(line,*) a,b
else
  read(*,*) a,b
endif
j=1
i=a/b
if ((i-int(i)).eq.0.0) j=0
write(line,'(i0)') j
write(*,'(a)') trim(adjustl(line))
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

