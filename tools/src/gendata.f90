!    GENDATA - Generates Gaussian Data Points
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

program gendata
implicit none
integer arg,narg,i,n
real m,s,rang
character line*1024

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Number of data points [1000]: ',narg,arg,line)
if (len_trim(line).eq.0) then
  n=1000
else
  read(line,*) n
endif

call readarg('Media: ',narg,arg,line)
read(line,*) m

call readarg('Standard Deviation: ',narg,arg,line)
read(line,*) s

call readarg('Output filename: ',narg,arg,line)
open(unit=1,file=trim(line))
do i=1,n
  write (1,*) i,(rang()*s+m)
enddo
end program

function rang()
!use ifport
!
! two gaussian random numbers generated from two uniform random
! numbers.
!
  implicit none
  real :: x!,y
  real :: pi,r1,r2,rang
!
  pi = 4.0*atan(1.0)
  r1 = -alog(1.0-rand())
  r2 = 2.0*pi*rand()
  r1 = sqrt(2.0*r1)
  x  = r1*cos(r2)
  rang=x
!  y  = r1*sin(r2)
end function rang
!
! read arguments
subroutine readarg(ques,narg,num,text)
implicit none
integer*4 narg,num
character text*(*),ques*(*)

num=num+1
write(*,'(/A$)') ques
if (narg.ge.num) then
  call GET_COMMAND_ARGUMENT(num,text)
  write(*,'(A)') trim(text)
else
  text=''
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='GENDATA'
prver='version 2.1'
prdesc='Generates Gaussian Data Points'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='11 Jan 2014'
lastdate='11 Jan 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

