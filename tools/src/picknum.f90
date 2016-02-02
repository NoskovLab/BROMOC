!    PICKNUM - Generates Random Natural Numbers. Uses CPU clock to generate randomness.
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

program picknum
!use ifport
implicit none
integer*4 arg,narg,nn,i,iseed, maxn
character line*256
logical*1 oput
arg=0
narg=COMMAND_ARGUMENT_COUNT()
oput=.true.
if (narg.eq.3) oput=.false. 
if (oput) call header()
iseed=0
call askclock(iseed)
write(line,'(I0)') iseed
call readarg('Seed ['//trim(line)//']: ',narg,arg,line,oput)
if (len_trim(line).gt.0.and.trim(line).ne.'0') read(line,*) iseed
call srand(iseed)

call readarg('Number of numbers [1]: ',narg,arg,line,oput)
if (len_trim(line).eq.0) line='1'
read(line,*) nn

call readarg('Max Integer Number [1000]: ',narg,arg,line,oput)
if (len_trim(line).eq.0) line='1000'
read(line,*) maxn

do i=1,nn
  write(*,*) int(rand()*maxn)+1
enddo

if (oput) write(*,'(/A/)') 'Normal Termination.'

end program

subroutine askclock(iseed)
implicit none
integer*4 iseed
if (iseed.le.0.or.iseed.ge.2147483647) call system_clock(count=iseed)
end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='PICKNUM'
prver='version 1.1'
prdesc='Generates Random Natural Numbers. Uses CPU clock to generate randomness.'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='20 Feb 2013'
lastdate='19 Mar 2014'

write(*,'(/A)') trim(prname)//' '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A/)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
end subroutine

! read arguments
subroutine readarg(ques,narg,num,text,oput)
implicit none
integer*4 narg,num
character text*(*),ques*(*)
logical*1 oput

num=num+1
if (oput) write(*,'(/A$)') ques
if (narg.ge.num) then
  call GET_COMMAND_ARGUMENT(num,text)
  if (oput) write(*,'(A)') trim(text)
else
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine


