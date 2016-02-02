!    STATS - Calculates media and standard deviation of a column data
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

program stats
implicit none
integer i,j,k,kode,c,narg,arg
double precision r,rm,rms,sr,ssr
character inputfile*128,line*128

call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()
call readarg('Input filename: ',narg,arg,inputfile)
call readarg('Column number: ',narg,arg,line)
read(line,*) c

open(unit=1,file=inputfile,IOSTAT=kode)
i=-1
do while (kode.eq.0)
  r=0
  read(1,*,IOSTAT=kode) (r,j=1,c)
  i=i+1
  sr=sr+r
  ssr=ssr+r**2
enddo
rm=sr/i
rms=(ssr/i-rm**2)
close(1)
write(*,*) 'Data points:',i
write(*,*) 'Media:',rm
write(*,*) '(True) Standard Deviation: ',sqrt(rms)
write(*,*) 'True Standard Deviation of the sample-mean: ',sqrt(rms/i)
write(*,*) 'Sample Standard Deviation: ',sqrt(rms*i/(i-1))
write(*,*) 'Standard Deviation/Error of the sample-mean: ',sqrt(rms/(i-1))
end program

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='STATS'
prver='version 1.5'
prdesc='Calculates media and standard deviation of a column data'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='03 Jul 2007'
lastdate='28 Sep 2012'

write(*,'(/A)') trim(prname)//' '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A/)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
end subroutine

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
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

