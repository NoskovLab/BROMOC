!    STT - Calculates media and standard deviation of a column data
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

program stt
implicit none
integer i,kode,narg,arg
double precision r,rm,rms,sr,ssr
character line*1024

arg=0
narg=COMMAND_ARGUMENT_COUNT()
do while (arg.lt.narg)
  call readarg(narg,arg,line)
  if (line.eq.'-h') then
    call header()
    write(*,*) 'USAGE:     cat file | stt [option]'
    write(*,*)
    write(*,*) '-h            Help'
    write(*,*) '-v            Version'
    write(*,*) '-n            Number of data points'
    write(*,*) '-m     -a     Media/Average'
    write(*,*) '-tsd   -sd    True Standard Deviation'
    write(*,*) '-tsdm         True Standard Deviation of the sample-mean'
    write(*,*) '-ssd          Sample Standard Deviation'
    write(*,*) '-ssdm  -se    Sample Standard Deviation of the sample-mean (Standard Error)'
    return
  else if (line.eq.'-v') then
    call header()
    return
  endif
enddo
arg=0
i=0
sr=0.0
ssr=0.0
read(5,*,IOSTAT=kode) r
do while (kode.eq.0)
  i=i+1
  sr=sr+r
  ssr=ssr+r**2
  read(5,*,IOSTAT=kode) r
enddo
rm=sr/i
rms=(ssr/i-rm**2)

do while (arg.lt.narg)
  call readarg(narg,arg,line)
  if (line.eq.'-n') then
   write(*,*) i
  else if (line.eq.'-m'.or.line.eq.'-a') then
    write(*,*) rm
  else if (line.eq.'-sd'.or.line.eq.'-tsd') then
    write(*,*) sqrt(rms)
  else if (line.eq.'-tsdm') then
    write(*,*) sqrt(rms/i)
  else if (line.eq.'-ssd') then
    write(*,*) sqrt(rms*i/(i-1))
  else if (line.eq.'-ssdm'.or.line.eq.'-se') then
    write(*,*) sqrt(rms/(i-1))
  endif
enddo

end program

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='STT'
prver='version 1.0'
prdesc='Calculates media and standard deviation of a column data'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='03 Jul 2007'
lastdate='17 Dec 2014'

write(*,'(/A)') trim(prname)//' '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A/)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
end subroutine

subroutine readarg(narg,num,text)
implicit none
integer*4 narg,num
character text*(*)

num=num+1
if (narg.ge.num) then
  call GET_COMMAND_ARGUMENT(num,text)
endif
text=adjustl(trim(text))
end subroutine

