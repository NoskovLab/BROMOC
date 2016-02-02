!    CORRECTBASELINE - Substract baseline defined in input parameters
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

program correctbaseline
implicit none
integer arg,narg,ftype,kode
character line*256,infile*256,outfile*256
real*8 a,b,c,y,x,bl

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Input filename: ',narg,arg,infile)
outfile=infile(1:index(infile,'.',back=.true.)-1)//'-bl'//infile(index(infile,'.',back=.true.):len_trim(infile))

call readarg('Function type (1/2): ',narg,arg,line)
read(line,*) ftype

call readarg('a = ',narg,arg,line)
read(line,*) a
call readarg('b = ',narg,arg,line)
read(line,*) b
call readarg('c = ',narg,arg,line)
read(line,*) c

open(unit=1,file=infile)
open(unit=2,file=outfile)
read(1,*,iostat=kode) x,y
do while (kode.eq.0)
  if (ftype.eq.1) then
    bl=a*exp(-b*x)+c
  elseif (ftype.eq.2) then
    bl=a*x**2+b*x+c
  else
    stop 'Incorrect function type'
  endif
  write(2,*) x,y-bl
  read(1,*,iostat=kode) x,y
enddo
close(1)
close(2)

write(*,*) 'Normal End'
end program

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

prname='CORRECTBASELINE'
prver='version 1.0'
prdesc='Substract baseline defined in input parameters'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='05 Jun 2013'
lastdate='05 Jun 2013'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

