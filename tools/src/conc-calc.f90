!    CONC-CALC - Calculates re-concentration from CONC-BROMOC output files
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

program conccalc
implicit none
integer i,j,narg,kode,x,y,arg
real*8 box(3),voltot,rad1,conc,cons
real*8, parameter :: pi=3.14159265358979323846264338327950288419716939937510d0
real*8,parameter :: pi43=4d0/3d0*pi
character*(256) pname,inpfile,line
logical*1 boxt

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Input number of particles filename: ',narg,arg,inpfile)

call readarg('Spherical or Box System (sph/box): ',narg,arg,line)
boxt=.false.
if (line(1:1).eq.'b'.or.line(1:1).eq.'B') boxt=.true.
  

if (boxt) then 
  call readarg('Box X: ',narg,arg,line)
  read(line,*) box(1)

  call readarg('Box Y: ',narg,arg,line)
  read(line,*) box(2)

  call readarg('Box Z: ',narg,arg,line)
  voltot=box(1)*box(2)*box(3)
else
  call readarg('Spherical Radius: ',narg,arg,line)
  read (line,*) rad1
  voltot=pi43*rad1**3
endif

write(*,*) 'Total Volume: ',voltot

pname=inpfile(1:index(inpfile,'.',back=.true.)-1)//'-out.dat'
cons=1d0/(0.00060221412927*voltot)

open(unit=1,file=inpfile)
open(unit=2,file=pname)

read(1,*,iostat=kode) x,conc,y

do while (kode.eq.0)
  write(2,*) x,y*cons,y
  read(1,*,iostat=kode) x,conc,y
enddo
close(1)
close(2)
write(*,'(/A)') 'Normal termination of DNACDF'

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

prname='CONC-CALC'
prver='version 1.0'
prdesc='Calculates re-concentration from CONC-BROMOC output files'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 May 2012'
lastdate='07 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

