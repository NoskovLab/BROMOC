!    TRANSCRD - Moves and rotates coordinates in CHARMM crd file
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

program transcrd
implicit none
real*8 x,y,z,dx,dy,dz,cx,cy,cz,rot(3,3),r(3)
integer narg,kode,n,arg
character inpfile*256,line*256

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

! Read input filename from argument or keyboard
call readarg('Input CHARMM CRD/COR filename: ',narg,arg,inpfile)
open(unit=1,file=inpfile,iostat=kode)

! Read output filename from argument or keyboard
call readarg('Output CHARMM CRD/COR filename: ',narg,arg,inpfile)
open(unit=2,file=inpfile)

! compute geometric center
n=0
cx=0d0
cy=0d0
cz=0d0
read(1,'(A)',iostat=kode) line
do while (kode.eq.0)
  if (line(1:1).ne.'*') then
    if (line(25:25).eq.'.') then
      n=n+1
      read(line(21:30),*) x
      read(line(31:40),*) y
      read(line(41:50),*) z
      cx=cx+x
      cy=cy+y
      cz=cz+z
    endif
  endif
  read(1,'(A)',iostat=kode) line
enddo
rewind(1)
cx=cx/n
cy=cy/n
cz=cz/n
write(line,*) 'Centroid: ',cx,cy,cz
write(*,'(a)') trim(adjustl(line))

! read displacement
call readarg('Displacement (x y z): ',narg,arg,line)
write(*,'(/A$)') 'Displacement (x y z): '
read(line,*) dx,dy,dz

call readarg('Rotational Matrix (xx yx zx xy yy zy xz yz zz): ',narg,arg,line)
read(line,*) rot(1,1),rot(2,1),rot(3,1),rot(1,2),rot(2,2),rot(3,2),rot(1,3),rot(2,3),rot(3,3)

read(1,'(A)',iostat=kode) line
do while (kode.eq.0)
  if (line(1:1).eq.'*') then
    write(2,'(A)') trim(line)
  else
    if (line(25:25).ne.'.') then 
      write(2,'(A)') trim(line)
    else
      read(line(21:30),*) x
      read(line(31:40),*) y
      read(line(41:50),*) z
      r(1)=x-dx
      r(2)=y-dy
      r(3)=z-dz
      r=matmul(r,transpose(rot))
      write(2,'(A,3F10.5,A)') line(1:20),r(1),r(2),r(3),line(51:len_trim(line))
    endif
  endif
  read(1,'(A)',iostat=kode) line
enddo
  
close(1)
close(2)

write(*,'(/A/)') 'Normal Termination.'

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

prname='TRANSCRD'
prver='version 1.1'
prdesc='Moves and rotates coordinates in crd file'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 Jan 2012'
lastdate='11 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine
