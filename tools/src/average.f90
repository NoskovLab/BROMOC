!    AVERAGE - Reads text files and computes average of the second column of each file
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

program average
implicit none
real*8,allocatable :: x(:),y(:),yy(:),s(:)
real*8 xt,yt,ina,ina1
integer arg,narg,i,j,n,kode,na,na1
character inpfile*256,line*256,bs*10

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

na=narg-1
na1=na-1
ina=1d0/dfloat(na)
ina1=1d0/dfloat(na1)
if (narg.ge.3) then
  ! open first file
  call readarg('',narg,arg,inpfile)
  open(unit=1,file=inpfile,iostat=kode)

  ! count lines
  read(1,'(A)',iostat=kode) line
  n=0
  do while (kode.eq.0)
    n=n+1
    read(1,'(A)',iostat=kode) line
  enddo
  
  ! start from first line
  rewind(1)

  ! allocate arrays
  allocate (x(n),y(n),yy(n),s(n))

  write(*,'(A,I0/)') 'Number of lines: ',n
  write(*,'(A,I0,A,A$)') 'File',1,': ',trim(inpfile)

  ! read first file data
  write(*,'(I10$)') 1
  do i=1,n
    read(1,*) x(i),y(i)
    yy(i)=y(i)**2
    write(*,'(A10,I10$)') bs,i
  enddo
  close(1)
  write(*,'(/)')

  ! read and sum the rest of the files
  do i=2,narg-1
    call readarg('',narg,arg,inpfile)
    open(unit=1,file=inpfile,iostat=kode)
    write(*,'(A,I0,A,A$)') 'File',i,': ',trim(inpfile)
    write(*,'(I10$)') 1
    do j=1,n
      read(1,*) xt,yt
      y(j)=y(j)+yt  ! addition
      yy(j)=yy(j)+yt**2
      write(*,'(A10,I10$)') bs,j
    enddo
    close(1)
    write(*,'(/)')
  enddo
else
  write(*,'(/A/)') 'USAGE: AVERAGE inpfile1 inpfile2 ... outfile'
  stop
endif

! average
s=yy
do i=1,n
  yy(i)=s(i)-y(i)**2*ina
  y(i)=y(i)*ina
enddo

! write output
call readarg('',narg,arg,inpfile)
open(unit=1,file=inpfile)
do i=1,n
  write(1,*) x(i),y(i),dsqrt(yy(i)*ina*ina1),dsqrt(yy(i)*ina1),dsqrt(s(i)*ina1)
enddo
close(1)

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

prname='AVERAGE'
prver='version 1.2'
prdesc='Reads text files and computes average of the second column of each file'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='24 Jan 2012'
lastdate='11 Jan 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
!write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

