!    HISTOGRAM - Generates Histograms and computes averages on second column (if present) for each bin
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

program histogram
implicit none
integer*4 i,narg,arg,nn,n,kode,col
real*8 step,maxi,mini
real*8,allocatable :: x(:),y(:),z(:),dat(:),daty(:)
character line*1024,inpfile*256,outfile*256
logical*1 norm
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
norm=.false.
! input data
call readarg('Input filename: ',narg,arg,inpfile) 
call readarg('Output filename: ',narg,arg,outfile)
call readarg('Set step: ',narg,arg,line)
read(line,*) step
call readarg('Normalize (y/n) [n]?: ',narg,arg,line)
if (line(1:1).eq.'y'.or.line(1:1).eq.'Y') norm=.true.

! process
open (unit=1,file=inpfile)
kode=0
n=0
do while(kode.eq.0)
  n=n+1
  read(1,'(A)',IOSTAT=kode) line
  if (n.eq.1.and.kode.eq.0) call countparm(line,col)
enddo
n=n-1
write(*,'(/A,I0,A,I0)') 'Columns detected: ',col,'     Lines detected: ',n
close(1)
allocate (dat(n))
if (col.gt.1) allocate (daty(n))
open (unit=1,file=inpfile)
do i=1,n
  read(1,'(A)') line
  if (col.eq.1) then 
    read(line,*) dat(i)
  else
    read(line,*) dat(i),daty(i)
  endif
enddo
close(1)
mini = int(minval(dat,dim=1)/step)*step
maxi = maxval(dat,dim=1)
nn=int((maxi-mini)/step + 0.5d0)+2
allocate (x(nn),y(nn))
if (col.gt.1) allocate (z(nn))
if (col.eq.1) then
  call histog(n,dat,step,nn,x,y,norm)
else
  call histog2(n,dat,daty,step,nn,x,y,z,norm)
endif
open (unit=1,file=outfile)
do i=1,nn
  if (col.eq.1) then 
    write(line,*) x(i),y(i)
  else
    write(line,*) x(i),y(i),z(i)
  endif
  write(1,'(A)') adjustl(trim(line))
enddo
close(1)
write(*,'(/A)') 'Normal termination.'
end program

! nn=int((max(dat)-min(dat))/step + 0.5d0)+1
subroutine histog(n,dat,step,nn,x,y,norm)
implicit none
integer*4 n,nn
real*8 dat(n),step,x(nn),y(nn)
real*8 a,b,c,e
integer*4 i,j,m
logical*1 norm

a=0d0
b=0d0
c=0d0
e=0d0
x=0d0
y=0d0
a = int(minval(dat,dim=1)/step)*step
b = maxval(dat,dim=1)
c = step
i = 0
e = a - c
do while (e.le.b) 
  i = i + 1
  e = e + c
  x(i) = e
enddo
do j = 1,n
  m = int((dat(j)-a)/c + 0.5d0) + 1
  y(m) = y(m) + 1d0
enddo
if (norm) y=y/(n*step)
end subroutine

subroutine histog2(n,dat,daty,step,nn,x,y,z,norm)
implicit none
integer*4 n,nn
real*8 dat(n),daty(n),step,x(nn),y(nn),z(nn)
real*8 a,b,c,e
integer*4 i,j,m
logical*1 norm

a=0d0
b=0d0
c=0d0
e=0d0
x=0d0
y=0d0
z=0d0
a = int(minval(dat,dim=1)/step)*step
b = maxval(dat,dim=1)
c = step
i = 0
e = a - c
do while (e.le.b)
  i = i + 1
  e = e + c
  x(i) = e
enddo
do j = 1,n
  m = int((dat(j)-a)/c + 0.5d0) + 1
  y(m) = y(m) + 1d0
  z(m) = z(m) + daty(j)
enddo
do j=1,nn
  if (y(j).ne.0d0) z(j)=z(j)/y(j)
enddo
if (norm) y=y/(n*step)
end subroutine

! count columns in a line
subroutine countparm(str,num)
implicit none
integer i,length,num
character str*(*)
logical chng
length=len_trim(str)
chng=.false.
num=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    chng=.false.
  else
    if (.not.chng) num=num+1
    chng=.true.
  endif
enddo
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
  text=''
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='HISTOGRAM'
prver='version 2.1'
prdesc='Generates Histograms and computes averages on second column (if present) for each bin'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='24 May 2008'
lastdate='07 Mar 2013'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

