!    CHGCNTRATE - Average First N Points
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

program chgcntrate
implicit none
real*8, allocatable :: col(:),ave(:)
integer*4 i,j,k,c,n,kode,t,tt,m,narg,arg,mm
character filename*64,outfile*64,line*1024,line2*1024
integer*4 u(1024),l(1024),nu

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Filename to read: ',narg,arg,filename)
open(unit=1,file=filename,IOSTAT=kode)
i=0
m=0
line=''
read(1,'(A)',IOSTAT=kode) line
do while (kode.eq.0)
  call findparm(line,1024,nu,l,u) 
  if (nu.gt.0) i=i+1
  if (nu.gt.m) m=nu
  line=''
  read(1,'(A)',IOSTAT=kode) line
enddo
t=i
write(*,*) 'Total number of data points: ',t
write(*,*) 'Maximum number of columns: ',m

call readarg('Choose column to consider from 2 to max: ',narg,arg,line)
if (len_trim(line).eq.0) then
  mm=2
else 
  read(line,*) mm
endif

write(line2,'(A,I0,A$)') 'Number of Data to Average [',10,']:'
call readarg(trim(line2)//' ',narg,arg,line)
if (len_trim(line).eq.0) then
  n=10
else
  read(line,*) n
endif

tt=t/n+1
allocate (col(t*m),ave(tt*m))
col=0d0
ave=0d0
rewind(1)
kode=0
i=0
line=''
read(1,'(A)',IOSTAT=kode) line
do while (kode.eq.0)
  call findparm(line,1024,nu,l,u)
  if (nu.gt.0) then 
    i=i+1
    do j=0,m-1
      if (j.le.nu) read(line(l(j+1):u(j+1)),*) col(i+j*t)
    enddo
  endif
  line=''
  read(1,'(A)',IOSTAT=kode) line
enddo
close(1)
k=0
c=1
do i=1,t
  if (k==n) then
    ave(c)=col(i-1)
    do j=1,mm-1
      ave(c+j*tt)=ave(c+j*tt)/n
    enddo
    do j=mm,m-1
      ave(c+j*tt)=col((i-1)+j*t)
    enddo
    c=c+1
    k=0
  endif
  do j=1,mm-1
    ave(c+j*tt)=ave(c+j*tt)+col(i+j*t)
  enddo
  k=k+1
enddo 
ave(c)=col(t)
do j=1,mm-1
  ave(c+j*tt)=ave(c+j*tt)/k
enddo
do j=mm,m-1
  ave(c+j*tt)=col(t+j*t)
enddo

call readarg('Output filename: ',narg,arg,outfile)
open(unit=2,file=outfile)
do i=1,c
  write(2,*) (ave(i+j*tt),j=0,m-1)
enddo
close(2)
end program

subroutine findparm(str,dimn,num,llim,ulim)
implicit none
integer dimn,i,length,num
integer ulim(dimn),llim(dimn)
character ( len = dimn ) str
logical chng
length=len_trim(str)
chng=.false.
num=0
ulim(1:dimn)=0
llim(1:dimn)=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
    endif
    chng=.true.
  endif
enddo
if (num.ne.0) then
  if (ulim(num).eq.0) ulim(num)=length
endif
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

prname='CHGCNTRATE'
prver='version 2.1'
prdesc='Average First N Points'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='11 Jan 2014'
lastdate='11 Jan 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine
