!    NOISEFILTER - Reduces noise from a given data set by averaging N contiguous points
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

program noisefilter
implicit none
real*8, allocatable :: col(:),ave(:)
real*8 av,tol,dif,mx,mn
integer*4 i,j,c,n,kode,t,ll,hl,d,m,narg,arg
character filename*64,outfile*64,line*1024
integer*4 u(1024),l(1024),nu
logical ig0,rtr

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Filename to read: ',narg,arg,filename)
open(unit=1,file=filename,IOSTAT=kode)
i=0
m=0
do while (kode.eq.0)
  line=''
  read(1,'(A)',IOSTAT=kode) line
  call findparm(line,1024,nu,l,u) 
  if (nu.gt.0) i=i+1
  if (nu.gt.m) m=nu
enddo
t=i
write(*,*) 'Total number of data points: ',t
write(*,*) 'Maximum number of columns: ',m
call readarg('Column Number to process: ',narg,arg,line)
read(line,*) c
if (c.gt.m.or.c.le.0) then 
  write(*,*) 'Selected Column not present using the last one'
  c=m
endif
allocate (col(t*m),ave(t))
col=0d0
rewind(1)
kode=0
i=0
do while (kode.eq.0)
  line=''
  read(1,'(A)',IOSTAT=kode) line
  call findparm(line,1024,nu,l,u)
  if (nu.gt.0) then 
    i=i+1
    do j=1,m
      if (j.le.nu) read(line(l(j):u(j)),*) col(i+(j-1)*t)
    enddo
  endif
enddo
close(1)
mn=col(1+(c-1)*t)
mx=col(1+(c-1)*t)
do i=2,t
  if (col(i+(c-1)*t).lt.mn) mn=col(i+(c-1)*t)
  if (col(i+(c-1)*t).gt.mx) mx=col(i+(c-1)*t)
enddo
tol=(mx-mn)*0.015d0

write(*,'(A,I0,A$)') 'Number of Data to Average [',2*int(t/80)+1,']: '
call readarg('',narg,arg,line)
if (len_trim(line).eq.0) then 
  n=2*int(t/80)+1
else
  read(line,*) n
endif

call readarg('Activate restrictions (y/n) [n]?: ',narg,arg,line)
rtr=.false.
line=trim(line)
if (line(1:1).eq.'y'.or.line(1:1).eq.'Y') rtr=.true.

if (rtr) then 
  write(*,'(A,F20.10,A$)') 'Tolerance [',tol,']: '
  call readarg('',narg,arg,line)
  if (len_trim(line).gt.0) then 
    read(line,*) tol
  endif
  call readarg('Ignore null data points (y/n) [y]? ',narg,arg,line)
  if (line(1:1)=='n'.or.line(1:1)=='N') then
    ig0=.false.
  else
    ig0=.true.
  endif
endif

! Convert to EVENs to the higher ODD neighbour
n=2*int(n/2)+1
! Do my thing
do i=1,t
  if (rtr.and.ig0.and.col(i+(c-1)*t).le.tol) then
    av=col(i+(c-1)*t)
  else
    ll=2*(i-1)+1
    hl=2*(t-i)+1
    d=n
    if (ll.lt.d) d=ll
    if (hl.lt.d) d=hl
    av=0d0
    do j=i-(d-1)/2,i+(d-1)/2
      av=av+col(j+(c-1)*t)
    enddo
    av=av/d
    if (rtr) then
      dif=dabs(av-col(i+(c-1)*t))
      do while (dif.gt.tol)
        if (d.gt.3) then  
          av=0
          d=d-2
          do j=i-(d-1)/2,i+(d-1)/2
            av=av+col(j+(c-1)*t)
          enddo
          av=av/d
          dif=dabs(av-col(i+(c-1)*t))
        else
          av=col(i+(c-1)*t)+dsign(tol,av-col(i+(c-1)*t))
          dif=tol
        endif
      enddo
    endif
  endif
  ave(i)=av
enddo
call readarg('Output filename: ',narg,arg,outfile)
open(unit=2,file=outfile)
do i=1,t
  write(2,*) (col(i+(j-1)*t),j=1,c-1),ave(i),(col(i+(j-1)*t),j=c+1,m)
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
if (ulim(num).eq.0) ulim(num)=length
end subroutine

! read arguments
subroutine readarg(ques,narg,num,text)
implicit none
integer*4 narg,num
character text*(*),ques*(*)

num=num+1
write(*,'(A$)') ques
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

prname='NOISEFILTER'
prver='version 2.2'
prdesc='Reduces noise from a given data set by averaging N contiguous points'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='05 Feb 2008'
lastdate='03 Dec 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine
