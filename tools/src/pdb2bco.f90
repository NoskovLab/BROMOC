!    PDB2BCO - Converts PDB to BCO (BROMOC Coordinates)
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

module pdb
implicit none
integer nn
integer,allocatable :: resls(:),atlsf(:),atlsi(:),ras(:)
real*8,allocatable :: w(:),e(:)
integer,allocatable ::  iatom(:),ires(:),ions(:),csoli(:),csolf(:)
character*4,allocatable :: typ(:),res(:),segid(:),resid(:)
real*8,allocatable :: rt(:,:)
real*8 :: delr
integer na,nsc
character*4,allocatable :: rtc(:)
logical*1 pdbopen
end module

program pdb2bco
use pdb
implicit none
integer tnf,arg,narg
character inpfile*256,outfile*256
character bs*8

bs=repeat(achar(8),len(bs))

call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Input  PDB filename: ',narg,arg,inpfile)
call readarg('Output BCO filename: ',narg,arg,outfile)


call readpdbna(inpfile)
write(*,'(A,I0)') 'Number of Particles: ',na
call readpdb(inpfile)
call readpdbtnf(inpfile,tnf)
write(*,'(A,I0)') 'Number of frames: ',tnf
call openpdbcoor(inpfile,2)
call readpdbcoor(2)

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc
do while (pdbopen.and.nsc.le.tnf)                  ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc             ! print frame number
  call writebco(outfile)    ! writeout
  call readpdbcoor(2)
enddo
close(2)

write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine writebco(filename)
use pdb
implicit none
integer*4 i
character numchar*3,filename*(*)
write(numchar,'(I3)') nsc
do i=1,3
  if (numchar(i:i).eq.' ') numchar(i:i)='0'
enddo
open(unit=77,file=trim(filename)//trim(numchar)//'.bco')
write(77,'(i5)') na
do i=1,na
  write(77,'(i1,1x,a1,1x,a2,3(1x,f15.8))') 1,'A',typ(i),rt(1,i),rt(2,i),rt(3,i)
enddo
write(77,'(i5)') 0
close(77)
end subroutine

subroutine readpdbna(input)
use pdb
implicit none
integer*4 j,kode
character input*256
character ln*6
logical once
open(unit=1,file=input,IOSTAT=kode)
j=0
once = .true.
read(1,'(A)',IOSTAT=kode) ln
do while (kode.eq.0.or.kode.eq.64.and.once)
  if (ln(1:3).eq.'END') then
    na=j
    once = .false.
  elseif ((ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM').and.once) then
    j=j+1
  endif
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
allocate (iatom(na),ires(na),typ(na),res(na),segid(na),resid(na),rt(3,na),w(na),e(na))
end subroutine

subroutine readpdbtnf(input,tnf)
use pdb
implicit none
integer*4 j,kode,tnf
character input*256
character ln*6
open(unit=1,file=input,IOSTAT=kode)
tnf=0
j=0
read(1,'(A)',IOSTAT=kode) ln
do while (kode.eq.0.or.kode.eq.64)
  if (ln(1:3).eq.'END') tnf=tnf+1
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
end subroutine

subroutine readpdb(input)
use pdb
implicit none
integer*4 j,kode,kk
character input*256
character ln*256,ann*5,rnn*4,frmt*64
logical once,rnlimit,inlimit
frmt='(A6,I5,x,A5,A5,I4,4x,3F8.3,2F6.2,6x,A4)'
open(unit=1,file=input,IOSTAT=kode)
j=0
!i=1
once = .true.
rnlimit=.false.
inlimit=.false.
read(1,'(A)',IOSTAT=kode) ln
do while ((kode.eq.0.or.kode.eq.64).and.once)
  if (ln(1:3).eq.'END') then
!    j=0
!    i=i+1
    once=.false.
  elseif (ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM') then
    j=j+1
!    if (i.eq.1) then
      read (ln,'(6x,A5,x,A4,x,A4,x,A4,4x,3F8.3,2F6.2,6x,A4)') ann,typ(j),res(j),rnn,rt(1,j),rt(2,j),rt(3,j),e(j),w(j),segid(j)
      typ(j)=adjustl(typ(j))
      res(j)=adjustl(res(j))
      iatom(j)=j
!      if (j.eq.1) then 
!        read(ann,'(I5)') iatom(j)
!      else 
!        if (iatom(j-1).lt.99999) then 
!          read(ann,'(I5)') iatom(j)
!        else
!          read(ann,'(Z5)') iatom(j)
!        endif
!      endif
      if (.not.rnlimit.and..not.inlimit) then 
        read(rnn,'(I4)',IOSTAT=kk) ires(j)
        if (ires(j).eq.9999) inlimit=.true.
      elseif (inlimit) then
        read(rnn,'(I4)',IOSTAT=kk) ires(j)
        if (ires(j).ne.9999) then 
          inlimit=.false.
          rnlimit=.true.
          read(rnn,'(Z4)',IOSTAT=kk) ires(j)
        endif
      elseif (rnlimit) then
        read(rnn,'(Z4)',IOSTAT=kk) ires(j)
      endif
      write(resid(j),'(I4)') ires(j)
!    else
!      read (ln,'(30x,3F8.3)') r(1,j,i),r(2,j,i),r(3,j,i)
!    endif
  endif
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
end subroutine

subroutine openpdbcoor(crdfile,u)
use pdb
implicit none
integer u
character crdfile*256

open(unit=u,file=crdfile)
nsc=0
pdbopen=.true.
end subroutine

subroutine readpdbcoor(unidad)
use pdb
implicit none
integer*4 i,kode,unidad
character ln*128
logical once
i=0
once = .true.
kode=0
do while (kode.eq.0.and.once)
  read(unidad,'(A)',IOSTAT=kode) ln
  if(kode.eq.0) then
    if (ln(1:3).eq.'END') then
      nsc=nsc+1
      once=.false.
    elseif (ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM') then
      i=i+1
      read (ln,'(30x,3F8.3)') rt(1,i),rt(2,i),rt(3,i)
    endif
  endif
enddo

if (kode.ne.0) then
  close(unidad)
  pdbopen=.false.
endif
end subroutine


! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='PDB2BCO'
prver='version 1.0'
prdesc='Converts PDB to BCO (BROMOC Coordinates)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 Nov 2012'
lastdate='01 Nov 2012'

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

