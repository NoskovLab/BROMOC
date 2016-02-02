!    BTREXTRACT - Extract a frame from BROMOC Trajectory (.btr) to another btr
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

module comun
implicit none
integer*4,allocatable :: typ(:),maxntyp(:)
integer*4 nsc,ntop,itype,nframe,maxntop,ir,lr
real*8 runtime
logical*1 inpopen
character*4, allocatable :: atnam2(:)
character bs*8
real*4,allocatable :: rt(:,:)
real*8 :: xtlabc6(6)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title

end module

program btrextract
use comun
implicit none
real*8 aa,bb
integer,allocatable :: fac(:)
integer h,i,j,k,l,m,a,b,c,d,narg,arg,ll(1000000),ul(1000000),n
character inpfile*256,line*1000000,outfile*256
logical*1 loopon

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)
call readarg('Output BROMOC trajectory (.btr) filename: ',narg,arg,outfile)
call readarg('Choose frame to save: ',narg,arg,line)
call findparm(trim(line),n,ll,ul)
allocate (fac(n))
do i=1,n
  fac(i)=chr2int(getparm(trim(line),n,ll,ul,i))
enddo

! Read gcmcbd input file
call readgcmcbdhead(inpfile,.false.)
call readgcmcbdbody()
call writebtrhead(trim(outfile),1,n,itype,atnam2)

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc
do while (inpopen)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  do i=1,n
    if (nsc.eq.fac(i)) then
      call writebtrbody(1,runtime,ntop,typ,rt)
    endif
  enddo
  if (inpopen) call readgcmcbdbody()          ! read next fram
enddo
close(1) 

write(*,'(/A)') 'Normal termination'
contains
  ! get parameter pn from line
  function getparm(str,num,llim,ulim,pn)
  implicit none
  integer num,pn
  integer llim(*),ulim(*)
  character getparm*(ulim(pn)-llim(pn)+1),str*(*)
  if (pn.gt.num.or.pn.lt.1.or.num.lt.1) then
    getparm=''
  else
    getparm=str(llim(pn):ulim(pn))
  endif
  end function
! char to integer
  function chr2int(str)
  implicit none
  integer chr2int,kode
  character str*(*)
  read(str,*,iostat=kode) chr2int
  if (kode.ne.0) then
    write(*,*) trim(str)
    stop 'Not an integer'
  endif
  end function

end program

! find parameters
subroutine findparm(str,num,llim,ulim)
implicit none
integer i,length,num
integer ulim(*),llim(*)
character str*(*)
logical chng
length=len_trim(str)
chng=.false.
num=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
      ulim(num)=0
    endif
    chng=.true.
  endif
enddo
if (ulim(num).eq.0) ulim(num)=length
end subroutine

subroutine readgcmcbdhead(inpfile,no)
use comun
implicit none
integer i, j, kode
integer*4 recordi,recordl
character*256 inpfile
logical no

open(14,file=inpfile,form='unformatted',iostat=kode,position='rewind',access='stream')
read(14,iostat=kode) recordi
read(14,iostat=kode) nframe                 ! number of frames
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

read(14,iostat=kode) recordi
read(14,iostat=kode) itype                  ! number of ions and nucleotides
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (.not.allocated(atnam2)) allocate (atnam2(itype))
if (.not.allocated(maxntyp)) allocate (maxntyp(itype))

read(14,iostat=kode) recordi
read(14,iostat=kode) (atnam2(j),j=1,itype)  ! ion and nucleotides types in char
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (.not.no) then
  write(*,'(A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(*,*) 'Fragment types: ',(atnam2(j),j=1,itype)
endif
nsc=0
if (kode.eq.0) then
  inpopen=.true.
else
  close(14)
endif
end subroutine

! check fortran binary record consistency
subroutine check(recordi,recordl,un,kode)
implicit none
integer*4 recordi,recordl,pos,un,kode
character*64 fname

if (kode.eq.0) then
  if (recordi.ne.recordl) then
    inquire(unit=un,pos=pos,name=fname)
    write (*,*)
    write (*,*) 'Unit: ',un
    write (*,*) 'Filename: ',trim(fname)
    write (*,*) 'Record not matching (init, last): ',recordi,recordl
    write (*,*) 'Position of last byte read: ',pos
    stop
  endif
endif
end subroutine

subroutine readgcmcbdbody()
use comun
implicit none
integer kode,i
integer*4 recordi,recordl

inquire(unit=14,pos=ir)

if (allocated(typ)) deallocate (typ)
if (allocated(rt)) deallocate (rt)
read(14,iostat=kode) recordi
read(14,iostat=kode) runtime !simulation time for each step (ns)
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

read(14,iostat=kode) recordi
read(14,iostat=kode) ntop ! number of particles
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (kode.eq.0) then 
  allocate (typ(ntop),rt(3,ntop))

  read(14,iostat=kode) recordi
  read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
  
  read(14,iostat=kode) recordi
  read(14,iostat=kode) (rt(1,i),i=1,ntop)
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
  
  read(14,iostat=kode) recordi
  read(14,iostat=kode) (rt(2,i),i=1,ntop)
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
  
  read(14,iostat=kode) recordi
  read(14,iostat=kode) (rt(3,i),i=1,ntop)
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
endif

inquire(unit=14,pos=lr)

if (kode.eq.0) then
  nsc=nsc+1
!  write(*,*) runtime,ntop,(typ(i),i=1,ntop)
else
  close(14)
  inpopen=.false.
endif
end subroutine

subroutine writebtrhead(outfile,unitw,nframew,itypew,atnamw)
implicit none
integer*4 j,kode,unitw
integer*4 nframew,itypew
character*4 atnamw(itypew)
character outfile*(*)

open(unit=unitw,file=trim(outfile),form='unformatted',iostat=kode,status='replace',action='write')
write(unitw) nframew
write(unitw) itypew
write(unitw) (atnamw(j),j=1,itypew)  ! ion and nucleotides types in char
end subroutine

subroutine writebtrbody(unitw,runtimew,ntopw,typw,rtw)
implicit none
real*8 runtimew
integer*4 i,ntopw,unitw
integer*4 typw(ntopw)
real*4 rtw(3,ntopw)

write(unitw) runtimew
write(unitw) ntopw
write(unitw) (typw(i),i=1,ntopw)
write(unitw) (rtw(1,i),i=1,ntopw)
write(unitw) (rtw(2,i),i=1,ntopw)
write(unitw) (rtw(3,i),i=1,ntopw)
!nscw=nscw+1
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

prname='BTREXTRACT'
prver='version 1.0'
prdesc='Extract from BROMOC Trajectory (.btr) a particular frame to another BROMOC trajectory file (.btr)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='28 Feb 2014'
lastdate='28 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

