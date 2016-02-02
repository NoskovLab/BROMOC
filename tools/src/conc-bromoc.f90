!    CONC-BROMOC - Reads BROMOC Trajectory and computes concentration for each particle type
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
integer*4,allocatable :: typ(:),ntyp(:,:)
integer*4 nsc,ntop,itype,nframe,maxntop,ir,lr
real*8 runtime
logical*1 inpopen
character*4, allocatable :: atnam(:)
character bs*8
real*4,allocatable :: rt(:,:)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title
real*8 rad1,delr,idelr,voltot,box(3)
real*8, parameter :: pi=3.14159265358979323846264338327950288419716939937510d0
real*8,parameter :: pi43=4d0/3d0*pi
end module

program concbromoc
use comun
implicit none
integer i,narg,arg
character inpfile*256,line*256
logical*1 boxt

do i=1,len(bs)
  bs(i:i)=achar(8)
enddo

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

call readarg('Spherical or Box System (sph/box): ',narg,arg,line)
boxt=.false.
if (line(1:1).eq.'b'.or.line(1:1).eq.'B') boxt=.true.

if (boxt) then 
  call readarg('Box X: ',narg,arg,line)
  read(line,*) box(1)

  call readarg('Box Y: ',narg,arg,line)
  read(line,*) box(2)

  call readarg('Box Z: ',narg,arg,line)
  read(line,*) box(3)
  voltot=box(1)*box(2)*box(3)
else
  call readarg('Spherical Radius: ',narg,arg,line)
  read (line,*) rad1
  voltot=pi43*rad1**3
endif

write(*,*) 'Total Volume: ',voltot

! Read gcmcbd input file
call readgcmcbdhead(inpfile,.false.)
call readgcmcbdbody()

allocate (ntyp(itype,nframe))
ntyp=0

! Print number of frames read
write(*,'(/A,I8$)') 'Frame: ',nsc
do while (inpopen)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  call countpart()
  call readgcmcbdbody()                       ! read next fram
enddo
close(3) ! close writedcd

call saveconc(inpfile(1:index(inpfile,'.',back=.true.)-1))

write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine countpart()
use comun
implicit none
integer i

do i=1,ntop
  ntyp(typ(i),nsc)=ntyp(typ(i),nsc)+1
enddo

end subroutine

subroutine saveconc(pname)
use comun
implicit none
character*(*) pname
integer i,j
real*8 ivoltot

ivoltot=1d0/(voltot*0.00060221412927d0)
do i=1,itype
  open(unit=1,file='conc-'//trim(atnam(i))//'-'//trim(adjustl(pname))//'.dat')
  do j=1,nsc
    write(1,*) j,ntyp(i,j)*ivoltot,ntyp(i,j)
  enddo
  close(unit=1)
enddo
end subroutine

subroutine readgcmcbdhead(inpfile,no)
use comun
implicit none
integer j,kode
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

if (.not.allocated(atnam)) allocate (atnam(itype))

read(14,iostat=kode) recordi
read(14,iostat=kode) (atnam(j),j=1,itype)  ! ion and nucleotides types in char
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (.not.no) then
  write(*,'(A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(*,*) achar(8)//'Fragment types: ',(atnam(j),j=1,itype)
endif
nsc=0
if (kode.eq.0) then
  inpopen=.true.
else
  close(14)
endif
end subroutine

subroutine check(recordi,recordl,unitn,kode)
implicit none
integer*4 recordi,recordl,pos,unitn,kode

if (kode.eq.0) then
  if (recordi.ne.recordl) then
    inquire(unit=unitn,pos=pos)
    write (*,*)
    write (*,*) 'Record not matching: ',recordi,recordl,' at byte: ',pos
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

prname='CONC-BROMOC'
prver='version 1.0'
prdesc='Reads BROMOC Trajectory and computes concentration for each particle type'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 May 2012'
lastdate='07 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

