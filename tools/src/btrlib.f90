!    BTRLIB - BTR Handler Library (BTR: Bromoc Trajectory File)
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

module btrlib
implicit none
integer*4 nsc,unitn
integer*4 ntop,itype,nframe,maxntop,ir,lr
integer*4,allocatable :: typ(:),maxntyp(:)
real*8 runtime
real*4,allocatable :: rt(:,:)
character*4, allocatable :: atnam(:)
logical*1 inpopen
!for writing
integer*4 unitw,nframew,itypew,ntopw,nscw
character*4, allocatable :: atnamw(:)
logical*1 outopen
real*8 runtimew
real*4,allocatable :: rtw(:,:)
integer*4,allocatable :: typw(:)

contains

subroutine checkmaxntop()
implicit none
integer*4 kode,i,cnttyp(itype),recordi,recordl

if (.not.allocated(maxntyp)) allocate (maxntyp(itype))
maxntyp=0
if (inpopen) kode=0
do while (kode.eq.0.and.nsc.lt.nframe)
  if (allocated(typ)) deallocate (typ)
  if (allocated(rt))  deallocate (rt)
  read(unitn,iostat=kode) recordi
  read(unitn,iostat=kode) runtime ! simulation time for each step (ns)
  read(unitn,iostat=kode) recordl
  call check(recordi,recordl,unitn,kode)
  read(unitn,iostat=kode) recordi
  read(unitn,iostat=kode) ntop ! number of particles
  read(unitn,iostat=kode) recordl
  call check(recordi,recordl,unitn,kode)
  if (kode.eq.0) then
    allocate (typ(ntop),rt(3,ntop))
    read(unitn,iostat=kode) recordi
    read(unitn,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
    read(unitn,iostat=kode) recordl
    call check(recordi,recordl,unitn,kode)
    cnttyp=0
    do i=1,ntop
      cnttyp(typ(i))=cnttyp(typ(i))+1
    enddo
    do i=1,itype
      if (cnttyp(i).gt.maxntyp(i)) maxntyp(i)=cnttyp(i)
    enddo
    read(unitn,iostat=kode) recordi
    read(unitn,iostat=kode) (rt(1,i),i=1,ntop)
    read(unitn,iostat=kode) recordl
    call check(recordi,recordl,unitn,kode)
    read(unitn,iostat=kode) recordi
    read(unitn,iostat=kode) (rt(2,i),i=1,ntop)
    read(unitn,iostat=kode) recordl
    call check(recordi,recordl,unitn,kode)
    read(unitn,iostat=kode) recordi
    read(unitn,iostat=kode) (rt(3,i),i=1,ntop)
    read(unitn,iostat=kode) recordl
    call check(recordi,recordl,unitn,kode)
    if (kode.eq.0) nsc=nsc+1
  endif
enddo
close(unitn)
maxntop=0
do i=1,itype
  maxntop=maxntop+maxntyp(i)
enddo
end subroutine

! read bromoc trj header
! if no is false prints info
subroutine readbtrhead(inpfile,un,no)
implicit none
integer*4 j,kode,un
integer*4 recordi,recordl
character inpfile*256,line*4096
logical no

unitn=un
open(unit=unitn,file=inpfile,form='unformatted',iostat=kode,position='rewind',access='stream')
read(unitn,iostat=kode) recordi
read(unitn,iostat=kode) nframe                 ! number of frames
read(unitn,iostat=kode) recordl
call check(recordi,recordl,unitn,kode)

read(unitn,iostat=kode) recordi
read(unitn,iostat=kode) itype                  ! number of ions and nucleotides
read(unitn,iostat=kode) recordl
call check(recordi,recordl,unitn,kode)

if (.not.allocated(atnam)) allocate (atnam(itype))

read(unitn,iostat=kode) recordi
read(unitn,iostat=kode) (atnam(j),j=1,itype)  ! ion and nucleotides types in char
read(unitn,iostat=kode) recordl
call check(recordi,recordl,unitn,kode)

if (.not.no) then
  write(*,'(/A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(line,*) 'Fragment types: ',(atnam(j),j=1,itype)
  write(*,'(A)') trim(adjustl(line))
endif

if (kode.eq.0) then
  inpopen=.true.
else
  close(unitn)
endif
end subroutine

subroutine writebtrhead(outfile,un)
implicit none
integer*4 j,kode,un
character outfile*256

unitw=un
open(unit=unitw,file=outfile,form='unformatted',iostat=kode,status='replace',action='write')
write(unitw) nframew 
write(unitw) itypew 
write(unitw) (atnamw(j),j=1,itypew)  ! ion and nucleotides types in char
outopen=.true.
end subroutine

subroutine closewritebtr()
implicit none
close(unitw)
outopen=.false.
end subroutine

subroutine closereadbtr()
implicit none
close(unitn)
inpopen=.false.
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

! read bromoc trj body
subroutine readbtrbody()
implicit none
integer*4 kode,i
integer*4 recordi,recordl

inquire(unit=unitn,pos=ir)

if (allocated(typ)) deallocate (typ)
if (allocated(rt)) deallocate (rt)
read(unitn,iostat=kode) recordi
read(unitn,iostat=kode) runtime !simulation time for each step (ns)
read(unitn,iostat=kode) recordl
call check(recordi,recordl,unitn,kode)

read(unitn,iostat=kode) recordi
read(unitn,iostat=kode) ntop ! number of particles
read(unitn,iostat=kode) recordl
call check(recordi,recordl,unitn,kode)

if (kode.eq.0) then 
  allocate (typ(ntop),rt(3,ntop))

  read(unitn,iostat=kode) recordi
  read(unitn,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(unitn,iostat=kode) recordl
  call check(recordi,recordl,unitn,kode)
  
  read(unitn,iostat=kode) recordi
  read(unitn,iostat=kode) (rt(1,i),i=1,ntop)
  read(unitn,iostat=kode) recordl
  call check(recordi,recordl,unitn,kode)
  
  read(unitn,iostat=kode) recordi
  read(unitn,iostat=kode) (rt(2,i),i=1,ntop)
  read(unitn,iostat=kode) recordl
  call check(recordi,recordl,unitn,kode)
  
  read(unitn,iostat=kode) recordi
  read(unitn,iostat=kode) (rt(3,i),i=1,ntop)
  read(unitn,iostat=kode) recordl
  call check(recordi,recordl,unitn,kode)
endif

inquire(unit=unitn,pos=lr)

if (kode.eq.0) then
  nsc=nsc+1
!  write(*,*) runtime,ntop,(typ(i),i=1,ntop)
else
  close(unitn)
  inpopen=.false.
endif
end subroutine

subroutine writebtrbody()
implicit none
integer*4 i
write(unitw) runtimew
write(unitw) ntopw
write(unitw) (typw(i),i=1,ntopw)
write(unitw) (rtw(1,i),i=1,ntopw)
write(unitw) (rtw(2,i),i=1,ntopw)
write(unitw) (rtw(3,i),i=1,ntopw)
nscw=nscw+1
end subroutine
end module
