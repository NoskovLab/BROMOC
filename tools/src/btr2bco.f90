!    BTR2BCO - Extract from BROMOC Trajectory (.btr) a BROMOC Coordinate file (.bco)
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

program btr2bco
use comun
implicit none
real*8 aa,bb
integer h,i,j,k,l,m,a,b,c,d,narg,arg
character inpfile*256,line*256,outfile*256
integer fac
logical*1 loopon

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)
call readarg('Input BROMOC coordinate (.bco) filename: ',narg,arg,outfile)

! Read gcmcbd input file
call readgcmcbdhead(inpfile,.false.)
call readgcmcbdbody()

if (nframe.gt.1) then 
  line=''
  do while (len_trim(line).eq.0)
    call readarg('Choose frame to save: ',narg,arg,line)
  enddo
  read(line,*) fac
else
  fac=1
endif

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

loopon=.true.
do while (inpopen.and.nsc.le.fac.and.loopon)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  if (nsc.eq.fac) then
    open(unit=1,file=outfile)
    call printframe(1)           ! print coordinates
    close(1) 
    loopon=.false.
  endif
  if (inpopen) call readgcmcbdbody()          ! read next fram
enddo

write(*,'(/A)') 'Normal termination'
end program

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

subroutine printframe(iunit)
use comun
implicit none
integer i,iunit,nsites,str,nions,j,typei,tp,mxtp
character ch*4,nn*1
write(*,'(//)') 
!write(*,'(A,I0,A,I0)') 'First and Last Byte: ',ir, '   ',lr
write(*,*) 'Simulation time for this frame (ns): ',runtime !simulation time for each step (ns)
write(*,*) 'Number of particles: ',ntop ! number of particles
!write(*,*) 'Type         Type#     x              y                z'

mxtp=0
do i=1,ntop
  ch=atnam2(typ(i))
  mxtp=max(typ(i),mxtp)
  if (ch.eq.'Ab'.or.ch.eq.'Cb'.or.ch.eq.'Gb'.or.ch.eq.'Tb'.or.ch.eq.'S'.or.ch.eq.'P') then
    nsites=i
  else
    exit
  endif
enddo

if (nsites.gt.0) then
  str=1
  write(iunit,'(i5)') nsites
  do i = 1, nsites
    ch=atnam2(typ(i))
    if (ch.eq.'Ab') then
      nn='A'
    elseif (ch.eq.'Tb') then
      nn='T'
    elseif (ch.eq.'Gb') then
      nn='G'
    elseif (ch.eq.'Cb') then
      nn='C'
    endif
    write(iunit,'(i1,1x,a1,1x,a2,3(1x,f15.8))') str,nn,adjustl(ch),rt(1,i),rt(2,i),rt(3,i)
  enddo
endif
tp=int(nsites/3)
nions = ntop - nsites
write(iunit,'(i5)') nions
if (nions.gt.0) then
  do i = nsites+1, ntop
    typei=typ(i)-mxtp+tp
    if (rt(3,i).lt.0) typei=-typei
    j=0
    write(iunit,'(2(i5,1x),3(f15.8,1x),i4)') i,typei,rt(1,i),rt(2,i),rt(3,i),j
  enddo
endif

inpopen=.false.
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

prname='BTR2BCO'
prver='version 1.1'
prdesc='Extract from BROMOC Trajectory (.btr) a BROMOC Coordinate file (.bco)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='08 Jan 2013'
lastdate='10 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

