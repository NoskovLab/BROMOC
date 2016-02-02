!    BTRCONV - Converts BROMOC trajectory to PDB,GRO,XTC,DCD
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
integer*4 nsc,ntop,itype,nframe,natoms,ndna,dtype
real*8 runtime
logical*1 inpopen,xtc,dcd
character*4, allocatable :: atnam2(:)
character bs*8
! for btr, dcd, xtc
real*4,allocatable :: rt(:,:)
! for xtc
real*4 :: box(9)
! for dcd
real*8 :: xtlabc6(6)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title
end module

program btrconv
use comun
implicit none
real*4 time,prec
integer*4 i,a,narg,rnsc,xd,ret,arg,iprec
character inpfile*256,line*256,num*64
integer*4 fac
logical*1 noions

bs=repeat(achar(8),len(bs))
call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()
! printout header
write(*,'(/A/)') 'USAGE: BTRCONV filename.btr savingfreq saveions (dcd/xtc)'

! read trajectory filename
call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

! Save frequency
call readarg('Saving Frequency: ',narg,arg,line)
read(line,*) fac

! Save ions?
noions=.true.
call readarg('Save ions (y/n)? [n] ',narg,arg,line)
if (line(1:1).eq.'Y'.or.line(1:1).eq.'y') noions=.false.

dcd=.true.
xtc=.true.
call readarg('Format (dcd/xtc) [both]:   ',narg,arg,line)
if (line(1:3).eq.'dcd'.or.line(1:3).eq.'DCD') then
  xtc=.false.
  dcd=.true.
elseif (line(1:3).eq.'xtc'.or.line(1:3).eq.'XTC') then
  xtc=.true.
  dcd=.false.
endif
    
! Read gcmcbd input file
call readbtrhead(inpfile,.false.)

! count particles and number of frames
call checknatoms()
rnsc=nsc
write(num,'(i0)') rnsc
call readarg('Final Frame ['//trim(num)//'] : ',narg,arg,line)
if (len_trim(line).gt.0) read(line,*) rnsc
if (allocated(rt)) deallocate (rt)
allocate (rt(3,natoms))

call readbtrhead(inpfile,.true.)
call readbtrbody()

! count number of dna particles
ndna=0
do i=1,itype
  if (atnam2(i).eq.'S'.or.atnam2(i).eq.'P'.or.atnam2(i).eq.'Ab'.or.     &
      atnam2(i).eq.'Cb'.or.atnam2(i).eq.'Gb'.or.atnam2(i).eq.'Tb') then
    ndna=ndna+maxntyp(i)
    dtype=i
  else
    exit
  endif
enddo
if (noions) natoms=ndna

a=index(inpfile,'.',back=.true.)
a=a-1
if (a.le.-1) then
  if (len_trim(inpfile).gt.0) then
    a=len_trim(inpfile)
  else
    stop 'Filename error'
  endif
endif

iprec=6
if (xtc) then 
  call readarg('Precision [6]:  ',narg,arg,line)
  if (len_trim(line).gt.0) read(line,*) iprec 
endif
prec=1.0*10**iprec

box=(/ 1e2, (0e0, i=1,3), 1e2, (0e0, i=1,3), 1e2 /)

if (dcd) then
  line=inpfile(1:a)//'.pdb'
  open(unit=2,file=trim(line))  ! open writepdb
  call writeout(2,1) ! write pdb
  close(2) ! close writegro
  ! Define dcd header
  icntrl=0
  hdr='CORD'
  icntrl(1)=int(rnsc/fac)
  icntrl(2)=fac
  icntrl(3)=fac
  icntrl(4)=rnsc
  icntrl(10)=1026003171
  icntrl(11)=1
  icntrl(20)=24
  itemp=2
  title='REMARKS FILENAME=bromoc.dcd CREATED BY BTRCONV v1.0         REMARKS DATE: 2012/10/29 CREATED BY: Pablo M. De Biase'
  xtlabc6=0d0
  ! create output filename
  line=inpfile(1:a)//'.dcd'
  ! write dcd header
  call writedcdhead(line,3)
endif

if (xtc) then 
  line=inpfile(1:a)//'.gro'
  open(unit=2,file=trim(line))  ! open writegro
  call writeout(2,2) ! write gro
  close(2) ! close writegro
  line=inpfile(1:a)//'.xtc'
  ! open xtc
  call xdrfopen(xd,trim(line),"w",ret)
endif

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

do while (inpopen.and.nsc.le.rnsc)          ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  if (mod(nsc,fac).eq.0) then
    time=runtime
    if (xtc) call writextc(xd,nsc,time,prec)
    if (dcd) call writedcdbody(3)
  endif
  if (nsc.le.rnsc) call readbtrbody()                       ! read next frame
enddo

if (dcd) close(3) ! close dcd
if (xtc) call xdrfclose(xd,ret)
write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine writedcdhead(dcdfile,un)
use comun
implicit none
integer*4 un
character dcdfile*256
open(unit=un,file=trim(dcdfile),form='unformatted')
write(un) hdr,icntrl
write(un) itemp,title(1:160)
write(un) natoms
end subroutine

subroutine writedcdbody(un)
use comun
implicit none
integer*4 un,i,j,k,l,m,rar(ndna+1:natoms)
real*4 x(3,natoms)

x(1:3,1:ndna)=rt(1:3,1:ndna)

if (natoms.gt.ndna) then
  k=ndna
  do i=dtype+1,itype
    l=0
    do j=1,ntop
      if (typ(j).eq.i) then
        l=l+1
        k=k+1
        rar(k)=j
      endif
    enddo
    m=rar(k)
    do j=l+1,maxntyp(i)
      k=k+1
      rar(k)=m
    enddo
  enddo
  if (k.ne.natoms) stop 'Unexpected error 3'

  do i=ndna+1,natoms
    x(1,i)=rt(1,rar(i))
    x(2,i)=rt(2,rar(i))
    x(3,i)=rt(3,rar(i))
  enddo
endif

write(un) xtlabc6
write(un) (x(1,i),i=1,natoms)
write(un) (x(2,i),i=1,natoms)
write(un) (x(3,i),i=1,natoms)
end subroutine

subroutine writextchead(natoms,step,time,xd)
implicit none
integer*4 magic,natoms,step,xd,ret
real*4 time
magic=1995
call xdrfint(xd,magic,ret)
call xdrfint(xd,natoms,ret)
call xdrfint(xd,step,ret)
call xdrffloat(xd,time,ret)

!open(unit=un,file=trim(xtcfile),form='unformatted',access='stream')
!write(un) magic
!write(un) natoms
!write(un) step
!write(un) time
end subroutine

subroutine writextc(xd,step,time,prec)
use comun
implicit none
integer*4 i,j,k,l,m,rar(ndna+1:natoms),xd,ret,step
real*4 x(3*natoms),prec,time

do i=1,ndna
  x(3*i-2)=1d-1*rt(1,i)
  x(3*i-1)=1d-1*rt(2,i)
  x(3*i)=1d-1*rt(3,i)
enddo
if (natoms.gt.ndna) then  
  k=ndna
  do i=dtype+1,itype
    l=0
    do j=1,ntop
      if (typ(j).eq.i) then
        l=l+1
        k=k+1
        rar(k)=j
      endif
    enddo
    m=rar(k)
    do j=l+1,maxntyp(i)
      k=k+1
      rar(k)=m
    enddo
  enddo
  if (k.ne.natoms) stop 'Unexpected error 3'
  
  do i=ndna+1,natoms
    x(3*i-2)=1d-1*rt(1,rar(i))
    x(3*i-1)=1d-1*rt(2,rar(i))
    x(3*i)=1d-1*rt(3,rar(i))
  enddo
endif

call writextchead(natoms,step,time,xd)
do i=1,9
   call xdrffloat(xd,box(i),ret)
end do
call xdrf3dfcoord(xd,x,natoms,prec,ret)

end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='BTRCONV'
prver='version 1.0'
prdesc='Converts BROMOC trajectory to PDB,GRO,XTC,DCD'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='29 Oct 2012'
lastdate='29 Oct 2012'

write(*,'(/A)') trim(prname)//' '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A/)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
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

subroutine readbtrhead(inpfile,no)
use comun
implicit none
integer j, kode
character*256 inpfile
logical no

open(14,file=inpfile,form='unformatted',iostat=kode,position='rewind')
read(14,iostat=kode) nframe                 ! number of frames
read(14,iostat=kode) itype                  ! number of ions and nucleotides
if (.not.allocated(atnam2)) allocate (atnam2(itype))
if (.not.allocated(maxntyp)) allocate (maxntyp(itype))
read(14,iostat=kode) (atnam2(j),j=1,itype)  ! ion and nucleotides types in char
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
  stop 'Cannot read filename'
endif
end subroutine

subroutine readbtrbody()
use comun
implicit none
integer kode,i

if (allocated(typ)) deallocate (typ)
read(14,iostat=kode) runtime !simulation time for each step (ns)
read(14,iostat=kode) ntop ! number of particles
if (kode.eq.0) then 
  allocate (typ(ntop))
  read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(14,iostat=kode) (rt(1,i),i=1,ntop)
  read(14,iostat=kode) (rt(2,i),i=1,ntop)
  read(14,iostat=kode) (rt(3,i),i=1,ntop)
endif
if (kode.eq.0) then
  nsc=nsc+1
else
  close(14)
  inpopen=.false.
endif
end subroutine

subroutine checknatoms()
use comun
implicit none
integer*4 kode,i,cnttyp(itype)

maxntyp=0
if (inpopen) kode=0
do while (kode.eq.0.and.nsc.lt.nframe)
  if (allocated(typ)) deallocate (typ)
  if (allocated(rt))  deallocate (rt)
  read(14,iostat=kode) runtime ! simulation time for each step (ns)
  read(14,iostat=kode) ntop ! number of particles
  if (kode.eq.0) then
    allocate (typ(ntop),rt(3,ntop))
    read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
    cnttyp=0
    do i=1,ntop
      cnttyp(typ(i))=cnttyp(typ(i))+1
    enddo
    do i=1,itype
      if (cnttyp(i).gt.maxntyp(i)) maxntyp(i)=cnttyp(i)
    enddo
    read(14,iostat=kode) (rt(1,i),i=1,ntop)
    read(14,iostat=kode) (rt(2,i),i=1,ntop)
    read(14,iostat=kode) (rt(3,i),i=1,ntop)
    if (kode.eq.0) nsc=nsc+1
  endif
enddo
close(14)
natoms=0
do i=1,itype
  natoms=natoms+maxntyp(i)
enddo
end subroutine

subroutine writeout(unitn,oformat)
use comun
implicit none
integer i,j,k,l,m,unitn,oformat,nawo,rar(ndna+1:natoms)
integer,allocatable :: rnwo(:)
real*4,allocatable :: rwo(:,:)
character,allocatable :: rtwo(:)*5,atwo(:)*5

nawo=natoms
allocate (rnwo(natoms),rwo(3,natoms),rtwo(natoms),atwo(natoms))

do j=1,ndna
  rwo(1:3,j)=rt(1:3,j)
  rtwo(j)=atnam2(typ(j))
  atwo(j)=atnam2(typ(j))
  rnwo(j)=typ(j)
enddo
if (natoms.gt.ndna) then
  k=ndna
  do i=dtype+1,itype
    do j=1,maxntyp(i)
      k=k+1
      rnwo(k)=i
      rtwo(k)=atnam2(i)
      atwo(k)=atnam2(i)
    enddo
  enddo
  if (k.ne.natoms) stop 'Unexpected error 1'

  k=ndna
  do i=dtype+1,itype
    l=0
    do j=1,ntop
      if (typ(j).eq.i) then
        l=l+1
        k=k+1
        rar(k)=j
      endif
    enddo
    m=rar(k)
    do j=l+1,maxntyp(i)
      k=k+1
      rar(k)=m
    enddo
  enddo
  if (k.ne.natoms) stop 'Unexpected error 2'

  do j=ndna+1,natoms
    rwo(1:3,j)=rt(1:3,rar(j))
  enddo
endif

if (oformat.eq.1) then 
  call writepdb(unitn,nawo,rwo,atwo,rtwo,rnwo)
else
  call writegro(unitn,nawo,rwo,atwo,rtwo,rnwo,box)
endif
deallocate (rnwo,rwo,rtwo,atwo)
end subroutine

subroutine writegro(unitn,na,r,at,rt,rn,box)
implicit none
integer i,na,rn(na),unitn
real*4 r(3,na),box(9)
character rt(na)*5,at(na)*5,title*80

title='GRO created by BTR2XTC  author: Pablo M. De Biase, October 2012'
write(unitn,'(a80)') title
write(unitn,'(i5)') na
do i=1,na
  write(unitn,'(i5,2a5,i5,3f8.3)') rn(i),rt(i),at(i),i,r(1,i)*1d-1,r(2,i)*1d-1,r(3,i)*1d-1
enddo
write (unitn,'(9f10.5)') box(1),box(5),box(9),box(2:4),box(6:8)
end subroutine

subroutine writepdb(unitn,na,r,at,rt,rn)
implicit none
integer i,na,rn(na),unitn
real*4 r(3,na)
character rt(na)*5,at(na)*5

do i=1,na
  write (unitn,'(A6,I5,x,A5,A5,I4,4x,3F8.3)') 'ATOM  ',i,at(i),rt(i),rn(i),r(1,i),r(2,i),r(3,i)
enddo
write (unitn,'(A)') 'END'
end subroutine

