!    RMSDBTR - Computes root mean square deviation of DNA with respect to the first structre from BROMOC Trajectory (.btr)
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
integer*4 nsc,ntop,itype,nframe,maxntop,ir,lr,ndna
real*8 runtime
logical*1 inpopen
character*4, allocatable :: atnam2(:)
character bs*8
real*4,allocatable :: rt(:,:),r0(:,:)
real*8 :: xtlabc6(6)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title

end module

program rmsdbtr
use comun
implicit none
real*8 rmsd,crmsd,rmsdf
integer narg,arg
character inpfile*256,outfile*256,line*2048

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)
call readarg('Input output filename: ',narg,arg,outfile)

! Read gcmcbd input file
call readgcmcbdhead(inpfile,.false.)
call readgcmcbdbody()
call checkndna()

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

open(unit=1,file=trim(outfile))
do while (inpopen)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  call compgyrad(rmsdf,rmsd,crmsd)           ! compute rmsd
  write(line,*) nsc,rmsdf,rmsd,crmsd                  ! write output
  write(1,'(A)') trim(line)
  if (inpopen) call readgcmcbdbody()          ! read next frame
enddo
close(1)

write(*,'(/A)') 'Normal termination'
end program

subroutine checkndna()
use comun
implicit none
integer i
character*4 chr

ndna=0
do i=1,ntop
  chr=atnam2(typ(i))
  if (chr=='Ab'.or.chr=='Gb'.or.chr=='Cb'.or.chr=='Tb'.or.chr=='S'.or.chr=='P') ndna=ndna+1
enddo
allocate (r0(3,ndna))
r0(1:3,1:ndna)=rt(1:3,1:ndna)
if (ndna==0) stop 'No DNA Particles'
end subroutine

subroutine compgyrad(rmsdf,rmsd,crmsd)
use comun
implicit none
real*8 rmsd,rx,ry,rz,rx2,ry2,rz2,x,y,z,x2,y2,z2,rgx,rgy,rgz,crmsd,rmsdf
integer i

rx=0.0
ry=0.0
rz=0.0
rx2=0.0
ry2=0.0
rz2=0.0
rmsdf=0.0
do i=1,ndna
  rx=rx+rt(1,i)
  ry=ry+rt(2,i)
  rz=rz+rt(3,i)
  rx2=rx2+rt(1,i)**2
  ry2=ry2+rt(2,i)**2
  rz2=rz2+rt(3,i)**2
  rmsdf=rmsdf+(rt(1,i)-r0(1,i))**2+(rt(2,i)-r0(2,i))**2+(rt(3,i)-r0(3,i))**2
enddo

x=(rx/ndna)**2
y=(ry/ndna)**2
z=(rz/ndna)**2
x2=rx2/ndna
y2=ry2/ndna
z2=rz2/ndna
rgx=x2-x
rgy=y2-y
rgz=z2-z
crmsd=sqrt(rgx+rgy+rgz)
rmsd=sqrt(x2+y2+z2)
rmsdf=sqrt(rmsdf/ndna)
end subroutine

subroutine readgcmcbdhead(inpfile,no)
use comun
implicit none
integer j, kode
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

prname='RMSDBTR'
prver='version 1.0'
prdesc='Computes root mean square deviation of DNA with respect to the first structre from BROMOC Trajectory (.btr)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='16 Jan 2014'
lastdate='16 Jan 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

