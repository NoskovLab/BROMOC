!    GYRADBTR - Computes radius of gyration (gyradius) of DNA BROMOC Trajectory (.btr)
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
integer*4 nsc,ntop,itype,nframe,maxntop,ir,lr,ndna,ndnaf,ndnal
real*8 runtime
logical*1 inpopen
character*4, allocatable :: atnam2(:)
character bs*8
real*4,allocatable :: rt(:,:)
real*8 :: xtlabc6(6),m(6)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title

end module

program gyradbtr
use comun
implicit none
real*8 x,y,z,x2,y2,z2,rgx,rgy,rgz,rgxy,rgc,xm,ym,zm,xm2,ym2,zm2,rgmx,rgmy,rgmz,rgmxy,rgmc
integer narg,arg,inif,finf
character inpfile*256,outfile*256,line*1024,str*10

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
write(*,*) 'Number of DNA particles: ',ndna
write(str,'(I0)') ndna
call readarg('Select Range of DNA particles to consider (eg: 10 40) [1 '//trim(str)//'] : ',narg,arg,line)
if (len_trim(line).eq.0) then
  ndnaf=1
  ndnal=ndna
else
  read(line,*) ndnaf, ndnal
  if (ndnaf.lt.1) ndnaf=1
  if (ndnal.lt.ndnaf) ndnal=ndna
endif

write(*,*) 'Total Frames: ',nframe
write(str,'(I0)') nframe
call readarg('Select Frame Range (eg: 10 100) [1 '//trim(str)//'] : ',narg,arg,line)
if (len_trim(line).eq.0) then
  inif=1
  finf=nframe
else
  read(line,*) inif, finf
  if (inif.lt.1) inif=1
  if (finf.lt.inif) finf=nframe
endif

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

open(unit=1,file=trim(outfile)//'raw.dat')
open(unit=2,file=trim(outfile)//'raw-m.dat')
open(unit=3,file=trim(outfile)//'.dat')
open(unit=4,file=trim(outfile)//'-m.dat')
open(unit=10,file=trim(outfile)//'-xy.dat')
open(unit=11,file=trim(outfile)//'-mxy.dat')
do while (inpopen.and.nsc.le.finf)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  if (nsc.ge.inif) then
    call compgyrad(x,y,z,x2,y2,z2,rgx,rgy,rgz,rgxy,rgc,xm,ym,zm,xm2,ym2,zm2,rgmx,rgmy,rgmz,rgmxy,rgmc)           ! compute gyradii
    write(1,*) nsc,x,y,z,x2,y2,z2                  ! write output
    write(2,*) nsc,xm,ym,zm,xm2,ym2,zm2            ! write output
    write(3,*) nsc,sqrt(rgx),sqrt(rgy),sqrt(rgz),sqrt(rgc)            ! write output
    write(4,*) nsc,sqrt(rgmx),sqrt(rgmy),sqrt(rgmz),sqrt(rgmc)       ! write output
    write(10,*) nsc,sqrt(rgxy),atan2(sqrt(rgy),sqrt(rgx))            ! write output
    write(11,*) nsc,sqrt(rgmxy),atan2(sqrt(rgmy),sqrt(rgmx))       ! write output
  endif
  if (inpopen) call readgcmcbdbody()          ! read next frame
enddo
close(1)
close(2)
close(3)
close(4)
close(10)
close(11)

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
if (ndna==0) stop 'No DNA Particles'
m(1)=94.9714822
m(2)=83.108959
m(3)=134.119288
m(4)=110.094498
m(5)=150.118718
m(6)=125.105935
end subroutine

function elemn(nn)
use comun
implicit none
character*4 chr
integer elemn,n,nn
chr=atnam2(nn)
n=0
if (chr=='P') then
  n=1
else if (chr=='S') then
  n=2
else if (chr=='Ab') then
  n=3
else if (chr=='Cb') then
  n=4
else if (chr=='Gb') then
  n=5
else if (chr=='Tb') then
  n=6
endif
if (n==0) stop 'Not DNA'
elemn=n
end function

subroutine compgyrad(x,y,z,x2,y2,z2,rgx,rgy,rgz,rgxy,rgc,xm,ym,zm,xm2,ym2,zm2,rgmx,rgmy,rgmz,rgmxy,rgmc)
use comun
implicit none
real*8 rgx,rgy,rgz,rgxy,rgc,rgmx,rgmy,rgmz,rgmxy,rgmc,rx2,ry2,rz2,rx,ry,rz,rmx,rmy,rmz,rmx2,rmy2,rmz2,x,y,z,x2,y2,z2,xm,ym,zm,xm2,ym2,zm2,mm
integer i,elemn,alldna

rx=0.0
ry=0.0
rz=0.0
rmx=0.0
rmy=0.0
rmz=0.0
rmx2=0.0
rmy2=0.0
rmz2=0.0
rx2=0.0
ry2=0.0
rz2=0.0
mm=0.0
alldna=ndnal-ndnaf+1

do i=ndnaf,ndnal
  rmx=rmx+rt(1,i)*m(elemn(typ(i)))
  rmy=rmy+rt(2,i)*m(elemn(typ(i)))
  rmz=rmz+rt(3,i)*m(elemn(typ(i)))
  rmx2=rmx2+m(elemn(typ(i)))*rt(1,i)**2
  rmy2=rmy2+m(elemn(typ(i)))*rt(2,i)**2
  rmz2=rmz2+m(elemn(typ(i)))*rt(3,i)**2
  rx=rx+rt(1,i)
  ry=ry+rt(2,i)
  rz=rz+rt(3,i)
  rx2=rx2+rt(1,i)**2
  ry2=ry2+rt(2,i)**2
  rz2=rz2+rt(3,i)**2
  mm=mm+m(elemn(typ(i)))
enddo
xm=(rmx/mm)**2
ym=(rmy/mm)**2
zm=(rmz/mm)**2
xm2=rmx2/mm
ym2=rmy2/mm
zm2=rmz2/mm
rgmx=xm2-xm
rgmy=ym2-ym
rgmz=zm2-zm
rgmxy=rgmx+rgmy
rgmc=rgmx+rgmy+rgmz

x=(rx/alldna)**2
y=(ry/alldna)**2
z=(rz/alldna)**2
x2=rx2/alldna
y2=ry2/alldna
z2=rz2/alldna
rgx=x2-x
rgy=y2-y
rgz=z2-z
rgxy=rgx+rgy
rgc=rgx+rgy+rgz
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

prname='GYRADBTR'
prver='version 1.0'
prdesc='Computes radius of gyration (gyradius) of DNA BROMOC Trajectory (.btr)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='16 Jan 2014'
lastdate='16 Jan 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

