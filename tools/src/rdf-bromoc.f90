!    RDF-BROMOC - Reads BROMOC trajectory (.btr) and computes RDF for spherical system
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
character*4, allocatable :: atnam(:)
character bs*8
real*4,allocatable :: rt(:,:)
real*8,allocatable :: g(:,:)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title
real*8 rad1,delr,idelr,voltot
integer maxbin,numpairs
real*8, parameter :: pi=3.14159265358979323846264338327950288419716939937510d0
real*8,parameter :: pi43=4d0/3d0*pi
end module

program rdfbromoc
use comun
implicit none
integer h,i,j,k,l,m,a,b,c,d,narg,arg
character inpfile*256,line*256
integer fac

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

call readarg('Spherical Radius: ',narg,arg,line)
read (line,*) rad1

voltot=pi43*rad1**3

call readarg('Resolution factor in Angstroms [0.1]: ',narg,arg,line)
if (len_trim(line).eq.0) then
  delr=0.1d0
else
  read(line,*) delr
endif
idelr=1d0/delr

write(*,*)

! Read gcmcbd input file
call readgcmcbdhead(inpfile,.false.)
call readgcmcbdbody()

numpairs=itype*(itype+1)/2
maxbin=2*rad1*idelr+1
allocate (g(maxbin,numpairs))

g=0.0d0

! Print number of frames read
write(*,'(/A,I8$)') 'Frame: ',nsc

do while (inpopen)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  call rdf()                                  ! calculate rdf for frame i
  call readgcmcbdbody()                       ! read next fram
enddo
close(3) ! close writedcd

g=g*(voltot/nsc)  ! Normalize
call saverdf()    ! write out results

write(*,'(/A)') 'Normal termination of DNACDF'
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

if (.not.allocated(atnam)) allocate (atnam(itype))
if (.not.allocated(maxntyp)) allocate (maxntyp(itype))

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

function indexi(itype,jtype)
implicit none
integer itype,jtype,indexi

indexi=MAX(itype,jtype)*(MAX(itype,jtype)-1)/2+MIN(itype,jtype)
end function

subroutine rdf()
use comun
implicit none
integer i,j,bin,numb(numpairs),indexi,typei
real*8 dv(3),rd,vola,volb,volc,vold,dist,vol,gt(maxbin,numpairs)

numb=0
do i=1,ntop
enddo

gt=0.0d0

do i=1,ntop-1
  do j=i+1,ntop
    dv=rt(1:3,j)-rt(1:3,i)
    rd=sqrt(dot_product(dv,dv))
    bin=rd*idelr+1
    typei=indexi(typ(i),typ(j))
    numb(typei)=numb(typei)+1
    dist=sqrt(dot_product(rt(1:3,i),rt(1:3,i)))
    call volume(rad1,bin*delr,dist,vola)
    call volume(rad1,(bin-1)*delr,dist,volb)
    dist=sqrt(dot_product(rt(1:3,j),rt(1:3,j)))
    call volume(rad1,bin*delr,dist,volc)
    call volume(rad1,(bin-1)*delr,dist,vold)
    if ((vola-volb).ne.0d0.and.(volc-vold).ne.0d0) then
      gt(bin,typei)=gt(bin,typei)+0.5d0/(vola-volb)+0.5d0/(volc-vold)
    else
      write(*,*) 'Warning: Division by zero'
    endif
  enddo
enddo
do i=1,itype
  do j=i,itype
    typei=indexi(i,j)
    g(:,typei)=g(:,typei)+gt(:,typei)/numb(typei)
  enddo
enddo
end subroutine

subroutine volume(r1,r2,d,v) ! volumen of a sphere or volume of intersection of two spheres
implicit none
real*8 r1,r2,d,v,rmin,rmax
real*8, parameter :: pi=3.14159265358979323846264338327950288419716939937510d0
real*8,parameter :: pi43=4d0/3d0*pi

rmin=min(r1,r2)
rmax=max(r1,r2)
if (d+rmin.gt.rmax.and.d.ne.0d0) then
  v=pi*(r1+r2-d)**2*(d**2+2d0*d*r2-3d0*r2**2+2d0*d*r1+6d0*r1*r2-3d0*r1**2)/(12d0*d)
else
  v=pi43*rmin**3
endif
end subroutine

subroutine saverdf()
use comun
implicit none
integer i,j,k,typei,indexi

do i=1,itype
  do j=i,itype
    typei=indexi(i,j)
    open(unit=1,file='rdf-'//trim(atnam(i))//'-'//trim(atnam(j))//'.dat')
    do k=1,maxbin
      write(1,*) k*delr,g(k,typei)
    enddo
    close(unit=1)
  enddo
enddo
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

prname='RDF-BROMOC'
prver='version 1.0'
prdesc='Reads BROMOC trajectory (.btr) and computes RDF for spherical system'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 Jan 2012'
lastdate='11 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

