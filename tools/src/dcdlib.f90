!    DCDLIB - DCD Handler Library
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

module dcdlib
implicit none

contains

subroutine xtl6toxtl12(xtlabc6,xtlabc12)
implicit none
real*8 xtlabc6(6),xtlabc12(12)

xtlabc6(1)=xtlabc12(1)
xtlabc6(2)=xtlabc12(2)
xtlabc6(3)=xtlabc12(5)
xtlabc6(4)=xtlabc12(3)
xtlabc6(5)=xtlabc12(6)
xtlabc6(6)=xtlabc12(9)

end subroutine

subroutine xtl12toxtl6(xtlabc12,xtlabc6)
implicit none
real*8 xtlabc6(6),xtlabc12(12)

xtlabc12(10:12)=0d0 ! box
xtlabc12(1)=xtlabc6(1)
xtlabc12(2)=xtlabc6(2)
xtlabc12(3)=xtlabc6(4)
xtlabc12(4)=xtlabc6(2)
xtlabc12(5)=xtlabc6(3)
xtlabc12(6)=xtlabc6(5)
xtlabc12(7)=xtlabc6(4)
xtlabc12(8)=xtlabc6(5)
xtlabc12(9)=xtlabc6(6)
end subroutine

subroutine readdcdebody(u,nsc,xtlabc12,na,x,charmm,dcdopen)
implicit none
integer kode,i,u,na,nsc
real*8 xtlabc12(12)
real*4 x(3,na)
logical*1 charmm,dcdopen

if (.not.charmm) read(u,iostat=kode) xtlabc12
read(u,iostat=kode) (x(1,i),i=1,na)
read(u,iostat=kode) (x(2,i),i=1,na)
read(u,iostat=kode) (x(3,i),i=1,na)
if (kode.eq.0) then
  nsc=nsc+1
else
  close(u)
  dcdopen=.false.
endif
end subroutine

subroutine readdcdhead(dcdfile,u,ntitle,tnf,hdr,icntrl,itemp,title,natom,charmm,dcdopen)
implicit none
integer*4 icntrl(20),itemp,kode,u,tnf
integer*4 nfile,npriv,nsavc,nstep,nfree
integer*4 natom,ntitle
character dcdfile*(*),hdr*4
character*1,allocatable :: title(:)
logical*1 charmm,dcdopen

open(unit=u,file=dcdfile,form='unformatted')
dcdopen=.true.
read(u) hdr,icntrl
ntitle=icntrl(20)/12*80
if (allocated(title)) deallocate (title)
allocate (title(ntitle))
title=''
read(u,iostat=kode) itemp,title
read(u) natom
nfile=icntrl(1)
npriv=icntrl(2)
nsavc=icntrl(3)
nstep=icntrl(4)
!if(icntrl(9).gt.0) print *, '# fixed atoms = ',icntrl(9)
nfree = natom-icntrl(9)
!print *, '# of free atoms = ',nfree
!print *, 'total # atom = ', natom,nstep,nsavc
charmm=.false.
if (icntrl(2).eq.0) charmm=.true.
if (nstep.le.0) nstep=1
if (nsavc.le.0) nsavc=1
tnf = nstep/nsavc

write(*,'(A,I0)') 'Total number of frames: ',tnf
end subroutine

subroutine readdcdbody(u,nsc,xtlabc6,na,x,charmm,dcdopen)
implicit none
integer kode,i,u,na,nsc
real*8 xtlabc6(6)
real*4 x(3,na)
logical*1 charmm,dcdopen

if (.not.charmm) read(u,iostat=kode) xtlabc6
read(u,iostat=kode) (x(1,i),i=1,na)
read(u,iostat=kode) (x(2,i),i=1,na)
read(u,iostat=kode) (x(3,i),i=1,na)
if (kode.eq.0) then  
  nsc=nsc+1
else
  close(u)
  dcdopen=.false.
endif
end subroutine

subroutine writedcdhead(dcdfile,un,tnf,hdr,icntrl,itemp,natoms)
implicit none
integer*4 un,tnf
character dcdfile*(*)
integer*4 icntrl(20),itemp,ntitle,natoms
character hdr*4
character*1,allocatable :: title(:)
open(unit=un,file=trim(dcdfile),form='unformatted')
icntrl=0
hdr='CORD'
icntrl(10)=1026003171
icntrl(11)=1
icntrl(20)=24
itemp=2
ntitle=icntrl(20)/12*80
allocate (title(ntitle))
call assgn(ntitle,title,'REMARKS CREATED BY DCDLIB v1.0 REMARKS DATE: 2013/05/02 CREATED BY: Pablo M. De Biase')
icntrl(1)=tnf
icntrl(2)=1
icntrl(3)=1
icntrl(4)=tnf
write(un) hdr,icntrl
write(un) itemp,title
write(un) natoms
deallocate (title)
end subroutine

subroutine writedcdebody(un,xtlabc12,nn,x)
implicit none
integer i,nn,un
real*4 :: x(3,nn)
real*8 :: xtlabc12(12)

write(un) xtlabc12
write(un) (x(1,i),i=1,nn)
write(un) (x(2,i),i=1,nn)
write(un) (x(3,i),i=1,nn)
end subroutine

subroutine writedcdbody(un,xtlabc6,nn,x)
implicit none
integer i,nn,un
real*4 :: x(3,nn)
real*8 :: xtlabc6(6)

write(un) xtlabc6
write(un) (x(1,i),i=1,nn)
write(un) (x(2,i),i=1,nn)
write(un) (x(3,i),i=1,nn)
end subroutine

subroutine assgn(m,var,text)
implicit none
integer i,n,m
character*1 var(m)
character text*(*)
n=len_trim(text)
if (n.gt.m) stop 'error in string length'
do i=1,n
  var(i)=text(i:i)
enddo
end subroutine

end module
