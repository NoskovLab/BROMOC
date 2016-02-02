!    FLUX - Computes the ions flux using BROMOC trajectory
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
integer*4,allocatable :: gfi(:),gff(:),fg(:),fp(:),ifp(:),typ(:),gft(:,:)
integer*4 gfn,fn,nsc,ion,ntop,itype,nframe,nn,maxntop,snd,fag,dna,sfreq,pn,ni
real*8 runtime
real*8 dt,uli,lli,cr,ori(3)
logical inpopen
character*4, allocatable :: atnam2(:)
character bs*8
integer*4,allocatable :: fti(:),ftf(:),ftl(:),pi(:),pni(:),pnf(:)
character*4,allocatable :: fl(:),gfl(:)
real*8,allocatable :: rt(:,:),flux(:,:),q(:),pz(:)
logical termon
end module

program flux_bromoc
use comun
implicit none
real*8 aa,bb
integer*4 h,i,j,k,l,m,a,b,c,d,nii,narg,arg
character inpfile*256,line*256
integer*4 ul(256),ll(256),num,fac

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

a=index(inpfile,'.',back=.true.)

! Read monk input file
call readmonkhead(inpfile,.true.)

! Build residue type list
call readarg('Input Sequence of Species Numbers: ',narg,arg,line)
call findparm(line,256,num,ll,ul)
fn=num
gfn=1
allocate (gfi(gfn),gff(gfn),gfl(gfn),fl(fn),fg(fn),fp(fn),ifp(itype))
do i=1,fn
  read(line(ll(i):ul(i)),*) fp(i)
  ifp(fp(i))=i  
  fl(i)=atnam2(fp(i))
enddo
gfl(1)='IONS'
gfi(1)=1
gff(1)=fn
fg(1:fn)=1

ion=1    ! the group fragment for ions
ni=gff(ion)-gfi(ion)+1

allocate (fti(ni),ftf(ni),flux(nframe,ni),q(ni),pni(ni),pnf(ni))

nii=0
do i=gfi(ion),gff(ion)
  nii=nii+1
  write(*,'(A,A,A$)') 'Set Charge for particle ',fl(i),': '
  call readarg('',narg,arg,line)
  read(line,*) q(nii)
enddo

call readarg('Saving Frequency, Time Step: ',narg,arg,line)
if(len_trim(line).eq.0) then
  stop 'integer real real'
else
  read(line,*) sfreq,dt
endif

call readarg('Lower Limit, Upper Limit, Cylinder Radius, Origin (x,y): ',narg,arg,line)
read(line,*) lli,uli,cr,ori(1),ori(2)
cr=cr**2

!Read first frame
call readmonkbody()

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

do while (inpopen)                        ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc              ! print frame number
  call fragilist()                        ! build fragment elements list for ions
  call calcflux()                         ! compute autocorrelation (absolute bead MSD and MSD of the chain center of mass)
  call readmonkbody()                     ! read next frame
enddo
inpfile=inpfile(1:a)
inpfile(a+1:)='flux'
! Print out result
call printout(inpfile)

write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine calcflux()
use comun
implicit none
integer i,j,k,l,m,c,d,a,pit(ntop),b
real*8 pzt(ntop),v,dv
logical go
m=0
b=0
do j=gfi(ion),gff(ion) ! for each ion type
  m=m+1
  flux(nsc,m)=0d0
  c=fti(m) ! first atom of the selected ion type
  d=ftf(m) ! last atom of the selected ion type
  if (nsc.gt.1) then
    do i=pni(m),pnf(m)
      a=c-1
      do while (a.lt.d)
        a=a+1
        if (pi(i).eq.ftl(a)) then
          flux(nsc,m)=flux(nsc,m)+q(m)*(rt(3,ftl(a))-pz(i))
          a=d   
        endif
      enddo
    enddo
  endif
  pni(m)=b+1
  do i=c,d
    a=ftl(i)
    dv=(rt(1,a)-ori(1))**2+(rt(2,a)-ori(2))**2
    v=rt(3,a)
    go=(v.le.uli.and.v.ge.lli.and.dv.le.cr)
    if (go) then
      b=b+1
      pit(b)=a
      pzt(b)=rt(3,a)
    endif
  enddo
  pnf(m)=b
enddo
if (allocated(pz))deallocate (pz)
if (allocated(pi))deallocate (pi)
allocate (pz(ntop),pi(ntop))
pi(1:b)=pit(1:b)
pz(1:b)=pzt(1:b)

end subroutine

subroutine fragilist()
use comun
implicit none
integer i,j,k,m

if (allocated(ftl)) deallocate (ftl)
allocate (ftl(ntop))

k=0
m=0
do i=gfi(ion),gff(ion)
  m=m+1
  fti(m)=k+1       ! initial number in ftl list for fragment type i
  do j=1,ntop
    if (fp(i).eq.typ(j)) then
      k=k+1
      ftl(k)=j     ! list that contain the position for each fragment type sorted by fragment type
    endif
  enddo
  ftf(m)=k         ! final number in ftl list for fragment type i
enddo
end subroutine

subroutine readmonkhead(inpfile,yes)
use comun
implicit none
integer i, j, kode
character*256 inpfile
character*2048 line
logical yes

open(14,file=inpfile,form='unformatted',iostat=kode,position='rewind')
read(14,iostat=kode) nframe                 ! number of frames
read(14,iostat=kode) itype                  ! number of ions and nucleotides
if (.not.allocated(atnam2)) allocate (atnam2(itype))
read(14,iostat=kode) (atnam2(j),j=1,itype)  ! ion and nucleotides types in char
if (yes) then
  write(*,'(A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(line,*)'Fragment types: ',(j,' '//atnam2(j),j=1,itype)
  write(*,'(A)') trim(line)
endif
nsc=0
if (kode.eq.0) then
  inpopen=.true.
else
  close(14)
endif
end subroutine

subroutine readmonkbody()
use comun
implicit none
integer kode,i,k
real*4, allocatable :: rtt(:,:)

if (allocated(typ)) deallocate (typ)
if (allocated(rt)) deallocate (rt)
if (allocated(rtt)) deallocate (rtt)
read(14,iostat=kode) runtime !simulation time for each step (ns)
read(14,iostat=kode) ntop ! number of particles
k=kode
if (k.eq.0) then 
  allocate (rtt(3,ntop),rt(3,ntop),typ(ntop))
  read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(14,iostat=kode) (rtt(1,i),i=1,ntop)
  read(14,iostat=kode) (rtt(2,i),i=1,ntop)
  read(14,iostat=kode) (rtt(3,i),i=1,ntop)
endif
if (kode.eq.0) then
  nsc=nsc+1
!  write(*,*) runtime,ntop,(typ(i),i=1,ntop)
  rt(1:3,1:ntop)=dble(rtt(1:3,1:ntop))
else
  close(14)
  inpopen=.false.
endif
if (k.eq.0) deallocate (rtt)
end subroutine

subroutine checkmaxntop()
use comun
implicit none
integer kode,i
real*4, allocatable :: rtt(:,:)

maxntop=0
if (inpopen) kode=0
do while (kode.eq.0)
  if (allocated(typ)) deallocate (typ)
  if (allocated(rtt)) deallocate (rtt)
  read(14,iostat=kode) runtime !simulation time for each step (ns)
  read(14,iostat=kode) ntop ! number of particles
  if (kode.eq.0) then
    if (ntop.gt.maxntop) maxntop=ntop
    allocate (rtt(3,ntop),typ(ntop))
    read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
    read(14,iostat=kode) (rtt(1,i),i=1,ntop)
    read(14,iostat=kode) (rtt(2,i),i=1,ntop)
    read(14,iostat=kode) (rtt(3,i),i=1,ntop)
  endif
  if (kode.ne.0) close(14)
enddo
end subroutine

subroutine printout(filename)
use comun
implicit none
integer i,j,k,h,m
real*8 cons,cum(ni+1),flx(ni+1)
character filename*256,line*2048

cons=1.60217656535d+5/(dabs(uli-lli)*dt*dfloat(sfreq))
! Print out result
open(unit=1,file=trim(filename))

write(line,*) 'frame   time  ',(fl(i),i=gfi(ion),gff(ion)),'  total |-> Cumulative current'
  write(1,'(A)') trim(line)
cum=0d0
do i=2,nsc
  flx(ni+1)=0d0
  do j=1,ni
    flx(j)=flux(i,j)*cons
    flx(ni+1)=flx(ni+1)+flux(i,j)
    cum(j)=cum(j)+flx(j)
  enddo
  flx(ni+1)=flx(ni+1)*cons
  cum(ni+1)=cum(ni+1)+flx(ni+1)
  
  write(line,*) i,(dfloat(i)-1.5d0)*dfloat(sfreq)*dt,(flx(j),j=1,ni+1),(cum(j)/dfloat(i-1),j=1,ni+1)
  write(1,'(A)') trim(line)
enddo
close(1)
end subroutine

subroutine findparm(str,dimn,num,llim,ulim)
implicit none
integer dimn,i,length,num
integer ulim(dimn),llim(dimn)
character ( len = dimn ) str
logical chng
length=len_trim(str)
chng=.false.
num=0
ulim(1:dimn)=0
llim(1:dimn)=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
    endif
    chng=.true.
  endif
enddo
if (ulim(num).eq.0) ulim(num)=length
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

prname='FLUX'
prver='version 2.3'
prdesc='Computes the ions flux using BROMOC trajectory (.btr)'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 Aug 2011'
lastdate='11 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

