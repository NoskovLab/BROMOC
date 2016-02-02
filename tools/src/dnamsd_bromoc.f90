!    DNAMSD-BROMOC - Computes mean square displacement of DNA fragments
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
integer,allocatable :: gfi(:),gff(:),fg(:),fp(:),ifp(:),typ(:),gft(:,:)
integer gfn,fn,nsc,ion,ntop,itype,nframe,nn,maxntop,snd,fag,dna
real*8 runtime,trv(3)
logical inpopen
character*4, allocatable :: atnam2(:)
real*8,allocatable :: msd(:),msdr(:),msdc(:),msdz(:),msdzr(:),msdzc(:),rdna(:,:,:),rcm0(:,:),trans(:)
character bs*8
integer,allocatable :: grf(:,:),gaf(:,:),sndi(:),sndf(:),pmdi(:,:),pmd(:),pmdf(:,:),pmdg(:)
integer,allocatable :: fti(:),ftf(:),ftl(:)
character*4,allocatable :: fl(:),gfl(:)
real*8,allocatable :: rt(:,:)
logical termon
end module

program dnamsd
use comun
implicit none
real*8 aa,bb,mdl
integer h,i,j,k,l,m,a,b,c,d
character inpfile*256,line*256,mskfile*256
integer ul(256),ll(256),num,fac,narg,arg

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

call readarg('Input Fragment Mask (.msk) filename [HELP and END]: ',narg,arg,mskfile)

if (len_trim(mskfile).eq.0) then
  write(*,'(A)') 'Order is important and must be 1)BASE 2)SUGAR 3)PHOSPATE 4)IONS'
  write(*,'(A)') 'Each record is composed of 4 characters'
  write(*,'(A)') 'Per line: 1st record -> group, 2nd record -> fragments'
  write(*,'(A)') 'Example of a mask file:'
  write(*,'(A)') 'BASEAb  Tb  Cb  Gb  '
  write(*,'(A)') 'SUGAS   '
  write(*,'(A)') 'PHOSP   '
  write(*,'(A)') 'IONSPOT CLA '
  stop
endif

! Read monk input file
call readmonkhead(inpfile,.false.)

allocate (msd(nframe),msdr(nframe),msdc(nframe),msdz(nframe),msdzr(nframe),msdzc(nframe),trans(nframe))

!Read mask file
call readmask(mskfile)

! Build residue type list
call fragtype

! Build fragment list and allocate cent var
write(*,'(A)') 'Read fragments:'
  write(*,'(A)') 'Fragment Group Number | Fragment Group Label | Fragment number | Fragment Label' 
do i=1,fn
  write(*,'(I5,3x,A4,3x,I5,3x,A4)') fg(i),gfl(fg(i)),i,fl(i)
enddo
write(*,'(A)') 'NOTE: In mask file group fragments order is BASE SUGA PHOS IONS'

ion=4    ! the group fragment for ions

allocate (fti(gff(ion)-gfi(ion)+1),ftf(gff(ion)-gfi(ion)+1))

call readarg('Choose Fragment number to compute the autocorrelation function [3]: ',narg,arg,line)
if(len_trim(line).eq.0) then
  fac=3
else
  read(line,*) fac
endif

call readarg('Input translocation reference vector (x y z): ',narg,arg,line)
trv=0d0
read(line,*) trv(1),trv(2),trv(3)
mdl=1d0/dsqrt(dot_product(trv,trv))
trv=trv*mdl

!Read first frame
call readmonkbody()

call fraglist()                         ! Rebuild fragment elements list for dna

allocate (rdna(3,snd,fag),rcm0(3,snd))

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

do while (inpopen)                        ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc              ! print frame number
  call fragilist()                        ! Rebuild fragment elements list for ions
  call calcautocorr(fac)                  ! compute autocorrelation (absolute bead MSD and MSD of the chain center of mass)
  call readmonkbody()                     ! read next frame
enddo

! Print out result
call printout()

write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine calcautocorr(fac)
use comun
implicit none
integer i,j,k,fac,a,b,nfac(snd),nbead
real*8 rrr(3,snd,fag),rcmt(3),rit(3),rcm(3),rdif(3)

do i=1,snd
  a=pmdi(fac,i)
  b=pmdf(fac,i)
  nfac(i)=0
  do j=a,b
    nfac(i)=nfac(i)+1
    if (nsc.eq.1) rdna(1:3,i,nfac(i))=rt(1:3,gaf(fg(fac),pmd(j)))
    rrr(1:3,i,nfac(i))=rt(1:3,gaf(fg(fac),pmd(j)))
  enddo
enddo
if (nsc.eq.1) then                            !
  do i=1,snd                                  !
    rcm0(1:3,i)=0d0                           !
    do j=1,nfac(i)                            !
      rcm0(1:3,i)=rcm0(1:3,i)+rdna(1:3,i,j)   !  compute center of fac for t=0
    enddo                                     !
    rcm0(1:3,i)=rcm0(1:3,i)/nfac(i)           !
  enddo                                       !
endif                                         !
msd(nsc)=0d0   !
msdr(nsc)=0d0  !  initialize mean-square displacement, msd relative to center & msd of the center
msdc(nsc)=0d0  !
msdz(nsc)=0d0   !
msdzr(nsc)=0d0  !  initialize z mean-square displacement, z msd relative to center & z msd of the center
msdzc(nsc)=0d0  !
trans(nsc)=0d0

nbead=0
do i=1,snd
  rcmt=0d0
  do j=1,nfac(i)                         !
    rcmt(1:3)=rcmt(1:3)+rrr(1:3,i,j)     ! 
  enddo                                  !  compute center of fac for t
  rcmt(1:3)=rcmt(1:3)/nfac(i)            !
  rcm(1:3)=rcmt(1:3)-rcm0(1:3,i)         !
  trans(nsc)=trans(nsc)+dot_product(rcm,trv)
  msdc(nsc)=msdc(nsc)+dot_product(rcm,rcm)      
  msdzc(nsc)=msdzc(nsc)+rcm(3)**2      
  do j=1,nfac(i)
    rit(1:3)=rrr(1:3,i,j)-rdna(1:3,i,j)
    rdif=rit-rcm
    msd(nsc)=msd(nsc)+dot_product(rit,rit)
    msdr(nsc)=msdr(nsc)+dot_product(rdif,rdif)
    msdz(nsc)=msdz(nsc)+rit(3)**2
    msdzr(nsc)=msdzr(nsc)+rdif(3)**2
  enddo
  nbead=nbead+nfac(i)
enddo
msd(nsc)=msd(nsc)/nbead     ! mean square displacement for each bead
msdr(nsc)=msdr(nsc)/nbead   ! mean square displacement relative to the center
msdc(nsc)=msdc(nsc)/snd     ! mean square displacement of the center
trans(nsc)=trans(nsc)/snd   ! translation along trv
msdz(nsc)=msdz(nsc)/nbead   ! z-mean square displacement for each bead
msdzr(nsc)=msdzr(nsc)/nbead ! z-mean square displacement relative to the center
msdzc(nsc)=msdzc(nsc)/snd   ! z-mean square displacement of the center
end subroutine

! fn is the number of fragments
! gfn is the number of group fragments
subroutine readmask(mskfile)
use comun
integer i,j,k,s,kode
character line*256,mskfile*256

open (unit=1,file=mskfile,iostat=kode)
gfn=0
fn=0
do while (kode.eq.0)
  line=''
  read(1,'(A)',iostat=kode) line
  k=len_trim(line)-4 ! because the first 8 chars correspond to the fragment group label and fragment label
  if (k.gt.0) then
    gfn=gfn+1
    if (mod(k,4).ne.0) then
      k=int(k/4)+1
    else
      k=int(k/4)
    endif
    fn=fn+k
  endif
enddo
rewind(1)
allocate (gfi(gfn),gff(gfn),gfl(gfn),fl(fn),fg(fn),fp(fn),ifp(fn))

kode=0
j=0
i=0
do while (kode.eq.0)
  line=''
  read(1,'(A)',iostat=kode) line
  k=len_trim(line)-4   ! substracting the group fragment label and fragment label
  if (k.gt.0) then
    j=j+1
    if (mod(k,4).ne.0) then
      k=int(k/4)+1
    else
      k=int(k/4)
    endif
    gfl(j)=line(1:4)
    gfi(j)=i+1
    do s=2,k+1
      i=i+1
      fl(i)=line(s*4-3:s*4)
      fg(i)=j    ! gives the group number as a function of the fragment number
    enddo
    gff(j)=i
  endif
enddo
close(1)

if (itype.gt.fn) stop 'Mask list is incomplete'
if (itype.lt.fn) stop 'Too many fragments in mask list'
end subroutine

subroutine fragtype()
use comun
implicit none
integer i,j
logical loopon,notfound

fp(1:fn)=0

do i=1,itype
  j=0
  loopon=.true.
  do while (j.lt.fn.and.loopon)
    j=j+1
    if (fl(j).eq.atnam2(i)) then
      fp(j)=i    ! < fragment type number in mask list order > fragment type number in input order
      ifp(i)=j   ! < fragment type number in input order > fragment type number in mask list order
      loopon=.false.
    endif  
  enddo
  if (loopon) then 
     write(*,*) 'Error: Fragment type ',atnam2(i),' was not found in the mask list'
     stop
  endif
enddo
notfound=.false.
do i=1,fn
  if (fp(i).eq.0) then 
    write(*,*) 'Could not find fragment ',fl(i),' in the input file'
    notfound=.true.
  endif
enddo
if (notfound) stop
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

subroutine fraglist()
use comun
implicit none
integer h,i,j,k,a,b,c
logical loopon,yes

i=0
dna=0
loopon=.true.
do while (i.lt.ntop.and.loopon)
  i=i+1
  if (ifp(typ(i)).le.gff(3)) then
    dna=dna+1
  else
    loopon=.false.
  endif
enddo


if (allocated(grf)) deallocate (grf)
if (allocated(gaf)) deallocate (gaf)
allocate (grf(3,dna),gaf(3,dna))

i=1
fag=0
snd=0
yes=.true.
do while (i.lt.dna)
  a=fg(ifp(typ(i)))
  b=fg(ifp(typ(i+1)))
  c=fg(ifp(typ(i+2)))
  yes=yes.and.a.le.3.and.b.le.3.and.c.le.3
  yes=yes.and.a.ge.1.and.b.ge.1.and.c.ge.1
  yes=yes.and.a.ne.b.and.b.ne.c.and.c.ne.a
  if (yes) then 
    fag=fag+1
    grf(a,fag)=ifp(typ(i))         !
    grf(b,fag)=ifp(typ(i+1))       !  Contains the fragment number in mask order for each group and complete nucleotide number
    grf(c,fag)=ifp(typ(i+2))       !
    gaf(a,fag)=i                   !
    gaf(b,fag)=i+1                 !  Contains the position in input order for each nucleotide fragment (base,sugar, phosp) and complete nucleotide number
    gaf(c,fag)=i+2                 !
  else
    i=i-1
    snd=snd+1
    yes=.true.
  endif
  i=i+3
enddo
if (allocated(sndi)) deallocate (sndi)
if (allocated(sndf)) deallocate (sndf)
allocate (sndi(snd),sndf(snd))

i=1
j=0
k=0
yes=.true.
sndi(1)=1                 ! initial position of strand in fag list
do while (i.lt.dna)
  a=fg(ifp(typ(i)))
  b=fg(ifp(typ(i+1)))
  c=fg(ifp(typ(i+2)))
  yes=yes.and.a.le.3.and.b.le.3.and.c.le.3
  yes=yes.and.a.ge.1.and.b.ge.1.and.c.ge.1
  yes=yes.and.a.ne.b.and.b.ne.c.and.c.ne.a
  if (yes) then
    j=j+1
  else
    k=k+1
    sndf(k)=j            ! final position of strand in fag list
    if (k+1.le.snd) sndi(k+1)=j+1  ! initial position of strand in fag list
    i=i-1
    yes=.true.
  endif
  i=i+3
enddo

b=gff(3)
a=gfi(1)
if (allocated(pmdi)) deallocate (pmdi)
if (allocated(pmdf)) deallocate (pmdf)
if (allocated(pmd)) deallocate (pmd)
if (allocated(pmdg)) deallocate (pmdg)
allocate (pmdi(a:b,snd),pmdf(a:b,snd),pmd(dna),pmdg(fag))

do i=1,fag
  pmdg(i)=i
enddo

h=0
do i=1,snd
  do j=gfi(1),gff(3)
    pmdi(j,i)=h+1   ! initial position in pmd list for each base and strand
    do k=sndi(i),sndf(i)
      if (j.eq.grf(fg(j),k)) then
        h=h+1
        pmd(h)=k ! pmd list contains the the position in fag list 
      endif
    enddo
    pmdf(j,i)=h    ! final position in pmd list for each base and strand
  enddo
enddo
end subroutine

subroutine readmonkhead(inpfile,no)
use comun
implicit none
integer i, j, kode
character*256 inpfile
logical no

open(14,file=inpfile,form='unformatted',iostat=kode,position='rewind')
read(14,iostat=kode) nframe                 ! number of frames
read(14,iostat=kode) itype                  ! number of ions and nucleotides
if (.not.allocated(atnam2)) allocate (atnam2(itype))
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

subroutine printout()
use comun
implicit none
integer i,j,k,h,m
character filename*256

! Print out result
open(unit=1,file='msd.dat')
open(unit=2,file='msdr.dat')
open(unit=3,file='msdc.dat')
open(unit=4,file='msdz.dat')
open(unit=7,file='msdzr.dat')
open(unit=8,file='msdzc.dat')
open(unit=9,file='transloc.dat')
do i=1,nsc
  write(1,*) i,msd(i),dsqrt(msd(i))
  write(2,*) i,msdr(i),dsqrt(msdr(i))
  write(3,*) i,msdc(i),dsqrt(msdc(i))
  write(4,*) i,msdz(i),dsqrt(msdz(i))
  write(7,*) i,msdzr(i),dsqrt(msdzr(i))
  write(8,*) i,msdzc(i),dsqrt(msdzc(i))
  write(9,*) i,trans(i)
enddo
close(1)
close(2)
close(3)
close(4)
close(7)
close(8)
close(9)
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

prname='DNAMSD-BROMOC'
prver='version 2.6'
prdesc='Computes mean square displacement of DNA fragments'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 Jul 2011'
lastdate='11 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine


