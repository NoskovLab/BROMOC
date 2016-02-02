!    MARKOVMACROTRJ - Produces Markovian Macrostates Trajectories
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
integer*4 nsc,ntop,itype,nframe,natoms,dtype,rnsc,tnf,ndna
real*8 runtime
logical*1 inpopen,charmm,inputbtr
character*4, allocatable :: atnam(:)
character bs*8
! for btr, dcd, xtc
real*4,allocatable :: rt(:,:)
! for xtc
real*4 :: box(9)
! for dcd/dcde
real*8 :: xtlabc6(6),xtlabc12(12)
integer*4 icntrl(20),itemp,ntitle
character hdr*4
character*1,allocatable :: title(:)
! for clust-id, micro, macro
integer*4 ncid,nmic,nmac
integer*4,allocatable :: cid(:),mic(:),mac(:),mic2mac(:),cid2mic(:),cid2mac(:),mic2cid(:),fmic(:),fmac(:),fcid(:)
end module

program markovmacrotrj
use comun
implicit none
real*4 prec,preci,time
integer*4 i,narg,xd,ret,arg,iprec,ipreci
integer*4,allocatable :: xda(:)
character*256 filename(5),line
integer*4 step
integer*1 trjin,trjout
logical*1 rdcde,wdcde

bs=repeat(achar(8),len(bs))
call header()
arg=0
inputbtr=.false.
narg=COMMAND_ARGUMENT_COUNT()
! printout header

! read trajectory filename
call readarg('Input DCD/DCDE/XTC/BTR trajectory (.dcde/.dcd/.xtc/.btr) filename: ',narg,arg,filename(1))
call checkext(filename(1))
if (ext(filename(1)).eq.'dcd') then 
  trjin=1
  rdcde=.false.
elseif (ext(filename(1)).eq.'dcde') then 
  trjin=2
  rdcde=.true.
elseif (ext(filename(1)).eq.'xtc') then 
  trjin=3
elseif (ext(filename(1)).eq.'btr') then 
  trjin=4
else
  write(*,'(//A//)') 'Extension not recognized'
  stop
endif

call readarg('Input cluster-ID for each frame output filename: ',narg,arg,filename(2))
call readarg('Input MSMBuilder Microstates ID output filename: ',narg,arg,filename(3))
call readarg('Input MSMBuilder Macrostates ID output filename: ',narg,arg,filename(4))
call readarg('Input Output Trajectory (.dcde/.dcd/.xtc/.btr) for Macrostates filename: ',narg,arg,filename(5))

call checkext(filename(5))
if (ext(filename(5)).eq.'dcd') then
  trjout=1
  wdcde=.false.
elseif (ext(filename(5)).eq.'dcde') then   
  trjout=2
  wdcde=.true.
elseif (ext(filename(5)).eq.'xtc') then
  trjout=3
elseif (ext(filename(5)).eq.'btr') then    
  trjout=4
  if (trjin.lt.4) stop 'Error: Writing in BROMOC format from other not BROMOC trajectory not allowed'
else
  write(*,'(//A//)') 'Extension not recognized'
  stop
endif

! read trj input files
if (trjin.eq.1) then
  call readdcdhead(filename(1),1)
  call readdcdbody(1,rdcde)
elseif (trjin.eq.2) then
  call readdcdhead(filename(1),1)
  call readdcdbody(1,rdcde)
elseif (trjin.eq.3) then
  call xdrfopen(xd,trim(filename(1)),"r",ret)
  nsc=0
  call readxtc(xd,step,time,preci)
  ipreci=int(log10(preci))
  write(*,'(/A,I0)') 'Reading precision at: ',ipreci
elseif (trjin.eq.4) then
  call readbtrhead(filename(1),1,.false.)
  call checknatoms(1)
  ! count particles and number of frames
  tnf=nsc
  if (allocated(rt)) deallocate (rt)
  allocate (rt(3,natoms))
  call readbtrhead(filename(1),1,.true.)
  call readbtrbody(1)
endif

! read clustering information
call readclustid(filename(2))
call readmicromsm(filename(3))
call readmacromsm(filename(4))

if (trjout.eq.3) allocate (xda(nmac))
! save pdb and gro
if (trjin.eq.4.and.trjout.lt.4) then
  ! count number of dna particles
  ndna=0
  do i=1,itype
    if (atnam(i).eq.'S'.or.atnam(i).eq.'P'.or.atnam(i).eq.'Ab'.or.     &
        atnam(i).eq.'Cb'.or.atnam(i).eq.'Gb'.or.atnam(i).eq.'Tb') then
      ndna=ndna+maxntyp(i)
      dtype=i
    else
      exit
    endif
  enddo
  call readarg('Save ions (y/n)? [n] ',narg,arg,line)
  if (line(1:1).eq.'Y'.or.line(1:1).eq.'y') natoms=ndna
  inputbtr=.true.
  line=fname(filename(5))//'.pdb'
  open(unit=2,file=trim(line))
  call writeout(2,1)
  close(2)
  line=fname(filename(5))//'.gro'
  open(unit=2,file=trim(line))  ! open writegro
  call writeout(2,2) ! write gro
  close(2) ! close writegro
endif

iprec=ipreci
if (trjout.eq.3) then 
  call readarg('Precision ['//num2str(iprec,1)//']:  ',narg,arg,line)
  if (len_trim(line).gt.0) read(line,*) iprec 
endif
prec=1.0*10**iprec

call convertdata(trjin,trjout,time)

! open files
if (trjout.le.2) then
  if (trjin.ge.3) then
    ! Define dcd header
    icntrl=0
    hdr='CORD'
    icntrl(10)=1026003171
    icntrl(11)=1
    icntrl(20)=24
    itemp=2
    allocate (title(ntitle))
    title='REMARKS FILENAME=unknown.dcd CREATED BY MARKOVMACROTRJ v1.0 REMARKS DATE: 2012/12/10 CREATED BY: Pablo M. De Biase'
  endif
  do i=1,nmac
    icntrl(1)=fmac(i)
    icntrl(3)=1
    icntrl(4)=fmac(i)
    if (trjout.eq.1) line=fname(filename(5))//num2str(i,3)//'.dcd'
    if (trjout.eq.2) line=fname(filename(5))//num2str(i,3)//'.dcde'
    call writedcdhead(trim(line),9+i)
  enddo
elseif (trjout.eq.3) then
  do i=1,nmac
    line=fname(filename(5))//num2str(i,3)//'.xtc'
    call xdrfopen(xda(i),trim(line),"w",ret)
  enddo
elseif (trjout.eq.4) then
  do i=1,nmac
    line=fname(filename(5))//num2str(i,3)//'.btr'
    call writebtrhead(trim(line),9+i,fmac(i),itype,atnam)
  enddo
endif
 
! Print number of frames read
write(*,'(/A,I8$)') 'Frame: ',nsc

do while (nsc.le.rnsc)          ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  call convertdata(trjin,trjout,time)
  if (mac(nsc).gt.0) then
    if (trjout.le.2) call writedcdbody(mac(nsc)+9,wdcde)
    if (trjout.eq.3) call writextc(xda(mac(nsc)),nsc,time,prec)
    if (trjout.eq.4) call writebtrbody(mac(nsc)+9,runtime,ntop,typ,rt)
  endif
  if (nsc.lt.rnsc) then
    if (trjin.le.2) call readdcdbody(1,rdcde)
    if (trjin.eq.3) call readxtc(xd,step,time,preci)
    if (trjin.eq.4) call readbtrbody(1)
  else
    nsc=nsc+1
  endif
enddo

if (trjin.eq.3) then
  call xdrfclose(xd,ret)
else
  close(1)
endif

do i=1,nmac
  if (trjout.eq.3) then 
    call xdrfclose(xda(i),ret)
  else
    close(9+i)
  endif
enddo
write(*,'(//A)') 'Normal termination of DNACDF'
contains
  function ext(filename)
  implicit none
  character*(*) filename
  character*(len_trim(filename)-index(filename,'.',back=.true.)) ext
  ext=filename(index(filename,'.',back=.true.)+1:len_trim(filename))
  end function
  function fname(filename)
  implicit none
  character*(*) filename
  character*(index(filename,'.',back=.true.)-1) fname
  fname=filename(1:index(filename,'.',back=.true.)-1)
  end function
  subroutine checkext(filename)
  implicit none
  character*(*) filename
  if (len_trim(filename).le.0) stop 'Filename not detected'
  if (index(filename,'.',back=.true.).le.0) stop 'Filename must have extension'
  end subroutine
  function num2str(n,mx)
  implicit none
  integer n,num,i,a,mx
  character num2str*(mx),numero*10
  numero='0123456789'
  num=n
  do i=mx,1,-1
    a=int(num/(10**(i-1)))
    num2str(mx-i+1:mx-i+1)=numero(a+1:a+1)
    num=num-a*10**(i-1)
  enddo
  end function
end program

subroutine convertdata(trjin,trjout,time)
use comun
implicit none
integer*1 trjin,trjout
integer*4 i
real*4 time

if (trjout.eq.1) then
  if (trjin.eq.2) then
    ! box
    xtlabc6(1)=xtlabc12(1)
    xtlabc6(2)=xtlabc12(2)
    xtlabc6(3)=xtlabc12(5)
    xtlabc6(4)=xtlabc12(3)
    xtlabc6(5)=xtlabc12(6)
    xtlabc6(6)=xtlabc12(9)
  elseif (trjin.eq.3) then
    !box
    xtlabc6(1)=1d1*box(1)
    xtlabc6(2)=1d1*box(2)
    xtlabc6(3)=1d1*box(5)
    xtlabc6(4)=1d1*box(3)
    xtlabc6(5)=1d1*box(6)
    xtlabc6(6)=1d1*box(9)
  elseif (trjin.eq.4) then
    ! box
    xtlabc6(1)=1d3
    xtlabc6(2)=0d0
    xtlabc6(3)=1d3
    xtlabc6(4)=0d0
    xtlabc6(5)=0d0
    xtlabc6(6)=1d3
  endif
elseif (trjout.eq.2) then
  if (trjin.eq.1) then
    ! box
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
  elseif (trjin.eq.3) then
    xtlabc12(10:12)=0d0 ! box
    xtlabc12(1:9)=1d1*box(1:9) ! box
  elseif (trjin.eq.4) then
    xtlabc12=(/ 1d3, (0d0,i=1,3), 1d3, (0d0,i=1,3), 1d3, (0d0,i=1,3) /)  ! box
  endif
elseif (trjout.eq.3) then
  if (trjin.eq.1) then
    !box
    box(1)=1e-1*sngl(xtlabc6(1))
    box(2)=1e-1*sngl(xtlabc6(2))
    box(3)=1e-1*sngl(xtlabc6(4))
    box(4)=1e-1*sngl(xtlabc6(2))
    box(5)=1e-1*sngl(xtlabc6(3))
    box(6)=1e-1*sngl(xtlabc6(5))
    box(7)=1e-1*sngl(xtlabc6(4))
    box(8)=1e-1*sngl(xtlabc6(5))
    box(9)=1e-1*sngl(xtlabc6(6))
  elseif (trjin.eq.2) then
    box=1e-1*sngl(xtlabc12(1:9))  ! box
  elseif (trjin.eq.4) then
    time=sngl(runtime)   ! time
    box=(/ 1e2, (0e0, i=1,3), 1e2, (0e0, i=1,3), 1e2 /)   ! box
  endif
endif

end subroutine

subroutine readclustid(filename)
use comun
implicit none
integer*4 kode,n,t
character filename*(*),ln*1024,fst*1
open(unit=123,file=trim(filename),iostat=kode)
read(123,'(A)',iostat=kode) ln
n=0
do while(kode.eq.0)
  fst=adjustl(ln)
  if(fst.ne.'@'.and.fst.ne.'#'.and.fst.ne.'!') n=n+1
  read(123,'(A)',iostat=kode) ln
enddo
close(123)
rnsc=n
allocate (cid(n),mic(n),mac(n))
cid=0
mic=0
mac=0
open(unit=123,file=trim(filename),iostat=kode)
read(123,'(A)',iostat=kode) ln
n=0
do while(kode.eq.0)
  fst=adjustl(ln)
  if(fst.ne.'@'.and.fst.ne.'#'.and.fst.ne.'!') then
    n=n+1
!    read(ln,*) t,cid(n)
    read(ln,*) cid(n)
    cid(n)=cid(n)+1
  endif
  read(123,'(A)',iostat=kode) ln
enddo
close(123)

end subroutine

subroutine readmicromsm(filename)
use comun
implicit none
integer*4 kode,i
character filename*(*),ln*1024
open(unit=123,file=trim(filename))
ncid=0
read(123,'(A)',iostat=kode) ln
do while(kode.eq.0.and.len_trim(ln).gt.0)
  ncid=ncid+1
  read(123,'(A)',iostat=kode) ln
enddo
close(123)
allocate (cid2mic(ncid),fcid(ncid),cid2mac(ncid))
cid2mic=0
fcid=0
cid2mac=0
open(unit=123,file=trim(filename))
nmic=0
do i=1,ncid
  read(123,*) cid2mic(i)
  cid2mic(i)=cid2mic(i)+1
  if (cid2mic(i).gt.0) nmic=nmic+1
enddo
close(123)
allocate (mic2cid(0:nmic),mic2mac(0:nmic),fmic(0:nmic))
mic2cid=0
mic2mac=0
fmic=0
do i=1,ncid
  mic2cid(cid2mic(i))=i
enddo
do i=1,rnsc
  mic(i)=cid2mic(cid(i))
enddo
return
end subroutine

subroutine readmacromsm(filename)
use comun
implicit none
integer*4 i
integer*4,allocatable :: cmac(:),mmac(:)
character*(*) filename

open(unit=123,file=trim(filename))
nmac=0
do i=1,nmic
  read(123,*) mic2mac(i)
  mic2mac(i)=mic2mac(i)+1
  if (mic2mac(i).gt.nmac) nmac=mic2mac(i)
enddo
close(123)
allocate(fmac(0:nmac),cmac(0:nmac),mmac(0:nmac))
fmac=0
cmac=0
mmac=0
do i=1,rnsc
  mac(i)=mic2mac(mic(i))
enddo
fcid=0
fmic=0
fmac=0
do i=1,rnsc
  fcid(cid(i))=fcid(cid(i))+1
  fmic(mic(i))=fmic(mic(i))+1
  fmac(mac(i))=fmac(mac(i))+1
enddo
cmac=0
mmac=0
do i=1,ncid
  cid2mac(i)=mic2mac(cid2mic(i))
enddo
do i=1,ncid
  cmac(cid2mac(i))=cmac(cid2mac(i))+1
enddo
do i=0,nmic
  mmac(mic2mac(i))=mmac(mic2mac(i))+1
enddo

write(*,'(/A)') '     Macrostate |       Frames |     Clusters |  Microstates'
do i=0,nmac
  write(*,'(4I15)') i,fmac(i),cmac(i),mmac(i)
enddo
write(*,'(A5,I10,3I15)') 'TOTAL',nmac,rnsc,ncid,nmic

deallocate (cmac,mmac)
end subroutine


subroutine writedcdhead(dcdfile,un)
use comun
implicit none
integer*4 un
character dcdfile*(*)
open(unit=un,file=trim(dcdfile),form='unformatted')
write(un) hdr,icntrl
write(un) itemp,title
write(un) natoms
end subroutine

subroutine writedcdbody(un,dcde)
use comun
implicit none
integer*4 un,i,j,k,l,m,rar(ndna+1:natoms)
real*4 x(3,natoms)
logical*1 dcde

if (inputbtr) then
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
endif
if (.not.charmm) then 
  if (dcde) then
    write(un) xtlabc12
  else
    write(un) xtlabc6
  endif
endif
if (inputbtr) then
  write(un) (x(1,i),i=1,natoms)
  write(un) (x(2,i),i=1,natoms)
  write(un) (x(3,i),i=1,natoms)
else
  write(un) (rt(1,i),i=1,natoms)
  write(un) (rt(2,i),i=1,natoms)
  write(un) (rt(3,i),i=1,natoms)
endif
end subroutine

subroutine readxtchead(natoms,step,time,xd)
implicit none
integer*4 magic,natoms,step,xd,ret
real*4 time
magic=1995
call xdrfint(xd,magic,ret)
if (magic.ne.1995) stop 'wrong magic number in xtc'
call xdrfint(xd,natoms,ret)
call xtccheck(ret,1)
call xdrfint(xd,step,ret)
call xtccheck(ret,2)
call xdrffloat(xd,time,ret)
call xtccheck(ret,3)
end subroutine

subroutine readxtc(xd,step,time,prec)
use comun
implicit none
integer*4 xd,ret,step,i
real*4 prec,time
real*4,allocatable :: x(:)

call readxtchead(natoms,step,time,xd)

allocate (x(3*natoms))
if (.not.allocated(rt)) allocate (rt(3,natoms))

do i=1,9
   call xdrffloat(xd,box(i),ret)
end do

do i=1,natoms
  x(3*i-2)=1e-1*rt(1,i)
  x(3*i-1)=1e-1*rt(2,i)
  x(3*i)=1e-1*rt(3,i)
enddo

call xdrf3dfcoord(xd,x,natoms,prec,ret)

do i=1,natoms
  rt(1,i)=x(3*i-2)*1e1
  rt(2,i)=x(3*i-1)*1e1
  rt(3,i)=x(3*i)*1e1
enddo

deallocate (x)
nsc=nsc+1
end subroutine

subroutine xtccheck(ret,ivar)
implicit none
! Passed variables
integer ret,ivar
if ( ret .eq. 0 ) then
   if (ivar.eq.1) print *,'> XTC-Error reading/writing natoms'
   if (ivar.eq.2) print *,'> XTC-Error reading/writing step'
   if (ivar.eq.3) print *,'> XTC-Error reading/writing time'
   if (ivar.eq.4) print *,'> XTC-Error reading/writing box'
   if (ivar.eq.5) print *,'> XTC-Error reading/writing x'
stop
endif
return
end

subroutine writextchead(natoms,step,time,xd)
implicit none
integer*4 magic,natoms,step,xd,ret
real*4 time
magic=1995
call xdrfint(xd,magic,ret)
call xdrfint(xd,natoms,ret)
call xdrfint(xd,step,ret)
call xdrffloat(xd,time,ret)
end subroutine

subroutine writextc(xd,step,time,prec)
use comun
implicit none
integer*4 i,j,k,l,m,rar(ndna+1:natoms),xd,ret,step
real*4 x(3*natoms),prec,time

if (inputbtr) then
  do i=1,ndna
    x(3*i-2)=1e-1*rt(1,i)
    x(3*i-1)=1e-1*rt(2,i)
    x(3*i)=1e-1*rt(3,i)
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
      x(3*i-2)=1e-1*rt(1,rar(i))
      x(3*i-1)=1e-1*rt(2,rar(i))
      x(3*i)=1e-1*rt(3,rar(i))
    enddo
  endif
else
  do i=1,natoms
    x(3*i-2)=1e-1*rt(1,i)
    x(3*i-1)=1e-1*rt(2,i)
    x(3*i)=1e-1*rt(3,i)
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

prname='MARKOVMACROTRJ'
prver='version 1.0'
prdesc='Produces Markovian Macrostates Trajectories'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='07 Dec 2012'
lastdate='07 Dec 2012'

write(*,'(/A)') trim(prname)//' '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A/)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
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

subroutine readbtrhead(inpfile,un,no)
use comun
implicit none
integer j, kode,un
character*(*) inpfile
logical no

open(un,file=trim(inpfile),form='unformatted',iostat=kode,position='rewind')
read(un,iostat=kode) nframe                 ! number of frames
read(un,iostat=kode) itype                  ! number of ions and nucleotides
if (.not.allocated(atnam)) allocate (atnam(itype))
if (.not.allocated(maxntyp)) allocate (maxntyp(itype))
read(un,iostat=kode) (atnam(j),j=1,itype)  ! ion and nucleotides types in char
if (.not.no) then
  write(*,'(A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(*,*) 'Fragment types: ',(atnam(j),j=1,itype)
endif
nsc=0
if (kode.eq.0) then
  inpopen=.true.
else
  close(un)
  stop 'Cannot read filename'
endif
end subroutine

subroutine readbtrbody(un)
use comun
implicit none
integer kode,i,un

if (allocated(typ)) deallocate (typ)
read(un,iostat=kode) runtime !simulation time for each step (ns)
read(un,iostat=kode) ntop ! number of particles
if (kode.eq.0) then 
  allocate (typ(ntop))
  read(un,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(un,iostat=kode) (rt(1,i),i=1,ntop)
  read(un,iostat=kode) (rt(2,i),i=1,ntop)
  read(un,iostat=kode) (rt(3,i),i=1,ntop)
endif
if (kode.eq.0) then
  nsc=nsc+1
else
  close(un)
  inpopen=.false.
endif
end subroutine

subroutine checknatoms(un)
use comun
implicit none
integer*4 kode,i,cnttyp(itype),un

maxntyp=0
if (inpopen) kode=0
do while (kode.eq.0.and.nsc.lt.nframe)
  if (allocated(typ)) deallocate (typ)
  if (allocated(rt))  deallocate (rt)
  read(un,iostat=kode) runtime ! simulation time for each step (ns)
  read(un,iostat=kode) ntop ! number of particles
  if (kode.eq.0) then
    allocate (typ(ntop),rt(3,ntop))
    read(un,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
    cnttyp=0
    do i=1,ntop
      cnttyp(typ(i))=cnttyp(typ(i))+1
    enddo
    do i=1,itype
      if (cnttyp(i).gt.maxntyp(i)) maxntyp(i)=cnttyp(i)
    enddo
    read(un,iostat=kode) (rt(1,i),i=1,ntop)
    read(un,iostat=kode) (rt(2,i),i=1,ntop)
    read(un,iostat=kode) (rt(3,i),i=1,ntop)
    if (kode.eq.0) nsc=nsc+1
  endif
enddo
close(un)
natoms=0
do i=1,itype
  natoms=natoms+maxntyp(i)
enddo
end subroutine

subroutine writeout(unitn,oformat)
use comun
implicit none
integer*4 i,j,k,l,m,unitn,oformat,rar(ndna+1:natoms)
integer*4 :: rnwo(natoms)
real*4 :: rwo(3,natoms)
character*5 :: rtwo(natoms),atwo(natoms)

do j=1,ndna
  rwo(1:3,j)=rt(1:3,j)
  rtwo(j)=atnam(typ(j))
  atwo(j)=atnam(typ(j))
  rnwo(j)=typ(j)
enddo
if (natoms.gt.ndna) then
  k=ndna
  do i=dtype+1,itype
    do j=1,maxntyp(i)
      k=k+1
      rnwo(k)=i
      rtwo(k)=atnam(i)
      atwo(k)=atnam(i)
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
  call writepdb(unitn,natoms,rwo,atwo,rtwo,rnwo)
else
  call writegro(unitn,natoms,rwo,atwo,rtwo,rnwo,box)
endif
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

subroutine readdcdhead(dcdfile,un)
use comun
implicit none
integer*4 kode,un
integer*4 nfile,npriv,nsavc,nstep,nfree
character dcdfile*256

title=''
open(unit=un,file=dcdfile,form='unformatted')
read(un) hdr,icntrl
ntitle=icntrl(20)/12*80
if (allocated(title)) deallocate (title)
allocate (title(ntitle))
read(un,iostat=kode) itemp,title
read(un) natoms
write(*,*) hdr,icntrl
write(*,*) itemp,title
write(*,*) natoms
nfile=icntrl(1)
npriv=icntrl(2)
nsavc=icntrl(3)
nstep=icntrl(4)
if(icntrl(9).gt.0) print *, '# fixed atoms = ',icntrl(9)
nfree = natoms-icntrl(9)
print *, '# of free atoms = ',nfree
print *, 'total # atom = ', natoms,nstep,nsavc
charmm=.false.
if (icntrl(2).eq.0) charmm=.true.
if (nstep.le.0) nstep=1
if (nsavc.le.0) nsavc=1
nsc = nstep/nsavc

allocate (rt(3,natoms))
write(*,'(A,I0)') 'Total number of frames: ',nsc
tnf=nsc
nsc=0
end subroutine

subroutine readdcdbody(un,dcde)
use comun
implicit none
integer*4 kode,i,un
logical*1 dcde
if (charmm) then
  xtlabc6=0d0
  xtlabc12=0d0
else
  if (dcde) then
    read(un,iostat=kode) xtlabc12
  else
    read(un,iostat=kode) xtlabc6
  endif
endif
read(un,iostat=kode) (rt(1,i),i=1,natoms)
read(un,iostat=kode) (rt(2,i),i=1,natoms)
read(un,iostat=kode) (rt(3,i),i=1,natoms)
if (kode.eq.0) then
  nsc=nsc+1
else
  close(un)
endif
end subroutine

subroutine writebtrhead(outfile,unitw,nframew,itypew,atnamw)
implicit none
integer*4 j,kode,unitw
integer*4 nframew,itypew
character*4 atnamw(itypew)
character outfile*(*)

open(unit=unitw,file=trim(outfile),form='unformatted',iostat=kode,status='replace',action='write')
write(unitw) nframew 
write(unitw) itypew 
write(unitw) (atnamw(j),j=1,itypew)  ! ion and nucleotides types in char
end subroutine

subroutine writebtrbody(unitw,runtimew,ntopw,typw,rtw)
implicit none
real*8 runtimew
integer*4 i,ntopw,unitw
integer*4 typw(ntopw)
real*4 rtw(3,ntopw)

write(unitw) runtimew
write(unitw) ntopw
write(unitw) (typw(i),i=1,ntopw)
write(unitw) (rtw(1,i),i=1,ntopw)
write(unitw) (rtw(2,i),i=1,ntopw)
write(unitw) (rtw(3,i),i=1,ntopw)
!nscw=nscw+1
end subroutine
