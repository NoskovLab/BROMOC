!    PDB2CRD - Converts CRD CHARMM/NAMD to NAMD/CHARMM CRD
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

module dcd
implicit none
integer*4 na,nsc,tnf
integer*4 icntrl(20),itemp
character hdr*4
integer ntitle
character*1,allocatable :: title(:)
real*4,allocatable :: rc(:,:)
real*8 :: xtlabc6(6)
logical*1 dcdopen,charmm
end module

program crd2crd
use dcd
implicit none
! crd and crd
real*8 :: w,e
integer*4 ::  nop,ires,ares
character*4 :: typ,res,segid,resid,iatom
real*8 :: rt(3),rs(3),rb(3)

integer*4 i,j,k,narg,arg,kode,ll(256),ul(256),num,cr,ca,pr,rbn,rsn
character inpfile*256,outfile*256,ln*256,a6*6,dcdfile*256
logical once,readexton,writeexton

call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()
dcdopen=.false.
call readarg('Input CRD filename: ',narg,arg,inpfile)
call readarg('Output CRD filename: ',narg,arg,outfile)
call readarg('Write EXTended (e) or Regular (r) CRD [r]? ',narg,arg,ln)
if (ln(1:1).eq.'e'.or.ln(1:1).eq.'E') then
  writeexton=.true.
else
  writeexton=.false.
endif

call readarg('Input DCD filename for new coordinates [ENTER for NONE]: ',narg,arg,dcdfile)
if (len_trim(dcdfile).gt.0) call readdcdhead(dcdfile,3)

open(unit=1,file=inpfile,IOSTAT=kode)
open(unit=2,file=outfile)
open(unit=4,file=trim(outfile)//'.xyz')
once=.true.
do while (once)
  read(1,'(A)') ln
  ln=adjustl(ln)
  if (ln(1:1).ne.'*') then
    call findparm(ln,num,ll,ul)
    if (num.eq.0) stop 'Error reading CRD'
    read(ln(ll(1):ul(1)),*) nop
    readexton=.false.
    if (num.eq.2) then 
      if (ln(ll(2):ul(2)).eq.'EXT') readexton=.true.
    endif
    once=.false.
    if (writeexton) then
      write(2,'(I10,A)') nop,'  EXT'
    else
      write(2,'(I5)') nop
    endif
  else
    write(2,'(A)') trim(ln)
  endif
enddo
if (readexton) then
  write(*,'(/A)') 'Input CRD: EXTended format detected.'
else
  write(*,'(/A)') 'Input CRD: Regular format detected.'
endif
if (writeexton) then
  write(*,'(/A/)') 'Output CRD: EXTended format.'
else
  write(*,'(/A/)') 'Ouput CRD: Regular format.'
endif

write(*,'(A,I0)') 'Number of atoms: ',nop

if (dcdopen) then
  if (nop.ne.na) stop 'not same number of atoms/particles'
  call readdcdbody(3) ! read first frame
  close(3)
  dcdopen=.true. 
endif

ca=0
cr=1
pr=-1
rs=0d0
rb=0d0
rbn=0
rsn=0
do j=1,nop
  if (readexton) then 
    read(1,'(I10,I10,2x,A4,6X,A4,4x,3F20.10,2X,A4,6X,A4,4x,F20.10)') k,ires,res,typ,(rt(i),i=1,3),segid,resid,w
  else
    read(1,'(I5,I5,1x,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)') k,ires,res,typ,(rt(i),i=1,3),segid,resid,w
  endif
  if (j==1) ares=ires
  if (dcdopen) rt=rc(:,j)
  if (ares.ne.ires) then 
    rbn=0
    rsn=0
    rs=0d0
    rb=0d0
    ares=ires
  endif
  if (typ=="C4'".or.typ=="O4'".or.typ=="C1'".or.typ=="C2'".or.typ=="C3'") then
    rsn=rsn+1
    rs=rs+rt
    if (rsn.eq.5) then 
      typ="S"
      rt=rs/5d0
    elseif (rsn.gt.5) then
      write(*,*) 'Possible error more elements for sugar than expected' 
    endif
  elseif ((typ(1:1)=='N'.or.typ(1:1)=='O'.or.typ(1:1)=='C').and.typ(3:4)=='  '.or.typ(3:4)=='M ') then
    rb=rb+rt
    rbn=rbn+1
    rt=rb/rbn
    if (res=='ADE'.and.rbn==10) then 
      typ="Ab"
    elseif (res=='CYT'.and.rbn==8) then 
      typ="Cb"
    elseif (res=='GUA'.and.rbn==11) then 
      typ="Gb"
    elseif (res=='THY'.and.rbn==9) then 
      typ="Tb"
    endif
  endif
  if (typ=="P".or.typ=="S".or.typ=="Ab".or.typ=="Cb".or.typ=="Gb".or.typ=="Tb") then
    ca=ca+1
    if (pr.ne.-1) then
      if (ires.ne.pr) cr=cr+1
    endif 
    if (writeexton) then 
      write(2,'(I10,I10,2x,A4,6X,A4,4x,3F20.10,2X,A4,6X,A4,4x,F20.10)') ca,cr,res,typ,(rt(i),i=1,3),segid,resid,w
    else  
      write(2,'(I5,I5,1x,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)') ca,cr,res,typ,(rt(i),i=1,3),segid,resid,w
    endif
    write(4,*) typ,rt(1:3)
    pr=ires
  endif
enddo
close(1)
close(2)
close(4)
write(*,'(/A)') 'Normal termination of DNACDF'
end program

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='PDB2CRD'
prver='version 1.0'
prdesc='Converts CRD CHARMM/NAMD to NAMD/CHARMM CRD'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='24 Oct 2012'
lastdate='24 Oct 2012'

write(*,'(/A)') trim(prname)//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
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

subroutine readdcdhead(dcdfile,un)
use dcd
implicit none
integer*4 kode,un
integer*4 nfile,npriv,nsavc,nstep,nfree
character dcdfile*256

title=''
open(unit=un,file=dcdfile,form='unformatted')
dcdopen=.true.
read(un) hdr,icntrl
ntitle=icntrl(20)/12*80
if (allocated(title)) deallocate (title)
allocate (title(ntitle))
read(un,iostat=kode) itemp,title
read(un) na
write(*,*) hdr,icntrl
write(*,*) itemp,title
write(*,*) na
nfile=icntrl(1)
npriv=icntrl(2)
nsavc=icntrl(3)
nstep=icntrl(4)
if(icntrl(9).gt.0) print *, '# fixed atoms = ',icntrl(9)
nfree = na-icntrl(9)
print *, '# of free atoms = ',nfree
print *, 'total # atom = ', na,nstep,nsavc
charmm=.false.
if (icntrl(2).eq.0) charmm=.true.
if (nstep.le.0) nstep=1
if (nsavc.le.0) nsavc=1
nsc = nstep/nsavc

allocate (rc(3,na))
write(*,'(A,I0)') 'Total number of frames: ',nsc
tnf=nsc
nsc=0
end subroutine

subroutine readdcdbody(un)
use dcd
implicit none
integer*4 kode,i,un
logical*1 dcde
if (charmm) then
  xtlabc6=0d0
else
  read(un,iostat=kode) xtlabc6
endif
read(un,iostat=kode) (rc(1,i),i=1,na)
read(un,iostat=kode) (rc(2,i),i=1,na)
read(un,iostat=kode) (rc(3,i),i=1,na)
if (kode.eq.0) then
  nsc=nsc+1
else
  close(un)
  dcdopen=.false.
endif
end subroutine

! identify words from line
subroutine findparm(str,num,llim,ulim)
implicit none
integer i,length,num
integer ulim(*),llim(*)
character str*(*)
logical chng
length=len_trim(str)
chng=.false.
num=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
      ulim(num)=0
    endif
    chng=.true.
  endif
enddo
if (ulim(num).eq.0) ulim(num)=length
end subroutine


