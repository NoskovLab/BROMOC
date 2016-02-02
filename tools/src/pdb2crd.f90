!    PDB2CRD - Converts PDB to CRD (CHARMM Coordinates)
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

program pdb2crd
implicit none
! pdb and crd
real*8 :: w,e
integer*4 ::  na,ires
character*4 :: typ,res,segid,resid,iatom
real*8 :: rt(3)

integer*4 i,j,k,narg,arg,kode
character inpfile*256,outfile*256,ln*256,a6*6
logical once,exton

call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()
call readarg('Input PDB filename: ',narg,arg,inpfile)
call readarg('Output CHARMM CRD filename: ',narg,arg,outfile)
call readarg('CHARMM Extended coordinates format (y/n) [n]? ',narg,arg,ln)
if (ln(1:1).eq.'y'.or.ln(1:1).eq.'Y') then
  exton=.true.
else
  exton=.false.
endif

open(unit=1,file=inpfile,IOSTAT=kode)
j=0
once = .true.
read(1,'(A)',IOSTAT=kode) a6
do while (kode.eq.0.or.kode.eq.64.and.once)
  if (a6(1:3).eq.'END') then
    na=j
    once = .false.
  elseif ((a6(1:4).eq.'ATOM'.or.a6(4:6).eq.'ATM').and.once) then
    j=j+1
  endif
  read(1,'(A)',IOSTAT=kode) a6
enddo
if (once) then
  write(*,*) 'WARNING: No END keyword'
  na=j
endif
write(*,*) 'Number of atoms: ',na
rewind(1)

open(unit=2,file=outfile)
if (exton) then
  write(2,'(I10,A)') na,'  EXT'
else
  write(2,'(I5)') na
endif

j=0
once = .true.
read(1,'(A)',IOSTAT=kode) ln
do while (kode.eq.0.or.kode.eq.64.and.once)
  if (ln(1:3).eq.'END') then
    once = .false.
  elseif ((ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM').and.once) then
    j=j+1
    read (ln,'(6x,4(x,A4),4x,3F8.3,2F6.2,6x,A4)') iatom,typ,res,resid,rt(1),rt(2),rt(3),e,w,segid
    read(resid,*) ires  
    if (exton) then 
      write(2,'(I10,I10,2x,A4,6X,A4,4x,3F20.10,2X,A4,6X,A4,4x,F20.10)') j,ires,res,typ,(rt(i),i=1,3),segid,resid,w
    else  
      write(2,'(I5,I5,1x,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5)') j,ires,res,typ,(rt(i),i=1,3),segid,resid,w
    endif
  endif
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
close(2)
write(*,'(/A)') 'Normal termination.'
end program

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='PDB2CRD'
prver='version 1.0'
prdesc='Converts PDB to CHARMM CRD'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='21 Aug 2012'
lastdate='21 Aug 2012'

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

