!    CATDCD - Concatenates several DCD trajectories into one
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

program catdcd
use dcdlib
use wordmod
implicit none
integer*4 narg,arg,nfiles,rfr,fac
character inpfile*256,outfile*256,bs*8,line*256,ext*64,ext2*64
! dcd/dcde
real*8 xtlabc12(12),xtlabc6(6)
real*4,allocatable :: x(:,:)
integer*4 icntrl(20),itemp,icntrl2(4)
integer*4 na,ntitle,nsc,tnf,nat
character hdr*4
character*1,allocatable :: title(:)
logical*1 charmm,dcdopen

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

write(*,'(/A)') 'Other Usage: CATDCD savingfrequency out.dcd in1.dcd in2.dcd'

!saving frequency
call readarg('Saving Frequency: ',narg,arg,line)
fac=chr2int(line)

! read output trajectory filename
call readarg('Output DCD/DCDE trajectory (.dcd/.dcde) filename: ',narg,arg,outfile)
call exten(outfile,ext)

! read input trajectory filename
call readarg('Input DCD/DCDE trajectory (.dcd/.dcde) filename: ',narg,arg,inpfile)
call exten(inpfile,ext2)

call readdcdhead(inpfile,2,ntitle,tnf,hdr,icntrl,itemp,title,na,charmm,dcdopen)
call writedcdhead(outfile,1,tnf,hdr,icntrl,itemp,na)
if (allocated(x)) deallocate (x)
allocate (x(3,na))
nfiles=0
rfr=0
do while (len_trim(inpfile).gt.0)
  if (dcdopen) then
    nfiles=nfiles+1
!
    nsc=0
    write(*,'(A,I8$)') 'Frame: ',nsc
    if (ext2.eq.'dcde') then
      call readdcdebody(2,nsc,xtlabc12,na,x,charmm,dcdopen)
    else
      call readdcdbody(2,nsc,xtlabc6,na,x,charmm,dcdopen)
    endif
    ! Print number of frames read
    do while (dcdopen)          ! open loop for each frame
      write(*,'(A8,I8$)') bs,nsc              ! print frame number
      if (mod(nsc,fac).eq.0) then
        rfr=rfr+1
        if (ext.eq.'dcde') then
          if (ext2.ne.'dcde') call xtl6toxtl12(xtlabc6,xtlabc12)
          call writedcdebody(1,xtlabc12,na,x)
        else
          if (ext2.eq.'dcde') call xtl12toxtl6(xtlabc12,xtlabc6)
          call writedcdbody(1,xtlabc6,na,x)
        endif
      endif
      if (ext.eq.'dcde') then
        call readdcdebody(2,nsc,xtlabc12,na,x,charmm,dcdopen)
      else
        call readdcdbody(2,nsc,xtlabc6,na,x,charmm,dcdopen)
      endif
    enddo
  endif
  inpfile=''
  if (narg.le.2.or.(narg.gt.2.and.arg.lt.narg)) call readarg('Input DCD/DCDE trajectory (.dcd/.dcde) filename: ',narg,arg,inpfile)
  call exten(inpfile,ext2)
  close(2)
  if (len_trim(inpfile).gt.0) then
    call readdcdhead(inpfile,2,ntitle,tnf,hdr,icntrl,itemp,title,nat,charmm,dcdopen)
    if(nat.ne.na) stop 'Number of atoms differ from previous file'
  endif
enddo
write(*,'(/A,I0)') 'Number of input files: ',nfiles 

close(1)

write(*,'(/A,I0)') 'Frames written to output: ',rfr
open(unit=15,file=outfile,form='unformatted',access='stream',action='readwrite',status='old')
icntrl2(1)=rfr
icntrl2(2)=1
icntrl2(3)=1
icntrl2(4)=rfr
write(unit=15,pos=9) icntrl2
close(unit=15)
write(*,'(/A)') 'Normal termination of CATDCD'
end program

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

prname='CATDCD'
prver='version 1.0'
prdesc='Concatenates several DCD trajectories into one'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 May 2013'
lastdate='01 May 2013'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

subroutine exten(filename,ext)
use wordmod
implicit none
character*(*) filename,ext
ext=lcase(filename(index(filename,'.',back=.true.)+1:len_trim(filename)))
end subroutine
