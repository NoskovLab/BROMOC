!    CATBTR - Concatenates several BROMOC trajectories into one
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

program catbtr
use btrlib
use wordmod
implicit none
real*8 runtimecum
integer*4 i,narg,arg,nfiles,rfr
character inpfile*256,outfile*256,bs*8,line*256
character*256,allocatable ::ifeven(:),ifodd(:),inpfiles(:)
integer*4 fac,counter

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

write(*,'(/A)') 'Other Usage: CATBTR savingfrequency out.btr in1.btr in2.btr'

!saving frequency
call readarg('Saving Frequency: ',narg,arg,line)
fac=chr2int(line)

! read output trajectory filename
call readarg('Output BROMOC trajectory (.btr) filename: ',narg,arg,outfile)

! read input trajectory filename
call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)
nframew=0
nfiles=0
do while (len_trim(inpfile).gt.0)
  call readbtrhead(inpfile,14,.false.)
  if (inpopen) then
    nfiles=nfiles+1
    if (nfiles.eq.1) then
      itypew=itype
      allocate (atnamw(itypew))
      atnamw=atnam
    endif
    if (mod(nfiles,2).eq.0) then
      allocate (ifeven(nfiles))
      ifeven(1:nfiles-1)=ifodd(1:nfiles-1)
      ifeven(nfiles)=inpfile
      deallocate (ifodd)
    else
      allocate (ifodd(nfiles))
      ifodd(nfiles)=inpfile
      if (nfiles.gt.1) then
        ifodd(1:nfiles-1)=ifeven(1:nfiles-1)
        deallocate (ifeven)
      endif
    endif
    nframew=nframew+nframe
    call closereadbtr
  endif
  inpfile=''
  if (narg.le.2.or.(narg.gt.2.and.arg.lt.narg)) call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)
enddo
write(*,'(A,I0)') 'Number of input files: ',nfiles 
allocate (inpfiles(nfiles))
if (mod(nfiles,2).eq.0) then
  inpfiles=ifeven
  deallocate (ifeven)
else
  inpfiles=ifodd
  deallocate (ifodd)
endif

counter=0
rfr=0
nscw=0
runtimew=0d0
runtimecum=0d0
call writebtrhead(outfile,15)
write(*,'(A,I8$)') 'Frame: ',counter
do i=1,nfiles
  call readbtrhead(inpfiles(i),14,.true.)
  ! Print number of frames read
  nsc=0  
  if (nsc.le.nframe) call readbtrbody()         ! read next frame
  do while (inpopen.and.nsc.le.nframe)          ! open loop for each frame
    counter=counter+1
    write(*,'(A8,I8$)') bs,counter              ! print frame number
    if (mod(nsc,fac).eq.0) then
      runtimew=runtimecum+runtime
      ntopw=ntop
      allocate (typw(ntopw),rtw(3,ntopw))
      typw=typ
      rtw=rt
      rfr=rfr+1
      call writebtrbody                         ! writebtr
      deallocate (typw,rtw)
    endif
    if (nsc.le.nframe) call readbtrbody()       ! read next frame
  enddo
  if (inpopen) call closereadbtr
  runtimecum=runtimew
enddo
call closewritebtr
write(*,'(/A,I0)') 'Frames written to output: ',rfr
open(unit=15,file=outfile,form='unformatted',access='stream',action='readwrite',status='old')
write(unit=15,pos=5) rfr
close(unit=15)
write(*,'(/A)') 'Normal termination of CATBTR'
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

prname='CATBTR'
prver='version 1.0'
prdesc='Concatenates several BROMOC trajectories into one'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='07 Aug 2012'
lastdate='07 Aug 2012'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

