!    CHECKBTR - Check BROMOC trajectories for inconsistencies
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

program checkbtr
use btrlib
use wordmod
implicit none
integer*4 narg,arg
character inpfile*256,bs*8

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

write(*,'(/A)') 'Other Usage: CHECKBTR in1.btr'

! read input trajectory filename
call readarg('Input BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

call readbtrhead(inpfile,14,.false.)
! Print number of frames read
nsc=0  
write(*,'(A,I8$)') 'Frame: ',nsc
if (nsc.le.nframe) call readbtrbody()         ! read next frame
do while (inpopen.and.nsc.le.nframe)          ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  if (nsc.le.nframe) call readbtrbody()       ! read next frame
enddo
if (inpopen) call closereadbtr
call closewritebtr
write(*,'(/A/A)') 'File Ok','Normal termination of CHECKBTR'
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

prname='CHECKBTR'
prver='version 1.0'
prdesc='Check BROMOC trajectories for inconsistencies'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='03 Oct 2012'
lastdate='03 Oct 2012'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

