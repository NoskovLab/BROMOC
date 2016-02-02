!    XTCREAD - Reads XTC GROMACS trajectories
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

program xtcread
use wordmod
implicit none
integer*4 narg,arg,nsc,nframe,frame,i,j,k
character inpfile*256,bs*8,line*256
integer*4 xd,nat,step,ret
real*4 time,box(9),prec
real*4,allocatable :: x(:)

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

write(*,'(/A)') 'Other Usage: XTCREAD in.xtc'

! read input trajectory filename
call readarg('Input XTC GROMACS trajectory (.xtc) filename: ',narg,arg,inpfile)

! open
call xtcopen(xd,trim(inpfile),'r',nat)
allocate (x(3*nat))

! count number of frames
nframe=0
write(*,'(/A,I8$)') 'Frame: ',nframe
call readxtc(xd,nat,step,time,box,x,prec,ret)
do while (ret.ne.0)
  nframe=nframe+1
  write(*,'(A8,I8$)') bs,nframe                  ! print frame number
  call readxtc(xd,nat,step,time,box,x,prec,ret)
enddo
call xdrfclose(xd,ret)
write(*,'(//A,I0)') 'Number of read frames: ',nframe
write(line,*) 'Number of particles: ',nat
write(*,'(A)') trim(line)
write(line,*) 'Last Time: ',time
write(*,'(A)') trim(line)
write(line,*) 'Last Step: ',step
write(*,'(A)') trim(line)
write(line,*) 'Last Prec: ',prec
write(*,'(A)') trim(line)

call readarg('Input Frame Number to print info: ',narg,arg,line)
read(line,*) frame

call xtcopen(xd,trim(inpfile),'r',nat)

! Print number of frames read
nsc=0  
write(*,'(/A,I8$)') 'Frame: ',nsc
do while (nsc.lt.frame.and.nsc.lt.nframe)     ! open loop for each frame
  nsc=nsc+1
  write(*,'(A8,I8$)') bs,nsc                  ! print frame number
  call readxtc(xd,nat,step,time,box,x,prec,ret)
  if (nsc.eq.frame) then
    write(line,*) 'Number of particles: ',nat
    write(*,'(//A)') trim(line)
    write(line,*) 'Time: ',time
    write(*,'(A)') trim(line)
    write(line,*) 'Step: ',step
    write(*,'(A)') trim(line)
    write(line,*) 'Prec: ',prec
    write(*,'(A)') trim(line)
    write(*,'(A)') 'In Angstroms: '
    write(line,*) 'Box: ',(box(i)*1e1,i=1,3)
    write(*,'(A)') trim(line)
    write(line,*) '     ',(box(i)*1e1,i=4,6)
    write(*,'(A)') trim(line)
    write(line,*) '     ',(box(i)*1e1,i=7,9)
    write(*,'(A)') trim(line)
    write(*,'(A)') 'Coordinates:'
    k=0
    do i=1,3*nat,3
      k=k+1
      write(line,*) k,(x(j)*1e1,j=i,i+2)
      write(*,'(A)') trim(line)
    enddo
  endif
enddo
call xdrfclose(xd,ret)
deallocate (x)
write(*,'(/A)') 'Normal termination of XTCREAD'
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

prname='XTCREAD'
prver='version 1.0'
prdesc='Reads XTC GROMACS trajectories'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='04 Oct 2012'
lastdate='04 Oct 2012'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

