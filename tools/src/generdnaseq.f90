!    GENERDNASEQ - Generates Random DNA Sequence. Uses CPU clock to generate randomness.
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

program generdnaseq
!use ifport
implicit none
integer*4 arg,narg,nucl,i,iseed,cont(4),tries
real*4 rn
character line*256,base*1
logical*1 ok
call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()
call readarg('Number of nucleotides [10]: ',narg,arg,line)
if (len_trim(line).eq.0) line='10'
read(line,*) nucl

iseed=0
call askclock(iseed)
write(line,'(I0)') iseed
call readarg('Seed ['//trim(line)//']: ',narg,arg,line)
if (len_trim(line).gt.0.and.trim(line).ne.'0') read(line,*) iseed
call srand(iseed)

write(*,'(//A/)') 'Sequence:'

tries=0
ok=.true.
do while (ok)
  tries=tries+1
  cont=0
  do i=1,nucl
    rn=4e0*rand()
    if (rn.lt.1e0) then
      base='A'
      cont(1)=cont(1)+1
    elseif (rn.lt.2e0) then
      base='C'
      cont(2)=cont(2)+1
    elseif (rn.lt.3e0) then
      base='G'
      cont(3)=cont(3)+1
    elseif (rn.lt.4e0) then 
      cont(4)=cont(4)+1
      base='T'
    endif
    write(*,'(A$)') base
  enddo
  write(*,'(I10,4(A,I0))') tries,' A=',cont(1),' C=',cont(2),' G=',cont(3),' T=',cont(4)
  if (mod(nucl,4).eq.0) then
    ok=.false.
    do i=2,4
      if (cont(i).ne.cont(1)) then
        ok=.true.
        exit
      endif
    enddo
  else
    ok=.false.
  endif
enddo
!write(*,'(//A,I0)') '# of A = ',cont(1)
!write(*,'(A,I0)') '# of C = ',cont(2)
!write(*,'(A,I0)') '# of G = ',cont(3)
!write(*,'(A,I0)') '# of T = ',cont(4)
!
!write(*,'(/A,I0)') 'Tries: ',tries
write(*,'(/A/)') 'Normal Termination.'

end program

subroutine askclock(iseed)
implicit none
integer*4 iseed
if (iseed.le.0.or.iseed.ge.2147483647) call system_clock(count=iseed)
end subroutine


! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='GENERDNASEQ'
prver='version 1.0'
prdesc='Generates Random DNA Sequence. Uses CPU clock to generate randomness.'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='19 Oct 2012'
lastdate='19 Oct 2012'

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


