!    BROMOC  -  CG-GCMC-BD
!    Electrodiffusion, Gran Canonical Monte Carlo, Brownian,Dynamics 
!    and Coarse Grain Model DNA Simulation Program.
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

! CHARACTER MANIPULATION LIBRARY

! identify words from line
subroutine findparm(str,num,llim,ulim)
implicit none
integer i,length,num
integer ulim(*),llim(*)
character str*(*)
logical*1 chng
length=len_trim(str)
if (length.le.0) then
  num=0
  return
endif
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
  
subroutine gtcpar(c1,c2,n,junk)
! This subroutine get the character string parameter "cpar"
! coming just after the word c2 in the chain c1.
! At the end the word and the parameter are removed
! from the c1.  If nothing was there, the default "idef"
! is used.
use charfuncmod
implicit none
character*(*) c1,c2
character*(*) junk
integer n,i,ln
logical*1 numon
junk=''
call getwrd(c1,c2,junk)
n=0
numon=.false.
ln=len_trim(junk)
if (ln.eq.0) then 
  write(*,*) 'Error: Nothing after '//trim(c2)
  stop
endif
do i=1,ln
  if (isnum(junk(i:i)).and..not.numon) then
    n=n+1
    numon=.true.
  elseif (.not.isnum(junk(i:i)).and.numon) then
    numon=.false.
  endif
enddo
if (n.eq.0) then 
  write(*,*) 'Error: No number after '//trim(c2)
  stop
endif
end subroutine

subroutine gtcrpar(c1,c2,n,junk)
! This subroutine get the character string parameter "cpar"
! coming just after the word c2 in the chain c1.
! At the end the word and the parameter are removed
! from the c1.  If nothing was there, the default "idef"
! is used.
use charfuncmod
implicit none
character*(*) c1,c2
character*(*) junk
integer n,i,ln
logical*1 numon
junk=''
call getwrd(c1,c2,junk)
n=0
numon=.false.
ln=len_trim(junk)
if (ln.eq.0) then
  write(*,*) 'Error: Nothing after '//trim(c2)
  stop
endif
do i=1,ln
  if (isrnum(junk(i:i)).and..not.numon) then
    n=n+1
    numon=.true.
  elseif (.not.isrnum(junk(i:i)).and.numon) then
    numon=.false.
  endif
enddo
if (n.eq.0) then
  write(*,*) 'Error: No number after '//trim(c2)
  stop
endif
end subroutine


subroutine getwrd(stfin,stin,stout)
use charfuncmod
implicit none
character*(*) stfin,stin,stout
character stfout*(len(stfin))
integer,dimension(len_trim(stfin)) :: ll,ul
integer n,i,j
stout=''
n=0
if (len_trim(stfin).gt.0) call findparm(stfin,n,ll,ul)
if (n.gt.1) then
  do i=1,n
    if (lcase(getparm(stfin,n,ll,ul,i)).eq.lcase(adjustl(stin))) then
      stout=getparm(stfin,n,ll,ul,i+1)
      stfout=''
      do j=1,n
        if (j.ne.i.and.j.ne.i+1) then 
          stfout=trim(stfout)//' '//getparm(stfin,n,ll,ul,j)
        endif
      enddo
      stfin=adjustl(stfout)
      return 
    endif
  enddo
endif
end subroutine

subroutine gtcipar(junk,iparv)
use charfuncmod
implicit none
character*(*) junk
integer i,j,ki,kf,ln
integer iparv(*)
logical*1 numon

j=0
numon=.false.
junk=trim(adjustl(junk))
ln=len_trim(junk)
ki=1
if (ln.eq.0) stop 'Error: character empty'
do i=1,ln
  if (isnum(junk(i:i)).and..not.numon) then
    j=j+1
    ki=i
    numon=.true.
  endif
  if ((.not.isnum(junk(i:i)).and.numon).or.i.eq.ln) then
    if (i.eq.ln) then
      kf=i
    else
      kf=i-1
    endif
    numon=.false.
    read(junk(ki:kf),*) iparv(j)
  endif
enddo

end subroutine

subroutine gtcdpar(junk,iparv)
use charfuncmod
implicit none
character*(*) junk
integer i,j,ki,kf,ln
real iparv(*)
logical*1 numon

j=0
numon=.false.
junk=trim(adjustl(junk))
ln=len_trim(junk)
ki=1
if (ln.eq.0) stop 'Error: character empty'
do i=1,ln
  if (isrnum(junk(i:i)).and..not.numon) then
    j=j+1
    ki=i
    numon=.true.
  endif
  if ((.not.isrnum(junk(i:i)).and.numon).or.i.eq.ln) then
    if (i.eq.ln) then
      kf=i
    else
      kf=i-1
    endif
    numon=.false.
    read(junk(ki:kf),*) iparv(j)
  endif
enddo
end subroutine

subroutine gtdpar(c1,c2,par,def)
! This subroutine get the real parameter "par"
! coming just after the word c2 in the chain c1.
! At the end the word and the parameter are removed
! from the c1.  If nothing was there, the default "def"
! is used.
use charfuncmod
implicit none
character*(*) c1,c2
character*(len_trim(c1)) junk
real par,def
junk=''
call getwrd(c1,c2,junk)
if(len_trim(junk).gt.0)then
  par=chr2real(junk)
else
  par=def
endif
return
end subroutine

subroutine gtipar(c1,c2,ipar,idef)
! This subroutine get the integer parameter "ipar"
! coming just after the word c2 in the chain c1.
! At the end the word and the parameter are removed
! from the c1.  If nothing was there, the default "idef"
! is used.
use charfuncmod
implicit none
character*(*) c1,c2
character*(len_trim(c1)) junk
integer ipar,idef
junk=''
call getwrd(c1,c2,junk)
if(len_trim(junk).gt.0)then
  ipar=chr2int(junk)
else
  ipar=idef
endif
return

end subroutine

subroutine gti8par(c1,c2,ipar,idef)
! This subroutine get the integer parameter "ipar"
! coming just after the word c2 in the chain c1.
! At the end the word and the parameter are removed
! from the c1.  If nothing was there, the default "idef"
! is used.
use charfuncmod
implicit none
character*(*) c1,c2
character*(len_trim(c1)) junk
integer*8 ipar,idef
junk=''
call getwrd(c1,c2,junk)
if(len_trim(junk).gt.0)then
  ipar=chr2int8(junk)
else
  ipar=idef
endif
return

end subroutine

subroutine clean(c1)
implicit none
!This subroutine removes any comments placed after !
character*(*) c1
integer i1 !,length
!length=len(c1)
i1=index(c1,'!')
if(i1.ne.0) c1=c1(1:i1-1) !c1(1:length)=c1(1:i1-1)
i1=index(c1,'#')
if(i1.ne.0) c1=c1(1:i1-1) !c1(1:length)=c1(1:i1-1)
i1=index(c1,'*')
if(i1.ne.0) c1=c1(1:i1-1) !c1(1:length)=c1(1:i1-1)
return
end subroutine

subroutine getlin(c1,inpu,outu)
use charfuncmod
implicit none
!This subroutine gets a  non-blank line of command in c1, if the 
!line terminates with a "-" then the next line is appended to it 
!and so on, until the maximum number of lines "nline" of len1gth 
!"len1" is read.
character*(*) c1
character*350 line
integer inpu,outu
integer begin,endlin,lenc1
logical*1 done

c1=''
line=''
begin=1
lenc1=len(c1)
endlin=0
done=.false.

1000 read(inpu,'(a)',end=2000) line

call clean(line)
line=adjustl(line)
if (len_trim(line).ne.0) then
  endlin=len_trim(line)
  if(check(line(endlin:),'-'))then
    done=.false.
  else
    done=.true.
  endif
  c1(begin:)=line(1:endlin)
  begin=begin+endlin
  if(begin.gt.lenc1)then
    write(abs(outu),'(a)')' * GETLIN * error, command line truncated'
    done=.true.
  endif
else
  done=.false.
endif

if(done)then
return !that was the last line
else
goto 1000 !append a new line
endif


! END-OF-FILE encountered
2000 c1='*END*'
write(abs(outu),'(a)') ' * GETLIN * error, end of file encountered'
stop
return

end subroutine

subroutine getfirst(stfin,stout)
! get first word in c1 and remove it and put it in c2
use charfuncmod
implicit none
integer n,j
character*(*) stfin,stout
character stfout*(len(stfin))
integer,dimension(len_trim(stfin)) :: ll,ul

stout=''
n=0
if (len_trim(stfin).gt.0) call findparm(stfin,n,ll,ul)
if (n.ge.1) then
  stout=getparm(stfin,n,ll,ul,1)
  stfout=''
  do j=2,n
    stfout=trim(stfout)//' '//getparm(stfin,n,ll,ul,j)
  enddo
  stfin=adjustl(stfout)
endif
end subroutine
