!    PB-PNP - Poisson-Boltzmann and Poisson-Nernst-Planck Equations Solver
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

module charfuncmod
implicit none
contains
  ! check if a single char is a number 
  function isnum(chr)
  implicit none
  integer*4 i
  character*1 chr
  character*10,parameter :: nums='0123456789'
  logical*1 isnum
  
  isnum=.false.
  
  do i=1,10
    if (chr.eq.nums(i:i)) then
      isnum=.true.
      return
    endif
  enddo
  end function isnum

  ! convert integer*4 to character
  function num2str(n,mx)
  implicit none
  integer*4 n,num,i,a,mx
  character num2str*(mx),numero*10
  numero='0123456789'
  num=n
  do i=mx,1,-1
    a=int(num/(10**(i-1)))
    num2str(mx-i+1:mx-i+1)=numero(a+1:a+1)
    num=num-a*10**(i-1)
  enddo
  end function num2str

  ! convert to lowercase
  function lcase(inchar)
  implicit none
  integer*4 i
  integer*1 s
  character ( len = * ) inchar
  character lcase*(len_trim(inchar))
  lcase=trim(inchar)
  do i=1,len_trim(inchar)
    s=iachar(lcase(i:i))
    if (s.ge.65.and.s.le.90) lcase(i:i)=achar(s+32)
  enddo
  end function lcase
  
  ! convert to uppercase
  function ucase(inchar)
  implicit none
  integer*4 i
  integer*1 s
  character ( len = * ) inchar
  character ucase*(len_trim(inchar))
  ucase=trim(inchar)
  do i=1,len_trim(inchar)
    s=iachar(ucase(i:i))
    if (s.ge.97.and.s.le.122) ucase(i:i)=achar(s-32)
  enddo
  end function ucase
  
  ! convert character to integer
  function chr2int(str)
  implicit none
  integer*4 chr2int,kode
  character str*(*)
  read(str,*,iostat=kode) chr2int
  if (kode.ne.0) stop 'Not an integer'
  end function chr2int
  
  function icnvrt(str,error)
  implicit none
  integer*4 icnvrt,kode
  character str*(*)
  logical*1 error
  read(str,*,iostat=kode) icnvrt
  error=kode.ne.0
  end function icnvrt
 
  ! convert character to real8
  function chr2real(str)
  implicit none
  integer*4 kode
  real*8 chr2real
  character str*(*)
  read(str,*,iostat=kode) chr2real
  if (kode.ne.0) stop 'Not a real'
  end function chr2real
 
  function fcnvrt(str)
  implicit none
  real*8 fcnvrt
  character str*(*)
  fcnvrt=chr2real(str)
  end function fcnvrt
 
  ! get parameter pn from line
  function getparm(str,num,llim,ulim,pn)
  implicit none
  integer*4 num,pn
  integer*4 llim(*),ulim(*)
  character getparm*(ulim(pn)-llim(pn)+1),str*(*)
  if (pn.gt.num.or.pn.lt.1.or.num.lt.1) then
    getparm=''
  else
    getparm=str(llim(pn):ulim(pn))
  endif
  end function getparm

  function isint(str)
  implicit none
  character str*(*)
  integer*4 test,kode
  logical*1 isint
  read(str,*,iostat=kode) test
  if (kode.eq.0) then
    isint=.true.
  else
    isint=.false.
  endif
  end function isint
  
  function isreal(str)
  implicit none
  character str*(*)
  integer*4 kode
  real*8 test
  logical*1 isreal
  read(str,*,iostat=kode) test
  if (kode.eq.0) then
    isreal=.true.
  else
    isreal=.false.
  endif
  end function isreal

  function isword(str)
  implicit none
  integer*4 i,s
  character str*(*)
  logical*1 isword
  isword=.true.
  if (len_trim(str).eq.0) isword=.false.
  do i=1,len_trim(str)
    s=iachar(str(i:i))
    if (.not.((s.ge.65.and.s.le.90).or.(s.ge.97.and.s.le.122))) then
      isword=.false.
      return
    endif
  enddo
  end function isword

  function iswordnum(str)
  implicit none
  integer*4 i,s
  character str*(*)
  logical*1 iswordnum
  iswordnum=.true.
  if (len_trim(str).eq.0) iswordnum=.false.
  do i=1,len_trim(str)
    s=iachar(str(i:i))
    if (.not.((s.ge.65.and.s.le.90).or.(s.ge.97.and.s.le.122).or.(s.ge.48.and.s.le.57))) then
      iswordnum=.false.
      return
    endif
  enddo
  end function iswordnum
  
  function isalfa(str)
  implicit none
  integer*4 i,s
  character str*(*)
  logical*1 isalfa
  isalfa=.true.
  if (len_trim(str).eq.0) isalfa=.false.
  do i=1,len_trim(str)
    s=iachar(str(i:i))
    if (s.le.32.or.s.ge.127) then
      isalfa=.false.
      return
    endif
  enddo
  end function isalfa

  function check(stfin,stin)
  implicit none
  character*(*) stfin,stin
  character stfout*(len(stfin))
  integer, dimension(len_trim(stfin)) :: ll,ul
  integer*4 n,i,j
  logical*1 check
  check=.false.
  stfout=''
  call findparm(stfin,n,ll,ul)
  do i=1,n
    if (lcase(getparm(stfin,n,ll,ul,i)).eq.lcase(adjustl(stin))) then
      check=.true.
      do j=1,n
        if (j.ne.i) stfout=trim(stfout)//' '//getparm(stfin,n,ll,ul,j)
      enddo
      stfin=adjustl(stfout)
      return 
    endif
  enddo
  end function check

  ! reads next line skipping blank lines. True if end is not reached
  function setline(u,line)
  implicit none
  integer*4 u,kode
  character*(*) line
  logical*1 setline
  read(u,'(a)',iostat=kode) line
  do while (len_trim(line).eq.0.and.kode.eq.0)
    read(u,'(a)',iostat=kode) line
  enddo
  setline=kode.eq.0
  end function setline

  ! reads the first word. If no word returns setword false
  function setword(word, line)
  implicit none
  character*(*) word,line
  logical*1 setword
  call getfirst(line,word)
  setword=len_trim(word).ne.0
  end function setword
  
  function setint(intvar,line)
  implicit none
  character*(*) line
  character*(len_trim(line)) word
  integer*4 intvar
  logical*1 setint
  call getfirst(line,word)
  setint=isint(word)
  if (setint) intvar=chr2int(word)
  end function setint

  function ptrim(str)
  implicit none
  character str*(*)
  logical*1 ptrim
  str=adjustl(str)
  ptrim=len_trim(str).gt.0
  end function ptrim

end module charfuncmod
