!    WORDMOD Word Handler Library/Module
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

module wordmod
implicit none
contains
  ! convert integer to character
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

  ! convert to lowercase
  function lcase(inchar)
  implicit none
  integer i
  integer*1 s
  character ( len = * ) inchar
  character ( len = len_trim(inchar) ) lcase
  lcase=trim(inchar)
  do i=1,len(lcase)
    s=iachar(lcase(i:i))
    if (s.ge.65.and.s.le.90) lcase(i:i)=achar(s+32)
  enddo
  end function
  
  ! convert to uppercase
  function ucase(inchar)
  implicit none
  integer i
  integer*1 s
  character ( len = * ) inchar
  character ( len = len_trim(inchar) ) ucase
  ucase=trim(inchar)
  do i=1,len(ucase)
    s=iachar(ucase(i:i))
    if (s.ge.97.and.s.le.122) ucase(i:i)=achar(s-32)
  enddo
  end function
  
  ! convert character to integer
  function chr2int(str)
  implicit none
  integer chr2int,kode
  character str*(*)
  read(str,*,iostat=kode) chr2int
  if (kode.ne.0) stop 'Not an integer'
  end function

  ! convert character to real8
  function chr2real(str)
  implicit none
  integer kode
  real*8 chr2real
  character str*(*)
  read(str,*,iostat=kode) chr2real
  if (kode.ne.0) stop 'Not a real'
  end function
  
  ! get parameter pn from line
  function getparm(str,num,llim,ulim,pn)
  implicit none
  integer num,pn
  integer llim(*),ulim(*)
  character getparm*(ulim(pn)-llim(pn)+1),str*(*)
  if (pn.gt.num.or.pn.lt.1.or.num.lt.1) then
    getparm=''
  else
    getparm=str(llim(pn):ulim(pn))
  endif
  end function

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
  end function
  
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
  end function

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
  end function

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
  end function
  
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
  end function

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
  
  subroutine countparm(str,num)
  implicit none
  integer i,length,num
  character str*(*)
  logical chng
  length=len_trim(str)
  chng=.false.
  num=0
  do i=1,length
    if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
      chng=.false.
    else
      if (.not.chng) num=num+1
      chng=.true.
    endif
  enddo
  end subroutine
 
end module


