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

      subroutine misc(com,strg,qstrg,wrd5)
c     This subroutine takes care of all the control statements
c     the internal parameters.
      use charfuncmod
      character*(*) com
      character*80 strg(10)
      character*12 label1,label2,wrd12
      character wrd4*4,wrd5*5,fmt*7,fnam*80
      logical*1 qstrg(10),test,error
      include 'mainio.fcm'
 
      integer*4 sign,istrg1,istrg2,istrg3,new_unit,lfnam
      real*8  rstrg1,rstrg2,rstrg3

c     Control statements
c     ------------------

  100 format(' Parameter #',i2,' set to "',a,'"',/)

      if    (wrd5.eq.'set  ') then
c           -----------------
        call getfirst(com,wrd4)
        if(ptrim(wrd4)) istrg=icnvrt(wrd4,error)
        if(.not.error.and.ptrim(com))then
          lstrg=len_trim(com)
          strg(istrg)=com
          com=com(lstrg+1:)
          qstrg(istrg)=isreal(strg(istrg))
          write(outu,100) istrg,trim(strg(istrg))
        endif
c
      elseif (wrd5.eq.'forma') then
c            ---------------
        if(ptrim(com))then
          frmt=com(1:7)
          write(outu,*) 'Get format ',frmt
          com=' '
        else
          frmt=' '
          write(outu,*) 'Format not understood'
        endif
c
      elseif (wrd5.eq.'incr'.or.wrd5.eq.'decr') then
c            ----------------      ----------------
        sign=1
        if(wrd5.eq.'decr') sign=-1
        call getfirst(com,wrd4)
        call getwrd(com,'by',wrd12)
        if(ptrim(wrd4)) istrg=icnvrt(wrd4,error)
        if(.not.error)then
  
          if(qstrg(istrg))then !check for real or integer?
            rstrg1=0.0
            rstrg2=0.0
            rstrg3=0.0
            if(ptrim(strg(istrg))) rstrg1=fcnvrt(strg(istrg))
            if(ptrim(wrd12)) rstrg2=fcnvrt(wrd12)
            rstrg3=rstrg1+sign*rstrg2
            if(frmt.eq.' ')then
               call fndfrmt(rstrg3,fmt) !find a proper format
            else
               fmt=frmt  !user given format
            endif   
            write(com,fmt) rstrg3
            if(ptrim(com)) strg(istrg)=com
            com=' '
          else
            istrg1=0
            istrg2=0
            istrg3=0
            if(ptrim(strg(istrg))) istrg1=icnvrt(strg(istrg),error)
            if(ptrim(wrd12)) istrg2=icnvrt(wrd12,error)
            istrg3=istrg1+sign*istrg2
            write(com,'(i25)') istrg3
            if (ptrim(com)) then 
              strg(istrg)=com
              com=''
            endif
            write(outu,100) istrg,trim(strg(istrg))
          endif
        endif
c
      elseif(wrd5.eq.'if   ')then
c           -----------------
      call getfirst(com,wrd4)
      if(ptrim(wrd4)) istrg=icnvrt(wrd4,error)
      if(.not.error)then
      call getfirst(com,wrd4)
c     Look in the chain c1 for the word coming 
      if(qstrg(istrg))then
      rstrg1=0.0
      rstrg2=0.0
      call getfirst(com,wrd12)
      if(ptrim(strg(istrg))) rstrg1=fcnvrt(strg(istrg))
      if(ptrim(wrd12)) rstrg2=fcnvrt(wrd12)
      else
      istrg1=0
      istrg2=0
      call getfirst(com,wrd12)
      if(ptrim(strg(istrg))) istrg1=icnvrt(strg(istrg),error)
      if(ptrim(wrd12)) istrg2=icnvrt(wrd12,error)
      rstrg1=istrg1
      rstrg2=istrg2
      endif
      endif

      test=.false.
      if(wrd4.eq.'.eq.')then
      test=rstrg1.eq.rstrg2
      elseif(wrd4.eq.'.ge.')then
      test=rstrg1.ge.rstrg2
      elseif(wrd4.eq.'.gt.')then
      test=rstrg1.gt.rstrg2
      elseif(wrd4.eq.'.le.')then
      test=rstrg1.le.rstrg2
      elseif(wrd4.eq.'.lt.')then
      test=rstrg1.lt.rstrg2
      endif
      test=test.and.(.not.error)

        if(test)then
        write(outu,102) 
  102   format(' If (test is .true.) then execute command',/)
        else
        write(outu,103) 
  103   format(' If (test is .false.) then ignore command',/)
        com=' ' !destroy the command
        endif


      else
c     ----
      return !no change mode to wrd5

      endif

c     Re-establish the command line
      wrd5=' '
      call getfirst(com,wrd5)

      return
      end subroutine misc

c============================================================================

!      subroutine subst(com,strg,outun)
!c     This subroutine takes care of the substitutions of internal parameters
!      use charfuncmod
!      character*(*) com
!      character*80 strg(10),wrd80
!      logical*1 error
!      integer*4 outun
!      include 'mainio.fcm'
!
!c     Make the substitutions first
!      istrg=0
! 1000 i1=index(com,'@')
!      if(i1.ne.0)then
!      istrg=icnvrt(com(i1+1:i1+1),error)
!      if(.not.error)then
!      lstrg=len_trim(strg(istrg))
!      lencom=len(com)
!      wrd80=strg(istrg)
!      com(i1:lencom)=trim(wrd80)//com(i1+2:lencom)
!      if(outun.gt.0) write(outun,100) istrg,trim(wrd80)
!  100 format(' Parameter #',i2,' substituted by "',a,'"')
!      goto 1000  !make sure all parameters have been changed
!      endif
!      endif
!      if(outun*istrg.gt.0) write(outun,*)
!      return
!      end subroutine subst

c============================================================================

!      subroutine gtlabel(com,wrd5)
!c     This subroutine control the movement to the labels
!      character*(*) com
!      character label1*12,label2*12,wrd5*5
!      include 'mainio.fcm'
!      call getfirst(com,label1)
!      rewind(unit=inpu,err=1000)
!
! 2000 read(inpu,100,end=3000) com
!  100 format(a)
!      call clean(com)
!      call getfirst(com,wrd5)
!      if(wrd5.eq.'label')then
!      call getfirst(com,label2)
!      if(label1.eq.label2)then
!      wrd5=' '
!      return
!      endif
!      endif
!      goto 2000
!
! 1000 write(outu,101) inpu
!  101 format(' error, could not rewind unit',i3)  
!      return
! 3000 com='*END*'
!      return
!      end subroutine gtlabel

      subroutine fndfrmt(x,fmt) 
C     This subroutine finds a proper format for the real number x
      real*8 x
      character*(*) fmt
      integer*4 logx
      logx=5

      if(x.ne.0.0d0) logx=int(dlog10(dabs(x)))

c     Medium range number
      if((logx.le.5).and.(logx.ge.-5))then
      fmt='(f12.5)'
      else
c     very large or very small numbers
      fmt='(e12.5)'

      endif
      return
      end subroutine fndfrmt
