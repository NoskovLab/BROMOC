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

program simul
!
!.....simul - main program unit.
!
!
use errormod
use stdiomod
use strtoolsmod
!
!.....Local variables:
!
implicit none
integer           marg, mflg
parameter         (marg = 2, mflg = 1)
character*80      arg(marg)
integer         i,la,nios
logical*1           getopt, larg(marg), lflg(mflg)
data (larg(i),i=1,marg) /marg*.false./
data (lflg(i),i=1,mflg) /mflg*.false./
!
!.....initialize strtools in gnu compiler
!
!*mdc*if gnu

!.....strtools.ini - explicit initialization of character constants.

null       = char(0)
newline    = char(10)
eof        = char(04)
eol        = newline//null
qualifier  = '-'
bell       = char(7)
backspace  = char(8)
tab        = char(9)
linefeed   = char(10)
formfeed   = char(12)
cr         = char(13)
esc        = char(27)
blank      = char(32)
apostrophe = char(39)
dquote     = char(34)
quote      = apostrophe

!*mdc*endif      
!
!.....parse command line
!
if (.not. getopt(arg, larg, marg, 'h', lflg, mflg, null, 0, null)) &
  call error ('simul','run program as: dna1 [input [output]] [-h]', faterr)
!
!.....help:
!
if (lflg(1)) then
   call header (stdout)
   call help (stdout)
   stop 0
endif
!
!.....open input and output files:
!
call ioinit ()
if (larg(1)) then
   call lualloc (inpu)
   la = len_trim(arg(1))
   open (inpu, file = arg(1)(1:la), status = 'old', iostat = nios)
   if (nios .ne. 0) call error ('simul', 'error while opening '//arg(1)(1:la), faterr)
else
   inpu = stdin
endif
if (larg(2)) then
   call lualloc(outu)
   la = len_trim(arg(2))
   open (outu, file = arg(2)(1:la), status='unknown', iostat=nios)
   if (nios .ne. 0) call error ('simul', 'error while opening '//arg(2)(1:la), faterr)
else
   outu = stdout
endif
luwrite = outu
!
!.....read & interpret input file:
!
call shell_simul
!
!.....end of run:
!
if (inpu .ne. stdin) then
   close (inpu)
   call lunalloc (inpu)
endif
if (outu .ne. stdout) then
   close (outu)
   call lunalloc (outu)
endif
end
