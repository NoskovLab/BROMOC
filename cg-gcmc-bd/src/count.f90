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

SUBROUTINE COUNT
!     Count all the ions and assign them to their appropriate buffer
!     nbuffer: the number of buffers in the system
!     ntot: number of particles in the system
!     buffer number ib concerns ions of type ibfftyp(ib)
!     particle i belongs to the buffer number ibuffer(i) 
!     (note that if ibuffer(i) is zero, then particle does not belong to any buffer
!      or is a fix ion)
use grandmod
use constamod
use stdiomod
use errormod
use nucleotmod
implicit none

integer ib,itype,i
real  radius2
logical*1 endlog

!     Initializations
do ib = 1, nbuffer
   nat(ib) = 0 
enddo
do i = nsites+1, ntot
! To which buffer does a particle of type (i) at this location belongs?
  ibuffer(i) = 0
  if (i.gt.(nsites+nfix)) then
    itype = abs(typei(i))
    endlog = .false.
    ib = 1
    do while (ib.le.nbuffer .and. .not.endlog)    
      if (itype.eq.ibfftyp(ib)) then
        if (z(i).ge.LZmin(ib) .and. z(i).le.LZmax(ib)) then
          if (Qsphere) then
            radius2 = (x(i)-cx)**2+(y(i)-cy)**2+(z(i)-cz)**2
            if (radius2.gt.Rsphe2) call error ('count', 'particles outside main system', faterr)      
            if ((radius2.ge.Rmin(ib)**2).and.(radius2.le.Rmax(ib)**2)) then
              ibuffer(i) = ib
              nat(ib) = nat(ib) + 1
              if (z(i).lt.cz) then
                typei(i) = -abs(typei(i))
              else
                typei(i) =  abs(typei(i))
              endif
               endlog = .true.
            endif
          else
            if (Qecyl) then
              radius2=((x(i)-cx)*iecx)**2+((y(i)-cy)*iecy)**2
              if (radius2.gt.1.0) call error ('count', 'particle/s outside main system', faterr)
            else
              if (x(i).lt.lx2m.or.x(i).gt.lx2p) call error ('count', 'particles outside main system', faterr)
              if (y(i).lt.ly2m.or.y(i).gt.ly2p) call error ('count', 'particles outside main system', faterr)
            endif
            if (z(i).lt.lz2m.or.z(i).gt.lz2p) call error ('count', 'particles outside main system', faterr)
            ibuffer(i) = ib
            nat(ib) = nat(ib) + 1
            if (z(i).lt.cz) then
              typei(i) = -abs(typei(i))
            else
              typei(i) =  abs(typei(i))
            endif
            endlog = .true.
          endif
        endif
      endif
      ib = ib + 1
    enddo
  endif   
enddo

do ib = 1, nbuffer
  ntotat(ib) = ntotat(ib) + nat(ib)
enddo

return
end
