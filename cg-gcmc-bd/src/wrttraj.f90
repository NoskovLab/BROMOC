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

subroutine wrttraj
use grandmod
use stdiomod
use nucleotmod
use charfuncmod, only: sng   !command parser

implicit none
integer i

if (Qpar) then
  do i = nsites+1, ntot
     typtyp(i)=nwtype(abs(typei(i)))
  enddo
endif
write(iuntrj) runtime                    ! simulation time in pico-second (before was nano)
write(iuntrj) ntot                       ! total number of ions and sites in motion
write(iuntrj) (typtyp(i),i=1,ntot)       ! ion and nucleotides types
write(iuntrj) (sng(x(i)),i=1,ntot)       ! x coordinates
write(iuntrj) (sng(y(i)),i=1,ntot)       ! y coordinates
write(iuntrj) (sng(z(i)),i=1,ntot)       ! z coordinates 
return
end
