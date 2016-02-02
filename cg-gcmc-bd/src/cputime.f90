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

function timer()
implicit none
integer*8 ta(8),timer
integer*8 dpm(12),dpy
dpm(1)=0
dpm(2)=31+dpm(1)
dpm(3)=28+dpm(2)
call date_and_time(values=ta)
if (mod(ta(1),4).eq.0) dpm(3)=dpm(3)+1
dpm(4)=31+dpm(3)
dpm(5)=30+dpm(4)
dpm(6)=31+dpm(5)
dpm(7)=30+dpm(6)
dpm(8)=31+dpm(7)
dpm(9)=31+dpm(8)
dpm(10)=30+dpm(9)
dpm(11)=31+dpm(10)
dpm(12)=30+dpm(11)
dpy=int((ta(1)-1)/4)*4*365.25+mod(ta(1)-1,4)*365
timer=((((dpy+dpm(ta(2))+(ta(3)-1))*24+ta(5))*60+ta(6))*60+ta(7))*1000+ta(8)
end function

