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

	real*8 function cputime(time0)
c
	dimension tarray(2)
c should return elapsed system time
c since time0 as real number of seconds
c
c	cputime = your_system_timer
c eg cray:
c	cputime = rtc()/1.e9 - time0
c vax/convex:
c	cputime = secnds(time0)
C g77	
	result=etime(tarray)
	cputime = tarray(1)
	end
