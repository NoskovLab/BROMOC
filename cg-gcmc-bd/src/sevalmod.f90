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

module sevalmod

contains

  function seval(n, u, x, y, b, c, d)
  !THIS FUNCTION EVALUATES CUBIC SPLINE      
  implicit none
  integer n
  real u, x(n), y(n), b(n), c(n), d(n)
  integer i, j, k
  real dx,seval
  data i/1/
  if ( i .ge. n ) i = 1
  if ( u .lt. x(i) ) go to 10
  if ( u .le. x(i+1) ) go to 30
  10 i = 1
  j = n+1
  20 k = (i+j)/2
  if ( u .lt. x(k) ) j = k
  if ( u .ge. x(k) ) i = k
  if ( j .gt. i+1 ) go to 20
  30 dx = u - x(i)
  seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
  return
  end function
  
  function sevald (n,u,x,b,c,d)
  !THIS FUNCTION EVALUATES 1-ST DERIVATIVE      
  implicit none
  integer n
  real u, x(n), b(n), c(n), d(n)
  integer i, j, k
  real dx,sevald
  data i/1/
  if ( i .ge. n ) i = 1
  if ( u .lt. x(i) ) go to 10
  if ( u .le. x(i+1) ) go to 30
  10 i = 1
  j = n+1
  20 k = (i+j)/2
  if ( u .lt. x(k) ) j = k
  if ( u .ge. x(k) ) i = k
  if ( j .gt. i+1 ) go to 20
  30 dx = u - x(i)
  sevald = b(i) + dx*(c(i) + dx*d(i))
  return
  end function

end module
