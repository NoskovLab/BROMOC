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

subroutine spline (n, x, y, b, c, d)
implicit none
integer n
real x(n), y(n), b(n), c(n), d(n)
!The coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!for a cubic interpolating spline
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) .le. x .le. x(i+1)
!input..
!  n = the number of data points or knots (n.ge.2)
!  x = the abscissas of the knots in strictly increasing order
!  y = the ordinates of the knots
!output..
!  b, c, d  = arrays of spline coefficients as defined above.
!using  p  to denote differentiation,
!  y(i) = s(x(i))
!  b(i) = sp(x(i))
!  c(i) = spp(x(i))/2
!  d(i) = sppp(x(i))/6  (derivative from the right)
!the accompanying function subprogram  seval  can be used
!to evaluate the spline.
integer nm1, ib, i
real t
nm1 = n-1
if ( n .lt. 2 ) return
if ( n .gt. 2 ) then
  !set up tridiagonal system
  !b = diagonal, d = offdiagonal, c = right hand side.
  d(1) = x(2) - x(1)
  c(2) = (y(2) - y(1))/d(1)
  do i = 2, nm1
    d(i) = x(i+1) - x(i)
    b(i) = 2.0*(d(i-1) + d(i))
    c(i+1) = (y(i+1) - y(i))/d(i)
    c(i) = c(i+1) - c(i)
  enddo
  !end conditions.  third derivatives at  x(1)  and  x(n)
  !obtained from divided differences
  b(1) = -d(1)
  b(n) = -d(n-1)
  c(1) = 0.0
  c(n) = 0.0
  if ( n .gt. 3 ) then
    c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
    c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
    c(1) = c(1)*d(1)**2/(x(4)-x(1))
    c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
  endif 
  !forward elimination
  do i = 2, n
    t = d(i-1)/b(i-1)
    b(i) = b(i) - t*d(i-1)
    c(i) = c(i) - t*c(i-1)
  enddo
  !back substitution
  c(n) = c(n)/b(n)
  do ib = 1, nm1
    i = n-ib
    c(i) = (c(i) - d(i)*c(i+1))/b(i)
  enddo
  !c(i) is now the sigma(i) of the text
  !compute polynomial coefficients
  b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0*c(n))
  do i = 1, nm1
    b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
    d(i) = (c(i+1) - c(i))/d(i)
    c(i) = 3.0*c(i)
  enddo
  c(n) = 3.0*c(n)
  d(n) = d(n-1)
else 
  b(1) = (y(2)-y(1))/(x(2)-x(1))
  c(1) = 0.0
  d(1) = 0.0
  b(2) = b(1)
  c(2) = 0.0
  d(2) = 0.0
endif
return
end
