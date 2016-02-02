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

SUBROUTINE ALPOL2(J,LMAX,MMAX,XC,AP)
!------------------------------------------------------------------------
!The Associate Legendre Polynomials 
!(Numerical recipes in Fortran 77: 6.8. Spherical Harmonics) 
!This is different from FUNCTION ALPOL because we don't consider xc > 1 here.

use grandmod
implicit none

integer lmax,mmax,J
real  xc,ap(ntot,0:LMAX,0:LMAX-2)
!local
integer l,m
real  fact,somx2

if(mmax.gt.lmax) then
   write(*,*) 'MMAX .GT. LMAX'
   stop 'BAD ARGUMENTS IN ALPOL'
elseif(mmax.lt.0) then
   stop 'BAD ARGUMENTS IN ALPOL'
endif

!compute p(m,m)
AP(J,0,0)=1.0
somx2=sqrt((1.0-xc)*(1.0+xc))
fact=1.0
do m=1,mmax
   AP(J,m,m)=AP(J,m-1,m-1)*fact*somx2
   fact=fact+2.0
enddo

!compute p(m+1,m)
do m=0,mmax
   AP(J,m+1,m)=xc*(2*m+1)*AP(J,m,m)
enddo

!compute p(l,m)
do l=2,lmax
   do m=0,l-2
      AP(J,l,m)=(xc*(2.0*l-1)*AP(J,l-1,m)-(l+m-1)*AP(J,l-2,m))/(l-m)
   enddo
enddo
 
RETURN
END SUBROUTINE
