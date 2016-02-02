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

SUBROUTINE DALPOL2(J,LMAX,MMAX,XC,AP,ADP)
!------------------------------------------------------------------------
!Derivatives of Associate Legendre Polynomials (From Smythe's BOOK)
!This is different from FUNCTION DALPOL because we don't consider xc > 1 here.

use grandmod
implicit none

integer lmax,mmax,j
real  xc,ap(ntot,0:lmax,0:mmax+1),adp(ntot,0:lmax,0:mmax)
!local
integer l,m
real  fact

IF(XC.EQ.1.0.OR.XC.EQ.-1.0) THEN
   DO L = 0,LMAX
      ADP(J,L,0) = 0.0
      IF(XC.EQ. 1.0) ADP(J,L,0) = L*(L+1.0)*0.5
      IF(XC.EQ.-1.0) ADP(J,L,0) = (-1.0)**(L+1)*L*(L+1.0)*0.5
      DO M = 1, MMAX
         ADP(J,L,M) = 0.0
      ENDDO
   ENDDO
   RETURN
ENDIF

DO L = 0, LMAX
   DO M = 0, MMAX
      FACT = 1.0/SQRT(1.0-XC*XC)
      ADP(J,L,M) = FACT*(AP(J,L,M+1)-FACT*M*XC*AP(J,L,M))
   ENDDO
ENDDO
 
RETURN
END SUBROUTINE
