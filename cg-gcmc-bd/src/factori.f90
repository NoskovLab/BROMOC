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

REAL FUNCTION FACTORI(N)
!------------------------------------------------------------------------
!N! calculation

implicit none
INTEGER   N,I

IF(N.EQ.0.or.N.EQ.1) THEN
   FACTORI=1.0
ELSE
   FACTORI=1.0
   DO I=2,N
      FACTORI=FACTORI*I
   ENDDO
ENDIF
 
RETURN
END FUNCTION
