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

SUBROUTINE SINMPHI2(J,M,CC0,S0,AS)
!------------------------------------------------------------------------
!sin(M*phi) calculation (M > 0)
use grandmod
implicit none

INTEGER   M,I,J
REAL    CC0,S0,AS(NTOT,0:M)

AS(J,0) = 0.0
AS(J,1) = S0
AS(J,2) = 2.0*CC0*S0
AS(J,3) = 2.0*CC0*AS(J,2)-S0
DO I = 4, M
  AS(J,I) = 2.0*CC0*AS(J,I-1)-AS(J,I-2)
ENDDO
 
RETURN
END
