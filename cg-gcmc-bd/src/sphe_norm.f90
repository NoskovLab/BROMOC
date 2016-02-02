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

SUBROUTINE SPHE_NORM(NMPOL,BNORM,SRDIST)
!-----------------------------------------------------------------------
!calculate the normalization constants of the spherical harmonics basis functions
!input: SRDIST -> radius sphere
!       NMPOL  -> number of multipoles
!output: BNORM -> Normalization constant

use constamod
implicit none
INTEGER NMPOL
REAL  BNORM(*),SRDIST
!local
INTEGER L,M,NORDER
REAL  SR2,LPART,RPART,UPFACTO,DNFACTO,FACTORI,ipi,itwopi

!always the same order in spherical harmonics
NORDER=1
SR2=SRDIST*SRDIST
RPART=SR2*SRDIST

itwopi=1.0/twopi
LPART=1.5*ITWOPI
ipi=0.5*ITWOPI
BNORM(NORDER)=SQRT(LPART/RPART)
DO L=1,NMPOL-1
   LPART=(2.0*L+1.0)*(2.0*L+3.0)*ipi
   RPART=RPART*SR2
   NORDER=NORDER+1
   BNORM(NORDER)=SQRT(LPART/RPART)
   LPART=(2.0*L+1.0)*(2.0*L+3.0)*ITWOPI       ! change for m > 0
   DO M=1,L
      NORDER=NORDER+1 
      UPFACTO=FACTORI(L-M)
      DNFACTO=FACTORI(L+M)
      BNORM(NORDER)=SQRT(LPART*UPFACTO/(RPART*DNFACTO))
      NORDER=NORDER+1 
      BNORM(NORDER)=BNORM(NORDER-1)
   ENDDO
ENDDO

RETURN
END
