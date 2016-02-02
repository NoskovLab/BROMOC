!    POTMAP-PBEQ-READ   Reads CHARMM PBEQ Potential Map
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

program potmap_pbeq_read
implicit none
INTEGER*4 IUNIT,NCLX,NCLY,NCLZ,ITYPE
REAL*4,allocatable ::  SMTHNG(:)
INTEGER*4 I,NC3,IFIR,ILAS,ini,fin
REAL*8  DCEL,XBCEN,YBCEN,ZBCEN
REAL*8  EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
character filename*256
write(*,'(A$)') 'Input CHARMM PBEQ Potential Map filename: '
read(*,*) filename

iunit=1
open(unit=iunit,file=filename,form='unformatted')
rewind(unit=iunit)

IF(IUNIT.GT.0)THEN
   READ(IUNIT) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
   READ(IUNIT) EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
   NC3=NCLx*NCLy*NCLz
   !IFIR=(itype-1)*NC3+1
   IFIR=1
   !ILAS=(itype-1)*NC3+NC3
   ILAS=NC3
   allocate (smthng(nc3))
   READ(IUNIT) (SMTHNG(I),I=IFIR,ILAS)
ENDIF

write(*,*) 'Range from ',ifir,' to ',ilas
write(*,'(A$)') 'Choose range to print: '
read(*,*) ini,fin

WRITE(*,*) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
WRITE(*,*) EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
write(*,*) (smthng(i),i=ini,fin)

close(iunit)
end program
