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

subroutine repwalls(walls)
!-----------------------------------------------------------------------
!This subroutine creates sucesive repulsion walls for repulsion map. This will allow particles to avoid entering into closed repulsion bodies (protein).

use gsbpmod
use nucleotmod

implicit none
integer*1,allocatable :: phivt(:)
integer*1 n,walls,one
integer a,b,c,i,j,k,pos,phisum,ck,ifir,ncel3,m,oo
one=1

ncel3=nclx2*ncly2*nclz2
oo=1

if (Qnmcden) oo=nttyp

allocate (phivt(ncel3*oo))
do m=1,oo
  ifir=(m-1)*ncel3
  do n=one,walls
    phivt=phiv
    ck=0
    do a = 1+n,nclx2-n
      do b = 1+n,ncly2-n
        do c = 1+n,nclz2-n
          pos=(a-1)*ncly2*nclz2 + (b-1)*nclz2 + c
          if (phiv(pos+ifir).eq.n) then 
            ck=ck+1
            phisum=0
            do i = a-1, a+1
              do j = b-1, b+1
                do k = c-1, c+1
                  phisum = phisum + phiv((i-1)*ncly2*nclz2 + (j-1)*nclz2 + k + ifir)
                enddo
              enddo
            enddo
            if (phisum.eq.27*n) phivt(pos)=n+one
          endif
        enddo
      enddo
    enddo
    phiv=phivt
    if (ck.eq.0) exit
  enddo
enddo
deallocate (phivt)
end subroutine
