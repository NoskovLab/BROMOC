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

subroutine countions(zold,zi,i,nf,nb)
use grandmod
implicit none
integer i,j,nf(1:ntype-nold,cntpts),nb(1:ntype-nold,cntpts)
real zold,zi

do j=1,cntpts
  if ((zold.lt.zcont(j)).and.(zi.ge.zcont(j))) then
    nf(i,j) = nf(i,j) + 1
  elseif((zold.ge.zcont(j)).and.(zi.lt.zcont(j)))then
    nb(i,j) = nb(i,j) + 1
  endif
enddo
end subroutine

subroutine ioncountout(icyst,pncross)
use grandmod
use nucleotmod
use stdiomod
use constamod
implicit none
integer*8 icyst
integer i,j,nt,ncount(ntype-nold),ncross(ntype-nold,cntpts),pncross(ntype-nold,cntpts) 
real curr(ntype-nold,cntpts),currave(ntype-nold,cntpts),cnst,cnst2
character*8192 ln

nt=ntype-nold
cnst=coulomb*ipico**2/(dt*float(icyst))
cnst2=coulomb*ipico**2/(dt*float(svcntfq))

ncount = 0
do i = nsites+nfix+1, ntot
  j = abs(typei(i))-nold
  ncount(j) = ncount(j) + 1
enddo
ncross=nforward-nbackward
currave=float(ncross)*cnst  ! in pico amperes
curr=float(ncross-pncross)*cnst2 ! in pico amperes
do j=1,nt 
  write(ln,*) atnam(j+nold),icyst,float(icyst)*dt,ncount(j),(zcont(i),nforward(j,i),nbackward(j,i),curr(j,i)*cg(j+nold),currave(j,i)*cg(j+nold),i=1,cntpts)
  write(iuncnt,'(a)') trim(ln)
enddo
write(ln,*) 'TOT ',icyst,float(icyst)*dt,(zcont(i),sum(curr(1:nt,i)*cg(nold+1:ntype)),sum(currave(1:nt,i)*cg(nold+1:ntype)),i=1,cntpts)
write(iuncnt,'(a)') trim(ln)
pncross=ncross
end subroutine

