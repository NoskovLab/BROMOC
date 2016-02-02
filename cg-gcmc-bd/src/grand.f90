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

SUBROUTINE GRAND(ngcmc,prob,icycle)
!THE GRAND CANONICAL MONTECARLO ALGORITHIM      
use grandmod
use constamod
use stdiomod
use nucleotmod
use errormod
implicit none

! NOTE: RANDOM IS A REA FUNCTION

integer ngcmc,prob(nbuffer,0:datom)
integer*8 icycle 
!LOCAL VARIABLES
integer igrand, ib, ip, iat, itype
!real*16 dener
real dener
real  rate, xnew, ynew, znew, nbufferr
nbufferr=float(nbuffer)

!Loop over number of steps for GCMC, ngcmc      
do igrand = 1, ngcmc

  if (rndm().lt.0.5) then 
    !Try to insert a particle into one of the buffers
    ib = int(nbufferr*rndm()) + 1 ! [1,nbuferr+1]
    if (ib.le.nbuffer .and. avnum(ib).gt.0.0) then
      itype = ibfftyp(ib) ! ion type
      call insert(ib,xnew,ynew,znew)
      call interact(dener,xnew,ynew,znew,itype,ntot+1,.false.)

      if (Qbufferbias(ib)) then
        rate = (avnum(ib)+kb(ib)*(avnum(ib)-real(ntotat(ib)/icycle)))
        rate = rate*exp(-(dener-mu(ib))*ikBT)/float(nat(ib)+1)
      else
        rate = (avnum(ib)/(nat(ib)+1))*exp(-(dener-mu(ib))*ikBT)
        ! Eq. 17 W. Im, and B. Roux Biophys. J. 79:188-801 (2000)          
      endif
      rate = rate/(1.0+rate) ! creation transition probability
      ! Eq. 5 W. Im, and B. Roux Biophys. J. 79:188-801 (2000) 
      if (rndm().le.rate) then
        ! The creation is accepted if a random number 'random(iseed)'
        ! between 0 and 1 is less than or equal to creation transition
        ! probability                   
        ener = ener + dener
        nat(ib) = nat(ib) + 1
        ntot = ntot + 1
        if (ntot.gt.datom) then
          call error ('grand', 'ntot is greater than datom',faterr)
        endif
        natom = natom + 1
        x(ntot) = xnew
        y(ntot) = ynew
        z(ntot) = znew
        typei(ntot) = itype
        if (znew.lt.cz) typei(ntot) = -itype
        ibuffer(ntot) = ib
        ninsert(ib) = ninsert(ib) + 1
      endif
    endif 
  else
  !  Try to remove a particle
    ib    = int(nbufferr*rndm())+1 ! number of buffer [1,nbuffer+1]
    ip    = int(float(nat(ib))*rndm())+1 ! number of particle [1,nat(ib)+1]
    if (ib.le.nbuffer .and. ip.le.nat(ib)) then
      itype = ibfftyp(ib) ! ion type
      call find(ibuffer,ib,nfix+nsites,ntot,ip,iat)
      call interact(dener,x(iat),y(iat),z(iat),itype,iat,.true.)

      if (Qbufferbias(ib)) then
        rate = (avnum(ib)+kb(ib)*(avnum(ib)-real(ntotat(ib)/icycle)))
        rate = rate*exp(-(dener-mu(ib))*ikBT)/float(nat(ib))
      else
        rate = (avnum(ib)/nat(ib))*exp(-(dener-mu(ib))*ikBT)
        ! Eq. 17 W. Im, and B. Roux Biophys. J. 79:188-801 (2000)
      endif
      rate = 1.0/(1.0+rate) ! destruction transition probability
      ! Eq. 6 W. Im, and B. Roux Biophys. J. 79:188-801 (2000)
      if (rndm().le.rate) then
       ! The destruction is accepted if a random number between 0 and 1 
       ! is less than or equal to destruction transition probability
        ener = ener - dener
        nremove(ib) = nremove(ib) + 1
        x(iat) = x(ntot)  ! copy coordinates of last particle
        y(iat) = y(ntot)
        z(iat) = z(ntot)
        typei(iat) = typei(ntot)
        ibuffer(iat) = ibuffer(ntot)
        nat(ib) = nat(ib) - 1
        natom = natom - 1
        ntot = ntot - 1
      endif
    endif  
  endif

  do ib = 1, nbuffer
    prob(ib,nat(ib)) = prob(ib,nat(ib)) + 1
  enddo

enddo

return
end
