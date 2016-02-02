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

SUBROUTINE MEMBRANE
!Paper => B. Roux, Biophys. J., 73:2980-2989 (1997)
use grandmod
use constamod
use stdiomod
use nucleotmod
use errormod
!LOCAL VARIABLES
implicit none
real sw1, dsw1, sw2, dsw2
real r, r2, charge
integer iat, itype, itype2
real  ir
logical*1 lgiat

ememb = 0.0 ! initilization

do iat = 1, ntot
  lgiat = iat.le.nsites .or. iat.ge.nsites+nfix+1
  lgiat = lgiat .and. Qforces
  if (iat.le.nsites) then  
    itype = typtyp(iat) 
    charge = 0.0 
    if (namsite(iat).eq.'P ') charge = cgnuc
  else
    itype2 = abs(typei(iat)) 
    charge = cg(itype2) 
    itype = nwtype(itype2) 
  endif

  ! Nernst transmembrane potential
  if (voltage.ne.0.0.and.charge.ne.0.0) then
    if (z(iat).lt.zmemb1) then ! REGION 1: z=z(iat)-zmemb1 < 0, lim{z->-inf} pot(1) = 0
      ememb = ememb + charge*afact*exp(ikappa*(z(iat)-zmemb1)) ! Eq. (31) paper
      if (lgiat) fz(iat) = fz(iat)- charge*afact*ikappa*exp(ikappa*(z(iat)-zmemb1))
    elseif ((z(iat).ge.zmemb1) .and. (z(iat).le.zmemb2)) then ! REGION 2
      ememb = ememb + charge*afact*(ceps*ikappa*(z(iat)-zmemb1) + 1.0) ! Eq. (31) paper
      if (lgiat) fz(iat) = fz(iat) - charge*afact*ceps*ikappa
    elseif (z(iat).gt.zmemb2) then ! REGION 3: z=z(iat)-zmemb2 > 0, lim{z->inf} pot(2) = voltage
      ememb = ememb + charge*(voltage-afact*exp(-ikappa*(z(iat)-zmemb2))) ! Eq. (31) paper
      if (lgiat) fz(iat) = fz(iat) - charge*afact*ikappa*exp(-ikappa*(z(iat)-zmemb2))
    endif
  endif ! voltage
  ! membrane boundary
  if (ampl1(itype).gt.0.0) then
    call switch1(sw1,dsw1,z(iat),p1(1,itype),p1(2,itype),zmemb1,zmemb2)
    if(sw1.ne.0.0) then ! if particle is not in the bulk
      if (Qpore) then ! cylindrical channel
        r2 = x(iat)**2+y(iat)**2
        call switch2(sw2,dsw2,r2,rcylinder(itype),p2(1,itype),p2(2,itype),r)
        if (sw2.ne.0.0) then ! if particle is not inside the pore
          if (sw1.eq.1.0) then ! if particle is inside membrane or in pore wall
            if (sw2.eq.1.0) then ! particle in membrane
              warn(itype)=warn(itype)+1
              if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine membrane :: particle inside membrane or protein - ',iat,'  '//atnam2(itype),x(iat),y(iat),z(iat)
              ememb=ememb+ampl2(itype)
            else        ! particle in pore wall
              if (lgiat) then
                ir=1.0/r
                fx(iat) = fx(iat) + ampl2(itype)*dsw2*x(iat)*ir
                fy(iat) = fy(iat) + ampl2(itype)*dsw2*y(iat)*ir
              endif
              ememb=ememb+ampl2(itype)*sw2
            endif
          else ! if particle is in membrane wall or in membrane+pore wall
            if (sw2.eq.1.0) then ! particle in membrane wall
              if (lgiat) fz(iat) = fz(iat) + ampl1(itype)*dsw1
              ememb=ememb+ampl1(itype)*sw1
            else ! particle in membrane and pore wall
              if (lgiat) then
                ir=1.0/r
                fx(iat) = fx(iat) + 0.5*ampl2(itype)*dsw2*x(iat)*ir
                fy(iat) = fy(iat) + 0.5*ampl2(itype)*dsw2*y(iat)*ir
                fz(iat) = fz(iat) + 0.5*ampl1(itype)*dsw1
              endif
              ememb = ememb + 0.5*(ampl2(itype)*sw2+ampl1(itype)*sw1)
            endif
          endif
        endif
      else ! no cylindrical channel
        if (sw1.eq.1.0) then
          warn(itype)=warn(itype)+1
          if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine membrane :: particle inside membrane or protein - ',iat,'  '//atnam2(itype),x(iat),y(iat),z(iat)
        endif
        ememb = ememb + ampl1(itype)*sw1
      endif
    endif
  endif
enddo 

return
end
