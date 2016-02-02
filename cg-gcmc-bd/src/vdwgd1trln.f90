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

subroutine vdwgd1trln
!-----------------------------------------------------------------------
!This subroutine computes the repulsive potental energy and forces
!based on the grid

!Repulsive forces - VDWF

use ioxmod
use constamod
use stdiomod
use grandmod
use nucleotmod
use gsbpmod     
!local variables
implicit none
integer ncyz,ncel3,i,ix,iy,iz,n1,n2,n3,in3,ifir,numb,itype2
real  vdwfx,vdwfy,vdwfz
real  xi,yi,zi,ai,bi,ci,fi
real  aisign,bisign,cisign,prefac
real  phisum,one,phis,esvdw
logical*1 ok
one=1.0

ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz
evdwgd = 0.0
ifir = 0
esvdw = svdw

!Main loop by atoms

do i = 1, ntot
  if (i.le.nsites .or. i.gt.(nsites+nfix)) then
    ok=x(i).le.xbcen2+tranx2.and.x(i).ge.xbcen2-tranx2.and. &
       y(i).le.ybcen2+trany2.and.y(i).ge.ybcen2-trany2.and. &
       z(i).le.vzmax.and.z(i).ge.vzmin
    if (ok) then
      !ion cartesian coordinates in the local grid system                 
      xi = x(i) + tranx2-xbcen2
      yi = y(i) + trany2-ybcen2
      zi = z(i) + tranz2-zbcen2
!      if (xi.ge.0.0.and.xi.le.2.0*tranx2 .and.yi.ge.0.0.and.yi.le.2.0*trany2 .and.zi.ge.0.0.and.zi.le.2.0*tranz2) then
      if (Qnmcden) then                                      
        if (Qnucl .and. i.le.nsites) then
          if (namsite(i).eq.'S ') then
            ifir = 0
          else if (namsite(i).eq.'P ') then
            ifir = ncel3 
          else if (namsite(i).eq.'Ab') then
            ifir = 2*ncel3 
          else if (namsite(i).eq.'Tb') then
            ifir = 3*ncel3 
          else if (namsite(i).eq.'Cb') then
            ifir = 4*ncel3 
          else ! namsite(i).eq.'Gb'
            ifir = 5*ncel3 
          endif
        else if (Qnucl .and. Qpar .and. i.gt.(nsites+nfix)) then
          if (istrs.eq.1) numb = abs(typei(i)) - (inuc+1)
          if (istrs.eq.2) numb = abs(typei(i)) - (2*inuc+1)
          ifir = (6 + numb)*ncel3 
        else if (.not.Qnucl .and. Qpar .and. i.gt.nfix) then
          ifir = (abs(typei(i))-1)*ncel3
        endif
      else
        if (Qsvdw) then 
          if (i.le.nsites) then 
            itype2=typtyp(i)
          else
            itype2=nwtype(abs(typei(i)))
          endif
          esvdw = svdw * scal(itype2)
        endif
      endif                                                     
      !initializations forces     
      vdwfx = 0.0
      vdwfy = 0.0
      vdwfz = 0.0
      !integer counter for ion cartesian coordinates             
      ix = int(xi*idcel2)
      iy = int(yi*idcel2)
      iz = int(zi*idcel2)
      if (ix.eq.nclx2-1) ix=nclx2-2
      if (iy.eq.ncly2-1) iy=ncly2-2
      if (iz.eq.nclz2-1) iz=nclz2-2

!Atom charge distribution by 8 adjacent grid points

      phisum = 0.0
      do n1 = ix, ix+1
        ai = xi - n1*dcel2
        aisign = sign(one,ai)
        ai = 1.0 - abs(ai)*idcel2
        do n2 = iy, iy+1
          bi = yi - n2*dcel2
          bisign = sign(one,bi)
          bi =1.0 - abs(bi)*idcel2
          do n3 = iz, iz+1
            ci = zi - n3*dcel2
            cisign = sign(one,ci)
            ci = 1.0 - abs(ci)*idcel2
            fi = ai*bi*ci ! fraction of the charge assigned to a
                          ! grid point
            in3 = n1*ncyz + n2*nclz2 + n3 + 1
            phis=phiv(in3+ifir)
            phisum = phisum + phis
 !Electrostatic forces
            if (Qforces) then
              prefac = phis*esvdw*idcel2
              if ((ai.lt.(1.0-rsmall)).and.(ai.gt.rsmall)) then
                vdwfx = vdwfx + aisign*bi*ci*prefac
              endif
              if ((bi.lt.(1.0-rsmall)).and.(bi.gt.rsmall)) then
                vdwfy = vdwfy + bisign*ai*ci*prefac
              endif
              if ((ci.lt.(1.0-rsmall)).and.(ci.gt.rsmall)) then
                vdwfz = vdwfz + cisign*ai*bi*prefac
              endif
            endif
! Electrostatic Energy 
            evdwgd = evdwgd + fi*esvdw*phis
          enddo ! n3
        enddo ! n2
      enddo ! n1
      if (phisum.ge.thold8) then
        if (.not.Qsvdw) then
          if (i.le.nsites) then
            itype2=typtyp(i)
          else
            itype2=nwtype(abs(typei(i)))
          endif
        endif
        warn(itype2)=warn(itype2)+1
        if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine vdwgd1trln :: particle inside membrane or protein - ',i,'  '//atnam2(itype2),x(i),y(i),z(i),phisum,thold8
      endif
      if (Qforces) then  
        fx(i) = fx(i) + vdwfx
        fy(i) = fy(i) + vdwfy
        fz(i) = fz(i) + vdwfz
      endif
    endif  
  endif   
enddo

return
end
