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

subroutine vdwgd1spln
! -----------------------------------------------------------------------
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
integer ncyz,ncel3,i,ix,iy,iz,ifir,numb,itype2
integer k,l,m,ipx,ipy,ipz
real  vdwfx,vdwfy,vdwfz,esvdw
real  xi,yi,zi,ai,bi,ci,fi,dai,dbi,dci,prefac
real  xc,yc,zc,m3,dm3
real  phisum,phivi,phivsv
logical*1 ok

ncyz = ncly2*nclz2
ncel3 = nclx2*ncyz
evdwgd = 0.0
ifir = 0
esvdw = svdw

!Main loop by atoms

do i = 1, ntot
  if (i.le.nsites .or. i.gt.(nsites+nfix)) then
    ok=x(i).le.xbcen2+tranx2-dcel2.and.x(i).ge.xbcen2-tranx2+dcel2.and. &
       y(i).le.ybcen2+trany2-dcel2.and.y(i).ge.ybcen2-trany2+dcel2.and. &
       z(i).le.vzmax.and.z(i).ge.vzmin
    if (ok) then
 !     ion cartesian coordinates in the local grid system
      xi = x(i) + tranx2-xbcen2
      yi = y(i) + trany2-ybcen2
      zi = z(i) + tranz2-zbcen2
!     if (xi.ge.0.0.and.xi.le.2.0*tranx2 .and.yi.ge.0.0.and.yi.le.2.0*trany2 .and.zi.ge.0.0.and.zi.le.2.0*tranz2) then
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
      ix = nint(xi*idcel2)
      iy = nint(yi*idcel2)
      iz = nint(zi*idcel2)

      if(ix.eq.0) ix=1
      if(iy.eq.0) iy=1
      if(iz.eq.0) iz=1
      if(ix.eq.nclx2-1)ix=nclx2-2
      if(iy.eq.ncly2-1)iy=ncly2-2
      if(iz.eq.nclz2-1)iz=nclz2-2
      
      !Atom charge distribution by 27 adjacent grid points
    
      phisum = 0.0
      do k = ix-1, ix+1
        ipx = k*ncyz
        xc = k*dcel2
        ai = 1.5 -(xi-xc)*idcel2
        dai = dm3(ai)
        ai = m3(ai)
        if (ai.ne.0.0.or.dai.ne.0.0) then
          do l = iy-1, iy+1
            ipy= l*nclz2 + ipx
            yc = l*dcel2
            bi = 1.5 -(yi-yc)*idcel2
            dbi = dm3(bi)
            bi = m3(bi)
            if (bi.ne.0.0.or.dbi.ne.0.0) then
              do m = iz-1, iz+1
                ipz = m + ipy + 1
                zc = m*dcel2
                ci = 1.5 -(zi-zc)*idcel2
                dci = dm3(ci)
                ci = m3(ci)
                if (ci.ne.0.0.or.dci.ne.0.0) then
                  fi = ai*bi*ci
                  phivi=phiv(ipz+ifir)
                  phivsv=phivi*esvdw
                  phisum = phisum + phivi
                  !Repulsive Forces
                  if (Qforces) then
                    prefac = phivsv*idcel2
                    vdwfx = vdwfx + dai*bi*ci*prefac
                    vdwfy = vdwfy + ai*dbi*ci*prefac
                    vdwfz = vdwfz + ai*bi*dci*prefac
                  endif
                  !Repulsive Energy 
                  evdwgd = evdwgd + fi*phivsv
                endif
              enddo
            endif
          enddo
        endif
      enddo
      if (phisum.ge.thold27) then
        if (.not.Qsvdw) then
          if (i.le.nsites) then
            itype2=typtyp(i)
          else
            itype2=nwtype(abs(typei(i)))
          endif
        endif
        warn(itype2)=warn(itype2)+1
        if (Qwarn) write(outu,'(a,i5,a,5f10.5)') 'Warning in routine vdwgd1spln :: particle inside membrane or protein - ',i,'  '//atnam2(itype2),x(i),y(i),z(i),phisum,thold27
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
