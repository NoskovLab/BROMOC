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

SUBROUTINE INTERACT(dener,xj,yj,zj,jtype,j,Qalert)
! calculate the interaction of particle "j" with the rest of the system
use apfmod
use ioxmod 
use grandmod
use gsbpmod
use constamod
use nucleotmod
use stdiomod
use errormod
use extramod
implicit none

!real*16 dener
real dener
real xj, yj, zj
integer i, j, itype, jtype, itype2, jtype2, is
real dist, dist2, dist6, idist, idist2
real sw1, dsw1, sw2, dsw2
real r,r2
real charge,eefp
logical*1 Qalert, Beion
integer isite1, isite2, isite3, isite4
real  dd, vard
real  ang, ang0, varang, vecR(3,3)
real  dihe, dihe0, vardihe
real  cte1, epsln, sgex2
real  angbond, angdih    
integer nindex

dener = 0.0 ! Initialization

Beion = Qpar .and. j.gt.nsites
if (Beion) then
  jtype2 = nwtype(jtype)
else 
  jtype2 = typtyp(j)
endif

if (Qenergy) then ! total energy
  if (Qmemb) then 
! Nernst transmembrane potential
    ememb  = 0.0 ! Initialization
    if (voltage.ne.0.0) then ! voltage -> transmembrane potential
      if (Beion) then
        charge = cg(jtype)
      else if (Qnucl) then
        charge = 0.0
        if (namsite(j).eq.'P ') charge = cgnuc
      endif 
      if (charge.ne.0.0) then
        if (zj.lt.zmemb1) then ! REGION 1: z=zj-zmemb1 < 0, lim{z->-inf} pot(1) = 0 
          ememb = ememb + charge*afact*exp(ikappa*(zj-zmemb1))
        else if ((zj.ge.zmemb1).and.(zj.le.zmemb2)) then ! REGION 2
          ememb = ememb + charge*afact*(ceps*ikappa*(zj-zmemb1)+1.0)
        else if (zj.gt.zmemb2) then ! REGION 3: z=zj-zmemb2 > 0, lim{z->inf} pot(2) = voltage
          ememb = ememb + charge*(voltage-afact*exp(-ikappa*(zj-zmemb2)))
        endif
      endif
    endif ! voltage
    ! Membrane boundary
    if (ampl1(jtype2).gt.0.0) then ! if membrane is not permeable to particle jtype2
      call switch1(sw1,dsw1,zj,p1(1,jtype2),p1(2,jtype2),zmemb1,zmemb2)
      if (sw1.ne.0.0) then ! particle is not in the bulk
        if (Qpore) then ! Cylindrical channel
          r2 = xj**2+yj**2
          call switch2(sw2,dsw2,r2,rcylinder(jtype2),p2(1,jtype2),p2(2,jtype2),r)
          if (sw2.ne.0.0) then ! particle is not inside pore
            if (sw1.eq.1.0) then ! particle is inside membrane or in pore wall
              if (sw2.eq.1.0) then ! particle is in membrane
                if (Qalert) then
                  warn(jtype2)=warn(jtype2)+1
                  if (Qwarn) write(outu,'(a,i5,a,5f10.5)')'Warning in routine interact :: particle inside membrane - ',j,'  '//atnam2(jtype2),xj,yj,zj
                endif
                ememb=ememb+1.0e10
              else  ! particle in pore wall
                ememb=ememb+ampl2(jtype2)*sw2
              endif
            else ! particle is in membrane wall or in membrane+pore wall
              if (sw2.eq.1.0) then ! particle in membrane wall
                ememb=ememb+ampl1(jtype2)*sw1
              else  ! particle in membrane and pore wall
                ememb = ememb + 0.5*(ampl2(jtype2)*sw2+ampl1(jtype2)*sw1)
              endif
            endif
          endif
        else ! no cylindrical channel
          if (sw1.eq.1.0) then
            if (Qalert) then
              warn(jtype2)=warn(jtype2)+1
              if (Qwarn) write(outu,'(a,i5,a,5f10.5)')'Warning in routine interact :: particle inside membrane - ',j,'  '//atnam2(jtype2),xj,yj,zj
            endif
            ememb=ememb+1.0e10
          else
            ememb = ememb + ampl1(jtype2)*sw1
          endif
        endif
      endif
    endif
    dener = ememb
  endif !Qmemb

  !reaction field contribution
  if (Qmmij) then
    if (shapes.EQ.'RECTBOX ') then
      if (j.eq.ntot+1) then ! in case of insert-try
        ntot = ntot+1
        call rect_rf0(xj,yj,zj,j,jtype)
        ntot = ntot-1
      else
       call rect_rf0(xj,yj,zj,j,jtype)
      endif
    elseif (shapes.eq.'SPHERE  ') then
      if (j.eq.ntot+1) then ! in case of insert-try
        ntot = ntot+1
        call sphe_rf0(xj,yj,zj,j,jtype)
        ntot = ntot-1
      else             
        call sphe_rf0(xj,yj,zj,j,jtype)
      endif  
    endif
    dener = dener + egsbpb
  endif ! Qmmij

  !reaction field parameter  
  if (Qrfpar) then
    if (j.eq.ntot+1) then ! insert
      ntot=ntot+1
      call rfparionj(xj,yj,zj,j,jtype)
      ntot=ntot-1
    else  ! remove
      call rfparionj(xj,yj,zj,j,jtype)
    endif
    dener=dener+erfpar
    egsbpb=erfpar
  endif ! Qrfpar

  !static external field contribution
  if (Qphix) then
    call staticf0(xj,yj,zj,j,jtype)
    dener = dener + egsbpa
  endif ! Qphix

  !grid-based repulsive potential
  if (Qphiv) then
    if (Qtrln) then
       call vdwgd0trln(xj,yj,zj,j,jtype,jtype2,Qalert)
    else
       call vdwgd0spln(xj,yj,zj,j,jtype,jtype2,Qalert)
    endif
    dener = dener + evdwgd
  endif ! Qphiv
  
  !nonbonded interaction between ions
  if (Beion) then
    if (Qnonbond) then
      do i = nsites+1, ntot
        if (i.ne.j) then
          itype = abs(typei(i))
          itype2=nwtype(itype)
          is=nindex(itype2,jtype2)
          dist2 = (x(i)-xj)**2+(y(i)-yj)**2+(z(i)-zj)**2
          if (Qefpot(is)) then
            call gety(is,dist2,eefp,dist)
            dener=dener+eefp
          else
            idist2 = 1.0/dist2
            if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
    !Lennard-Jones 6-12 potential + electrostatic interaction    
            if (Qchr(is)) then
              dener = dener + fct(is)*idist
            endif
            if (Qlj(is)) then
              dist6 = (sgp2(is)*idist2)**3
              dener  = dener + epp4(is)*dist6*(dist6-1.0)
            endif
     !water-mediated short-range ion-ion interaction  
            if (Qsrpmfi(is)) then
              if (dist2.le.rth) then
                dist=1.0/idist
                esrpmf = c0(is)*exp((c1(is)-dist)*c2(is))*cos(c3(is)*pi*(c1(is)-dist))+c4(is)*(c1(is)*idist)**6
                if (dist.ge.srpx) esrpmf=esrpmf*exp(-srpk*(dist-srpx))-srpy  ! smoothly fix discontinuity 
                ! Eq. 9 W. Im,and B. Roux J. Mol. Biol. 322:851-869
                ! (2002)
                dener = dener + esrpmf
              endif
            endif  
          endif 
        endif ! j  
      enddo ! i=nsites+1,...,ntot  
    endif !Qnonbond
  endif ! Beion 
 !interaction between interaction sites (nucleotides)
  if (Qnucl .and. .not.Beion) then
    if (Qnobond) then
   !  Bonded interactions
   !  -------------------
   !  The strech energy is calculated by means of a
   !  Taylor expansion around natural bond lenght.
   !  This contribution contains harmonic and anharmonic
   !  interactions.
      do i = 1, nbond
        isite1 = sitebond(i,1)
        isite2 = sitebond(i,2)
        if (isite1.eq.j .or. isite2.eq.j) then
          dd = distbond(i) ! natural bond length
          dist = sqrt((x(isite1)-x(isite2))**2 + (y(isite1)-y(isite2))**2 + (z(isite1)-z(isite2))**2)
          vard = dist - dd
          dener = dener + epsnuc*vard**2*(1.0+100.0*vard**2) 
        endif        
      enddo
   !  The bending energy is calculated
   !  by means of a Taylor expansion around natural bond
   !  angle. This contribution contains harmonic interactions.
      do i = 1, nangle
        isite1 = siteangle(i,1)
        isite2 = siteangle(i,2) ! central site
        isite3 = siteangle(i,3)
        if (isite1.eq.j .or. isite2.eq.j .or. isite3.eq.j) then
          ang0 = valangle(i) ! natural bond angle
          vecR(1,1) = x(isite1) - x(isite2)
          vecR(2,1) = y(isite1) - y(isite2)
          vecR(3,1) = z(isite1) - z(isite2)
          vecR(1,2) = x(isite3) - x(isite2)
          vecR(2,2) = y(isite3) - y(isite2)
          vecR(3,2) = z(isite3) - z(isite2)
          ang = angbond(vecR)
          varang = ang - ang0
          dener = dener + 700.0*epsnuc*varang**2
        endif  
      enddo  
   !  The torsional energy is calculated
   !  using Fourier series for natural torsional angle.
      do i = 1, ndihe
        isite1 = sitedihe(i,1)
        isite2 = sitedihe(i,2)
        isite3 = sitedihe(i,3)
        isite4 = sitedihe(i,4)
        if (isite1.eq.j .or. isite2.eq.j .or. isite3.eq.j .or. isite4.eq.j) then
          dihe0 = valdihe(i) ! natural torsion angle
          vecR(1,1) = x(isite2) - x(isite1)
          vecR(2,1) = y(isite2) - y(isite1)
          vecR(3,1) = z(isite2) - z(isite1)
          vecR(1,2) = x(isite3) - x(isite2)
          vecR(2,2) = y(isite3) - y(isite2)
          vecR(3,2) = z(isite3) - z(isite2)
          vecR(1,3) = x(isite4) - x(isite3)
          vecR(2,3) = y(isite4) - y(isite3)
          vecR(3,3) = z(isite4) - z(isite3)
          dihe = angdih(vecR)
          vardihe = dihe - dihe0
          dener = dener + 28.0*epsnuc*(1.0-cos(vardihe))
        endif 
      enddo
    endif ! Qnobond
    if (Qnonbond) then   
   !   Non-bonded interactions
   !   -----------------------
   !   Native contacts
      do i = 1, nstack
        isite1 = sitestack(i,1)
        isite2 = sitestack(i,2)
        if (isite1.eq.j .or. isite2.eq.j) then
          dist2 = (x(isite2)-x(isite1))**2 + (y(isite2)-y(isite1))**2 + (z(isite2)-z(isite1))**2
          cte1 = (sgstack(i)*sgstack(i)/dist2)**3
          dener = dener + 4.0*epsnuc*cte1*(cte1-1.0)
        endif 
      enddo  
   !  Hydrogen bonding
      do i = 1, nbp
        isite1 = sitebp(i,1)
        isite2 = sitebp(i,2)
        if (isite1.eq.j .or. isite2.eq.j) then
          if (namsite(isite1).eq.'Gb' .or. namsite(isite1).eq.'Cb') then
            epsln = 2.532*epsnuc
          else
            epsln = 2.0*epsnuc
          endif
          dist2 = (x(isite2)-x(isite1))**2 + (y(isite2)-y(isite1))**2 + (z(isite2)-z(isite1))**2
          cte1 = sgbp(i)*sgbp(i)/dist2
          dener = dener + epsln*cte1**5*(20.0*cte1-24.0)
        endif
      enddo   
   !  Excluded volume
      do i = 1, nex
        isite1 = siteex(i,1)
        isite2 = siteex(i,2)
        if (isite1.eq.j .or. isite2.eq.j) then
          dist2 = (x(isite2)-x(isite1))**2 + (y(isite2)-y(isite1))**2 + (z(isite2)-z(isite1))**2
          sgex2=sgex(i)**2
          if (dist2.lt.sgex2) then
            cte1 = (sgex2/dist2)**3
            dener = dener + 4.0*epsnuc*cte1*(cte1-1.0) + epsnuc
          endif  
        endif 
      enddo  
   !  Coulomb interaction
      do i = 1, nqq
        isite1 = siteqq(i,1)
        isite2 = siteqq(i,2)
        if (isite1.eq.j .or. isite2.eq.j) then
          dist2 = (x(isite2)-x(isite1))**2 + (y(isite2)-y(isite1))**2 + (z(isite2)-z(isite1))**2
          dist = sqrt(dist2)
          idist=1.0/dist
          if (Qdebyhyb) then
            if(outbox(x(isite1),y(isite1),z(isite1)).and.outbox(x(isite2),y(isite2),z(isite2))) Qdeby=.true.
          endif
          if (.not.Qdeby) then
            dener = dener + fctn*idist
          else  
            dener = dener + fctn*idist*exp(-dist*ikappa) 
          endif
          if (Qdebyhyb) Qdeby=.false.  
        endif
      enddo 
   ! Solvent-induced contribution
      if (Qsolv) then
        do i = 1, nsolv
          isite1 = siteslv(i,1)
          isite2 = siteslv(i,2)
          if (isite1.eq.j .or. isite2.eq.j) then
            dist2 = (x(isite2)-x(isite1))**2 + (y(isite2)-y(isite1))**2 + (z(isite2)-z(isite1))**2
            dist = sqrt(dist2)
            cte1 = exp((13.38-dist)*0.1875)
            dener = dener + epsolv*(1.0-cte1)**2 - epsolv
          endif
        enddo
      endif 
    endif ! Qnonbond 
  endif ! Qnucl .and. .not.Beion  
  ! nonbonded interactions between interaction sites and ions       
  if (Qpar .and. Qnucl .and. Qnonbond) then
    if (Beion) then
      do i = 1, nsites
        itype2 = typtyp(i)
        is=nindex(itype2,jtype2)
        dist2 = (x(i)-xj)**2 + (y(i)-yj)**2 + (z(i)-zj)**2    
        if (Qefpot(is)) then
          call gety(is,dist2,eefp,dist)
          dener=dener+eefp
        else
          idist2 = 1.0/dist2
          if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
          if (Qchr(is)) then 
            dener = dener + fct(is)*idist
          endif
          if (Qlj(is)) then
            dist6 = (sgp2(is)*idist2)**3
            dener = dener + epp4(is)*dist6*(dist6-1.0)
          endif
   ! Compute Short-Range correction contribution
          if (Qsrpmfi(is)) then
            if (dist2.le.rth) then
              dist=1.0/idist
              esrpmfmx = c0(is)*exp((c1(is)-dist)*c2(is))*cos(c3(is)*pi*(c1(is)-dist))+c4(is)*(c1(is)*idist)**6
              if (dist.ge.srpx) esrpmfmx=esrpmfmx*exp(-srpk*(dist-srpx))-srpy  ! smoothly fix discontinuity 
                ! Eq. 9 W. Im,and B. Roux J. Mol. Biol. 322:851-869
                ! (2002)
              dener = dener + esrpmfmx
            endif
          endif
        endif
      enddo ! i=1,...,nsites
    else
      do i = nsites+1, ntot
        itype = abs(typei(i))
        itype2 = nwtype(itype)
        is=nindex(itype2,jtype2)
        dist2 = (x(i)-xj)**2 + (y(i)-yj)**2 + (z(i)-zj)**2
        if (Qefpot(is)) then
          call gety(is,dist2,eefp,dist)
          dener=dener+eefp
        else
          idist2 = 1.0/dist2
          if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
          if (Qchr(is)) then
            dener = dener + fct(is)*idist
          endif
          if (Qlj(is)) then
            dist6 = (sgp2(is)*idist2)**3
            dener = dener + epp4(is)*dist6*(dist6-1.0)
          endif
 ! Compute Short-Range correction contribution
          if (Qsrpmfi(is)) then 
            if (dist2.le.rth) then
              dist=1.0/idist
              esrpmfmx = c0(is)*exp((c1(is)-dist)*c2(is))*cos(c3(is)*pi*(c1(is)-dist))+c4(is)*(c1(is)*idist)**6
              if (dist.ge.srpx) esrpmfmx=esrpmfmx*exp(-srpk*(dist-srpx))-srpy  ! smoothly fix discontinuity 
                ! Eq. 9 W. Im,and B. Roux J. Mol. Biol. 322:851-869
                ! (2002)
              dener = dener + esrpmfmx
            endif
          endif
        endif
      enddo ! i=nsites+1,...,ntot
    endif ! Beion
  endif ! Qpar .and. Qnucl .and. Qnonbond
endif !Qenergy

return
end
