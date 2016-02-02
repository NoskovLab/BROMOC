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

SUBROUTINE ENERGY
use apfmod
use ioxmod
use grandmod
use gsbpmod
use constamod
use nucleotmod
use extramod
use efpmod
implicit none

integer i, j, itype, itype2, jtype, jtype2, is
real dist, dist2, dist6, idist, idist2, idistkd
real de, dc

real  esrpmf0,esrpmf1,esrpmf2,esrpmf3
real  cofo
integer isite1, isite2, isite3, isite4
real  dd,vard,eefp,fdf,fdv,vard2
real  ang, ang0, varang, cte1
real  dihe, dihe0, vardihe,didstmp
real  epsln,sgex2,f1(3),f2(3),f3(3),f4(3),v1(3),v2(3),v3(3),m1,m2,m3,modval
real  dist12,esolve,n1(3),n2(3),n1v,n2v,v22,v2v,av,bv,iv22,v12,v23
logical*1 ok,ok2
integer nindex

! Initializations
ener     = 0.0
eelec    = 0.0
eefpot   = 0.0
eefpotmx = 0.0
evdw     = 0.0
esrpmf   = 0.0
esrpmfmx = 0.0
ememb    = 0.0
esolve   = 0.0
egsbpa   = 0.0
egsbpb   = 0.0
evdwgd   = 0.0
ebond    = 0.0
eang     = 0.0
edihe    = 0.0
estack   = 0.0 
ebp      = 0.0
eex      = 0.0
eqq      = 0.0
esolv    = 0.0
eqqmx    = 0.0
evdwmx   = 0.0
econ     = 0.0


if (Qenergy) then
  ! Initializations
  if (Qforces) then
     do i = 1, ntot
        fx(i) = 0.0
        fy(i) = 0.0
        fz(i) = 0.0
     enddo
  endif

  ! membrane contribution
  if (Qmemb) call membrane
  ener = ememb

  ! reaction field contribution
  if (Qmmij) then
     if(shapes.EQ.'RECTBOX ') then
        call rect_rf1
     else if(shapes.eq.'SPHERE  ') then
        call sphe_rf1
     endif
     ener = ener + egsbpb
  endif

  ! static external field contribution
  if (Qphix) then
     call staticf1
     ener = ener + egsbpa
  endif

  ! grid-based repulsive potential
  if (Qphiv) then
     if (Qtrln) then
        call vdwgd1trln
     else
        call vdwgd1spln
     endif
     ener = ener + evdwgd
  endif
  
  if (Qrfpar) then
    call rfparion
    ener = ener + erfpar
    egsbpb=erfpar
  endif

  ! nonbonded interaction between ions
  if (Qpar) then                 
    if (Qnonbond) then
      do j = nsites+nfix+1, ntot
        jtype  = abs(typei(j))
        jtype2 = nwtype(jtype) ! convert atnam to atnam2
        do i = nsites+1, j-1
          itype = abs(typei(i))
          itype2=nwtype(itype) ! convert atnam to atnam2
          is=nindex(itype2,jtype2)
          dist2 = ((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
          if (Qefpot(is)) then
            if (Qforces) then
              call getyd(is,dist2,eefp,de,dist)
            else
              call gety(is,dist2,eefp,dist)
            endif
            eefpot = eefpot + eefp
          else 
            idist2 = 1.0/dist2
            if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
            if (Qchr(is)) then
              cofo=fct(is)*idist
              eelec = eelec + cofo
            endif
            if (Qlj(is)) then
              dist6 =(sgp2(is)*idist2)**3
              dist12=dist6**2
              evdw = evdw + epp4(is)*(dist12-dist6) ! van der waals potential
            endif
            if (Qsrpmfi(is)) then 
              if (dist2.le.rth) then
                dist=1.0/idist
                esrpmf1 = exp((c1(is)-dist)*c2(is))
                esrpmf2 = cos(c3(is)*pi*(c1(is)-dist))
                esrpmf3 = (c1(is)*idist)**6
                esrpmf0 = c0(is)*esrpmf1*esrpmf2+c4(is)*esrpmf3 
                if (dist.ge.srpx) then ! smoothly fix discontinuity 
                  fdf=exp(-srpk*(dist-srpx))-srpy
                  fdv=esrpmf0 
                  esrpmf0=fdv*fdf
                endif
                esrpmf = esrpmf + esrpmf0
              endif
            endif
          endif  
          if (Qforces) then
            if (.not.Qefpot(is)) then
              de=0.0
              if (Qchr(is)) then
                de=cofo*idist2
              endif
              if (Qlj(is)) then
                de=de+epp4(is)*(2.0*dist12-dist6)*6.0*idist2 ! van der waals forces
              endif
              if (Qsrpmfi(is)) then 
                if (dist2.le.rth) then 
                  dc=(-c2(is)*esrpmf2+c3(is)*pi*sin(c3(is)*pi*(c1(is)-dist)))*c0(is)*esrpmf1-6.0*c4(is)*esrpmf3*idist! forces
                  if (dist.ge.srpx) dc=dc*fdf-fdv*srpk*(fdf+srpy)  ! smoothly fix discontinuity 
                  de=de-dc*idist
                endif
              endif
            endif
            if (de.ne.0.0) then 
              if (i.gt.(nsites+nfix)) then
                fx(i) = fx(i) + de*(x(i)-x(j))
                fy(i) = fy(i) + de*(y(i)-y(j))
                fz(i) = fz(i) + de*(z(i)-z(j))
              endif
              fx(j) = fx(j) - de*(x(i)-x(j))
              fy(j) = fy(j) - de*(y(i)-y(j))
              fz(j) = fz(j) - de*(z(i)-z(j))
            endif
          endif
        enddo ! i = nsites+1,...,(j-1)
      enddo ! j = nsites+nfix+1,...,ntot 
    endif !Qnonbond
  endif  ! Qpar

  ! interaction between interaction sites (nucleotides)         
  if (Qnucl) then
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
        dd = distbond(i) ! natural bond length
        v1(1)=x(isite2)-x(isite1)
        v1(2)=y(isite2)-y(isite1)
        v1(3)=z(isite2)-z(isite1)
        dist = sqrt(dot_product(v1,v1))
        vard = dist - dd
        vard2 = vard**2 
        ebond = ebond + epsnuc*vard2*(1.0+100.0*vard2)
        ok = Qforces .and. stfree(isite1).or.stfree(isite2)
        if (ok) then
          de = epsnuc*vard*(2.0+400.0*vard2)/dist
          f1=de*v1
          if (stfree(isite1)) then   
            fx(isite1) = fx(isite1) + f1(1)
            fy(isite1) = fy(isite1) + f1(2)
            fz(isite1) = fz(isite1) + f1(3)
          endif
          if (stfree(isite2)) then
            fx(isite2) = fx(isite2) - f1(1)
            fy(isite2) = fy(isite2) - f1(2)
            fz(isite2) = fz(isite2) - f1(3)
          endif 
        endif
      enddo       
  !  The bending energy is calculated 
  !  by means of a Taylor expansion around natural bond
  !  angle. This contribution contains harmonic interactions. 
      do i = 1, nangle
        isite1 = siteangle(i,1)
        isite2 = siteangle(i,2) ! central site
        isite3 = siteangle(i,3) 
        ang0 = valangle(i) ! natural bond angle
        v1(1) = x(isite1) - x(isite2)
        v2(1) = x(isite3) - x(isite2)
        v1(2) = y(isite1) - y(isite2)
        v2(2) = y(isite3) - y(isite2)
        v1(3) = z(isite1) - z(isite2)
        v2(3) = z(isite3) - z(isite2)
        m1=dot_product(v1,v1) 
        m2=dot_product(v2,v2)
        m3=dot_product(v1,v2)
        modval=1.0/sqrt(m1*m2)
        ang = acos(m3*modval)
        varang = ang - ang0
        eang = eang + 700.0*epsnuc*varang**2
        ok = Qforces .and. stfree(isite1).or.stfree(isite2).or.stfree(isite3)
        if (ok) then
          de = 1400.0*epsnuc*varang*modval/sin(ang)
          f1=de*v1
          f2=de*v2
          if (stfree(isite1)) then
            fx(isite1) = fx(isite1) + f2(1)
            fy(isite1) = fy(isite1) + f2(2)
            fz(isite1) = fz(isite1) + f2(3)
          endif
          if (stfree(isite2)) then
            fx(isite2) = fx(isite2) - (f1(1)+f2(1))
            fy(isite2) = fy(isite2) - (f1(2)+f2(2))
            fz(isite2) = fz(isite2) - (f1(3)+f2(3))
          endif
          if (stfree(isite3)) then
            fx(isite3) = fx(isite3) + f1(1)
            fy(isite3) = fy(isite3) + f1(2)
            fz(isite3) = fz(isite3) + f1(3)
          endif
        endif
      enddo  
   ! The torsional energy is calculated 
   ! using Fourier series for natural torsional angle. 
      do i = 1, ndihe
        isite1 = sitedihe(i,1)
        isite2 = sitedihe(i,2)
        isite3 = sitedihe(i,3)
        isite4 = sitedihe(i,4)
        dihe0 = valdihe(i) ! natural torsion angle
        v1(1) = x(isite1) - x(isite2)
        v2(1) = x(isite3) - x(isite2)
        v3(1) = x(isite3) - x(isite4)
        v1(2) = y(isite1) - y(isite2)
        v2(2) = y(isite3) - y(isite2)
        v3(2) = y(isite3) - y(isite4)
        v1(3) = z(isite1) - z(isite2)
        v2(3) = z(isite3) - z(isite2)
        v3(3) = z(isite3) - z(isite4)
        call cross_product(v1,v2,n1)
        call cross_product(v2,v3,n2)
        v22=dot_product(v2,v2)
        v2v=sqrt(v22)
        av=dot_product(v1,n2)
        bv=dot_product(n1,n2)
        dihe = atan2(v2v*av,bv)
        vardihe = dihe - dihe0
        edihe = edihe + 28.0*epsnuc*(1.0-cos(vardihe))
        ok = Qforces .and. stfree(isite1).or.stfree(isite2).or.stfree(isite3).or.stfree(isite4)
        if (ok) then
          n1v=dot_product(n1,n1)
          n2v=dot_product(n2,n2)
          v12=dot_product(v1,v2)
          v23=dot_product(v2,v3)
          de=28.0*epsnuc*sin(vardihe)*v2v
          f1=-de/n1v*n1          
          f4=de/n2v*n2          
          iv22=1.0/v22
          v12=v12*iv22
          v23=v23*iv22
          f2=-f1+v12*f1-v23*f4
          f3=-f4-v12*f1+v23*f4
          if (stfree(isite1)) then
            fx(isite1)=fx(isite1) + f1(1)
            fy(isite1)=fy(isite1) + f1(2)
            fz(isite1)=fz(isite1) + f1(3)
          endif
          if (stfree(isite2)) then
            fx(isite2) = fx(isite2) + f2(1)
            fy(isite2) = fy(isite2) + f2(2)
            fz(isite2) = fz(isite2) + f2(3)
          endif
          if (stfree(isite3)) then
            fx(isite3) = fx(isite3) + f3(1)
            fy(isite3) = fy(isite3) + f3(2)
            fz(isite3) = fz(isite3) + f3(3)
          endif
          if (stfree(isite4)) then
            fx(isite4) = fx(isite4) + f4(1) 
            fy(isite4) = fy(isite4) + f4(2)
            fz(isite4) = fz(isite4) + f4(3)
          endif
        endif
      enddo
    endif  ! Qnobond 
    if (Qnonbond) then
    !Non-bonded interactions
    !-----------------------
    !Native contacts
      do i = 1, nstack
        isite1 = sitestack(i,1)
        isite2 = sitestack(i,2)
        v1(1)=x(isite2)-x(isite1)
        v1(2)=y(isite2)-y(isite1)
        v1(3)=z(isite2)-z(isite1)
        dist2 = dot_product(v1,v1)
        idist2=1.0/dist2
        cte1 = (sgstack(i)**2*idist2)**3   
        estack = estack + 4.0*epsnuc*cte1*(cte1-1.0)
        ok = Qforces .and. stfree(isite1).or.stfree(isite2)
        if (ok) then
          de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
          f1=de*v1 
          if (stfree(isite1)) then
            fx(isite1) = fx(isite1) - f1(1)  
            fy(isite1) = fy(isite1) - f1(2)
            fz(isite1) = fz(isite1) - f1(3)
          endif                            
          if (stfree(isite2)) then         
            fx(isite2) = fx(isite2) + f1(1)
            fy(isite2) = fy(isite2) + f1(2)
            fz(isite2) = fz(isite2) + f1(3)
          endif
        endif  
      enddo  
    ! Hydrogen bonding
      do i = 1, nbp
        isite1 = sitebp(i,1)
        isite2 = sitebp(i,2)
        if (namsite(isite1).eq.'Gb' .or. namsite(isite1).eq.'Cb') then
          epsln = 2.532*epsnuc*scalepairing
        else  
          epsln = 2.0*epsnuc*scalepairing
        endif
        v1(1)=x(isite2)-x(isite1)
        v1(2)=y(isite2)-y(isite1)
        v1(3)=z(isite2)-z(isite1)
        dist2 = dot_product(v1,v1)
        idist2=1.0/dist2
        cte1 = sgbp(i)*sgbp(i)*idist2
        ebp = ebp + epsln*cte1**5*(20.0*cte1-24.0) 
        ok = Qforces .and. stfree(isite1).or.stfree(isite2)
        if (ok) then
          de = 240.0*epsln*cte1**5*(cte1-1.0)*idist2
          f1=de*v1 
          if (stfree(isite1)) then
            fx(isite1) = fx(isite1) - f1(1)
            fy(isite1) = fy(isite1) - f1(2)
            fz(isite1) = fz(isite1) - f1(3)
          endif                            
          if (stfree(isite2)) then         
            fx(isite2) = fx(isite2) + f1(1)
            fy(isite2) = fy(isite2) + f1(2)
            fz(isite2) = fz(isite2) + f1(3)   
          endif
        endif                
      enddo
    ! Excluded volume
      do i = 1, nex
        isite1 = siteex(i,1)
        isite2 = siteex(i,2)
        v1(1)=x(isite2)-x(isite1)
        v1(2)=y(isite2)-y(isite1)
        v1(3)=z(isite2)-z(isite1)
        dist2 = dot_product(v1,v1)
        idist2=1.0/dist2
        sgex2=sgex(i)**2
        if (dist2.lt.sgex2) then
          cte1 = (sgex2*idist2)**3      
          eex = eex + 4.0*epsnuc*cte1*(cte1-1.0) + epsnuc
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then           
            de = 24.0*epsnuc*cte1*(2.0*cte1-1.0)*idist2
            f1=de*v1 
            if (stfree(isite1)) then
              fx(isite1) = fx(isite1) - f1(1)
              fy(isite1) = fy(isite1) - f1(2)
              fz(isite1) = fz(isite1) - f1(3)
            endif                            
            if (stfree(isite2)) then         
              fx(isite2) = fx(isite2) + f1(1)
              fy(isite2) = fy(isite2) + f1(2)
              fz(isite2) = fz(isite2) + f1(3)
            endif
          endif
        endif
      enddo 
    ! Coulomb interaction
      do i = 1, nqq
        isite1 = siteqq(i,1)
        isite2 = siteqq(i,2)
        v1(1)=x(isite2)-x(isite1)
        v1(2)=y(isite2)-y(isite1)
        v1(3)=z(isite2)-z(isite1)
        dist=sqrt(dot_product(v1,v1)) 
        idist=1.0/dist
        if (Qdebyhyb) then
          if(outbox(x(isite1),y(isite1),z(isite1)).and.outbox(x(isite2),y(isite2),z(isite2))) Qdeby=.true.
        endif
        if (Qdeby) then
          idistkd=exp(-dist*ikappa)
          eqq = eqq + fctn*idist*idistkd
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then
            de = fctn*(dist+kappa)*idist**3*ikappa*idistkd
            f1=de*v1
            if (stfree(isite1)) then
              fx(isite1) = fx(isite1) - f1(1)
              fy(isite1) = fy(isite1) - f1(2)
              fz(isite1) = fz(isite1) - f1(3)
            endif                            
            if (stfree(isite2)) then         
              fx(isite2) = fx(isite2) + f1(1)
              fy(isite2) = fy(isite2) + f1(2)
              fz(isite2) = fz(isite2) + f1(3) 
            endif
          endif      
        else
          eqq = eqq + fctn*idist
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then
            de = fctn*idist**3
            f1=de*v1
            if (stfree(isite1)) then
              fx(isite1) = fx(isite1) - f1(1)
              fy(isite1) = fy(isite1) - f1(2)
              fz(isite1) = fz(isite1) - f1(3)
            endif
            if (stfree(isite2)) then
              fx(isite2) = fx(isite2) + f1(1)
              fy(isite2) = fy(isite2) + f1(2)
              fz(isite2) = fz(isite2) + f1(3)
            endif
          endif
        endif
        if (Qdebyhyb) Qdeby=.false.
      enddo 
      !     Solvent-induced contribution
      if (Qsolv) then
        do i = 1, nsolv
          isite1 = siteslv(i,1)
          isite2 = siteslv(i,2)
          v1(1)=x(isite2)-x(isite1)
          v1(2)=y(isite2)-y(isite1)
          v1(3)=z(isite2)-z(isite1)
          dist=sqrt(dot_product(v1,v1)) 
          cte1 = exp((13.38-dist)*0.1875)
          esolv = esolv + epsolv*(1.0-cte1)**2 - epsolv
          ok = Qforces .and. stfree(isite1).or.stfree(isite2)
          if (ok) then
            de = 0.375*epsolv*cte1*(cte1-1.0)/dist
            f1=de*v1
            if (stfree(isite1)) then
              fx(isite1) = fx(isite1) - f1(1)
              fy(isite1) = fy(isite1) - f1(2)
              fz(isite1) = fz(isite1) - f1(3)
            endif                            
            if (stfree(isite2)) then         
              fx(isite2) = fx(isite2) + f1(1)
              fy(isite2) = fy(isite2) + f1(2)
              fz(isite2) = fz(isite2) + f1(3)
            endif
          endif
        enddo
      endif 
    endif ! Qnonbond
    ! APFOR    
    if (Qapfor) then
      do i=1,afn
        j=sn(i)
        if (Qforces.and.stfree(j)) then
          fx(j)=fx(j)+af(1,i)
          fy(j)=fy(j)+af(2,i)
          fz(j)=fz(j)+af(3,i)
        endif
      enddo
    endif
    if (Qcontrans) then
      do i=1,ctn
        j=csn(i)
        if (kx(i).ne.0.0) then
          if (j.eq.0) then
            if (Qunsplit) then
              xcon=sum(x(1:nsites1st))*insites*2.0-contrx(i)
              fx(1:nsites1st)=fx(1:nsites1st)-kx(i)*xcon*insites*2.0
              econ = econ + 0.5*kx(i)*xcon**2
              xcon=sum(x(1+nsites1st:nsites))*insites*2.0-contrx(ctn+1)
              fx(1+nsites1st:nsites)=fx(1+nsites1st:nsites)-kx(i)*xcon*insites*2.0
              econ = econ + 0.5*kx(i)*xcon**2
            else
              xcon=sum(x(1:nsites))*insites-contrx(i)
              fx(1:nsites)=fx(1:nsites)-kx(i)*xcon*insites
              econ = econ + 0.5*kx(i)*xcon**2
            endif
          else
            fx(j)=fx(j)-kx(i)*(x(j)-contrx(i))
            econ = econ + 0.5*kx(i)*(x(j)-contrx(i))**2
          endif 
        endif
        if (ky(i).ne.0.0) then
          if (j.eq.0) then
            if (Qunsplit) then
              ycon=sum(y(1:nsites1st))*insites*2.0-contry(i)
              fy(1:nsites1st)=fy(1:nsites1st)-ky(i)*ycon*insites*2.0
              econ = econ + 0.5*ky(i)*ycon**2
              ycon=sum(y(1+nsites1st:nsites))*insites*2.0-contry(ctn+1)
              fy(1+nsites1st:nsites)=fy(1+nsites1st:nsites)-ky(i)*ycon*insites*2.0
              econ = econ + 0.5*ky(i)*ycon**2
            else
              ycon=sum(y(1:nsites))*insites-contry(i)
              fy(1:nsites)=fy(1:nsites)-ky(i)*ycon*insites
              econ = econ + 0.5*ky(i)*ycon**2
            endif
          else
            fy(j)=fy(j)-ky(i)*(y(j)-contry(i))
            econ = econ + 0.5*ky(i)*(y(j)-contry(i))**2
          endif
        endif
        if (kz(i).ne.0.0) then
          if (j.eq.0) then
            if (Qunsplit) then
              zcon=sum(z(1:nsites1st))*insites*2.0-contrz(i)
              fz(1:nsites1st)=fz(1:nsites1st)-kz(i)*zcon*insites*2.0
              econ = econ + 0.5*kz(i)*zcon**2
              zcon=sum(z(1+nsites1st:nsites))*insites*2.0-contrz(ctn+1)
              fz(1+nsites1st:nsites)=fz(1+nsites1st:nsites)-kz(i)*zcon*insites*2.0
              econ = econ + 0.5*kz(i)*zcon**2
            else
              zcon=sum(z(1:nsites))*insites-contrz(i)
              fz(1:nsites)=fz(1:nsites)-kz(i)*zcon*insites
              econ = econ + 0.5*kz(i)*zcon**2
            endif
          else
            fz(j)=fz(j)-kz(i)*(z(j)-contrz(i))
            econ = econ + 0.5*kz(i)*(z(j)-contrz(i))**2
          endif
        endif
      enddo
    endif
  endif ! Qnuc
  dids(1:5,nsites+1:ntot)=0.0
  !  nonbonded interactions between interaction sites and ions
  if (Qpar .and. Qnucl .and. Qnonbond) then
    do j = 1, nsites
      jtype2 = typtyp(j)
      do i = nsites+1, ntot
        itype = abs(typei(i))
        itype2 = nwtype(itype)
        is=nindex(itype2,jtype2)
  ! Compute distances dna fragment-ion
        dist2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
        ok=.false.
        ok2=Qforces .and.(i.gt.(nsites+nfix).or.stfree(j))
        if (Qefpot(is)) then
          if (ok2) then
            call getyd(is,dist2,eefp,de,dist)
          else
            call gety(is,dist2,eefp,dist)
          endif
          eefpotmx=eefpotmx+eefp
          if (Qproxdiff) then
            if (dist.gt.0.0) then
              if (j.eq.1) then 
                dids(1,i)=dist-dmi(is) 
                dids(2,i)=dist 
                dids(3,i)=x(i)-x(j)
                dids(4,i)=y(i)-y(j)
                dids(5,i)=z(i)-z(j)
              else
                didstmp=dist-dmi(is)
                if (didstmp.lt.dids(1,i)) then
                  dids(1,i)=didstmp
                  dids(2,i)=dist
                  dids(3,i)=x(i)-x(j)
                  dids(4,i)=y(i)-y(j)
                  dids(5,i)=z(i)-z(j)
                endif
              endif
            endif
          endif
        else
          idist2 = 1.0/dist2
          if (Qchr(is).or.Qsrpmfi(is)) idist = sqrt(idist2)
  ! Compute Coulomb contribution
          if (Qchr(is)) then
            cofo=fct(is)*idist
            eqqmx = eqqmx + cofo
          endif
  ! Compute Lennard Jones contribution 
          if (Qlj(is)) then 
            dist6 = (sgp2(is)*idist2)**3
            dist12 = dist6**2
            evdwmx = evdwmx + epp4(is)*(dist12-dist6)
          endif
  ! Compute Short-Range correction contribution
          if (Qsrpmfi(is)) then
            if (dist2.le.rth) then
              dist=1.0/idist
              esrpmf1 = exp((c1(is)-dist)*c2(is))
              esrpmf2 = cos(c3(is)*pi*(c1(is)-dist))
              esrpmf3 = (c1(is)*idist)**6
              esrpmf0 = c0(is)*esrpmf1*esrpmf2+c4(is)*esrpmf3
              if (dist.ge.srpx) then ! smoothly fix discontinuity 
                fdf=exp(-srpk*(dist-srpx))-srpy
                fdv=esrpmf0
                esrpmf0=fdv*fdf
              endif
              esrpmfmx = esrpmfmx + esrpmf0
            endif
          endif
        endif
  ! Compute forces
        if (ok2) then    
          if (.not.Qefpot(is)) then
            de=0.0
            if (Qchr(is)) then
              de = cofo*idist2  ! Coulomb force
            endif
            if (Qlj(is)) then
              de = de+epp4(is)*(2.0*dist12-dist6)*6.0*idist2 ! vdw force
            endif
            if (Qsrpmfi(is)) then
              if (dist2.le.rth) then
                dc=(-c2(is)*esrpmf2 + c3(is)*pi*sin(c3(is)*pi*(c1(is)-dist)))*c0(is)*esrpmf1-6.0*c4(is)*esrpmf3*idist! forces
                if (dist.ge.srpx) dc=dc*fdf-fdv*srpk*(fdf+srpy)  ! smoothly fix discontinuity 
                de = de - dc*idist
              endif
            endif
          endif
          if (de.ne.0.0) then 
            if (i.gt.(nsites+nfix)) then
              fx(i) = fx(i) + de*(x(i)-x(j))
              fy(i) = fy(i) + de*(y(i)-y(j))
              fz(i) = fz(i) + de*(z(i)-z(j))
            endif
            if (stfree(j)) then
              fx(j) = fx(j) - de*(x(i)-x(j))
              fy(j) = fy(j) - de*(y(i)-y(j))
              fz(j) = fz(j) - de*(z(i)-z(j))
            endif
          endif
        endif
      enddo ! i=nsites+1,...,ntot 
    enddo   ! j=1,...,nsites
  endif
!write(*,*) 'Energy: ',egsbpb,egsbpa,evdwgd,eelec,evdw,esrpmf,esrpmfmx,ebond,eang,edihe,estack,ebp,eex,eqq,esolv,eqqmx,evdwmx,eefpot,eefpotmx,econ   ! debug
  ener = ener + eelec + evdw + esrpmf + esrpmfmx + ebond + eang + edihe + estack + ebp + eex + eqq + esolv + eqqmx + evdwmx + eefpot + eefpotmx + econ
endif                     !Qenergy
return
end
