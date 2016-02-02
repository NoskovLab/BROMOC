!    PB-PNP - Poisson-Boltzmann and Poisson-Nernst-Planck Equations Solver
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

      SUBROUTINE PROFILE(Ntype,Iontype,Temp,Nclx,Ncly,Nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           rsphe,xsphe,ysphe,zsphe,
     $           PHI,Cion,Zion,Diffusion,Qrho,Qflux,Qcurrent)
c------------------------------------------------------------------------
c Calculate density or flux profiles along Z-axis for each ion
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'
      integer*4 Ntype,NCLX,NCLY,NCLZ
      real*8  Temp,dcel,Zion(*),Diffusion(*)
      real*8  tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8  rsphe,xsphe,ysphe,zsphe
      real*4  PHI(*),Cion(*)
      CHARACTER*4 Iontype(*)
      LOGICAL Qrho,Qflux,Qcurrent
c local
      real*4  Zprofile(8)
      real*8  volume(8),area(8),NUMION(8),Isum(8),Isum2(8)
      real*8  Ic,Ic2,Icsum,Icsum2
      real*8  dcel2,dcel3,zz,factor,factor2
      real*8  phi0,phi6,avephi,c0,c6
      real*8  zmaxsphe,zminsphe,rsphe2,xx,yy,r2
      integer*4 ncyz,nc3,i,kg,ig,jg,iii,ip0,ipi,ipd
c

      if(Ntype.GT.8) then
         stop 'Increase dimension of Zprofile ... in PROFILE'
      endif

c

      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz
      dcel2=dcel*dcel
      dcel3=dcel2*dcel
      factor  = celec  / ( kboltz * Temp / kcalmol)
      zmaxsphe= rsphe+zsphe
      zminsphe=-rsphe+zsphe
      rsphe2  = rsphe*rsphe

c
      write(outu,'(6x,a)')
      write(outu,'(6x,a)') 
     $     'Ion accessible cross-sectional area along Z : [A^2]'
      write(outu,'(9x,a,7x,8(4x,a4,4x))') 'Z',(Iontype(i),i=1,ntype)

      do i=1,ntype
         NUMION(i)=0.0d0
         volume(i)=0.0d0
      enddo
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(rsphe.gt.0.0d0) then
            if(zz.gt.zmaxsphe.or.zz.lt.zminsphe) goto 1001
         endif
         do i=1,Ntype
            Zprofile(i)=0.0d0
            area(i)=0.0d0
         enddo
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz

               if(rsphe.gt.0.0d0) then
                  yy=(jg-1)*dcel-trany+ybcen
                  r2=xx*xx+yy*yy+zz*zz
                  if(r2.le.rsphe2) then
                     do i=1,ntype
                        ipi=ip0+nc3*(i-1)
                        if(Cion(ipi).ne.0.0d0) then
                           area(i)=area(i)+1.0d0
                           Zprofile(i)=Zprofile(i)+Cion(ipi)
                        endif
                     enddo
                  endif
               else
                  do i=1,ntype
                     ipi=ip0+nc3*(i-1)
                     if(Cion(ipi).ne.0.0d0) then
                        area(i)=area(i)+1.0d0
                        Zprofile(i)=Zprofile(i)+Cion(ipi)
                     endif
                  enddo
               endif
            enddo
         enddo
         do i=1,ntype
            volume(i)=volume(i)+area(i)*dcel3
            NUMION(i)=NUMION(i)+Zprofile(i)*dcel3
         enddo
      write(outu,'(3x,f10.3,2x,8f12.3)') zz,(area(i)*dcel2,i=1,ntype)
 1001 enddo

      write(outu,'(6x,a)')
      write(outu,'(40x,8(4x,a4,4x))') (Iontype(i),i=1,ntype)
      write(outu,'(6x,a,8f12.3)')
     $     'Ion accessible volume [Angs**3] :',(volume(i),i=1,ntype)
      write(outu,'(6x,a,8f12.3)')
     $     'Total number of ions            :',(NUMION(i),i=1,ntype)


c

      IF(Qrho) THEN

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Integrated number of ions along Z : ',
     $                      '[unit charge]'
      write(outu,'(6x,a)')
      write(outu,'(9x,a,8x,8(a4,8x))') 'Z',(Iontype(i),i=1,ntype)

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(rsphe.gt.0.0d0) then
            if(zz.gt.zmaxsphe.or.zz.lt.zminsphe) goto 1002
         endif
         do i=1,Ntype
            Zprofile(i)=0.0d0
         enddo
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz

               if(rsphe.gt.0.0d0) then
               yy=(jg-1)*dcel-trany+ybcen
               r2=xx*xx+yy*yy+zz*zz
                  if(r2.le.rsphe2) then
                     do i=1,ntype
                        ipi=ip0+nc3*(i-1)
                        if(Cion(ipi).ne.0.0d0) then
                           Zprofile(i)=Zprofile(i)+Cion(ipi)
                        endif
                     enddo
                  endif
               else
                  do i=1,ntype
                     ipi=ip0+nc3*(i-1)
                     if(Cion(ipi).ne.0.0d0) then
                        Zprofile(i)=Zprofile(i)+Cion(ipi)
                     endif
                  enddo
               endif
            enddo
         enddo
         write(outu,'(3x,f10.3,2x,8(e11.4,2x))') 
     $        zz,(Zprofile(i)*dcel3,i=1,ntype)
 1002 enddo


      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 
     $     'Density Profile based on ion-accessible area along Z : ',
     $     '[unit charge]/[A^3]'
      write(outu,'(6x,a)')
      write(outu,'(9x,a,8x,8(a4,8x))') 'Z',(Iontype(i),i=1,ntype)

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(rsphe.gt.0.0d0) then
            if(zz.gt.zmaxsphe.or.zz.lt.zminsphe) goto 1003
         endif
         do i=1,Ntype
            area(i)=0.0d0
            Zprofile(i)=0.0d0
         enddo
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz

               if(rsphe.gt.0.0d0) then 
                  yy=(jg-1)*dcel-trany+ybcen
                  r2=xx*xx+yy*yy+zz*zz
                  if(r2.le.rsphe2) then
                     do i=1,ntype
                        ipi=ip0+nc3*(i-1)
                        if(Cion(ipi).ne.0.0d0) then
                           area(i)=area(i)+1.0d0
                           Zprofile(i)=Zprofile(i)+Cion(ipi)
                        endif
                     enddo
                  endif
               else
                  do i=1,ntype
                     ipi=ip0+nc3*(i-1)
                     if(Cion(ipi).ne.0.0d0) then
                        area(i)=area(i)+1.0d0
                        Zprofile(i)=Zprofile(i)+Cion(ipi)
                     endif
                  enddo
               endif
            enddo
         enddo
         write(outu,'(3x,f10.3,2x,8(e11.4,2x))') 
     $        zz,(Zprofile(i)/(area(i)+rsmall),i=1,ntype)
 1003 enddo


      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 
     $     'Density Profile based on system area along Z : ',
     $     '[unit charge]/[A^3]'
      write(outu,'(6x,a)')
      write(outu,'(9x,a,8x,8(a4,8x))') 'Z',(Iontype(i),i=1,ntype)

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         do i=1,Ntype
            area(i)=0.0d0 
            Zprofile(i)=0.0d0
            if(rsphe.eq.0.0d0) area(i)=nclx*ncly
         enddo
         if(rsphe.gt.0.0d0) then
            if(zz.gt.zmaxsphe.or.zz.lt.zminsphe) goto 1004
         endif
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz

               if(rsphe.gt.0.0d0) then 
                  yy=(jg-1)*dcel-trany+ybcen
                  r2=xx*xx+yy*yy+zz*zz
                  if(r2.le.rsphe2) then
                     do i=1,ntype
                        ipi=ip0+nc3*(i-1)
                        if(Cion(ipi).ne.0.0d0) then
                           area(i)=area(i)+1.0d0
                           Zprofile(i)=Zprofile(i)+Cion(ipi)
                        endif
                     enddo
                  endif
               else
                  do i=1,ntype
                     ipi=ip0+nc3*(i-1)
                     if(Cion(ipi).ne.0.0d0) then
                        Zprofile(i)=Zprofile(i)+Cion(ipi)
                     endif
                  enddo
               endif
            enddo
         enddo

         write(outu,'(3x,f10.3,2x,8(e11.4,2x))') 
     $        zz,(Zprofile(i)/(area(i)+rsmall),i=1,ntype)
 1004 enddo

      ENDIF      ! Qrho

c

      IF(Qflux) THEN

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Flux (Jz) Profile along Z : ',
     $                      '[unit charge]/[ps][Angs^2]'
      write(outu,'(6x,a)')
      write(outu,'(9x,a,8x,8(a4,8x))') 'Z',(Iontype(i),i=1,ntype)

      do kg=1,nclz-1                         ! exclude upper boundary points
         zz=(kg-1)*dcel-tranz+zbcen
         do i=1,Ntype
            area(i)=0.0d0
            Zprofile(i)=0.0d0
         enddo
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               
               do i=1,ntype
                  ipi=ip0+(i-1)*nc3
                  if(Cion(ipi).ne.0.0d0.and.Cion(ipi+1).ne.0.0d0) then
                     area(i)=area(i)+1.0d0
                     phi0=phi(ip0)  *factor
                     phi6=phi(ip0+1)*factor
                     avephi=(phi0+phi6)*0.5d0
                     c0=Cion(ipi)  *exp(phi0*zion(i)) !conc. -> effective conc.
                     c6=Cion(ipi+1)*exp(phi6*zion(i)) !conc. -> effective conc.
                     ipd=kg+(i-1)*nclz
                     Zprofile(i)=Zprofile(i)-
     $                  Diffusion(ipd)*exp(-avephi*zion(i))*(c6-c0)/dcel
                  endif
               enddo

            enddo
         enddo
         write(outu,'(3x,f10.3,2x,8(e11.4,2x))') 
     $        zz+dcel*0.5d0,(Zprofile(i)/(area(i)+rsmall),i=1,ntype)
      enddo

      ENDIF   ! Qflux

c

      IF(Qcurrent) THEN

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Current Profile along Z : ','[pA]'
      write(outu,'(6x,a)')
      write(outu,'(9x,a,8x,8(a4,8x))') 'Z',(Iontype(i),i=1,ntype)

      factor2=coulomb/psec/pico*dcel2
      Icsum2=0.0d0
      do i=1,ntype
         Isum(i)=0.0d0
         Isum2(i)=0.0d0
      enddo
      do kg=1,nclz-1                         ! exclude upper boundary points
         zz=(kg-1)*dcel-tranz+zbcen
         do i=1,Ntype
            Zprofile(i)=0.0d0
         enddo
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               
               do i=1,ntype
                  ipi=ip0+(i-1)*nc3
                  if(Cion(ipi).ne.0.0d0.and.Cion(ipi+1).ne.0.0d0) then
                     phi0=phi(ip0)  *factor
                     phi6=phi(ip0+1)*factor
                     avephi=(phi0+phi6)*0.5d0
                     c0=Cion(ipi)  *exp(phi0*zion(i)) !conc. -> effective conc.
                     c6=Cion(ipi+1)*exp(phi6*zion(i)) !conc. -> effective conc.
                     ipd=kg+(i-1)*nclz
                     Zprofile(i)=Zprofile(i)-
     $                  Diffusion(ipd)*exp(-avephi*zion(i))*(c6-c0)/dcel
                  endif
               enddo

            enddo
         enddo
         write(outu,'(3x,f10.3,2x,8(e11.4,2x))') 
     $        zz+dcel*0.5d0,(Zprofile(i)*Zion(i)*factor2,i=1,ntype)
         Icsum=0.0d0
         do i=1,Ntype
            Ic=Zprofile(i)*Zion(i)*factor2
            Isum(i)=Isum(i)+Ic
            Isum2(i)=Isum2(i)+Ic*Ic
            Icsum=Icsum+Ic
         enddo
         Icsum2=Icsum2+Icsum*Icsum
      enddo

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Current Statistics : '
      Icsum=0.0d0
      do i=1,ntype
         Ic=Isum(i)/(nclz-1)
         Ic2=Isum2(i)/(nclz-1)
         write(outu,'(6x,a4,a,f10.4,a,f10.4,a)') 
     $        iontype(i),'  : ',Ic,' +/-',sqrt(Ic2-Ic*Ic),' [pA]'
         Icsum=Icsum+Isum(i)
      enddo
      Icsum=Icsum/(nclz-1)
      Icsum2=Icsum2/(nclz-1)
      write(outu,'(6x,a,f10.4,a,f10.4,a)') 
     $        'TOTAL : ',Icsum,' +/-',sqrt(Icsum2-Icsum*Icsum),' [pA]'

      ENDIF   ! Qcurrent
c
      return
      end


      SUBROUTINE COUNTERION1(NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $           VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN,CONC,TEMP,
     $           Qnonlinear,Qpartlinear)
c-----------------------------------------------------------------------
c This subroutine computes the number of counter ions
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'
      integer*4 NCLX,NCLY,NCLZ
      real*8  dcel,tranz,zbcen,conc,temp
      real*8  vmemb,tmemb,zmemb 
      real*4  phi(*),mcden(*)
      logical*1 Qnonlinear,Qpartlinear
c local
      real*8  factor1
      real*8  bulk_num,volume,nion_num,pion_num,bulk_rho
      real*8  nion_numzp,pion_numzp,nion_numzn,pion_numzn
      real*8  nion_nummb,pion_nummb,nfactor,pfactor,area
      real*8  dcel2,dcel3,zz,zc,zmemb2,phif
      integer*4 ip0,nc3,ncyz,ig,jg,kg,iii
c
      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz
      dcel2=dcel*dcel
      dcel3=dcel2*dcel
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall

c conversion factoer from 1/(kcal/(mol*e)) to 1/(e/A)
      factor1 = celec  / ( kboltz * Temp / kcalmol)

c calculate ion accessible volume
      volume=0.0D0
      do ip0=1,nc3
         if(mcden(ip0).ne.0.0d0) volume=volume+dcel3
      enddo

c bulk density and number of counter ions
      bulk_rho=conc*(avogadro/liter)*(angstrom**3) ! [unit charge]/A^3
      bulk_num=bulk_rho*volume

c deviation from bulk number of count ions
      pion_num=0.0D0
      nion_num=0.0D0
      pion_numzp=0.0D0
      nion_numzp=0.0D0
      pion_numzn=0.0D0
      nion_numzn=0.0D0
      pion_nummb=0.0D0
      nion_nummb=0.0D0
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            do 101 jg=1,ncly
               ip0=iii+(jg-1)*nclz
               if(mcden(ip0).ne.0.0d0) then
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
                  if(zz.ge.0.0d0) then
                     pion_numzp=pion_numzp+pfactor
                     nion_numzp=nion_numzp+nfactor
                  else
                     pion_numzn=pion_numzn+pfactor
                     nion_numzn=nion_numzn+nfactor
                  endif
               endif
 101        enddo
         enddo
      enddo
c
      write(outu,'(6x,a)')
      write(outu,'(6x,a,f13.5,a)')
     $     'Ion accessible volume                   :',volume,
     $     ' [Angs**3]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Bulk density                            :',bulk_rho,
     $     ' [unit charge]/[Angs**3]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Bulk number of counter ions             :',bulk_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions                 :',pion_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions (Z > 0)         :',pion_numzp,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions (Z < 0)         :',pion_numzn,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions                 :',nion_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions (Z > 0)         :',nion_numzp,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions (Z < 0)         :',nion_numzn,
     $     ' [unit charge]'


c  Counter Ion Distributions
      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Integrated number of ions along Z : ',
     $                      '[unit charge]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         pion_num=0.0D0
         nion_num=0.0D0
         area=0.0d0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               if(mcden(ip0).ne.0.0d0) then
                  area=area+1.0d0
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
            enddo
         enddo
         write(outu,'(6x,2f10.3,2x,2(e11.4,2x))') 
     $        zz,area*dcel2,pion_num,nion_num
      enddo

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 
     $     'Density Profile based on ion-accessible area along Z : ',
     $     '[unit charge]/[A^3]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         pion_num=0.0D0
         nion_num=0.0D0
         area=0.0d0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               if(mcden(ip0).ne.0.0d0) then
                  area=area+1.0d0
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
            enddo
         enddo
         write(outu,'(6x,2f10.3,2x,2(e11.4,2x))') 
     $        zz,area*dcel2,
     $        pion_num/(area*dcel3+rsmall),nion_num/(area*dcel3+rsmall)
      enddo

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 
     $     'Density Profile based on system area along Z : ',
     $     '[unit charge]/[A^3]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      area=nclx*ncly
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         pion_num=0.0D0
         nion_num=0.0D0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            zz=(kg-1)*dcel-tranz+zbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               if(mcden(ip0).ne.0.0d0) then
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
            enddo
         enddo
         write(outu,'(6x,2f10.3,2x,2(e11.4,2x))') 
     $        zz,area*dcel2,
     $        pion_num/(area*dcel3+rsmall),nion_num/(area*dcel3+rsmall)
      enddo
c
      RETURN
      END


      SUBROUTINE COUNTERION2(NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $           VMEMB,TMEMB,ZMEMB,CONC,TEMP,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           rsphe,xsphe,ysphe,zsphe,rdist,
     $           Qnonlinear,Qpartlinear)
c-----------------------------------------------------------------------
c This subroutine computes the number of counter ions
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'
      integer*4 NCLX,NCLY,NCLZ
      real*8  dcel,conc,temp,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8  rsphe,xsphe,ysphe,zsphe,rdist
      real*8  vmemb,tmemb,zmemb 
      real*4  phi(*),mcden(*)
      logical*1 Qnonlinear,Qpartlinear
c local
      real*8  factor1
      real*8  bulk_num,volume,nion_num,pion_num,bulk_rho
      real*8  nion_numzp,pion_numzp,nion_numzn,pion_numzn
      real*8  nion_nummb,pion_nummb,nfactor,pfactor,area
      real*8  dcel2,dcel3,zz,zc,zmemb2,phif
      integer*4 ip0,nc3,ncyz,ig,jg,kg,iii
      real*8  zmaxsphe,zminsphe,rsphe2,xx,yy,r2,rdist2
c
      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz
      dcel2=dcel*dcel
      dcel3=dcel2*dcel
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall
      zmaxsphe= rsphe+zsphe
      zminsphe=-rsphe+zsphe
      rsphe2  = rsphe*rsphe
      rdist2  = rdist*rdist

c conversion factoer from 1/(kcal/(mol*e)) to 1/(e/A)
      factor1 = celec  / ( kboltz * Temp / kcalmol)

c calculate ion accessible volume
      volume=0.0D0
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if((zz.gt.zmaxsphe.or.zz.lt.zminsphe).and.rdist.eq.0.0d0)
     $        goto 1001
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
                  if(rsphe.gt.0.0d0) then
                     r2=xx*xx+yy*yy+zz*zz
                     if(r2.le.rsphe2) volume=volume+dcel3
                     goto 101
                  endif
                  if(rdist.gt.0.0d0) then
                     r2=xx*xx+yy*yy
                     if(r2.le.rdist2) volume=volume+dcel3
                     goto 101
                  endif
               endif

 101        enddo
         enddo
 1001 enddo

c bulk density and number of counter ions
      bulk_rho=conc*(avogadro/liter)*(angstrom**3) ! [unit charge]/A^3
      bulk_num=bulk_rho*volume

c deviation from bulk number of count ions
      pion_num=0.0D0
      nion_num=0.0D0
      pion_numzp=0.0D0
      nion_numzp=0.0D0
      pion_numzn=0.0D0
      nion_numzn=0.0D0
      pion_nummb=0.0D0
      nion_nummb=0.0D0
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if((zz.gt.zmaxsphe.or.zz.lt.zminsphe).and.rdist.eq.0.0d0)
     $        goto 1002
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
                  if(rsphe.gt.0.0d0) then
                     r2=xx*xx+yy*yy+zz*zz
                     if(r2.gt.rsphe2) goto 102
                  endif
                  if(rdist.gt.0.0d0) then
                     r2=xx*xx+yy*yy
                     if(r2.gt.rdist2) goto 102
                  endif
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
                  if(zz.ge.0.0d0) then
                     pion_numzp=pion_numzp+pfactor
                     nion_numzp=nion_numzp+nfactor
                  else
                     pion_numzn=pion_numzn+pfactor
                     nion_numzn=nion_numzn+nfactor
                  endif
               endif
 102        enddo
         enddo
 1002 enddo
c
      write(outu,'(6x,a)')
      write(outu,'(6x,a,f13.5,a)')
     $     'Ion accessible volume                   :',volume,
     $     ' [Angs**3]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Bulk density                            :',bulk_rho,
     $     ' [unit charge]/[Angs**3]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Bulk number of counter ions             :',bulk_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions                 :',pion_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions (Z > 0)         :',pion_numzp,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions (Z < 0)         :',pion_numzn,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions                 :',nion_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions (Z > 0)         :',nion_numzp,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions (Z < 0)         :',nion_numzn,
     $     ' [unit charge]'


c  Counter Ion Distributions
      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Integrated number of ions along Z : ',
     $                      '[unit charge]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if((zz.gt.zmaxsphe.or.zz.lt.zminsphe).and.rdist.eq.0.0d0)
     $        goto 1003
         pion_num=0.0D0
         nion_num=0.0D0
         area=0.0d0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
                  if(rsphe.gt.0.0d0) then
                     r2=xx*xx+yy*yy+zz*zz
                     if(r2.gt.rsphe2) goto 103
                  endif
                  if(rdist.gt.0.0d0) then
                     r2=xx*xx+yy*yy
                     if(r2.gt.rdist2) goto 103
                  endif
                  area=area+1.0d0
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
 103        enddo
         enddo
         write(outu,'(6x,2f10.3,2x,2(e11.4,2x))') 
     $        zz,area*dcel2,pion_num,nion_num
 1003 enddo

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 
     $     'Density Profile based on ion-accessible area along Z : ',
     $     '[unit charge]/[A^3]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if((zz.gt.zmaxsphe.or.zz.lt.zminsphe).and.rdist.eq.0.0d0)
     $        goto 1004
         pion_num=0.0D0
         nion_num=0.0D0
         area=0.0d0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
                  if(rsphe.gt.0.0d0) then
                     r2=xx*xx+yy*yy+zz*zz
                     if(r2.gt.rsphe2) goto 104
                  endif
                  if(rdist.gt.0.0d0) then
                     r2=xx*xx+yy*yy
                     if(r2.gt.rdist2) goto 104
                  endif
                  area=area+1.0d0
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
 104        enddo
         enddo
         write(outu,'(6x,2f10.3,2x,2(e11.4,2x))') 
     $        zz,area*dcel2,
     $        pion_num/(area*dcel3+rsmall),nion_num/(area*dcel3+rsmall)
 1004 enddo

      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 
     $     'Density Profile based on system area along Z : ',
     $     '[unit charge]/[A^3]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if((zz.gt.zmaxsphe.or.zz.lt.zminsphe).and.rdist.eq.0.0d0)
     $        goto 1005
         area=0.0d0
         pion_num=0.0D0
         nion_num=0.0D0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
                  if(rsphe.gt.0.0d0) then
                     r2=xx*xx+yy*yy+zz*zz
                     if(r2.gt.rsphe2) goto 105
                  endif
                  if(rdist.gt.0.0d0) then
                     r2=xx*xx+yy*yy
                     if(r2.gt.rdist2) goto 105
                  endif
                  area=area+1.0d0
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
 105        enddo
         enddo
         write(outu,'(6x,2f10.3,2x,2(e11.4,2x))') 
     $        zz,area*dcel2,
     $        pion_num/(area*dcel3+rsmall),nion_num/(area*dcel3+rsmall)
 1005 enddo
c
      RETURN
      END


      SUBROUTINE PORIN(NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $           VMEMB,TMEMB,ZMEMB,CONC,TEMP,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           Qnonlinear,Qpartlinear)
c-----------------------------------------------------------------------
c This subroutine computes the number of counter ions
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'
      integer*4 NCLX,NCLY,NCLZ
      real*8  dcel,conc,temp,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8  vmemb,tmemb,zmemb 
      real*4  phi(*),mcden(*)
      logical*1 Qnonlinear,Qpartlinear
c local
      real*8  factor1
      real*8  bulk_num,volume,nion_num,pion_num,bulk_rho
      real*8  nion_numzp,pion_numzp,nion_numzn,pion_numzn
      real*8  nion_nummb,pion_nummb,nfactor,pfactor,area
      real*8  dcel2,dcel3,zz,zc,zmemb2,phif
      integer*4 ip0,nc3,ncyz,ig,jg,kg,iii
      real*8  xx,yy

      real*8  area1,area2,area3
      real*8  minarea1,minarea2,minarea3
      real*8  minzz1,minzz2,minzz3

      real*8  xave1,yave1,rave1,cost1,sint1
      real*8  xave2,yave2,rave2,cost2,sint2
      real*8  xave3,yave3,rave3,cost3,sint3
      real*8  r2,xnew,ynew,ellip1,ellip2,ellip3
c
      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz
      dcel2=dcel*dcel
      dcel3=dcel2*dcel
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall

c
c pore1         
      xave1 = 21.600829d0
      yave1 = -7.021551d0
      rave1 = sqrt(xave1*xave1+yave1*yave1)
      cost1 = xave1/rave1
      sint1 = yave1/rave1
      rave1 = rave1 - 3.0d0
      xave1 = rave1*cost1
      yave1 = rave1*sint1
c pore2
      xave2 =-16.881256d0
      yave2 =-15.196091d0
      rave2 = sqrt(xave2*xave2+yave2*yave2)
      cost2 = xave2/rave2
      sint2 = yave2/rave2
      rave2 = rave2 - 3.0d0
      xave2 = rave2*cost2
      yave2 = rave2*sint2
c pore3
      xave3 = -4.719573d0
      yave3 = 22.217642d0
      rave3 = sqrt(xave3*xave3+yave3*yave3)
      cost3 = xave3/rave3
      sint3 = yave3/rave3
      rave3 = rave3 - 3.0d0
      xave3 = rave3*cost3
      yave3 = rave3*sint3

c conversion factoer from 1/(kcal/(mol*e)) to 1/(e/A)
      factor1 = celec  / ( kboltz * Temp / kcalmol)

c calculate ion accessible volume
      volume=0.0D0
      minarea1=10000.0d0
      minarea2=10000.0d0
      minarea3=10000.0d0
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(zz.gt.16.0d0.or.zz.lt.-16.0d0) goto 1001
         area1=0.0D0
         area2=0.0D0
         area3=0.0D0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
c     pore1
                  xnew  = (xx-xave1)*cost1+(yy-yave1)*sint1
                  ynew  = (yy-yave1)*cost1-(xx-xave1)*sint1
                  ellip1= xnew*xnew/225d0 + ynew*ynew/324d0
c     pore2
                  xnew  = (xx-xave2)*cost2+(yy-yave2)*sint2
                  ynew  = (yy-yave2)*cost2-(xx-xave2)*sint2
                  ellip2= xnew*xnew/225d0 + ynew*ynew/324d0
c     pore3
                  xnew  = (xx-xave3)*cost3+(yy-yave3)*sint3
                  ynew  = (yy-yave3)*cost3-(xx-xave3)*sint3
                  ellip3= xnew*xnew/225d0 + ynew*ynew/324d0
c     common pore
                  r2=xx*xx+yy*yy

                  if(ellip1.le.1.0d0) area1=area1+1.0d0
                  if(ellip2.le.1.0d0) area2=area2+1.0d0
                  if(ellip3.le.1.0d0) area3=area3+1.0d0
                  
                  if(ellip1.gt.1.0d0.and.
     $               ellip2.gt.1.0d0.and.
     $               ellip3.gt.1.0d0.and.r2.gt.20.0d0) then
                     mcden(ip0)=0.0d0
                  else
                     volume=volume+dcel3
                  endif
               endif

            enddo
         enddo
         if(area1.lt.minarea1) then
            minarea1=area1
            minzz1=zz
         endif
         if(area2.lt.minarea2) then 
            minarea2=area2
            minzz2=zz
         endif
         if(area3.lt.minarea3) then
            minarea3=area3
            minzz3=zz
         endif
         write(71,'(6x,4f10.3)') 
     $        zz,area1*dcel2,area2*dcel2,area3*dcel2
 1001 enddo

      write(72,'(6x,6f10.3)') minzz1,minarea1*dcel2,
     $     minzz2,minarea2*dcel2,minzz3,minarea3*dcel2


c bulk density and number of counter ions
      bulk_rho=conc*(avogadro/liter)*(angstrom**3) ! [unit charge]/A^3
      bulk_num=bulk_rho*volume

c deviation from bulk number of count ions
      pion_num=0.0D0
      nion_num=0.0D0
      pion_numzp=0.0D0
      nion_numzp=0.0D0
      pion_numzn=0.0D0
      nion_numzn=0.0D0
      pion_nummb=0.0D0
      nion_nummb=0.0D0
      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(zz.gt.16.0d0.or.zz.lt.-16.0d0) goto 1002
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen

               if(mcden(ip0).ne.0.0d0) then
c     pore1
                  xnew  = (xx-xave1)*cost1+(yy-yave1)*sint1
                  ynew  = (yy-yave1)*cost1-(xx-xave1)*sint1
                  ellip1= xnew*xnew/225d0 + ynew*ynew/324d0
c     pore2
                  xnew  = (xx-xave2)*cost2+(yy-yave2)*sint2
                  ynew  = (yy-yave2)*cost2-(xx-xave2)*sint2
                  ellip2= xnew*xnew/225d0 + ynew*ynew/324d0
c     pore3
                  xnew  = (xx-xave3)*cost3+(yy-yave3)*sint3
                  ynew  = (yy-yave3)*cost3-(xx-xave3)*sint3
                  ellip3= xnew*xnew/225d0 + ynew*ynew/324d0
c     common pore
                  r2=xx*xx+yy*yy
                  
                  if(ellip1.gt.1.0d0.and.
     $               ellip2.gt.1.0d0.and.
     $               ellip3.gt.1.0d0.and.r2.gt.20.0d0) goto 102

                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
                  if(zz.ge.0.0d0) then
                     pion_numzp=pion_numzp+pfactor
                     nion_numzp=nion_numzp+nfactor
                  else
                     pion_numzn=pion_numzn+pfactor
                     nion_numzn=nion_numzn+nfactor
                  endif
               endif
 102        enddo
         enddo
 1002 enddo
c
      write(outu,'(6x,a)')
      write(outu,'(6x,a,f13.5,a)')
     $     'Ion accessible volume                   :',volume,
     $     ' [Angs**3]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Bulk density                            :',bulk_rho,
     $     ' [unit charge]/[Angs**3]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Bulk number of counter ions             :',bulk_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions                 :',pion_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions (Z > 0)         :',pion_numzp,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of positive ions (Z < 0)         :',pion_numzn,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions                 :',nion_num,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions (Z > 0)         :',nion_numzp,
     $     ' [unit charge]'
      write(outu,'(6x,a,f13.5,a)')
     $     'Number of negative ions (Z < 0)         :',nion_numzn,
     $     ' [unit charge]'


c  Counter Ion Distributions
      write(outu,'(6x,a)')
      write(outu,'(6x,2a)') 'Integrated number of ions along Z : ',
     $                      '[unit charge]'
      write(outu,'(6x,a)')
      write(outu,'(11x,a,7x,a,7x,a,7x,a)') 
     $        'Z','AREA','+ ION','- ION'

      do kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(zz.gt.16.0d0.or.zz.lt.-16.0d0) goto 1003
         pion_num=0.0D0
         nion_num=0.0D0
         area=0.0d0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen
               r2=xx*xx+yy*yy+zz*zz
               if(mcden(ip0).ne.0.0d0) then
c     pore1
                  xnew  = (xx-xave1)*cost1+(yy-yave1)*sint1
                  ynew  = (yy-yave1)*cost1-(xx-xave1)*sint1
                  ellip1= xnew*xnew/225d0 + ynew*ynew/324d0
c     pore2
                  xnew  = (xx-xave2)*cost2+(yy-yave2)*sint2
                  ynew  = (yy-yave2)*cost2-(xx-xave2)*sint2
                  ellip2= xnew*xnew/225d0 + ynew*ynew/324d0
c     pore3
                  xnew  = (xx-xave3)*cost3+(yy-yave3)*sint3
                  ynew  = (yy-yave3)*cost3-(xx-xave3)*sint3
                  ellip3= xnew*xnew/225d0 + ynew*ynew/324d0
c     common pore
                  r2=xx*xx+yy*yy
                  
                  if(ellip1.gt.1.0d0.and.
     $               ellip2.gt.1.0d0.and.
     $               ellip3.gt.1.0d0.and.r2.gt.20.0d0) goto 103

                  area=area+1.0d0
                  phif=phi(ip0)*factor1
                  if(vmemb.ne.0.0d0)then
                     zc=(kg-1)*dcel
                     if(zc.gt.zmemb2) phif=phif-vmemb*factor1
                  endif
                  if(Qnonlinear) then
                     pfactor=bulk_rho*dcel3*exp(-phif)
                     nfactor=bulk_rho*dcel3*exp(+phif)
                  elseif(Qpartlinear) then
                     if(phif.gt.0.0d0) then
                        pfactor=bulk_rho*dcel3*exp(-phif)
                        nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                     else
                        pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                        nfactor=bulk_rho*dcel3*exp(+phif)
                     endif
                  else
                     pfactor=bulk_rho*dcel3*(1.0d0 - phif)
                     nfactor=bulk_rho*dcel3*(1.0d0 + phif)
                  endif
                  pion_num=pion_num+pfactor
                  nion_num=nion_num+nfactor
               endif
 103        enddo
         enddo
         write(outu,'(6x,2f10.3,2x,4(e11.4,2x))') 
     $     zz,area*dcel2,pion_num,nion_num,pion_num/3.0d0,nion_num/3.0d0
 1003 enddo
c
      RETURN
      END


