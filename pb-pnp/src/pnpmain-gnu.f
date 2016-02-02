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

c     ==============
      PROGRAM PNP0
c     ==============
      use charfuncmod
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'pnp.fcm'
      include 'consta.fcm'
      include 'misc.fcm'     !command parser
      include 'mainio.fcm'

      character com*2048,fnam*80,wrd5*5,wrd4*4,rtype*5,fnam1*80
      character*128 line
      integer*4   i,j,iunit,itype,lfnam,punit,new_unit,lfnam1
c for GBRAD      
      real*8    rmax,volume
c for diffusion_constant
      real*8    zz,zmin,zmax,diffus,aaa
      integer*4   kg,ipd,n12
c for time
      real*8       start,finish,cputime
!      character*8  hour ! for ifort
      character*24  hour ! for gfortran
      logical*1      EOF

      write(outu,*)
      write(outu,101)
     $   'POISSON-NERNST-PLANCK AND POISSON-BOLTZMANN EQUATIONS SOLVER'
      write(outu,101)
     $   '============================================================'
      start = cputime(0.0D0)
      call time(hour)
      write(outu,107) 'started at: ',hour
      write(outu,*)
      write(outu,*)

 101  format(6x,a)

c initialize all parameters
      call pnpinit

      data inpu,outu/5,6/
      frmt = ' '

c read input
      call rdtitl(inpu)
 1000 call getlin('CMD> ',com,inpu,outu)
      call getfirst(com,wrd5)
      call misc(com,strg,qstrg,wrd5)
      call upper(wrd5)

c.......................................................................
      if(wrd5.eq.'     ')then
c        ---------------
         goto 1000

c.......................................................................
      elseif(wrd5.eq.'IONTY')then                              ! IONTYPE
c            ---------------
 1001    call getlin('---> ',com,inpu,outu)
         call misc(com,strg,qstrg,wrd5)
         if(.not.check(com,'end'))then
            ntype = ntype + 1
            if(ntype.gt.dtype) then
               write(outu,101) 
     $              'Number of ion types exceeds dtype in pnp.fcm'
               write(outu,'(1x,a,i10,a,i10)') 
     $              'dtype :',dtype,'  Number of ion types :',ntype
               stop
            endif
            call getfirst(com,iontype(ntype))
            call gtdpar(com,'charge',Zion(ntype),0.0d0)
            call gtdpar(com,'radius',Rion(ntype),0.0d0)            
            goto 1001
         endif

         write(outu,*)
         write(outu,'(6x,a,i3,a)') 'There are ',ntype,' ion types'
         write(outu,'(16x,a,2x,a)') 'charge [e] ','radius [Ang.] '
         do i=1,ntype 
            write(outu,'(6x,i3,1x,a,1x,f8.3,5x,f8.3)')
     $           i,iontype(i),Zion(i),Rion(i)
         enddo
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'DIFFU')then                   ! DIFFUSION_CONSTANT
c            ---------------
         call gtipar(com,'nclz',nclz,65)
         call gtdpar(com,'dcel',dcel,0.5D0) 
         call gtdpar(com,'zbcen',zbcen,0.0D0) 
         tranz = 0.5d0*(nclz-1)*dcel
         call gtipar(com,'dunit',iunit,0)
         if(iunit.eq.0) then
 1002       call getlin('---> ',com,inpu,outu)
            call misc(com,strg,qstrg,wrd5)
            if(.not.check(com,'end'))then
               call getfirst(com,wrd4)
               call fatnam(iontype,ntype,wrd4,itype,.false.)
               call gtdpar(com,'zmax',zmax,0.0D0) 
               call gtdpar(com,'zmin',zmin,0.0D0) 
               call gtdpar(com,'diffusion',diffus,0.1d0) ! [A^2/ps]
               do kg=1,nclz
                  zz=(kg-1)*dcel-tranz+zbcen
                  ipd=(itype-1)*nclz+kg
                  if(zz.ge.zmin.and.zz.le.zmax) Diffusion(ipd)=diffus
               enddo
               goto 1002
            endif
         else
            write(outu,102) 
     $           'Reading DIFFUSION_CONSTANT from unit ',iunit
            do kg=1,nclz
               read(iunit,*) zz,(Diffusion((i-1)*nclz+kg),i=1,ntype)
            enddo
         endif

         write(outu,*)
         write(outu,101)
     $        'Diffusion constant in [A^2/ps] along Z-axis'
         write(outu,'(9x,a,9x,8(a4,8x))') 'Z',(Iontype(i),i=1,ntype)
         do kg=1,nclz
            zz=(kg-1)*dcel-tranz+zbcen
            write(outu,'(3x,f10.3,2x,8(f10.4,2x))') 
     $           zz,(Diffusion((i-1)*nclz+kg),i=1,ntype)
         enddo

         write(outu,*)
 102     format(6x,a,i3)

c.......................................................................
      elseif(wrd5.eq.'BOUND')then               ! BOUNDARY_CONCENTRATION
c            ---------------
 1003    call getlin('---> ',com,inpu,outu)
         call misc(com,strg,qstrg,wrd5)
         if(.not.check(com,'end'))then
            call getfirst(com,wrd4)
            call fatnam(iontype,ntype,wrd4,itype,.false.)
            call gtdpar(com,'ctop',Ctop(itype),0.0d0) !concentration in mol/L
            call gtdpar(com,'cbot',Cbot(itype),0.0d0) !concentration in mol/L
            goto 1003
         endif

         write(outu,*)
         write(outu,101)
     $        'Boundary concentrations in top and bottom along Z-axis'
         write(outu,'(18x,a)') 'Ctop    Cbot'
         do i=1,ntype
            write(outu,'(6x,i3,1x,a,2f8.3)')
     $           i,iontype(i),ctop(i),cbot(i)
         enddo
c
         do i=1,ntype
            Ctop(i)=Ctop(i)*(avogadro/liter)*(angstrom**3) ! [unit charge]/A^3
            Cbot(i)=Cbot(i)*(avogadro/liter)*(angstrom**3)
         enddo

         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'OPEN ')then                                 ! OPEN
c            ---------------
         call gtipar(com,'unit',iunit,-1)
         call getwrd(com,'name',fnam)
         lfnam=len_trim(fnam)
         write(outu,105) trim(fnam),iunit
 105     format(6x,'open file ',a,' as unit ',i3)
         if(check(com,'write'))then
            if(check(com,'file '))then
            open(unit=iunit,file=fnam,form='unformatted')
            else
            open(unit=iunit,file=fnam,form='formatted')
            endif
         elseif(check(com,'read'))then
            if(check(com,'file '))then
            open(unit=iunit,file=fnam,form='unformatted')
            rewind(unit=iunit)
            else
            open(unit=iunit,file=fnam,form='formatted')
            rewind(unit=iunit)
            endif
         endif
         write(outu,*)
         
c.......................................................................
      elseif(wrd5.eq.'PROTE')then                    ! PROTEIN_STRUCTURE
c            ---------------
         if(check(com,'charmm'))then
         call gtipar(com,'cunit',iunit,-1)
         call getwrd(com,'type',rtype)
         if(iunit.gt.0) then
         n12 = 1
        write(outu,102)
     $              'Reading CRD & Charges from unit ',iunit
               call CRDRDR(iunit,datom,natom,iresid,resid,aname,
     $              x,y,z,segid,isegid,cg)
            endif
            call gtipar(com,'runit',iunit,-1)
            if(iunit.gt.0) then
               write(outu,102) 
     $              'Reading CRD & Radii from unit ',iunit
               call CRDRDR(iunit,datom,natom,iresid,resid,aname,
     $              x,y,z,segid,isegid,radius)
            call gtipar(com,'punit',punit,-1)
            if (punit .gt. 0) then
            write(outu,102)
     $     'Reading Parameters from unit ',punit
            if (n12.gt.0) then
      write(outu,*) 'WARNING! You try to re-define parameters!' 
            endif
            n11 =-1.0d0
1111        call getlin('---> ',com,punit,-1)
            call misc(com,strg,qstrg,wrd5)
            if(.not.(check(com,'!').or.check(com,'end'))) then
            n11 = n11 +1
            goto 1111
            endif
         rewind (punit)
300      continue
         read(punit,'(a)') line
         if (line(1:1).eq.'!') goto 300
         do i=1,n11
         read(punit,111) resid1(i),aname1(i),pb(i),
     >   vdw(i),stern(i),wmain(i)
         enddo
111      FORMAT(a4,2x,a4,2x,4f10.5)
         do i =1,natom
         radius(i) = 0.0d0  
         cg(i) = 0.0d0  
         enddo
         write(*,*) 'Radii type is ', rtype
         do i =1,natom
         do j =1,n11
         if ((resid(i) .eq. resid1(j)) .and. 
     >   (aname(i) .eq. aname1(j))) then
         cg(i) = wmain(j)
         if (lcase(rtype).eq.'pb') then
         radius(i) = pb(j)  
         endif
         if (lcase(rtype).eq.'stern') then
         radius(i) = stern(j)  
         endif
         if (lcase(rtype) .eq. 'vdw') then
         radius(i) = vdw(j)  
         endif
         endif
         enddo
         enddo
         endif
         endif
         else
            call gtipar(com,'unit',iunit,5)
            write(outu,102) 
     $           'Reading PROTEIN_STRUCTURE from unit ',iunit
            read(iunit,*) natom
            do i=1,natom
               read(iunit,*) j,x(i),y(i),z(i),cg(i),radius(i)
            enddo
         endif
         write(outu,*)
         write(outu,'(6x,a,i5,a)') 
     $        'PROTEIN_STRUCTURE has been read for ',natom,' atoms'
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'ITERA')then                  !ITERATION PARAMETERS
c            ---------------
c maxpnp  : maximum number of iterations for both PBEQ and PNP 
c maxphi  : maximum number of iterations for PBEQ
c maxcion : maximum number of iterations for NP
c tolcion : tolerance value for PNP solver
c tolphi  : tolerance value for PBEQ solver
c lambda1 : mixing factor for potential
c lambda2 : mixing factor for concentration

         call gtipar(com,'maxpnp',maxpnp,100) 
         call gtipar(com,'maxphi',maxphi,2000)
         call gtipar(com,'maxcion',maxcion,2000)
         call gtdpar(com,'tolphi',tolphi,2.0D-6)
         call gtdpar(com,'tolcion',tolcion,2.0D-10) 
         call gtdpar(com,'lambda',Lambda1,1.0d0)
         call gtdpar(com,'lambdaphi',Lambda1,lambda1)
         call gtdpar(com,'lambdacion',Lambda2,lambda1)

         write(outu,*)
         write(outu,101) 
     $        'ITERATION PARAMETERS has been read'
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'GRID ')then                      ! GRID PARAMETERS
c            ---------------

c ncel, nclx, ncly, nclz : number of grid points
c xbcen, ybcen, zbcen : center of box
c dcel : grid spacing

         call gtipar(com,'ncel',ncel,65) 
         call gtipar(com,'nclx',nclx,ncel) 
         call gtipar(com,'ncly',ncly,ncel) 
         call gtipar(com,'nclz',nclz,ncel) 
         call gtdpar(com,'xbcen',xbcen,0.0D0) 
         call gtdpar(com,'ybcen',ybcen,0.0D0) 
         call gtdpar(com,'zbcen',zbcen,0.0D0) 
         call gtdpar(com,'dcel',dcel,0.5D0) 

c number of grid points should be odd
         if(mod(nclx,2).eq.0) nclx=nclx+1
         if(mod(ncly,2).eq.0) ncly=ncly+1
         if(mod(nclz,2).eq.0) nclz=nclz+1

         write(outu,*)
         write(outu,101) 
     $        'GRID PARAMETERS has been read'
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'PHYSI')then                  ! PHYSICAL_PARAMETERS
c            ---------------
c watr  : solvent probe radius
c ionr  : ion exclusion radius (Stern layer)
c epsw, epsp : dielectric constants of water and protein
c temp  : temperature in [K]
c conc  : concentration in [mol/L]

         call gtdpar(com,'watr',WATR,0.0D0) 
         call gtdpar(com,'ionr',IONR,0.0D0) 
         call gtdpar(com,'epsw',EPSW,80.0D0) 
         call gtdpar(com,'epsp',EPSP,1.0D0) 
         call gtdpar(com,'temp',Temp,300.0D0) 
         call gtdpar(com,'conc',CONC,0.0D0) 

         write(outu,*)
         write(outu,101) 
     $        'PHYSICAL PARAMETERS has been read'
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'MEMBR')then                              !MEMBRANE
c            ---------------

c epsm  : dielectric constants of membrane
c epsh  : dielectric constants of head groups
c tmemb : thickness of membrane including head groups
c htmemb: thickness of head groups
c zmemb : membrane position along Z-axis
c vmemb : transmembrane potentail in [volts]

         call gtdpar(com,'epsm',EPSM,1.0D0)
         call gtdpar(com,'epsh',EPSH,1.0D0) 
         call gtdpar(com,'tmemb',TMEMB,0.0D0) 
         call gtdpar(com,'zmemb',ZMEMB,0.0D0) 
         call gtdpar(com,'vmemb',VMEMB,0.0D0) 
         call gtdpar(com,'htmemb',HTMEMB,0.0D0) 
 
c     convert to [unit charge]/[Angstroms] for PBEQ1 and PNP1 subroutine
         VMEMB=VMEMB/14.40145905

         write(outu,*)
         write(outu,101) 
     $        'MEMBRANE PARAMETERS has been read'
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'OBJEC')then                               !OBJECTS
c            ---------------
c SPHERE PARAMETERS
c rsphe,xyzsphe : radius and positions of a sphere
c epss : dielectric constant of a spherical cavity

         call gtdpar(com,'rsphe',rsphe,0.0D0) 
         call gtdpar(com,'xsphe',xsphe,0.0D0) 
         call gtdpar(com,'ysphe',ysphe,0.0D0) 
         call gtdpar(com,'zsphe',zsphe,0.0D0) 
         call gtdpar(com,'epss',epss,1.0D0) 
         Qskap=check(com,'skappa')

c CYLINDER PARAMETERS
c rcyln, hcyln, xyzcyln : radius, height, and positions of a cylinder
c epsc : dielectric constant inside cylinder

         call gtdpar(com,'rcyln',rcyln,0.0D0) 
         call gtdpar(com,'hcyln',hcyln,0.0D0) 
         call gtdpar(com,'xcyln',xcyln,0.0D0) 
         call gtdpar(com,'ycyln',ycyln,0.0D0) 
         call gtdpar(com,'zcyln',zcyln,0.0D0) 
         call gtdpar(com,'epsc',epsc,1.0D0) 
         Qckap=check(com,'ckappa')

c BOX PARAMETERS
c b{x,y,z}max and b{x,y,z} : X, Y, and Z dimensions of an orthorhombic box
c epsb : dielectric constant of a box
         call gtdpar(com,'bxmax',bxmax,0.0D0) 
         call gtdpar(com,'bxmin',bxmin,0.0D0) 
         call gtdpar(com,'bymax',bymax,0.0D0) 
         call gtdpar(com,'bymin',bymin,0.0D0) 
         call gtdpar(com,'bzmax',bzmax,0.0D0) 
         call gtdpar(com,'bzmin',bzmin,0.0D0) 
         call gtdpar(com,'epsb',epsb,1.0D0) 
         Qbkap=check(com,'bkappa')

c OVERLAY CYLINDER PARAMETERS
c rocyl, hocyl, xyzocyl : radius, height, and positions of a cylinder
c epso : dielectric constant inside cylinder

         call gtdpar(com,'rocyl',rocyl,0.0D0) 
         call gtdpar(com,'hocyl',hocyl,0.0D0) 
         call gtdpar(com,'xocyl',xocyl,0.0D0) 
         call gtdpar(com,'yocyl',yocyl,0.0D0) 
         call gtdpar(com,'zocyl',zocyl,0.0D0) 
         call gtdpar(com,'epso',epso,1.0D0) 
         Qokap=check(com,'okappa')

c Elliptic CYLINDER PARAMETERS
c Ax, By, tmin, hecyln, xyzecyln : radius, height, and positions of a cylinder
c epsec : dielectric constant inside cylinder

         call gtdpar(com,'ax',ax,0.0D0)
         call gtdpar(com,'by',by,0.0D0)
         call gtdpar(com,'tmin',tmin,1.0D0)
         call gtdpar(com,'hecyln',hecyln,0.0D0)
         call gtdpar(com,'xecyln',xecyln,0.0D0)
         call gtdpar(com,'yecyln',yecyln,0.0D0)
         call gtdpar(com,'zecyln',zecyln,0.0D0)
         call gtdpar(com,'epsec',epsec,1.0D0)
c         write(*, *) 'ax', ax 
c         write(*, *) 'by', by 
c         write(*, *) 'tmin', tmin 
c         write(*, *) 'hecyln',hecyln 
c         write(*, *) 'Qeckap', Qeckap
         Qeckap=check(com,'eckap')
c         write(*, *) 'Qeckap', Qeckap

         write(outu,*)
         write(outu,101) 
     $        'OBJECTS PARAMETERS has been read'
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'PBEQ ')then                           !PBEQ SOLVER
c            ---------------
         Qiterate=check(com,'skipgridsetup')  ! skip prep1

         IF(Qiterate) THEN
            Qnonlinear =check(com,'nonlinear')  ! non-linear PBEQ solver
            Qpartlinear=check(com,'partlinear') ! partially linearized PBEQ solver
            Qunder     =check(com,'underrelax') ! under-relaxation for non-linear PB
         ELSE
! if we don't have this, we may have some problems with the setup of
! transmembrane potential (see prep.f (mayer_vmemb))
            QPNP=.false.

! prepare grid parameters, physical constants, membrane, and so on.
            CALL PREP1(com)
         ENDIF

! PBEQ solvers
         IF(Qnonlinear) THEN
            CALL PBEQ2(MAXPHI,TOLPHI,KAPPA2,Temp,Lambda1,
     $           NCLX,NCLY,NCLZ,DCEL,VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,
     $           QphiXYPBC,QphiXYZPBC,Qunder)
         ELSEIF(Qpartlinear) THEN
            CALL PBEQ3(MAXPHI,TOLPHI,KAPPA2,Temp,Lambda1,
     $           NCLX,NCLY,NCLZ,DCEL,VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,
     $           QphiXYPBC,QphiXYZPBC,Qunder)
         ELSE
            CALL PBEQ1(MAXPHI,TOLPHI,KAPPA2,KAP2TOP,KAP2BOT,
     $           NCLX,NCLY,NCLZ,DCEL,TMEMB,ZMEMB,TRANZ,ZBCEN,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,QphiXYPBC,QphiXYZPBC)
         ENDIF

c calculate electrostatic solvation energy
         CALL ENPB1(NATOM,X,Y,Z,CG,NCLX,NCLY,NCLZ,DCEL,PHI,
     $        TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,0.5D0)

         QPBEQ =.true.
         QMCDEN=.false.
         QPHI  =.false.
         write(outu,*)

!.......................................................................
      elseif(wrd5.eq.'PNP  ')then                           ! PNP SOLVER  
!            ---------------
         QPNP=.true.
         IF(QPBEQ) THEN
            QPHI=.true.
         ELSE
            write(outu,101)
            write(outu,101)
     $      'Warning: PBEQ should be solved first for PNP solutions'
         ENDIF

! prepare grid parameters, physical constants, membrane, and so on.
         CALL PREP1(com)

! 3d-PNP solver
         CALL PNP1(MAXPNP,MAXPHI,MAXCION,TOLCION,TOLPHI,
     $        Lambda1,Lambda2,
     $        TEMP,Ntype,NCLX,NCLY,NCLZ,DCEL,
     $        PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,EFFEPS,WAR,
     $        Cion,Zion,Iontype,Diffusion,
     $        QCionXYPBC,QphiXYPBC,QPHI)

c calculate electrostatic solvation energy
         CALL ENPB1(NATOM,X,Y,Z,CG,NCLX,NCLY,NCLZ,DCEL,PHI,
     $        TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,0.5D0)

         QPBEQ =.true.
         QMCDEN=.false.
         QPHI  =.false.
         write(outu,*)

!.......................................................................
      elseif(wrd5.eq.'PROFI')then                              ! PROFILE
!            ---------------
         IF(QPNP) THEN

! density, flux, or current profiles along Z-axis
            Qrho=check(com,'rho')         ! density profile alogn Z-axis
            Qflux=check(com,'flux')       ! flux profile alogn Z-axis
            Qcurrent=check(com,'current') ! current profile alogn Z-axis
            call gtdpar(com,'rsphe',rsphe,0.0D0) 
            call gtdpar(com,'xsphe',xsphe,0.0D0) 
            call gtdpar(com,'ysphe',ysphe,0.0D0) 
            call gtdpar(com,'zsphe',zsphe,0.0D0) 

            IF(Qrho.or.Qflux.or.Qcurrent) THEN
               CALL PROFILE(Ntype,Iontype,Temp,Nclx,Ncly,Nclz,dcel,
     $              tranx,trany,tranz,xbcen,ybcen,zbcen,
     $              rsphe,xsphe,ysphe,zsphe,
     $              PHI,Cion,Zion,Diffusion,Qrho,Qflux,Qcurrent)
            ENDIF

            rsphe   =0.0d0
         ELSE
            write(outu,101) 'Grid is not prepared !'
         ENDIF

         write(outu,*)

!.......................................................................
      elseif(wrd5.eq.'COUNT')then                           ! COUNTERION
!            ---------------
! counter ion profile along Z

         Qnonlinear =check(com,'nonlinear')
         Qpartlinear =check(com,'partlinear')
         call gtdpar(com,'rsphe',rsphe,0.0D0) 
         call gtdpar(com,'xsphe',xsphe,0.0D0) 
         call gtdpar(com,'ysphe',ysphe,0.0D0) 
         call gtdpar(com,'zsphe',zsphe,0.0D0) 
         call gtdpar(com,'rdist',rdist,0.0D0) 

         Qporin =check(com,'porin')

         IF(QPBEQ) THEN
            TRANX = 0.5D0*(NCLX-1)*DCEL
            TRANY = 0.5D0*(NCLY-1)*DCEL
            TRANZ = 0.5D0*(NCLZ-1)*DCEL
            IF(Qporin)then
               CALL PORIN(NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $              VMEMB,TMEMB,ZMEMB,CONC,TEMP,
     $              TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $              Qnonlinear,Qpartlinear)
            ELSEIF(rsphe.eq.0.0d0.and.rdist.eq.0.0d0) THEN
               CALL COUNTERION1(NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $              VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN,CONC,TEMP,
     $              Qnonlinear,Qpartlinear)
            ELSE
               CALL COUNTERION2(NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $              VMEMB,TMEMB,ZMEMB,CONC,TEMP,
     $              TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $              rsphe,xsphe,ysphe,zsphe,rdist,
     $              Qnonlinear,Qpartlinear)
            ENDIF
         ELSE
            write(outu,101) 'Grid is not prepared !'
         ENDIF

         rsphe=0.0d0
         write(outu,*)

!.......................................................................
      elseif(wrd5.eq.'RXNFL')then                               ! RXNFLD
!            ---------------

         QMMIJ=.true.

         CALL PREP1(com)
         CALL RXNFLD(com)

         write(outu,*)

!.......................................................................
      elseif(wrd5.eq.'RFPAR')then                              ! GBRADII
!            ---------------
         
         QRFPAR=.true.
         Qbox=check(com,'box')                  
c maximum radius (cut-off) for numerical integration          
         call gtdpar(com,'rcut',rmax,50.0D0)
c         rmax=50.0D0 
         
c prepare grid parameters, physical constants, membrane, and so on
         CALL PREP1(com) 
         
c check shape and size of the GCMC/BD simulation system        
         IF (Qbox) THEN
            volume=(bxmax-bxmin)*(bymax-bymin)*(bzmax-bzmin)
            write(outu,101)
     $      'Orthorhombic GCMC/BD simulation system'
            write(outu,101)
     $      'Size and location defined by OBJECTS BOX PARAMETERS'     
         ELSE
            volume=rsphe    
            write(outu,101)
     $      'Spherical GCMC/BD simulation system'
            write(outu,101)
     $      'Size and location defined by OBJECTS SPHERE PARAMETERS'    
         END IF 
         IF (volume.LT.1.0D-10) THEN
            write(outu,*)
            write(outu,101)
     $      'Warning: Simulation system has no volume !'
            write(outu,101) 
     $      'Reaction field parameters are not calculated'
! calculate reaction field parameters 
         ELSE          
            write(outu,*)
            write(outu,101) 'Calculating reaction field parameters'
            CALL rfpar(ntype,nclx,ncly,nclz,dcel,
     $           xbcen,ybcen,zbcen,
     $           epsx,epsy,epsz,epsw,tranx,trany,tranz,
     $           tmemb,htmemb,zmemb,epsm,epsh,     
     $           rsphe,xsphe,ysphe,zsphe,
     $           bxmax,bymax,bzmax,bxmin,bymin,bzmin,rmax,rion,Qbox,
     $           gnowat,gwater,gsrfen,greff)
         END IF            

!         CALL gridz(2,nclx,ncly,nclz,gsrfen)
!         CALL gridz(2,nclx,ncly,nclz,greff)
         
         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'WRITE')then                                ! WRITE
c            ---------------
         if(QPBEQ.or.QPNP.or.QMMIJ.or.QRFPAR) then
            CALL WRIGD0(com)
         else
            write(outu,101) 'Grid is not prepared !'
         endif

         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'READ ')then                                 ! READ
c            ---------------
         CALL REAGD0(com)

         write(outu,*)

c.......................................................................
      elseif(wrd5.eq.'STOP ')then                                 ! STOP 
c            ---------------
         write(outu,*)
         finish = cputime(start)
         IF(finish.le.60.0d0) THEN
            write(outu,106) 'cpu time   :',finish,' seconds'
         ELSEIF(finish.gt.3600.0d0) THEN
            finish=finish/3600.0d0
            write(outu,106) 'cpu time   :',finish,' hours'
         ELSE
            finish=finish/60.0d0
            write(outu,106) 'cpu time   :',finish,' minutes'
         ENDIF
         call time(hour)
         write(outu,107)    'finished at: ',hour
         stop

 106     format(6x,a,f9.2,a)
 107     format(6x,2a)

c.......................................................................
      elseif(wrd5.eq.'*END*')then
c            ---------------
         write(outu,103)
 103     format(6x,'*ERROR*  END-OF-FILE ENCOUNTERED IN UNIT',i3)
         goto 2000
      else
         write(outu,104)
 104     format(6x,'*ERROR*  Unrecognized command')

c.......................................................................
      endif

      goto 1000
c
 2000 STOP
      END

! for GFORTRAN
      subroutine time(hour)
      character*24  hour
      hour=ctime(time8())
      end subroutine
