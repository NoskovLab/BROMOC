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

      SUBROUTINE PREP1(com)
c------------------------------------------------------------------------
c     Prepare grid parameters, physical parameters, and so on.
c
      use charfuncmod
      implicit none
      include 'pnp.fcm'
      include 'consta.fcm'
      include 'mainio.fcm'

      character*(*) com

c local 
      integer*4   i
      integer*4   pncel,pnclx,pncly,pnclz
      real*8    pxbcen,pybcen,pzbcen,pdcel,ptranx,ptrany,ptranz
      real*8    factor0


c LINEAR or NON-LINEAR PBEQ solver?
c =================================
      Qnonlinear =check(com,'nonlinear')   ! non-linear PBEQ solver
      Qpartlinear=check(com,'partlinear')  ! partially linearized PBEQ solver
      Qunder     =check(com,'underrelax')  ! under-relaxation for non-linear PBEQ solver
      if(qnonlinear.or.qpartlinear.or.qunder) then
         write(outu,102)
         IF(Qnonlinear) WRITE(OUTU,102) 'NON-LINEAR PBEQ SOLVER'
         if(qpartlinear) write(outu,102) 
     $        'PARTIALLY LINEARIZED PBEQ SOLVER'
         if(qunder) write(outu,102) 
     $        'Under-relaxation will be done with fixed Lambda'
      endif


c ITERATION PARAMETERS for the SOR (Successive OverRelaxation) method
c ====================
c maxpnp  : maximum number of iterations for both PBEQ and PNP 
c maxphi  : maximum number of iterations for PBEQ
c maxcion : maximum number of iterations for NP
c tolcion : tolerance value for PNP solver
c tolphi  : tolerance value for PBEQ solver
c lambda1 : mixing factor for potential
c lambda2 : mixing factor for concentration

      write(outu,102)
      write(outu,102) 'iteration parameters'
      write(outu,104) 
     $     'Maximum # iterations for PNP             (MAXPNP) =',MAXPNP
      write(outu,104) 
     $     'Maximum # iterations for potentials      (MAXPHI) =',MAXPHI
      write(outu,104) 
     $     'Maximum # iterations for concentrations (MAXCION) =',MAXCION
      write(outu,105) 
     $     'Tolerance value for potentials           (TOLPHI) =',TOLPHI
      write(outu,105) 
     $     'Tolerance value for concentrations      (TOLCION) =',TOLCION
      write(outu,105)
     $     'Global mixing factor for potentials   (LAMBDAPHI) =',LAMBDA1
      write(outu,105)
     $     'Global mixing factor for conc.       (LAMBDACION) =',LAMBDA2

 102  format(6x,a,f8.3,a,f8.3)
 104  format(6x,a,1x,i8)
 103  format(6x,a,1x,f8.3,a,2x)
 105  format(6x,a,1x,e10.3)

c GRID PARAMETERS
c ===============
c ncel, nclx, ncly, nclz : number of grid points
c xbcen, ybcen, zbcen : center of box
c dcel : grid spacing

      TRANX = 0.5D0*(NCLX-1)*DCEL
      TRANY = 0.5D0*(NCLY-1)*DCEL
      TRANZ = 0.5D0*(NCLZ-1)*DCEL

      write(outu,102)
      WRITE(OUTU,102) 'GRID PARAMETERS & BOX SIZE'
      write(outu,'(6x,a,3i6)')
     $     'Number of grid points:',NCLX,NCLY,NCLZ
      write(outu,102)
     $     'Box in X from ',XBCEN-TRANX,' to ',XBCEN+TRANX
      write(outu,102)
     $     'Box in Y from ',YBCEN-TRANY,' to ',YBCEN+TRANY
      write(outu,102)
     $     'Box in Z from ',ZBCEN-TRANZ,' to ',ZBCEN+TRANZ


c PHYSICAL PARAMETERS
c ===================
c watr  : solvent probe radius
c ionr  : ion exclusion radius (Stern layer)
c epsw, epsp, epsm : dielectric constants of water, protein, and membrane
c epsh  : dielectric constants of head groups
c temp  : temperature in [K]
c conc  : concentration in [mol/L]
c tmemb : thickness of membrane including head groups
c htmemb: thickness of head groups
c zmemb : membrane position along Z-axis
c vmemb : transmembrane potentail in [volts]

      write(outu,102)
      WRITE(OUTU,103) 'PHYSICAL PARAMETERS'
      write(outu,103) 
     $     'Solvent probe radius               (WATR) =',WATR,' [Angs]'
      write(outu,103) 
     $     'Ion exclusion radius (Stern layer) (IONR) =',IONR,' [Angs]'
      write(outu,103) 
     $     'Solvent dielectric constant        (EPSW) =',EPSW
      write(outu,103) 
     $     'Protein dielectric constant        (EPSP) =',EPSP
      write(outu,103)
     $     'Temperature                        (TEMP) =',TEMP,' [K]'
      write(outu,103)
     $     'Salt concentration                 (CONC) =',CONC,' [mol/L]'

c                 4*pi*Coulomb**2                I
c     kappa^2 = --------------------------- * ------ , I=Sum_i qi*qi*Ci
c                 4*pi*eps0*angstrom*kboltz    Temp
      kap2bot=0.0d0
      kap2top=0.0d0
      factor0=(Coulomb*Coulomb)/eps0/angstrom/kboltz
      do i=1,ntype
         kap2bot=kap2bot+factor0*Zion(i)*Zion(i)*Cbot(i)/Temp
         kap2top=kap2top+factor0*Zion(i)*Zion(i)*Ctop(i)/Temp
      enddo
      kappa2=(kap2top+kap2bot)/2.0d0

      if(conc.gt.0.0d0) then
         kappa2=2530.362733d0*conc/temp
         if(kap2bot.eq.0.0d0) kap2bot=kappa2
         if(kap2top.eq.0.0d0) kap2top=kappa2
      endif

      if(kappa2.gt.0.0d0) write(outu,103) 
     $    'Debye screening length                    =',
     $    sqrt(epsw/kappa2),' [angs]'

      if(tmemb.gt.0.0d0)then
         write(outu,103)
         WRITE(OUTU,103) 'MEMBRANE'
         write(outu,103)
     $    'Membrane thickness along Z        (TMEMB) =',TMEMB,' [Angs]'
         write(outu,103)
     $    'Membrane position along Z         (ZMEMB) =',ZMEMB,' [Angs]'
         write(outu,103)
     $    'Membrane dielectric constant       (EPSM) =',EPSM
         if(htmemb.gt.0.0d0) then
           write(outu,103)
     $    'Membrane headgroup thickness     (HTMEMB) =',HTMEMB,' [Angs]'
           write(outu,103)
     $    'Headgroup dielectric constant      (EPSH) =',EPSH
         endif
         write(outu,'(6x,2(a,f8.3),a,2x)')
     $    'The membrane goes from                    =',
     $     ZMEMB-TMEMB/2,' to ',ZMEMB+TMEMB/2,' [Angs]'
         write(outu,103)
     $    'Transmembrane voltage along Z     (VMEMB) =',
     $     VMEMB*14.40145905d0,' [Volts]'
      endif

c SPHERE PARAMETERS
c ===================
c rsphe,xyzsphe : radius and positions of a sphere
c epss : dielectric constant of a spherical cavity

      IF(rsphe.gt.0.0d0)THEN
         WRITE(OUTU,103)
         WRITE(OUTU,103) 'SPHERE'
         WRITE(OUTU,103)
     $     'Sphere radius                     (RSPHE) =',rsphe,' [Angs]'
         WRITE(OUTU,103)
     $     'Sphere dielectric constant         (EPSS) =',epss
         WRITE(OUTU,'(6x,a,3f8.3)')
     $     'Sphere position                           =',
     $     xsphe,ysphe,zsphe
         IF(Qskap) THEN
            WRITE(OUTU,103) 
     $           'Mobile ions are accessible inside the sphere'
         ELSE
            WRITE(OUTU,103) 
     $           'Mobile ions are not accessible inside the sphere'
         ENDIF
      ENDIF

c CYLINDER PARAMETERS
c ===================
c rcyln, hcyln, xyzcyln : radius, height, and positions of a cylinder
c epsc : dielectric constant inside cylinder

      IF(rcyln.gt.0.0d0)THEN
         WRITE(OUTU,103)
         WRITE(OUTU,103) 'CYLINDER'
         WRITE(OUTU,103)
     $     'Cylinder height along Z           (HCYLN) =',HCYLN,' [Angs]'
         WRITE(OUTU,103)
     $     'Cylinder radius in XY plane       (RCYLN) =',RCYLN,' [Angs]'
         WRITE(OUTU,103)
     $     'Cylinder dielectric constant       (EPSC) =',EPSC
         WRITE(OUTU,'(6x,a,3f8.3)')
     $     'Cylinder position                         =',
     $     XCYLN,YCYLN,ZCYLN
         IF(qckap) THEN
            WRITE(OUTU,103) 
     $           'Mobile ions are accessible inside the cylinder'
         ELSE
            WRITE(OUTU,103) 
     $           'Mobile ions are not accessible inside the cylinder'
         ENDIF
      ENDIF

c ELLIPTIC CYLINDER PARAMETERS
c ===================
c ax, by, hecyln, xyzecyln : radius, height, and positions of an elliptic cylinder
c epsec : dielectric constant inside cylinder

      IF(ax.gt.0.0)THEN
         WRITE(OUTU,103)
         WRITE(OUTU,103) 'ELLIPTIC CYLINDER'
         WRITE(OUTU,103)
     $     'Elliptic cylinder height along Z       (HECYLN) ='
     $     ,HECYLN,' [Angs]'
         WRITE(OUTU,103)
     $     'Elliptic cylinder radius in X axis         (AX) ='
     $     ,ax,' [Angs]'
         WRITE(OUTU,103)
     $     'Elliptic cylinder radius in Y axis         (BY) ='
     $     ,by,' [Angs]'
         WRITE(OUTU,103)
     $     'Area scaling factor in Z axis              (tmin) ='
     $     ,tmin,' []'
         WRITE(OUTU,103)
     $     'Elliptic cylinder dielectric constant   (EPSEC) =',EPSEC
         WRITE(OUTU,'(6x,a,3f8.3)')
     $     'Elliptic cylinder position                      =',
     $     XECYLN,YECYLN,ZECYLN
         IF(Qeckap) THEN
            WRITE(OUTU,103)
     $           'Ions are accessible inside the elliptic cylinder'
         ELSE
            WRITE(OUTU,103)
     $           'Ions are not accessible inside the elliptic cylinder'
         ENDIF
      ENDIF


c BOX PARAMETERS
c ===================
c b{x,y,z}max and b{x,y,z} : X, Y, and Z dimensions of an orthorhombic box
c epsb : dielectric constant of a box

      IF(bxmax.ne.bxmin)THEN
         WRITE(OUTU,103)
         WRITE(OUTU,103) 'BOX'
         WRITE(OUTU,103)
     $     'Box dielectric constant            (EPSB) =',epsb
         WRITE(OUTU,102)
     $     'Box in X from ',bxmin,' to ',bxmax
         WRITE(OUTU,102)
     $     'Box in Y from ',bymin,' to ',bymax
         WRITE(OUTU,102)
     $     'Box in Z from ',bzmin,' to ',bzmax
         IF(Qbkap) THEN
            WRITE(OUTU,103) 
     $           'Mobile ions are accessible inside the box'
         ELSE
            WRITE(OUTU,103) 
     $           'Mobile ions are not accessible inside the box'
         ENDIF
      ENDIF

c OVERLAY CYLINDER PARAMETERS
c ===================
c rocyl, hocyl, xyzocyl : radius, height, and positions of a cylinder
c epso : dielectric constant inside cylinder

      IF(rocyl.gt.0.0d0)THEN
         WRITE(OUTU,103)
         WRITE(OUTU,103) 'OVERLAY CYLINDER'
         WRITE(OUTU,103)
     $     'Cylinder height along Z           (HCYLN) =',HOCYL,' [Angs]'
         WRITE(OUTU,103)
     $     'Cylinder radius in XY plane       (RCYLN) =',ROCYL,' [Angs]'
         WRITE(OUTU,103)
     $     'Cylinder dielectric constant       (EPSC) =',EPSO
         WRITE(OUTU,'(6x,a,3f8.3)')
     $     'Cylinder position                         =',
     $     XOCYL,YOCYL,ZOCYL
         IF(Qokap) THEN
            WRITE(OUTU,103) 
     $           'Mobile ions are accessible inside the cylinder'
         ELSE
            WRITE(OUTU,103) 
     $           'Mobile ions are not accessible inside the cylinder'
         ENDIF
      ENDIF

c MOLECULAR SURFACE
c ==================
c molecular surface = ( contact + reentrant ) surface
c MAPT : max number of points on the sphere surface (set to 10000 in pnp.fcm)
c POX,POY,POZ : random points on the surface

      QREEN =check(com,'reentrant') .or. check(com,'reen')
      IF(QREEN) THEN
         WRITE(OUTU,103)
         WRITE(OUTU,103)
     $  'Molecular (contact+reentrant) surface will be considered'
         IF(WATR.le.0.0d0) WRITE(OUTU,103)
     $  'WATR should be greater than 0.0 to consider molecular surface'
      ENDIF

c BOUNDARY CONDITIONS
c ===================
c Zero boundary potentials
      QZEROBP=check(com,'zerobp')
      IF(QZEROBP) then
         WRITE(OUTU,103) 
         WRITE(OUTU,103) 
     $     'Zero potentials will be used at the boundary points'
      ENDIF

c XY Periodic boundary condition for potentials or ionic concentrations
      QphiXYPBC =check(com,'phixypbc')
      QphiXYZPBC=check(com,'phixyzpbc')
      QCionXYPBC=check(com,'cionxypbc')

      IF(QphiXYPBC.or.QCionXYPBC.or.QphiXYZPBC) then
         WRITE(OUTU,103)
         IF(QphiXYZPBC) WRITE(OUTU,103)
     $ 'Periodic boundary conditions in XYZ will be used for potentials'
         IF(QphiXYPBC) WRITE(OUTU,103)
     $  'Periodic boundary conditions in XY will be used for potentials'
         IF(QCionXYPBC) WRITE(OUTU,'(6x,2A)')
     $   'Periodic boundary conditions in XY will be used for ',
     $   'ionic concentrations'
      ENDIF

      IF(QMMIJ) return

c Focussing method for potentials or ionic concentrations
      QphiFOCUS=check(com,'phifocus')
      QCionFOCUS=check(com,'cionfocus')

      IF(QphiFOCUS.or.QCionFOCUS) then
         WRITE(OUTU,103)
         IF(QphiFOCUS) WRITE(OUTU,103)
     $   'Focussing method will be used for potentials'
         IF(QCionFOCUS) WRITE(OUTU,103)
     $   'Focussing method will be used for ionic concentrations'

c for previous grid parameters
         call gtipar(com,'pncel',pncel,ncel) 
         call gtipar(com,'pnclx',pnclx,nclx) 
         call gtipar(com,'pncly',pncly,ncly) 
         call gtipar(com,'pnclz',pnclz,nclz) 
         call gtdpar(com,'pxbcen',pxbcen,xbcen) 
         call gtdpar(com,'pybcen',pybcen,ybcen) 
         call gtdpar(com,'pzbcen',pzbcen,zbcen) 
         call gtdpar(com,'pdcel',pdcel,dcel) 
         ptranx=0.5D0*(pnclx-1)*pdcel
         ptrany=0.5D0*(pncly-1)*pdcel
         ptranz=0.5D0*(pnclz-1)*pdcel

         WRITE(OUTU,'(6x,A,3I6)')
     $        'Number of coarse grid points:',PNCLX,PNCLY,PNCLZ
         WRITE(OUTU,102)
     $        'Coarse box in X from ',PXBCEN-PTRANX,' to ',PXBCEN+PTRANX
         WRITE(OUTU,102)
     $        'Coarse box in Y from ',PYBCEN-PTRANY,' to ',PYBCEN+PTRANY
         WRITE(OUTU,102)
     $        'Coarse box in Z from ',PZBCEN-PTRANZ,' to ',PZBCEN+PTRANZ

         IF(PXBCEN-PTRANX.GT.XBCEN-TRANX .OR.  
     $      PXBCEN+PTRANX.LT.XBCEN+TRANX .OR.  
     $      PYBCEN-PTRANY.GT.YBCEN-TRANY .OR.  
     $      PYBCEN+PTRANY.LT.YBCEN+TRANY .OR.  
     $      PZBCEN-PTRANZ.GT.ZBCEN-TRANZ .OR.  
     $      PZBCEN+PTRANZ.LT.ZBCEN+TRANZ) THEN 
            WRITE(OUTU,'(/,6X,2A)')
     $           'PREP1 Warning: Focussed system should be ',
     $           'inside previous one.'
            STOP
         ENDIF

c construct boundary potentials
         IF(QphiFOCUS) THEN
            CALL FOCUSPHI(PNCLX,PNCLY,PNCLZ,PDCEL,
     $           PTRANX,PTRANY,PTRANZ,PXBCEN,PYBCEN,PZBCEN,
     $           NCLX,NCLY,NCLZ,DCEL,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           PHI,PHIB)
         ENDIF

c construct boundary concentrations
         IF(QcionFOCUS) THEN
            CALL FOCUSCion(PNCLX,PNCLY,PNCLZ,PDCEL,
     $           PTRANX,PTRANY,PTRANZ,PXBCEN,PYBCEN,PZBCEN,
     $           NCLX,NCLY,NCLZ,DCEL,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           Cion,BCion,Ntype)
         ENDIF
      ENDIF

c PREVIOUS SPACE-DEPENDENT FUNCTIONS
c ==================================
c different MCDEN array for ions
      Qnmcden=check(com,'nmcden')  
c ion-accessible space function
      IF(.not.QMCDEN) QMCDEN =check(com,'mcden')
c use PHI for initial guess for NP equation
      IF(.not.QPHI)   QPHI   =check(com,'phi')  

      IF(QMCDEN.or.QPHI) THEN
         WRITE(OUTU,103) 
         IF(QMCDEN) THEN
            WRITE(OUTU,'(6x,2A)')
     $      'Previous mobile-ion charge density will be used for ',
     $      'MCDEN array (no modification in MAYER)'
         ENDIF
         IF(QPHI) THEN
            IF(Qnonlinear) WRITE(OUTU,'(6x,2A)')
     $    'Previous potentials will be used for the initial guess for ',
     $    'nonlinear PBEQ solver'
            IF(Qpartlinear) WRITE(OUTU,'(6x,2A)')
     $    'Previous potentials will be used for the initial guess for ',
     $    'partially linearized PBEQ solver'
            IF(QPNP) WRITE(OUTU,'(6x,2A)')
     $    'Previous potentials will be used for the initial guess for ',
     $    'Nernst-Planck solver in PNP solver'
         ENDIF
         IF( (QMCDEN.or.QPHI) .and. .not. QPBEQ ) THEN
            WRITE(OUTU,'(6x,2A)')
     $    'There is no previous potentials or mobile-ion charge density'
            stop
         ENDIF
      ENDIF


c CONSTRUCT ALL SPACE-DEPENDENT FUNCTIONS : PHI,Cion,FCDEN,MCDEN,EPSX,EPSY,EPSZ
c =======================================

c      CALL MAYER(NATOM,X,Y,Z,CG,RADIUS,EPSW,EPSP,KAPPA2,Temp,
c     $     WATR,IONR,NCLX,NCLY,NCLZ,DCEL,
c     $     PHI,FCDEN,MCDEN,EPSX,EPSY,EPSZ,
c     $     Ntype,Cion,Zion,CTOP,CBOT,kap2top,kap2bot,
c     $     VMEMB,TMEMB,ZMEMB,EPSM,HTMEMB,EPSH,
c     $     epso,xocyl,yocyl,zocyl,rocyl,hocyl,Qokap,
c     $     epsb,bxmax,bxmin,bymax,bymin,bzmax,bzmin,Qbkap,
c     $     epsc,xcyln,ycyln,zcyln,rcyln,hcyln,Qckap,
c     $     epss,xsphe,ysphe,zsphe,rsphe,Qskap,
c     $     TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
c     $     QREEN,MAPT,POX,POY,POZ,LISTR,FISTR,
c     $     QPNP,QphiXYPBC,QphiXYZPBC,QPHI,QMCDEN,Qnmcden,
c     $     QZEROBP,QphiFOCUS,PHIB,QCionFOCUS,BCion,
c     $     Qnonlinear,Qpartlinear)

      call mayer(natom,x,y,z,cg,radius,epsw,epsp,kappa2,temp,
     $           watr,ionr,nclx,ncly,nclz,dcel,
     $           phi,fcden,mcden,epsx,epsy,epsz,
     $           ntype,cion,zion,ctop,cbot,kap2top,kap2bot,
     $           vmemb,tmemb,zmemb,epsm,htmemb,epsh,
     $           epso,xocyl,yocyl,zocyl,rocyl,hocyl,qokap,
     $           epsb,bxmax,bxmin,bymax,bymin,bzmax,bzmin,qbkap,
     $           epsc,xcyln,ycyln,zcyln,rcyln,hcyln,qckap,
     $           epsec,xecyln,yecyln,zecyln,ax,by,tmin,hecyln,qeckap,
     $           epss,xsphe,ysphe,zsphe,rsphe,qskap,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           qreen,mapt,pox,poy,poz,listr,fistr,
     $           qpnp,qphixypbc,qphixyzpbc,qphi,qmcden,qnmcden,
     $           qzerobp,qphifocus,phib,qcionfocus,bcion,
     $           qnonlinear,qpartlinear)
      write(outu,103) 
c
      return
      end

