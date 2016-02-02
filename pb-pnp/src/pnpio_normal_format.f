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

      SUBROUTINE WRIGD0(com)
c-----------------------------------------------------------------------
c     write various grid parameters
c
      implicit none
      include 'pnp.fcm'
      include 'mmij.fcm'
      include 'consta.fcm'
      include 'mainio.fcm'

      character*(*) com
      character*4 wrd4
      real*4  PRFTR
      integer*4 iunit,itype
c
      call gtipar(com,'unit',iunit,outu)
      WRITE(OUTU,'(6X,A,I3)') 'Writing to unit ',iunit

      PRFTR=1.0d0

 100  FORMAT(6X,A)

      IF(check(com,'phi'))THEN
         IF(check(com,'kcal'))THEN
            PRFTR=CELEC
            WRITE(OUTU,100)
     $           'Electrostatic potential in [KCAL/MOL]/[UNIT CHARGE]'
         ELSEIF (check(com,'volts'))THEN
            PRFTR=14.40145905D0
            WRITE(OUTU,100) 'Electrostatic potential in [VOLTS]'
         ELSE
            PRFTR=1.0d0
            WRITE(OUTU,100)
     $           'Electrostatic potential in [UNIT CHARGE]/[ANGS]'
         ENDIF
         IF(check(com,'card'))THEN
            CALL WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  PHI,PRFTR)
         ELSE
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,PHI,PRFTR,1,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'kappa2').or.check(com,'mcden').or.
     $        check(com,'access'))THEN
         PRFTR=1.0d0/(DCEL*DCEL)
         WRITE(OUTU,100) 'Debye screening factor KAPPA2 '
         IF(check(com,'card'))THEN
            IF(check(com,'pdb')) THEN
               CALL WRIGDPDB(IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $              TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,MCDEN)
            ELSE
               CALL WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $              TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $              MCDEN,PRFTR)
            ENDIF
         ELSE
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,MCDEN,PRFTR,1,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'charge'))THEN
         PRFTR=DCEL/(2.0d0*TWOPI)
         WRITE(OUTU,100) 'charge distributions on the grid '
         IF(check(com,'card'))THEN
            CALL WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  FCDEN,PRFTR)
         ELSE
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,FCDEN,PRFTR,1,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'epsx'))THEN
         WRITE(OUTU,100) 'X set of dielectric function'
         IF(check(com,'card'))THEN
            CALL WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  EPSX,PRFTR)
         ELSE
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,EPSX,PRFTR,1,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'epsy'))THEN
         WRITE(OUTU,100) 'Y set of dielectric function'
         IF(check(com,'card'))THEN
            CALL WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  EPSY,PRFTR)
         ELSE
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,EPSY,PRFTR,1,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'epsz'))THEN
         WRITE(OUTU,100) 'Z set of dielectric function'
         IF(check(com,'card'))THEN
            CALL WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  EPSZ,PRFTR)
         ELSE
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,EPSZ,PRFTR,1,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'cion'))THEN
         call getfirst(com,wrd4)
         call fatnam(iontype,ntype,wrd4,itype,.false.)
         WRITE(OUTU,'(6x,a,1x,a4,a)')
     $        'Concentration in [unit charge]/[Angs^3] for',
     $        iontype(itype),' ion'
         IF(check(com,'card'))THEN
            CALL WRIGD2(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  Cion,itype)
         ELSE
            PRFTR = 1.0d0
            CALL SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,Cion,PRFTR,itype,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
         ENDIF

      ELSEIF (check(com,'flux'))THEN
         PRFTR = celec  / ( kboltz * Temp / kcalmol)
         call getfirst(com,wrd4)
         call fatnam(iontype,ntype,wrd4,itype,.false.)
         WRITE(OUTU,'(6x,a,1x,a4,a)')
     $        'Fluxes (Jx,Jy,Jz) in [unit charge]/[ps][Angs^2] for',
     $        iontype(itype),' ion'
         IF(check(com,'card'))THEN
            CALL WRIGD3(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $                  TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $                  PHI,Cion,Zion,Diffusion,PRFTR,itype)
         ENDIF

      ELSEIF (check(com,'cg'))THEN
         IF(check(com,'card'))THEN
            IF(QPNP) THEN
               CALL WRITECG1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,Cion,Zion,
     $              TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,Ntype)
            ELSE
               CALL WRITECG2(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $              TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,CONC,TEMP)
            ENDIF
         ENDIF

      ELSEIF (check(com,'mmij'))THEN
         WRITE(OUTU,100) 'The MMIJ matrix'
         CALL SAVGD2(IUNIT,XNPOL,YNPOL,ZNPOL,'RECTBOX ',
     $        LSTPX,LSTPY,LSTPZ,MMIJ)

      ELSEIF (check(com,'rfpar'))THEN
         call getfirst(com,wrd4)
         call fatnam(iontype,ntype,wrd4,itype,.false.)
         WRITE(OUTU,'(6x,a,1x,a4,a)')
     $        'Reaction field parameter for',
     $        iontype(itype),' ion'
         IF(check(com,'card'))THEN
            WRITE(OUTU,'(6x,a)') 
     $         'Formatted output is not available for rfpar!'        
            WRITE(OUTU,'(6x,a)') 
     $         'Writing binary file!'        
         ENDIF
         CALL SAVGD3(IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $      XBCEN,YBCEN,ZBCEN,ITYPE,RION(ITYPE),IONTYPE(ITYPE),
     $      GSRFEN,GREFF)

      ELSEIF (check(com,'help'))THEN
         WRITE(OUTU,100)
     $   'Syntax: WRITE <quantity> <fact> <ion> [{CARD}{FILE} <range>]'
         WRITE(OUTU,100) 'where <quantity> can be: '
         WRITE(OUTU,100) 'phi    - electrostatic potential'
         WRITE(OUTU,100) 'kappa2 - Debye screening factor'
         WRITE(OUTU,100) 'charge - charges on the lattice'
         WRITE(OUTU,100) 'epsx   - X sets of dielectric function'
         WRITE(OUTU,100) 'epsy   - Y sets of dielectric function'
         WRITE(OUTU,100) 'epsz   - Z sets of dielectric function'
         WRITE(OUTU,100) 'cion   - concentration of <iontype>'
         WRITE(OUTU,100) 'flux   - fluexes (Jx, Jy, Jz)'
         WRITE(OUTU,100) 'cg     - ionic charges'
         WRITE(OUTU,100) 'mmij   - reaction field matrix MIJ'
         WRITE(OUTU,100) 'rfpar  - reaction field parameter'
         WRITE(OUTU,100) 
     $      'where <ion> for phi, flux, and rfpar can be: iontype'
         WRITE(OUTU,100) 'where <fact> for phi can be:'
         WRITE(OUTU,100) 'kcal   - write phi in [Kcal/mol][unit charge]'
         WRITE(OUTU,100) 'volts  - write phi in [Volts]'
         WRITE(OUTU,100) '       - write phi in [unit charge]/[Angs]'
         WRITE(OUTU,100) 'where <range> can be:'
         WRITE(OUTU,100) 'xfirst [REAL] xlast [REAL] '
         WRITE(OUTU,100) 'yfirst [REAL] ylast [REAL] '
         WRITE(OUTU,100) 'zfirst [REAL] zlast [REAL] '
         WRITE(OUTU,100)
      ELSE
         WRITE(OUTU,100) 'Set the quantity to write (try write help)'

      ENDIF

      RETURN
      END


      SUBROUTINE WRIGD1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           SMTHNG,PRFTR)
c-----------------------------------------------------------------------
c
      implicit none
      include 'consta.fcm'
      character*(*) com
      real*4 smthng(*),prftr
      real*8 dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      integer*4 iunit,nclx,ncly,nclz
      integer*4 NC3,IXS,IXF,IYS,IYF,IZS,IZF,I,J,K,L
      REAL*8  XFIR,XLAS,YFIR,YLAS,ZFIR,ZLAS,XTEM,YTEM,ZTEM

      call gtdpar(com,'xfirst',xfir,0.0d0)
      call gtdpar(com,'yfirst',yfir,0.0d0)
      call gtdpar(com,'zfirst',zfir,0.0d0)
      call gtdpar(com,'xlast',xlas,0.0d0)
      call gtdpar(com,'ylast',ylas,0.0d0)
      call gtdpar(com,'zlast',zlas,0.0d0)

      nc3=nclx*ncly*nclz

      xfir   = xfir+tranx-xbcen
      ixs=int(xfir/dcel+rsmall)+1
      if(ixs.le.0)ixs=1
      if(ixs.gt.nclx)ixs=nclx

      yfir   = yfir+trany-ybcen
      iys=int(yfir/dcel+rsmall)+1
      if(iys.le.0)iys=1
      if(iys.gt.ncly)iys=ncly

      zfir   = zfir+tranz-zbcen
      izs=int(zfir/dcel+rsmall)+1
      if(izs.le.0)izs=1
      if(izs.gt.nclz)izs=nclz

      xlas   = xlas+tranx-xbcen
      ixf=int(xlas/dcel+rsmall)+1
      if(ixf.le.0)ixf=1
      if(ixf.gt.nclx)ixf=nclx

      ylas   = ylas+trany-ybcen
      iyf=int(ylas/dcel+rsmall)+1
      if(iyf.le.0)iyf=1
      if(iyf.gt.ncly)iyf=ncly

      zlas   = zlas+tranz-zbcen
      izf=int(zlas/dcel+rsmall)+1
      if(izf.le.0)izf=1
      if(izf.gt.nclz)izf=nclz

      do i=ixs,ixf
         xtem=(i-1)*dcel-tranx+xbcen
         do j=iys,iyf
            ytem=(j-1)*dcel-trany+ybcen
            do k=izs,izf
               l=(i-1)*ncly*nclz+(j-1)*nclz+k
               ztem=(k-1)*dcel-tranz+zbcen
               write(IUNIT,20) xtem,ytem,ztem,SMTHNG(l)*PRFTR
            enddo
         enddo
      enddo
 20   format(3f10.5,2x,e16.7)

      RETURN
      END


      SUBROUTINE WRIGD2(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           Cion,itype)
c-----------------------------------------------------------------------
c
      implicit none
      include 'consta.fcm'
      character*(*) com
      real*4  cion(*)
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      integer*4 iunit,nclx,ncly,nclz,itype
c local
      integer*4 NC3,IXS,IXF,IYS,IYF,IZS,IZF,I,J,K
      integer*4 NCyz,IP0X,IP0Y,IP0,IPI
      REAL*8  XFIR,XLAS,YFIR,YLAS,ZFIR,ZLAS,XTEM,YTEM,ZTEM

      call gtdpar(com,'xfirst',xfir,0.0d0)
      call gtdpar(com,'yfirst',yfir,0.0d0)
      call gtdpar(com,'zfirst',zfir,0.0d0)
      call gtdpar(com,'xlast',xlas,0.0d0)
      call gtdpar(com,'ylast',ylas,0.0d0)
      call gtdpar(com,'zlast',zlas,0.0d0)

      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz

      xfir   = xfir+tranx-xbcen
      ixs=int(xfir/dcel+rsmall)+1
      if(ixs.le.0)ixs=1
      if(ixs.gt.nclx)ixs=nclx

      yfir   = yfir+trany-ybcen
      iys=int(yfir/dcel+rsmall)+1
      if(iys.le.0)iys=1
      if(iys.gt.ncly)iys=ncly

      zfir   = zfir+tranz-zbcen
      izs=int(zfir/dcel+rsmall)+1
      if(izs.le.0)izs=1
      if(izs.gt.nclz)izs=nclz

      xlas   = xlas+tranx-xbcen
      ixf=int(xlas/dcel+rsmall)+1
      if(ixf.le.0)ixf=1
      if(ixf.gt.nclx)ixf=nclx

      ylas   = ylas+trany-ybcen
      iyf=int(ylas/dcel+rsmall)+1
      if(iyf.le.0)iyf=1
      if(iyf.gt.ncly)iyf=ncly

      zlas   = zlas+tranz-zbcen
      izf=int(zlas/dcel+rsmall)+1
      if(izf.le.0)izf=1
      if(izf.gt.nclz)izf=nclz

      DO I=IXS,IXF
         XTEM=(I-1)*DCEL-TRANX+XBCEN
         IP0X=(I-1)*NCyz
         DO J=IYS,IYF
            YTEM=(J-1)*DCEL-TRANY+YBCEN
            IP0Y=(J-1)*NCLz+IP0X
            DO K=IZS,IZF
               ZTEM=(K-1)*DCEL-TRANZ+ZBCEN
               IP0=IP0Y+K
               IPI=IP0+(itype-1)*NC3

               WRITE(IUNIT,20) XTEM,YTEM,ZTEM,Cion(IPI)

            ENDDO
         ENDDO
      ENDDO
 20   FORMAT(3f10.5,2x,e12.4)

      RETURN
      END


      SUBROUTINE WRIGD3(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,
     $           PHI,Cion,Zion,Diffusion,PRFTR,itype)
c-----------------------------------------------------------------------
c
      implicit none
      include 'consta.fcm'
      character*(*) com
      real*4  phi(*),cion(*),prftr
      real*8  zion(*),diffusion(*)
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      integer*4 iunit,nclx,ncly,nclz,itype
c local
      integer*4 NC3,IXS,IXF,IYS,IYF,IZS,IZF,I,J,K
      integer*4 NCyz,IP0X,IP0Y,IP0,IPI,IPD
      REAL*8  XFIR,XLAS,YFIR,YLAS,ZFIR,ZLAS,XTEM,YTEM,ZTEM
      REAL*4  PHI0,PHI2,PHI4,PHI6,avePHI,C0,C2,C4,C6
      REAL*8  FACTOR3,Jx,Jy,Jz
      

      call gtdpar(com,'xfirst',xfir,0.0d0)
      call gtdpar(com,'yfirst',yfir,0.0d0)
      call gtdpar(com,'zfirst',zfir,0.0d0)
      call gtdpar(com,'xlast',xlas,0.0d0)
      call gtdpar(com,'ylast',ylas,0.0d0)
      call gtdpar(com,'zlast',zlas,0.0d0)

      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz

      xfir   = xfir+tranx-xbcen
      ixs=int(xfir/dcel+rsmall)+1
      if(ixs.le.0)ixs=1
      if(ixs.gt.nclx)ixs=nclx

      yfir   = yfir+trany-ybcen
      iys=int(yfir/dcel+rsmall)+1
      if(iys.le.0)iys=1
      if(iys.gt.ncly)iys=ncly

      zfir   = zfir+tranz-zbcen
      izs=int(zfir/dcel+rsmall)+1
      if(izs.le.0)izs=1
      if(izs.gt.nclz)izs=nclz

      xlas   = xlas+tranx-xbcen
      ixf=int(xlas/dcel+rsmall)+1
      if(ixf.le.0)ixf=1
      if(ixf.gt.nclx)ixf=nclx

      ylas   = ylas+trany-ybcen
      iyf=int(ylas/dcel+rsmall)+1
      if(iyf.le.0)iyf=1
      if(iyf.gt.ncly)iyf=ncly

      zlas   = zlas+tranz-zbcen
      izf=int(zlas/dcel+rsmall)+1
      if(izf.le.0)izf=1
      if(izf.gt.nclz)izf=nclz

      FACTOR3=Zion(itype)*PRFTR

      DO I=IXS,IXF
         XTEM=(I-1)*DCEL-TRANX+XBCEN
         IP0X=(I-1)*NCyz
         DO J=IYS,IYF
            YTEM=(J-1)*DCEL-TRANY+YBCEN
            IP0Y=(J-1)*NCLz+IP0X
            DO K=IZS,IZF
               ZTEM=(K-1)*DCEL-TRANZ+ZBCEN
               IP0=IP0Y+K
               IPI=IP0+(itype-1)*NC3
               IPD=K+(itype-1)*NCLz

               Jx=0.0D0
               Jy=0.0D0
               Jz=0.0D0

               C0=Cion(IPI)
               PHI0=PHI(IP0)*FACTOR3

               IF(I.EQ.NCLx) GOTO 111
               IF(J.EQ.NCLy) GOTO 111
               IF(K.EQ.NCLz) GOTO 111
               IF(C0.EQ.0.0D0) GOTO 110
               C0=C0*EXP(PHI0)             ! conc. -> effective conc.

               C2=Cion(IPI+NCyz)
               IF(C2.NE.0.0D0) THEN
                  PHI2=PHI(IP0+NCyz)*FACTOR3
                  C2=C2*EXP(PHI2)          ! conc. -> effective conc.
                  avePHI=(PHI0+PHI2)*0.5D0
                  Jx=-Diffusion(ipd)*EXP(-avePHI)*(C2-C0)/DCEL
               ENDIF
               C4=Cion(IPI+NCLz)
               IF(C4.NE.0.0D0) THEN
                  PHI4=PHI(IP0+NCLz)*FACTOR3
                  C4=C4*EXP(PHI4)          ! conc. -> effective conc.
                  avePHI=(PHI0+PHI4)*0.5D0
                  Jy=-Diffusion(ipd)*EXP(-avePHI)*(C4-C0)/DCEL
               ENDIF
               C6=Cion(IPI+1)
               IF(C6.NE.0.0D0) THEN
                  PHI6=PHI(IP0+1)*FACTOR3
                  C6=C6*EXP(PHI6)          ! conc. -> effective conc.
                  avePHI=(PHI0+PHI6)*0.5D0
                  Jz=-Diffusion(ipd)*EXP(-avePHI)*(C6-C0)/DCEL
               ENDIF

 110           WRITE(IUNIT,20) XTEM,YTEM,ZTEM,Jx,Jy,Jz

 111        ENDDO
         ENDDO
      ENDDO
 20   FORMAT(3f10.5,2x,3e12.4)

      RETURN
      END


      SUBROUTINE SAVGD1(IUNIT,NCLX,NCLY,NCLZ,DCEL,SMTHNG,PRFTR,itype,
     $           XBCEN,YBCEN,ZBCEN,EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM)
C-----------------------------------------------------------------------
c
      implicit none
      integer*4 IUNIT,NCLX,NCLY,NCLZ,ITYPE
      REAL*4  SMTHNG(*),PRFTR
      integer*4 I,NC3,IFIR,ILAS
      REAL*8  DCEL,XBCEN,YBCEN,ZBCEN
      REAL*8  EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
c
      NC3=NCLx*NCLy*NCLz
      IFIR=(itype-1)*NC3+1
      ILAS=(itype-1)*NC3+NC3

      IF(IUNIT.GT.0)THEN
         WRITE(IUNIT) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
         WRITE(IUNIT) EPSW,EPSP,CONC,TMEMB,ZMEMB,EPSM
         WRITE(IUNIT) (PRFTR*SMTHNG(I),I=IFIR,ILAS)
      ENDIF
c
      RETURN
      END


      SUBROUTINE SAVGD2(UNIT,POL1,POL2,POL3,SHAPE,
     $           LSTP1,LSTP2,LSTP3,MIJ)
C-----------------------------------------------------------------------
      implicit none
      include 'mainio.fcm'

      REAL*8  MIJ(*)
      integer*4 UNIT,POL1,POL2,POL3
      integer*4 LSTP1(*),LSTP2(*),LSTP3(*)
      CHARACTER*8  SHAPE
c local
      integer*4 NTPOL,NTPOL2,I,J
c
      NTPOL =POL1*POL2*POL3
      NTPOL2=NTPOL*NTPOL
      WRITE(OUTU,'(6x,2A)') 
     $ 'FIRST RECORD IS: XNPOL,YNPOL,ZNPOL,LSTPX,LSTPY,LSTPZ,MIJ FOR ',
     $  SHAPE     
      WRITE(UNIT) SHAPE
      WRITE(UNIT) POL1,POL2,POL3
      WRITE(UNIT) (LSTP1(I),I=1,NTPOL)
      WRITE(UNIT) (LSTP2(I),I=1,NTPOL)
      WRITE(UNIT) (LSTP3(I),I=1,NTPOL)
      WRITE(UNIT) (MIJ(I),I=1,NTPOL2)

c      write(6,'(10f10.5)') ((MIJ((I-1)*NTPOL+J),J=1,10),I=1,10)
c      write(6,*)
c
      RETURN
      END

      
      SUBROUTINE SAVGD3(IUNIT,NCLX,NCLY,NCLZ,DCEL,
     $           XBCEN,YBCEN,ZBCEN,TYPION,RADION,NAMION,
     $           SMTHN1,SMTHN2)
c-----------------------------------------------------------------------
c     
      IMPLICIT NONE
      
      integer*4 IUNIT,NCLX,NCLY,NCLZ,TYPION
      integer*4 I,NC3,IFIR,ILAS
      REAL*4  SMTHN1(*),SMTHN2(*)
      REAL*8  DCEL,XBCEN,YBCEN,ZBCEN,RADION     
      CHARACTER*4 NAMION 
c
      NC3=NCLX*NCLY*NCLZ
      IFIR=(TYPION-1)*NC3+1
      ILAS=TYPION*NC3
c
      IF(IUNIT.GT.0)THEN
         WRITE(IUNIT) NCLX,NCLY,NCLZ,DCEL,XBCEN,YBCEN,ZBCEN
         WRITE(IUNIT) RADION,RADION,RADION,RADION,RADION,RADION
         WRITE(IUNIT) (SMTHN1(I),I=IFIR,ILAS)
         WRITE(IUNIT) (SMTHN2(I),I=IFIR,ILAS)        
      ENDIF      
c
      RETURN
      END

            
      SUBROUTINE WRIGDPDB(UNIT,NCLX,NCLY,NCLZ,DCEL,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,SMTHNG)
c-----------------------------------------------------------------------
c
      implicit none
      integer*4 UNIT,NCLX,NCLY,NCLZ
      REAL*8  DCEL,TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN
      REAL*4  SMTHNG(*)
c local
      integer*4  il,i,j,k,ncyz,cnt
      real*8   xi,yi,zi
c
      cnt=0
      ncyz=ncly*nclz
      do il=1,nclx*ncyz
         if(smthng(il).ne.0.0d0) then
            cnt=cnt+1
            i=int((il-1)/ncyz)+1
            j=int(mod((il-1),ncyz)/nclz)+1
            k=mod(mod((il-1),ncyz),nclz)+1
            xi=dcel*(i-1)-tranx+xbcen
            yi=dcel*(j-1)-trany+ybcen
            zi=dcel*(k-1)-tranz+zbcen
            write(UNIT,101)
     $           'ATOM',CNT,' POL POL     1',
     $           XI*10d0,YI*10d0,ZI*10d0,SMTHNG(IL),'  1.00      LPOL'
         ENDIF
      enddo
 101  format(a,2x,i5,1x,a,4x,3f8.3,f6.2,a)
c
      RETURN
      END


      SUBROUTINE REAGD0(com)
c-----------------------------------------------------------------------
c     read various grid parameters
c
      implicit none
      include 'pnp.fcm'
      include 'consta.fcm'
      include 'mainio.fcm'

      character*(*) com
      character*4 wrd4
      integer*4 iunit,itype,nc3,ifir,ilas,i,ip0
      integer*4 RNCLX,RNCLY,RNCLZ
      real*8  RDCEL,RXBCEN,RYBCEN,RZBCEN
      real*8  REPSW,REPSP,RCONC,RTMEMB,RZMEMB,REPSM
      real*8  RTRANX,RTRANY,RTRANZ
c
      call gtipar(com,'unit',iunit,-1)
      WRITE(OUTU,'(6X,A,I3)') 'Reading from unit ',iunit

 100  FORMAT(6X,A)


      READ(IUNIT) RNCLX,RNCLY,RNCLZ,RDCEL,RXBCEN,RYBCEN,RZBCEN
      READ(IUNIT) REPSW,REPSP,RCONC,RTMEMB,RZMEMB,REPSM
      WRITE(OUTU,101)
      WRITE(OUTU,101) 'Number of grid points in X   (NCLX)= ',RNCLX 
      WRITE(OUTU,101) 'Number of grid points in Y   (NCLY)= ',RNCLY 
      WRITE(OUTU,101) 'Number of grid points in Z   (NCLZ)= ',RNCLZ 
      WRITE(OUTU,102) 'Center of box in X          (XBCEN)= ',RXBCEN
      WRITE(OUTU,102) 'Center of box in Y          (YBCEN)= ',RYBCEN
      WRITE(OUTU,102) 'Center of box in Z          (ZBCEN)= ',RZBCEN
      WRITE(OUTU,102) 'Grid spacing                (DCEL) = ',RDCEL
      WRITE(OUTU,102) 
      WRITE(OUTU,102) 'Solvent dielectric constant (EPSW) = ',REPSW
      WRITE(OUTU,102) 'Protein dielectric constant (EPSP) = ',REPSP
      WRITE(OUTU,102) 'Salt concentration          (CONC) = ',RCONC
      IF(RTMEMB.GT.0.0d0)THEN
         WRITE(OUTU,102)
         WRITE(OUTU,102)
     $        'Membrane thickness along Z  (TMEMB)= ',RTMEMB
         WRITE(OUTU,102)
     $        'Membrane position along Z   (ZMEMB)= ',RZMEMB
         WRITE(OUTU,102)
     $        'Membrane dielectric constant(EPSM) = ',REPSM
      ENDIF
      RTRANX = 0.5D0*(RNCLX-1)*RDCEL
      RTRANY = 0.5D0*(RNCLY-1)*RDCEL
      RTRANZ = 0.5D0*(RNCLZ-1)*RDCEL
      WRITE(OUTU,102)
      WRITE(OUTU,102)
     $     'Box in X from ',RXBCEN-RTRANX,' to ',RXBCEN+RTRANX
      WRITE(OUTU,102)
     $     'Box in Y from ',RYBCEN-RTRANY,' to ',RYBCEN+RTRANY
      WRITE(OUTU,102)
     $     'Box in Z from ',RZBCEN-RTRANZ,' to ',RZBCEN+RTRANZ
 101  FORMAT(6X,A,I6,A)
 102  FORMAT(6X,A,F8.3,A,F8.3)


      NC3=RNCLx*RNCLy*RNCLz
      IF(check(com,'phi'))THEN
         READ(IUNIT) (PHI(I),I=1,NC3)
         QPBEQ=.true.
         QPHI =.true. 

      ELSEIF(check(com,'mcden').or.check(com,'kappa2').or.
     ,       check(com,'access'))THEN
         call getfirst(com,wrd4)
         call fatnam(iontype,ntype,wrd4,itype,.true.)
         IFIR=(itype-1)*NC3+1
         ILAS=(itype-1)*NC3+NC3
         READ(IUNIT) (MCDEN(I),I=IFIR,ILAS)
         DO I=1,NC3
            IP0=(itype-1)*NC3+I
            IF(MCDEN(IP0).ne.0.0d0) MCDEN(IP0)=1.0d0
         ENDDO
         QPBEQ =.true.
         QMCDEN=.true.

      ELSEIF (check(com,'cion'))THEN
         call getfirst(com,wrd4)
         call fatnam(iontype,ntype,wrd4,itype,.false.)
         IFIR=(itype-1)*NC3+1
         ILAS=(itype-1)*NC3+NC3
         READ(IUNIT) (Cion(I),I=IFIR,ILAS)
         QPNP =.true.

      ENDIF

      RETURN
      END


      SUBROUTINE CRDRDR(iunit,datom,
     $           natom,iresid,resid,aname,x,y,z,segid,isegid,wmain)
c----------------------------------------------------------------------
c     CHARMM CRDfile READER
c
      implicit none
      include 'mainio.fcm'
      integer*4       iunit,datom,natom
      integer*4       iresid(*),isegid(*)
      real*8        wmain(*),x(*),y(*),z(*)
      character*4   resid(*),segid(*),aname(*)
c local
      character*128 line
      integer*4       anum,i

 300  continue
      read(iunit,'(a)') line
      if (line(1:1).eq.'*') goto 300

      read(line,'(i5)') natom
      if(natom.gt.datom) then
         write(outu,'(6x,a)') 
     $        'Number of atoms exceeds DATOM in pnp.fcm.'
         write(outu,'(6x,a,i10,a,i10)') 
     $        'DATOM :',datom,'   Number of Atoms :',natom
         stop
      endif

      do i=1,natom
         read(iunit,101) anum,iresid(i),resid(i),aname(i),
     $        x(i),y(i),z(i),segid(i),isegid(i),wmain(i)
      enddo
 101    format(i5,1x,i4,1x,a4,1x,a4,3f10.5,1x,a4,1x,i4,f10.5)
c
      return
      end


      SUBROUTINE WRITECG2(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,PHI,MCDEN,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,CONC,TEMP)
c-----------------------------------------------------------------------
c This subroutine computes the number of counter ions from PB
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'
      integer*4 NCLX,NCLY,NCLZ,IUNIT
      character*(*) com
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen,conc,temp
      real*4  phi(*),mcden(*)
c local
      real*8  factor1
      real*8  ioncg,bulk_rho,sumioncg,phif
      real*8  dcel2,dcel3,xx,yy,zz,zc
      real*8  zmin,zmax
      integer*4 ip0,nc3,ncyz,ig,jg,kg,iii,anum
c
      call gtdpar(com,'zmin',zmin,-100000.0d0)
      call gtdpar(com,'zmax',zmax,+100000.0d0)

      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz
      dcel2=dcel*dcel
      dcel3=dcel2*dcel

c conversion factoer from 1/(kcal/(mol*e)) to 1/(e/A)
      factor1 = celec  / ( kboltz * Temp / kcalmol)

c bulk density and number of counter ions
      bulk_rho=conc*(avogadro/liter)*(angstrom**3) ! [unit charge]/A^3

      anum=0
      do 110 kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(zz.gt.zmax.or.zz.lt.zmin) goto 110
         sumioncg=0.0D0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do 111 jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen
               if(xx.le.22.25D0.and.xx.ge.-22.25D0.and.
     $            yy.le.22.25D0.and.yy.ge.-22.25D0) goto 111
               if(mcden(ip0).ne.0.0d0) then
                  anum=anum+1
                  phif=phi(ip0)*factor1
                  ioncg=bulk_rho*dcel3*(exp(-phif)-exp(+phif))
                  sumioncg=sumioncg+ioncg
                  if(ioncg.ge.0.0D0) then
                     write(iunit,101) anum,anum,'POT ','POT ',
     $                    xx,yy,zz,'POT ',anum,ioncg
                  else
                     write(iunit,101) anum,anum,'CLA ','CLA ',
     $                    xx,yy,zz,'CLA ',anum,ioncg
                  endif
               endif
 111        enddo
         enddo
         write(*,*) zz,sumioncg
 110  enddo
 101  format(i5,1x,i4,1x,a4,1x,a4,3f10.5,1x,a4,1x,i4,f10.5)
c
      RETURN
      END


      SUBROUTINE WRITECG1(COM,IUNIT,NCLX,NCLY,NCLZ,DCEL,CION,Zion,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,Ntype)
c-----------------------------------------------------------------------
c This subroutine computes the number of counter ions from PNP
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'
      integer*4 NCLX,NCLY,NCLZ,IUNIT,Ntype
      character*(*) com
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen,zion(*)
      real*4  cion(*)
c local
      real*8  ioncg,sumioncg
      real*8  dcel2,dcel3,xx,yy,zz,zc
      real*8  zmin,zmax
      integer*4 ip0,nc3,ncyz,ig,jg,kg,i,iii,anum,ipi
c
      call gtdpar(com,'zmin',zmin,-100000.0d0)
      call gtdpar(com,'zmax',zmax,+100000.0d0)

      ncyz=ncly*nclz
      nc3=nclx*ncly*nclz
      dcel2=dcel*dcel
      dcel3=dcel2*dcel

      do 110 kg=1,nclz
         zz=(kg-1)*dcel-tranz+zbcen
         if(zz.gt.zmax.or.zz.lt.zmin) goto 110
         sumioncg=0.0D0
         do ig=1,nclx
            iii=(ig-1)*ncyz+kg
            xx=(ig-1)*dcel-tranx+xbcen
            do 111 jg=1,ncly
               ip0=iii+(jg-1)*nclz
               yy=(jg-1)*dcel-trany+ybcen
               if(xx.le.22.25D0.and.xx.ge.-22.25D0.and.
     $            yy.le.22.25D0.and.yy.ge.-22.25D0) goto 111
               ioncg=0.0D0
               do i=1,ntype
                  ipi=ip0+nc3*(i-1)
                  if(Cion(ipi).ne.0.0d0) then
                     ioncg=ioncg+Zion(i)*Cion(ipi)*dcel3
                  endif
               enddo
               if(ioncg.ne.0.0D0)then
                  anum=anum+1
                  sumioncg=sumioncg+ioncg
                  if(ioncg.ge.0.0D0) then
                     write(iunit,101) anum,anum,'POT ','POT ',
     $                    xx,yy,zz,'POT ',anum,ioncg
                  else
                     write(iunit,101) anum,anum,'CLA ','CLA ',
     $                    xx,yy,zz,'CLA ',anum,ioncg
                  endif
               endif
 111        enddo
         enddo
         write(*,*) zz,sumioncg
 110  enddo
 101  format(i5,1x,i4,1x,a4,1x,a4,3f10.5,1x,a4,1x,i4,f10.5)
c
      return
      end

