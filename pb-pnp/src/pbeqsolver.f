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

      SUBROUTINE PBEQ1(MAXPHI,TOLPHI,KAPPA2,KAP2TOP,KAP2BOT,
     $           NCLX,NCLY,NCLZ,DCEL,TMEMB,ZMEMB,TRANZ,ZBCEN,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,QphiXYPBC,QphiXYZPBC)
c-----------------------------------------------------------------------
c     (Linearized) Poisson-Boltzmann solver for discretized lattice.  
c     It is assumed that all arrays have been prepared elsewhere.
c
c                        (in Z)
c                        ip0+1=ip6
c                               | ip3
c                               | /
c                               |/
c                    ip1 ----- ip0 ----- ip2   (in X)
c                              /| 
c                             / |
c                           ip4 |
c                    (in Y)    ip5=ip0-1 
c
c NOTE: FCDEN already contained 4*pi/h
c
      implicit none
      include 'consta.fcm'
      include 'mainio.fcm'

      integer*4 MAXPHI,NCLX,NCLY,NCLZ
      REAL*4  PHI(*),EPSX(*),EPSY(*),EPSZ(*),FCDEN(*),MCDEN(*)
      REAL*8  TOLPHI,DCEL,KAPPA2,KAP2TOP,KAP2BOT
      REAL*8  TMEMB,ZMEMB,TRANZ,ZBCEN
      LOGICAL QphiXYPBC,QphiXYZPBC


c Local variables 
      REAL*8  DCEL2,TOP,BOT,RESID,RJAC,Omega,ANORM
      REAL*8  ZC,ZMEMBCENTER
      REAL*4  PHI0,EPS1,EPS2,EPS3,EPS4,EPS5,EPS6
      integer*4 I,NCYZ,NC3,Nphi
      integer*4 IG,JG,KG,IP0,IP0X,IP0Y,IP1,IP2,IP3,IP4,IP5,IP6
      integer*4 KSW
c
      NCyz=NCLy*NCLz
      NC3=NCLx*NCyz
      DCEL2=DCEL*DCEL

c initialize some parameters
      RJAC=(COS(PI/(NCLx-1))+COS(PI/(NCLy-1))+COS(PI/(NCLz-1)))/3.0D0
      Omega=1.0D0  ! initial mixing factors
      KSW=1

c OUTPUT format
      WRITE(OUTU,'(6x,a)')'Linearized PBEQ ITERATIONs'

c convert MCDEN array to kappa2
      zmembcenter=tranz-zbcen+zmemb-rsmall
      do kg=1,nclz
         zc=(kg-1)*dcel
         IF(zc.lt.zmembcenter) THEN
c            write(*,*) 'bot',zc,kap2bot
            do ig=1,nclx
               ip0x=(ig-1)*ncyz
               do jg=1,ncly
                  ip0=ip0x+(jg-1)*nclz+kg
                  mcden(ip0)=mcden(ip0)*kap2bot*dcel2
               enddo
            enddo
         ELSEIF(zc.ge.zmembcenter) THEN
c            write(*,*) 'top',zc,kap2top
            do ig=1,nclx
               ip0x=(ig-1)*ncyz
               do jg=1,ncly
                  ip0=ip0x+(jg-1)*nclz+kg
                  mcden(ip0)=mcden(ip0)*kap2top*dcel2
               enddo
            enddo
c         ELSE
c            write(*,*) 'middle',zc,kappa2
c            do ig=1,nclx
c               ip0x=(ig-1)*ncyz
c               do jg=1,ncly
c                  ip0=ip0x+(jg-1)*nclz+kg
c                  mcden(ip0)=mcden(ip0)*kappa2*dcel2
c               enddo
c            enddo
         ENDIF
      enddo

c========================================================================c
c                    Linearized PBEQ MAIN CYCLE                          c
c========================================================================c

      IF(QphiXYZPBC) THEN
         DO Nphi=1,MAXPHI
            ANORM=0.0D0
            DO IG=1,NCLx
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=1,NCLy
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW,NCLZ-KSW+1,2
                     IP0= IP0Y + KG
                     IP5=-1
                     IF(KG.EQ.1) IP5=NCLz-1
                     IP6=1
                     IF(KG.EQ.NCLz) IP6=-(NCLz-1)
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0+IP5)
                     EPS6=EPSZ(IP0)
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MCDEN(IP0)
                     
                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0+IP5)*EPS5+PHI(IP0+IP6)*EPS6 +
     $                   FCDEN(IP0)

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO
            IF(Nphi.EQ.1)THEN
               OMEGA=1d0/(1d0-0.5d0*RJAC**2)
            ELSE
               OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
            ENDIF

            ANORM=ANORM/NC3

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ELSEIF(QphiXYPBC) THEN
         DO Nphi=1,MAXPHI
            ANORM=0.0D0
            DO IG=1,NCLx
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=1,NCLy
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW+1,NCLZ-KSW,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MCDEN(IP0)
                     
                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FCDEN(IP0)

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO
            IF(Nphi.EQ.1)THEN
               OMEGA=1d0/(1d0-0.5d0*RJAC**2)
            ELSE
               OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
            ENDIF

            ANORM=ANORM/NC3

c            write(*,*) nphi,omega,anorm
            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ELSE

         DO Nphi=1,MAXPHI
            ANORM=0.0D0
            DO IG=2,NCLx-1
               IP0X=(IG-1)*NCyz
               DO JG=2,NCLy-1
                  IP0Y=(JG-1)*NCLz+IP0X
                  DO KG=KSW+1,NCLZ-KSW,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0-NCyz)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0-NCLz)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MCDEN(IP0)
                     
                     TOP=PHI(IP0-NCyz)*EPS1+PHI(IP0+NCyz)*EPS2+
     $                   PHI(IP0-NCLz)*EPS3+PHI(IP0+NCLz)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FCDEN(IP0)

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO
            IF(Nphi.EQ.1)THEN
               OMEGA=1d0/(1d0-0.5d0*RJAC**2)
            ELSE
               OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
            ENDIF

            ANORM=ANORM/NC3

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ENDIF

      WRITE(OUTU,'(6x,A,2E13.6)')
     $     'Calculation not converged for the Poisson-Boltzmann EQ.',
     $     ANORM,TOLPHI

 101  CONTINUE
      WRITE(OUTU,'(6X,A,I5)') 'Number of iterations: ',Nphi
c
      RETURN
      END


      SUBROUTINE PBEQ2(MAXPHI,TOLPHI,KAPPA2,Temp,Omega,
     $           NCLX,NCLY,NCLZ,DCEL,VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,
     $           QphiXYPBC,QphiXYZPBC,Qunder)
c-----------------------------------------------------------------------
c     (Non-linear) Poisson-Boltzmann solver for discretized lattice.  
c     It is assumed that all arrays have been prepared elsewhere.
c
c                        (in Z)
c                        ip0+1=ip6
c                               | ip3
c                               | /
c                               |/
c                    ip1 ----- ip0 ----- ip2   (in X)
c                              /| 
c                             / |
c                           ip4 |
c                    (in Y)    ip5=ip0-1 
c
c NOTE: FCDEN already contained 4*pi/h
c
      implicit none
      include 'consta.fcm'
      include 'mainio.fcm'

      integer*4 MAXPHI,NCLX,NCLY,NCLZ
      REAL*4  PHI(*),EPSX(*),EPSY(*),EPSZ(*),FCDEN(*),MCDEN(*)
      REAL*8  TOLPHI,DCEL,KAPPA2,Temp,Omega
      REAL*8  VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN
      LOGICAL QphiXYPBC,QphiXYZPBC,Qunder

c Local variables 
      REAL*8  DCEL2,TOP,BOT,RESID,RJAC,ANORM
      REAL*8  FACTOR1,PHIF,RATIO,MIONCD
      REAL*8  ZC,ZMEMB2,MKAPPA2,FPROTCD
      REAL*4  PHI0,EPS1,EPS2,EPS3,EPS4,EPS5,EPS6
      integer*4 I,NCYZ,NC3,Nphi
      integer*4 IG,JG,KG,IP0,IP0X,IP0Y,IP1,IP2,IP3,IP4,IP5,IP6
      integer*4 KSW

c some factors
      NCyz=NCLy*NCLz
      NC3=NCLx*NCyz
      DCEL2=DCEL*DCEL
      FACTOR1=celec /(kboltz*Temp/kcalmol)       ! 1/(kcal/(mol*e))->1/(e/A)
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall
      MKAPPA2=KAPPA2*DCEL2

c initialize some parameters
      RJAC=(COS(PI/(NCLx-1))+COS(PI/(NCLy-1))+COS(PI/(NCLz-1)))/3.0D0
      KSW=1

c OUTPUT format
      WRITE(OUTU,'(6x,a)')
     $     'Non-linear PBEQ ITERATIONs for 1:1 charge-paired salt '

c========================================================================c
c                    Non-linear PBEQ MAIN CYCLE                          c
c========================================================================c

      IF(QphiXYZPBC) THEN
         DO Nphi=1,MAXPHI
            ANORM=0.0D0

            DO IG=1,NCLx
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=1,NCLy
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW,NCLZ-KSW+1,2
                     IP0= IP0Y + KG
                     IP5=-1
                     IF(KG.EQ.1) IP5=NCLz-1
                     IP6=1
                     IF(KG.EQ.NCLz) IP6=-(NCLz-1)
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0+IP5)
                     EPS6=EPSZ(IP0)

c convert MCDEN array to (kappa2 * nonlinearity ratio)
c add the background charge (kappa2 * nonlinearity ratio * Vmemb) to FCDEN
                     MIONCD=0.0d0
                     FPROTCD=FCDEN(IP0)
                     IF(MCDEN(IP0).ne.0.0d0) THEN
                        PHIF=PHI0*FACTOR1
                        IF(VMEMB.NE.0.0d0)THEN
                           ZC=(KG-1)*DCEL
                           IF(ZC.GT.ZMEMB2) PHIF=(PHI0-VMEMB)*FACTOR1
                        ENDIF
                        RATIO=EXP(PHIF)
                        RATIO=(RATIO-1.0d0/RATIO)*0.5d0/PHIF
                        MIONCD=MKAPPA2*RATIO
                        IF(ZC.GT.ZMEMB2) FPROTCD=FPROTCD+MIONCD*VMEMB
                     ENDIF
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MIONCD
                     
                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0+IP5)*EPS5+PHI(IP0+IP6)*EPS6 +
     $                   FPROTCD

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO

            IF(.not.Qunder) THEN
               IF(Nphi.EQ.1)THEN
                  OMEGA=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
               ENDIF
            ENDIF

            ANORM=ANORM/NC3

            write(*,*) nphi,omega,anorm

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ELSEIF(QphiXYPBC) THEN
         DO Nphi=1,MAXPHI
            ANORM=0.0D0

            DO IG=1,NCLx
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=1,NCLy
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW+1,NCLZ-KSW,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)

c convert MCDEN array to ( kappa2 * nonlinearity ratio)
c add the background charge (kappa2 * nonlinearity ratio * Vmemb) to FCDEN
                     MIONCD=0.0d0
                     FPROTCD=FCDEN(IP0)
                     IF(MCDEN(IP0).ne.0.0d0) THEN
                        PHIF=PHI0*FACTOR1
                        IF(VMEMB.NE.0.0d0)THEN
                           ZC=(KG-1)*DCEL
                           IF(ZC.GT.ZMEMB2) PHIF=(PHI0-VMEMB)*FACTOR1
                        ENDIF
                        RATIO=EXP(PHIF)
                        RATIO=(RATIO-1.0d0/RATIO)*0.5d0/PHIF
                        MIONCD=MKAPPA2*RATIO
                        IF(ZC.GT.ZMEMB2) FPROTCD=FPROTCD+MIONCD*VMEMB
                     ENDIF
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MIONCD
                     
                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FPROTCD

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO

            IF(.not.Qunder) THEN
               IF(Nphi.EQ.1)THEN
                  OMEGA=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
               ENDIF
            ENDIF

            ANORM=ANORM/NC3

            write(*,*) nphi,omega,anorm

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ELSE

         DO Nphi=1,MAXPHI
            ANORM=0.0D0

            DO IG=2,NCLx-1
               IP0X=(IG-1)*NCyz
               DO JG=2,NCLy-1
                  IP0Y=(JG-1)*NCLz+IP0X
                  DO KG=KSW+1,NCLZ-KSW,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0-NCyz)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0-NCLz)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)

c convert MCDEN array to ( kappa2 * nonlinearity ratio)
c add the background charge (kappa2 * nonlinearity ratio * Vmemb) to FCDEN
                     MIONCD=0.0d0
                     FPROTCD=FCDEN(IP0)
                     IF(MCDEN(IP0).ne.0.0d0) THEN
                        PHIF=PHI0*FACTOR1
                        IF(VMEMB.NE.0.0d0)THEN ! for focussing with vmemb
                           ZC=(KG-1)*DCEL
                           IF(ZC.GT.ZMEMB2) PHIF=(PHI0-VMEMB)*FACTOR1
                        ENDIF
                        RATIO=EXP(PHIF)
                        RATIO=(RATIO-1.0d0/RATIO)*0.5d0/PHIF
                        MIONCD=MKAPPA2*RATIO
                        IF(ZC.GT.ZMEMB2) FPROTCD=FPROTCD+MIONCD*VMEMB
                     ENDIF
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MIONCD
                     
                     TOP=PHI(IP0-NCyz)*EPS1+PHI(IP0+NCyz)*EPS2+
     $                   PHI(IP0-NCLz)*EPS3+PHI(IP0+NCLz)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FPROTCD

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO

            IF(.not.Qunder) THEN
               IF(Nphi.EQ.1)THEN
                  OMEGA=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
               ENDIF
            ENDIF

            ANORM=ANORM/NC3

            write(*,*) nphi,omega,anorm

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ENDIF

      WRITE(OUTU,'(6x,A,2E13.6)')
     $     'Calculation not converged for the Poisson-Boltzmann EQ.',
     $     ANORM,TOLPHI

 101  CONTINUE
      WRITE(OUTU,'(6X,A,I5)') 'Number of iterations: ',Nphi
c
      RETURN
      END


      SUBROUTINE PBEQ3(MAXPHI,TOLPHI,KAPPA2,Temp,Omega,
     $           NCLX,NCLY,NCLZ,DCEL,VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,
     $           QphiXYPBC,QphiXYZPBC,Qunder)
c-----------------------------------------------------------------------
c     (Partially linearized) Poisson-Boltzmann solver for discretized lattice.  
c     It is assumed that all arrays have been prepared elsewhere.
c
c                        (in Z)
c                        ip0+1=ip6
c                               | ip3
c                               | /
c                               |/
c                    ip1 ----- ip0 ----- ip2   (in X)
c                              /| 
c                             / |
c                           ip4 |
c                    (in Y)    ip5=ip0-1 
c
c NOTE: FCDEN already contained 4*pi/h
c
      implicit none
      include 'consta.fcm'
      include 'mainio.fcm'

      integer*4 MAXPHI,NCLX,NCLY,NCLZ
      REAL*4  PHI(*),EPSX(*),EPSY(*),EPSZ(*),FCDEN(*),MCDEN(*)
      REAL*8  TOLPHI,DCEL,KAPPA2,Temp,Omega
      REAL*8  VMEMB,TMEMB,ZMEMB,TRANZ,ZBCEN
      LOGICAL QphiXYPBC,QphiXYZPBC,Qunder

c Local variables 
      REAL*8  DCEL2,TOP,BOT,RESID,RJAC,ANORM
      REAL*8  FACTOR1,PHIF,RATIO,MIONCD
      REAL*8  ZC,ZMEMB2,MKAPPA2,FPROTCD
      REAL*4  PHI0,EPS1,EPS2,EPS3,EPS4,EPS5,EPS6
      integer*4 I,NCYZ,NC3,Nphi
      integer*4 IG,JG,KG,IP0,IP0X,IP0Y,IP1,IP2,IP3,IP4,IP5,IP6
      integer*4 KSW

c some factors
      NCyz=NCLy*NCLz
      NC3=NCLx*NCyz
      DCEL2=DCEL*DCEL
      FACTOR1=celec /(kboltz*Temp/kcalmol)       ! 1/(kcal/(mol*e))->1/(e/A)
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall
      MKAPPA2=KAPPA2*DCEL2

c initialize some parameters
      RJAC=(COS(PI/(NCLx-1))+COS(PI/(NCLy-1))+COS(PI/(NCLz-1)))/3.0D0
      KSW=1

c OUTPUT format
      WRITE(OUTU,'(6x,a)')
     $'Partially linearized PBEQ ITERATIONs for 1:1 charge-paired salt'

c========================================================================c
c                    Non-linear PBEQ MAIN CYCLE                          c
c========================================================================c

      IF(QphiXYZPBC) THEN
         DO Nphi=1,MAXPHI
            ANORM=0.0D0

            DO IG=1,NCLx
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=1,NCLy
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW,NCLZ-KSW+1,2
                     IP0= IP0Y + KG
                     IP5=-1
                     IF(KG.EQ.1) IP5=NCLz-1
                     IP6=1
                     IF(KG.EQ.NCLz) IP6=-(NCLz-1)
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0+IP5)
                     EPS6=EPSZ(IP0)

c convert MCDEN array to ( kappa2 * nonlinearity ratio)
c add the background charge (kappa2 * nonlinearity ratio * Vmemb) to FCDEN
                     MIONCD=0.0d0
                     FPROTCD=FCDEN(IP0)
                     IF(MCDEN(IP0).ne.0.0d0) THEN
                        PHIF=PHI0*FACTOR1
                        IF(VMEMB.NE.0.0d0)THEN
                           ZC=(KG-1)*DCEL
                           IF(ZC.GT.ZMEMB2) PHIF=(PHI0-VMEMB)*FACTOR1
                        ENDIF
                        IF(PHIF.GT.0.0d0) THEN
                           RATIO=(1.0d0+PHIF-EXP(-PHIF))*0.5d0/PHIF
                           MIONCD=MKAPPA2*RATIO
                        ELSE
                           RATIO=(EXP(PHIF)-1.0d0+PHIF)*0.5d0/PHIF
                           MIONCD=MKAPPA2*RATIO
                        ENDIF
                        IF(ZC.GT.ZMEMB2) FPROTCD=FPROTCD+MIONCD*VMEMB
                     ENDIF
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MIONCD
                     
                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0+IP5)*EPS5+PHI(IP0+IP6)*EPS6 +
     $                   FPROTCD

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO

            IF(.not.Qunder) THEN
               IF(Nphi.EQ.1)THEN
                  OMEGA=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
               ENDIF
            ENDIF

            ANORM=ANORM/NC3

            write(*,*) nphi,omega,anorm

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ELSEIF(QphiXYPBC) THEN
         DO Nphi=1,MAXPHI
            ANORM=0.0D0

            DO IG=1,NCLx
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=1,NCLy
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW+1,NCLZ-KSW,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)

c convert MCDEN array to ( kappa2 * nonlinearity ratio)
c add the background charge (kappa2 * nonlinearity ratio * Vmemb) to FCDEN
                     MIONCD=0.0d0
                     FPROTCD=FCDEN(IP0)
                     IF(MCDEN(IP0).ne.0.0d0) THEN
                        PHIF=PHI0*FACTOR1
                        IF(VMEMB.NE.0.0d0)THEN
                           ZC=(KG-1)*DCEL
                           IF(ZC.GT.ZMEMB2) PHIF=(PHI0-VMEMB)*FACTOR1
                        ENDIF
                        IF(PHIF.GT.0.0d0) THEN
                           RATIO=(1.0d0+PHIF-EXP(-PHIF))*0.5d0/PHIF
                           MIONCD=MKAPPA2*RATIO
                        ELSE
                           RATIO=(EXP(PHIF)-1.0d0+PHIF)*0.5d0/PHIF
                           MIONCD=MKAPPA2*RATIO
                        ENDIF
                        IF(ZC.GT.ZMEMB2) FPROTCD=FPROTCD+MIONCD*VMEMB
                     ENDIF
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MIONCD
                     
                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FPROTCD

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO

            IF(.not.Qunder) THEN
               IF(Nphi.EQ.1)THEN
                  OMEGA=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
               ENDIF
            ENDIF

            ANORM=ANORM/NC3

            write(*,*) nphi,omega,anorm

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ELSE

         DO Nphi=1,MAXPHI
            ANORM=0.0D0

            DO IG=2,NCLx-1
               IP0X=(IG-1)*NCyz
               DO JG=2,NCLy-1
                  IP0Y=(JG-1)*NCLz+IP0X
                  DO KG=KSW+1,NCLZ-KSW,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0-NCyz)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0-NCLz)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)

c convert MCDEN array to ( kappa2 * nonlinearity ratio)
c add the background charge (kappa2 * nonlinearity ratio * Vmemb) to FCDEN
                     MIONCD=0.0d0
                     FPROTCD=FCDEN(IP0)
                     IF(MCDEN(IP0).ne.0.0d0) THEN
                        PHIF=PHI0*FACTOR1
                        IF(VMEMB.NE.0.0d0)THEN ! for focussing with vmemb
                           ZC=(KG-1)*DCEL
                           IF(ZC.GT.ZMEMB2) PHIF=(PHI0-VMEMB)*FACTOR1
                        ENDIF
                        IF(PHIF.GT.0.0d0) THEN
                           RATIO=(1.0d0+PHIF-EXP(-PHIF))*0.5d0/PHIF
                           MIONCD=MKAPPA2*RATIO
                        ELSE
                           RATIO=(EXP(PHIF)-1.0d0+PHIF)*0.5d0/PHIF
                           MIONCD=MKAPPA2*RATIO
                        ENDIF
                        IF(ZC.GT.ZMEMB2) FPROTCD=FPROTCD+MIONCD*VMEMB
                     ENDIF
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6+MIONCD
                     
                     TOP=PHI(IP0-NCyz)*EPS1+PHI(IP0+NCyz)*EPS2+
     $                   PHI(IP0-NCLz)*EPS3+PHI(IP0+NCLz)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FPROTCD

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM=ANORM+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA*RESID/BOT
                  ENDDO
                  KSW=3-KSW
               ENDDO
            ENDDO

            IF(.not.Qunder) THEN
               IF(Nphi.EQ.1)THEN
                  OMEGA=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
               ENDIF
            ENDIF

            ANORM=ANORM/NC3

            write(*,*) nphi,omega,anorm

            IF((Nphi.GT.1).AND.(ANORM.LT.TOLPHI)) GOTO 101 
         ENDDO

      ENDIF

      WRITE(OUTU,'(6x,A,2E13.6)')
     $     'Calculation not converged for the Poisson-Boltzmann EQ.',
     $     ANORM,TOLPHI

 101  CONTINUE
      WRITE(OUTU,'(6X,A,I5)') 'Number of iterations: ',Nphi
c
      RETURN
      END


      SUBROUTINE PBEQ4(MAXITS,TOL,EPSW,NCLX,NCLY,NCLZ,
     $           PHI,CDEN,QFOCUS,QphiXYPBC,QphiXYZPBC)
c-----------------------------------------------------------------------
c     calculation in bulk solution for non-orthogonal basis set in GSBP.
c
      implicit none
      include 'consta.fcm'
      include 'mainio.fcm'

      REAL*4  PHI(*),CDEN(*)
      REAL*8  OMEGA,TOL,EPSW
      integer*4 MAXITS,NCLX,NCLY,NCLZ
      LOGICAL QFOCUS,QphiXYPBC,QphiXYZPBC
c Local variables 
      REAL*8  RESID,ANORM,E123,RJAC
      integer*4 NCYZ,NC3,N
      integer*4 IG,JG,KG,IP0,IP0X,IP1X,IP2X,IP0Y,IP3Y,IP4Y,IP5Z,IP6Z
      integer*4 KSW
c
      NCyz=NCLy*NCLz
      NC3=NCLx*NCyz

c initialize some parameters
      RJAC=(COS(PI/(NCLx-1))+COS(PI/(NCLy-1))+COS(PI/(NCLz-1)))/3.0D0
      Omega=1.0D0  ! initial mixing factors
      KSW=1

      IF(QphiXYZPBC) THEN
c periodic in XYZ

      DO N=1,MAXITS
         ANORM=0.0D0

         E123=6.0D0*EPSW
         DO IG=1,NCLx
            IP0X=(IG-1)*NCyz
            IP1X=-NCyz
            IF(IG.EQ.1) IP1X=NC3-NCyz     ! (NCLx-1)*NCyz
            IP2X=NCyz
            IF(IG.EQ.NCLx) IP2X=-(NC3-NCyz)
            DO JG=1,NCLy
               IP0Y=(JG-1)*NCLz+IP0X
               IP3Y=-NCLz
               IF(JG.EQ.1) IP3Y=NCyz-NCLz ! (NCLy-1)*NCLz
               IP4Y=NCLz
               IF(JG.EQ.NCLy) IP4Y=-(NCyz-NCLz)
cwi               DO KG=1,NCLz
               DO KG=KSW,NCLZ-KSW+1,2
                  IP0 = IP0Y + KG
                  IP5Z=-1
                  IF(KG.EQ.1) IP5Z=NCLz-1
                  IP6Z=1
                  IF(KG.EQ.NCLz) IP6Z=-(NCLz-1)
                  RESID=(PHI(IP0+IP1X)+PHI(IP0+IP2X)+
     $                   PHI(IP0+IP3Y)+PHI(IP0+IP4Y)+
     $                   PHI(IP0+IP5Z)+PHI(IP0+IP6Z))*EPSW-
     $                   PHI(IP0)*E123+CDEN(IP0)
                  ANORM=ANORM+ABS(RESID)
                  PHI(IP0)=PHI(IP0)+OMEGA*RESID/E123
               ENDDO
               KSW=3-KSW
            ENDDO
         ENDDO

         ANORM=ANORM/NC3

         IF(N.EQ.1)THEN
            OMEGA=1d0/(1d0-0.5d0*RJAC**2)
         ELSE
            OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
         ENDIF

         IF((N.GT.1).AND.(ANORM.LT.TOL)) GOTO 19

      ENDDO

c     ---------------------------------------
      ELSEIF(QphiXYPBC .AND. .NOT.QFOCUS)THEN
c     ---------------------------------------
c periodic in XY, Z fixed edge-value (in membrane) 

      DO N=1,MAXITS
         ANORM=0.0D0

         E123=6.0D0*EPSW
         DO IG=1,NCLx
            IP0X=(IG-1)*NCyz
            IP1X=-NCyz
            IF(IG.EQ.1) IP1X=NC3-NCyz     ! (NCLx-1)*NCyz
            IP2X=NCyz
            IF(IG.EQ.NCLx) IP2X=-(NC3-NCyz)
            DO JG=1,NCLy
               IP0Y=(JG-1)*NCLz+IP0X
               IP3Y=-NCLz
               IF(JG.EQ.1) IP3Y=NCyz-NCLz ! (NCLy-1)*NCLz
               IP4Y=NCLz
               IF(JG.EQ.NCLy) IP4Y=-(NCyz-NCLz)
               DO KG=KSW+1,NCLZ-KSW,2
                  IP0 = IP0Y + KG
                  RESID=(PHI(IP0+IP1X)+PHI(IP0+IP2X)+
     $                   PHI(IP0+IP3Y)+PHI(IP0+IP4Y)+
     $                   PHI(IP0-1)+PHI(IP0+1))*EPSW-
     $                   PHI(IP0)*E123+CDEN(IP0)
                  ANORM=ANORM+ABS(RESID)
                  PHI(IP0)=PHI(IP0)+OMEGA*RESID/E123
               ENDDO
               KSW=3-KSW
            ENDDO
         ENDDO

         ANORM=ANORM/NC3

         IF(N.EQ.1)THEN
            OMEGA=1d0/(1d0-0.5d0*RJAC**2)
         ELSE
            OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
         ENDIF

         IF((N.GT.1).AND.(ANORM.LT.TOL)) GOTO 19

      ENDDO

c     ----
      ELSE
c     ----

      DO N=1,MAXITS
         ANORM=0.0D0

         DO IG=2,NCLx-1
            IP0X=(IG-1)*NCyz
            DO JG=2,NCLy-1
               IP0Y=(JG-1)*NCLz+IP0X
               DO KG=KSW+1,NCLZ-KSW,2
                  IP0 = IP0Y + KG
                  E123=6.0D0*EPSW
                  RESID=(PHI(IP0-NCyz)+PHI(IP0+NCyz)+
     $                   PHI(IP0-NCLz)+PHI(IP0+NCLz)+
     $                   PHI(IP0-1)+PHI(IP0+1))*EPSW -
     $                   PHI(IP0)*E123+CDEN(IP0)
                  ANORM=ANORM+ABS(RESID)
                  PHI(IP0)=PHI(IP0)+OMEGA*RESID/E123
               ENDDO
               KSW=3-KSW
            ENDDO
         ENDDO

         ANORM=ANORM/NC3

         IF(N.EQ.1)THEN
            OMEGA=1d0/(1d0-0.5d0*RJAC**2)
         ELSE
            OMEGA=1d0/(1d0-0.25d0*RJAC**2*OMEGA)
         ENDIF

         IF((N.GT.1).AND.(ANORM.LT.TOL)) GOTO 19

      ENDDO

      ENDIF
c     -----

      WRITE(OUTU,'(6x,A,2X,2E13.6)')
     $     'Calculation not converged, MAXITS is too small ',ANORM,TOL

 19   CONTINUE
c
      RETURN
      END


      SUBROUTINE ENPB1(NATOM,X,Y,Z,CG,NCLX,NCLY,NCLZ,DCEL,PHI,
     $           TRANX,TRANY,TRANZ,XBCEN,YBCEN,ZBCEN,FACTOR)
c-----------------------------------------------------------------------
c This subroutine computes the total electrostatic free energy of solvation.
c
      implicit none
      include 'consta.fcm'
      include 'mainio.fcm'
      integer*4 NATOM,NCLX,NCLY,NCLZ
      real*8  x(*),y(*),z(*),cg(*)
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen,FACTOR
      real*4  phi(*)
c local
      integer*4 ncyz,i,ix,iy,iz,n1,in1,n2,in2,n3,in3
      real*8  chi,xi,yi,zi,ai,bi,ci,fi,enpb
c
      ncyz=ncly*nclz
      enpb=0.0D0

c Main loop by atoms
      do 101 i=1,natom
         chi=cg(i)
         IF(CHI.EQ.0.0d0) GOTO 101
         xi=x(i)+tranx-xbcen
         yi=y(i)+trany-ybcen
         zi=z(i)+tranz-zbcen
         if(xi.lt.0.0D0 .or. xi.gt.2*tranx .or.
     $      yi.lt.0.0D0 .or. yi.gt.2*trany .or.
     $      zi.lt.0.0D0 .or. zi.gt.2*tranz) GOTO 101
         
         ix=int(xi/dcel+rsmall)+1
         iy=int(yi/dcel+rsmall)+1
         iz=int(zi/dcel+rsmall)+1

c Atom charge distribution by 8 adjacent grid points
         do n1=1,2
            in1=ix+n1-1
            ai=abs(xi-(in1-1)*dcel)/dcel
            in1=(in1-1)*ncyz
            ai=1d0-ai
               
            do n2=1,2
               in2=iy+n2-1
               bi=abs(yi-(in2-1)*dcel)/dcel
               in2=(in2-1)*nclz
               bi=ai*(1d0-bi)
                  
               do n3=1,2
                  in3=iz+n3-1
                  ci=abs(zi-(in3-1)*dcel)/dcel
                  fi=bi*(1d0-ci)
                  in3=in1+in2+in3

                  enpb=enpb+FACTOR*fi*chi*phi(in3)*CELEC

               enddo
            enddo
         enddo
 101  enddo
c
      WRITE(OUTU,'(6X,A,F13.5)')
     $     'Electrostatic energy [KCAL/MOL] = ',ENPB

      RETURN
      END
