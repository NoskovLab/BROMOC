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

      SUBROUTINE PNP1(MAXPNP,MAXPHI,MAXCION,TOLCION,TOLPHI,
     $           Lambda1,Lambda2,
     $           TEMP,Ntype,NCLX,NCLY,NCLZ,DCEL,
     $           PHI,EPSX,EPSY,EPSZ,FCDEN,MCDEN,EFFEPS,WA,
     $           Cion,Zion,Iontype,Diffusion,
     $           QCionXYPBC,QphiXYPBC,QPHI)
c-----------------------------------------------------------------------
c     Nernst-Planck and Poisson solver for discretized lattice.  
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
c NOTE: (1) FCDEN already contained 4*pi/h for the Poisson Eq. solver
c       (2) The Nernst-Planck Eq. solver uses the Slotboom transformation, so
c           the effective dielectric profile and effective concentration.
c       (3) INPUT and OUTPUT of Cion are not effective concentrations
c                                        but just concentrations.
c
      implicit none
      include 'consta.fcm'
      include 'mainio.fcm'

      integer*4 MAXPNP,MAXPHI,MAXCION,Ntype,NCLX,NCLY,NCLZ
      REAL*4  PHI(*),EPSX(*),EPSY(*),EPSZ(*)
      REAL*4  FCDEN(*),MCDEN(*),Cion(*),EFFEPS(*),WA(*)
      REAL*8  TOLCION,TOLPHI,DCEL,TEMP,Zion(*),Diffusion(*)
      REAL*8  Lambda1,Lambda2
      LOGICAL QCionXYPBC,QphiXYPBC,QPHI
      CHARACTER*4 IONTYPE(*)
c Local variables 
      REAL*8  TOP,BOT,RESID,RJAC,Omega1,Omega2,ANORM(10)
      REAL*8  FACTOR1,FACTOR2,FACTOR3,MIONCD,DCEL2
      REAL*4  PHI0
      REAL*4  C0,C1,C2,C3,C4,C5,C6,CSUM
      REAL*4  A0,A1,A2,A3,A4,A5,A6,HA1,HA2,HA3,HA4,HA5,HA6
      REAL*4  EPS1,EPS2,EPS3,EPS4,EPS5,EPS6
      integer*4 I,NCYZ,NC3,Nion,Nphi,NPNP,NC3NP(10),IONCNT(10)
      integer*4 IG,JG,KG,IP0,IPI,IP0X,IP0Y,IP1,IP2,IP3,IP4,IPD
      integer*4 XFIR1,XLAS1,YFIR1,YLAS1,KSW1
      integer*4 XFIR2,XLAS2,YFIR2,YLAS2,KSW2
      integer*4 SUMNion
c
      if(Ntype.GT.9) then
         stop 'Increase dimension of Anorm and Ioncnt in PNP1'
      endif

c some factors
      NCyz=NCLy*NCLz
      NC3=NCLx*NCyz
      DCEL2=DCEL*DCEL
      FACTOR1=celec /(kboltz*Temp/kcalmol)       ! 1/(kcal/(mol*e))->1/(e/A)
      FACTOR2=2.0D0*TWOPI*DCEL2                  ! for MCDEN

c initialize some parameters
      Omega1=1.0D0            ! initial mixing factors for PHI
      Omega2=1.0D0            ! initial mixing factors for Cion
      KSW1=1
      KSW2=1
      RJAC=(COS(PI/(NCLx-1))+COS(PI/(NCLy-1))+COS(PI/(NCLz-1)))/3.0D0
      DO I=1,Ntype+1
         ANORM(I)=0.0d0

      ENDDO

c check boundary conditions for potentials
      IF(QphiXYPBC) THEN
         XFIR1=1
         XLAS1=NCLx
         YFIR1=1
         YLAS1=NCLy
      ELSE
         XFIR1=2
         XLAS1=NCLx-1
         YFIR1=2
         YLAS1=NCLy-1
      ENDIF

c check boundary conditions for concentrations
      IF(QCionXYPBC) THEN
         XFIR2=1
         XLAS2=NCLx
         YFIR2=1
         YLAS2=NCLy
      ELSE
         XFIR2=2
         XLAS2=NCLx-1
         YFIR2=2
         YLAS2=NCLy-1
      ENDIF

c calculate number of ion-accessible grid points and
c initialize EFFEPS (effective dielectric profile) array
      DO I=1,Ntype
         NC3NP(I)=0
      ENDDO
      DO I=1,Ntype
         DO IP0=1,NC3
            IPI=IP0+(I-1)*NC3
            IF(MCDEN(IPI).NE.0.0D0) NC3NP(I)=NC3NP(I)+1
         ENDDO
      ENDDO
      DO I=1,NC3
         EFFEPS(I)=0.0D0
      ENDDO

c check non-zero concentration grid points which are surrounded by zero ones
c if there is, set zero concentration
      DO I=1,Ntype
         DO IG=XFIR2,XLAS2
            IP0X=(IG-1)*NCyz
            IP1=-NCyz
            IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
            IP2=NCyz
            IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
            DO JG=YFIR2,YLAS2
               IP0Y=(JG-1)*NCLz+IP0X
               IP3=-NCLz
               IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
               IP4=NCLz
               IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
               DO KG=2,NCLZ-1
                  IP0=IP0Y+KG
                  IPI=(I-1)*NC3+IP0
                  C0=MCDEN(IPI)
                  IF(C0.EQ.0.0D0) GOTO 100
                  C1=MCDEN(IPI+IP1)
                  C2=MCDEN(IPI+IP2)
                  C3=MCDEN(IPI+IP3)
                  C4=MCDEN(IPI+IP4)
                  C5=MCDEN(IPI-1)
                  C6=MCDEN(IPI+1)
                  CSUM=C1+C2+C3+C4+C5+C6
                  IF(CSUM.EQ.0.0D0) THEN
                     Cion(IPI)=0.0D0
                     WRITE(OUTU,'(6x,a,a4,a)') 
     $                   'ZERO concentration neighbors for ',iontype(I),
     $                   '... set to zero'
                  ENDIF
 100           ENDDO
            ENDDO
         ENDDO
      ENDDO

c concentration -> effective concentration [to modify boundary concentrations]
      IF(QPHI)THEN
         DO I=1,Ntype
            FACTOR3=Zion(I)*FACTOR1
            DO IP0=1,NC3
               IPI=IP0+(I-1)*NC3
               C0=Cion(IPI)
               IF(C0.ne.0.0d0) Cion(IPI)=C0*exp(FACTOR3*PHI(IP0))
            ENDDO
         ENDDO
      ENDIF

c OUTPUT format
      WRITE(OUTU,'(6x,a)')'PNP ITERATIONs'
      WRITE(OUTU,'(6x,a,3x,a,3x,8(a4,4x))') 
     $     'NPNP','Poisson',(Iontype(i),i=1,ntype)

c========================================================================c
c                              PNP MAIN CYCLE                            c
c========================================================================c

      DO NPNP=1,MAXPNP

c use the PBEQ potential for the first guess of the concentrations
         IF(NPNP.EQ.1.AND.QPHI) THEN
            Nphi = 0
            GOTO 101
         ENDIF

c update the charge density of mobile ions (MCDEN array)
c =====================================================

         IF(NPNP.EQ.1) THEN
c concentration
            DO IP0=1,NC3
               MIONCD=0.0D0
               PHI0=PHI(IP0)
               DO I=1,Ntype
                  IPI=(I-1)*NC3+IP0
                  C0=Cion(IPI)
                  IF(C0.ne.0.0d0) THEN
                     MIONCD=MIONCD+Zion(I)*C0
                  ENDIF
               ENDDO
               MCDEN(IP0)=MIONCD*FACTOR2
            ENDDO
         ELSE
c effective concentration
            DO IP0=1,NC3
               MIONCD=0.0D0
               PHI0=PHI(IP0)
               DO I=1,Ntype
                  IPI=(I-1)*NC3+IP0
                  C0=Cion(IPI)
                  IF(C0.ne.0.0d0) THEN
                     MIONCD=MIONCD+
     $                    Zion(I)*C0*exp(-Zion(I)*FACTOR1*PHI0)
                  ENDIF
               ENDDO
               MCDEN(IP0)=MIONCD*FACTOR2
            ENDDO
         ENDIF


c solve Poisson equation and update the electrostatic potentials
c =============================================================
         DO IP0=1,NC3
            WA(IP0)=PHI(IP0)
         ENDDO

         DO Nphi=1,MAXPHI
            ANORM(Ntype+1)=0.0D0
            DO IG=XFIR1,XLAS1
               IP0X=(IG-1)*NCyz
               IP1=-NCyz
               IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
               IP2=NCyz
               IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
               DO JG=YFIR1,YLAS1
                  IP0Y=(JG-1)*NCLz+IP0X
                  IP3=-NCLz
                  IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                  IP4=NCLz
                  IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                  DO KG=KSW1+1,NCLZ-KSW1,2
                     IP0 = IP0Y + KG
                  
                     PHI0=PHI(IP0)
                     EPS1=EPSX(IP0+IP1)
                     EPS2=EPSX(IP0)
                     EPS3=EPSY(IP0+IP3)
                     EPS4=EPSY(IP0)
                     EPS5=EPSZ(IP0-1)
                     EPS6=EPSZ(IP0)
                  
                     BOT=EPS1+EPS2+EPS3+EPS4+EPS5+EPS6

                     TOP=PHI(IP0+IP1)*EPS1+PHI(IP0+IP2)*EPS2+
     $                   PHI(IP0+IP3)*EPS3+PHI(IP0+IP4)*EPS4+
     $                   PHI(IP0-1)*EPS5+PHI(IP0+1)*EPS6 +
     $                   FCDEN(IP0)+MCDEN(IP0)

c phi = phi_old + omega*(phi_new - phi_old)
                     RESID=TOP-PHI0*BOT
                     ANORM(Ntype+1)=ANORM(Ntype+1)+ABS(RESID)
                     PHI(IP0)=PHI0+OMEGA1*RESID/BOT
                  ENDDO
                  KSW1=3-KSW1
               ENDDO
            ENDDO
            IF(NPNP.EQ.1.and.Nphi.EQ.1)THEN
               OMEGA1=1d0/(1d0-0.5d0*RJAC**2)
            ELSE
               OMEGA1=1d0/(1d0-0.25d0*RJAC**2*OMEGA1)
            ENDIF

            ANORM(Ntype+1)=ANORM(Ntype+1)/NC3

            if(nphi.le.2.or.npnp.eq.1) 
     $           write(*,*) 'PHI ', nphi,anorm(ntype+1)

c check the convergence of PHI for the current PNP step
            IF(Nphi.GT.1.and.ANORM(Ntype+1).LT.TOLPHI) GOTO 101 
         ENDDO

         WRITE(OUTU,'(6x,A,2E13.6,A,I3)')
     $        'Calculation not converged for the Poisson EQ.',
     $        ANORM(Ntype+1),TOLPHI,' -- the PNP step :',NPNP

 101     CONTINUE

c mix PHI(new) with PHI(old)
         IF(NPNP.NE.1) THEN
            DO IP0=1,NC3
               PHI(IP0)=Lambda1*PHI(IP0)+(1.D0-Lambda1)*WA(IP0)
            ENDDO
         ENDIF

c concentration -> effective concentration [to modify boundary concentrations]
         IF(NPNP.EQ.1.AND..not.QPHI) THEN
            DO I=1,Ntype
               FACTOR3=Zion(I)*FACTOR1
               DO IP0=1,NC3
                  IPI=IP0+(I-1)*NC3
                  C0=Cion(IPI)
                  IF(C0.ne.0.0D0) Cion(IPI)=C0*exp(FACTOR3*PHI(IP0))
               ENDDO
            ENDDO
         ENDIF

c solve Nernst-Planck equation and update the concentration of each ion
c =====================================================================
         SumNion=0
         DO I=1,Ntype
c setup working array WA and effective dielectric profile array 
            FACTOR3=Zion(I)*FACTOR1
            DO IP0=1,NC3
               IPI=IP0+(I-1)*NC3
               WA(IP0)=Cion(IPI)
               EFFEPS(IP0)=0.0D0
               IF(Cion(IPI).NE.0.0D0) THEN
                  kg=MOD(MOD((IP0-1),NCyz),NCLz)+1
                  IPD=(I-1)*NCLz+kg
                  EFFEPS(IP0)=Diffusion(IPD)*exp(-FACTOR3*PHI(IP0))
               ENDIF
            ENDDO

            Omega2=1.0D0        ! initial mixing factors for Cion
c solve transfromed Nernst-Planck equation
            DO Nion=1,MAXCION
               ANORM(I)=0.0D0
               DO IG=XFIR2,XLAS2
                  IP0X=(IG-1)*NCyz
                  IP1=-NCyz
                  IF(IG.EQ.1) IP1=NC3-NCyz ! (NCLx-1)*NCyz
                  IP2=NCyz
                  IF(IG.EQ.NCLx) IP2=-(NC3-NCyz)
                  DO JG=YFIR2,YLAS2
                     IP0Y=(JG-1)*NCLz+IP0X
                     IP3=-NCLz
                     IF(JG.EQ.1) IP3=NCyz-NCLz ! (NCLy-1)*NCLz
                     IP4=NCLz
                     IF(JG.EQ.NCLy) IP4=-(NCyz-NCLz)
                     DO KG=KSW2+1,NCLZ-KSW2,2
                        IP0=IP0Y+KG
                        IPI=IP0+(I-1)*NC3
                        C0=Cion(IPI)
                        IF(C0.EQ.0.0D0) GOTO 102

                        A0=EFFEPS(IP0)
                        A1=EFFEPS(IP0+IP1)
                        A2=EFFEPS(IP0+IP2)
                        A3=EFFEPS(IP0+IP3)
                        A4=EFFEPS(IP0+IP4)
                        A5=EFFEPS(IP0-1)
                        A6=EFFEPS(IP0+1)
                        HA1=2.0d0*A0*A1/(A0+A1)
                        HA2=2.0d0*A0*A2/(A0+A2)
                        HA3=2.0d0*A0*A3/(A0+A3)
                        HA4=2.0d0*A0*A4/(A0+A4)
                        HA5=2.0d0*A0*A5/(A0+A5)
                        HA6=2.0d0*A0*A6/(A0+A6)
                     
                        BOT=HA1+HA2+HA3+HA4+HA5+HA6

                        TOP=Cion(IPI+IP1)*HA1+Cion(IPI+IP2)*HA2+
     $                      Cion(IPI+IP3)*HA3+Cion(IPI+IP4)*HA4+
     $                      Cion(IPI-1)*HA5+Cion(IPI+1)*HA6

c c = c_old + omega*(c_new - c_old)
                        RESID=TOP-C0*BOT
                        ANORM(I)=ANORM(I)+ABS(RESID)
                        Cion(IPI)=C0+OMEGA2*RESID/BOT

 102                 ENDDO
                     KSW2=3-KSW2
                  ENDDO
               ENDDO

               IF(NPNP.EQ.1.and.Nion.EQ.1)THEN
                  OMEGA2=1d0/(1d0-0.5d0*RJAC**2)
               ELSE
                  OMEGA2=1d0/(1d0-0.25d0*RJAC**2*OMEGA2)
               ENDIF

               ANORM(I)=ANORM(I)/NC3NP(I)

               if(nion.le.2.or.npnp.eq.1) 
     $              write(*,*) iontype(i),nion,anorm(i)

               IF(Nion.GT.1.and.ANORM(I).LT.TOLCION) GOTO 103
            ENDDO
            WRITE(OUTU,'(6x,A,A4,2E13.6,A,I3)')
     $           'Calculation not converged for ',Iontype(I),
     $           ANORM(I),TOLCION,' -- the PNP step :',NPNP

 103        CONTINUE

c mix Cion(new) with Cion(old)
            IF(NPNP.NE.1) THEN
               DO IP0=1,NC3
                  IPI=IP0+(I-1)*NC3
                  Cion(IPI)=Lambda2*Cion(IPI)+(1.D0-Lambda2)*WA(IP0)
               ENDDO
            ENDIF

            SumNion=SumNion+Nion
            IONCNT(I)=Nion
         ENDDO

c check the divergence or convergence of Cion(i) and PHI 
c ======================================================
         DO I=1,Ntype+1
            IF(ANORM(I).GT.1.0d0) THEN
               WRITE(OUTU,'(6X,A)') 'PNP diverged'
               GOTO 104
            ENDIF
         ENDDO

         IF(SumNion.EQ.Ntype*2.and.Nphi.EQ.2) GOTO 104

         WRITE(OUTU,'(6X,i4,6x,i4,2x,8(i4,4x))') 
     $        NPNP,NPHI,(IONCNT(i),i=1,ntype)

      ENDDO     ! end of PNP cycle

      WRITE(OUTU,'(6X,A)') 'PNP does not converged'
      WRITE(OUTU,'(6X,i5,3x,i5,4x,8(i5,3x))') 
     $     NPNP,NPHI,(IONCNT(i),i=1,ntype)

 104  CONTINUE
      WRITE(OUTU,'(6X,i4,3x,i5,4x,8(i5,3x))') 
     $     NPNP,NPHI,(IONCNT(i),i=1,ntype)

      WRITE(OUTU,'(6X,A,I5)') 'Number of iterations: ',NPNP

c effective concentration -> concentration
      DO I=1,Ntype
         FACTOR3=Zion(I)*FACTOR1
         DO IP0=1,NC3
            IPI=IP0+(I-1)*NC3
            C0=Cion(IPI)
            IF(C0.ne.0.0d0) Cion(IPI)=C0*exp(-FACTOR3*PHI(IP0))
         ENDDO
      ENDDO
c
      RETURN
      END

c
c      IG=int(Nclx/2+rsmall)+1
c      JG=int(Ncly/2+rsmall)+1
c      DO KG=1,NCLz
c         IP0=(IG-1)*NCyz+(JG-1)*NCLz+KG
c         write(6,'(f10.5,4e20.10)') (KG-1)*Dcel-(NCLz-1)*DCEL*0.5,
c     $        PHI(IP0),
c     $        Cion(IP0)*exp(-Zion(1)*FACTOR1*PHI(IP0)),
c     $        Cion(IP0+NC3)*exp(-Zion(2)*FACTOR1*PHI(IP0)),
c     $        Cion(IP0)*exp(-Zion(1)*FACTOR1*PHI(IP0)) -
c     $        Cion(IP0+NC3)*exp(-Zion(2)*FACTOR1*PHI(IP0))
c      ENDDO
c

