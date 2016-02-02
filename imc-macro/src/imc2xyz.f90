!    IMC-MACRO - Inverse Monte Carlo for Macromolecules
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

module invmc
!  Declarations
!   NTOT - max number of particles
!   NTYPES - max number of particles types
!   NAMAX - max number of grid points
!   NKVM - max number of k-vectors in reciprocal Ewald
!
implicit none
real*8,parameter :: PI=3.1415926535897932d0 
logical LMU,LXMOL
integer,parameter :: NTYPES=1000
integer NPOTM,NKVM,NSPEC(NTYPES),NTOT
real*8 CH(NTYPES),rnd,inint
real*8,allocatable :: X(:),Y(:),Z(:),Q(:)
real*8,allocatable :: RDF(:),POT(:,:,:),RAS(:) 
integer*4,allocatable :: IADR(:),ITYPE(:),IPOT(:,:,:),IND(:,:,:),INA(:),ITYP1(:),ITYP2(:)
real*8 BOXL,HBOXL,VOL,ENER,DR,DE,COULF,AF,FQ,ALPHA,FK,PIAL,FK2,ENERS,EFUR,RCUT,RMIN,RMAX,AVEXP,AVE2,ELIM,VIR,VIRS,VIRE
integer NOP,NPOT,NA,ISTEP,NMKS0,NKV,NAV,IMU,NPAIR,IPRINT
real*8,allocatable :: CORS(:),CORP(:,:)
integer*4,allocatable :: KX(:),KY(:),KZ(:)
real*8,allocatable :: RKV(:),SSIN(:),SCOS(:),DDSIN(:),DDCOS(:)
real*8,allocatable :: STAB(:,:),CTAB(:,:)
external rnd,inint 
end module

program IMC_30
!
!  Construction of effective potentials for monoatomic systems from
!  given RDF
!
!  The program performs one iteration of the inverse 
!  Monte Carlo scheme, which was originally suggested in:
!  A.P.Lyubartsev, A.Laaksonen, Phys.Rev.E, v.52(4) 3730 (1995)
!
!  Relevant reference is also:
!  A.P.Lyubartsev, A.Laaksonen, Comp.Phys.Comm,v.121-122, 57 (1999) 
!  
!  Version 3.0 
!    
use invmc
implicit none
real*8,allocatable :: CROSS(:,:),DIFF(:)
real*8,allocatable :: COR(:),SHELV(:)
integer*4,allocatable :: IUCMP(:),IPVT(:)
character,allocatable :: LNAME(:)*8,NMS(:)*4
integer*4 kode
character*20 FILRDF,FILPOT,FOUT,FDMP
character LAB*12,LABEL*12,DATE*8,xyz*64
logical*1 LPOT,LDMP,LRST
!real*8,parameter :: EPS0=8.854187817620d-12,ELCH=1.60217656535d-19,AVAG=6.0221417930d23,BOLTZ=1.380648813d-23
real*8,parameter :: EPS0=8.854d-12,ELCH=1.602d-19,AVAG=6.02252d23,BOLTZ=1.38054d-23 ! debug
integer*4 NTYP,NMKS,IOUT,IAV,I,IAC,IC,IMKU,INFO,IP,IPT,ISTP,IT,IT1,IT2,ITYP,IUD,J,JC,JT,JTYP,K,NAP,NMKSF,NR,NTPP,NUR
real*8 REGP,DPOTM,RTM,EPS,TEMP,CHI,CP,CRC,CRR,DEE,DIFRC,DLR,DX,DY,DZ,FAC,FELC,FELR,FNR,OSM,POTEN,POTNEW,PRES,RDFC,RDFP
real*8 RDFINC,RDFREF,RR,RRN,RRO,SHIFT,X1,Y1,Z1,W
!  input
namelist /INPUT/ NTYP,NA,NSPEC,NMKS,NMKS0,LPOT,FILRDF,FILPOT,FOUT,FDMP,AF,FQ,BOXL,DR,IOUT,CH,IAV,IPRINT,REGP,DPOTM,LDMP,LRST,RTM,EPS,TEMP
LABEL = 'IMC v3.00'
DATE  = '27-06-2011'
DR    = 5d0
IAV   = 10 
REGP  = 1d0
DPOTM = 5d0 
RTM   = 10d0 
AF    = 3d0      !  erfc(AF)=0.
FQ    = 9d0      !  exp(-FQ**2)=0.
EPS   = 78.3d0
TEMP  = 300d0

! Read INPUT
read(*,INPUT)
write(*,*)' Starting... read data...'

if(NTYP.gt.NTYPES) stop ' increase NTYPES'

NPOTM=NA*NTYP*(NTYP+1)/2
!NAMAX=NA

! allocate
allocate (RDF(NPOTM),POT(NA,NTYP,NTYP),RAS(NA)) ! RDFs
allocate (IADR(0:NTYP),IPOT(NA,NTYP,NTYP),IND(NA,NTYP,NTYP),INA(NPOTM),ITYP1(NPOTM),ITYP2(NPOTM)) ! INT
allocate (CORS(NPOTM),CORP(NPOTM,NPOTM)) ! CORR
allocate (CROSS(NPOTM,NPOTM),DIFF(NPOTM))
allocate (COR(NPOTM),SHELV(NA))
allocate (IUCMP(NPOTM),IPVT(NPOTM))
allocate (LNAME(NTYP),NMS(NTYP))

! calc NTOT
NTOT=0
do ITYP=1,NTYP
  NTOT=NTOT+NSPEC(ITYP)
enddo

!allocate the rest
allocate (ITYPE(NTOT)) !INT
allocate (X(NTOT),Y(NTOT),Z(NTOT),Q(NTOT)) ! COORD

! randomize seed using clock
call init_random_seed() 

!  calculation of electrostatic constant 
VOL=BOXL**3
HBOXL=0.5d0*BOXL
if(DR.gt.HBOXL)DR=HBOXL
COULF = 1d10*ELCH**2/(4d0*PI*EPS0*EPS*BOLTZ*TEMP)
write(*,*)' COULF=',COULF

!  addresses
IADR(0)=0
NOP=0
I=0 
do ITYP=1,NTYP
  do J=1,NSPEC(ITYP)
    I=I+1
    ITYPE(I)=ITYP
    Q(I)=CH(ITYP) 
  enddo
  NOP=NOP+NSPEC(ITYP)
  IADR(ITYP)=IADR(ITYP-1)+NSPEC(ITYP)
enddo

J=0
do I=1,NA
  do ITYP=1,NTYP
    do JTYP=1,ITYP
      J=J+1
      IPOT(I,ITYP,JTYP)=J
      IPOT(I,JTYP,ITYP)=J
      IND(I,ITYP,JTYP)=-1
      IND(I,JTYP,ITYP)=-1
      INA(J)=I
      ITYP1(J)=ITYP
      ITYP2(J)=JTYP
    enddo
  enddo
enddo
NPOT=J

!    Input of potential
if(LPOT)then
  open(unit=3,file=FILPOT,status='old')
  read(3,*)NTPP,NAP,RMIN,RMAX
  if(NAP.ne.NA)stop 'wrong NA in pot.file'
  if(NTPP.ne.NTYP)stop 'wrong NTYP in pot.file'
  do IT1=1,NTYP
    do IT2=IT1,NTYP
      do NR=1,NA 
        read(3,*)RAS(NR),POT(NR,IT1,IT2)
        POT(NR,IT2,IT1)=POT(NR,IT1,IT2)
        if(IPRINT.ge.7) write(*,*) RAS(NR),POT(NR,IT1,IT2),NR,IT1,IT2
      enddo
    enddo
  enddo
  close(3)   
  write(*,*)' potential is taken from ',FILPOT
endif

!    First approximation - from RDF 
open(unit=4,file=FILRDF,status='old')
read(4,*)NTPP,NAP,RMIN,RMAX 
RDFINC=(RMAX-RMIN)/NAP
do I=1,NTPP
  read(4,*)LNAME(I)
  NMS(I)=LNAME(I)(1:4)
enddo
if(NAP.ne.NA)stop 'wrong NA in pot.file'
if(NTPP.ne.NTYP)stop 'wrong NTYP in pot.file'

read(4,*,IOSTAT=kode)RR,RDFP,IT1,IT2
do while(kode.eq.0)
  NR=(RR-RMIN)*NA/(RMAX-RMIN)+1
  if (NR.ge.1.and.NR.le.NA) then 
    IP=IPOT(NR,IT1,IT2)
    RDF(IP)=RDFP
    SHIFT=CH(IT1)*CH(IT2)*COULF/RMAX
    DLR=-dlog(RDFP)                
    if(RDFP.gt.0d0)then
      if(.NOT.LPOT)POT(NR,IT1,IT2)=SHIFT+DLR
      IND(NR,IT1,IT2)=1
      if(.not.LPOT)POT(NR,IT2,IT1)=SHIFT+DLR
      IND(NR,IT2,IT1)=1
    else
      POT(NR,IT1,IT2)=500d0+100d0*(NA-NR)
      IND(NR,IT1,IT2)=0
      POT(NR,IT2,IT1)=500d0+100d0*(NA-NR)
      IND(NR,IT2,IT1)=0
    endif
  endif
  read(4,*,IOSTAT=kode)RR,RDFP,IT1,IT2
enddo
write(*,*)' RDFs were loaded from ',FILRDF 
if(.not.LPOT)write(*,*)' Mean force potential is used'
do NR=1,NA
  RAS(NR)=RMIN+(dfloat(NR)-1d0)*(RMAX-RMIN)/dfloat(NA) 
  SHELV(NR)=4d0*PI*RDFINC*(RAS(NR)**2+RDFINC**2/12d0) 
enddo
if(IPRINT.ge.5)write(*,*) ' N1  N2      r       rdf          pot       ind'
do IT=1,NTYP
  do JT=1,IT
    do NR=1,NA
      if(IND(NR,IT,JT).eq.-1)then
        IP=IPOT(NR,IT,JT)
        POT(NR,IT,JT)=500d0+100d0*(NA-NR)
        IND(NR,IT,JT)=0
        POT(NR,JT,IT)=500d0+100d0*(NA-NR)
        IND(NR,JT,IT)=0 
        RDF(IP)=0d0
      endif
    IP=IPOT(NR,IT,JT)
    if(IPRINT.ge.6) write(*,'(2i5,3f12.4,I4)')IT,JT,RAS(NR),RDF(IP),POT(NR,IT,JT),IND(NR,IT,JT)
    enddo
  enddo
enddo  
RCUT=RMAX-RDFINC 
write(*,*)' Cut off ',RCUT

! Allocate more
call calcnkvm()
allocate (KX(NKVM),KY(NKVM),KZ(NKVM),RKV(NKVM),SSIN(NKVM),SCOS(NKVM),DDSIN(NKVM),DDCOS(NKVM)) ! FUR
allocate (STAB(NTOT,NKVM),CTAB(NTOT,NKVM))

 
! Read Restart if any
if(LRST)then  
  open(unit=45,file=FDMP,status='old',form='unformatted',IOSTAT=kode)
  if (kode.ne.0) then 
    write(*,*) ' file ',FDMP,' not found '
    stop
  endif
  read(45)LAB
  if(LAB.ne.LABEL)then
    write(*,*)'!!! Wrong restart file !!! Label =',LAB
  endif
  read(45)NAV,NMKSF,NMKS0,IMKU,NPOT,IMU,IUD,ENERS,AVEXP,AVE2,VIRS,VIRE
  read(45)(X(I),Y(I),Z(I),I=1,NOP)
  read(45)(CORS(I),I=1,NPOT)
  do I=1,NPOT
    read(45)(CORP(I,J),J=I,NPOT)
    if(IPRINT.ge.8)write(*,'(601I4)')I,(CORP(K,I),K=1,I)
  enddo                          
  close(45)
  write(*,*)'restart file read from ',FDMP
  call ENERGY
  write(*,*)' Restart energy = ',ENER
else
  write(*,*) 'No restart file'
endif

!write(*,'(A$)') 'output xyz filename: '
!read(*,*) xyz
xyz='restart.xyz'
open(unit=77,file=xyz)
write(77,*) NOP
write(77,*) 
do J=1,NOP
  X(J) = X(J) - BOXL * ININT(X(J)/BOXL)
  Y(J) = Y(J) - BOXL * ININT(Y(J)/BOXL)
  Z(J) = Z(J) - BOXL * ININT(Z(J)/BOXL)
  write(77,*) NMS(ITYPE(J)),X(J),Y(J),Z(J)
enddo
close(77)
write(*,*) 'Normal termination'

end program

subroutine calcnkvm()
use invmc
implicit none
integer*4 IKX,IKY,IKZ,IKV,KY1,KY2,KMAXX
real*8    CUTK2,RK2

ALPHA=AF/RCUT

if(ALPHA.gt.0.001)then
  PIAL=(PI/ALPHA)**2
  CUTK2=FQ/PIAL
  KMAXX=dsqrt(CUTK2)*BOXL+1
else
  KMAXX=0
endif

IKV=0
do IKZ=-KMAXX,KMAXX
  do IKX=0,KMAXX
    if(IKX.eq.0)then
      if(IKZ.gt.0)then
        KY1=0
      else
        KY1=1
      endif
      KY2=KMAXX
    else
      KY1=-KMAXX
      KY2=KMAXX-IKX
    endif
    do IKY=KY1,KY2
      RK2=(IKX/BOXL)**2+(IKY/BOXL)**2+(IKZ/BOXL)**2
      if(RK2.le.CUTK2) IKV=IKV+1
    enddo
  enddo
enddo
NKVM=IKV
end subroutine
!   
!=============== EWALD =======================================
!
SUBROUTINE EWALD(QPE)
use invmc
implicit none
integer*4,save :: INIT,IKX,IKY,IKZ,IKV,KY1,KY2,KMAXX,I,KXV,KYV,KZV
real*8,save    :: ESF,CUTK2,RK2,CFUR,AFF,SCS,SSN,SC,SS1,CC1
real*8         :: QPE
data INIT/0/
! save 
! Initialisation 
if(INIT.eq.0)then
  INIT=1
  HBOXL = 0.5d0*BOXL
  VOL = BOXL**3
  ALPHA=AF/RCUT
  if(ALPHA.gt.0.001)then
    PIAL=(PI/ALPHA)**2
    CUTK2=FQ/PIAL
    KMAXX=dsqrt(CUTK2)*BOXL+1
  else
    PIAL=1d0
    CUTK2=0d0
    KMAXX=0
    write(*,*)' no Ewald summation '
  endif
  FK2=PI/BOXL
  FK=2d0*FK2
  IKV=0
  CFUR=COULF/(VOL*PI)
  do IKZ=-KMAXX,KMAXX
    do IKX=0,KMAXX
      if(IKX.eq.0)then
        if(IKZ.gt.0)then
          KY1=0
        else
          KY1=1
        endif
        KY2=KMAXX
      else
        KY1=-KMAXX
        KY2=KMAXX-IKX
      endif
      do IKY=KY1,KY2
        RK2=(IKX/BOXL)**2+(IKY/BOXL)**2+(IKZ/BOXL)**2
        if(RK2.le.CUTK2)then
          IKV=IKV+1
          if(IKV.gt.NKVM)stop 'increase NKVM'
          KX(IKV)=IKX
          KY(IKV)=IKY
          KZ(IKV)=IKZ
          RKV(IKV)=CFUR*dexp(-RK2*PIAL)/RK2
        endif
      enddo
    enddo
  enddo      
  NKV=IKV
  write(*,*)'Kmax=',KMAXX,'    Num. of k-vectors ',NKV
  ESF=0d0
  AFF=ALPHA/dsqrt(PI)
  do I=1,NOP
    ESF=ESF-AFF*COULF*Q(I)**2
  enddo 
endif
QPE = ESF
!   reciprocal space Evald
do IKV=1,NKV
  KXV=KX(IKV)
  KYV=KY(IKV)
  KZV=KZ(IKV)
  SCS=0d0
  SSN=0d0
!$DIR NO_RECURRENCE
  do I=1,NOP
    SC=(KXV*X(I)+KYV*Y(I)+KZV*Z(I))*FK 
    SS1=Q(I)*dsin(SC)
    CC1=Q(I)*dcos(SC)
    SSN=SSN+SS1
    SCS=SCS+CC1 
    STAB(I,IKV)=SS1
    CTAB(I,IKV)=CC1
  enddo
  SSIN(IKV)=SSN
  SCOS(IKV)=SCS
  QPE=QPE+RKV(IKV)*(SSN**2+SCS**2)
enddo
end subroutine
!
!========================= DIFEW ================================
!       
subroutine DIFEW(X1,Y1,Z1,I,CHI,DEE)
use invmc
implicit none
integer*4 IKV,IKX,IKY,IKZ,I
real*8 DEE,CHI,X1,Y1,Z1,DDD,DCS,SCN,SCO,SCP,SINR,DSN
!   reciprocal space Evald
DEE=0d0
CHI=2d0*CHI 
do IKV=1,NKV
  IKX=KX(IKV)
  IKY=KY(IKV)
  IKZ=KZ(IKV)
  SCO=(IKX*X(I)+IKY*Y(I)+IKZ*Z(I))*FK2
  SCN=(IKX*X1  +IKY*Y1  +IKZ*Z1  )*FK2 
  SCP=SCO+SCN
  SINR=dsin(SCN-SCO)*CHI
  DSN=dcos(SCP)*SINR
  DCS=-dsin(SCP)*SINR
  DDSIN(IKV)=DSN
  DDCOS(IKV)=DCS
  DDD=(2d0*SSIN(IKV)+DSN)*DSN+(2d0*SCOS(IKV)+DCS)*DCS  
  DEE=DEE+DDD*RKV(IKV)
enddo
end subroutine 
!
!=================== RECSIN ========================================
!
subroutine RECSIN
use invmc
implicit none
integer*4 IKV

do IKV=1,NKV
  SSIN(IKV)=SSIN(IKV)+DDSIN(IKV)
  SCOS(IKV)=SCOS(IKV)+DDCOS(IKV)
enddo      
end subroutine
! 
!===================== ENERGY ========================================
!     
subroutine ENERGY
use invmc
implicit none
integer*4 I,J,IT,JT,NR
real*8 DEE,DDE,DX,DY,DZ,ETT,RR
call EWALD(EFUR) 
write(*,*)' Efur=',EFUR
DE=0d0
do I=2,NOP   
  IT=ITYPE(I)
  do J=1,I-1 
    JT=ITYPE(J) 
    DX=X(I)-X(J)
    DY=Y(I)-Y(J)
    DZ=Z(I)-Z(J)      
    call PBC(DX,DY,DZ)
    RR=dsqrt(DX**2+DY**2+DZ**2)
    if(RR.le.RCUT)then
      NR=(RR-RMIN)*NA/(RMAX-RMIN)+1
      DDE=POT(NR,IT,JT)
      DEE=CH(IT)*CH(JT)*COULF*(derfc(ALPHA*RR)-1d0)/RR
      EFUR=EFUR+DEE 
      ETT=DDE+DEE
      DE =DE +ETT
    endif
  enddo
enddo 
ENER=EFUR+DE
end subroutine
! 
!===================== CORREL ========================================
!     
subroutine CORREL                                                                
use invmc
implicit none
real*8 CCOR(NPOTM),EIT(NTOT),ENEW,CCC,DDE,DEE,DX,DY,DZ,EDL,ETT,RR
integer*4 I,J,IT,JT,NR,NR1,NR2,IPT

call EWALD(ENEW) 
EFUR=ENEW
DE=0d0     
VIR=0d0
do I=1,NPOT
  CCOR(I)=0d0
enddo
do I=1,NOP
  EIT(I)=0d0
enddo  
do I=2,NOP   
  IT=ITYPE(I) 
  EDL=0d0
  do J=1,I-1 
    JT=ITYPE(J)
    DX=X(I)-X(J)
    DY=Y(I)-Y(J)
    DZ=Z(I)-Z(J)
    call PBC(DX,DY,DZ)
    RR=dsqrt(DX**2+DY**2+DZ**2)
    if(RR.le.RCUT)then
      NR=(RR-RMIN)*NA/(RMAX-RMIN)+1
      if(IND(NR,IT,JT).eq.0.and.ISTEP.gt.NMKS0)then
        if (IPRINT.gt.5) then
          write(*,*)' !!! Wrong configuration:'
          write(*,'(I5,3f10.3)')I,X(I),Y(I),Z(I)
          write(*,'(I5,3f10.3)')J,X(J),Y(J),Z(J)
          write(*,'(I5,a,f10.3)')NR,' R= ',RR
        endif
      endif
      IPT=IPOT(NR,IT,JT)
      DDE=POT(NR,IT,JT)
      NR1=(RR-RMIN)*NA/(RMAX-RMIN)+0.5d0
      if(NR1.le.0)NR1=1 
      if(IND(NR1,IT,JT).eq.0)NR1=NR1+1
      NR2=NR1+1         
      if(NR2.le.NA)then
        VIR=VIR+RR*(POT(NR2,IT,JT)-POT(NR1,IT,JT))/(NR2-NR1)        
      else
        if (IPRINT.gt.5) write(*,*)' abnormal virial',I,J,IT,JT,NR1,NR2,RR 
      endif
      CCOR(IPT)=CCOR(IPT)+1d0
      DEE=CH(IT)*CH(JT)*COULF*(derfc(ALPHA*RR)-1d0)/RR
      EFUR=EFUR+DEE 
      ETT=DDE+DEE
      DE =DE +ETT
      EIT(I)=EIT(I)+ETT
      EIT(J)=EIT(J)+ETT
    endif
  enddo
enddo  
VIR=VIR*NA/((RMAX-RMIN)*3d0)
ENEW=ENEW+DE
if(abs(ENER-ENEW).gt.0.01)write(*,'(I7,a,f12.4,a,f12.4)') ISTEP,' En.error.  New:',ENEW,'    old:',ENER
ENER=ENEW
!   correlators                       
do I=1,NPOT
  CCC=CCOR(I)
  if(CCC.ne.0d0)then
    CORS(I)=CORS(I)+CCC
    do J=I,NPOT
      CORP(I,J)=CORP(I,J)+CCC*CCOR(J)
    enddo
  endif
enddo 
ENERS=ENERS+ENER
NAV=NAV+1
end subroutine
!
!======================= PBC ========================================
!
subroutine PBC(XX,YY,ZZ)
use invmc
implicit none
real*8 XX,YY,ZZ

if(XX.gt. HBOXL)XX=XX-BOXL
if(XX.lt.-HBOXL)XX=XX+BOXL
if(YY.gt. HBOXL)YY=YY-BOXL
if(YY.lt.-HBOXL)YY=YY+BOXL
if(ZZ.gt. HBOXL)ZZ=ZZ-BOXL
if(ZZ.lt.-HBOXL)ZZ=ZZ+BOXL
return
end subroutine
!
!================= RANDOM =======================
! 
function rnd()
implicit none
real*8 a,rnd
call random_number(a)
rnd=a
end function

SUBROUTINE init_random_seed()
implicit none
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE

function iint(num)
implicit none
integer iint
real*8 num
if (num.ge.0D0) then
  iint=int(num)
else
  if ((num-int(num)).eq.0D0) then
    iint=int(num)
  else
    iint=int(num)-1
  endif
endif
end function

function inint(num)
implicit none
integer iint
real*8 inint,num
inint=dble(iint(num+0.5))
end function

