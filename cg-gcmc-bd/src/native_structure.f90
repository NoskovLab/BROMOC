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

subroutine native_structure
!This routine obtain the site interaction coordinates in
!(B,10_1,3.38)-DNA conformation as well as the natural bond
!lenghts, bond angles and dihedral angles. (xnat,ynat,znat) are
!cartesian coordinates and (rnat,phinat,znat) are cylindrical
!polar coordinates. The xnat, ynat, znat, and rnat coordinates
!are in Angstroms and phinat is in degrees.
use constamod
use stdiomod
use grandmod
use nucleotmod
use errormod
!Local variables
implicit none
integer isite, nuc, i,j,k
integer jbase, jS, jP, jPE, jPF
real  cte,icte
real  angAb, angTb, angCb, angGb
real  angS, angP
real  xAb(3), yAb(3), zAb(3)
real  xTb(3), yTb(3), zTb(3)
real  xCb(3), yCb(3), zCb(3)
real  xGb(3), yGb(3), zGb(3)
real  xxS(3), yyS(3), zzS(3)
real  xP(3), yP(3), zP(3)
real  vecR(3,3)
real  angbond, angdih
real  phPSAbd, phPSTbd, phPSCbd, phPSGbd, phPSAb2d, phPSTb2d,phPSCb2d, phPSGb2d, phSPSd, phPSPd
real  dhAbSPSd, dhTbSPSd, dhCbSPSd, dhGbSPSd, dhSPSAbd,dhSPSTbd, dhSPSCbd, dhSPSGbd, dhSPSPd, dhPSPSd
real  cylAb(3,2), cylTb(3,2), cylCb(3,2), cylGb(3,2)
real  cylS(3,2), cylP(3,2)

if (istrs.eq.0) return

cte = pi/180.0
icte = 1.0/cte
nuc=3*inuc
jPE=1
jPF=1
if (.not.QfirstP) then
  nuc=3*inuc-1
  if (Qinvstr) then
    jPF=0
  else
    jPE=0
  endif
endif

!Cylindrical polar coordinates
!First strand
!Bases
cylAb(1,1) = cylall(1,1)
cylAb(2,1) = cylall(1,2)
cylAb(3,1) = cylall(1,3)
cylTb(1,1) = cylall(2,1)
cylTb(2,1) = cylall(2,2)
cylTb(3,1) = cylall(2,3)
cylCb(1,1) = cylall(3,1)
cylCb(2,1) = cylall(3,2)
cylCb(3,1) = cylall(3,3)
cylGb(1,1) = cylall(4,1)
cylGb(2,1) = cylall(4,2)
cylGb(3,1) = cylall(4,3)
!Sugar
cylS(1,1) = cylall(5,1)
cylS(2,1) = cylall(5,2)
cylS(3,1) = cylall(5,3)
!Phosphate
cylP(1,1) = cylall(6,1) 
cylP(2,1) = cylall(6,2)
cylP(3,1) = cylall(6,3)
if (Qinvstr) then
  !Bases
  cylAb(2:3,1) = -cylAb(2:3,1)
  cylTb(2:3,1) = -cylTb(2:3,1)
  cylCb(2:3,1) = -cylCb(2:3,1)
  cylGb(2:3,1) = -cylGb(2:3,1)
  !Sugar
  cylS(2:3,1) = -cylS(2:3,1)
  !Phosphate
  cylP(2:3,1) = -cylP(2:3,1)
endif 
!Second strand
if (istrs.eq.2) then
  cylAb(1,2)   = cylAb(1,1)
  cylTb(1,2)   = cylTb(1,1)
  cylCb(1,2)   = cylCb(1,1)
  cylGb(1,2)   = cylGb(1,1)
  cylS (1,2)   = cylS(1,1)
  cylP (1,2)   = cylP(1,1)
  cylAb(2:3,2) = -cylAb(2:3,1)
  cylTb(2:3,2) = -cylTb(2:3,1)
  cylCb(2:3,2) = -cylCb(2:3,1)
  cylGb(2:3,2) = -cylGb(2:3,1)
  cylS (2:3,2) = -cylS(2:3,1)
  cylP (2:3,2) = -cylP(2:3,1)
endif

!OBTENTION OF NATURAL BONDS, BOND ANGLES AND DIHEDRAL ANGLES
do i = 1, 3
  angAb  = (cylAb(2,1)+(i-1)*ain)*cte
  angTb  = (cylTb(2,1)+(i-1)*ain)*cte
  angCb  = (cylCb(2,1)+(i-1)*ain)*cte
  angGb  = (cylGb(2,1)+(i-1)*ain)*cte
  angS   = (cylS(2,1)+(i-1)*ain)*cte
  angP   = (cylP(2,1)+(i-1)*ain)*cte
  xAb(i) = cylAb(1,1)*cos(angAb)
  yAb(i) = cylAb(1,1)*sin(angAb)
  zAb(i) = cylAb(3,1) + (i-1)*din
  xTb(i) = cylTb(1,1)*cos(angTb)
  yTb(i) = cylTb(1,1)*sin(angTb)
  zTb(i) = cylTb(3,1) + (i-1)*din
  xCb(i) = cylCb(1,1)*cos(angCb)
  yCb(i) = cylCb(1,1)*sin(angCb)
  zCb(i) = cylCb(3,1) + (i-1)*din
  xGb(i) = cylGb(1,1)*cos(angGb)
  yGb(i) = cylGb(1,1)*sin(angGb)
  zGb(i) = cylGb(3,1) + (i-1)*din
  xxS(i) = cylS(1,1)*cos(angS)
  yyS(i) = cylS(1,1)*sin(angS)
  zzS(i) = cylS(3,1) + (i-1)*din
  xP(i)  = cylP(1,1)*cos(angP)
  yP(i)  = cylP(1,1)*sin(angP)
  zP(i)  = cylP(3,1) + (i-1)*din
enddo
if (Qinvstr) then
  i=1
  j=2
  k=3
else
  i=3
  j=2
  k=1
endif 
!Natural bond lenghts
dSAb = sqrt((xAb(i)-xxS(i))**2+(yAb(i)-yyS(i))**2+(zAb(i)-zzS(i))**2)
dSTb = sqrt((xTb(i)-xxS(i))**2+(yTb(i)-yyS(i))**2+(zTb(i)-zzS(i))**2)
dSCb = sqrt((xCb(i)-xxS(i))**2+(yCb(i)-yyS(i))**2+(zCb(i)-zzS(i))**2)
dSGb = sqrt((xGb(i)-xxS(i))**2+(yGb(i)-yyS(i))**2+(zGb(i)-zzS(i))**2)
dPS5 = sqrt(( xP(i)-xxS(i))**2+ (yP(i)-yyS(i))**2+( zP(i)-zzS(i))**2)  
dPS3 = sqrt(( xP(i)-xxS(j))**2+ (yP(i)-yyS(j))**2+( zP(i)-zzS(j))**2)
!Natural bond angles
vecR(1,1) = xAb(i)-xxS(i)
vecR(2,1) = yAb(i)-yyS(i)
vecR(3,1) = zAb(i)-zzS(i)
vecR(1,2) =  xP(i)-xxS(i)
vecR(2,2) =  yP(i)-yyS(i)
vecR(3,2) =  zP(i)-zzS(i)
phPSAb  = angbond(vecR)
vecR(1,1) = xTb(i)-xxS(i)
vecR(2,1) = yTb(i)-yyS(i)
vecR(3,1) = zTb(i)-zzS(i)
phPSTb  = angbond(vecR)
vecR(1,1) = xCb(i)-xxS(i)
vecR(2,1) = yCb(i)-yyS(i)
vecR(3,1) = zCb(i)-zzS(i)
phPSCb  = angbond(vecR)
vecR(1,1) = xGb(i)-xxS(i)
vecR(2,1) = yGb(i)-yyS(i)
vecR(3,1) = zGb(i)-zzS(i)
phPSGb  = angbond(vecR)
vecR(1,1) = xAb(j)-xxS(j)
vecR(2,1) = yAb(j)-yyS(j)
vecR(3,1) = zAb(j)-zzS(j)
vecR(1,2) = xP(i)-xxS(j)
vecR(2,2) = yP(i)-yyS(j)
vecR(3,2) = zP(i)-zzS(j)
phPSAb2 = angbond(vecR)
vecR(1,1) = xTb(j)-xxS(j)
vecR(2,1) = yTb(j)-yyS(j)
vecR(3,1) = zTb(j)-zzS(j)
phPSTb2 = angbond(vecR)
vecR(1,1) = xCb(j)-xxS(j)
vecR(2,1) = yCb(j)-yyS(j)
vecR(3,1) = zCb(j)-zzS(j)
phPSCb2 = angbond(vecR)
vecR(1,1) = xGb(j)-xxS(j)
vecR(2,1) = yGb(j)-yyS(j)
vecR(3,1) = zGb(j)-zzS(j)
phPSGb2 = angbond(vecR)
vecR(1,1) = xxS(i)-xP(i)
vecR(2,1) = yyS(i)-yP(i) 
vecR(3,1) = zzS(i)-zP(i)
vecR(1,2) = xxS(j)-xP(i)
vecR(2,2) = yyS(j)-yP(i)
vecR(3,2) = zzS(j)-zP(i)
phSPS   = angbond(vecR)
vecR(1,1) = xP(i)-xxS(j)
vecR(2,1) = yP(i)-yyS(j)
vecR(3,1) = zP(i)-zzS(j)
vecR(1,2) = xP(j)-xxS(j)
vecR(2,2) = yP(j)-yyS(j)
vecR(3,2) = zP(j)-zzS(j)
phPSP   = angbond(vecR)
!Natural dihedral angles
vecR(1,1) = xxS(i)-xAb(i)
vecR(2,1) = yyS(i)-yAb(i)
vecR(3,1) = zzS(i)-zAb(i)
vecR(1,2) =  xP(i)-xxS(i)
vecR(2,2) =  yP(i)-yyS(i)
vecR(3,2) =  zP(i)-zzS(i)
vecR(1,3) = xxS(j)- xP(i)
vecR(2,3) = yyS(j)- yP(i)
vecR(3,3) = zzS(j)- zP(i)
dhAbSPS = angdih(vecR)
vecR(1,1) = xxS(i)-xTb(i)
vecR(2,1) = yyS(i)-yTb(i)
vecR(3,1) = zzS(i)-zTb(i)
dhTbSPS = angdih(vecR)
vecR(1,1) = xxS(i)-xCb(i)
vecR(2,1) = yyS(i)-yCb(i)
vecR(3,1) = zzS(i)-zCb(i)
dhCbSPS = angdih(vecR)
vecR(1,1) = xxS(i)-xGb(i)
vecR(2,1) = yyS(i)-yGb(i)
vecR(3,1) = zzS(i)-zGb(i)
dhGbSPS = angdih(vecR)
vecR(1,1) = vecR(1,2)
vecR(2,1) = vecR(2,2)
vecR(3,1) = vecR(3,2)
vecR(1,2) = vecR(1,3)
vecR(2,2) = vecR(2,3)
vecR(3,2) = vecR(3,3)
vecR(1,3) = xAb(j)-xxS(j)
vecR(2,3) = yAb(j)-yyS(j)  
vecR(3,3) = zAb(j)-zzS(j)
dhSPSAb = angdih(vecR)
vecR(1,3) = xTb(j)-xxS(j)
vecR(2,3) = yTb(j)-yyS(j)
vecR(3,3) = zTb(j)-zzS(j)
dhSPSTb = angdih(vecR)
vecR(1,3) = xCb(j)-xxS(j)
vecR(2,3) = yCb(j)-yyS(j)
vecR(3,3) = zCb(j)-zzS(j)
dhSPSCb = angdih(vecR)
vecR(1,3) = xGb(j)-xxS(j)
vecR(2,3) = yGb(j)-yyS(j)
vecR(3,3) = zGb(j)-zzS(j)
dhSPSGb = angdih(vecR)
vecR(1,3) = xP(j)-xxS(j)
vecR(2,3) = yP(j)-yyS(j)
vecR(3,3) = zP(j)-zzS(j)
dhSPSP = angdih(vecR)
vecR(1,1) = vecR(1,2)
vecR(2,1) = vecR(2,2)
vecR(3,1) = vecR(3,2)
vecR(1,2) = vecR(1,3)
vecR(2,2) = vecR(2,3)
vecR(3,2) = vecR(3,3)
vecR(1,3) = xxS(k)-xP(j)
vecR(2,3) = yyS(k)-yP(j)      
vecR(3,3) = zzS(k)-zP(j)
dhPSPS = angdih(vecR)

!NATIVE STRUCTURE FOR B ISOFORM OF DNA
!First strand  (5'-3' direction)
jbase = 0
jS    = 0
jP    = 0
do isite = 1, nuc
  if (namsite(isite).eq.'Ab') then
    jbase = jbase + 1
    rnat(isite)   = cylAb(1,1)
    phinat(isite) = cylAb(2,1) + (jbase-1)*ain
    znat(isite)   = cylAb(3,1) + (jbase-1)*din
  else if (namsite(isite).eq.'Tb') then
    jbase = jbase + 1
    rnat(isite)   = cylTb(1,1)
    phinat(isite) = cylTb(2,1) + (jbase-1)*ain
    znat(isite)   = cylTb(3,1) + (jbase-1)*din
  else if (namsite(isite).eq.'Cb') then
    jbase = jbase + 1
    rnat(isite)   = cylCb(1,1)
    phinat(isite) = cylCb(2,1) + (jbase-1)*ain
    znat(isite)   = cylCb(3,1) + (jbase-1)*din
  else if (namsite(isite).eq.'Gb') then
    jbase = jbase + 1
    rnat(isite)   = cylGb(1,1)
    phinat(isite) = cylGb(2,1) + (jbase-1)*ain
    znat(isite)   = cylGb(3,1) + (jbase-1)*din
  else if (namsite(isite).eq.'S ') then
    jS = jS + 1
    rnat(isite)   = cylS(1,1)
    phinat(isite) = cylS(2,1) + (jS-1)*ain
    znat(isite)   = cylS(3,1) + (jS-1)*din
  else if (namsite(isite).eq.'P ') then
    jP = jP + 1
    rnat(isite)   = cylP(1,1)
    phinat(isite) = cylP(2,1) + (jP-jPE)*ain
    znat(isite)   = cylP(3,1) + (jP-jPE)*din
  else
    call error ('native_structure', 'wrong name for sites',faterr)
  endif
  xnat(isite) = rnat(isite)*cos(phinat(isite)*cte)
  ynat(isite) = rnat(isite)*sin(phinat(isite)*cte)
enddo
!Second strand (3'-5' direction)
if (istrs.eq.2) then
  jbase = 0
  jS    = 0
  jP    = 0
  do isite = nsites, nuc+1, -1
    if (namsite(isite).eq.'Ab') then
      jbase = jbase + 1
      rnat(isite)   = cylAb(1,2)
      phinat(isite) = cylAb(2,2) + (jbase-1)*ain
      znat(isite)   = cylAb(3,2) + (jbase-1)*din
    else if (namsite(isite).eq.'Tb') then
      jbase = jbase + 1
      rnat(isite)   = cylTb(1,2)
      phinat(isite) = cylTb(2,2) + (jbase-1)*ain
      znat(isite)   = cylTb(3,2) + (jbase-1)*din
    else if (namsite(isite).eq.'Cb') then
      jbase = jbase + 1
      rnat(isite)   = cylCb(1,2)
      phinat(isite) = cylCb(2,2) + (jbase-1)*ain
      znat(isite)   = cylCb(3,2) + (jbase-1)*din
    else if (namsite(isite).eq.'Gb') then
      jbase = jbase + 1
      rnat(isite)   = cylGb(1,2)
      phinat(isite) = cylGb(2,2) + (jbase-1)*ain
      znat(isite)   = cylGb(3,2) + (jbase-1)*din
    else if (namsite(isite).eq.'S ') then
      jS = jS + 1
      rnat(isite)   = cylS(1,2)
      phinat(isite) = cylS(2,2) + (jS-1)*ain
      znat(isite)   = cylS(3,2) + (jS-1)*din
    else if (namsite(isite).eq.'P ') then
      jP = jP + 1
      rnat(isite)   = cylP(1,2)
      phinat(isite) = cylP(2,2) + (jP-jPF)*ain
      znat(isite)   = cylP(3,2) + (jP-jPF)*din
    else
      call error ('native_structure', 'wrong name for sites',faterr)
    endif
    xnat(isite) = rnat(isite)*cos(phinat(isite)*cte)
    ynat(isite) = rnat(isite)*sin(phinat(isite)*cte)
  enddo
endif

!OUTPUT
phPSAbd = phPSAb*icte
phPSTbd = phPSTb*icte
phPSCbd = phPSCb*icte
phPSGbd = phPSGb*icte
phPSAb2d = phPSAb2*icte
phPSTb2d = phPSTb2*icte
phPSCb2d = phPSCb2*icte
phPSGb2d = phPSGb2*icte
phSPSd   = phSPS*icte
phPSPd   = phPSP*icte
dhAbSPSd = dhAbSPS*icte
dhTbSPSd = dhTbSPS*icte
dhCbSPSd = dhCbSPS*icte
dhGbSPSd = dhGbSPS*icte
dhSPSAbd = dhSPSAb*icte
dhSPSTbd = dhSPSTb*icte
dhSPSCbd = dhSPSCb*icte
dhSPSGbd = dhSPSGb*icte
dhSPSPd  = dhSPSP*icte
dhPSPSd  = dhPSPS*icte
if (Qninfo) then
  write(outu,'(/6x,a)') 'Natural bond lenghts, bond angles and dihedral angles'
  write(outu,'(6x,a)') 'NOTE: distances in Ang. and angles in'//' degrees'
  write(outu,'(/6x,a,3x,a)') 'BOND','NATURAL BOND LENGHT'
  write(outu,'(6x,a,3x,a)') '----','-------------------'
  write(outu,'(6x,a,3x,f6.3)') 'S-Ab',dSAb
  write(outu,'(6x,a,3x,f6.3)') 'S-Tb',dSTb
  write(outu,'(6x,a,3x,f6.3)') 'S-Cb',dSCb
  write(outu,'(6x,a,3x,f6.3)') 'S-Gb',dSGb
  write(outu,'(6x,a,1x,f6.3)') 'S(5)-P',dPS5
  write(outu,'(6x,a,1x,f6.3)') 'S(3)-P',dPS3
  write(outu,'(/6x,a,3x,a)') 'BOND ANGLE','NATURAL BOND ANGLE'
  write(outu,'(6x,a,3x,a)') '----------','------------------'
  write(outu,'(6x,a,3x,f7.3)') 'P-(5)S-Ab',phPSAbd
  write(outu,'(6x,a,3x,f7.3)') 'P-(5)S-Tb',phPSTbd
  write(outu,'(6x,a,3x,f7.3)') 'P-(5)S-Cb',phPSCbd
  write(outu,'(6x,a,3x,f7.3)') 'P-(5)S-Gb',phPSGbd
  write(outu,'(6x,a,3x,f7.3)') 'P-(3)S-Ab',phPSAb2d
  write(outu,'(6x,a,3x,f7.3)') 'P-(3)S-Tb',phPSTb2d
  write(outu,'(6x,a,3x,f7.3)') 'P-(3)S-Cb',phPSCb2d
  write(outu,'(6x,a,3x,f7.3)') 'P-(3)S-Gb',phPSGb2d
  write(outu,'(6x,a,3x,f7.3)') 'S(5)-P-(3)S',phSPSd
  write(outu,'(6x,a,3x,f7.3)') 'P-(3)S(5)-P',phPSPd
  write(outu,'(/6x,a,3x,a)') 'DIHEDRAL ANGLE','NATURAL DIHEDRAL ANGLE'
  write(outu,'(6x,a,3x,a)') '--------------','----------------------'
  write(outu,'(6x,a,3x,f7.3)') 'Ab-S(5)-P-(3)S',dhAbSPSd
  write(outu,'(6x,a,3x,f7.3)') 'Tb-S(5)-P-(3)S',dhTbSPSd
  write(outu,'(6x,a,3x,f7.3)') 'Cb-S(5)-P-(3)S',dhCbSPSd
  write(outu,'(6x,a,3x,f7.3)') 'Gb-S(5)-P-(3)S',dhGbSPSd
  write(outu,'(6x,a,3x,f7.3)') 'S(5)-P-(3)S-Ab',dhSPSAbd
  write(outu,'(6x,a,3x,f7.3)') 'S(5)-P-(3)S-Tb',dhSPSTbd
  write(outu,'(6x,a,3x,f7.3)') 'S(5)-P-(3)S-Cb',dhSPSCbd
  write(outu,'(6x,a,3x,f7.3)') 'S(5)-P-(3)S-Gb',dhSPSGbd
  write(outu,'(6x,a,3x,f8.3)') 'S(5)-P-(3)S(5)-P',dhSPSPd
  write(outu,'(6x,a,3x,f8.3)') 'P-(3)S(5)-P-(3)S',dhPSPSd
endif
write(outu,'(/6x,a)') 'Native structure for B isoform of DNA'
write(outu,'(/6x,a)') 'STRAND SITE# NUCLEOT SITE CARTESIAN COORDINATES (X,Y,Z)'
write(outu,'(6x,a)') '------ ----- ------- ---- -----------------------------'
do i = 1, nsites
  write(outu,'(6x,i2,5x,i4,2x,a1,5x,a2,3x,f8.3,2x,f8.3,2x,f8.3)') strand(i),i,namnucl(i),namsite(i),xnat(i),ynat(i),znat(i)
enddo
return
end subroutine
