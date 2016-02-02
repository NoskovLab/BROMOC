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

! Fortran common block for Grand Canonical Monte Carlo
! Dimension of arrays
module grandmod
implicit none
integer   dtype, datom, mxcnt
parameter (datom = 23000, dtype = 1000) ! if datom bigger than 23000 start giving segm fault in ifort, solution -heap-arrays 0
integer   dbuff
parameter (dbuff = 2*dtype)
!integer   dindx
!parameter (dindx = dtype*(dtype+1)/2)
real,external :: rndm,rndm1
real    ikbt,ikbtdna
real    temp,kbt,fact0n1,fact0n2,kbtdt,dnatemp,kbtdna
real*8 runtime
integer ntype, nold, nbuffer, ntot, natom, nfix

integer typei(datom)
real    x(datom),y(datom),z(datom)
real    fx(datom),fy(datom),fz(datom)
integer typtyp(datom),ibuffer(datom)

integer,allocatable :: warn(:)
integer     ibfftyp(dbuff),nat(dtype),idtype(dtype)
integer     nremove(dtype),ninsert(dtype),ntotat(dtype),nwtype(dtype)
real        eps(dtype),sigma(dtype),cg(dtype),diffusion(dtype),ampl1(dtype),p1(2,dtype)
real        ampl2(dtype),p2(2,dtype),rcylinder(dtype),deltaz,ampl3(dtype),p3(dtype)
real        volume(dbuff),avnum(dbuff),mu(dbuff),density(dbuff),kb(dbuff),LZmin(dbuff),LZmax(dbuff),Rmin(dbuff),Rmax(dbuff)
character*4 atnam(dtype), atnam2(dtype)
logical*1   Qbufferbias(dtype)

integer,allocatable :: nforward(:,:),nbackward(:,:)  
real,allocatable :: zcont(:)
integer nprint,nanal,cntpts,svcntfq
integer igr, nframe
real    dt,mcmax,bdmax
real dids(5,datom),fact2a(dtype),fact1a(dtype),fact2pd,beta,diff0,ibeta,diffcutoff
!real*16 ener,eelec,evdw,ememb,esrpmf,esrpmfmx,eefpot,eefpotmx,econ,erfpar
real ener,eelec,evdw,ememb,esrpmf,esrpmfmx,eefpot,eefpotmx,econ,erfpar
real    rsphe,rsphe2,lx,ly,lz,lx2p,ly2p,lz2p,lx2m,ly2m,lz2m,cx,cy,cz,iecx,iecy
real, allocatable :: epp4(:),sgp2(:)
real, allocatable :: c0(:),c1(:),c2(:),c3(:),c4(:)
real    cdie,rth,srpx,srpk,srpy,tvol
real    afact,kappa,ikappa
!integer typat(datom) 
logical*1 Qsphere, Qecyl, Qbox
logical*1 Qenergy, Qforces, Qnobond, Qnonbond, Qmemb, Qgr, Qrho, Qrdna, Qprob, Qdiffuse, Qsrpmf, Qionpair, Qenerprofile, Qprofile, Qsec, Qpres, Qpore, Qwarn, Qcountion, Qproxdiff 
logical*1,allocatable :: Qlj(:),Qsrpmfi(:),Qefpot(:)

! membrane parameters
real    voltage, thick2, zmemb, plength2, pcenter
real    czmax, czmin
real  tmemb, epsm, zmemb1,zmemb2, ceps

contains

function outbox(xi,yi,zi)
implicit none
real xi,yi,zi
logical*1 outbox

if (Qsphere) then
  if (((xi-cx)**2+(yi-cy)**2+(zi-cz)**2).gt.Rsphe2) then
    outbox=.true.
  else
    outbox=.false.
  endif 
elseif (Qecyl) then
  if ((((xi-cx)*iecx)**2+((yi-cy)*iecy)**2).gt.1.0.or.zi.lt.lz2m.or.zi.gt.lz2p) then  
    outbox=.true.
  else
    outbox=.false.
  endif
else
  if (xi.lt.lx2m.or.xi.gt.lx2p.or.yi.lt.ly2m.or.yi.gt.ly2p.or.zi.lt.lz2m.or.zi.gt.lz2p) then
    outbox=.true.
  else
    outbox=.false.
  endif
endif
end function
end module



