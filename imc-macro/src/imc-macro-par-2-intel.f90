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

module math
implicit none
contains
  function inintv(vec)
  real vec(3),inintv(3)
  inintv(1)=inint(vec(1))
  inintv(2)=inint(vec(2))
  inintv(3)=inint(vec(3))
  end function

  function invmat(mat)
  real invmat(3,3),mat(3,3)
  invmat(1,1)=mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
  invmat(2,1)=mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
  invmat(3,1)=mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
  invmat(1,2)=mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
  invmat(2,2)=mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
  invmat(3,2)=mat(3,1)*mat(1,2)-mat(1,1)*mat(3,2)
  invmat(1,3)=mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
  invmat(2,3)=mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
  invmat(3,3)=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
  invmat=invmat/dot_product(mat(1,1:3),invmat(1:3,1))
  end function

  function cross_product(u,v)
  implicit none
  real u(3),v(3),w(3),cross_product(3)
  w(1)=u(2)*v(3)-u(3)*v(2)
  w(2)=u(3)*v(1)-u(1)*v(3)
  w(3)=u(1)*v(2)-u(2)*v(1)
  cross_product=w
  end function

  function inint(num)
  implicit none
  real inint,num
  inint=iint(num+0.5)
  end function

  function iint(num)
  implicit none
  integer iint
  real num
  if (num.ge.0.0) then
    iint=int(num)
  else
    if ((num-int(num)).eq.0.0) then
      iint=int(num)
    else
      iint=int(num)-1
    endif
  endif
  end function

end module

module invmc
!  Declarations
!   nop - max number of particles
!   ntypes - max number of particles types
!   namax - max number of grid points
!   nkvm - max number of k-vectors in reciprocal Ewald
!
implicit none
real,parameter :: pi=3.14159265358979323846264338327950288419716939937510
real,parameter :: pi2=2.0*pi 
real,parameter :: i3=1e0/3e0 
integer nkvm,nfix,nfree,nfxfr
real,allocatable :: x(:),y(:),z(:),q(:)
real,allocatable :: rdf(:),pot(:,:,:),ras(:),ch(:)
integer,allocatable :: itype(:),ipot(:,:,:),ina(:),ityp1(:),ityp2(:),nspec(:),nspecf(:),nspecfr(:)
logical*1,allocatable :: ind(:,:,:)
integer*1,allocatable :: potres(:,:)
real vol,dr,coulf,af,fq,alpha,rcut,rcut2,rmin,rmax,avexp,ave2,vir,virs,vire,iri
real boxlx,boxly,boxlz,iboxlx,iboxly,iboxlz,bv(3,3),invbv(3,3),inop
real*8 de,fce,fcee,ener,eners,efur,esf
integer nop,npot,na,istep,nmks0,nkv,nav,npair,iprint
integer,allocatable :: cors(:),corp(:,:)
integer,allocatable :: kx(:),ky(:),kz(:)
real,allocatable :: rkv(:),ssin(:),scos(:),ddsin(:),ddcos(:)
character label*22,datestamp*10
logical*1 lelec,lpot,box
end module

program imcmacro
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
!  Version 3.8
!    
use ifport  ! use only for intel compiler, comment for gnu compiler
use invmc
use math
implicit none
real,allocatable :: cross(:,:),diff(:)
real,allocatable :: cor(:),shelv(:)
real,allocatable :: ssinl(:),scosl(:),ddsinl(:),ddcosl(:)
integer,allocatable :: iucmp(:),ipvt(:)
character,allocatable :: nms(:)*4,fpn(:)*4
character :: nmst1*4,nmst2*4,keyword*5
integer*8 timer,wall
integer kode,iseed,wxyzfq
integer rstfq
character*256 fdmp
logical*1 ldmp,lrst
character*256 filrdf,filpot,fout,wxyznm,rxyznm,respotnm
character*1024 line
logical*1 loopon,wxyz,rxyz,ldmppot,lzm,lrespot,lseppot,lseprdf,latvec
real,parameter :: eps0=8.854187817620e-12,elch=1.60217656535e-19,avag=6.0221412927e23,boltz=1.380648813e-23
integer ntyp,nmks,iout,iav,i,ic,info,ip,ipt,it,it1,it2,ityp,j,jc,jt,jtyp,nmksf,nr,nur,nap,ntpp
real regp,dpotm,rtm,eps,temp,chi,crc,crr,dee,difrc,dlr,dx,dy,dz,felc,felr,fnr,osm,poten,potnew,pres,rdfc,rdfp
real rdfinc,rdfref,rr,rrn,rrn2,rro,rro2,shift,x1,y1,z1,cent(3),aone,zeromove
real vrvs,vr,vs,vs2,a,b,b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,pos(3)
real*8 deloc,efurl,del,enerl,enersl
integer jobs,tid,nth,isd,iud,iudl,cova,aun,omp_get_num_threads,omp_get_thread_num,istepl,navl
real virel,virsl,virl
integer,allocatable :: corsl(:),corpl(:,:)
real,allocatable :: xl(:),yl(:),zl(:)

integer*1 restyp
!  input
!namelist /input/ nmks,nmks0,lpot,filrdf,filpot,fout,af,fq,b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,dr,iout,iav,iprint,regp,dpotm,rtm,eps,temp,iseed,wxyz,rxyz,wxyznm,rxyznm,wxyzfq,ldmppot,zeromove,lzm,lelec,lrespot,respotnm,lseppot,lseprdf
namelist /input/ nmks,nmks0,lpot,filrdf,filpot,fout,fdmp,af,fq,b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,dr,iout,iav,iprint,regp,dpotm,ldmp,lrst,rtm,eps,temp,iseed,wxyz,rxyz,wxyznm,rxyznm,wxyzfq,rstfq,ldmppot,zeromove,lzm,lelec,lrespot,respotnm,lseppot,lseprdf

label    = 'IMC-MACRO-PAR-2 v3.90'
datestamp = '19-02-2014'
dr       = 5e0                  ! max particle displacement at each step
iav      = 0                    ! how often < SaSb > evaluated (if iav=0 => iav=1.5*number of particles)
regp     = 1e0                  ! regularization parameter - between 0 and 1
dpotm    = 20e0                 ! maximum change of the potential at this iteration
rtm      = 10e0                 ! keep this 
af       = 3e0                  ! erfc(af)=0. Ewald parameter, keep this. (erfc(AF) must be small) Put zero if no electrostatics forces out cut-off
fq       = 9e0                  ! exp(-fq**2)=0. Defines k-cut-off in reciprocal Ewald (exp(-FQ) must be small)
eps      = 78.3e0               ! dielectric permittivity (default water) 
temp     = 300e0                ! temperature Note: Ewald parameters, Eps and temperature define only out cut-off corrections of electrostatic interactions
iseed    = 0                    ! random seed integer number. Omit or use number lower-equal zero to use cpu clock as seeder
iprint   = 5                    ! level of output ( 5 is a good value)
nmks     = 6100000              ! num of MC-steps
nmks0    = 100000               ! num of MC steps for equilibration
lpot     = .false.              ! .t. - trial potential from the input file   .f. - otherwise mean force potential used
filrdf   = 'imc-macro.rdf'      ! file with reference RDF
filpot   = 'imc-macro-in.pot'   ! input file with "trial" potential (not needed if LPOT=.f.)
fout     = 'imc-macro-out.pot'  ! improved potential (output), to be used as input at next iteration
b1x      = 0e0                  ! lattice vector 1 x
b1y      = 0e0                  ! lattice vector 1 y 
b1z      = 0e0                  ! lattice vector 1 z
b2x      = 0e0                  ! lattice vector 2 x
b2y      = 0e0                  ! lattice vector 2 y
b2z      = 0e0                  ! lattice vector 2 x
b3x      = 0e0                  ! lattice vector 3 z
b3y      = 0e0                  ! lattice vector 3 y
b3z      = 0e0                  ! lattice vector 3 z
iout     = 1000                 ! frequency parameter for writing output
wxyz     = .false.              ! If .true. writes an xyz file
wxyznm   = 'imc-macro-out.xyz'  ! .xyz output filename 
wxyzfq   = 1000                 ! frequency of saving for xyz 
rxyz     = .false.              ! If .true. reads from an xyz file fixed and free particle coordinates and names
rxyznm   = 'imc-macro-in.xyz'   ! .xyz output filename
ldmppot  = .false.              ! dump first guess of potential or the readed potential and stop
lzm      = .true.               ! if .true., if S is zero move potential to limit
zeromove = 1.0                  ! if S is zero move potential up to this limit
lelec    = .true.               ! if .true. compute electrostatics and ewald
lrespot  = .false.              ! if .true. several restriction to the potential for pairs may be specified in respotnm
respotnm = 'imc-macro-fix.pot'  ! Format: particle type 1 particle type 2 free,fix,sas,scale,shift,fixce 
lseppot  = .false.              ! Dump separated potential files for each pair (can be used combined with ldmppot)
lseprdf  = .false.              ! Dump separated RDF and S files for each pair 
ldmp     = .false.              ! ignored, not used

! Time Stamp
call timestamp()
wall=timer()

! Header
call printheader()

! Read input
write(*,*)'Starting... read data...'
write(*,*)
read(*,input)
write(*,input)

latvec=b1y.ne.0e0.or.b1z.ne.0e0.or.b2x.ne.0e0.or.b2z.ne.0e0.or.b3x.ne.0e0.or.b3y.ne.0e0
boxlx=b1x
boxly=b2y
boxlz=b3z 
bv(1:3,1)=(/b1x, b1y, b1z/)
bv(1:3,2)=(/b2x, b2y, b2z/)
bv(1:3,3)=(/b3x, b3y, b3z/)

invbv=matmul(transpose(bv),bv)
invbv=invmat(invbv)

if (latvec) then
  lelec=.false.
  write(*,'(a)') 'Lattice Vectors: '
  write(*,'(a,3g26.18)') 'Vector 1: ',bv(1,1),bv(2,1),bv(3,1)
  write(*,'(a,3g26.18)') 'Vector 2: ',bv(1,2),bv(2,2),bv(3,2)
  write(*,'(a,3g26.18)') 'Vector 3: ',bv(1,3),bv(2,3),bv(3,3)
  write(*,'(a)') 'If lelec on, switching off'
  call volume(bv(1:3,1),bv(1:3,2),bv(1:3,3),vol)
  box=.false.
else
  write(*,'(a,3g26.18)') 'BOX Size (x y z): ',boxlx,boxly,boxlz
  iboxlx=1.0/boxlx
  iboxly=1.0/boxly
  iboxlz=1.0/boxlz
  vol=boxlx*boxly*boxlz
  box=.true.
endif
if (vol.le.0.0) stop 'Error: Box size/lattice vectors wrong'

! Get Header from RDF file
open(unit=4,file=filrdf,status='old')
read(4,*) ntyp,na,rmin,rmax  ! Reads number of particle types, number of data points, min and max distance
allocate (nms(ntyp),ch(ntyp),nspec(ntyp),nspecf(ntyp),nspecfr(ntyp),potres(ntyp,ntyp))
allocate (shelv(na),pot(na,ntyp,ntyp),ras(na),ind(na,ntyp,ntyp),ipot(na,ntyp,ntyp))
rdfinc=(rmax-rmin)/float(na)
iri=float(na)/(rmax-rmin)
aone=1-rmax/na/10e0
nms=''
do i=1,ntyp
  read(4,*) nms(i),ch(i),nspec(i)  ! Reads particle type label, particle charge, number of particles of this type
enddo

! calc nop
nop=0
do ityp=1,ntyp
  nop=nop+nspec(ityp)
enddo
inop=1e0/nop

!allocate nop dependent arrays
allocate (itype(nop)) ! INT
allocate (x(nop),y(nop),z(nop),q(nop),fpn(nop)) ! COORD
allocate (xl(nop),yl(nop),zl(nop)) ! COORD

npot=na*ntyp*(ntyp+1)/2
write(*,*) 'NPOT = ',npot

! allocate npot dependent arrays
allocate (rdf(npot)) ! RDFs
allocate (ina(npot),ityp1(npot),ityp2(npot)) ! INT
allocate (cors(npot),corp(npot,npot)) ! CORR
allocate (corsl(npot),corpl(npot,npot)) ! CORR
allocate (cross(npot,npot),diff(npot),cor(npot),iucmp(npot),ipvt(npot))

nfix=0
nspecf=0
nspecfr=0
nfree=0
cent=0e0
nfxfr=0
if (rxyz) then ! if fixed particle system activated
! Reads xyz
  open(unit=77,file=rxyznm,status='old')
  read(77,*) nfxfr
  read(77,*) nfix
  if (nfxfr.lt.nfix) stop 'Number of fixed particles is bigger than number of total particles in xyz'
  nfree=nfxfr-nfix
  if (nfxfr.gt.nop) then
    close(77)
    write(*,*) 'nfree+nfix',nfxfr
    write(*,*) 'Number of Particles in rdf file: ',nop
    stop ' nfree + nfix is bigger than number of particles declared in rdf file'
  endif
  write(*,*) 
  do i=1,nfxfr
    read(77,*) fpn(i),x(i),y(i),z(i)
  enddo
  close(77)
  ! Compute Center
  if (nfix.gt.0) then
    do i=1,nfix
      cent(1)=cent(1)+x(i)
      cent(2)=cent(2)+y(i)
      cent(3)=cent(3)+z(i)
    enddo
    cent(1)=cent(1)/nfix
    cent(2)=cent(2)/nfix
    cent(3)=cent(3)/nfix
  endif
  ! Translate system to the center
  do i=1,nfxfr
    x(i)=x(i)-cent(1)
    y(i)=y(i)-cent(2)
    z(i)=z(i)-cent(3)
  enddo
  ! Compare size of box and size of fixed particle system
  if (nfix.gt.1) call checkboxsize()
  ! Put everything in the box
  do i=1+nfix,nfxfr
    call pbc(x(i),y(i),z(i))
  enddo 
  ! compute number of fixed and free species
  do i=1,nfxfr
    j=0
    loopon=.true.
    do while (j.lt.ntyp.and.loopon)
      j=j+1
      if (nms(j).eq.fpn(i)) then 
        itype(i)=j
        q(i)=ch(j)
        if (i.le.nfix) then
          nspecf(j)=nspecf(j)+1
        else
          nspecfr(j)=nspecfr(j)+1
        endif
        loopon=.false.
      endif
    enddo
  enddo
  ! Check for species inconsistencies
  do i=1,ntyp
    if (nspecf(i).gt.nspec(i)) stop 'There are more fixed than total particles of one kind'
    if (nspecf(i)+nspecfr(i).gt.nspec(i)) stop 'There are more fixed+free than total particles of one kind'
  enddo
endif

! iav 
if (iav.le.0) then
  iav=int(1.5*(nop-nfix))
  write(*,*) 'Setting iav to ',iav
  write(*,*)
endif

write(*,*) 'Number of Particles: ',nop
write(*,*) 'Number of Fixed Part.: ',nfix
write(*,*) 'Number of Free Located Part.: ',nfree
write(*,*) 'Number of Free Unlocated Part.: ',nop-nfxfr

!  addresses remaining particles if any
i=nfxfr
do ityp=1,ntyp
  do j=1,nspec(ityp)-nspecf(ityp)-nspecfr(ityp),1
    i=i+1
    itype(i)=ityp
    fpn(i)=nms(ityp)
    q(i)=ch(ityp)
  enddo
enddo

! locate remaining free particles
! be aware that it is not avoiding overlaps
do i=nfxfr+1,nop,1
  pos=bv(:,1)*(rand()-0.5e0)+bv(:,2)*(rand()-0.5e0)+bv(:,3)*(rand()-0.5e0)
  x(i)=pos(1)
  y(i)=pos(2)
  z(i)=pos(3)
enddo

! put all free particles inside box
do i=nfix+1,nop
  call pbc(x(i),y(i),z(i))
enddo

! seed randomizer using given seed or using clock if seed lower-equal zero
call init_rand_seed(iseed) 

!  calculation of electrostatic constant 
if (lelec) then 
  if(dr.gt.min(boxlx,boxly,boxlz)/2) dr=min(boxlx,boxly,boxlz)/2
  coulf = 1e10*elch**2/(4.0*pi*eps0*eps*boltz*temp)
  write(*,*)
  write(*,*)'CoulF = ',coulf
endif

j=0
ind=.false.
do i=1,na
  do ityp=1,ntyp
    do jtyp=1,ityp
      j=j+1
      ipot(i,ityp,jtyp)=j
      ipot(i,jtyp,ityp)=j
      ina(j)=i
      ityp1(j)=ityp
      ityp2(j)=jtyp
    enddo
  enddo
enddo

! read fixed potential
potres=1 ! free all by default
if (lrespot) then
  write(*,*)
  write(*,*) 'Restricting potential for the following pairs:'
  open(unit=79,file=respotnm,status='old')
  read(79,*,iostat=kode) it1,it2,keyword
  do while (kode.eq.0)
    if (keyword.eq.'fix') then 
      restyp=0
    elseif (keyword.eq.'scale') then 
      restyp=2
    elseif (keyword.eq.'shift') then 
      restyp=3
    elseif (keyword.eq.'sas') then 
      restyp=4
    elseif (keyword.eq.'fixce') then
      restyp=5
    else ! free
      restyp=1
      keyword='free'
    endif
    potres(it1,it2)=restyp
    potres(it2,it1)=restyp
    write(*,*) it1,it2,'(',nms(it1),' ',nms(it2),')     ',keyword
    read(79,*,iostat=kode) it1,it2,keyword
  enddo 
  write(*,*)
  close(79)
endif

! input of potential if present
if(lpot)then
  open(unit=3,file=filpot,status='old')
  read(3,*)ntpp,nap,rmin,rmax
  if(nap.ne.na)stop 'Wrong na in pot.file'
  if(ntpp.ne.ntyp)stop 'Wrong ntyp in pot.file'
  do it1=1,ntyp
    do it2=it1,ntyp
      do nr=1,na 
        read(3,*)ras(nr),pot(nr,it1,it2)
        pot(nr,it2,it1)=pot(nr,it1,it2)
        if(iprint.ge.7) write(*,*) ras(nr),pot(nr,it1,it2),nr,it1,it2
      enddo
    enddo
  enddo
  close(3)   
  write(*,*)'Potential is taken from ',filpot
endif

rcut=rmax
rcut2=rcut**2
alpha=af/rcut
write(*,*)'Cut off ',rcut
write(*,*)'Alpha ',alpha
write(*,*)


! Get RDF and make first approximation
read(4,*,iostat=kode)rr,rdfp,it1,it2
do while(kode.eq.0)
  nr=(rr-rmin)*iri+aone
  if (nr.ge.1.and.nr.le.na) then 
    ip=ipot(nr,it1,it2)
    rdf(ip)=rdfp
    shift=0e0
    if (lelec) shift=ch(it1)*ch(it2)*coulf/rmax
    dlr=-log(rdfp)                
    if(rdfp.gt.0e0)then
      if(.not.lpot)pot(nr,it1,it2)=shift+dlr
      ind(nr,it1,it2)=.true.
      if(.not.lpot)pot(nr,it2,it1)=shift+dlr
      ind(nr,it2,it1)=.true.
    endif
  endif
  read(4,*,iostat=kode) rr,rdfp,it1,it2
enddo
close(4)
write(*,*)'RDFs were loaded from ',filrdf 

if(.not.lpot)write(*,*)'Mean force potential is used'

do nr=1,na
  ras(nr)=rmin+nr*rdfinc
!  shelv(nr)=4e0*pi*rdfinc*(ras(nr)**2+rdfinc**2/12e0) 
  shelv(nr)=4e0*pi*rdfinc*(ras(nr)**2+rdfinc**2*i3-rdfinc*ras(nr))
enddo

if (.not.lpot) call completepotential(ntyp)

! output potential and stop if required
if (ldmppot) then 
  !   Output of potential
  open(unit=2,file=fout,status='unknown')
  write(2,*)ntyp,na,rmin,rmax
  do it=1,ntyp
    do jt=it,ntyp
      do nr=1,na
        write(2,*)ras(nr),pot(nr,it,jt),it,jt
      enddo
    enddo
  enddo
  close(2)
! dump separated potential files for each pair
  if(lseppot) then 
    do it=1,ntyp
      do jt=it,ntyp
        call lcase(nms(it),nmst1)
        call lcase(nms(jt),nmst2)
        open(unit=13,file=trim(nmst1)//'-'//trim(nmst2)//'.pot',status='replace',form='formatted')
        do nr=1,na
          if(ind(nr,it,jt)) then 
            write(13,'(f8.3$)') ras(nr)
            write(13,*) pot(nr,it,jt)
          endif
        enddo
        close(13)
      enddo
    enddo
  endif
  write(*,*) 'Potential Files Dumped'
  stop
endif

if(iprint.ge.6) then 
  write(*,*) 'n1  n2      r       rdf          pot       ind'
  do it=1,ntyp
    do jt=1,it
      do nr=1,na
        write(*,'(2i5,3f12.4,i4)') it,jt,ras(nr),rdf(ipot(nr,it,jt)),pot(nr,it,jt),ind(nr,it,jt)
      enddo
    enddo
  enddo  
endif

if (lelec) then 
  ! Allocate more
  call calcnkvm()
  allocate (kx(nkvm),ky(nkvm),kz(nkvm),rkv(nkvm),ssin(nkvm),scos(nkvm),ddsin(nkvm),ddcos(nkvm)) ! FUR
  allocate (ssinl(nkvm),scosl(nkvm),ddsinl(nkvm),ddcosl(nkvm)) ! FUR
else
  allocate (ssinl(1),scosl(1),ddsinl(1),ddcosl(1)) ! FUR
endif

efur=0e0 
! Set initial variables if no restart
!   summators=0
nmksf=0
iud=0
nav=0 
eners=0e0 
avexp=0e0 
ave2=0e0
virs=0e0
vire=0e0
cors=0
corp=0
call fixedenergy
call energy
write(*,*) 'Fixed Particles Energy = ',fce
write(*,*) 'Fixed Particles Coulombic Energy = ',fcee
write(*,*) 'Init. energy = ',ener

if (wxyz) then 
  open(unit=88,file=wxyznm,status='unknown') 
  if (mod(istep,wxyzfq).eq.0) then
    write(88,*) nop
    write(line,*) nfix,'  Step #',nmksf,' LatVec: ',bv
    write(88,'(A)') trim(line)
    do i=1,nop
      write(88,*) fpn(i),x(i),y(i),z(i)
    enddo
  endif
endif

!$omp parallel private(tid) 
tid = omp_get_thread_num()
if (tid .eq. 0) then
  nth = omp_get_num_threads()
  print *, 'Number of threads =', nth
end if
!$omp end parallel

!   Monte Carlo  
do istep=1,nmks0
  i=(nop-nfix)*rand()+1+nfix
  it=itype(i)
  chi=ch(it)
  x1=x(i)+dr*(rand()-0.5e0)
  y1=y(i)+dr*(rand()-0.5e0)
  z1=z(i)+dr*(rand()-0.5e0)
  call pbc(x1,y1,z1)
  !  Energy difference
  de=0d0
!$omp parallel private(j,jt,dx,dy,dz,rrn2,rrn,deloc,dee,rro2,rro,nr)
  deloc=0d0
!$omp do
  do j=1,nop
    if(i.ne.j)then
      jt=itype(j)
      dx=x(j)-x1
      dy=y(j)-y1
      dz=z(j)-z1
      call pbc(dx,dy,dz)
      rrn2=dx**2+dy**2+dz**2
      if(rrn2.lt.rcut2)then
        rrn=sqrt(rrn2)
        nr=(rrn-rmin)*iri+1
        deloc=deloc+pot(nr,it,jt)
        if (lelec) then
          dee=chi*ch(jt)*coulf*(erfc(alpha*rrn)-1e0)/rrn
          deloc=deloc+dee
        endif
      endif
      dx=x(j)-x(i)
      dy=y(j)-y(i)
      dz=z(j)-z(i)
      call pbc(dx,dy,dz)
      rro2=dx**2+dy**2+dz**2
      if(rro2.lt.rcut2)then
        rro=sqrt(rro2)
        nr=(rro-rmin)*iri+1
        deloc=deloc-pot(nr,it,jt)
        if (lelec) then
          dee=chi*ch(jt)*coulf*(erfc(alpha*rro)-1e0)/rro
          deloc = deloc-dee
        endif
      endif
    endif
  enddo
!$omp end do
!$omp critical
de=de+deloc
!$omp end critical
!$omp end parallel
  !  Corrections due to electrostatics
  if (lelec) then 
    call difew(x1,y1,z1,i,chi,dee)
    de=de+dee
  endif
  !  MC transition
  if (de.le.0.0.or.(de.le.18e0.and.exp(-de).ge.rand())) then 
  !  New configuration
    x(i)=x1
    y(i)=y1
    z(i)=z1
    ener=ener+de
    iud=iud+1
    if (lelec) call recsin
  endif 
  ! writes step, energy & pressure
  if(mod(istep,iout).eq.0)then
    pres=nop-vir-efur*i3
    write(*,'(i10,a5,3f12.4,a7,f12.4)') istep,'  en:',ener*inop,vir*inop,efur*inop*i3,'  pres=',pres 
  endif
enddo   
cova=(nmks-nmks0)/(nth*iav)+1
nmksf=cova*nth*iav+nmks0
write(*,*) nmks,nmks0,nth,iav,nmksf,cova
!$omp parallel private(jobs,x1,y1,z1,i,it,chi,del,aun,istepl,j,jt,dx,dy,dz,rrn2,rrn,deloc,dee,rro2,rro,nr,xl,yl,zl,tid,isd,iudl,ssinl,scosl,ddsinl,ddcosl,efurl,virl,enerl,enersl,corsl,corpl,navl,virsl,virel,pres)
xl=x
yl=y
zl=z
iudl=0
navl=0
if (lelec) then 
  ssinl=ssin
  scosl=scos
  ddsinl=ddsin
  ddcosl=ddcos
endif
efurl=efur
virl=0.0
virel=0.0
virsl=0.0
enerl=ener
enersl=0.0
corsl=0
corpl=0
tid = omp_get_thread_num()
call system_clock(count=isd)
call srand(isd*(tid+1))
!$omp do
do jobs=1,nth
  do aun=1,cova
    do istepl=1,iav
      i=(nop-nfix)*rand()+1+nfix
      it=itype(i)
      chi=ch(it)
      x1=xl(i)+dr*(rand()-0.5e0)
      y1=yl(i)+dr*(rand()-0.5e0)
      z1=zl(i)+dr*(rand()-0.5e0)
      call pbc(x1,y1,z1)
      !  Energy difference
      del=0d0
      do j=1,nop
        if(i.ne.j)then
          jt=itype(j)
          dx=xl(j)-x1
          dy=yl(j)-y1
          dz=zl(j)-z1
          call pbc(dx,dy,dz)
          rrn2=dx**2+dy**2+dz**2
          if(rrn2.lt.rcut2)then
            rrn=sqrt(rrn2)
            nr=(rrn-rmin)*iri+1
            del=del+pot(nr,it,jt)
            if (lelec) then
              dee=chi*ch(jt)*coulf*(erfc(alpha*rrn)-1e0)/rrn
              del = del+dee
            endif
          endif
          dx=xl(j)-xl(i)
          dy=yl(j)-yl(i)
          dz=zl(j)-zl(i)
          call pbc(dx,dy,dz)
          rro2=dx**2+dy**2+dz**2
          if(rro2.lt.rcut2)then
            rro=sqrt(rro2)
            nr=(rro-rmin)*iri+1
            del=del-pot(nr,it,jt)
            if (lelec) then
              dee=chi*ch(jt)*coulf*(erfc(alpha*rro)-1e0)/rro
              del = del-dee
            endif
          endif
        endif
      enddo
      !  Corrections due to electrostatics
      if (lelec) then
        call difewpar(x1,y1,z1,i,chi,dee,ssinl,scosl,ddsinl,ddcosl)
        del=del+dee
      endif
      !  MC transition
      if (del.le.0.0.or.(del.le.18e0.and.exp(-del).ge.rand())) then
      !  New configuration
        xl(i)=x1
        yl(i)=y1
        zl(i)=z1
        enerl=enerl+del
        iudl=iudl+1
        if (lelec) call recsinpar(ssinl,scosl,ddsinl,ddcosl)
      endif
    enddo
    !   Averaging
    call correlpar(ssinl,scosl,efurl,del,virl,enerl,enersl,corsl,corpl,navl,xl,yl,zl)
    virsl=virsl+virl
    virel=virel+efurl*i3
    ! writes step, energy & pressure
    pres=nop-virl-efurl*i3
    write(*,'(i0,x,i0,x,a5,3f12.4,a7,f12.4)') tid,aun,'  en:',enerl*inop,virl*inop,efurl*inop*i3,'  pres=',pres
  enddo
enddo
!$omp end do
!$omp critical
nav=nav+navl
iud=iud+iudl
eners=eners+enersl
cors=cors+corsl
corp=corp+corpl
virs=virs+virsl
vire=vire+virel
! writes xyzq
if (wxyz) then
!  if (mod(istep,wxyzfq).eq.0) then
  write(88,*) nop
  write(line,*) nfix,'  Step #',istep,' LatVec: ',bv,' Proc# ',tid
  write(88,'(A)') trim(line)
  do i=1,nop
    x1=xl(i)
    y1=yl(i)
    z1=zl(i)
    if (i.gt.nfix) call pbc(x1,y1,z1)
    write(88,*) fpn(i),x1,y1,z1
  enddo
!  endif
endif
!$omp end critical
!$omp end parallel
if (wxyz) close(88) 

fnr=1e0/float(nav)

!time stamp and wall time
call timestamp()
write(*,*) 'Wall Time (ms): ',timer()-wall


!   Output of MC
write(*,*)'Calculation of effective potentials'
write(*,*)'Program ',label,' from ',datestamp
write(*,*)'-----------------------------------'
write(*,*)'File with original RDF:',filrdf
write(*,*)'Number of particles:   ',nop
write(*,*)'Number of species:     ',ntyp  
write(*,*)'Sampled configurations:',nav
write(*,*)'Number of particles of each species:'
do i=1,ntyp
  write(*,'(i4,i6,3x,a6,f6.2)')i,nspec(i),nms(i),ch(i) 
enddo 
write(*,*)
write(*,*)'Total MC steps:        ',nmksf
write(*,*)'MC steps with averaging',nmksf-nmks0
write(*,*)'Acceptance ratio       ',float(iud)/float(nmksf)
write(*,*)'Max. displacemet       ',dr
write(*,*)'Temperature            ',temp
write(*,*)'Dielectric permitivity ',eps
write(*,*)'------------------------------------'
write(*,*)'Average energy:        ',eners*fnr*inop 
pres=nop-(virs+vire)*fnr 
osm =pres*inop
write(*,*)'Presure                ',pres
write(*,*)'Osmotic coefficient    ',osm
write(*,*)'Pressure components    ',nop,-virs*fnr,-vire*fnr 
write(*,*)'------------------------------------'

! Radial distribution function
if(iprint.ge.6)write(*,*)'Radial distribution functions:'
felc=0e0
felr=0e0
do it=1,ntyp
  do jt=it,ntyp
    if(iprint.ge.6)then
      write(*,*)' This pair:',nms(it),' - ',nms(jt) 
      write(*,*)'   R       RDF   RDFref     <S(R)> ','   <S(R)ref>    delta           typ '
    endif
    do nr=1,na 
      ipt=ipot(nr,it,jt)
      crc=float(cors(ipt))*fnr 
      rdfref=rdf(ipt)
      crr=rdfref*shelv(nr)*nspec(it)*nspec(jt)/vol
      rdfc=crc*vol/(nspec(it)*nspec(jt)*shelv(nr))
      if(it.eq.jt)then
!   may be important!
!   How do you normalize RDF between likewise particles?
!   The same should be done ~ 100 lines below
!   Normalize on N*(N-1)/2
!	        fac=0.5*(dfloat(nspec(it))-1.d0)/dfloat(nspec(it))
!   Normalize on N**2/2
        crr=crr*(nspec(it)-1)/(nspec(it)*2e0)
        rdfc=rdfc*2e0*dfloat(nspec(it))/dfloat(nspec(it)-1)
      endif 
      felc=felc+(crc-crr)**2                         ! total error
      felr=felr+(rdfref-rdfc)**2
      difrc=(crc-crr)*regp 
      if(crc.ne.0e0)then
        if(difrc/crc.gt. rtm)difrc= crc*rtm
        if(difrc/crc.lt.-rtm)difrc=-crc*rtm 
      endif
      cor(ipt)=difrc
      if(ind(nr,it,jt).and.iprint.ge.6) write(*,'(f9.4,2f10.5,3f10.4,5x,3a4)') ras(nr),rdfc,rdfref,crc,crr,difrc,'rdf:',nms(it),nms(jt)
    enddo
  enddo
enddo
! check here 
felc=sqrt(felc/npot)
felr=sqrt(felr/npot)
write(*,'(2(A,F16.10))') 'Standard Deviation from reference: RDF => ',felr,'  S => ',felc

!   Analysis
!   Zeros elimination
ic=0
do i=1,npot
  nr=ina(i)
  it1=ityp1(i)
  it2=ityp2(i)
  if(ind(nr,it1,it2).and.cors(i).ne.0) then
    if (potres(it1,it2).ge.1)then ! free, shift, scale or sas
      ic=ic+1
      iucmp(ic)=i
      diff(ic)=cor(i)
!   elseif (potres(it1,it2).eq.0) then ! fix
!     scale, shift and sas continue below
    endif 
  endif
enddo

nur=ic 

write(*,*)' Number of equations ',nur 
write(*,*)' regularization parameter ',regp
write(*,*)' max. relative change rdf ',rtm
write(*,*)' max. change of potential ',dpotm 

do ic=1,nur
  i=iucmp(ic)
  do jc=ic,nur 
    j=iucmp(jc)
    if(i.gt.j)write(*,*)'!!!  Violation in uncompressing' 
    cross(ic,jc) = float(corp(i,j))*fnr-float(cors(i))*float(cors(j))*fnr**2
    cross(jc,ic)=cross(ic,jc)
  enddo
enddo 

if(iprint.ge.8)then
  do ic=1,nur    
    i=iucmp(ic)
    write(*,'(i4,10(200f10.6/4x))') i,(cross(ic,jc),jc=1,nur)
  enddo
endif
    
! CROSS * X = DIFF
! nur -> number of linear equations. order of cross matrix
! 1 -> number of columns of matrix DIFF
! npot -> leading dimension of cross
! CROSS -> input CROSS(npot,nur) / output factors L and U of factorization of CROSS=P*L*U
! ipvt -> integer array with nur dimension corresponds to permutation matrix P
! DIFF -> input DIFF(NUR,1) / output X(NUR,1)
! INFO if 0 solution obtained successfully
call sgesv(nur,1,cross,npot,ipvt,diff,npot,info)

if(info.ne.0)then
   write(*,*)'STOP singular correlation matrix. Error #',info
   stop
endif

cor=0e0
do ic=1,nur
  i=iucmp(ic)
  cor(i)=diff(ic)
enddo
do i=1,npot
  nr=ina(i)
  it1=ityp1(i)
  it2=ityp2(i)
  if(lzm.and.ind(nr,it1,it2).and.cors(i).eq.0)cor(i)=-zeromove
  if(cor(i).lt.-dpotm)cor(i)=-dpotm
  if(cor(i).gt.dpotm)cor(i)=dpotm
  if(potres(it1,it2).eq.0)cor(i)=0e0
enddo

if (lrespot) then 
  do it=1,ntyp
    do jt=it,ntyp
      vrvs=0e0
      vr=0e0
      vs=0e0
      vs2=0e0
      it1=0
      do nr=1,na
        ipt=ipot(nr,it,jt)
        if (ind(nr,it,jt).and.cors(ipt).ne.0) then
          it1=it1+1
          ! SUM(VR*VS)
          vrvs=vrvs+(pot(nr,it,jt)+cor(ipt))*pot(nr,it,jt) ! SUM(VR*VS)
          vr=vr+pot(nr,it,jt)+cor(ipt) ! SUM(VR) (reference vector)
          vs=vs+pot(nr,it,jt) ! SUM(VS) (original vector to modify)
          vs2=vs2+pot(nr,it,jt)**2! SUM(VS^2)
        endif
      enddo
      if (potres(it,it).eq.2) then ! scale
        a=vrvs/vs2
      elseif (potres(it,it).eq.3.or.potres(it,it).eq.5) then ! shift or fixcenter
        b=(vr-vs)/it1
      elseif (potres(it,it).eq.4) then ! scale and shift
        a=(it1*vrvs-vr*vs)/(it1*vs2-vs**2)
        b=(vr-a*vs)/it1
      endif
      do nr=1,na
        if (ind(nr,it,jt)) then
          ipt=ipot(nr,it,jt)
          if (potres(it,it).eq.2) then ! scale
            cor(ipt)=(a-1)*pot(nr,it,jt)
          elseif (potres(it,it).eq.3) then ! shift
            cor(ipt)=b
          elseif (potres(it,it).eq.4) then ! scale and shift
            cor(ipt)=(a-1)*pot(nr,it,jt)+b
          elseif (potres(it,it).eq.5) then ! fix center
            cor(ipt)=cor(ipt)-b
          endif
        endif
      enddo
    enddo
  enddo
endif

!   Output of potential
open(unit=2,file=fout,status='unknown')
write(2,*)ntyp,na,rmin,rmax
write(*,*)'------------------------------------'
write(*,*)'Corrected effective potential'
do it=1,ntyp
  do jt=it,ntyp
    if(nspec(it).ne.nspecf(it).or.nspec(jt).ne.nspecf(jt)) then
      write(*,*)' This pair:',nms(it),' - ',nms(jt) 
      write(*,*)'   R       RDF    RDFref   Potential     init.Pot.'
    endif
    if (lseppot.or.lseprdf) then
      call lcase(nms(it),nmst1)
      call lcase(nms(jt),nmst2)
      if (lseppot) then 
        open(unit=13,file=trim(nmst1)//'-'//trim(nmst2)//'.pot',status='replace',form='formatted')
      endif
      if (lseprdf) then
        open(unit=14,file=trim(nmst1)//'-'//trim(nmst2)//'.rdf',status='replace',form='formatted')
      endif
    endif
    do nr=1,na 
      ipt=ipot(nr,it,jt)
      crc=float(cors(ipt))*fnr
      rdfc=crc*vol/(nspec(it)*nspec(jt)*shelv(nr))
      if(it.eq.jt)then
!   may be important!
!   How do you normalize RDF between likewise particles?
!   The same should be done ~ 100 lines below
!   Normalize on N*(N-1)/2
!	        fac=0.5*(dfloat(nspec(it))-1.d0)/dfloat(nspec(it))
!   Normalize on N**2/2
        rdfc=rdfc*2e0*dfloat(nspec(it))/dfloat(nspec(it)-1)
      endif
      poten=pot(nr,it,jt)
      potnew=poten+cor(ipt) 
      if(ind(nr,it,jt)) write(*,'(f9.4,5f10.5,7x,3a4) ') ras(nr),rdfc,rdf(ipt),potnew,poten,cor(ipt),'pot:',nms(it),nms(jt)
      write(2,*)ras(nr),potnew,it,jt
      if(lseppot.and.ind(nr,it,jt)) then
        write(13,'(f8.3$)') ras(nr)
        write(13,*) potnew
      endif
      if(lseprdf.and.ind(nr,it,jt)) then
        write(14,'(f8.3$)') ras(nr)
        write(14,*) rdfc,cors(ipt)
      endif
    enddo
    if(lseppot) close(13)
    if(lseprdf) close(14)
  enddo
enddo
close(2)
write(*,*) 'Normal termination'
! Time Stamp
call timestamp()
write(*,*) 'Wall Time (ms): ',timer()-wall
end program

subroutine calcnkvm()
use invmc
implicit none
integer ikx,iky,ikz,ikv,ky1,ky2,kmaxx,kmaxy,kmaxz
real    cutk2,rk2,pial

if(alpha.gt.0.001)then
  pial=(pi/alpha)**2
  cutk2=fq/pial
  kmaxx=sqrt(cutk2)*boxlx+1
  kmaxy=sqrt(cutk2)*boxly+1
  kmaxz=sqrt(cutk2)*boxlz+1
else
  kmaxx=0
  kmaxy=0
  kmaxz=0
endif

ikv=0
do ikz=-kmaxz,kmaxz
  do ikx=0,kmaxx
    if(ikx.eq.0)then
      if(ikz.gt.0)then
        ky1=0
      else
        ky1=1
      endif
      ky2=kmaxy
    else
      ky1=-kmaxy
      ky2=kmaxy-ikx
    endif
    do iky=ky1,ky2
      rk2=(ikx*iboxlx)**2+(iky*iboxly)**2+(ikz*iboxlz)**2
      if(rk2.le.cutk2) ikv=ikv+1
    enddo
  enddo
enddo
nkvm=ikv
end subroutine
!   
!=============== ewald =======================================
!
subroutine ewald(qpe)
use invmc
implicit none
logical,save   :: init
integer      :: ikx,iky,ikz,ikv,ky1,ky2,kmaxx,kmaxy,kmaxz,i,kxv,kyv,kzv
real         :: pial,rk2,cutk2,cfur,aff,scs,ssn,sc,ss1,cc1
real*8       :: qpe
data init/.true./
! save 
! initialisation 
if(init)then
  init=.false.
  if(alpha.gt.0.001)then
    pial=(pi/alpha)**2
    cutk2=fq/pial
    kmaxx=sqrt(cutk2)*boxlx+1
    kmaxy=sqrt(cutk2)*boxly+1
    kmaxz=sqrt(cutk2)*boxlz+1
  else
    pial=1e0
    cutk2=0e0
    kmaxx=0
    kmaxy=0
    kmaxz=0
    write(*,*)'No Ewald summation '
  endif
  ikv=0
  cfur=coulf/(vol*pi)
  do ikz=-kmaxz,kmaxz
    do ikx=0,kmaxx
      if(ikx.eq.0)then
        if(ikz.gt.0)then
          ky1=0
        else
          ky1=1
        endif
        ky2=kmaxy
      else
        ky1=-kmaxy
        ky2=kmaxy-ikx
      endif
      do iky=ky1,ky2
        rk2=(ikx*iboxlx)**2+(iky*iboxly)**2+(ikz*iboxlz)**2
        if(rk2.le.cutk2)then
          ikv=ikv+1
!          if(ikv.gt.nkvm)stop 'Increase nkvm'
          kx(ikv)=ikx
          ky(ikv)=iky
          kz(ikv)=ikz
          rkv(ikv)=cfur*exp(-rk2*pial)/rk2
        endif
      enddo
    enddo
  enddo      
  nkv=ikv
  write(*,*)'kmaxx=',kmaxx,' kmaxy=',kmaxy,'kmaxz=',kmaxz,'    Num. of k-vectors ',nkv
  esf=0e0
  aff=alpha/sqrt(pi)
  do i=1,nop
    esf=esf-aff*coulf*q(i)**2
  enddo 
endif
qpe = esf
!   reciprocal space Evald
do ikv=1,nkv
  kxv=kx(ikv)
  kyv=ky(ikv)
  kzv=kz(ikv)
  scs=0e0
  ssn=0e0
!$DIR NO_RECURRENCE
  do i=1,nop
    sc=(kxv*x(i)*iboxlx + kyv*y(i)*iboxly + kzv*z(i)*iboxlz)*pi2 
    ss1=q(i)*sin(sc)
    cc1=q(i)*cos(sc)
    ssn=ssn+ss1
    scs=scs+cc1 
  enddo
  ssin(ikv)=ssn
  scos(ikv)=scs
  qpe=qpe+rkv(ikv)*(ssn**2+scs**2)
enddo
end subroutine
!   
!=============== ewaldpar =======================================
!
subroutine ewaldpar(qpe,ssinl,scosl)
use invmc
implicit none
logical,save :: init
integer      :: ikx,iky,ikz,ikv,ky1,ky2,kmaxx,kmaxy,kmaxz,i,kxv,kyv,kzv
real         :: pial,rk2,cutk2,cfur,aff,scs,ssn,sc,ss1,cc1,ssinl(*),scosl(*)
real*8       :: qpe
data init/.true./
! save 
! initialisation 
if(init)then
  init=.false.
  if(alpha.gt.0.001)then
    pial=(pi/alpha)**2
    cutk2=fq/pial
    kmaxx=sqrt(cutk2)*boxlx+1
    kmaxy=sqrt(cutk2)*boxly+1
    kmaxz=sqrt(cutk2)*boxlz+1
  else
    pial=1e0
    cutk2=0e0
    kmaxx=0
    kmaxy=0
    kmaxz=0
    write(*,*)'No Ewald summation '
  endif
  ikv=0
  cfur=coulf/(vol*pi)
  do ikz=-kmaxz,kmaxz
    do ikx=0,kmaxx
      if(ikx.eq.0)then
        if(ikz.gt.0)then
          ky1=0
        else
          ky1=1
        endif
        ky2=kmaxy
      else
        ky1=-kmaxy
        ky2=kmaxy-ikx
      endif
      do iky=ky1,ky2
        rk2=(ikx*iboxlx)**2+(iky*iboxly)**2+(ikz*iboxlz)**2
        if(rk2.le.cutk2)then
          ikv=ikv+1
!          if(ikv.gt.nkvm)stop 'Increase nkvm'
          kx(ikv)=ikx
          ky(ikv)=iky
          kz(ikv)=ikz
          rkv(ikv)=cfur*exp(-rk2*pial)/rk2
        endif
      enddo
    enddo
  enddo
  nkv=ikv
  write(*,*)'kmaxx=',kmaxx,' kmaxy=',kmaxy,'kmaxz=',kmaxz,'    Num. of k-vectors ',nkv
  esf=0e0
  aff=alpha/sqrt(pi)
  do i=1,nop
    esf=esf-aff*coulf*q(i)**2
  enddo
endif
qpe = esf
!   reciprocal space Evald
do ikv=1,nkv
  kxv=kx(ikv)
  kyv=ky(ikv)
  kzv=kz(ikv)
  scs=0e0
  ssn=0e0
!$DIR NO_RECURRENCE
  do i=1,nop
    sc=(kxv*x(i)*iboxlx + kyv*y(i)*iboxly + kzv*z(i)*iboxlz)*pi2
    ss1=q(i)*sin(sc)
    cc1=q(i)*cos(sc)
    ssn=ssn+ss1
    scs=scs+cc1
  enddo
  ssinl(ikv)=ssn
  scosl(ikv)=scs
  qpe=qpe+rkv(ikv)*(ssn**2+scs**2)
enddo
end subroutine
!
!========================= difew ================================
!       
subroutine difew(x1,y1,z1,i,chi,dee)
use invmc
implicit none
integer ikv,ikx,iky,ikz,i
real dee,chi,x1,y1,z1,ddd,dcs,scn,sco,scp,sinr,dsn
!   reciprocal space Evald
dee=0e0
chi=2e0*chi 
do ikv=1,nkv
  ikx=kx(ikv)
  iky=ky(ikv)
  ikz=kz(ikv)
  sco=(ikx*x(i)*iboxlx + iky*y(i)*iboxly + ikz*z(i)*iboxlz)*pi2
  scn=(ikx*x1  *iboxlx + iky*y1  *iboxly + ikz*z1  *iboxlz)*pi2 
  scp=sco+scn
  sinr=sin(scn-sco)*chi
  dsn=cos(scp)*sinr
  dcs=-sin(scp)*sinr
  ddsin(ikv)=dsn
  ddcos(ikv)=dcs
  ddd=(2e0*ssin(ikv)+dsn)*dsn+(2e0*scos(ikv)+dcs)*dcs  
  dee=dee+ddd*rkv(ikv)
enddo
end subroutine

!
!========================= difew-par ================================
!       
subroutine difewpar(x1,y1,z1,i,chi,dee,ssinl,scosl,ddsinl,ddcosl)
use invmc
implicit none
integer ikv,ikx,iky,ikz,i
real dee,chi,x1,y1,z1,ddd,dcs,scn,sco,scp,sinr,dsn
real ssinl(*),scosl(*),ddsinl(*),ddcosl(*)

!   reciprocal space Evald
dee=0e0
chi=2e0*chi
do ikv=1,nkv
  ikx=kx(ikv)
  iky=ky(ikv)
  ikz=kz(ikv)
  sco=(ikx*x(i)*iboxlx + iky*y(i)*iboxly + ikz*z(i)*iboxlz)*pi2
  scn=(ikx*x1  *iboxlx + iky*y1  *iboxly + ikz*z1  *iboxlz)*pi2
  scp=sco+scn
  sinr=sin(scn-sco)*chi
  dsn=cos(scp)*sinr
  dcs=-sin(scp)*sinr
  ddsinl(ikv)=dsn
  ddcosl(ikv)=dcs
  ddd=(2e0*ssinl(ikv)+dsn)*dsn+(2e0*scosl(ikv)+dcs)*dcs
  dee=dee+ddd*rkv(ikv)
enddo
end subroutine
 
!
!=================== recsin ========================================
!
subroutine recsin
use invmc
implicit none
integer ikv

do ikv=1,nkv
  ssin(ikv)=ssin(ikv)+ddsin(ikv)
  scos(ikv)=scos(ikv)+ddcos(ikv)
enddo      
end subroutine

!
!=================== recsin-par ========================================
!
subroutine recsinpar(ssinl,scosl,ddsinl,ddcosl)
use invmc
implicit none
integer ikv
real ssinl(*),scosl(*),ddsinl(*),ddcosl(*)

do ikv=1,nkv
  ssinl(ikv)=ssinl(ikv)+ddsinl(ikv)
  scosl(ikv)=scosl(ikv)+ddcosl(ikv)
enddo
end subroutine

!
!==================== fixedenergy =======================
!
subroutine fixedenergy
use invmc
implicit none
integer i,j,it,jt,nr
real dde,dee,dx,dy,dz,rr,rr2

! compute fixed coulombic energy
fce=0d0
fcee=0d0
do i=2,nfix,1
  it=itype(i)
  do j=1,i-1,1
    jt=itype(j)
    dx=x(i)-x(j)
    dy=y(i)-y(j)
    dz=z(i)-z(j)
    call pbc(dx,dy,dz)
    rr2=dx**2+dy**2+dz**2
    if(rr2.lt.rcut2)then
      rr=sqrt(rr2)
      nr=(rr-rmin)*iri+1
      dde=pot(nr,it,jt)
      if (lelec) then 
        dee=ch(it)*ch(jt)*coulf*(erfc(alpha*rr)-1e0)/rr
        fce =fce + dee
        fcee=fcee+ dee
      endif
      fce =fce + dde
    endif
  enddo
enddo
end subroutine
! 
!===================== energy ========================================
!     
subroutine energy
use invmc
implicit none
integer i,j,it,jt,nr
real dee,dde,dx,dy,dz,rr
real*8 enew

if (lelec) then 
  call ewald(efur) 
  write(*,*)' Efur=',efur
  enew=efur
else
  enew=0d0
endif
de=0d0
do i=1+nfix,nop   
  it=itype(i)
  do j=1,i-1,1 
    jt=itype(j) 
    dx=x(i)-x(j)
    dy=y(i)-y(j)
    dz=z(i)-z(j)      
    call pbc(dx,dy,dz)
    rr=sqrt(dx**2+dy**2+dz**2)
    if(rr.lt.rcut)then
      nr=(rr-rmin)*iri+1
      dde=pot(nr,it,jt)
      if (lelec) then 
        dee=ch(it)*ch(jt)*coulf*(erfc(alpha*rr)-1e0)/rr
        efur=efur+dee 
        de =de + dee
      endif
      de =de + dde 
    endif
  enddo
enddo 
ener=enew+de+fce
if (lelec) efur=efur+fcee
end subroutine
! 
!===================== correl ========================================
!     
subroutine correl                                                                
use invmc
implicit none
real dde,dee,dx,dy,dz,rr,rr2 !,eit(nop),ett
real*8 enew
integer i,j,it,jt,nr,nr1,nr2,ipt,ccor(npot),ccc

if (lelec) then 
  call ewald(enew) 
  efur=enew
else
  enew=0d0
endif
de=0d0     
vir=0e0
ccor=0
!eit=0e0
do i=1+nfix,nop   
  it=itype(i) 
  do j=1,i-1,1
    jt=itype(j)
    dx=x(i)-x(j)
    dy=y(i)-y(j)
    dz=z(i)-z(j)
    call pbc(dx,dy,dz)
    rr2=dx**2+dy**2+dz**2
    if(rr2.lt.rcut2)then
      rr=sqrt(rr2)
      nr=(rr-rmin)*iri+1
      dde=pot(nr,it,jt)
      if (lelec) then 
        dee=ch(it)*ch(jt)*coulf*(erfc(alpha*rr)-1e0)/rr
        efur=efur+dee 
        de=de+dee
      endif
      de =de + dde 
      nr1=(rr-rmin)*iri+0.5e0
      if(nr1.lt.1)then
        nr1=1 
      elseif(nr1+1.gt.na) then
        nr1=na-1
      endif
      nr2=nr1+1         
      vir=vir+rr*(pot(nr2,it,jt)-pot(nr1,it,jt)) !/(nr2-nr1)        
      ipt=ipot(nr,it,jt)
      ccor(ipt)=ccor(ipt)+1
    endif
  enddo
enddo  
vir=vir*iri*i3
enew=enew+de+fce
if (lelec) efur=efur+fcee
if(abs(ener-enew).gt.0.05)write(*,'(i10,a,f12.4,a,f12.4)') istep,' En.error.  New:',enew,'    old:',ener
ener=enew
!   correlators                       
do i=1,npot
  ccc=ccor(i)
  if(ccc.ne.0)then
    cors(i)=cors(i)+ccc
    do j=i,npot
      corp(i,j)=corp(i,j)+ccc*ccor(j)
    enddo
  endif
enddo 
! cors = S
! corp = SaSb
eners=eners+ener
nav=nav+1
end subroutine
! 
!===================== correlpar ========================================
!     
subroutine correlpar(ssinl,scosl,efurl,del,virl,enerl,enersl,corsl,corpl,navl,xl,yl,zl)
use invmc
implicit none
real dde,dee,dx,dy,dz,rr,rr2,virl !,eit(nop),ett
real*8 enew,efurl,del,enerl,enersl
integer i,j,it,jt,nr,nr1,nr2,ipt,ccor(npot),ccc,navl
real ssinl(*),scosl(*),xl(*),yl(*),zl(*)
integer corsl(npot),corpl(npot,npot)

if (lelec) then
  call ewaldpar(enew,ssinl,scosl)
  efurl=enew
else
  enew=0d0
endif
del=0d0
virl=0e0
ccor=0
!eit=0e0
do i=1+nfix,nop
  it=itype(i)
  do j=1,i-1,1
    jt=itype(j)
    dx=xl(i)-xl(j)
    dy=yl(i)-yl(j)
    dz=zl(i)-zl(j)
    call pbc(dx,dy,dz)
    rr2=dx**2+dy**2+dz**2
    if(rr2.lt.rcut2)then
      rr=sqrt(rr2)
      nr=(rr-rmin)*iri+1
      dde=pot(nr,it,jt)
      if (lelec) then
        dee=ch(it)*ch(jt)*coulf*(erfc(alpha*rr)-1e0)/rr
        efurl=efurl+dee
        del=del+dee
      endif
      del =del + dde
      nr1=(rr-rmin)*iri+0.5e0
      if(nr1.lt.1)then
        nr1=1
      elseif(nr1+1.gt.na) then
        nr1=na-1
      endif
      nr2=nr1+1
      virl=virl+rr*(pot(nr2,it,jt)-pot(nr1,it,jt)) !/(nr2-nr1)        
      ipt=ipot(nr,it,jt)
      ccor(ipt)=ccor(ipt)+1
    endif
  enddo
enddo
virl=virl*iri*i3
enew=enew+del+fce
if (lelec) efurl=efurl+fcee
if(abs(enerl-enew).gt.0.05)write(*,'(a,f12.4,a,f12.4)') ' En.error.  New:',enew,'    old:',enerl
enerl=enew
!   correlators                       
do i=1,npot
  ccc=ccor(i)
  if(ccc.ne.0)then
    corsl(i)=corsl(i)+ccc
    do j=i,npot
      corpl(i,j)=corpl(i,j)+ccc*ccor(j)
    enddo
  endif
enddo
! cors = S
! corp = SaSb
enersl=enersl+enerl
navl=navl+1
end subroutine

!======================= pbc ========================================
!

subroutine pbc(xx,yy,zz)
use invmc
implicit none
real xx,yy,zz
if (box) then
  call pbcbox(xx,yy,zz)
else
  call pbcgener(xx,yy,zz)
endif
end subroutine

subroutine pbcgener(xx,yy,zz)
use invmc
use math
implicit none
real pos(3),xx,yy,zz,vec(3)
pos=(/xx,yy,zz/)
vec=matmul(invbv,matmul(pos,bv))
pos=pos-matmul(bv,inintv(vec))
xx=pos(1)
yy=pos(2)
zz=pos(3)
end subroutine

subroutine pbcbox(xx,yy,zz)
use invmc
use math
implicit none
real xx,yy,zz
xx = xx - boxlx * inint(xx*iboxlx)
yy = yy - boxly * inint(yy*iboxly)
zz = zz - boxlz * inint(zz*iboxlz)
end subroutine

subroutine init_rand_seed(iseed)
implicit none
integer iseed
if (iseed.le.0) call system_clock(count=iseed)
call srand(iseed)
write(*,*) 'Using seed: ',iseed
end subroutine

!*****************************************************************************80
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!  Example:
!    31 May 2001   9:45:54.872 AM
!  Modified:
!    06 August 2005
!  Author:
!    john Burkardt
!  Parameters:
!    None
subroutine timestamp ( )
  implicit none
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm
end subroutine

function timer()
implicit none
integer*8 ta(8),timer
integer*8 dpm(12),dpy
dpm(1)=0
dpm(2)=31+dpm(1)
dpm(3)=28+dpm(2)
call date_and_time(values=ta)
if (mod(ta(1),4).eq.0) dpm(3)=dpm(3)+1
dpm(4)=31+dpm(3)
dpm(5)=30+dpm(4)
dpm(6)=31+dpm(5)
dpm(7)=30+dpm(6)
dpm(8)=31+dpm(7)
dpm(9)=31+dpm(8)
dpm(10)=30+dpm(9)
dpm(11)=31+dpm(10)
dpm(12)=30+dpm(11)
dpy=int((ta(1)-1)/4)*(4*365.25)+mod(ta(1)-1,4)*365
timer=((((dpy+dpm(ta(2))+(ta(3)-1))*24+ta(5))*60+ta(6))*60+ta(7))*1000+ta(8)
end function

subroutine checkboxsize()
use invmc
implicit none
integer i
real minr(3),maxr(3),diff(3)
minr(1)=x(1)
maxr(1)=x(1)
minr(2)=y(1)
maxr(2)=y(1)
minr(3)=z(1)
maxr(3)=z(1)
do i=2,nfix
  if (x(i).lt.minr(1)) minr(1)=x(i)
  if (x(i).gt.maxr(1)) maxr(1)=x(i)
  if (y(i).lt.minr(2)) minr(2)=y(i)
  if (y(i).gt.maxr(2)) maxr(2)=y(i)
  if (z(i).lt.minr(3)) minr(3)=z(i)
  if (z(i).gt.maxr(3)) maxr(3)=z(i)
enddo
diff(1)=maxr(1)-minr(1)
diff(2)=maxr(2)-minr(2)
diff(3)=maxr(3)-minr(3)
if (diff(1).ge.boxlx) then
  write(*,*) 'Fixed particle system is too big, increase the box size in x more than ',diff(1)
  stop
elseif (diff(2).ge.boxly) then
  write(*,*) 'Fixed particle system is too big, increase the box size in y more than ',diff(2)
  stop
elseif (diff(3).ge.boxlz) then
  write(*,*) 'Fixed particle system is too big, increase the box size in z more than ',diff(3)
  stop
else 
  write(*,*) 'Fixed particle system fits into periodic box.'
endif
end subroutine

subroutine printheader()
use invmc
implicit none
write(*,*) ' ___________________________________________________ '
write(*,*) '|     IMC for Macromolecules                        |'
write(*,*) '|     (Inverse Monte Carlo for Macromolecules)      |'
write(*,*) '|     ',label,       '     ',datestamp, '         |'
write(*,*) '|___________________________________________________|'
write(*,*) '|                                                   |'
write(*,*) '|      Author:  Pablo M. De Biase                   |'
write(*,*) '|               Biological Sciences                 |'
write(*,*) '|               University of Calgary               |'
write(*,*) '|               Calgary, Alberta, Canada            |'
write(*,*) '|               e-mail: pablodebiase@gmail.com      |'
write(*,*) '|                                                   |'
write(*,*) '|      Based on IMC 2.1 from:                       |'
write(*,*) '|               Alexander Lyubartsev                |' 
write(*,*) '|               Division of Physical Chemistry      |'
write(*,*) '|               Arrhenius Lab.                      |'
write(*,*) '|               Stockholm University                |'
write(*,*) '|               S-10691 Stockholm Sweden            |'
write(*,*) '|               e-mail: sasha@physc.su.se           |'
write(*,*) '|___________________________________________________|'
write(*,*) 
end subroutine

subroutine completepotential(ntyp)
use invmc
implicit none
integer i,it,jt,ntyp,ppi,ppf,nr,ip
real kc,bc,shift

do it=1,ntyp
  do jt=1,it
    ppi=0
    ppf=0
    shift=0e0
    if (lelec) shift=ch(it)*ch(jt)*coulf/rmax
    do nr=1,na
      ip=ipot(nr,it,jt)
      if(.not.ind(nr,it,jt))then
        if (nspec(it)==nspecf(it).and.nspec(jt)==nspecf(jt)) then
          if(.not.lpot) then 
            pot(nr,it,jt)=0e0
            pot(nr,jt,it)=0e0
          endif
          rdf(ip)=1e0
        else
          ppf=nr+1
!          pot(nr,it,jt)=500e0+100e0*(na-nr)
!          pot(nr,jt,it)=500e0+100e0*(na-nr)
!          rdf(ip)=0e0
        endif
      else
        if (nspec(it).ne.nspecf(it).or.nspec(jt).ne.nspecf(jt)) then
          if (ppf.gt.ppi) then
            if (ppf.eq.nr) then
              if (ppi.eq.0) then
                if (lpot) then 
                  kc=(pot(nr,it,jt))*ras(ppf)**12
                else
                  kc=(shift-log(rdf(ip)))*ras(ppf)**12
                endif
                do i=ppi+1,ppf-1
                  pot(i,it,jt)=kc/ras(i)**12
                  rdf(ipot(i,it,jt))=exp(-pot(i,it,jt))
                  pot(i,jt,it)=pot(i,it,jt)
                enddo
              else
                kc=(rdf(ipot(ppf,it,jt))-rdf(ipot(ppi,it,jt)))/(ras(ppf)-ras(ppi))
                bc=(rdf(ipot(ppi,it,jt))*ras(ppf)-rdf(ipot(ppf,it,jt))*ras(ppi))/(ras(ppf)-ras(ppi))
                do i=ppi+1,ppf-1
                  rdf(ipot(i,it,jt))=kc*ras(i)+bc
                  if (.not.lpot) then 
                    pot(i,it,jt)=shift-log(rdf(ipot(i,it,jt)))
                    pot(i,jt,it)=pot(i,it,jt)
                  endif
                enddo
                if (lpot) then
                  kc=(pot(ppf,it,jt)-pot(ppi,it,jt))/(ras(ppf)-ras(ppi))
                  bc=(pot(ppi,it,jt)*ras(ppf)-pot(ppf,it,jt)*ras(ppi))/(ras(ppf)-ras(ppi))
                  do i=ppi+1,ppf-1
                    pot(i,it,jt)=kc*ras(i)+bc
                    pot(i,jt,it)=pot(i,it,jt)
                  enddo
                endif
              endif
            endif
          endif
          ppi=nr
        endif
      endif
    enddo
    if (ppf.eq.na+1) then
      do i=ppi+1,na
        rdf(ipot(i,it,jt))=rdf(ipot(ppi,it,jt))
        if (lpot) then
          pot(i,it,jt)=pot(ppi,it,jt) 
          pot(i,jt,it)=pot(ppi,it,jt) 
        else
          pot(i,it,jt)=shift-log(rdf(ipot(i,it,jt)))
          pot(i,jt,it)=pot(i,it,jt)
        endif
      enddo
    endif
  enddo
enddo
end subroutine

subroutine lcase(inchar,outchar)
implicit none
integer i
character s*1
character*(*) inchar,outchar

do i=1,len_trim(inchar)
  s=inchar(i:i)
  if (iachar(s).ge.65.and.iachar(s).le.90) then
    outchar(i:i)=achar(iachar(s)+32)
  else
    outchar(i:i)=inchar(i:i)
  endif
enddo
outchar=outchar(1:len_trim(inchar))
end subroutine

subroutine volume(boxv1,boxv2,boxv3,vvv)
use math
implicit none
real boxv1(3),boxv2(3),boxv3(3),vvv

vvv=dot_product(cross_product(boxv1,boxv2),boxv3)

end subroutine

