!    XSSURF - Reads BROMOC Trajectory files and Repulsion maps to computes the 
!             cross-section free surface along z-axis considering the pore and DNA
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

module btr
implicit none
integer*4,allocatable :: typ(:),maxntyp(:)
integer*4 ntop,itype,nframe,maxntop,ir,lr
real*8 runtime
character*4, allocatable :: atnam2(:)
character bs*8
real*4,allocatable :: rt(:,:)
real*8 :: xtlabc6(6)
integer*4 icntrl(20),itemp
character hdr*4
integer*4,parameter :: ntitle=160
character ( len = ntitle ) title
integer ndna,dtype
real*8,allocatable :: rad(:)
real*8,parameter :: ab=4.15,tb=4.05,cb=3.4,gb=3.4,sug=3.45,pho=3.55 !radii
! d: effective radius of the species
!dCLA is the minimum distance in efpot table from CLA to the other species
!dPOT is the minimum distance in efpot table from POT to the other species
!rPOT is the Lennard Jones sigma from MD parameters for POT (1.76375)
!rCLA is the Lennard Jones sigma from MD params for CLA (2.27)
!dCLA values
!aden-cla.pot 4.600
!cyto-cla.pot 3.600
!guan-cla.pot 3.400
!pho-cla.pot 4.000
!suga-cla.pot 3.900
!thym-cla.pot 4.100
!dPOT values
!aden-pot.pot 3.700
!cyto-pot.pot 3.200
!guan-pot.pot 3.400
!pho-pot.pot 3.100
!suga-pot.pot 3.000
!thym-pot.pot 4.000
!Effective Radius is computed as an average of dPOT and dCLA: d=(rPOT+rCLA)/2
! aden	4.15
! cyto	3.4
! guan	3.4
! pho	3.55
! suga	3.45
! thym	4.05

end module

module gsbp
implicit none
integer*4 nclx, ncly, nclz, ncyz, mini
real*8  dcel, xbcen, ybcen, zbcen
real*8  epsw, epspp, conc, tmemb, zmemb1, epsm
real*8 tranx, trany, tranz, idcel
integer*4 ncel3
integer*4,allocatable :: lv(:)
integer*1,allocatable :: phiv(:)
integer*1,allocatable :: phivo(:)
real*8,allocatable :: surfav(:),surf(:),surfdna(:)
end module

program xssurf
implicit none
integer*4 nsc
integer narg,arg
character inpfile*256,outfile*256,line*256,bs*8
real*8 low,up
logical*1 inpopen

bs=repeat(achar(8),len(bs))

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
call readarg('BROMOC repulsion map (.pbeq) filename: ',narg,arg,inpfile)
open(1,file=inpfile,form='unformatted')
call readphi(1,6)
close(1)
call areapot(0.0,0.0)

! 5 stdin, 6 stdout, 0 stderr
call readarg('BROMOC trajectory (.btr) filename: ',narg,arg,inpfile)

call readarg('Output filename: ',narg,arg,outfile)

call readarg('Lower Limit for indexes: ',narg,arg,line)
read(line,*) low
call readarg('Upper Limit for indexes: ',narg,arg,line)
read(line,*) up

! Read btr input file
call readbtrhead(inpfile,.false.,nsc,inpopen)
call readbtrbody(nsc,inpopen)

call typeradius()

do while (inpopen)                            ! open loop for each frame
  write(*,'(A8,I8$)') bs,nsc
  call mapdna()
  call compsurfdna(0.0,0.0) 
  call compsurf(0.0,0.0) 
  if (inpopen) call readbtrbody(nsc,inpopen)          ! read next fram
enddo
call aver(nsc)
open(1,file=outfile)
call writeout(1,low,up)
close(1)
write(*,'(/A)') 'Normal termination'
end program

subroutine aver(n)
use gsbp
integer*4 n
surfav=surfav/n
surfdna=surfdna/n
return
end subroutine 

subroutine writeout(iun,low,up)
use gsbp
implicit none
integer iun,i
real*8 idx1,idx2,low,up,z
character line*1024

idx1=0.0
idx2=0.0
do i=1,nclz
  z=(i-1)*dcel-tranz+zbcen
  if (z.ge.low.and.z.le.up) then
    idx1=idx1+1.0/surf(i)
    idx2=idx2+1.0/surfav(i)
  endif
  write(line,*) z,surf(i),surfav(i),surfdna(i)
  write(1,'(A)') trim(line)
enddo
idx1=idcel/idx1
idx2=idcel/idx2
write(*,'(//)')
write(*,*) 'Lower and Upper boundaries: ',low,up
write(*,*) 'Conductance index of the pore: ',idx1
write(*,*) 'Conductance index of the pore with DNA: ',idx2
return
end subroutine

subroutine compsurf(xi,yi)
use gsbp
implicit none
integer i,j,k,ix,iy
real xi,yi
integer*4 surfz

phivo=phivo+phiv
ix=int((xi+tranx-xbcen)*idcel)
iy=int((yi+trany-ybcen)*idcel)
do k=0,nclz-1
  surfz=0
  do i=ix-lv(k+1),ix+lv(k+1)
    do j=iy-lv(k+1),iy+lv(k+1)
      if (phivo(i*ncyz+j*nclz+k+1).eq.0) surfz=surfz+1
    enddo
  enddo
  surfav(k+1)=surfav(k+1)+surfz*dcel**2
enddo
return
end subroutine

subroutine compsurfdna(xi,yi)
use gsbp
implicit none
integer i,j,k,ix,iy
real xi,yi
integer*4 surfz

ix=int((xi+tranx-xbcen)*idcel)
iy=int((yi+trany-ybcen)*idcel)
do k=0,nclz-1
  surfz=0
  do i=ix-mini,ix+mini
    do j=iy-mini,iy+mini
      if (i*ncyz+j*nclz+k+1.gt.ncel3) stop 'out of boundary' ! debug
      if (phivo(i*ncyz+j*nclz+k+1).gt.0) surfz=surfz+1
    enddo
  enddo
  surfdna(k+1)=surfdna(k+1)+surfz*dcel**2
enddo
return
end subroutine


subroutine mapdna()
use btr
use gsbp
implicit none
integer i,j,k,l
integer ix,iy,iz,md
real*4 xi,yi,zi,d2,md2

phivo=0  

do l=1,ndna
  md=rad(typ(l))*idcel+1
  md2=rad(typ(l))**2
  ix=int((rt(1,l)+tranx-xbcen)*idcel)
  iy=int((rt(2,l)+trany-ybcen)*idcel)
  iz=int((rt(3,l)+tranz-zbcen)*idcel)
  do i=ix-md,ix+md
    if (i.ge.0.and.i.lt.nclx) then
      do j=iy-md,iy+md
        if (j.ge.0.and.j.lt.ncly) then
          do k=iz-md,iz+md
            if (k.ge.0.and.k.lt.nclz) then
              xi=i*dcel-tranx+xbcen
              yi=j*dcel-trany+ybcen
              zi=k*dcel-tranz+zbcen
              d2=(xi-rt(1,l))**2+(yi-rt(2,l))**2+(zi-rt(3,l))**2
              if (d2.le.md2) phivo(i*ncyz+j*nclz+k+1)=1 !+phivo(i*ncyz+j*nclz+k+1)
            endif
          enddo
        endif
      enddo
    endif
  enddo
enddo
return
end subroutine

subroutine typeradius()
use btr
implicit none
integer i,j
ndna=0
do j=1,ntop
  i=typ(j)
  if (atnam2(i).eq.'S'.or.atnam2(i).eq.'P'.or.atnam2(i).eq.'Ab'.or.     &
      atnam2(i).eq.'Cb'.or.atnam2(i).eq.'Gb'.or.atnam2(i).eq.'Tb') then
    ndna=j
  else
    exit
  endif
enddo
do i=1,itype
  if (atnam2(i).eq.'S'.or.atnam2(i).eq.'P'.or.atnam2(i).eq.'Ab'.or.     &
      atnam2(i).eq.'Cb'.or.atnam2(i).eq.'Gb'.or.atnam2(i).eq.'Tb') then
    dtype=i
  else
    exit
  endif
enddo
allocate (rad(dtype))
do i=1,dtype
  if (atnam2(i).eq.'S') then
    rad(i)=sug  
  elseif (atnam2(i).eq.'P') then
    rad(i)=pho 
  elseif (atnam2(i).eq.'Ab') then
    rad(i)=ab 
  elseif (atnam2(i).eq.'Cb') then
    rad(i)=cb  
  elseif (atnam2(i).eq.'Gb') then
    rad(i)=gb
  elseif (atnam2(i).eq.'Tb') then
    rad(i)=tb  
  endif 
enddo
end subroutine

subroutine readbtrhead(inpfile,no,nsc,inpopen)
use btr
implicit none
integer i, j, kode
integer*4 recordi,recordl,nsc
character*256 inpfile
logical*1 no,inpopen

open(14,file=inpfile,form='unformatted',iostat=kode,position='rewind',access='stream')
read(14,iostat=kode) recordi
read(14,iostat=kode) nframe                 ! number of frames
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

read(14,iostat=kode) recordi
read(14,iostat=kode) itype                  ! number of ions and nucleotides
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (.not.allocated(atnam2)) allocate (atnam2(itype))
if (.not.allocated(maxntyp)) allocate (maxntyp(itype))

read(14,iostat=kode) recordi
read(14,iostat=kode) (atnam2(j),j=1,itype)  ! ion and nucleotides types in char
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (.not.no) then
  write(*,'(A,I0)') 'Number of frames: ',nframe
  write(*,'(A,I0)') 'Number of fragment types: ',itype
  write(*,*) 'Fragment types: ',(atnam2(j),j=1,itype)
endif
nsc=0
if (kode.eq.0) then
  inpopen=.true.
else
  close(14)
endif
end subroutine

! check fortran binary record consistency
subroutine check(recordi,recordl,un,kode)
implicit none
integer*4 recordi,recordl,pos,un,kode
character*64 fname

if (kode.eq.0) then
  if (recordi.ne.recordl) then
    inquire(unit=un,pos=pos,name=fname)
    write (*,*)
    write (*,*) 'Unit: ',un
    write (*,*) 'Filename: ',trim(fname)
    write (*,*) 'Record not matching (init, last): ',recordi,recordl
    write (*,*) 'Position of last byte read: ',pos
    stop
  endif
endif
end subroutine

subroutine readbtrbody(nsc,inpopen)
use btr
implicit none
integer kode,i
integer*4 recordi,recordl,nsc
logical*1 inpopen

inquire(unit=14,pos=ir)

if (allocated(typ)) deallocate (typ)
if (allocated(rt)) deallocate (rt)
read(14,iostat=kode) recordi
read(14,iostat=kode) runtime !simulation time for each step (ns)
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

read(14,iostat=kode) recordi
read(14,iostat=kode) ntop ! number of particles
read(14,iostat=kode) recordl
call check(recordi,recordl,14,kode)

if (kode.eq.0) then 
  allocate (typ(ntop),rt(3,ntop))

  read(14,iostat=kode) recordi
  read(14,iostat=kode) (typ(i),i=1,ntop) ! ion and DNA sites types
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
  
  read(14,iostat=kode) recordi
  read(14,iostat=kode) (rt(1,i),i=1,ntop)
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
  
  read(14,iostat=kode) recordi
  read(14,iostat=kode) (rt(2,i),i=1,ntop)
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
  
  read(14,iostat=kode) recordi
  read(14,iostat=kode) (rt(3,i),i=1,ntop)
  read(14,iostat=kode) recordl
  call check(recordi,recordl,14,kode)
endif

inquire(unit=14,pos=lr)

if (kode.eq.0) then
  nsc=nsc+1
!  write(*,*) runtime,ntop,(typ(i),i=1,ntop)
else
  close(14)
  inpopen=.false.
endif
end subroutine

! read arguments
subroutine readarg(ques,narg,num,text)
implicit none
integer*4 narg,num
character text*(*),ques*(*)

num=num+1
write(*,'(/A$)') ques
if (narg.ge.num) then
  call GET_COMMAND_ARGUMENT(num,text)
  write(*,'(A)') trim(text)
else
  text=''
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='XS-SURF'
prver='version 1.0'
prdesc='Computes the cross-section free surface along z-axis considering the pore and DNA'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='03 Dec 2014'
lastdate='03 Dec 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

subroutine readphi(unit,outu)
!-----------------------------------------------------------------------
! read INPUT static external field PHIX or grid-based repulsion potential PHIV
!            and miscelaneous parameters
use gsbp
!Input
implicit none
integer*4 unit, outu
integer*4 i
real*4,allocatable ::  phivt(:)
integer*4 ifir, ilas

read(unit) nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen
read(unit) epsw,epspp,conc,tmemb,zmemb1,epsm
tranx = 0.5*(nclx-1)*dcel
trany = 0.5*(ncly-1)*dcel
tranz = 0.5*(nclz-1)*dcel
ncel3  = nclx*ncly*nclz
idcel = 1.0/dcel
ncyz=ncly*nclz

!Writting in output file            
write(outu,*)
write(outu,*) 'Number of grid point in X   (nclx) = ',nclx 
write(outu,*) 'Number of grid point in Y   (ncly) = ',ncly 
write(outu,*) 'Number of grid point in Z   (nclz) = ',nclz 
write(outu,*) 'Grid spacing                (dcel) = ',dcel
write(outu,*) 'Center of box in X          (xbcen)= ',xbcen
write(outu,*) 'Center of box in Y          (ybcen)= ',ybcen
write(outu,*) 'Center of box in Z          (zbcen)= ',zbcen
write(outu,*) 
write(outu,*) 'Solvent dielectric constant (epsw) = ',epsw
write(outu,*) 'Protein dielectric constant (epsp) = ',epspp
write(outu,*) 'Salt concentration          (conc) = ',conc
if (tmemb.gt.0.0) then
  write(outu,*)
  write(outu,*) 'Membrane thickness along Z  (tmemb)= ',tmemb
  write(outu,*) 'Membrane position along Z   (zmemb)= ',zmemb1
  write(outu,*) 'Membrane dielectric constant(epsm) = ',epsm
endif
write(outu,*)
write(outu,*) 'Box in X from ',xbcen-tranx,' to ',xbcen+tranx
write(outu,*) 'Box in Y from ',ybcen-trany,' to ',ybcen+trany
write(outu,*) 'Box in Z from ',zbcen-tranz,' to ',zbcen+tranz

!Grid-based repulsion potential        
ifir = 1
ilas = ncel3
if (allocated(phiv)) deallocate (phiv)
allocate(phivt(ifir:ilas))
allocate(phiv(ifir:ilas))
allocate(phivo(ifir:ilas))
read(unit) (phivt(i),i=ifir,ilas)
do i = ifir, ilas
  if (phivt(i).ne.0.0) then
    phiv(i) = 0
  else
    phiv(i) = 1
  endif
enddo
deallocate (phivt)
return
end subroutine

subroutine areapot(xi,yi)
use gsbp
implicit none
integer i,j,k,l,ix,iy,jjmp
real xi,yi
logical*1 ety

allocate (lv(nclz),surfav(nclz),surf(nclz),surfdna(nclz))
lv=0
surf=0.0
surfdna=0.0
surfav=0.0
ix=int((xi+tranx-xbcen)*idcel)
iy=int((yi+trany-ybcen)*idcel)
mini=min(nclx-1,ncly-1)/2
do k=0,nclz-1
  l=0
  ety=.true.
  surf(k+1)=0.0
  do while (l.le.mini.and.ety)
    ety=.false.
    do i=ix-l,ix+l
      if (i.eq.ix-l.or.i.eq.ix+l) then
        jjmp=1
      else
        jjmp=2*l
        if (l.eq.0) jjmp=1
      endif
      do j=iy-l,iy+l,jjmp
        if (phiv(i*ncyz+j*nclz+k+1).eq.0) then
          ety=.true.
          surf(k+1)=surf(k+1)+1.0
        endif
      enddo
    enddo
    l=l+1
  enddo
  lv(k+1)=l-1
  surf(k+1)=surf(k+1)*dcel**2
enddo
return
end subroutine


