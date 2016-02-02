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

subroutine readphi(unit,outu,wrd4,adjust)
!-----------------------------------------------------------------------
! read INPUT static external field PHIX or grid-based repulsion potential PHIV
!            and miscelaneous parameters
use ioxmod
use gsbpmod
use grandmod
use nucleotmod
use errormod
!Input
implicit none
integer unit, outu
character*4  wrd4
!Local variables
integer i, ncel3, numb, j, k, x1, x2, y1, y2, z1, z2
integer*4 nclx, ncly, nclz, nxo, nyo, nzo, pos, posn
integer ifir, ilas, ifir0, isite, itype
real tranx, trany, tranz
real*8  xbcen, ybcen, zbcen, dcel,xmmn,xmmm,ymmn,ymmm,zmmn,zmmm
real*8  epsw, epspp, conc
real*8  tmembr, epsmr, zmemb1r
logical*1 adjust
real*4,allocatable ::  phixt(:),phivt(:)

read(unit) nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen
read(unit) epsw,epspp,conc,tmembr,zmemb1r,epsmr
tmemb=tmembr
zmemb1=zmemb1r
epsm=epsmr

tranx = 0.5*(nclx-1)*dcel
trany = 0.5*(ncly-1)*dcel
tranz = 0.5*(nclz-1)*dcel
ncel3  = nclx*ncly*nclz

write(outu,'(6x,/A)') '         SYSTEM-SIZE                        MAP-DIMENSION'
write(outu,'(6x,A,2(F12.5,A,F12.5,5x))') 'x ',lx2m,' - ',lx2p,xbcen-tranx,' - ',xbcen+tranx
write(outu,'(6x,A,2(F12.5,A,F12.5,5x))') 'y ',ly2m,' - ',ly2p,ybcen-trany,' - ',ybcen+trany
write(outu,'(6x,A,2(F12.5,A,F12.5,5x))') 'z ',lz2m,' - ',lz2p,zbcen-tranz,' - ',zbcen+tranz
if (lx2m.lt.xbcen-tranx.or.lx2p.gt.xbcen+tranx.or.ly2m.lt.ybcen-trany.or.ly2p.gt.ybcen+trany.or.lz2m.lt.zbcen-tranz.or.lz2p.gt.zbcen+tranz) then
  write(outu,'(6x,a)') 'SYSTEM-DIMENSION does not fit into MAP-DIMENSION'
else 
  write(outu,'(6x,a)') 'SYSTEM-DIMENSION fits into MAP-DIMENSION'
endif
if (lx2m.gt.xbcen-tranx.or.lx2p.lt.xbcen+tranx.or.ly2m.gt.ybcen-trany.or.ly2p.lt.ybcen+trany.or.lz2m.gt.zbcen-tranz.or.lz2p.lt.zbcen+tranz) then
  write(outu,'(6x,a)') 'MAP-DIMENSION does not fit into SYSTEM-DIMENSION'
else 
  write(outu,'(6x,a)') 'MAP-DIMENSION fits into SYSTEM-DIMENSION'
endif


!Writting in output file            
write(outu,*)
write(outu,101) 'Number of grid point in X   (nclx) = ',nclx 
write(outu,101) 'Number of grid point in Y   (ncly) = ',ncly 
write(outu,101) 'Number of grid point in Z   (nclz) = ',nclz 
write(outu,102) 'Grid spacing                (dcel) = ',dcel
write(outu,102) 'Center of box in X          (xbcen)= ',xbcen
write(outu,102) 'Center of box in Y          (ybcen)= ',ybcen
write(outu,102) 'Center of box in Z          (zbcen)= ',zbcen
write(outu,*) 
write(outu,102) 'Solvent dielectric constant (epsw) = ',epsw
write(outu,102) 'Protein dielectric constant (epsp) = ',epspp
write(outu,102) 'Salt concentration          (conc) = ',conc
if (tmemb.gt.0.0) then
  write(outu,*)
  write(outu,102) 'Membrane thickness along Z  (tmemb)= ',tmemb
  write(outu,102) 'Membrane position along Z   (zmemb)= ',zmemb1
  write(outu,102) 'Membrane dielectric constant(epsm) = ',epsm
endif
write(outu,*)
write(outu,103) 'Box in X from ',xbcen-tranx,' to ',xbcen+tranx
write(outu,103) 'Box in Y from ',ybcen-trany,' to ',ybcen+trany
write(outu,103) 'Box in Z from ',zbcen-tranz,' to ',zbcen+tranz

!Static external field          
if (wrd4.eq.'PHIX') then
  if (allocated(phix)) deallocate (phix)
! adjust not actually working when center is not zero for system and maps
  if (adjust) then
    allocate (phixt(ncel3))
    read(unit) (phixt(i),i=1,ncel3)
    write(outu,'(/6x,a)') 'Resizing map to the system limits to optimize memory usage'
    write(outu,'(6x,a)') 'New boundaries:' 
    x1=int((lx2m+tranx-xbcen)/dcel)
    y1=int((ly2m+trany-ybcen)/dcel)
    z1=int((lz2m+tranz-zbcen)/dcel)
    x2=int((lx2p+tranx-xbcen)/dcel)+2
    y2=int((ly2p+trany-ybcen)/dcel)+2
    z2=int((lz2p+tranz-zbcen)/dcel)+2
    if (x1.lt.1) x1=1
    if (y1.lt.1) y1=1
    if (z1.lt.1) z1=1
    if (x2.gt.nclx) x2=nclx
    if (y2.gt.ncly) y2=ncly
    if (z2.gt.nclz) z2=nclz
    xmmn=(x1-1)*dcel-tranx+xbcen
    ymmn=(y1-1)*dcel-trany+ybcen
    zmmn=(z1-1)*dcel-tranz+zbcen
    xmmm=(x2-1)*dcel-tranx+xbcen
    ymmm=(y2-1)*dcel-trany+ybcen
    zmmm=(z2-1)*dcel-tranz+zbcen
    nxo=nclx
    nyo=ncly
    nzo=nclz
    nclx=x2-x1+1
    ncly=y2-y1+1
    nclz=z2-z1+1
    ncel3 = nclx*ncly*nclz
    tranx = 0.5*(nclx-1)*dcel
    trany = 0.5*(ncly-1)*dcel
    tranz = 0.5*(nclz-1)*dcel
    xbcen = 0.5*(xmmm+xmmn) 
    ybcen = 0.5*(ymmm+ymmn) 
    zbcen = 0.5*(zmmm+zmmn) 
    write(outu,101) 'Number of grid point in X   (nclx) = ',nclx
    write(outu,101) 'Number of grid point in Y   (ncly) = ',ncly
    write(outu,101) 'Number of grid point in Z   (nclz) = ',nclz
    write(outu,*)
    write(outu,102) 'Center of box in X          (xbcen)= ',xbcen
    write(outu,102) 'Center of box in Y          (ybcen)= ',ybcen
    write(outu,102) 'Center of box in Z          (zbcen)= ',zbcen
    write(outu,*)
    write(outu,103) 'Box in X from ',xbcen-tranx,' to ',xbcen+tranx
    write(outu,103) 'Box in Y from ',ybcen-trany,' to ',ybcen+trany
    write(outu,103) 'Box in Z from ',zbcen-tranz,' to ',zbcen+tranz
    allocate (phix(ncel3))
    posn=0
    do i=x1,x2
      do j=y1,y2
        do k=z1,z2
          pos=(i-1)*nyo*nzo + (j-1)*nzo + k
          posn=posn+1
          phix(posn)=phixt(pos)
        enddo
      enddo
    enddo
    deallocate(phixt)
  else
    allocate(phix(ncel3))
    read(unit) (phix(i),i=1,ncel3)
  endif
  nclx1 = nclx
  ncly1 = ncly
  nclz1 = nclz
  dcel1 = dcel
  idcel1 = 1.0/dcel1
  tranx1 = tranx
  trany1 = trany
  tranz1 = tranz
  xbcen1 = xbcen
  ybcen1 = ybcen
  zbcen1 = zbcen
endif

!Grid-based repulsion potential        
if (wrd4.eq.'PHIV') then
  if (Qnmcden) then
    if (adjust) write(outu,'(6x,a)') 'Warning: Cannot resize map, not implemented for nmcden'
    ifir=1
    ifir0=0
    numb=0
    if (Qnucl) then
      ilas = (6-1)*ncel3 + ncel3
      ifir0 = ilas
    endif
    if (Qpar) then
      if (Qnucl .and. istrs.eq.1) numb = inuc + 1
      if (Qnucl .and. istrs.eq.2) numb = 2*inuc + 1
      ifir = ifir0 + (ntype-numb)*ncel3 + 1
      ilas = ifir0 + (ntype-numb)*ncel3 + ncel3 
    endif
    if (allocated(phiv)) deallocate (phiv)
    allocate (phiv(ifir:ilas),phivt(ifir:ilas))
    if (Qnucl) then 
      do isite = 1, 6 ! sites (S,P,Ab,Tb,Cb,Gb)
        ifir = (isite-1)*ncel3 + 1
        ilas = (isite-1)*ncel3 + ncel3
        unit = vecphiv(isite) 
        if (ifir.eq.1) then
          read(unit) (phivt(i),i=ifir,ilas)
        else
          read(unit) nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen
          read(unit) epsw,epspp,conc,tmembr,zmemb1r,epsmr
          tmemb=tmembr
          zmemb1=zmemb1r
          epsm=epsmr
          read(unit) (phivt(i),i=ifir,ilas)
        endif      
        do i = ifir, ilas
          if (phivt(i).ne.0.0) then ! BAJO CHISTERA!!!
            phiv(i) = 0
          else
            phiv(i) = 1
          endif
        enddo                
      enddo
      ifir0 = ilas
    endif ! Qnucl
    if (Qpar) then
      numb = 1   
      if (Qnucl .and. istrs.eq.1) numb = inuc + 1
      if (Qnucl .and. istrs.eq.2) numb = 2*inuc + 1  
      j = 0
      if (Qnucl) j = 6
      do itype = numb, ntype
        j = j + 1
        ifir = ifir0 + (itype-numb)*ncel3 + 1
        ilas = ifir0 + (itype-numb)*ncel3 + ncel3
        unit = vecphiv(j)  
        if (ifir.eq.1) then
          read(unit) (phivt(i),i=ifir,ilas)
        else
          read(unit) nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen
          read(unit) epsw,epspp,conc,tmembr,zmemb1r,epsmr
          tmemb=tmembr
          zmemb1=zmemb1r
          epsm=epsmr
          read(unit) (phivt(i),i=ifir,ilas)
        endif
        ! BAJO CHISTERA
        do i = ifir, ilas
          if (phivt(i).ne.0.0) then
            phiv(i) = 0
          else
            phiv(i) = 1
          endif
        enddo
      enddo
    endif ! Qpar
    deallocate (phivt)
  else
    ifir = 1
    ilas = ncel3
    if (allocated(phiv)) deallocate (phiv)
    allocate (phivt(ifir:ilas))
    if (adjust) then
      read(unit) (phivt(i),i=ifir,ilas)
      write(outu,'(/6x,a)') 'Resizing map to the system limits to optimize memory usage'
      write(outu,'(6x,a)') 'New boundaries:'
      x1=int((lx2m+tranx-xbcen)/dcel)
      y1=int((ly2m+trany-ybcen)/dcel)
      z1=int((lz2m+tranz-zbcen)/dcel)
      x2=int((lx2p+tranx-xbcen)/dcel)+2
      y2=int((ly2p+trany-ybcen)/dcel)+2
      z2=int((lz2p+tranz-zbcen)/dcel)+2
      if (x1.lt.1) x1=1
      if (y1.lt.1) y1=1
      if (z1.lt.1) z1=1
      if (x2.gt.nclx) x2=nclx
      if (y2.gt.ncly) y2=ncly
      if (z2.gt.nclz) z2=nclz
      xmmn=(x1-1)*dcel-tranx+xbcen
      ymmn=(y1-1)*dcel-trany+ybcen
      zmmn=(z1-1)*dcel-tranz+zbcen
      xmmm=(x2-1)*dcel-tranx+xbcen
      ymmm=(y2-1)*dcel-trany+ybcen
      zmmm=(z2-1)*dcel-tranz+zbcen
      nxo=nclx
      nyo=ncly
      nzo=nclz
      nclx=x2-x1+1
      ncly=y2-y1+1
      nclz=z2-z1+1
      ncel3 = nclx*ncly*nclz
      tranx = 0.5*(nclx-1)*dcel
      trany = 0.5*(ncly-1)*dcel
      tranz = 0.5*(nclz-1)*dcel
      xbcen = 0.5*(xmmm+xmmn)
      ybcen = 0.5*(ymmm+ymmn)
      zbcen = 0.5*(zmmm+zmmn)
      write(outu,101) 'Number of grid point in X   (nclx) = ',nclx
      write(outu,101) 'Number of grid point in Y   (ncly) = ',ncly
      write(outu,101) 'Number of grid point in Z   (nclz) = ',nclz
      write(outu,*)
      write(outu,102) 'Center of box in X          (xbcen)= ',xbcen
      write(outu,102) 'Center of box in Y          (ybcen)= ',ybcen
      write(outu,102) 'Center of box in Z          (zbcen)= ',zbcen
      write(outu,*)
      write(outu,103) 'Box in X from ',xbcen-tranx,' to ',xbcen+tranx
      write(outu,103) 'Box in Y from ',ybcen-trany,' to ',ybcen+trany
      write(outu,103) 'Box in Z from ',zbcen-tranz,' to ',zbcen+tranz
      ifir=1
      ilas=ncel3
      allocate (phiv(ifir:ilas))
      posn=0
      do i=x1,x2
        do j=y1,y2
          do k=z1,z2
            pos=(i-1)*nyo*nzo + (j-1)*nzo + k
            posn=posn+1
            if (phivt(pos).ne.0.0) then 
              phiv(posn)= 0
            else
              phiv(posn)= 1
            endif
          enddo
        enddo
      enddo
    else
      allocate(phiv(ifir:ilas))
      read(unit) (phivt(i),i=ifir,ilas)
      do i = ifir, ilas
        if (phivt(i).ne.0.0) then
          phiv(i) = 0
        else
          phiv(i) = 1
        endif
      enddo
    endif
    deallocate (phivt)
  endif
  nclx2 = nclx
  ncly2 = ncly
  nclz2 = nclz
  dcel2 = dcel
  idcel2 = 1.0/dcel2
  tranx2 = tranx
  trany2 = trany
  tranz2 = tranz
  xbcen2 = xbcen
  ybcen2 = ybcen
  zbcen2 = zbcen
endif
!Formats
101  format(6x,a,i6)
102  format(6x,a,f8.3)       
103  format(6x,a,f8.3,a,f8.3)

return
end
