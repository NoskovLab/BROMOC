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

subroutine staticplot(x1,y1,z1,x2,y2,z2,res,unitn,Qnohead)
!-----------------------------------------------------------------------
!This subroutine prints the electrostatic potential resulting from the static field
use gsbpmod

!local variables
implicit none
integer ncyz,i,ix,iy,iz,n1,in1,n2,in2,n3,in3
real  xi,yi,zi,xii,yii,zii
real x1,y1,z1,x2,y2,z2,res,pos,dist,idist
real stat,st(0:1,0:1,0:1),dx,dy,dz,st2(0:1,0:1),st3(0:1)
integer unitn,didi
logical*1 ok,Qnohead
character line*1024

ncyz = ncly1*nclz1

!Main loop by atoms
dist=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
idist=1.0/dist
didi=int(dist/res)+1

if (.not.Qnohead) then
  write(unitn,*) 'Writing Static Field ...'
  write(unitn,*) 'From: ',x1,y1,z1
  write(unitn,*) 'To: ',x2,y2,z2
  write(unitn,*) 'Resolution: ',res
  write(unitn,'(a)') '       Position                 StaticField                x                       y                       z'
endif

do i = 1, didi
  pos=(i-1)*res ! position along vector
  xii=(x2-x1)*idist*pos+x1  ! x position
  yii=(y2-y1)*idist*pos+y1  ! y position
  zii=(z2-z1)*idist*pos+z1  ! z position
  ok=xii.le.xbcen1+tranx1.and.xii.ge.xbcen1-tranx1.and. &
     yii.le.ybcen1+trany1.and.yii.ge.ybcen1-trany1.and. &
     zii.le.zbcen1+tranz1.and.zii.ge.zbcen1-tranz1
  if (ok) then
    xi = xii + tranx1-xbcen1  !
    yi = yii + trany1-ybcen1  ! moving center
    zi = zii + tranz1-zbcen1  !
    ix = int(xi*idcel1) 
    iy = int(yi*idcel1) 
    iz = int(zi*idcel1)
    if (ix.eq.nclx1-1) ix=nclx1-2
    if (iy.eq.ncly1-1) iy=ncly1-2
    if (iz.eq.nclz1-1) iz=nclz1-2
    dx=xi*idcel1-ix ! save (x-x0)/dcel
    dy=yi*idcel1-iy ! save (y-y0)/dcel
    dz=zi*idcel1-iz ! save (z-z0)/dcel
    do n1=0,1
      in1=(ix+n1)*ncyz
      do n2=0,1
        in2=(iy+n2)*nclz1
        do n3=0,1
          in3=in1+in2+iz+n3+1
          st(n1,n2,n3)=phix(in3) ! save 8 neighbours static field values
        enddo
      enddo
    enddo
    ! linearly extrapolate sf values for x
    do n2=0,1
      do n3=0,1
        st2(n2,n3)=(st(1,n2,n3)-st(0,n2,n3))*dx+st(0,n2,n3)
      enddo
    enddo
    ! linearly extrapolate sf values for y
    do n3=0,1
      st3(n3)=(st2(1,n3)-st2(0,n3))*dy+st2(0,n3)
    enddo
    ! linearly extrapolate sf value for z 
    stat=(st3(1)-st3(0))*dz+st3(0)  
    write(line,*) pos,stat,xii,yii,zii
    write(unitn,'(a)') trim(line)
  endif
enddo 

return
end subroutine

subroutine repulplot(x1,y1,z1,x2,y2,z2,res,unitn,Qnohead)
!-----------------------------------------------------------------------
!This subroutine prints the electrostatic potential resulting from the static field

use gsbpmod

!local variables
implicit none
integer ncyz,i,ix,iy,iz,n1,in1,n2,in2,n3,in3,k,l,m,ipx,ipy,ipz
real  xi,yi,zi,ai,bi,ci,fi,xii,yii,zii,m3,xc,yc,zc
real x1,y1,z1,x2,y2,z2,res,pos,dist,idist
real stat,st(0:1,0:1,0:1),dx,dy,dz,st2(0:1,0:1),st3(0:1),e3
integer unitn,didi
logical*1 ok,Qnohead
character line*1024

ncyz = ncly2*nclz2

!Main loop by atoms
dist=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
idist=1.0/dist
didi=int(dist/res)+1
if (.not.Qnohead) then
  write(unitn,*) 'Writing Repulsion Field ...'
  write(unitn,*) 'From: ',x1,y1,z1
  write(unitn,*) 'To: ',x2,y2,z2
  write(unitn,*) 'Resolution: ',res
  write(unitn,'(a)') 'Position                 RepulsionField-Trilinear RepulsionField-BSpline          x                        y                         z'
endif

do i = 1, didi
  pos=(i-1)*res ! position along vector
  xii=(x2-x1)*idist*pos+x1  ! x position
  yii=(y2-y1)*idist*pos+y1  ! y position
  zii=(z2-z1)*idist*pos+z1  ! z position
  ok=xii.le.xbcen2+tranx2-dcel2.and.xii.ge.xbcen2-tranx2+dcel2.and. &
     yii.le.ybcen2+trany2-dcel2.and.yii.ge.ybcen2-trany2+dcel2.and. &
     zii.le.zbcen2+tranz2-dcel2.and.zii.ge.zbcen2-tranz2+dcel2
  if (ok) then
    xi = xii + tranx2-xbcen2  !
    yi = yii + trany2-ybcen2  ! moving center
    zi = zii + tranz2-zbcen2  !
    ix = int(xi*idcel2)
    iy = int(yi*idcel2)
    iz = int(zi*idcel2)
    if (ix.eq.nclx2-1) ix=nclx2-2
    if (iy.eq.ncly2-1) iy=ncly2-2
    if (iz.eq.nclz2-1) iz=nclz2-2
    ! Method 1 (trilinear but Pablo's algorithm)
    dx=xi*idcel2-ix ! save (x-x0)/dcel
    dy=yi*idcel2-iy ! save (y-y0)/dcel
    dz=zi*idcel2-iz ! save (z-z0)/dcel
    do n1=0,1
      in1=(ix+n1)*ncyz
      do n2=0,1
        in2=(iy+n2)*nclz2
        do n3=0,1
          in3=in1+in2+iz+n3+1
          st(n1,n2,n3)=phiv(in3) ! save 8 neighbours static field values
        enddo
      enddo
    enddo
    ! linearly extrapolate sf values for x
    do n2=0,1
      do n3=0,1
        st2(n2,n3)=(st(1,n2,n3)-st(0,n2,n3))*dx+st(0,n2,n3)
      enddo
    enddo
    ! linearly extrapolate sf values for y
    do n3=0,1
      st3(n3)=(st2(1,n3)-st2(0,n3))*dy+st2(0,n3)
    enddo
    ! linearly extrapolate sf value for z 
    stat=(st3(1)-st3(0))*dz+st3(0)

!    ! Method 2 = Method 1 (trilinear, original algorithm) 
!    e2=0.0
!    do n1 = 1, 2
!      in1 = ix + n1
!      ai = xi -(in1-1)*dcel2
!      ai = 1.0 - abs(ai)*idcel2
!      in1 = (in1-1)*ncyz
!      do n2 = 1, 2
!        in2 = iy + n2
!        bi = yi - (in2-1)*dcel2
!        bi = 1.0 - abs(bi)*idcel2
!        in2 = (in2-1)*nclz2
!        do n3=1,2
!          in3 = iz + n3
!          ci = zi - (in3-1)*dcel2
!          ci = 1.0 - abs(ci)*idcel2
!          fi = ai*bi*ci
!          in3 = in1 + in2 + in3
!          e2 = e2 + fi*phiv(in3)
!        enddo ! n3
!      enddo ! n2
!    enddo ! n1
    
    ! Method 3 (3rd order b-spline)
    ix=nint(xi*idcel2)
    iy=nint(yi*idcel2)
    iz=nint(zi*idcel2)
    if(ix.eq.0) ix=1
    if(iy.eq.0) iy=1
    if(iz.eq.0) iz=1
    if(ix.eq.nclx2-1)ix=nclx2-2
    if(iy.eq.ncly2-1)iy=ncly2-2
    if(iz.eq.nclz2-1)iz=nclz2-2
    e3=0.0
    do k = ix-1, ix+1
      ipx = k*ncyz
      xc = k*dcel2
      ai = 1.5 - (xi-xc)*idcel2
      ai = m3(ai)
      if (ai.ne.0.0) then
        do l = iy-1, iy+1
          ipy = l*nclz2 + ipx
          yc = l*dcel2
          bi = 1.5 - (yi-yc)*idcel2
          bi = m3(bi)
          if (bi.ne.0.0) then
            do m = iz-1, iz+1
              ipz = m + ipy + 1
              zc = m*dcel2
              ci = 1.5 - (zi-zc)*idcel2
              fi = ai*bi*m3(ci)
              if (fi.ne.0.0) then
                e3 = e3 + fi*phiv(ipz)
              endif
            enddo
          endif
        enddo
      endif
    enddo
    write(line,*) pos,stat,e3,xii,yii,zii
    write(unitn,'(a)') trim(line)
  endif
enddo 

return
end subroutine

subroutine rfparplot(x1,y1,z1,x2,y2,z2,res,unitn,Qnohead,nion)
!-----------------------------------------------------------------------
!This subroutine prints the electrostatic potential resulting from the static field

use gsbpmod

!local variables
implicit none
integer ncyz,i,ix,iy,iz,n1,in1,n2,in2,n3,in3,itype,ncel3,nion
real  xi,yi,zi,xii,yii,zii
real x1,y1,z1,x2,y2,z2,res,pos,dist,idist
real stat(nion,2),st(0:1,0:1,0:1,nion,2),dx,dy,dz,st2(0:1,0:1,nion,2),st3(0:1,nion,2)
integer unitn,didi
logical*1 ok,Qnohead
character line*1024

ncyz = ncly3*nclz3
ncel3= nclx3*ncyz

!Main loop by atoms
dist=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
idist=1.0/dist
didi=int(dist/res)+1

if (.not.Qnohead) then
  write(unitn,*) 'Writing Reaction Field ...'
  write(unitn,*) 'From: ',x1,y1,z1
  write(unitn,*) 'To: ',x2,y2,z2
  write(unitn,*) 'Resolution: ',res
  write(unitn,'(a)') '       Position              SRFEN ion 1 SRFEN ion 2 ... REFF ion1 REFF ion2 ...           x                       y                       z'
endif

do i = 1, didi
  pos=(i-1)*res ! position along vector
  xii=(x2-x1)*idist*pos+x1  ! x position
  yii=(y2-y1)*idist*pos+y1  ! y position
  zii=(z2-z1)*idist*pos+z1  ! z position
  ok=xii.le.xbcen3+tranx3.and.xii.ge.xbcen3-tranx3.and. &
     yii.le.ybcen3+trany3.and.yii.ge.ybcen3-trany3.and. &
     zii.le.zbcen3+tranz3.and.zii.ge.zbcen3-tranz3
  if (ok) then
    xi = xii + tranx3-xbcen3  !
    yi = yii + trany3-ybcen3  ! moving center
    zi = zii + tranz3-zbcen3  !
    ix = int(xi*idcel3) 
    iy = int(yi*idcel3) 
    iz = int(zi*idcel3)
    if(ix.eq.nclx3-1)ix=nclx3-2
    if(iy.eq.ncly3-1)iy=ncly3-2
    if(iz.eq.nclz3-1)iz=nclz3-2
    dx=xi*idcel3-ix ! save (x-x0)/dcel
    dy=yi*idcel3-iy ! save (y-y0)/dcel
    dz=zi*idcel3-iz ! save (z-z0)/dcel
    do n1=0,1
      in1=(ix+n1)*ncyz
      do n2=0,1
        in2=(iy+n2)*nclz3
        do n3=0,1
          in3=in1+in2+iz+n3+1
          do itype=1,nion
            st(n1,n2,n3,itype,1)=gsrfen(in3+(itype-1)*ncel3)
            st(n1,n2,n3,itype,2)=greff(in3+(itype-1)*ncel3)
          enddo
        enddo
      enddo
    enddo
    ! linearly extrapolate sf values for x
    do n2=0,1
      do n3=0,1
        do itype=1,nion
          st2(n2,n3,itype,1)=(st(1,n2,n3,itype,1)-st(0,n2,n3,itype,1))*dx+st(0,n2,n3,itype,1)
          st2(n2,n3,itype,2)=(st(1,n2,n3,itype,2)-st(0,n2,n3,itype,2))*dx+st(0,n2,n3,itype,2)
        enddo
      enddo
    enddo
    ! linearly extrapolate sf values for y
    do n3=0,1
      do itype=1,nion
        st3(n3,itype,1)=(st2(1,n3,itype,1)-st2(0,n3,itype,1))*dy+st2(0,n3,itype,1)
        st3(n3,itype,2)=(st2(1,n3,itype,2)-st2(0,n3,itype,2))*dy+st2(0,n3,itype,2)
      enddo
    enddo
    ! linearly extrapolate sf value for z 
    do itype=1,nion
      stat(itype,1)=(st3(1,itype,1)-st3(0,itype,1))*dz+st3(0,itype,1)  
      stat(itype,2)=(st3(1,itype,2)-st3(0,itype,2))*dz+st3(0,itype,2)  
    enddo
    write(line,*) pos,stat(1:nion,1:2),xii,yii,zii
    write(unitn,'(a)') trim(line)
  endif
enddo 

return
end subroutine

subroutine areapot(xi,yi,unitn)
use gsbpmod
use grandmod

!local variables
implicit none
integer ncyz,i,j,k,l,ix,iy,jjmp
real xi,yi,stat,sts(nclz2),ens(nclz2),ene1,ene2
integer unitn,mini,surfz(nclz2)
character line*1024
logical*1 ety

ncyz = ncly2*nclz2
ix=int((xi+tranx2-xbcen2)*idcel2)
iy=int((yi+trany2-ybcen2)*idcel2)
mini=min(nclx2-1,ncly2-1)/2
do k=0,nclz2-1
  l=0
  ety=.true.
  surfz(k+1)=0
  sts(k+1)=0.0
  ens(k+1)=0.0
  do while (l<mini.and.ety)
    ety=.false.
    do i=ix-l,ix+l
      if (i.eq.ix-l.or.i.eq.ix+l) then
        jjmp=1
      else
        jjmp=2*l
        if (l.eq.0) jjmp=1
      endif
      do j=iy-l,iy+l,jjmp
        if (phiv(i*ncyz+j*nclz2+k+1).eq.0) then
          ety=.true.
          surfz(k+1)=surfz(k+1)+1
          call getstatic(i*dcel2-tranx2+xbcen2,j*dcel2-trany2+ybcen2,k*dcel2-tranz2+zbcen2,stat)
          sts(k+1)=sts(k+1)+stat
          call staticener(i*dcel2-tranx2+xbcen2,j*dcel2-trany2+ybcen2,k*dcel2-tranz2+zbcen2-dcel2/2,ene1)
          call staticener(i*dcel2-tranx2+xbcen2,j*dcel2-trany2+ybcen2,k*dcel2-tranz2+zbcen2+dcel2/2,ene2)
          ens(k+1)=ens(k+1)+exp(-((ene2-ene1)*ikBT)**2)
        endif
      enddo
    enddo 
    l=l+1
  enddo
  write(line,*) k+1,surfz(k+1),surfz(k+1)*dcel2**2,sts(k+1)/surfz(k+1),ens(k+1)*dcel2,xi,yi,k*dcel2-tranz2+zbcen2
  write(unitn,'(a)') trim(line)
enddo 
return
end subroutine

subroutine getstatic(xii,yii,zii,stat)
use gsbpmod
integer ncyz,ix,iy,iz,n1,in1,n2,in2,n3,in3
real xi,yi,zi,xii,yii,zii
real stat,st(0:1,0:1,0:1),dx,dy,dz,st2(0:1,0:1),st3(0:1)
logical*1 ok

stat=0.0
ncyz = ncly1*nclz1
ok=xii.le.xbcen1+tranx1.and.xii.ge.xbcen1-tranx1.and. &
   yii.le.ybcen1+trany1.and.yii.ge.ybcen1-trany1.and. &
   zii.le.zbcen1+tranz1.and.zii.ge.zbcen1-tranz1
if (ok) then
  xi = xii + tranx1-xbcen1  !
  yi = yii + trany1-ybcen1  ! moving center
  zi = zii + tranz1-zbcen1  !
  ix = int(xi*idcel1)
  iy = int(yi*idcel1)
  iz = int(zi*idcel1)
  if (ix.eq.nclx1-1) ix=nclx1-2
  if (iy.eq.ncly1-1) iy=ncly1-2
  if (iz.eq.nclz1-1) iz=nclz1-2
  dx=xi*idcel1-ix ! save (x-x0)/dcel
  dy=yi*idcel1-iy ! save (y-y0)/dcel
  dz=zi*idcel1-iz ! save (z-z0)/dcel
  do n1=0,1
    in1=(ix+n1)*ncyz
    do n2=0,1
      in2=(iy+n2)*nclz1
      do n3=0,1
        in3=in1+in2+iz+n3+1
        st(n1,n2,n3)=phix(in3) ! save 8 neighbours static field values
      enddo
    enddo
  enddo
  ! linearly extrapolate sf values for x
  do n2=0,1
    do n3=0,1
      st2(n2,n3)=(st(1,n2,n3)-st(0,n2,n3))*dx+st(0,n2,n3)
    enddo
  enddo
  ! linearly extrapolate sf values for y
  do n3=0,1
    st3(n3)=(st2(1,n3)-st2(0,n3))*dy+st2(0,n3)
  enddo
  ! linearly extrapolate sf value for z 
  stat=(st3(1)-st3(0))*dz+st3(0)
endif
return
end subroutine

subroutine staticener(xj,yj,zj,ener)
!-----------------------------------------------------------------------
!This subroutine computes only the external static field energy
!for one particle used in subroutine INTERACT in simul.f   

use constamod
use gsbpmod
implicit none
!input variables
integer j, jtype
real  xj, yj, zj
!local variables
integer ncyz,ix,iy,iz,n1,n2,n3,in3
real  xi,yi,zi,ai,bi,ci,fi,ener
logical*1 ok

ener = 0.0

ok=xj.le.xbcen1+tranx1.and.xj.ge.xbcen1-tranx1.and. &
   yj.le.ybcen1+trany1.and.yj.ge.ybcen1-trany1.and. &
   zj.le.zbcen1+tranz1.and.zj.ge.zbcen1-tranz1
if (ok) then
  ncyz = ncly1*nclz1
!  ion cartesian coordinates in the local grid system              
  xi = xj + tranx1-xbcen1
  yi = yj + trany1-ybcen1
  zi = zj + tranz1-zbcen1
!  integer*4 counter for ion cartesian coordinates        
  ix = int(xi*idcel1)
  iy = int(yi*idcel1)
  iz = int(zi*idcel1)
  if (ix.eq.nclx1-1) ix=nclx1-2
  if (iy.eq.ncly1-1) iy=ncly1-2
  if (iz.eq.nclz1-1) iz=nclz1-2

!Atom charge distribution by 8 adjacent grid points

  do n1 = ix, ix+1
    ai = xi - n1*dcel1
    ai = 1.0 - abs(ai)*idcel1
    do n2 = iy, iy+1
      bi = yi - n2*dcel1
      bi = 1.0 - abs(bi)*idcel1
      do n3 = iz, iz+1
        ci = zi - n3*dcel1
        ci = 1.0 - abs(ci)*idcel1
        fi = ai*bi*ci
        in3 = n1*ncyz + n2*nclz1 + n3 + 1
!Electrostatic Energy 
        ener = ener + fi*phix(in3)*celec
      enddo ! n3
    enddo ! n2
  enddo ! n1
endif ! ok

return
end subroutine

subroutine statxd(ix1,iy1,iz1,ix2,iy2,iz2,unitn,Qnohead)
!-----------------------------------------------------------------------
!This subroutine prints the electrostatic potential resulting from the static field in the format x,y,z,pot 
use gsbpmod

!local variables
implicit none
integer ncyz,i,j,k
real  xii,yii,zii
integer ix1,iy1,iz1,ix2,iy2,iz2
real res,pos,dist,idist
integer unitn
logical*1 Qnohead
character line*1024

ncyz = ncly1*nclz1
if (ix1.lt.0) ix1=0
if (iy1.lt.0) iy1=0
if (iy1.lt.0) iz1=0
if (ix2.ge.nclx1) ix2=nclx1-1
if (iy2.ge.ncly1) iy2=ncly1-1
if (iz2.ge.nclz1) iz2=nclz1-1

if (.not.Qnohead) then
  write(unitn,*) 'Writing Static Field ASCII File ...'
  write(unitn,*) 'X Range: ',ix1,' - ',ix2
  write(unitn,*) 'Y Range: ',iy1,' - ',iy2
  write(unitn,*) 'Z Range: ',iz1,' - ',iz2
  write(unitn,'(a)') '   x        y        z        pot'
endif

do i=ix1,ix2
  do j=iy1,iy2    
    do k=iz1,iz2
      xii=(i-1)*dcel1-tranx1+xbcen1
      yii=(j-1)*dcel1-trany1+ybcen1
      zii=(k-1)*dcel1-tranz1+zbcen1
      write(line,*) xii,yii,zii,phix(i*ncyz+j*nclz1+k+1)
      write(unitn,'(a)') trim(line)
    enddo
  enddo
enddo

return
end subroutine

subroutine repxd(ix1,iy1,iz1,ix2,iy2,iz2,unitn,Qnohead)
!-----------------------------------------------------------------------
!This subroutine prints the electrostatic potential resulting from the repulsion field in the format x,y,z,rep
use gsbpmod

!local variables
implicit none
integer ncyz,i,j,k
real  xii,yii,zii
integer ix1,iy1,iz1,ix2,iy2,iz2
real res,pos,dist,idist
integer unitn
logical*1 Qnohead
character line*1024

ncyz = ncly2*nclz2
if (ix1.lt.0) ix1=0
if (iy1.lt.0) iy1=0
if (iy1.lt.0) iz1=0
if (ix2.ge.nclx2) ix2=nclx2-1
if (iy2.ge.ncly2) iy2=ncly2-1
if (iz2.ge.nclz2) iz2=nclz2-1

if (.not.Qnohead) then
  write(unitn,*) 'Writing Repulsion Field ASCII File ...'
  write(unitn,*) 'X Range: ',ix1,' - ',ix2
  write(unitn,*) 'Y Range: ',iy1,' - ',iy2
  write(unitn,*) 'Z Range: ',iz1,' - ',iz2
  write(unitn,'(a)') '   x        y        z        pot'
endif

do i=ix1,ix2
  do j=iy1,iy2
    do k=iz1,iz2
      xii=(i-1)*dcel2-tranx2+xbcen2
      yii=(j-1)*dcel2-trany2+ybcen2
      zii=(k-1)*dcel2-tranz2+zbcen2
      write(line,*) xii,yii,zii,phiv(i*ncyz+j*nclz2+k+1)
      write(unitn,'(a)') trim(line)
    enddo
  enddo
enddo

return
end subroutine

