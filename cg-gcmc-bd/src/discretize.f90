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

subroutine discretize(is,nnp,mnp,nxf)
! output variables: ep, nxi, dmi ,dm2  (all shared)
use efpmod
use constamod
use grandmod
use extramod

implicit none
integer k,is,nn,mnp,nnp,nxf
real xx(1-nnp+mnp),yy(1-nnp+mnp)
real scal,scald,x1,x2,d1,d2,dd2,y0,kf
real cc0,cc1,cc2,cc3,cc4,dist12,dist6,fdf,fdv

if (Qchr(is)) then
  kf=fct(is)
else
  kf=0.0
endif
do k=nnp,mnp
  x1=k*res
  x2=1.0/x1
  scal=0.0
  scald=0.0
!  scaldd=0.0
  scal=kf*x2     !coulomb potential
  scald=-kf*x2**2
!  if (k.eq.mnp) then
!    y1=2.0*kf*x2**3
!    scaldd=scal
!  endif
  if (Qlj(is)) then
    dist6=sgp2(is)**3*x2**6
    dist12=dist6**2
    scal=scal+epp4(is)*(dist12-dist6) ! van der waals potential
    scald=scald+epp4(is)*(dist6-2.0*dist12)*6.0*x2! van der waals force
!!  scaldd=scaldd+epp4(is)*(156.0*dist12-42.0*dist6)*x2**2
  endif
  if (Qsrpmfi(is)) then
    if (x1**2.le.rth) then
      cc1   = exp((c1(is)-x1)*c2(is))
      cc2   = cos(c3(is)*pi*(c1(is)-x1))
      cc3   = (c1(is)*x2)**6
      cc0   = c0(is)*cc1*cc2+c4(is)*cc3
      if (x1.ge.srpx) then ! smoothly fix discontinuity 
        fdf=exp(-srpk*(x1-srpx))-srpy
        fdv=cc0
        cc0=fdv*fdf
      endif
      scal  = scal + cc0
      cc4   = (-c2(is)*cc2+c3(is)*pi*sin(c3(is)*pi*(c1(is)-x1)))*c0(is)*cc1-6.0*c4(is)*cc3*x2! forces
      if (x1.ge.srpx) cc4=cc4*fdf-fdv*srpk*(fdf+srpy)  ! smoothly fix discontinuity 
      scald = scald + cc4
      !! cc0=cc1*cc2
      !! cc4=cc1*sin(c3(is)*pi*(c1(is)-x1))
      !! scaldd=scaldd+(c2(is)**2-c3(is)**2*pi**2)*c0(is)*cc0-2.0*c0(is)*c2(is)*c3(is)*pi*cc4+42.0*c4(is)*cc3*x2**2 !check derivative
    endif
  endif
  xx(k-nnp+1)=x1
  yy(k-nnp+1)=scal
!  yp(k-nnp+1)=-scald/x1 ! compute -(derivative of potential)/x = -force/x, force=derivative of potential
  if (k.eq.nnp) d1=scald
!  if (k.eq.mnp) d2=scald ! for a/r6 + b/r12 + k/r
enddo
nn=1-nnp+mnp
! for a/r6 + b/r12 + k/r
!sc(3,is)=xx(nn)**5*((d2*xx(nn)**2 - 11.0*kf)/6.0 + 2.0*yy(nn)*xx(nn))
!sc(4,is)=xx(nn)**11*((5.0*kf - d2*xx(nn)**2)/6.0 - yy(nn)*xx(nn))

! for a/r6 + k/r
sc(3,is)=(yy(nn)*xx(nn)-kf)*xx(nn)**5

! for a/r6
!sc(3,is)=yy(nn)*xx(nn)**6

! for a/r6 + b/r12 + k/r
!d2=-kf/xx(nn)**2-6.0*sc(3,is)/xx(nn)**7-12.0*sc(4,is)/xx(nn)**13
!dd2=2.0*kf/xx(nn)**3+42.0*sc(3,is)/xx(nn)**8+156.0*sc(4,is)/xx(nn)**14

! for a/r6 + k/r
d2=-kf/xx(nn)**2-6.0*sc(3,is)/xx(nn)**7
dd2=2.0*kf/xx(nn)**3+42.0*sc(3,is)/xx(nn)**8

! for a/r6
!d2=-6.0*sc(3,is)/xx(nn)**7
!dd2=42.0*sc(3,is)/xx(nn)**8


call squarespline(nn,xx,yy,ep(1:3,nxi(is):nxf),d1,d2,dd2)

sc(1,is)=-d1*xx(1)**13/12.0            ! compute a for y=a/x**12+b
y0=ep(1,nxi(is))+xx(1)*(ep(2,nxi(is))+xx(1)*ep(3,nxi(is)))
sc(2,is)=y0-sc(1,is)/xx(1)**12         ! compute b
!dd1=156.0*sc(1,is)/xx(1)**14           ! compute second derivative

!yp(nn)=-d2/xx(nn)
!yp(nn-1)=(yp(nn)+yp(nn-1))*0.5
!yp(nn-2)=(yp(nn-1)+yp(nn-2))*0.5
!
!dd1=-(yp(1)+dd1)/xx(1)   ! transform aceleration (derivative of force)
!dd2=-(yp(nn)+dd2)/xx(nn) !        to -(derivative of force)/x
!call squarespline(nn,xx,yp,ep(4:6,nxi(is):nxf),dd1,dd2,0.0)
end subroutine
