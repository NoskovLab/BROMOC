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

subroutine splinepot(is,nn,xx,yy,nxf)
use efpmod
use extramod
use grandmod
implicit none
integer nn,is,nxf
!integer,parameter :: nda=1
real xx(nn),yy(nn),d1,d2,dd1,dd2,kf!,a !,b,c,d,yda,xda,den

dmi(is)=xx(1)
dm2(1,is)=xx(1)**2
dm2(2,is)=xx(nn)**2
! Shifting potential
if(xx(nn).eq.0.0.or.xx(1).eq.0.0) stop 'Division by zero'
if (Qchr(is)) then
  kf=fct(is)
else
  kf=0.0
endif

! for k/x + a/x6 + b/x12
!a=0.0
!b=0.0
!c=0.0
!d=0.0
!do i=nn-nda+1,nn
!  yda=xx(i)**5*(xx(i)*yy(i)-kf)
!  xda=1.0/xx(i)**6
!  a=a+yda
!  b=b+xda
!  c=c+xda**2
!  d=d+xda*yda
!enddo
!den=1.0/(nda*c-b**2)
!sc(3,is)=(a*c-d*b)*den
!sc(4,is)=(d*nda-a*b)*den
!d2=-kf/xx(nn)**2-6.0*sc(3,is)/xx(nn)**7-12.0*sc(4,is)/xx(nn)**13
!dd2=2.0*kf/xx(nn)**3+42.0*sc(3,is)/xx(nn)**8+156.0*sc(4,is)/xx(nn)**14

! for k/x + a/x6
! for nda!=1
!a=0.0
!do i=nn-nda+1,nn
!  a = a + (xx(i)*yy(i)-kf)*xx(i)**5
!enddo
!sc(3,is)=a/nda
! for nda=1
sc(3,is)=(xx(nn)*yy(nn)-kf)*xx(nn)**5

d2=-kf/xx(nn)**2-6.0*sc(3,is)/xx(nn)**7
dd2=2.0*kf/xx(nn)**3+42.0*sc(3,is)/xx(nn)**8

!for nda!=1
!yy(nn)=kf/xx(nn)+sc(3,is)/xx(nn)**6

! for a/x6
!sc(3,is)=yy(nn)*xx(nn)**6
!d2=-6.0*sc(3,is)/xx(nn)**7
!dd2=42.0*sc(3,is)/xx(nn)**8

dd1=1.0/(xx(1)**12-xx(2)**12)
sc(1,is)=(yy(2)-yy(1))*xx(2)**12*xx(1)**12*dd1
sc(2,is)=(yy(1)*xx(1)**12-yy(2)*xx(2)**12)*dd1
d1=-12.0*sc(1,is)/xx(1)**13
dd1=156.0*sc(1,is)/xx(1)**14
call squarespline(nn,xx,yy,ep(1:3,nxi(is):nxf),d1,d2,dd2)

!ct2=0
!do l=nxi(is),nxf
!  ct2=ct2+1
!  yy(ct2)=-ep(2,l)/xx(ct2)-2.0*ep(3,l) ! compute -(derivative of potential)/x = -force/x, force=derivative of potential
!enddo

!dd1=-(yy(1)+dd1)/xx(1)   ! transform aceleration (derivative of force)
!dd2=-(yy(nn)+dd2)/xx(nn) !        to derivative of -force/x
!call squarespline(nn,xx,yy,ep(4:6,nxi(is):nxf),dd1,dd2,0.0)

end subroutine
