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

! Computes Splining Square Polynom Table
subroutine squarespline(n,x,y,p,d1,d2,dd2)
implicit none
integer n
real x(n),y(n),p(3,n),yp(n),pp(2,n),d1,d2,y0,yn,dd2

if (n.lt.1) then
  stop 'Error: n is lower than 1'
elseif (n.eq.1) then
  p(3,1)=0.0
  p(2,1)=0.0
  p(1,1)=y(1)
elseif (n.eq.2) then
  p(3,1)=0.0
  p(2,1)=(y(2)-y(1))/(x(2)-x(1))
  p(1,1)=y(1)-p(2,1)*x(1)
  p(3,2)=p(3,1)
  p(2,2)=p(2,1)
  p(1,2)=p(1,1)
else
  y0=y(1)
  yn=y(n)
  yp=y
  call sqfit(n,x,yp,p,d1,d2)  ! square interpolation averaged
  call fixder(n,x,p,yp)       ! derive and fix continuity of derivates and put it in yp
  call linipol(n,x,yp,pp,dd2) ! linear interpolation of derivatives yp to construct a lineal equation table (conserves second derivative at yn)
  call intlin(n,x,pp,p,y0)    ! integration of lineal equation table to get the square polynom
  call fixdy(n,x,p,yn)        ! fix yn-y0 preserving derivatives
endif
end subroutine 

! Integrates linear interpolation
subroutine intlin(n,x,pp,p,y0)
implicit none
integer n,i
real p(3,n),pp(2,n),x(n),y0

p(3,1)=pp(2,1)*0.5
p(2,1)=pp(1,1)
p(1,1)=y0-x(1)*(p(2,1)+p(3,1)*x(1))
do i=2,n
  p(3,i)=pp(2,i)*0.5
  p(2,i)=pp(1,i)
  p(1,i)=p(1,i-1)+x(i)*(p(2,i-1)-p(2,i)+x(i)*(p(3,i-1)-p(3,i)))
enddo
end subroutine

! Linear interpolation
subroutine linipol(n,x,y,p,dd)
implicit none
integer n,i
real x(n),y(n),p(2,n),tt,dd

do i=1,n-1
  tt=1.0/(x(i+1)-x(i))
  p(2,i)=(y(i+1)-y(i))*tt
  p(1,i)=(y(i)*x(i+1)-y(i+1)*x(i))*tt
enddo
!! Setting up p(2,n) so that preserves dy
!p(2,n)=dy
!do i=2,n
!  p(2,n)=p(2,n)-p(2,i-1)*0.5*(x(i)**2-x(i-1)**2)-p(1,i-1)*(x(i)-x(i-1))
!enddo
!p(2,n)=2.0*(p(2,n)+y(n)*x(n))/x(n)**2

! preserves second derivative at yn 
p(2,n)=dd
p(1,n)=y(n)-p(2,n)*x(n)
end subroutine

! Fix function so that yn is at the same place that originally and preserving derivatives
! by adding a corrective square equation 
subroutine fixdy(n,x,p,yn)
implicit none
integer n,i
real x(n),p(3,n),yn,y0,ynn,a,b,c,h
!Save first y
i=1
y0=p(1,i)+x(i)*(p(2,i)+x(i)*p(3,i))
!Save last y
i=n
ynn=p(1,i)+x(i)*(p(2,i)+x(i)*p(3,i))
! Compute correction preserving yn, y0 and first & last derivatives
b=2.0*(yn-ynn)/(x(n)-x(1))
c=-b*x(1)*0.5
h=(x(2)-x(1))/(n-1)
! Apply correction
do i=1,n-1
  p(1,i)=p(1,i)+c
  p(2,i)=p(2,i)+b
  a=-b*0.5/(x(i)+h*(i-1))
  p(3,i)=p(3,i)+a
enddo
i=n
p(1,i)=p(1,i)+c
p(2,i)=p(2,i)+b
a=-b*0.5/x(i)
p(3,i)=p(3,i)+a

end subroutine

! Derivate and Fix derivatives
subroutine fixder(n,x,p,yp)
implicit none
integer n,i
real yp(n),x(n),p(3,n)

yp(1)=p(2,1)+2.0*p(3,1)*x(1)
do i=2,n-1
  yp(i)=0.5*(p(2,i)+2.0*p(3,i)*x(i)+p(2,i-1)+2.0*p(3,i-1)*x(i))
enddo
i=n
yp(i)=p(2,i)+2.0*p(3,i)*x(i)
end subroutine

! Square interpolation averaged
subroutine sqfit(n,x,y,p,d1,d2)
implicit none
integer n,i
real x(n),y(n),p(3,n),tt,d1,d2

i=1
tt=1.0/(x(i+1)-x(i))                   !
p(3,i)=((y(i+1)-y(i))*tt-d1)*tt        !
p(2,i)=d1-2.0*p(3,i)*x(i)              ! Create first square equation conserving the first given derivative
p(1,i)=y(i)-x(i)*(p(2,i)+p(3,i)*x(i))  !

i=2
tt=(y(i+1)-y(i))/(x(i+1)-x(i))                             !
p(3,i)=(tt-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))    ! Compute next square equation
p(2,i)=tt-p(3,i)*(x(i+1)+x(i))                             !
p(1,i)=y(i)-x(i)*(p(2,i)+p(3,i)*x(i))                      !

do i=3,n-1
  tt=(y(i+1)-y(i))/(x(i+1)-x(i))                             !
  p(3,i)=(tt-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))    ! Compute next square equation
  p(2,i)=tt-p(3,i)*(x(i+1)+x(i))                             !
  p(1,i)=y(i)-x(i)*(p(2,i)+p(3,i)*x(i))                      !
  p(3,i-1)=0.5*(p(3,i-1)+p(3,i))          !
  p(2,i-1)=0.5*(p(2,i-1)+p(2,i))          ! Do average with the previous
  p(1,i-1)=0.5*(p(1,i-1)+p(1,i))          !
enddo

tt=1.0/(x(n)-x(n-1))                      !
p(3,n)=(d2-(y(n)-y(n-1))*tt)*tt           !
p(2,n)=d2-2.0*p(3,n)*x(n)                 !  Compute the last square equation conserving the last given derivative
p(1,n)=y(n)-x(n)*(p(2,n)+p(3,n)*x(n))     !
p(3,n-1)=0.5*(p(3,n-1)+p(3,n))       !
p(2,n-1)=0.5*(p(2,n-1)+p(2,n))       ! Do average with the previous
p(1,n-1)=0.5*(p(1,n-1)+p(1,n))       !
end subroutine

! Computes potential
subroutine gety(is,x2,y,x)
use efpmod
use extramod
implicit none
integer i,is
real x,x2,y,ix2
x=0.0
if (x2.lt.dm2(1,is)) then
  ix2=1.0/x2
  y=sc(1,is)*ix2**6 + sc(2,is)
elseif (x2.gt.dm2(2,is)) then
  ix2=1.0/x2
  if (Qchr(is)) then
    y=sc(3,is)*ix2**3+fct(is)*sqrt(ix2) ! coulomb term for tail
  else
    y=sc(3,is)*ix2**3 !+sc(4,is)*ix2**6
  endif
else
  x=sqrt(x2)
  i=int((x-dmi(is))*ires)+nxi(is)
  y=ep(1,i)+x*(ep(2,i)+x*ep(3,i))
endif
end subroutine

! Computes potential and forces/x
subroutine getyd(is,x2,y,yp,x)
use efpmod
use extramod
implicit none
integer i,is
real x2,x,y,yp,ix2,c,f

x=0.0
if (x2.lt.dm2(1,is)) then
  ix2=1.0/x2
  c=sc(1,is)*ix2**6
  y=c + sc(2,is)
  yp=12.0*c*ix2
elseif (x2.gt.dm2(2,is)) then
  ix2=1.0/x2
  c=sc(3,is)*ix2**3
!  d=sc(4,is)*ix2**6
  if (Qchr(is)) then ! coulomb term for tail
    f=fct(is)*sqrt(ix2)
    y=c + f !+ d
!    yp=(12.0*d + 6.0*c + f)*ix2
    yp=(6.0*c + f)*ix2
  else
    y=c !+ d
!    yp=(6.0*c + 12.0*d)*ix2 
    yp=6.0*c*ix2 
  endif  
else
  x=sqrt(x2)
  i=int((x-dmi(is))*ires)+nxi(is)
  y=ep(1,i)+x*(ep(2,i)+x*ep(3,i))
  yp=-ep(2,i)/x-2.0*ep(3,i)
endif
end subroutine

! Computes force/x
subroutine getyp(is,x2,yp)
use efpmod
use extramod
implicit none
integer i,is
real x,x2,yp,ix2

if (x2.lt.dm2(1,is)) then
  yp=12.0*sc(1,is)/x2**7 
elseif (x2.gt.dm2(2,is)) then
  ix2=1.0/x2
  if (Qchr(is)) then 
    yp=6.0*sc(3,is)*ix2**4+fct(is)*sqrt(ix2)*ix2    ! coulomb term for tail
  else
    yp=6.0*sc(3,is)*ix2**4 !+ 12.0*sc(4,is)*ix2**7
  endif
else
  x=sqrt(x2)
  i=int((x-dmi(is))*ires)+nxi(is)
  yp=-ep(2,i)/x-2.0*ep(3,i)
endif
end subroutine


