!    COMPESSCOOR - Computes Essential Coordinates for DNA.
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

module xyz
implicit none
integer*4 nop
real*8,allocatable :: x(:,:)
character*4,allocatable :: pn(:)

end module

program compesscoor
use mathmod
use xyz
implicit none
integer*4 arg,narg,prin,n
real*8,allocatable :: cx(:,:),cxa(:,:),xx(:,:)
real*8,external :: funk
real*8 ad,aa,asd,asa,ave(3,6),cave(3,6),g,praxis,t0,h0,dx(20),pr,rot(3,3) !,cent(3),ccent(3)
character line*256

n=20
t0=0.0001
h0=10.0
prin=1

call header()
arg=0
narg=COMMAND_ARGUMENT_COUNT()

write(*,'(/A)') 'Warning: XYZ DNA must have the specified format: '
write(*,'(A)') '         First strand order from -z to z from 5p to 3p direction' 
write(*,'(A)') '         Second strand order from z to -z in 5p to 3p direction' 
write(*,'(A/)') '         Phos (5p) Base Sugar sucesive order' 

call readarg('XYZ Filename: ',narg,arg,line) ! read filename from argument or ask
call readxyz(trim(line)) ! read xyz
call rmgcent(nop,x) ! move to geom center to 0,0,0
call aligntoz(nop,x) ! align to z axis
call writexyz2('struct.xyz',nop,pn,x) 
allocate (cx(3,nop),cxa(3,nop),xx(3,nop)) ! allocate cylindrical coordinates vector
xx=x
call rotsecstrnd(nop,xx) ! rotates second strand
call cart2cyl(nop,xx,cx) ! calculate cylindrical coordinates
!call writexyz2('cylin.xyz',nop,pn,cx) ! debug
!call cyl2cart(nop,cx,xx) !  debug
!call writexyz2('cart.xyz',nop,pn,xx) ! debug
cx(2,:)=cx(2,:)+180d0 ! convert from -180/180 to 0/360
call averdistnuc(nop,cx,ad,asd)  ! compute average distance inter-nucleotides
call averrotnuc(nop,cx,aa,asa)  ! compute average rotation inter-nucleotides
call alignnucleot(nop,cx,cxa,ad,aa,asd,asa) ! compute superimpose cyl coordinates for each type
cxa(2,:)=cxa(2,:)-180d0 ! convert from -180/180 to 0/360
call cyl2cart(nop,cxa,xx)
call avecartcoor(nop,xx,pn,ave)
call cart2cyl(6,ave,cave) ! convert to cylindrical
g=(sum(cave(3,3:6))/4d0+sum(cave(3,1:2)))/3d0 ! get geometric center
cave(3,1:6)=cave(3,1:6)-g ! centering 
!cent=sum(ave,dim=2)/6d0
!call cart2cyl(1,cent,ccent)
cave(2,1:6)=cave(2,1:6)-asa/2d0+180d0

call rebuilddna(nop,pn,ave,ad,aa,asd,asa,xx)

call rotationmatrix(nop,x,xx,rot)
xx=matmul(transpose(rot),xx)

call writexyz2('rebuilt.xyz',nop,pn,xx) 

call cyl2cart(6,cave,ave)

write(*,*)
write(*,*)
write(*,*) 'P  ',ave(1:3,1),'     ',cave(1:3,1)
write(*,*) 'S  ',ave(1:3,2),'     ',cave(1:3,2)
write(*,*) 'Ab ',ave(1:3,3),'     ',cave(1:3,3)
write(*,*) 'Cb ',ave(1:3,4),'     ',cave(1:3,4)
write(*,*) 'Gb ',ave(1:3,5),'     ',cave(1:3,5)
write(*,*) 'Tb ',ave(1:3,6),'     ',cave(1:3,6)
write(*,*)
write(*,*) 'Internuc dist = ',ad,'  Internuc ang = ',aa
write(*,*)
write(*,*)

call rebuilddna2(nop,pn,cave,ad,aa,xx)

call rotationmatrix(nop,x,xx,rot)
xx=matmul(transpose(rot),xx)

call writexyz2('rebuilt2.xyz',nop,pn,xx) 

call cyl2cart(6,cave,ave)
dx(1:3)=ave(1:3,1)
dx(4:6)=ave(1:3,2)
dx(7:9)=ave(1:3,3)
dx(10:12)=ave(1:3,4)
dx(13:15)=ave(1:3,5)
dx(16:18)=ave(1:3,6)
dx(19)=ad
dx(20)=aa


write(*,*) 'Alignment function Initial Value =', dsqrt(funk(dx,n)/dfloat(nop))

pr = praxis(t0,h0,n,prin,dx,funk)

write(*,*) 'Aligment function Final Value =', dsqrt(funk(dx,n)/dfloat(nop))

ave(1:3,1) = dx(1:3)
ave(1:3,2) = dx(4:6)
ave(1:3,3) = dx(7:9)
ave(1:3,4) = dx(10:12)
ave(1:3,5) = dx(13:15)
ave(1:3,6) = dx(16:18)
ad=dx(19)
aa=dx(20)

call cart2cyl(6,ave,cave)
write(*,*) 'Readjusted parameters'
write(*,*)
write(*,*) 'P  ',ave(1:3,1),'     ',cave(1:3,1)
write(*,*) 'S  ',ave(1:3,2),'     ',cave(1:3,2)
write(*,*) 'Ab ',ave(1:3,3),'     ',cave(1:3,3)
write(*,*) 'Cb ',ave(1:3,4),'     ',cave(1:3,4)
write(*,*) 'Gb ',ave(1:3,5),'     ',cave(1:3,5)
write(*,*) 'Tb ',ave(1:3,6),'     ',cave(1:3,6)
write(*,*)
write(*,*) 'Internuc dist = ',ad,'  Internuc ang = ',aa
write(*,*)
write(*,*)

call rebuilddna2(nop,pn,cave,ad,aa,xx)

call rotationmatrix(nop,x,xx,rot)
xx=matmul(transpose(rot),xx)

call writexyz2('rebuilt3.xyz',nop,pn,xx)

! de pablo DNA
!ave(1:3,1)=(/-0.628,8.896,2.186/)
!ave(1:3,2)=(/ 2.365,6.568,1.280/)
!ave(1:3,3)=(/ 0.575,0.516,0.051/)
!ave(1:3,4)=(/ 0.199,2.287,0.187/)
!ave(1:3,5)=(/ 0.628,0.540,0.053/)
!ave(1:3,6)=(/ 0.159,2.344,0.191/)
!asd=0d0
!aa=36d0
!ad=3.38d0
!asa=0d0
!call rebuilddna(nop,pn,ave,ad,aa,asd,asa,xx)
!
!
!call writexyz2('depablo.xyz',nop,pn,xx)

write(*,'(/A/)') 'Normal Termination.'

end program

function funk(dx,n)
use xyz
use mathmod
implicit none
integer*4 n
real*8 dx(n),funk,ad,aa,ave(3,6),xx(3,nop),cave(3,6),rot(3,3)

ave(1:3,1) = dx(1:3)  
ave(1:3,2) = dx(4:6)  
ave(1:3,3) = dx(7:9)  
ave(1:3,4) = dx(10:12)
ave(1:3,5) = dx(13:15)
ave(1:3,6) = dx(16:18)
ad=dx(19)
aa=dx(20)
call cart2cyl(6,ave,cave)
call rebuilddna2(nop,pn,cave,ad,aa,xx)
call rotationmatrix(nop,x,xx,rot)
funk=sum((matmul(transpose(rot),xx)-x)**2)
end function

subroutine rebuilddna(nop,pn,ave,ad,aa,asd,asa,xx)
implicit none
integer*4 nop,i,j
real*8 ave(3,6),ad,aa,asd,asa,xx(3,nop),avcl(3,3),avcr(3,3),cx(3,nop),cenz
character*4 pn(nop)

! PBS

! build dna
j=0
do i=1,nop,3
  if (i.eq.nop/2+1) j=0
  avcr(1:3,1)=ave(1:3,1)
  if (pn(i+1).eq.'Ab') then
    avcr(1:3,2)=ave(1:3,3)
  elseif (pn(i+1).eq.'Cb') then
    avcr(1:3,2)=ave(1:3,4)
  elseif (pn(i+1).eq.'Gb') then
    avcr(1:3,2)=ave(1:3,5)
  elseif (pn(i+1).eq.'Tb') then
    avcr(1:3,2)=ave(1:3,6)
  endif
  avcr(1:3,3)=ave(1:3,2)
  call cart2cyl(3,avcr,avcl)
  xx(3,i:i+2)=avcl(3,1:3)+dfloat(j)*ad
  xx(2,i:i+2)=avcl(2,1:3)+dfloat(j)*aa
  xx(1,i:i+2)=avcl(1,1:3)
  j=j+1
enddo

! rotates 2nd strand
xx(2,nop/2+1:nop)=xx(2,nop/2+1:nop)-asa

! correct angles
!call correctang(nop,xx(2,1:nop)) ! not needed

cx=xx
! convert cylinder to cartesian
call cyl2cart(nop,cx,xx)

! center
! center 1st strand
cenz=2d0*sum(xx(3,1:nop/2))/dfloat(nop)
xx(3,1:nop/2)=xx(3,1:nop/2)-cenz
! center 2nd strand
cenz=2d0*sum(xx(3,nop/2+1:nop))/dfloat(nop)
xx(3,1:nop/2)=xx(3,nop/2+1:nop)-cenz

! translate 2nd strand
xx(3,nop/2+1:nop)=xx(3,nop/2+1:nop)-asd

!invert second strand
xx(2:3,nop/2+1:nop)=-xx(2:3,nop/2+1:nop)

!center whole dna
call rmgcent(nop,xx)

end subroutine

subroutine rebuilddna2(nop,pn,cave,ad,aa,xx)
implicit none
integer*4 nop,i,j,sgni
real*8 cave(3,6),ad,aa,xx(3,nop),cx(3,nop),avcl(3,3),sgn
character*4 pn(nop)

! build dna
j=0
sgn=1d0
sgni=1
do i=1,nop,3
  if (i.eq.nop/2+1) then
    j=j-1
    sgn=-1d0
    sgni=-1
  endif
  avcl(1:3,1)=cave(1:3,1)
  if (pn(i+1).eq.'Ab') then
    avcl(1:3,2)=cave(1:3,3)
  elseif (pn(i+1).eq.'Cb') then
    avcl(1:3,2)=cave(1:3,4)
  elseif (pn(i+1).eq.'Gb') then
    avcl(1:3,2)=cave(1:3,5)
  elseif (pn(i+1).eq.'Tb') then
    avcl(1:3,2)=cave(1:3,6)
  endif
  avcl(1:3,3)=cave(1:3,2)
  cx(3,i:i+2)=sgn*avcl(3,1:3)+dfloat(j)*ad
  cx(2,i:i+2)=sgn*avcl(2,1:3)+dfloat(j)*aa
  cx(1,i:i+2)=avcl(1,1:3)
  j=j+sgni
enddo

! convert cyl to cart
call cyl2cart(nop,cx,xx)

! center both
call rmgcent(nop,xx)
end subroutine

subroutine avecartcoor(nop,x,pn,ave)
implicit none
integer*4 nop,i,ele(6)
real*8 x(3,nop),ave(3,6)
character*4 pn(nop)

! 1 P / 2 S / 3 Ab / 4 Cb / 5 Gb / 6 Tb
ele=0
ele(1)=nop/3
ele(2)=nop/3

ave=0d0
do i=1,nop,3
  ave(1:3,1)=ave(1:3,1)+x(1:3,i)
  if (pn(i+1).eq.'Ab') then
    ele(3)=ele(3)+1
    ave(1:3,3)=ave(1:3,3)+x(1:3,i+1)
  elseif (pn(i+1).eq.'Cb') then
    ele(4)=ele(4)+1
    ave(1:3,4)=ave(1:3,4)+x(1:3,i+1)
  elseif (pn(i+1).eq.'Gb') then
    ele(5)=ele(5)+1
    ave(1:3,5)=ave(1:3,5)+x(1:3,i+1)
  elseif (pn(i+1).eq.'Tb') then
    ele(6)=ele(6)+1
    ave(1:3,6)=ave(1:3,6)+x(1:3,i+1)
  endif
  ave(1:3,2)=ave(1:3,2)+x(1:3,i+2)
enddo

do i=1,6
  ave(1:3,i)=ave(1:3,i)/dfloat(ele(i))
enddo
end subroutine

subroutine cyl2cart(nop,cx,x)
implicit none
real*8,parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679d0
real*8,parameter :: pi180inv = pi/180d0
integer*4 nop,i
real*8 cx(3,nop),x(3,nop)

do i=1,nop
  x(1,i)=cx(1,i)*dcos(cx(2,i)*pi180inv)
  x(2,i)=cx(1,i)*dsin(cx(2,i)*pi180inv)
  x(3,i)=cx(3,i)
enddo

end subroutine

subroutine alignnucleot(nop,cx,cxa,ad,aa,asd,asa)
implicit none
integer*4 nop,i,j
real*8 ad,aa,asd,asa,cx(3,nop),cxa(3,nop)

! copy
cxa=cx
do i=nop/2+1,nop
  ! translate second strand
  cxa(3,i)=cxa(3,i)+asd
  ! rotates second strand
  cxa(2,i)=cxa(2,i)+asa
enddo
! superimpose z position
  ! first strand
j=0
do i=1,nop/2,3
  cxa(3,i)=cxa(3,i)-dfloat(j)*ad
  cxa(3,i+1)=cxa(3,i+1)-dfloat(j)*ad
  cxa(3,i+2)=cxa(3,i+2)-dfloat(j)*ad
  j=j+1
enddo
  ! second strand
j=0
do i=nop/2+1,nop,3
  cxa(3,i)=cxa(3,i)-dfloat(j)*ad
  cxa(3,i+1)=cxa(3,i+1)-dfloat(j)*ad
  cxa(3,i+2)=cxa(3,i+2)-dfloat(j)*ad
  j=j+1
enddo
! superimpose angular position
  !first strand
j=0
do i=1,nop/2,3
  cxa(2,i)=cxa(2,i)-dfloat(j)*aa
  cxa(2,i+1)=cxa(2,i+1)-dfloat(j)*aa
  cxa(2,i+2)=cxa(2,i+2)-dfloat(j)*aa
  j=j+1
enddo
  !second strand
j=0
do i=nop/2+1,nop,3
  cxa(2,i)=cxa(2,i)-dfloat(j)*aa
  cxa(2,i+1)=cxa(2,i+1)-dfloat(j)*aa
  cxa(2,i+2)=cxa(2,i+2)-dfloat(j)*aa
  j=j+1
enddo
! correct angles
!call correctang(nop,cxa(2,1:nop))

end subroutine

subroutine correctang(n,ang)
use mathmod
implicit none
integer*4 n,i
real*8 ang(n)

do i=1,n
  ang(i)=ang(i)-360d0*inint(ang(i)/360d0)
enddo
end subroutine

! rotates second strand
subroutine rotsecstrnd(nop,x)
implicit none
integer*4 nop
real*8 x(3,nop)

x(2:3,nop/2+1:nop)=-x(2:3,nop/2+1:nop)

write(*,*) '2nd Strand has been rotated 180 degress  along y-z plane'
write(*,*)

end subroutine

! compute average rot inter-nucleotides
subroutine averrotnuc(nop,cx,aa,asa)
implicit none
integer*4 nop,i,j,nnop
real*8 cx(3,nop),aa,ad1,ad2,angcnt,ang(nop/3),angprev,sum1,sum2,asa


! first strand
j=0
do i=1,nop/2,3
  j=j+1
  ang(j)=cx(2,i)
  j=j+1
  ang(j)=cx(2,i+2)
enddo
nnop=j

angprev=ang(1)
angcnt=0d0
do i=2,nnop
  if (angprev.gt.180d0.and.ang(i).le.180d0) angcnt=angcnt+360d0
  angprev=ang(i)
  ang(i)=ang(i)+angcnt
enddo

sum1=sum(ang(1:nnop))

ad1=(ang(nnop)+ang(nnop-1)-ang(1)-ang(2))

! second strand
j=0
do i=nop/2+1,nop,3
  j=j+1
  ang(j)=cx(2,i)
  j=j+1
  ang(j)=cx(2,i+2)
enddo
nnop=j

angprev=ang(1)
angcnt=0d0
do i=2,nnop
  if (angprev.gt.180d0.and.ang(i).le.180d0) angcnt=angcnt+360d0
  angprev=ang(i)
  ang(i)=ang(i)+angcnt
enddo

sum2=sum(ang(1:nnop))

ad2=(ang(nnop)+ang(nnop-1)-ang(1)-ang(2))

aa=(ad1+ad2)/dfloat(2*nop/3-4)

write(*,*) 'Average angle internucleotides = ',aa
write(*,*) 
asa=(sum1-sum2)/dfloat(nnop)
write(*,*) 'Average angle strand separation (1-2)= ',asa
write(*,*) 
end subroutine

! compute average distance inter-nucleotides
subroutine averdistnuc(nop,cx,ad,asd)
implicit none
integer*4 nop
!integer*4 nop,i,hnop
real*8 cx(3,nop),ad,asd
!real*8 ad1,ad2,adn,zp,inop

!hnop=nop/2
!inop=2*(hnop/3-1)
!ad1=0d0
!zp=0d0
!do i=1,hnop,3
!  adn=cx(3,i)+cx(3,i+2) !+cx(3,i+1)
!  if (i.gt.3) ad1=ad1+dabs(adn-zp)
!  zp=adn
!enddo
!ad1=ad1/inop
!ad1=(cx(3,hnop)+cx(3,hnop-2)-cx(3,1)-cx(3,3))/inop  ! this replaces the prev commented algorithm

!ad2=0d0
!zp=0d0
!do i=hnop+1,nop,3
!  adn=cx(3,i)+cx(3,i+2) !+cx(3,i+2)
!  if (i.gt.hnop+3) ad2=ad2+dabs(zp-adn)
!  zp=adn
!enddo
!ad2=ad2/inop
!write(*,*) ad2
!write(*,*) (-cx(3,nop)-cx(3,nop-2)+cx(3,hnop+1)+cx(3,hnop+3))/inop
!
!ad=(ad1+ad2)/2d0
ad=(cx(3,nop)+cx(3,nop-2)-cx(3,nop/2+1)-cx(3,nop/2+3)+cx(3,nop/2)+cx(3,nop/2-2)-cx(3,1)-cx(3,3))/dfloat(2*nop/3-4)

write(*,*) 'Average distance internucleotides = ',ad
write(*,*)
asd=(sum(cx(3,1:nop/2))-sum(cx(3,nop/2+1:nop)))/dfloat(nop/2)
write(*,*) 'Displacement between strands (1-2) = ',asd
write(*,*)

end subroutine

! calculate cylindrical coordinates
subroutine cart2cyl(nop,x,cx)
implicit none
integer*4 nop,i
real*8,parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679d0
real*8,parameter :: pi180 = 180d0/pi
real*8 x(3,nop),cx(3,nop)
do i=1,nop
  cx(1,i)=dsqrt(dot_product(x(1:2,i),x(1:2,i))) ! comp radio
enddo
cx(2,:)=datan2(x(2,:),x(1,:))*pi180 ! angle
cx(3,:)=x(3,:)  ! z 

end subroutine

! align to z axis
subroutine aligntoz(nop,x)
use mathmod
implicit none
integer*4 np,order(3),nop
real*8 bfv(3),a(3,3),s(3),v(3,3),rot(3,3),ref(3),rsum,rv1,rv2,x(3,nop) ! best fit vector
np=3
a=matmul(x,transpose(x))
call diag(a,np,s,v)
call msort(s,order,3)
s=s(order(1:3))
v=v(1:3,order(1:3))
bfv=v(1:3,1)
bfv=bfv/dsqrt(dot_product(bfv,bfv))
write(*,*) 
write(*,*) 'Alignment info:' 
write(*,*) 'Vector: ',bfv
if (dot_product(bfv,(/ 1d0, 0d0, 0d0 /)).gt.0.996d0) then ! if aligned to 1,0,0 use 0,1,0
  ref=(/ 0d0, 1d0, 0d0 /)
else
  ref=(/ 1d0, 0d0, 0d0 /)
endif
rot(2,:)=cross_product(bfv,ref) ! set y
rot(1,:)=cross_product((rot(2,:)),bfv) ! set x
rot(3,:)=bfv ! set z 
rsum=sum(matmul(transpose(x),x))
rv1=sum(matmul(transpose(x),(/0d0,0d0,1d0/)))
rv2=sum(matmul(transpose(x),bfv))
write(*,*) 'Alignment to z=',rsum-rv1
write(*,*) 'Alignment to vector=',rsum-rv2
! align to vector
x=matmul(rot,x)
rsum=sum(matmul(transpose(x),x))
rv1=sum(matmul(transpose(x),(/0d0,0d0,1d0/)))
write(*,*) 'After alignment, Alignment to z=',rsum-rv1

a=matmul(x,transpose(x))
call diag(a,np,s,v)
call msort(s,order,3)
s=s(order(1:3))
v=v(1:3,order(1:3))
bfv=v(1:3,1)
bfv=bfv/dsqrt(dot_product(bfv,bfv))
write(*,*) 'Vector after alignment: ',bfv
write(*,*) 

end subroutine

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='COMPESSCOOR'
prver='version 1.0'
prdesc='Computes Essential Coordinates for DNA.'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='22 Oct 2012'
lastdate='22 Oct 2012'

write(*,'(/A)') trim(prname)//' '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A/)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
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
  read(*,'(A)') text
endif
text=adjustl(trim(text))
end subroutine

! remove geometric center
subroutine rmgcent(nop,x)
implicit none
integer*4 nop
real*8 x(3,nop)
x=x-spread(sum(x,dim=2)/dfloat(nop),dim=2,ncopies=nop)
end subroutine

subroutine readxyz(filename)
use xyz
implicit none
character*(*) filename
integer*4 un,i

un=1

open(unit=un,file=filename) !open file 
read(un,*) nop ! read number of particles
read(un,*)  ! read comment line
if (.not.allocated(x)) allocate (x(3,nop))  ! allocate x
if (.not.allocated(pn)) allocate (pn(nop))  ! allocate pn

do i=1,nop
  read(un,*) pn(i),x(1,i),x(2,i),x(3,i)  ! read pn and x
enddo
close(un) ! close file
end subroutine

subroutine writexyz(filename)
use xyz
implicit none
character*(*) filename
integer*4 un,i

un=2

open(unit=un,file=filename) !open file 
write(un,*) nop ! write number of particles
write(un,*)  ! write comment line

do i=1,nop
  write(un,*) pn(i),x(1,i),x(2,i),x(3,i)  !write pn and x
enddo
close(un) ! close file
end subroutine

subroutine writexyz2(filename,nop,pn,x)
implicit none
integer*4 nop
real*8 x(3,nop)
character*4 pn(nop) 
character*(*) filename
integer*4 un,i

un=2

open(unit=un,file=filename) !open file 
write(un,*) nop ! write number of particles
write(un,*)  ! write comment line

do i=1,nop
  write(un,*) pn(i),x(1,i),x(2,i),x(3,i)  !write pn and x
enddo
close(un) ! close file
end subroutine

