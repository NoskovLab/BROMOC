!    MATH Library
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

!program test
!implicit none
!real*8 a(5,4),b(5),x(4) 
!a(1,1)=1.0
!a(2,1)=2.0
!a(3,1)=3.0
!a(4,1)=1.5
!a(5,1)=4.0
!a(1,2)=-1.0
!a(2,2)=2.0
!a(3,2)=-3.0
!a(4,2)=1.5
!a(5,2)=4.0
!a(1,3)=-1.0
!a(2,3)=-2.0
!a(3,3)=-3.0
!a(4,3)=1.5
!a(5,3)=4.0
!a(1,4)=1.0
!a(2,4)=-2.0
!a(3,4)=3.0
!a(4,4)=-1.5
!a(5,4)=-4.0
!x(1)=2.0
!x(2)=0.25
!x(3)=-2.0
!x(4)=-0.25
!!   2.0000000000000004       0.25000000000000189        2.0000000000000013      -0.24999999999999697     
!write(*,*) x
!b=matmul(a,x)
!write(*,*) b
!x=0d0
!call bax(5,4,b,a,x)
!write(*,*) x
!
!end program

! Solve the matricial equation system B=Ax
subroutine bax(n,m,b,a,x)
implicit none
integer*4 n,m
real*8 b(n),x(m),a(n,m),ainv(m,n)
call pinv(n,m,a,ainv)
x=matmul(ainv,b)
end subroutine

! matrix diagonalization 
subroutine diag(a,np,s,v)
implicit none
integer*4 np
real*8 a(np,np),s(np),v(np,np),e(np)
v=a
s=0d0
call tred2(v,np,s,e)
call tqli(s,e,np,v)
end subroutine

subroutine pinv(np,mp,mat,imat)
implicit none
integer*4 i,np,mp
real*8 mat(np,mp),imat(mp,np),ss(mp),vv(mp,mp),s(mp,mp)
real*8,parameter :: low=1d-8
call diag(matmul(transpose(mat),mat),mp,ss,vv)         ! A = V*S2*Vt
!  write(*,*) ss
s=0d0
do i=1,mp
  if (abs(ss(i)).le.low) then
    s(i,i)=0d0
  else
    s(i,i)=1d0/ss(i)     ! S^2 -> S^-2
  endif
enddo
! U = R*V*S-1
! Ut= S-1*Vt*Rt
! R-1 = V*S-1*Ut
! R-1 = V*S-2*Vt*Et
imat=matmul(matmul(matmul(vv,s),transpose(vv)),transpose(mat))
end subroutine

subroutine tred2(a,np,d,e)
implicit none
integer*4 np,n,i,l,k,j
real*8 a(np,np),d(np),e(np),h,f,g,hh,scale
n=np
if(n.gt.1)then
  do 18 i=n,2,-1  
    l=i-1
    h=0d0
    scale=0d0
    if(l.gt.1)then
      do 11 k=1,l
        scale=scale+abs(a(i,k))
11    continue
      if(scale.eq.0d0)then
        e(i)=a(i,l)
      else
        do 12 k=1,l
          a(i,k)=a(i,k)/scale
          h=h+a(i,k)**2
12      continue
        f=a(i,l)
        g=-sign(dsqrt(h),f)
        e(i)=scale*g
        h=h-f*g
        a(i,l)=f-g
        f=0d0
        do 15 j=1,l
          a(j,i)=a(i,j)/h
          g=0d0
          do 13 k=1,j
            g=g+a(j,k)*a(i,k)
13        continue
          if(l.gt.j)then
            do 14 k=j+1,l
              g=g+a(k,j)*a(i,k)
14          continue
          endif
          e(j)=g/h
          f=f+e(j)*a(i,j)
15      continue
        hh=f/(h+h)
        do 17 j=1,l
          f=a(i,j)
          g=e(j)-hh*f
          e(j)=g
          do 16 k=1,j
            a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16        continue
17      continue
      endif
    else
      e(i)=a(i,l)
    endif
    d(i)=h
18 continue
endif
d(1)=0d0
e(1)=0d0
do 23 i=1,n
  l=i-1
  if(d(i).ne.0d0)then
    do 21 j=1,l
      g=0d0
      do 19 k=1,l
        g=g+a(i,k)*a(k,j)
19    continue
      do 20 k=1,l
        a(k,j)=a(k,j)-g*a(k,i)
20    continue
21  continue
  endif
  d(i)=a(i,i)
  a(i,i)=1d0
  if(l.ge.1)then
    do 22 j=1,l
      a(i,j)=0d0
      a(j,i)=0d0
22  continue
  endif
23 continue
return
end subroutine

subroutine tqli(d,e,np,z)
implicit none
integer*4 np,n,i,l,iter,m,k
real*8 d(np),e(np),z(np,np),dd,g,r,s,c,p,f,b
n=np
if (n.gt.1) then
  do 11 i=2,n
    e(i-1)=e(i)
11  continue
  e(n)=0d0
  do 15 l=1,n
    iter=0
1    do 12 m=l,n-1
      dd=abs(d(m))+abs(d(m+1))
      if (abs(e(m))+dd.eq.dd) go to 2
12    continue
    m=n
2    if(m.ne.l)then
      if(iter.eq.30) stop 'too many iterations'
      iter=iter+1
      g=(d(l+1)-d(l))/(2d0*e(l))
      r=dsqrt(g**2+1d0)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1d0
      c=1d0
      p=0d0
      do 14 i=m-1,l,-1
        f=s*e(i)
        b=c*e(i)
        if(abs(f).ge.abs(g))then
          c=g/f
          r=dsqrt(c**2+1d0)
          e(i+1)=f*r
          s=1d0/r
          c=c*s
        else
          s=f/g
          r=dsqrt(s**2+1d0)
          e(i+1)=g*r
          c=1d0/r  
          s=s*c
        endif
        g=d(i+1)-p
        r=(d(i)-g)*s+2d0*c*b
        p=s*r
        d(i+1)=g+p
        g=c*r-b
        do 13 k=1,n
          f=z(k,i+1)
          z(k,i+1)=s*z(k,i)+c*f
          z(k,i)=c*z(k,i)-s*f
13        continue
14      continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=0d0
      go to 1
    endif
15  continue
endif
return
end subroutine

