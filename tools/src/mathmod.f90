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

module mathmod
implicit none
contains
  function inint(num)
  implicit none
  real*8 inint,num
  inint=dfloat(iint(num))
  end function

  function iint(num)
  implicit none
  integer*4 iint
  real*8 num
  if (num.ge.0.0d0) then
    iint=int(num)
  else
    if ((num-dfloat(int(num))).eq.0.0d0) then
      iint=int(num)
    else
      iint=int(num)-1
    endif
  endif
  end function

  function invmat(mat)
  implicit none
  real*8 invmat(3,3),mat(3,3)
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
  real*8 u(3),v(3),w(3),cross_product(3)
  w(1)=u(2)*v(3)-u(3)*v(2)
  w(2)=u(3)*v(1)-u(1)*v(3)
  w(3)=u(1)*v(2)-u(2)*v(1)
  cross_product=w
  end function
  
  subroutine svd(np,mp,a,u,s,v,is)
  implicit none
  integer*4 np,mp,i
  real*8 a(np,mp),u(np,mp),v(mp,mp),ss(mp),s(mp,mp),is(mp,mp)
  call diag(matmul(transpose(a),a),mp,ss,v)
  is=0d0
  s=0d0
  do i=1,mp
    is(i,i)=1d0/dsqrt(ss(i))
    s(i,i)=ss(i)
  enddo
  u=matmul(matmul(a,v),is)
  end subroutine

  subroutine svd2(np,mp,a,u,ss,v) !,is)
  implicit none
  integer*4 np,mp,i
  real*8 a(np,mp),u(np,mp),v(mp,mp),ss(mp),is(mp,mp)!,s(mp,mp)
  is=matmul(transpose(a),a)
  call diag(is,mp,ss,v)
  is=0d0
  !s=0d0
  do i=1,mp
    is(i,i)=1d0/dsqrt(ss(i))
  !  s(i,i)=ss(i)
  enddo
  u=matmul(matmul(a,v),is)
  end subroutine
  
  ! Solve the matricial equation system B=Ax
  subroutine bax(n,m,b,a,x)
  implicit none
  integer*4 n,m
  real*8 b(n),x(n),a(n,m),ainv(m,n)
  call pinv(n,m,a,ainv)
  x=matmul(b,transpose(ainv))
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

  subroutine Merge(A,NA,B,NB,C,NC,D,N)
  integer*4, intent(in) :: NA,NB,NC,N         ! Normal usage: NA+NB = NC
  integer*4, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
  integer*4, intent(in)     :: B(NB)
  integer*4, intent(in out) :: C(NC)
  real*8, intent(in)        :: D(N)
  integer*4 :: I,J,K
  I = 1; J = 1; K = 1;
  do while(I <= NA .and. J <= NB)
     if (D(A(I)) >= D(B(J))) then  ! exchange > by < to invert order
        C(K) = A(I)
        I = I+1
     else
        C(K) = B(J)
        J = J+1
     endif
     K = K + 1
  enddo
  do while (I <= NA)
     C(K) = A(I)
     I = I + 1
     K = K + 1
  enddo
  return
  end subroutine
  
  recursive subroutine MergeSort(A,M,R,N,T)
  implicit none
  integer*4, intent(in) :: N,M
  integer*4, dimension(N), intent(in out) :: R
  real*8, dimension(M), intent(in) :: A
  integer*4, dimension((N+1)/2), intent (out) :: T
  integer*4 :: NA,NB,V
  if (N < 2) return
  if (N == 2) then
     if (A(R(1)) < A(R(2))) then   ! exchange < by > to invert order
        V = R(1)
        R(1) = R(2)
        R(2) = V
     endif
     return
  endif
  NA=(N+1)/2
  NB=N-NA
  call MergeSort(A,M,R,NA,T)
  call MergeSort(A,M,R(NA+1),NB,T)
  if (A(R(NA)) < A(R(NA+1))) then  ! exchange < by > to invert order
     T(1:NA)=R(1:NA)
     call Merge(T,NA,R(NA+1),NB,R,N,A,M)
  endif
  return
  end subroutine
  
  subroutine msort(a,b,n)
  integer*4 n,i
  real*8, dimension(n) :: a
  integer*4, dimension(n) :: b
  integer*4, dimension ((n+1)/2) :: t
  do i=1,n
    b(i)=i
  enddo
  call MergeSort(a,n,b,n,t) ! order from higher to lower
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

  function det(m) ! determinant 
  implicit none
  real*8 det, m(3,3)
  det=m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2))+m(1,2)*(m(2,3)*m(3,1)-m(2,1)*m(3,3))+m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
  end function

  subroutine rotationmatrix(n,mat1,mat2,rot)
  implicit none
  integer*4 n
  real*8 mat1(3,n),mat2(3,n),rot(3,3),a(3,3),s(3),u(3,3),v(3,3),one
  if (n.lt.4) stop 'number of particles cannot be lower than 4'
  a=matmul(mat1,transpose(mat2))
  one=sign(1d0,det(a))
  call svd2(3,3,a,u,s,v)
  v(1:3,3)=v(1:3,3)*one
  rot=matmul(v,transpose(u))
  end subroutine

end module
