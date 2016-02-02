!    FIT-IT Fit mu (chemical potential) vs concentration and extrapolate mu for given concentration
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

program fitit
implicit none
integer narg,n,i,na,p,arg
character concc*100,filename*256,pname*256
real*8 conc,mu,x(1000),y1(1000),y2(1000) !,conc2,mu2
real*8,allocatable :: aa1(:),aa2(:),bb1(:),bb2(:)

! printout header
call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()

call readarg(' Data filename: ',narg,arg,filename)

call readarg(' Number of parameters: ',narg,arg,concc)
read(concc,*) na
write(*,*)

allocate (aa1(na),aa2(na)) !,bb1(na),bb2(na))

call readdata(trim(filename),n,x,y1,y2)

call calcabc(n,x,y1,na,aa1)
write(*,*)
call calcabc(n,x,y2,na,aa2)
!call calcabc2(n,x,y1,na,bb1)
!call calcabc2(n,x,y2,na,bb2)

call readarg(' Molar Concentration: ',narg,arg,concc)
read(concc,*) conc

call readarg(' Write output (yes/no)? [no]: ',narg,arg,concc)
write(*,'(/A$)') ' Write output (yes/no)?: '
if (concc(1:1).eq.'y'.or.concc(1:1).eq.'Y') then
  concc='yes'
else
  concc='no'
endif
write(*,'(x,A/)') trim(concc)
if (concc.eq.'yes') then
  p=index(filename,'.',back=.true.)-1
  if (p.le.0) p=len_trim(filename)
  pname=filename(1:p)//'.out'
  write(*,*) 'Output written to: ',trim(pname)
  call writedata(pname,n,x,y1,y2,na,aa1,aa2)
endif

call calcmu(na,aa1,conc,mu)
!call calcmu2(na,bb1,conc,mu2)
write(*,*) 'Chemical Potential (MU) for PA1: ',mu !,mu2
call calcmu(na,aa2,conc,mu)
!call calcmu2(na,bb2,conc,mu2)
write(*,*) 'Chemical Potential (MU) for PA2: ',mu !,mu2
end program

subroutine writedata(filename,n,x,y1,y2,na,aa1,aa2)
implicit none
integer na,n,i
real*8 x(n),y1(n),y2(n),aa1(na),aa2(na),mu,mu2
character filename*(*),line*4096

open (unit=1,file=filename)
do i=1,n
  call calcmu(na,aa1,x(i),mu)
  call calcmu(na,aa2,x(i),mu2)
  write(line,*) x(i),y1(i),mu,y2(i),mu2
  write(1,'(A)') trim(line)
enddo
close(1)

end subroutine

subroutine calcabc(n,x,y,na,a)
implicit none
integer n,i,j,k,ifault,nullty,na
real*8 x(n),y(n),a(na),px(na,na),py(na),xx,yy,imat(na,na),r2,sse,sst
real*8 aa(na*(na+1)/2),ww(na),cc(na*(na+1)/2)

px=0d0
py=0d0
sst=0d0
sse=0d0
do i=1,n
  xx=dsqrt(x(i))
  yy=y(i)
  sst=sst+y(i)**2
  sse=sse+y(i)
  do j=1,na
    py(j)=py(j)+yy*xx**j
    do k=j,na
      px(k,j)=px(k,j)+xx**j*xx**k
      if (j.ne.k) px(j,k)=px(k,j)
    enddo
  enddo
enddo

k=0
do i=1,na
  do j=1,i
    k=k+1
    aa(k)=px(i,j)
  enddo
enddo

call syminv ( aa, na, cc, ww, nullty, ifault )
if (ifault.ne.0) stop 'error inverting matrix'

k=0
do i=1,na
  do j=1,i
    k=k+1
    imat(i,j)=cc(k)
    if (i.ne.j) imat(j,i)=cc(k)
  enddo
enddo

a=matmul(py,transpose(imat))
sst=sst-sse**2/n
sse=0d0
do i=1,n
  call calcmu(na,a,x(i),yy)
  sse=sse+(yy-y(i))**2
enddo
! Determine the Coefficient of Determination (R2)
r2=1-sse/sst

do i=1,na
  write(*,'(A,I0,A$)') ' a(',i,')= '
  write(*,*) a(i)
enddo
write(*,*) 'R=',dsqrt(r2)
end subroutine

subroutine calcabc2(n,x,y,na,a)
implicit none
integer n,i,j,k,ifault,nullty,na
real*8 x(n),y(n),a(na),px(na,na),py(na),xx,yy,imat(na,na),r2,sse,sst
real*8 aa(na*(na+1)/2),ww(na),cc(na*(na+1)/2)

px=0d0
py=0d0
sst=0d0
sse=0d0
do i=1,n
  xx=dlog(1d0+x(i))
  yy=dlog(1d0-y(i))
  sst=sst+y(i)**2
  sse=sse+y(i)
  do j=1,na
    py(j)=py(j)+yy*xx**j
    do k=j,na
      px(k,j)=px(k,j)+xx**j*xx**k
      if (j.ne.k) px(j,k)=px(k,j)
    enddo
  enddo
enddo

k=0
do i=1,na
  do j=1,i
    k=k+1
    aa(k)=px(i,j)
  enddo
enddo

call syminv ( aa, na, cc, ww, nullty, ifault )
if (ifault.ne.0) stop 'error inverting matrix'

k=0
do i=1,na
  do j=1,i
    k=k+1
    imat(i,j)=cc(k)
    if (i.ne.j) imat(j,i)=cc(k)
  enddo
enddo

a=matmul(py,transpose(imat))
sst=sst-sse**2/n
sse=0d0
do i=1,n
  call calcmu2(na,a,x(i),yy)
  sse=sse+(yy-y(i))**2
enddo
! Determine the Coefficient of Determination (R2)
r2=1-sse/sst

do i=1,na
  write(*,'(A,I0,A$)') ' a(',i,')= '
  write(*,*) a(i)
enddo
write(*,*) 'R=',dsqrt(r2)
end subroutine

subroutine readdata(filename,n,x1,x2,x3)
implicit none
integer kode,n
character filename*(*)
real*8 x1(*),x2(*),x3(*)

open(unit=1,file=filename)
n=1
read(1,*,iostat=kode) x1(n),x2(n),x3(n)
do while (kode.eq.0)
  n=n+1
  read(1,*,iostat=kode) x1(n),x2(n),x3(n)
enddo
n=n-1
close(1)
end subroutine

subroutine calcmu(n,a,conc,mu)
implicit none
integer n,i
real*8 conc,mu,a(n),x
x=dsqrt(conc)
mu=0d0
do i=1,n
  mu=mu+a(i)*x**i
enddo
end subroutine

subroutine calcmu2(n,a,conc,mu)
implicit none
integer n,i
real*8 a(n),conc,mu,x,y
x=dlog(1d0+conc)
y=0d0
do i=1,n
  y=y+a(i)*x**i
enddo
mu=1d0-dexp(y)
end subroutine

subroutine cholesky ( a, n, nn, u, nullty, ifault )

!*****************************************************************************80
!
!! CHOLESKY computes the Cholesky factorization of a PDS matrix.
!
!  Discussion:
!
!    For a positive definite symmetric matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!    The missing initialization "II = 0" has been added to the code.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix
!    stored by rows in lower triangular form as a one dimensional array,
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, integer ( kind = 4 ) NN, the dimension of A, (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, if NN < (N*(N+1))/2.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(nn)
  real ( kind = 8 ), parameter :: eta = 1.0D-09
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) rsq
  real ( kind = 8 ) u(nn)
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
end
subroutine syminv ( a, n, c, w, nullty, ifault )

!*****************************************************************************80
!
!! SYMINV computes the inverse of a symmetric matrix.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    FORTRAN77 version by Michael Healy
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 7:
!    Inversion of a Positive Semi-Definite Symmetric Matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 198-199.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix stored
!    by rows in lower triangular form as a one dimensional array, in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) C((N*(N+1))/2), the inverse of A, or generalized
!    inverse if A is singular, stored using the same storage scheme employed
!    for A.  The program is written in such a way that A and U can share storage.
!
!    Workspace, real ( kind = 8 ) W(N).
!
!    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY is zero,
!    the matrix is judged to have full rank.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no error detected.
!    1, N < 1.
!    2, A is not positive semi-definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a((n*(n+1))/2)
  real ( kind = 8 ) c((n*(n+1))/2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mdiag
  integer ( kind = 4 ) ndiag
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x

  ifault = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  nrow = n
!
!  Compute the Cholesky factorization of A.
!  The result is stored in C.
!
  nn = ( n * ( n + 1 ) ) / 2

  call cholesky ( a, n, nn, c, nullty, ifault )

  if ( ifault /= 0 ) then
    return
  end if
!
!  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
!  of C, row by row starting with the last row.
!  IROW = the row number,
!  NDIAG = location of last element in the row.
!
  irow = nrow
  ndiag = nn

  do
!
!  Special case, zero diagonal element.
!
    if ( c(ndiag) == 0.0D+00 ) then

      l = ndiag
      do j = irow, nrow
        c(l) = 0.0D+00
        l = l + j
      end do

    else

      l = ndiag
      do i = irow, nrow
        w(i) = c(l)
        l = l + i
      end do

      icol = nrow
      jcol = nn
      mdiag = nn

      do

        l = jcol

        if ( icol == irow ) then
          x = 1.0D+00 / w(irow)
        else
          x = 0.0D+00
        end if

        k = nrow

        do while ( irow < k )

          x = x - w(k) * c(l)
          k = k - 1
          l = l - 1

          if ( mdiag < l ) then
            l = l - k + 1
          end if

        end do

        c(l) = x / w(irow)

        if ( icol <= irow ) then
          exit
        end if

        mdiag = mdiag - icol
        icol = icol - 1
        jcol = jcol - 1

      end do

    end if

    ndiag = ndiag - irow
    irow = irow - 1

    if ( irow <= 0 ) then
      exit
    end if

  end do

  return
end

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

prname='FIT-IT'
prver='version 1.0'
prdesc='Fit mu (chemical potential) vs concentration and extrapolate mu for given concentration'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='01 May 2012'
lastdate='11 Feb 2014'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

