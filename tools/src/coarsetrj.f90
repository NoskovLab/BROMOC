!    COARSETRJ - Convert atomistic DNA to coarse grained DNA
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

module comun
implicit none
integer nn,nion,tion,tnf
integer,allocatable :: resls(:),atlsf(:),atlsi(:),ras(:)
real*8,allocatable :: w(:),e(:)
integer,allocatable ::  iatom(:),ires(:),ions(:),csoli(:),csolf(:)
character*4,allocatable :: typ(:),res(:),segid(:),resid(:)
real*8 :: bv(3,3),bv2(3,3),invbv(3,3),invbv2(3,3)
real*8,allocatable :: rt(:,:)
real*8,allocatable :: cent(:,:,:),csol(:,:)
real*8 :: delr
character*4,allocatable :: fel(:),fl(:),gfl(:)
integer,allocatable :: gfi(:),gff(:),fi(:),ff(:),fg(:),fp(:)
integer gfn,fn,efn
integer na,nsc
integer rtn,dp(4),soln
integer,allocatable :: rti(:),rtf(:),rlrt(:)
integer agfpn,snd,afrn,frgn
integer, allocatable :: agfpl(:),grf(:,:),frg(:),frgg(:),pmd(:),pmdi(:,:),pmdf(:,:),sndi(:),sndf(:),afr(:),afri(:,:),afrf(:,:)
character*4,allocatable :: rtc(:)
logical*1 dcdopen,depablo,bvon
logical*1 cnvbox,center,masson,reass,charmm
contains

  function inintv(vec)
  implicit none
  real*8 vec(3),inintv(3)
  inintv(1)=inint(vec(1))
  inintv(2)=inint(vec(2))
  inintv(3)=inint(vec(3))
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

  function inint(num)
  implicit none
  real*8 inint,num
  inint=dfloat(iint(num+0.5d0))
  end function

  function iint(num)
  implicit none
  integer iint
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

  function pbc(pos)
  implicit none
  real*8 pbc(3),pos(3)
  pbc=pos-matmul(bv,inintv(matmul(invbv,matmul(pos,bv))))
  end function
  
  function inbox(pos)
  implicit none
  real*8 pos(3)
  logical*1 inbox
  inbox=all(int(inintv(matmul(invbv2,matmul(pos,bv2)))).eq.0)
  end function
end module

program coarsetrj
use comun
implicit none
integer i,j,k,fnt
character crdfile*256,dcdfile*256,mskfile*256,line*256,frmt*4
integer ul(256),ll(256),num,narg,arg
integer*1 stf,cf,wof,woft
character bs*8

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

bvon=.false.

call readarg('Input System Topology (Charmm .crd or .pdb) filename (filename format): ',narg,arg,line)
read(line,*) crdfile,frmt
if (frmt.eq.'pdb'.or.frmt.eq.'PDB') then
  stf=1
elseif (frmt.eq.'crd'.or.frmt.eq.'CRD') then
  stf=2
endif

call readarg('Input Coordinates .crd .pdb .dcd or .dcde filename (filename format): ',narg,arg,line)
read(line,*) dcdfile,frmt
if (frmt.eq.'pdb'.or.frmt.eq.'PDB') then
  cf=1
elseif (frmt.eq.'crd'.or.frmt.eq.'CRD') then
  cf=2
elseif (frmt.eq.'dcd'.or.frmt.eq.'DCD') then
  cf=3
elseif (frmt.eq.'dcde'.or.frmt.eq.'DCDE') then
  cf=4
endif

call readarg('Input Fragment Mask (.msk) filename [HELP and END]: ',narg,arg,mskfile)
if (len_trim(mskfile).eq.0) then
  write(*,'(A)') 'Order is important and must be 1)BASE 2)SUGAR 3)PHOSPATE'
  write(*,'(A)') 'Each record is composed of 4 characters'
  write(*,'(A)') 'Per line: 1st record -> group, 2nd record -> fragment, 3rd and higher record -> fragments elements'
  write(*,'(A)') 'Example of a mask file:'
  write(*,'(A)') "BASEADENN9  C5  N7  C8  H8  N1  C2  H2  N3  C4  C6  N6  H61 H62 "
  write(*,'(A)') "BASECYTON1  C6  H6  C5  H5  C2  O2  N3  C4  N4  H41 H42 "
  write(*,'(A)') "BASEGUANN9  C4  N2  H21 H22 N3  C2  N1  H1  C6  O6  C5  N7  C8  H8  "
  write(*,'(A)') "BASETHYMN1  C6  H6  C2  O2  N3  H3  C4  O4  C5  C5M H51 H52 H53 "
  write(*,'(A)') "SUGASUGAC5' H5' H5''C4' H4' O4' C1' H1' C2' H2''H2' C3' H3' "
  write(*,'(A)') "PHOSPHO P   "  
  write(*,'(A)') "PHOST5P H5T "
  stop
endif

write(*,'(A$)') 'Reading topology from '
! Read crd
if (stf.eq.1) then
  write(*,'(A)') 'pdb ...'
  call readpdbna(crdfile)
  write(*,'(A,I0)') 'Number of atoms: ',na
  call readpdb(crdfile)
elseif (stf.eq.2) then
  write(*,'(A)') 'CHARMM crd ...'
  call readcharmmcrd(crdfile)
  write(*,'(A,I0)') 'Number of atoms: ',na
endif

write(*,'(A)') 'Reading mask file ...'
!Read mask file
call readmask(mskfile)

write(*,'(A)') 'Building residue list ...'
! Build residue list
call reslist()

! reassing O3' to next residue
call readarg('Reassign O3 (y/n) [n] ? ',narg,arg,line)
reass=.false.
if (line(1:1).eq.'y'.or.line(1:1).eq.'Y') then
  reass=.true.
  call reassigno3()
  call reslist()
endif

write(*,'(A,I0)') 'Read residues: ',nn

call readarg('Locate coarse bases at geometric/mass center (y/n) [n]? ',narg,arg,line)
if (line(1:1).eq.'Y'.or.line(1:1).eq.'y') then
  depablo=.false.
else
  depablo=.true.
  write(*,'(A)') 'Turning on de Pablo et al coarse grained DNA'
  write(*,'(A)') 'Select the atomtype where to center the base (usually N1 for G/A & N3 for C/T)'
  do j=1,4
    k=ff(j)-fi(j)+1
    write(line,'(A,I0,A)') '(4x,',k,'(I4,x))' 
    write(*,line) (i,i=1,k) 
    write(*,*) fl(j),' ',(fel(i),' ',i=fi(j),ff(j))
    if (fl(j)(1:1).eq.'A') dp(j)=5
    if (fl(j)(1:1).eq.'C') dp(j)=6
    if (fl(j)(1:1).eq.'G') dp(j)=6
    if (fl(j)(1:1).eq.'T') dp(j)=5
    call readarg('',narg,arg,line)
    if (len_trim(line).gt.0) read(line,*) dp(j)
    write(*,'(I0,A,A)') dp(j),' ',fel(fi(j)-1+dp(j))
  enddo
endif

call readarg('Use center of mass to define center of coarse particles (otherwise will use geometric center) (y/n) [y]? ',narg,arg,line)
if (line(1:1).eq.'N'.or.line(1:1).eq.'n') then
  masson=.false.
else
  masson=.true.
endif

! Buid residue type list
call restyplist()

! Inform residue types read and enquire more info'
write(*,'(A)') 'Residue types read: '
write(*,'(A)') 'Number Symbol Residues'
do i=1,rtn
  write(*,'(I5,3x,A4,3x,I5)') i,rtc(i),rtf(i)-rti(i)+1
enddo

call readarg('Define residue type number sequence of ions/solvent to consider: ',narg,arg,line)
if (len_trim(line).eq.0) then
  nion=0
else
  call findparm(line,256,num,ll,ul)
  nion=num
  allocate (ions(nion))
endif
tion=0
do i=1,nion
  read(line(ll(i):ul(i)),*) ions(i)
  tion=tion+rtf(ions(i))-rti(ions(i))+1
enddo
if (tion.gt.0.and.nion.gt.0) allocate (csol(1:3,tion),csoli(nion),csolf(nion))

! Build fragment list and allocate cent var
call fraglist()
write(*,'(A,I0)') 'Total Nucleotides fragments: ',agfpn
  write(*,'(A)') 'Fragment Group Number | Fragment Group Label | Fragment Label | Number of Elements | Found Fragments' 
do i=1,fn
  write(*,'(I5,3x,A4,3x,A4,3x,I5,3x,I5)') fg(i),gfl(fg(i)),fl(i),ff(i)-fi(i)+1,fp(i)
enddo

!Read dcd header and first frame
write(*,'(A$)') 'Reading coordinates from '
if (cf.eq.1) then
  write(*,'(A)') 'pdb ...'
  call readpdbtnf(crdfile)
  call openpdbcoor(crdfile,2)
  call readpdbcoor(2)
elseif (cf.eq.2) then
  write(*,'(A)') 'CHARMM crd ...'
  tnf=1
  call opencharmmcoor(crdfile,2)
  call readcharmmcoor(2)
elseif (cf.eq.3) then
  write(*,'(A)') 'dcd ...'
  call readdcdhead(dcdfile,2)
  call readdcdbody(2)
elseif (cf.eq.4) then
  write(*,'(A)') 'dcde ...'
  call readdcdhead(dcdfile,2)
  call readdcdebody(2)
endif
if (reass) rt=rt(1:3,ras(1:na)) 

call readarg('Center DNA to centroid (y/n) [n]? ',narg,arg,line)
if (line(1:1).eq.'Y'.or.line(1:1).eq.'y') then
  center=.true.
else
  center=.false.
endif

call readarg('Write PDB, CRD or XYZ format [xyz]? ',narg,arg,line)
if (line(1:3).eq.'pdb'.or.line(1:3).eq.'PDB') then
  wof=1
elseif (line(1:3).eq.'crd'.or.line(1:3).eq.'CRD') then
  wof=3
else 
  wof=2
endif

call readarg('Write DCD [none]? ',narg,arg,line)
if (line(1:3).eq.'dcd'.or.line(1:3).eq.'DCD') then
  woft=1
else
  woft=0
endif

if (wof==2) then
  call readarg('Convert Box to different box size (y/n) [n]? ',narg,arg,line)
  if (line(1:1)=='y'.or.line(1:1)=='Y') then
    cnvbox=.true.
    write(*,*) 'Actual Periodic Box Vectors: '
    write(*,*) bv
    call readarg('Box Vectors (ax ay az bx by bz cx cy cz) [ENTER to use actual]: ',narg,arg,line)
    if (len_trim(line).ne.0) then 
      bvon=.true.
      read(line,*) bv(1,1),bv(2,1),bv(3,1),bv(1,2),bv(2,2),bv(3,2),bv(1,3),bv(2,3),bv(3,3)
    endif
    if (bvon) invbv=invmat(matmul(transpose(bv),bv))
    call readarg('New Box Vectors (ax ay az bx by bz cx cy cz): ',narg,arg,line)
    read(line,*) bv2(1,1),bv2(2,1),bv2(3,1),bv2(1,2),bv2(2,2),bv2(3,2),bv2(1,3),bv2(2,3),bv2(3,3)
    invbv2=invmat(matmul(transpose(bv2),bv2))
  else
    cnvbox=.false.
  endif
else
  cnvbox=.false.
endif
write(*,'(A,I0)') 'Number of frames: ',tnf
call readarg('Choose frame number: ',narg,arg,line)
read(line,*) fnt

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc
if (wof.eq.1) open(unit=3,file='coarse.pdb') 
if (wof.eq.2) open(unit=3,file='coarse.xyz')
if (woft.eq.1) call writedcdhead('coarse.dcd',4)
do while (dcdopen.and.nsc.le.tnf)                  ! open loop for each frame
  call centall()                                   ! Compute geometric center (centroid)
  write(*,'(A8,I8$)') bs,nsc             ! print frame number
  if (wof.eq.1.and.nsc.eq.fnt) call writeout(3)    ! writeout
  if (wof.eq.2.and.nsc.eq.fnt) call writexyz(3)    ! writexyz
  if (wof.eq.3.and.nsc.eq.fnt) call writecharmmcrd('notcoarsed.cor',3)
  if (woft.eq.1) call writedcdbody(4)
  if (cf.eq.1) then
    call readpdbcoor(2)
  elseif (cf.eq.2) then
    call readcharmmcoor(2)
  elseif (cf.eq.3) then
    call readdcdbody(2)
  elseif (cf.eq.4) then
    call readdcdebody(2)
  endif
  if (reass) rt=rt(1:3,ras(1:na))
enddo
close(2)
close(3) ! close writeout
if (woft.eq.1) close(4)
write(*,'(/A)') 'Normal termination of DNACDF'
end program

subroutine writexyz(unitn)
use comun
implicit none
integer i,j,k,l,nnn,unitn,nq,nqq
character*4,allocatable ::  nam(:),namq(:)
real*8,allocatable :: rq(:,:),rqq(:,:)
real*8 cnt(3),dv(3)

nnn=gfn*agfpn
do i=1,nion
  nnn=nnn+atlsf(rtf(ions(i)))-atlsi(rti(ions(i)))+1
enddo
allocate (rq(3,na),rqq(3,27*na),nam(na),namq(27*na))

nq=0
cnt=0d0
do i=1,gfn
  do j=1,agfpn
    nq=nq+1
    nam(nq)=fl(grf(i,j))
    rq(1:3,nq)=cent(1:3,i,j)
    cnt(1:3)=cnt(1:3)+rq(1:3,nq)
  enddo
enddo
cnt=cnt/nq
do i=1,nion
  do j=atlsi(rti(ions(i))),atlsf(rtf(ions(i)))
    nq=nq+1
    nam(nq)=res(j)
    rq(1:3,nq)=rt(1:3,j)
  enddo
enddo

if (center) then
  do i=1,nq
    rq(1:3,i)=rq(1:3,i)-cnt(1:3)
  enddo
endif

if (bvon) then
  do i=1,nq
    rq(1:3,i)=pbc(rq(1:3,i))
  enddo
endif
if (cnvbox) then
  nqq=0
  do i=-1,1,1
    do j=-1,1,1
      do k=-1,1,1
        do l=1,nq
          dv=rq(1:3,l)+i*bv(1:3,1)+j*bv(1:3,2)+k*bv(1:3,3)
          if (inbox(dv)) then
            nqq=nqq+1
            namq(nqq)=nam(l)
            rqq(1:3,nqq)=dv(1:3)
          endif
        enddo
      enddo
    enddo
  enddo
else
  nqq=nq
  rqq(1:3,1:nqq)=rq(1:3,1:nq)
  namq(1:nqq)=nam(1:nq)
endif

! write xyz
write(unitn,*) nqq
write(unitn,*) 
do i=1,nqq
  write(unitn,*) namq(i),rqq(1:3,i)
enddo

deallocate (rq,rqq,nam,namq)

end subroutine

subroutine writeout(unitn)
use comun
implicit none
integer i,j,k,unitn,nawo,c,d,a,b
integer,allocatable :: rnwo(:)
real*8,allocatable :: rwo(:,:)
character,allocatable :: rtwo(:)*5,atwo(:)*5
integer,allocatable :: conwo(:,:)

if (nion.gt.0) then
  nawo=gfn*agfpn+csolf(nion)
else
  nawo=gfn*agfpn
endif
allocate (rnwo(nawo),rwo(3,nawo),rtwo(nawo),atwo(nawo),conwo(nawo,2))

k=0
do i=1,gfn
  do j=1,agfpn
    k=k+1
    rnwo(k)=i
    rtwo(k)=gfl(i)
    atwo(k)=fl(grf(i,j))
    rwo(1:3,k)=cent(1:3,i,j)
  enddo
enddo
b=gfn
do j=1,nion
  b=b+1
  c=csoli(j) ! first atom of the selected ion type
  d=csolf(j) ! last atom of the selected ion type
  do i=c,d
    k=k+1
    a=atlsi(rti(ions(j)))
    rnwo(k)=b
    rtwo(k)=res(a)
    atwo(k)=res(a)
    rwo(1:3,k)=csol(1:3,i)
  enddo
enddo
call connect(conwo,nawo)
call writepdb(unitn,nawo,rwo,atwo,rtwo,rnwo,conwo)
deallocate (rnwo,rwo,rtwo,atwo,conwo)
end subroutine

subroutine connect(con,nawo)
use comun
implicit none
integer i,j,nawo
integer con(nawo,2)

con=0

do i=1,agfpn 
  con(i,1)=i+agfpn
  con(i+agfpn,1)=i+2*agfpn
enddo
do i=1,snd
  do j=sndi(i),sndf(i)-1
    con(j+agfpn,2)=j+agfpn*2+1
  enddo
enddo
end subroutine

subroutine writepdb(unitn,na,r,at,rt,rn,con)
implicit none
integer i,na,rn(na),unitn
real*8 r(3,na)
character rt(na)*5,at(na)*5
integer con(na,2)

do i=1,na 
  write (unitn,'(A6,I5,x,A5,A5,I4,4x,3F8.3)') 'ATOM  ',i,at(i),rt(i),rn(i),r(1,i),r(2,i),r(3,i)
enddo
do i=1,na
  if (con(i,2).ne.0) then
    write(unitn,'(a6,14i5)') 'CONECT',i,con(i,1),con(i,2)
  elseif (con(i,1).ne.0) then
    write(unitn,'(a6,14i5)') 'CONECT',i,con(i,1)
  endif
enddo
write (unitn,'(A)') 'END'
end subroutine

! fn is the number of fragments
! gfn is the number of group fragments
! efn is the number of elements 
subroutine readmask(mskfile)
use comun
implicit none
integer i,j,k,kk,s,kode
character line*256,mskfile*256,gflt*4

open (unit=1,file=mskfile,iostat=kode)
j=0
efn=0
gfn=0
gflt=''
do while (kode.eq.0)
  line=''
  read(1,'(A)',iostat=kode) line
  k=len_trim(line)-8 ! because the first 8 chars correspond to the fragment group label and fragment label
  if (k.gt.0) then 
    j=j+1
    if (mod(k,4).ne.0) then 
      k=int(k/4)+1
    else
      k=int(k/4)
    endif
    efn=efn+k
    if (gflt.ne.line(1:4)) gfn=gfn+1
    gflt=line(1:4)
  endif
enddo
fn=j
rewind(1)
allocate (gfi(gfn),gff(gfn),fi(fn),ff(fn),fel(efn),fl(fn),gfl(gfn),fg(fn),fp(fn))
kode=0
j=0
do while (kode.eq.0)
  line=''
  read(1,'(A)',iostat=kode) line
  k=len_trim(line)-8   ! substracting the group fragment label and fragment label
  if (k.gt.0) then
    j=j+1
    if (mod(k,4).ne.0) then
      k=int(k/4)+1
    else
      k=int(k/4)
    endif
    if (j.eq.1) then 
      i=1
      gfl(i)=line(1:4)
      gfi(i)=j
      fl(j)=line(5:8)
      fi(j)=1
      ff(j)=k
      fg(j)=i
      gff(i)=j
    else
      if (gfl(i).ne.line(1:4)) then
        gff(i)=j-1
        i=i+1
        gfl(i)=line(1:4)
        gfi(i)=j
      endif
      fg(j)=i
      fl(j)=line(5:8)
      fi(j)=ff(j-1)+1
      ff(j)=ff(j-1)+k
    endif
    gff(i)=j
    do s=fi(j),ff(j)
      kk=s-fi(j)+1+2
      fel(s)=line(kk*4-3:kk*4) ! reading fragment element label
    enddo
  endif
enddo
close(1)
end subroutine

subroutine findparm(str,dimn,num,llim,ulim)
implicit none
integer dimn,i,length,num
integer ulim(dimn),llim(dimn)
character ( len = dimn ) str
logical chng
length=len_trim(str)
chng=.false.
num=0
ulim(1:dimn)=0
llim(1:dimn)=0
do i=1,length
  if (iachar(str(i:i)).le.32.or.iachar(str(i:i)).ge.127) then
    if (chng) ulim(num)=i-1
    chng=.false.
  else
    if (.not.chng) then
      num=num+1
      llim(num)=i
    endif
    chng=.true.
  endif
enddo
if (ulim(num).eq.0) ulim(num)=length
end subroutine

subroutine readcharmmcrd(crdfile)
use comun
implicit none
integer i,num
character line*1024 ! ,frmt*64
integer ul(1024),ll(1024)
character crdfile*256
logical loopon
loopon=.true.
!frmt='(I10,I10,1X,A4,1X,A4,3F10.5,1X,A4,1X,a4,F10.5)'
!         1         1  ADE       H5T           -42.6356481721        8.6779119410       -0.3240608554  DNAA      1               0.0000000000
open(unit=1,file=crdfile)

do while (loopon)
  read(1,'(A)') line
  if (line(1:1).ne.'*') then
    call findparm(line,1024,num,ll,ul)
    read(line(ll(1):ul(1)),*) na
    loopon=.false.
  endif
enddo
allocate (iatom(na),ires(na),typ(na),res(na),segid(na),resid(na),rt(3,na),w(na))
do i=1,na
  read(1,'(A)') line
  call findparm(line,1024,num,ll,ul)
!  read(line(ll(1):ul(1)),*) iatom(i)
  iatom(i)=i
  read(line(ll(2):ul(2)),*) ires(i)
  res(i)=line(ll(3):ul(3))
  typ(i)=line(ll(4):ul(4))
  read(line(ll(5):ul(5)),*) rt(1,i)
  read(line(ll(6):ul(6)),*) rt(2,i)
  read(line(ll(7):ul(7)),*) rt(3,i)
  segid(i)=line(ll(8):ul(8))
  resid(i)=line(ll(9):ul(9))
  read(line(ll(10):ul(10)),*) w(i)
enddo
close(1)
end subroutine

subroutine reassigno3()
use comun
implicit none
integer i,j,k,l,m,rea(na),rsn
real*8 rttmp(3),wtmp,etmp
integer irestmp
character*4 typtmp,restmp,segidtmp,residtmp

allocate (ras(na))
rsn=0
do i=1,na
  if (typ(i).eq."O3' ") then
    rsn=rsn+1
    rea(rsn)=i
  endif
  ras(i)=i
enddo

do i=1,rsn
  l=rea(i)
  do j=1,nn
    if (ires(l).eq.resls(j)) then
      k=atlsf(j)
      ras(l)=k
      ras(k)=l
      if (j.eq.nn) then 
        irestmp=resls(nn)+1
        restmp=res(atlsf(nn))
      else
        m=atlsi(j+1)
        irestmp=resls(j+1)
        restmp=res(m)
      endif
      typtmp=typ(l)
      rttmp=rt(1:3,l)
      segidtmp=segid(l)
      residtmp=resid(l)
      wtmp=w(l)
      etmp=e(l)

      ires(l)=ires(k)
      res(l)=res(k)
      typ(l)=typ(k)
      rt(1:3,l)=rt(1:3,k)
      segid(l)=segid(k)
      resid(l)=resid(k)
      w(l)=w(k)
      e(l)=e(k)

      ires(k)=irestmp
      res(k)=restmp
      typ(k)=typtmp
      rt(1:3,k)=rttmp
      segid(k)=segidtmp
      resid(k)=residtmp
      w(k)=wtmp
      e(k)=etmp
      exit
    endif
  enddo
enddo
end subroutine

subroutine opencharmmcoor(crdfile,u)
use comun
implicit none
integer u
character crdfile*256

open(unit=u,file=crdfile)
nsc=0
dcdopen=.true.
end subroutine

subroutine readcharmmcoor(unidad)
use comun
implicit none
integer i,num,nan,kode,unidad
character line*1024 ! ,frmt*64
integer ul(1024),ll(1024)
logical loopon

kode=0
loopon=.true.
do while (loopon.and.kode.eq.0)
  read(unidad,'(A)',iostat=kode) line
  if (kode.eq.0.and.line(1:1).ne.'*') then
    call findparm(line,1024,num,ll,ul)
    read(line(ll(1):ul(1)),*) nan
    loopon=.false.
  endif
enddo
if (kode.eq.0) then
  do i=1,na
    read(unidad,'(A)',iostat=kode) line
    call findparm(line,1024,num,ll,ul)
    read(line(ll(5):ul(5)),*) rt(1,i)
    read(line(ll(6):ul(6)),*) rt(2,i)
    read(line(ll(7):ul(7)),*) rt(3,i)
  enddo
endif
if (kode.eq.0) then
  nsc=nsc+1
else 
  close(unidad)
  dcdopen=.false.
endif
end subroutine

subroutine writecharmmcrd(crdfile,u)
use comun
implicit none
integer i,j,u
character crdfile*(*),frmt*64
frmt='(I10,I10,2x,A4,6X,A4,4x,3F20.10,2X,A4,6X,A4,4x,F20.10)'
!         1         1  ADE       H5T           -42.6356481721        8.6779119410       -0.3240608554  DNAA      1               0.0000000000
open(unit=u,file=trim(adjustl(crdfile)))
write(u,'(I10,A)') na,'  EXT'
do i=1,na
  write(u,frmt) iatom(i),ires(i),res(i),typ(i),(rt(j,i),j=1,3),segid(i),resid(i),w(i)
enddo
close(u)
end subroutine

subroutine writepdbnew(crdfile,con)
use comun
implicit none
integer i,j,con(na,2)
character crdfile*256,frmt*64
frmt='(A6,I5,x,A5,A5,I4,4x,3F8.3,2F6.2,6x,A4)'
!ATOM      5 H5'' CYT     1     -25.455  -1.135  -2.002  1.00  0.00      DNAA
open(unit=67,file=crdfile)
do i=1,na
  write(67,frmt) 'ATOM  ',iatom(i),typ(i),res(i),ires(i),(rt(j,i),j=1,3),1.00,w(i),segid(i)
enddo
do i=1,na
  if (con(i,2).ne.0) then
    write(67,'(a6,14i5)') 'CONECT',i,con(i,1),con(i,2)
  elseif (con(i,1).ne.0) then
    write(67,'(a6,14i5)') 'CONECT',i,con(i,1)
  endif
enddo
write (67,'(A)') 'END'

close(67)
end subroutine

subroutine readdcdebody(u)
use comun
implicit none
integer kode,i,u
real*8 xtlabc(12)
real*4 rtt(3,na)

if (.not.charmm) read(u,iostat=kode) xtlabc
read(u,iostat=kode) (rtt(1,i),i=1,na)
read(u,iostat=kode) (rtt(2,i),i=1,na)
read(u,iostat=kode) (rtt(3,i),i=1,na)
if (kode.eq.0) then
  nsc=nsc+1
  rt(1:3,1:na)=dble(rtt(1:3,1:na))
  !Build Box vectors
  bv(1,1)=xtlabc(1)
  bv(2,1)=xtlabc(2)
  bv(3,1)=xtlabc(3)
  bv(1,2)=xtlabc(4)
  bv(2,2)=xtlabc(5)
  bv(3,2)=xtlabc(6)
  bv(1,3)=xtlabc(7)
  bv(2,3)=xtlabc(8)
  bv(3,3)=xtlabc(9)
  if (.not.charmm) bvon=.true.
  invbv=invmat(matmul(transpose(bv),bv))
else
  close(u)
  dcdopen=.false.
endif
end subroutine

subroutine readdcdhead(dcdfile,u)
use comun
implicit none
integer icntrl(20),itemp,kode,u
integer nfile,npriv,nsavc,nstep,nfree
integer natom,ntitle
character dcdfile*256,hdr*4
character*1,allocatable :: title(:)

open(unit=u,file=dcdfile,form='unformatted')
dcdopen=.true.
read(u) hdr,icntrl
ntitle=icntrl(20)/12*80
if (allocated(title)) deallocate (title)
allocate (title(ntitle))
title=''
read(u,iostat=kode) itemp,title
read(u) natom
if (na.ne.natom) stop 'Number of atoms differ between .crd and .dcd'
nfile=icntrl(1)
npriv=icntrl(2)
nsavc=icntrl(3)
nstep=icntrl(4)
!if(icntrl(9).gt.0) print *, '# fixed atoms = ',icntrl(9)
nfree = natom-icntrl(9)
!print *, '# of free atoms = ',nfree
!print *, 'total # atom = ', natom,nstep,nsavc
charmm=.false.
if (icntrl(2).eq.0) charmm=.true.
if (nstep.le.0) nstep=1
if (nsavc.le.0) nsavc=1
tnf = nstep/nsavc

write(*,'(A,I0)') 'Total number of frames: ',tnf
nsc=0
end subroutine

subroutine readdcdbody(u)
use comun
implicit none
integer kode,i,u
real*8 xtlabc(6)
real*4 rtt(3,na)

if (.not.charmm) read(u,iostat=kode) xtlabc
read(u,iostat=kode) (rtt(1,i),i=1,na)
read(u,iostat=kode) (rtt(2,i),i=1,na)
read(u,iostat=kode) (rtt(3,i),i=1,na)
if (kode.eq.0) then  
  nsc=nsc+1
  rt(1:3,1:na)=dble(rtt(1:3,1:na))
  !Build Box vectors
  bv(1,1)=xtlabc(1)
  bv(2,1)=xtlabc(2)
  bv(3,1)=xtlabc(4)
  bv(1,2)=xtlabc(2)
  bv(2,2)=xtlabc(3)
  bv(3,2)=xtlabc(5)
  bv(1,3)=xtlabc(4)
  bv(2,3)=xtlabc(5)
  bv(3,3)=xtlabc(6)
  if (.not.charmm) bvon=.true.
  invbv=invmat(matmul(transpose(bv),bv))
else
  close(u)
  dcdopen=.false.
endif
end subroutine

!resl(nn) list that contains the residue number in each element (per residue)
!nn number of residues
!atlsi initial atom position per residue
!atlsf final atom position per residue
subroutine reslist()
use comun
implicit none
integer rl(na),alf(na),i,ali(na)
if (allocated(atlsi)) deallocate (atlsi)
if (allocated(resls)) deallocate (resls)
if (allocated(atlsf)) deallocate (atlsf)
if (allocated(rlrt)) deallocate (rlrt)
rl(1)=ires(1)
ali(1)=1
alf(1)=1
nn=1
do i=2,na
  if (ires(i).ne.rl(nn)) then
    alf(nn)=i-1
    nn=nn+1
    ali(nn)=i
    rl(nn)=ires(i)
    alf(nn)=i
  endif
enddo
alf(nn)=na

allocate (atlsi(nn),resls(nn),atlsf(nn),rlrt(nn)) 

resls(1:nn)=rl(1:nn)
atlsi(1:nn)=ali(1:nn)
atlsf(1:nn)=alf(1:nn)
end subroutine

!rtn is the residue type number (1 = DNA, 2 = Na, 3 = Cl, 4 = WAT)
!rlrt is a residue list sorted by residue type
!rti(rtn) is the initial position of rtn in the residue list sorted by residue type
!rtf(rtn) is the final position of rtn in the residue list sorted by residue type
!rtc(rtn) is the residue type character (ex: DNA )
subroutine restyplist()
use comun
implicit none 
integer i,j,k,mat(nn,nn),din(nn)
character*4 rtct(nn)
logical sta

rtn=1
din(rtn)=1
rtct(rtn)=res(1)
mat(rtn,din(rtn))=1
do i=2,nn
   sta=.false.
   j=rtn
   do while (j.gt.0.and..not.sta)
     if (res(atlsi(i)).eq.rtct(j)) then
       sta=.true.
       din(j)=din(j)+1
       mat(j,din(j))=i
     else
       j=j-1
     endif
   enddo
   if (.not.sta) then
     rtn=rtn+1
     din(rtn)=1
     rtct(rtn)=res(atlsi(i))
     mat(rtn,din(rtn))=i
   endif
enddo

allocate (rtf(rtn),rtc(rtn),rti(rtn)) 
rtc(1:rtn)=rtct(1:rtn)
k=0
do i=1,rtn
  rti(i)=k+1
  do j=1,din(i)
    k=k+1
    rlrt(k)=mat(i,j)
  enddo
  rtf(i)=k   
enddo
end subroutine

subroutine fraglist()
use comun
implicit none
integer h,i,j,k,l,ii
character sgid*4
logical frgpres(fn,nn),gfpres(gfn,nn),agfpres(nn),pres,apres

do i=1,nn
  do k=1,fn
    apres=.true.
    l=fi(k)-1
    do while (apres.and.l.lt.ff(k))
      l=l+1
      pres=.false.
      j=atlsi(i)-1
      do while (.not.pres.and.j.lt.atlsf(i))
        j=j+1
        if (typ(j).eq.fel(l)) pres=.true.
      enddo
      apres=apres.and.pres
    enddo
    frgpres(k,i)=apres ! logical list ordered by residue number and fragment number that is true when the fragment number is present in the residue number
  enddo
enddo 

do i=1,nn
  do j=1,gfn
    pres=.false.
    k=gfi(j)-1
    do while (.not.pres.and.k.lt.gff(j))
      k=k+1
      pres=pres.or.frgpres(k,i)
    enddo
    gfpres(j,i)=pres ! logical list ordered by residue number and group number that is true when the group number is present in the residue number
  enddo
enddo

agfpn=0
do i=1,nn
  apres=.true.
  j=0
  do while (apres.and.j.lt.gfn)
    j=j+1
    apres=apres.and.gfpres(j,i)
  enddo
  agfpres(i)=apres              ! logical list ordered by residue number that gives true when all group fragments are present
  if (agfpres(i)) agfpn=agfpn+1 ! number of residues with all the group fragments present (active residues)
enddo

allocate (agfpl(agfpn),grf(gfn,agfpn))

j=0
do i=1,nn
  if (agfpres(i)) then
    j=j+1
    agfpl(j)=i  ! agfpl(1:agfpn)->residue number: is a list that gives the residue numbers which have all the group fragment present (active residues)
  endif
enddo

sgid=segid(atlsi(agfpl(1)))
snd=1
do i=2,agfpn
  if (sgid.ne.segid(atlsi(agfpl(i)))) snd=snd+1 ! snd is the number of strands
  sgid=segid(atlsi(agfpl(i)))
enddo

allocate (sndi(snd),sndf(snd))

! set the begining and end of each strand in agfpl order
sgid=segid(atlsi(agfpl(1)))
snd=1
sndi=1
do i=2,agfpn
  if (sgid.ne.segid(atlsi(agfpl(i)))) then 
    sndf(snd)=i-1 ! end
    snd=snd+1
    sndi(snd)=i   ! begining
  endif
  sgid=segid(atlsi(agfpl(i)))
enddo
sndf(snd)=agfpn

fp=0
do i=1,gfn
  do j=1,agfpn
    k=gfi(i)-1
    apres=.true.
    do while (apres.and.k.lt.gff(i))
      k=k+1
      if (frgpres(k,agfpl(j))) then 
        grf(i,j)=k     ! matrix containing the fragment number based on the group number and the active residues
        fp(k)=fp(k)+1  ! number of elements found for each fragment
        apres=.false.
      endif
    enddo
  enddo
enddo

frgn=0
do i=1,gfn
  h=0
  do j=gfi(i),gff(i)
    if (fp(j).gt.0) h=h+1 
  enddo
  if (h.gt.1) frgn=frgn+h ! number of fragments belonging to groups of more that 1 fragment
enddo

allocate (frg(frgn),frgg(frgn),pmd(2*gfn*agfpn),pmdi(gfn+frgn,snd),pmdf(gfn+frgn,snd))

k=0
do i=1,gfn
  h=0
  do j=gfi(i),gff(i)
    if (fp(j).gt.0) h=h+1 
  enddo
  if (h.gt.1) then
    do j=gfi(i),gff(i)
      k=k+1
      frg(k)=j  ! list with fragments belonging to groups of more that 1 present fragment -> fragment number
      frgg(k)=i ! list with fragments belonging to groups of more that 1 present fragment -> group fragment number
    enddo
  endif
enddo

!pmd(l)=j contains the position in the agfpl list of each fragment>1 per strand
!pmdi(i,k) says the initial position in pmd where is located the fragment>1(i) and strand (k)
!pmdf(i,k) says the final position in pmd where is located the fragment>1(i) and strand (k)

l=0
do k=1,snd
h=0
  do i=1,gfn
    h=h+1
    pmdi(h,k)=l+1
    do j=sndi(k),sndf(k)
      l=l+1
      pmd(l)=j
    enddo
    pmdf(h,k)=l
  enddo
  do i=1,frgn
    h=h+1
    pmdi(h,k)=l+1
    do j=sndi(k),sndf(k)
      if (grf(frgg(i),j).eq.frg(i)) then
        l=l+1
        pmd(l)=j
      endif
    enddo
    pmdf(h,k)=l
  enddo
enddo

allocate (afr(na),afri(1:gfn,1:agfpn),afrf(1:gfn,1:agfpn))

afrn=0
do i=1,gfn
  do j=1,agfpn
    k=grf(i,j) ! gives the corresponding fragment number
    ii=agfpl(j) ! gives the corresponding residue number
    apres=.true.
    l=fi(k)-1
    afri(i,j)=afrn+1 ! points to the begining of gfn and agfpn in afr list
    do while (apres.and.l.lt.ff(k))
      l=l+1
      pres=.false.
      h=atlsi(ii)-1
      do while (.not.pres.and.h.lt.atlsf(ii))
        h=h+1
        if (typ(h).eq.fel(l)) then 
          pres=.true.
          afrn=afrn+1
          afr(afrn)=h ! contains the atom number list ordered by gfn, agfn
        endif
      enddo
      apres=apres.and.pres
    enddo
    afrf(i,j)=afrn  ! points to the end of the gfn and agfpn in afr list
  enddo
enddo

allocate (cent(1:3,1:gfn,1:agfpn))
end subroutine

subroutine centmass(r,n,c,at)
implicit none
integer i,n
real*8 r(3,n),c(3),m(n),mass
character*(*) at(n)

do i=1,n
  if (at(i)(1:3).eq.'CLA') then
    m(i)=35.4527d0
  elseif (at(i)(1:3).eq.'POT') then
    m(i)=39.0983d0
  elseif (at(i)(1:1).eq.'C') then
    m(i)=12.011d0
  elseif (at(i)(1:1).eq.'O') then
    m(i)=15.9994d0
  elseif (at(i)(1:1).eq.'H') then
    m(i)=1.00794d0
  elseif (at(i)(1:1).eq.'N') then
    m(i)=14.00674d0
  elseif (at(i)(1:1).eq.'P') then
    m(i)=30.973762d0
  elseif (at(i)(1:1).eq.'S') then
    m(i)=32.066d0
  else 
    m(i)=1.0d0
  endif
enddo
c(1:3)=0d0
mass=0d0
do i=1,n
  c(1)=c(1)+r(1,i)*m(i)
  c(2)=c(2)+r(2,i)*m(i)
  c(3)=c(3)+r(3,i)*m(i)
  mass=mass+m(i)
enddo
c(1)=c(1)/mass
c(2)=c(2)/mass
c(3)=c(3)/mass

end subroutine

subroutine centroid(r,n,c)
implicit none
integer i,n
real*8 r(3,n),c(3)

c(1:3)=0.0
do i=1,n
  c(1)=c(1)+r(1,i)
  c(2)=c(2)+r(2,i)
  c(3)=c(3)+r(3,i)
enddo
c(1)=c(1)/dfloat(n)
c(2)=c(2)/dfloat(n)
c(3)=c(3)/dfloat(n)

end subroutine

subroutine centall()
use comun
implicit none
integer i,j,k,c,d
if (depablo) then
  i=1
  do j=1,agfpn
    c=grf(i,j)
    k=dp(c)-1+afri(i,j)
    cent(1:3,i,j)=rt(1:3,afr(k))
  enddo
else
  i=1
  do j=1,agfpn ! fragment number
    c=afri(i,j)
    d=afrf(i,j)
    if (masson) then 
      call centmass(rt(1:3,afr(c:d)),d-c+1,cent(1:3,i,j),typ(afr(c:d)))
    else
      call centroid(rt(1:3,afr(c:d)),d-c+1,cent(1:3,i,j))
    endif
  enddo
endif
do i=2,gfn                ! residue number per residue type
  do j=1,agfpn ! fragment number
    c=afri(i,j)
    d=afrf(i,j)
    if (masson) then
      call centmass(rt(1:3,afr(c:d)),d-c+1,cent(1:3,i,j),typ(afr(c:d)))
    else
      call centroid(rt(1:3,afr(c:d)),d-c+1,cent(1:3,i,j))
    endif
  enddo
enddo

k=0
do i=1,nion ! each solvent/ion type
  csoli(i)=k+1
  do j=rti(ions(i)),rtf(ions(i)) ! all ions/solvent for each ion/solvent type
    k=k+1
    c=atlsi(j)
    d=atlsf(j)
    if (masson) then
      call centmass(rt(1:3,c:d),d-c+1,csol(1:3,k),typ(c:d))
    else
      call centroid(rt(1:3,c:d),d-c+1,csol(1:3,k))
    endif
  enddo
  csolf(i)=k
enddo

end subroutine

! check if a character is a letter
function chklet(s)
implicit none
character s*1
logical chklet

chklet=.true.
if ((iachar(s).ge.65.and.iachar(s).le.90).or.(iachar(s).ge.97.and.iachar(s).le.122)) then
  chklet=chklet.and..true.
else
  chklet=chklet.and..false.
endif
end function

! get atomtype
subroutine getat(s,d,at)
implicit none
integer d,k
logical chklet
character at*2
character ( len = d ) s
k=0
do while(k.lt.d)
  k=k+1
  if (chklet(s(k:k))) then
    at(1:1)=s(k:k)
    if (chklet(s(k+1:k+1))) then
      at(2:2)=s(k+1:k+1)
    else
      at(2:2)=' '
    endif
    k=d
  endif
enddo

end subroutine

subroutine readpdbna(input)
use comun
implicit none
integer*4 j,kode
character input*256
character ln*6
logical once
open(unit=1,file=input,IOSTAT=kode)
j=0
once = .true.
read(1,'(A)',IOSTAT=kode) ln
do while (kode.eq.0.or.kode.eq.64.and.once)
  if (ln(1:3).eq.'END') then
    na=j
    once = .false.
  elseif ((ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM').and.once) then
    j=j+1
  endif
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
allocate (iatom(na),ires(na),typ(na),res(na),segid(na),resid(na),rt(3,na),w(na),e(na))
end subroutine

subroutine readpdbtnf(input)
use comun
implicit none
integer*4 j,kode
character input*256
character ln*6
open(unit=1,file=input,IOSTAT=kode)
tnf=0
j=0
read(1,'(A)',IOSTAT=kode) ln
do while (kode.eq.0.or.kode.eq.64)
  if (ln(1:3).eq.'END') tnf=tnf+1
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
end subroutine

subroutine readpdb(input)
use comun
implicit none
integer*4 j,kode,kk
character input*256
character ln*256,ann*5,rnn*4,frmt*64
logical once,rnlimit,inlimit
frmt='(A6,I5,x,A5,A5,I4,4x,3F8.3,2F6.2,6x,A4)'
open(unit=1,file=input,IOSTAT=kode)
j=0
!i=1
once = .true.
rnlimit=.false.
inlimit=.false.
read(1,'(A)',IOSTAT=kode) ln
do while (kode.eq.0.or.kode.eq.64.and.once)
  if (ln(1:3).eq.'END') then
!    j=0
!    i=i+1
    once=.false.
  elseif (ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM') then
    j=j+1
!    if (i.eq.1) then
      read (ln,'(6x,A5,x,A4,x,A4,x,A4,4x,3F8.3,2F6.2,6x,A4)') ann,typ(j),res(j),rnn,rt(1,j),rt(2,j),rt(3,j),e(j),w(j),segid(j)
      typ(j)=adjustl(typ(j))
      res(j)=adjustl(res(j))
      iatom(j)=j
!      if (j.eq.1) then 
!        read(ann,'(I5)') iatom(j)
!      else 
!        if (iatom(j-1).lt.99999) then 
!          read(ann,'(I5)') iatom(j)
!        else
!          read(ann,'(Z5)') iatom(j)
!        endif
!      endif
      if (.not.rnlimit.and..not.inlimit) then 
        read(rnn,'(I4)',IOSTAT=kk) ires(j)
        if (ires(j).eq.9999) inlimit=.true.
      elseif (inlimit) then
        read(rnn,'(I4)',IOSTAT=kk) ires(j)
        if (ires(j).ne.9999) then 
          inlimit=.false.
          rnlimit=.true.
          read(rnn,'(Z4)',IOSTAT=kk) ires(j)
        endif
      elseif (rnlimit) then
        read(rnn,'(Z4)',IOSTAT=kk) ires(j)
      endif
      write(resid(j),'(I4)') ires(j)
!    else
!      read (ln,'(30x,3F8.3)') r(1,j,i),r(2,j,i),r(3,j,i)
!    endif
  endif
  read(1,'(A)',IOSTAT=kode) ln
enddo
close(1)
end subroutine

subroutine openpdbcoor(crdfile,u)
use comun
implicit none
integer u
character crdfile*256

open(unit=u,file=crdfile)
nsc=0
dcdopen=.true.
end subroutine

subroutine readpdbcoor(unidad)
use comun
implicit none
integer*4 i,kode,unidad
character ln*128
logical once
i=0
once = .true.
kode=0
do while (kode.eq.0.and.once)
  read(unidad,'(A)',IOSTAT=kode) ln
  if(kode.eq.0) then
    if (ln(1:3).eq.'END') then
      nsc=nsc+1
      once=.false.
    elseif (ln(1:4).eq.'ATOM'.or.ln(4:6).eq.'ATM') then
      i=i+1
      read (ln,'(30x,3F8.3)') rt(1,i),rt(2,i),rt(3,i)
    endif
  endif
enddo

if (kode.ne.0) then
  close(unidad)
  dcdopen=.false.
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

prname='COARSETRJ'
prver='version 1.0'
prdesc='Convert atomistic DNA to coarse grained DNA'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='18 Apr 2013'
lastdate='18 Apr 2013'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
end subroutine

subroutine writedcdhead(dcdfile,un)
use comun
implicit none
integer*4 un
character dcdfile*(*)
integer*4 icntrl(20),itemp,ntitle,natoms
character hdr*4
character*1,allocatable :: title(:)
open(unit=un,file=trim(dcdfile),form='unformatted')
if (nion.gt.0) then
  natoms=gfn*agfpn+csolf(nion)
else
  natoms=gfn*agfpn
endif
icntrl=0
hdr='CORD'
icntrl(10)=1026003171
icntrl(11)=1
icntrl(20)=24
itemp=2
ntitle=icntrl(20)/12*80
allocate (title(ntitle))
call assgn(ntitle,title,'REMARKS FILENAME=coarse.dcd CREATED BY COARSETRJ v1.0 REMARKS DATE: 2013/04/18 CREATED BY: Pablo M. De Biase')
icntrl(1)=tnf
icntrl(2)=1
icntrl(3)=1
icntrl(4)=tnf
write(un) hdr,icntrl
write(un) itemp,title
write(un) natoms
end subroutine

subroutine writedcdbody(un)
use comun
implicit none
integer i,j,k,l,nnn,un,nq,nqq
real*8,allocatable :: rq(:,:),rqq(:,:)
real*8 cnt(3),dv(3)
real*4,allocatable :: x(:,:)
real*8 :: xtlabc6(6)

nnn=gfn*agfpn
do i=1,nion
  nnn=nnn+atlsf(rtf(ions(i)))-atlsi(rti(ions(i)))+1
enddo
allocate (rq(3,na),rqq(3,27*na))

nq=0
cnt=0d0
do i=1,gfn
  do j=1,agfpn
    nq=nq+1
    rq(1:3,nq)=cent(1:3,i,j)
    cnt(1:3)=cnt(1:3)+rq(1:3,nq)
  enddo
enddo
cnt=cnt/nq
do i=1,nion
  do j=atlsi(rti(ions(i))),atlsf(rtf(ions(i)))
    nq=nq+1
    rq(1:3,nq)=rt(1:3,j)
  enddo
enddo

if (center) then
  do i=1,nq
    rq(1:3,i)=rq(1:3,i)-cnt(1:3)
  enddo
endif

if (bvon) then
  do i=1,nq
    rq(1:3,i)=pbc(rq(1:3,i))
  enddo
endif
if (cnvbox) then
  nqq=0
  do i=-1,1,1
    do j=-1,1,1
      do k=-1,1,1
        do l=1,nq
          dv=rq(1:3,l)+i*bv(1:3,1)+j*bv(1:3,2)+k*bv(1:3,3)
          if (inbox(dv)) then
            nqq=nqq+1
            rqq(1:3,nqq)=dv(1:3)
          endif
        enddo
      enddo
    enddo
  enddo
else
  nqq=nq
  rqq(1:3,1:nqq)=rq(1:3,1:nq)
endif

allocate (x(3,nqq))
xtlabc6=0d0
x=rqq(1:3,1:nqq)
! write 
write(un) xtlabc6
write(un) (x(1,i),i=1,nqq)
write(un) (x(2,i),i=1,nqq)
write(un) (x(3,i),i=1,nqq)

deallocate (rq,rqq,x)

end subroutine

subroutine assgn(m,var,text)
implicit none
integer i,n,m
character*1 var(m)
character text*(*)
n=len_trim(text)
if (n.gt.m) stop 'error in string length'
do i=1,n
  var(i)=text(i:i)
enddo
end subroutine
