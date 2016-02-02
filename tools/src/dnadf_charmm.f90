!    DNADF - Computes the RDF and CDF for DNA fragments and free particles.
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
integer nn,nion,tion,ionpairs
integer,allocatable :: resls(:),atlsf(:),atlsi(:)
real*8,allocatable :: rc(:,:),w(:),g(:,:,:),gt(:),gtt(:),gr(:,:,:),gtr(:),gttr(:),gions(:,:),gdna(:,:)
integer,allocatable ::  iatom(:),ires(:),ions(:),csoli(:),csolf(:)
character*4,allocatable :: typ(:),res(:),segid(:),resid(:)
real*8,allocatable :: rt(:,:)
real*8,allocatable :: cent(:,:,:),csol(:,:),boxvol(:),bvall(:,:,:)
real*8 :: bv(3,3),invbv(3,3)
real*8 :: delr,dnavol
character*4,allocatable :: fel(:),fl(:),gfl(:)
integer,allocatable :: gfi(:),gff(:),fi(:),ff(:),fg(:),fp(:)
integer gfn,fn,efn
integer maxbin
integer na,nsc
integer rtn,dp(4),soln
integer,allocatable :: rti(:),rtf(:),rlrt(:),solute(:)
integer agfpn,snd,afrn,frgn
integer, allocatable :: agfpl(:),grf(:,:),frg(:),frgg(:),pmd(:),pmdi(:,:),pmdf(:,:),sndi(:),sndf(:),afr(:),afri(:,:),afrf(:,:)
character*4,allocatable :: rtc(:)
logical dcdopen,depablo
logical termon,charmm
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
  inint=iint(num+0.5d0)
  end function

  function iint(num)
  implicit none
  integer iint
  real*8 num
  if (num.ge.0.0d0) then
    iint=int(num)
  else
    if ((num-int(num)).eq.0.0d0) then
      iint=int(num)
    else
      iint=int(num)-1
    endif
  endif
  end function

  ! convert to lowercase
  function lcase(inchar)
  implicit none
  integer i
  integer*1 s
  character ( len = * ) inchar
  character ( len = len_trim(inchar) ) lcase
  lcase=trim(inchar)
  do i=1,len(lcase)
    s=iachar(lcase(i:i))
    if (s.ge.65.and.s.le.90) lcase(i:i)=achar(s+32)
  enddo
  end function

end module

program dnadf
use comun
implicit none
real*8 aa,bb,vvv,rw
integer h,i,j,k,l,fra,tnf
character crdfile*256,dcdfile*256,mskfile*256,line*256,exten*64
integer ul(256),ll(256),num,narg,arg
logical wout,cdv,ext,docdf,dordf,dordfion
character bs*8

! DNA Volume Correction Factor
! Coefficient obtained for 14 Poly THY from Volume computed using Gaussian 09 PM3 and from Volume computed using ARVO algorithm
! 5582.0727567655500177 PM3
! 3178.3411032982212 ARVO
! volcorfac=1.7562849849479445888d0 
rw=0.6827609835262d0  ! Water VDW Radii considered that gives the PM3 volume

call header
arg=0
narg=COMMAND_ARGUMENT_COUNT()
bs=repeat(achar(8),len(bs))

call readarg('Input CHARMM coordinates (.crd/.cor) or PDB (.pdb) filename: ',narg,arg,crdfile)
exten=lcase(crdfile(index(crdfile,'.',back=.true.)+1:len_trim(crdfile)))

! Read crd
if (exten.eq.'pdb') then
  call readpdbna(crdfile,na)
  allocate (iatom(na),ires(na),typ(na),res(na),segid(na),resid(na),rc(3,na),w(na))
  call readpdb(crdfile)
else
  call readcharmmcrd(crdfile)
endif

call readarg('Input CHARMM trajectory dcd or dcde filename: ',narg,arg,dcdfile)
exten=lcase(dcdfile(index(dcdfile,'.',back=.true.)+1:len_trim(dcdfile)))

if (exten.eq.'dcde') then
  ext=.true.
else
  ext=.false.
endif

call readarg('Compute CDF (y/n) [n]? ',narg,arg,line)
if (line(1:1).eq.'y'.or.line(1:1).eq.'Y') then
  docdf=.true.
else
  docdf=.false.
endif

call readarg('Compute RDF for DNA-IONS (y/n) [y]? ',narg,arg,line)
if (line(1:1).eq.'n'.or.line(1:1).eq.'N') then
  dordf=.false.
else
  dordf=.true.
endif

call readarg('Compute RDF for IONS (y/n) [y]? ',narg,arg,line)
if (line(1:1).eq.'n'.or.line(1:1).eq.'N') then
  dordfion=.false.
else
  dordfion=.true.
endif


if (dordf.or.docdf) then
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
endif

call readarg('Input resolution factor in Angstroms [0.1]: ',narg,arg,line)
if (len_trim(line).eq.0) then
  delr=0.1d0
else
  read(line,*) delr
endif

! Build residue list
call reslist()

! Buid residue type list
call restyplist()

if (dordf.or.docdf) then 
  !Read mask file
  call readmask(mskfile)
  
  write(*,'(A,I0)') 'Read residues: ',nn
  
  call readarg('Define center of coarse particles at geometric center (y/n) [n] ? ',narg,arg,line)
  if (.not.(line(1:1).eq.'y'.or.line(1:1).eq.'Y')) then
    depablo=.true.
    write(*,'(A)') 'Turning on de Pablo et al coarse grained DNA'
    write(*,'(A)') 'Select the atomtype where to center the base (usually N1 for G/A & N3 for C/T)'
    do j=1,4
      k=ff(j)-fi(j)+1
      write(line,'(A,I0,A)') '(4x,',k,'(I4,x))' 
      write(*,line) (i,i=1,k) 
      write(*,*) fl(j),' ',(fel(i),' ',i=fi(j),ff(j))
      read(*,*) dp(j)
    enddo
  else
    depablo=.false.
  endif
endif 

! Inform residue types read and enquire more info'
write(*,'(A)') 'Residue types read: '
write(*,'(A)') 'Number  Symbol Residues'
do i=1,rtn
  write(*,'(I5,3x,A4,3x,I5)') i,rtc(i),rtf(i)-rti(i)+1
enddo

call readarg('Residue Type Number sequence for ions/solvent: ',narg,arg,line)
call findparm(line,256,num,ll,ul)
nion=num
tion=0
allocate (ions(nion))
do i=1,nion
  read(line(ll(i):ul(i)),*) ions(i)
  tion=tion+rtf(ions(i))-rti(ions(i))+1
enddo
allocate (csol(1:3,tion),csoli(nion),csolf(nion))
  
if (dordf.or.docdf) then 
  write(*,'(A)') 'Relevant for normalization: '
  call readarg('Substract DNA volume to the periodic cell volume (y/n)? [n]: ',narg,arg,line)
  if (line(1:1).eq.'y'.or.line(1:1).eq.'Y') then
    cdv=.true.
  else
    cdv=.false.
  endif
  
  ! Define solute
  if (cdv) then 
    write(*,'(A$)') 'Residue Type Number sequence that correspond to the solute (DNA): '
    line=''
    do while(len_trim(line).eq.0)
      read(*,'(A)') line
    enddo
    call findparm(line,256,num,ll,ul)
    allocate (solute(na))
    soln=0
    l=0
    do i=1,num
      read(line(ll(i):ul(i)),*) h
      do j=rti(h),rtf(h)
        do k=atlsi(j),atlsf(j)
          l=l+1
          solute(l)=k
        enddo
      enddo
    enddo
    soln=l
  endif
  ! Build fragment list and allocate cent var
  call fraglist()
  write(*,'(A)') 'Read fragments:'
    write(*,'(A)') 'Fragment Group Number | Fragment Group Label | Fragment Label | Number of Elements | Found Fragments' 
  do i=1,fn
    write(*,'(I5,3x,A4,3x,A4,3x,I5,3x,I5)') fg(i),gfl(fg(i)),fl(i),ff(i)-fi(i)+1,fp(i)
  enddo
endif
if (docdf) then
  call readarg('Select the fragment group number to be used as a fixed Reference Axis (RA) [ENTER to use each fragment group as a RA] : ',narg,arg,line)
  if (len_trim(line).eq.0) then
    fra=0
  else 
    read(line,*) fra 
  endif
endif

!Read dcd header and first frame
call readdcdhead(dcdfile,tnf)
if (ext) then 
  call readdcdebody()
else
  call readdcdbody()
endif
invbv=invmat(matmul(transpose(bv),bv))
allocate (boxvol(tnf),bvall(3,3,tnf))

! Compute Total Volume
call volume(bv(1:3,1),bv(1:3,2),bv(1:3,3),vvv)
write (*,'(A)') 'Box Vector for first frame: '
write (*,*) bv(1:3,1)
write (*,*) bv(1:3,2)
write (*,*) bv(1:3,3)

write (*,'(/A$)') 'Total Box Volume for first frame: '
write (*,*) vvv

! Compute DNA volume
dnavol=0d0
if (cdv) call dnavolume(rw)

!Determine the smallest side of the box
bb=dsqrt(dot_product(bv(1:3,1),bv(1:3,1)))
do i=2,3
  aa=dsqrt(dot_product(bv(1:3,i),bv(1:3,i)))
  if (aa.lt.bb) bb=aa
enddo

write(exten,*) bb/2+delr
exten=adjustl(exten)
call readarg('Choose max distance limit in Ang to compute distribution functions [ '//trim(exten)//' ]: ',narg,arg,line)
if (len_trim(line).eq.0) then
  maxbin = int(bb/delr/2)+1 ! Estimate the optimum number of bins for the length of the rdf/cdf
else
  read(line,*) aa
  maxbin = int(aa/delr)
endif

if (docdf) allocate (g(maxbin,nion,gfn+frgn),gt(maxbin),gtt(maxbin))
if (dordf) allocate (gr(maxbin,nion,gfn+frgn),gtr(maxbin),gttr(maxbin)) ! allocate g of r 
if (dordfion) then
  ionpairs=nion*(nion+1)/2
  allocate (gions(maxbin,ionpairs))
endif

if (docdf) then
  call readarg('Consider terminal residues for CDF (y/n) [y]? ',narg,arg,line)
  if (line(1:1).eq.'n'.or.line(1:1).eq.'N') then
    termon=.false.
  else
    termon=.true.
  endif
endif

call readarg('Write PDB with all elements (y/n) [n]? ',narg,arg,line)
wout=.false.
if (line(1:1).eq.'y'.or.line(1:1).eq.'Y') wout=.true.

! Print number of frames read
write(*,'(A,I8$)') 'Frame: ',nsc

! initialize g and gions
if (docdf) g=0d0
if (dordf) gr=0d0
if (dordfion) gions=0d0 

if (wout) open(unit=2,file='dnacdf.pdb') ! open writexyz
do while (dcdopen)                       ! open loop for each frame
  call volume(bv(1:3,1),bv(1:3,2),bv(1:3,3),vvv)
  boxvol(nsc)=vvv
  bvall(:,:,nsc)=bv
  if (docdf.or.dordf) call centsolute()  ! Compute geometric center (centroid)
  call centsolvent()                     ! Compute geometric center (centroid)
  if (wout)     call writeout(2)         ! writeout
  if (docdf)    call calccdf(fra)        ! compute cdf
  if (dordf)    call calcrdf()           ! compute rdf
  if (dordfion) call calcrdfions()       ! compute rdf for ions
  if (ext) then 
    call readdcdebody()
  else
    call readdcdbody()
  endif
  invbv=invmat(matmul(transpose(bv),bv))
  write(*,'(A8,I8$)') bs,nsc             ! print frame number
enddo
if (wout) close(2) ! close writeout

! NORMALIZATION by cylinder shell surface/spherical shell volume and by number of frames
if (docdf) call normcdf()
if (dordf) call normrdf()
if (dordfion) call normrdfions()

! Print out result
if (docdf) call printcdf()
if (dordf) call printrdf()
if (dordfion) call printrdfion()
call printvol()

write(*,'(/A)') 'Normal termination of DNACDF'
end program

! program output header
subroutine header()
implicit none
character prname*64,prver*32,prdesc*256,startdate*32,lastdate*32,author*64

prname='DNADF'
prver='version 2.6'
prdesc='Computes the RDF and CDF for DNA fragments and free particles.'
author='Pablo M. De Biase (pablodebiase@gmail.com)'
startdate='18 Apr 2011'
lastdate='30 Apr 2013'

write(*,'(/A)') trim(prname)//'  '//trim(prver)
write(*,'(A)') trim(prdesc)
write(*,'(A)') '      created by '//trim(author)//' since '//trim(startdate)//' and last modified '//trim(lastdate)
write(*,'(/A/)') 'NOTE: Press ENTER to accept the default option in between []'
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

subroutine printcdf()
use comun
implicit none
integer i,j,k
character filename*256

! Print out result cdf between dna and ions
do j=1,nion ! for each ion type
  do k=1,gfn ! for each fragment
    write(filename,'(5A)') 'cdf-',trim(gfl(k)),'-',trim(rtc(ions(j))),'.dat'
    open(unit=1,file=filename)
    do i=1,maxbin
      write(1,*) dble(i)*delr,'   ',g(i,j,k)
    enddo
    close(1)
  enddo
  do k=1,frgn ! for each fragment
    write(filename,'(5A)') 'cdf-',trim(fl(frg(k))),'-',trim(rtc(ions(j))),'.dat'
    open(unit=1,file=filename)
    do i=1,maxbin
      write(1,*) dble(i)*delr,'   ',g(i,j,k+gfn)
    enddo
    close(1)
  enddo
enddo
end subroutine

subroutine printrdf()
use comun
implicit none
integer i,j,k
character filename*256

! Print out result rdf between dna and ions
do j=1,nion ! for each ion type
  do k=1,gfn ! for each fragment
    write(filename,'(5A)') 'rdf-',trim(gfl(k)),'-',trim(rtc(ions(j))),'.dat'
    open(unit=2,file=filename)
    do i=1,maxbin
      write(2,*) dble(i)*delr,'   ',gr(i,j,k)
    enddo
    close(2)
  enddo
  do k=1,frgn ! for each fragment
    write(filename,'(5A)') 'rdf-',trim(fl(frg(k))),'-',trim(rtc(ions(j))),'.dat'
    open(unit=2,file=filename)
    do i=1,maxbin
      write(2,*) dble(i)*delr,'   ',gr(i,j,k+gfn)
    enddo
    close(2)
  enddo
enddo
end subroutine

subroutine printrdfion()
use comun
implicit none
integer i,j,k,h
character filename*256

!print rdf for ions
h=0
do i=1,nion
  h=h+1
  write(filename,'(5A)') 'rdf-',trim(rtc(ions(i))),'-',trim(rtc(ions(i))),'.dat'
  open(unit=1,file=filename)
  do j=1,maxbin
    write(1,*) dble(j)*delr,'   ',gions(j,h)
  enddo
  close(1)
enddo
do i=1,nion
  do j=i+1,nion
    h=h+1
    write(filename,'(5A)') 'rdf-',trim(rtc(ions(i))),'-',trim(rtc(ions(j))),'.dat'
    open(unit=1,file=filename)
    do k=1,maxbin
      write(1,*) dble(k)*delr,'   ',gions(k,h)
    enddo
    close(1)
  enddo
enddo
end subroutine 

subroutine printvol()
use comun
implicit none
integer i,j
character filename*256,line*512

!print boxvolume
filename='boxvolume.dat'
open(unit=1,file=filename)
do i=1,nsc
  write(line,*) i,boxvol(i),boxvol(i)-dnavol,dnavol
  write(1,'(a)') trim(line)
enddo
close(1)

!print vectors
filename='vectors.dat'
open(unit=1,file=filename)
do i=1,nsc
  write(line,*) i,(bvall(1:3,j,i),j=1,3)
  write(1,'(a)') trim(line)
enddo
close(1)

end subroutine

subroutine normcdf()
use comun
implicit none
integer i
real*8 rp,rg
real*8,parameter :: pi = 3.14159265358979323846264d0

! NORMALIZATION by cylinder shell surface and by number of frames
rp=0d0
do i=1,maxbin
  rg=dble(i)*delr
  g(i,1:nion,1:gfn+frgn) = g(i,1:nion,1:gfn+frgn)/(rg**2 - rp**2)/dble(nsc)/pi
  rp=rg
enddo
end subroutine

subroutine normrdf()
use comun
implicit none
integer i
real*8 rp,rg
real*8,parameter :: pi133 = 4d0/3d0*3.14159265358979323846264d0

! NORMALIZATION by cylinder shell surface and by number of frames
rp=0d0
do i=1,maxbin
  rg=dble(i)*delr
  gr(i,1:nion,1:gfn+frgn) = gr(i,1:nion,1:gfn+frgn)/(rg**3 - rp**3)/dble(nsc)/pi133
  rp=rg
enddo
end subroutine

subroutine normrdfions()
use comun
implicit none
integer i
real*8 rp,rg
real*8,parameter :: pi133 = 4d0/3d0*3.14159265358979323846264d0

! NORMALIZATION by cylinder shell surface and by number of frames
rp=0d0
do i=1,maxbin
  rg=dble(i)*delr
  gions(i,1:ionpairs) = gions(i,1:ionpairs)/(rg**3 - rp**3)/dble(nsc)/pi133
  rp=rg
enddo
end subroutine

subroutine calccdf(fra)
use comun
implicit none
integer h,j,k,l,fra,a,b,c,d,nobp,nob

do j=1,nion ! for each ion type
  h=0
  do k=1,gfn
    h=h+1
    gtt=0d0
    nob=0
    do l=1,snd
      a=pmdi(h,l)
      b=pmdf(h,l)
      c=csoli(j) ! first atom of the selected ion type
      d=csolf(j) ! last atom of the selected ion type
      if (fra.lt.1.or.fra.gt.gfn) then
        call cdf(cent(1:3,k,1:agfpn),cent(1:3,k,1:agfpn),agfpn,sndi(l),sndf(l),csol(1:3,c:d),d-c+1,pmd(a:b),b-a+1,nobp)
      else
        call cdf(cent(1:3,fra,1:agfpn),cent(1:3,k,1:agfpn),agfpn,sndi(l),sndf(l),csol(1:3,c:d),d-c+1,pmd(a:b),b-a+1,nobp)
      endif
      nob=nob+nobp
      gtt(1:maxbin)=gtt(1:maxbin)+gt(1:maxbin)
    enddo
    g(1:maxbin,j,h)=g(1:maxbin,j,h)+gtt(1:maxbin)/dble(nob)
  enddo
  do k=1,frgn
    h=h+1
    gtt=0d0
    nob=0
    do l=1,snd
      a=pmdi(h,l)
      b=pmdf(h,l)
      c=csoli(j) ! first atom of the selected ion type
      d=csolf(j) ! last atom of the selected ion type
      if (fra.lt.1.or.fra.gt.gfn) then
        call cdf(cent(1:3,frgg(k),1:agfpn),cent(1:3,frgg(k),1:agfpn),agfpn,sndi(l),sndf(l),csol(1:3,c:d),d-c+1,pmd(a:b),b-a+1,nobp)
      else
        call cdf(cent(1:3,fra,1:agfpn),cent(1:3,frgg(k),1:agfpn),agfpn,sndi(l),sndf(l),csol(1:3,c:d),d-c+1,pmd(a:b),b-a+1,nobp)
      endif
      nob=nob+nobp
      gtt(1:maxbin)=gtt(1:maxbin)+gt(1:maxbin)
    enddo
    g(1:maxbin,j,h)=g(1:maxbin,j,h)+gtt(1:maxbin)/dble(nob)
  enddo
enddo
end subroutine

subroutine calcrdf()
use comun
implicit none
integer h,j,k,l,a,b,c,d,nobp,nob

do j=1,nion ! for each ion type
  h=0
  do k=1,gfn
    h=h+1
    gttr=0d0
    nob=0
    do l=1,snd
      a=pmdi(h,l)
      b=pmdf(h,l)
      c=csoli(j) ! first atom of the selected ion type
      d=csolf(j) ! last atom of the selected ion type
      call rdf(cent(1:3,k,1:agfpn),agfpn,csol(1:3,c:d),d-c+1,pmd(a:b),b-a+1,nobp)
      nob=nob+nobp
      gttr(1:maxbin)=gttr(1:maxbin)+gtr(1:maxbin)
    enddo
    gr(1:maxbin,j,h)=gr(1:maxbin,j,h)+gttr(1:maxbin)/dble(nob)
  enddo
  do k=1,frgn
    h=h+1
    gttr=0d0
    nob=0
    do l=1,snd
      a=pmdi(h,l)
      b=pmdf(h,l)
      c=csoli(j) ! first atom of the selected ion type
      d=csolf(j) ! last atom of the selected ion type
      call rdf(cent(1:3,frgg(k),1:agfpn),agfpn,csol(1:3,c:d),d-c+1,pmd(a:b),b-a+1,nobp)
      nob=nob+nobp
      gttr(1:maxbin)=gttr(1:maxbin)+gtr(1:maxbin)
    enddo
    gr(1:maxbin,j,h)=gr(1:maxbin,j,h)+gttr(1:maxbin)/dble(nob)
  enddo
enddo
end subroutine

subroutine calcrdfions()
use comun
implicit none
integer i,j,a,b,c,d,h

h=0
do i=1,nion
  a=csoli(i)
  b=csolf(i)
  h=h+1
  call rdfions(csol(1:3,a:b),b-a+1,csol(1:3,a:b),b-a+1,.true.,h)
enddo
do i=1,nion
  a=csoli(i)
  b=csolf(i)
  do j=i+1,nion
    c=csoli(j) ! first atom of the selected ion type
    d=csolf(j) ! last atom of the selected ion type
    h=h+1
    call rdfions(csol(1:3,a:b),b-a+1,csol(1:3,c:d),d-c+1,.false.,h)
  enddo
enddo
end subroutine

subroutine rdfions(x1,n1,x2,n2,same,h)
use comun
implicit none
logical same
integer n1,n2,i,j,bin,n,h
real*8 x1(3,n1),x2(3,n2)
real*8 cons
real*8 dv(3),rcd,vvv,ggion(maxbin)

ggion=0d0

if (same) then
  do i=1,n1-1
    do j=i+1,n1
      dv=x1(1:3,j)-x1(1:3,i)
      call pbc(dv)
      rcd =dsqrt(dot_product(dv,dv)) ! Computing the distance between the ion and center of the fragment
      bin = int(rcd/delr) + 1 ! Define the position vector based on the distance
      if(bin.le.maxbin) ggion(bin)=ggion(bin)+2d0 ! add odh in the position bin of the g(r)
    enddo
  enddo
  n=n1-1
else
  do i=1,n1
    do j=1,n2
      dv=x2(1:3,j)-x1(1:3,i)
      call pbc(dv)
      rcd =dsqrt(dot_product(dv,dv)) ! Computing the distance between the ion and center of the fragment
      bin = int(rcd/delr) + 1 ! Define the position vector based on the distance
      if(bin.le.maxbin) ggion(bin)=ggion(bin)+1d0 ! add odh in the position bin of the g(r)
    enddo
  enddo
  n=n2
endif

!cons depends on the box volume so it depends on each frame
call volume(bv(1:3,1),bv(1:3,2),bv(1:3,3),vvv)
cons=(vvv-dnavol)/(dble(n1)*dble(n))
gions(1:maxbin,h)=gions(1:maxbin,h)+ggion(1:maxbin)*cons
end subroutine

subroutine pbc(pos)
use comun
implicit none
real*8 pos(3)
pos=pos - matmul(bv,inintv(matmul(invbv,matmul(pos,bv))))
end subroutine

subroutine writexyz(unitn)
use comun
implicit none
integer i,j,nnn,unitn
character nam(16)*2,ionc(5)*2
character strng*32,ionst*10
strng='H B C N O F P S BrI MnFeCoNiCuZn'
ionst='ClK NaLiMg'
do i=1,31,2
  j=(i+1)/2
  nam(j)=strng(i:i+1)
enddo
do i=1,9,2
  j=(i+1)/2
  ionc(j)=ionst(i:i+1)
enddo
nnn=gfn*agfpn
do i=1,nion
  nnn=nnn+atlsf(rtf(ions(i)))-atlsi(rti(ions(i)))+1
enddo
write(unitn,*) nnn
write(unitn,*)
do i=1,gfn
  do j=1,agfpn
    write(unitn,*) nam(grf(i,j)),cent(1:3,i,j)
  enddo
enddo
do i=1,nion
  do j=atlsi(rti(ions(i))),atlsf(rtf(ions(i)))
    write(unitn,*) ionc(i),rt(1:3,j)
  enddo
enddo
end subroutine

subroutine writeout(unitn)
use comun
implicit none
integer i,j,k,unitn,nawo,c,d,a,b
integer,allocatable :: rnwo(:)
real*8,allocatable :: rwo(:,:)
character,allocatable :: rtwo(:)*5,atwo(:)*5
integer,allocatable :: conwo(:,:)

nawo=gfn*agfpn+csolf(nion)
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
allocate (iatom(na),ires(na),typ(na),res(na),segid(na),resid(na),rc(3,na),w(na))
do i=1,na
  read(1,'(A)') line
  call findparm(line,1024,num,ll,ul)
  read(line(ll(1):ul(1)),*) iatom(i)
  read(line(ll(2):ul(2)),*) ires(i)
  res(i)=line(ll(3):ul(3))
  typ(i)=line(ll(4):ul(4))
  read(line(ll(5):ul(5)),*) rc(1,i)
  read(line(ll(6):ul(6)),*) rc(2,i)
  read(line(ll(7):ul(7)),*) rc(3,i)
  segid(i)=line(ll(8):ul(8))
  resid(i)=line(ll(9):ul(9))
  read(line(ll(10):ul(10)),*) w(i)
enddo
close(1)
end subroutine

subroutine writecharmmcrd(crdfile)
use comun
implicit none
integer i,j
character crdfile*256,frmt*64
frmt='(I10,I10,2x,A4,6X,A4,4x,3F20.10,2X,A4,6X,A4,4x,F20.10)'
!         1         1  ADE       H5T           -42.6356481721        8.6779119410       -0.3240608554  DNAA      1               0.0000000000
open(unit=67,file=crdfile)
write(67,'(I10,A)') na,'  EXT'
do i=1,na
  write(67,frmt) iatom(i),ires(i),res(i),typ(i),(rt(j,i),j=1,3),segid(i),resid(i),w(i)
enddo
close(67)
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

subroutine volume(boxv1,boxv2,boxv3,vvv)
implicit none
real*8 boxv1(3),boxv2(3),boxv3(3),vvv

vvv=dot_product(crossprod(boxv1,boxv2),boxv3)

contains 
! Cross Product: w = u x v
  function crossprod(u,v)
  implicit none
  real*8 u(3),v(3),w(3),crossprod(3)
  w(1)=u(2)*v(3)-u(3)*v(2)
  w(2)=u(3)*v(1)-u(1)*v(3)
  w(3)=u(1)*v(2)-u(2)*v(1)
  crossprod=w
  end function
end subroutine

! Radial Distribution Function
! fcc: fragment centers coordinates
! nc: number of centers
! itc: ion type coordinates
! ni: number of ions
! fltc: fragment list to compute
! nfltc: number of elements in fltc
! ils: inferior limit for each strand
! sls: superior limit for each strand
! nob number of bases considered
subroutine rdf(fcc,nc,itc,ni,fltc,nfltc,nob)
use comun
implicit none
integer nc,ni,nfltc,nob
real*8 fcc(3,nc),itc(3,ni)
integer i,j,bin,fltc(1:nfltc)
real*8 cons
real*8 dv(3),rcd,vvv
real*8,parameter :: pi = 3.14159265358979323846264

gtr=0d0
do i=1,nfltc,1
  do j=1,ni
    dv=itc(1:3,j)-fcc(1:3,fltc(i)) ! Computing vector between ion and center of fragment
    call pbc(dv)
    rcd =dsqrt(dot_product(dv,dv)) ! Computing the distance between the ion and center of the fragment
    bin = int(rcd/delr) + 1 ! Define the position vector based on the distance
    if(bin.le.maxbin) gtr(bin)=gtr(bin)+1d0 ! add odh in the position bin of the g(r)
  enddo
enddo
!cons depends on the box volume so it depends on each frame
if (nfltc.ge.1) then
  nob=nfltc
else
  nob=1
endif
call volume(bv(1:3,1),bv(1:3,2),bv(1:3,3),vvv)
cons=(vvv-dnavol)/dble(ni)
gtr(1:maxbin)=gtr(1:maxbin)*cons
end subroutine

! Cylindrical Distribution Function
! rac: reference axis centers coordinates
! fcc: fragment centers coordinates
! nc: number of centers
! itc: ion type coordinates
! ni: number of ions
! fltc: fragment list to compute
! nfltc: number of elements in fltc
! ils: inferior limit for each strand
! sls: superior limit for each strand
subroutine cdf(rac,fcc,nc,ils,sls,itc,ni,fltc,nfltc,nob)
use comun
implicit none
integer nc,ni,ini,fin,ils,sls,nfltc,nob
real*8 rac(3,nc),fcc(3,nc),itc(3,ni)
integer i,j,k,bin,fltc(1:nfltc)
real*8 cons,vvv
real*8 dv(3),h,odh,ax(3),ca(3),hsc,hic,rcd,him
real*8,parameter :: pi = 3.14159265358979323846264
logical doit

gt=0d0
if (fltc(1).eq.ils) then
  ini=2
else
  ini=1
endif
if (fltc(nfltc).eq.sls) then
  fin=nfltc-1
else
  fin=nfltc
endif
do i=ini,fin,1
  ax=rac(1:3,fltc(i)-1)-rac(1:3,fltc(i)+1)  ! vector between prev and next rac element
  call pbc(ax)
  h=dsqrt(dot_product(ax,ax)) ! Compute the module of ax in h
  ax=ax/h ! normalizing ax
  h=h/2 ! Compute height of the cylinder that is half of the module of ax in h
  odh=1d0/h ! dividing unit by the height of the cylinder
  ca=fcc(1:3,fltc(i))-rac(1:3,fltc(i)+1) ! vector between inferior element of rac and the corresponding element of fcc that is in between
  hic=-dot_product(ca,ax)/2 ! Computing the cylinder inferior height from the center of the fragment
  hsc=h+hic ! Computing the cylinder superior height from the center
  do j=1,ni
    dv=itc(1:3,j)-fcc(1:3,fltc(i)) ! Computing vector between ion and center of fragment
    call pbc(dv)
    him=dot_product(dv,ax) ! dot product between the ref axis and the ion vector that is the projection of the ion vector over ax
    if (him.lt.hsc.and.him.ge.hic) then ! if the him is in between hsc and hic then consider that the ion is inside the cylinder
      rcd =dsqrt(dot_product(dv,dv)-him**2) ! Computing the distance between the ion and the reference axis of the cylinder centered in the fragment
      bin = int(rcd/delr) + 1 ! Define the position vector based on the distance
      if(bin.le.maxbin) gt(bin)=gt(bin)+odh ! add odh in the position bin of the g(r)
    endif
  enddo
enddo
if (termon) then
  do k=1,2
    doit=.false.
    if (k.eq.1.and.fltc(1).eq.ils) then 
      ini=1
      i=ini
      ax=rac(1:3,fltc(i))-rac(1:3,fltc(i)+1)  ! vector between actual and next fcc element
      doit=.true.
    elseif (k.eq.2.and.fltc(nfltc).eq.sls) then
      fin=nfltc
      i=fin
      ax=rac(1:3,fltc(i)-1)-rac(1:3,fltc(i))  ! vector between prev and actual fcc element
      doit=.true.
    endif
    if (doit) then
      call pbc(ax)
      h=dsqrt(dot_product(ax,ax)) ! Compute the module of ax in h that is the height of the cylinder
      ax=ax/h ! normalizing ax
      odh=1d0/h ! dividing unit by the height of the cylinder
      hic=-h/2 ! Computing the cylinder inferior height from the center of the fragment
      hsc=h/2 ! Computing the cylinder superior height from the center
      do j=1,ni
        dv=itc(1:3,j)-fcc(1:3,fltc(i)) ! Computing vector between ion and center of fragment
        call pbc(dv)
        him=dot_product(dv,ax) ! dot product between the ref axis and the ion vector that is the projection of the ion vector over ax
        if (him.lt.hsc.and.him.ge.hic) then ! if the him is in between hsc and hic then consider that the ion is inside the cylinder
          rcd =dsqrt(dot_product(dv,dv)-him**2) ! Computing the distance between the ion and the reference axis of the cylinder centered in the fragment
          bin = int(rcd/delr) + 1 ! Define the position vector based on the distance
          if(bin.le.maxbin) gt(bin)=gt(bin)+odh ! add odh in the position bin of the g(r)
        endif
      enddo
    endif
  enddo
endif
!cons depends on the box volume so it depends on each frame
if (fin.ge.ini) then
  nob=fin-ini+1
else
  nob=1
endif
call volume(bv(1:3,1),bv(1:3,2),bv(1:3,3),vvv)
cons=(vvv-dnavol)/dble(ni)
gt(1:maxbin)=gt(1:maxbin)*cons
end subroutine

subroutine readdcdebody()
use comun
implicit none
integer kode,i
real*8 xtlabc(12)
real*4 rtt(3,na)

read(1,iostat=kode) xtlabc
read(1,iostat=kode) (rtt(1,i),i=1,na)
read(1,iostat=kode) (rtt(2,i),i=1,na)
read(1,iostat=kode) (rtt(3,i),i=1,na)
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
else
  close(1)
  dcdopen=.false.
endif
end subroutine

subroutine readdcdhead(dcdfile,tnf)
use comun
implicit none
integer icntrl(20),itemp,kode,tnf
!integer nfile,npriv,nsavc,nstep,nfree
integer natom
character dcdfile*256,hdr*4
integer ntitle
character*1,allocatable ::  title(:)

open(unit=1,file=dcdfile,form='unformatted')
dcdopen=.true.
read(1) hdr,icntrl
ntitle=icntrl(20)/12*80
if (allocated(title)) deallocate (title)
allocate (title(ntitle))
read(1,iostat=kode) itemp,title
read(1) natom
if (na.ne.natom) stop 'Number of atoms differ between .crd and .dcd'
if (icntrl(2).eq.0) then 
  charmm=.true.
else
  charmm=.false.
endif
!nfile=icntrl(1)
!npriv=icntrl(2)
!nsavc=icntrl(3)
!nstep=icntrl(4)
!if(icntrl(9).gt.0) print *, '# fixed atoms = ',icntrl(9)
!nfree = natom-icntrl(9)
!print *, '# of free atoms = ',nfree
!print *, 'total # atom = ', natom,nstep,nsavc
!nsc = nstep/nsavc

allocate (rt(3,na))
tnf=icntrl(4)/icntrl(3)
write(*,'(A,I0)') 'Total number of frames: ',tnf
nsc=0
end subroutine

subroutine readdcdbody()
use comun
implicit none
integer kode,i
real*8 xtlabc(6)
real*4 rtt(3,na)

if (.not.charmm) read(1,iostat=kode) xtlabc
read(1,iostat=kode) (rtt(1,i),i=1,na)
read(1,iostat=kode) (rtt(2,i),i=1,na)
read(1,iostat=kode) (rtt(3,i),i=1,na)
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
else
  close(1)
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
integer i,j,k
integer :: din(nn)
character*4 :: rtct(nn)
integer,allocatable :: mat(:)
logical sta

allocate (mat(nn**2))
rtn=1
din(rtn)=1
rtct(rtn)=res(1)
mat(rtn+(din(rtn)-1)*nn)=1
do i=2,nn
   sta=.false.
   j=rtn
   do while (j.gt.0.and..not.sta)
     if (res(atlsi(i)).eq.rtct(j)) then
       sta=.true.
       din(j)=din(j)+1
       mat(j+(din(j)-1)*nn)=i
     else
       j=j-1
     endif
   enddo
   if (.not.sta) then
     rtn=rtn+1
     din(rtn)=1
     rtct(rtn)=res(atlsi(i))
     mat(rtn+(din(rtn)-1)*nn)=i
   endif
enddo

allocate (rtf(rtn),rtc(rtn),rti(rtn)) 
rtc(1:rtn)=rtct(1:rtn)
k=0
do i=1,rtn
  rti(i)=k+1
  do j=1,din(i)
    k=k+1
    rlrt(k)=mat(i+(j-1)*nn)
  enddo
  rtf(i)=k   
enddo
deallocate (mat)
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
c(1)=c(1)/dble(n)
c(2)=c(2)/dble(n)
c(3)=c(3)/dble(n)

end subroutine

subroutine centsolute()
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
    call centroid(rt(1:3,afr(c:d)),d-c+1,cent(1:3,i,j))
  enddo
endif
do i=2,gfn                ! residue number per residue type
  do j=1,agfpn ! fragment number
    c=afri(i,j)
    d=afrf(i,j)
    call centroid(rt(1:3,afr(c:d)),d-c+1,cent(1:3,i,j))
  enddo
enddo
end subroutine

subroutine centsolvent()
use comun
implicit none
integer i,j,k,c,d

k=0
do i=1,nion ! each solvent/ion type
  csoli(i)=k+1
  do j=rti(ions(i)),rtf(ions(i)) ! all ions/solvent for each ion/solvent type
    k=k+1
    c=atlsi(j)
    d=atlsf(j)
    call centroid(rt(1:3,c:d),d-c+1,csol(1:3,k))
  enddo
  csolf(i)=k
enddo

end subroutine

! Uses ARVO Algorithm
! The authors ask to acknowledge the use of ARVO in all resulting publications (by referring to [Jan Busa, Jozef Dzurina,
! Edik Hayryan, Shura Hayryan, Chin-Kun Hu, Jan Plavka, Imrich Pokorny, Jaroslav Skrivanek, and Ming-Chya Wu,
! ARVO: A Fortran package for computing the solvent accessible surface area and the excluded volume of overlapping
! spheres via analytic equations, Comp. Phys. Comm. (Computer Physics Communications) 165 (2005) 59-96)
subroutine dnavolume(rw)
use comun
implicit none
integer i,j,ns
real*8, parameter :: pi=3.14159265358979323846264d0
integer, parameter :: ks=1000,kl=1000,ka=2000,ki=100000
character at*2
! ks = soln
!                ks - maximal spheres' number
!                kl - maximal neighbors' number of one sphere (local spheres' number)     
!                ka - maximal angles' or arcs' number
!                ki - maximal neighbors' relations' number = cca.
!                     spheres' number * maximal neighbors' number
!        eps_north_pole - accuracy level in the function North_Pole_test
!          eps_deltat - accuracy level in the subroutine circles_intersection
!        eps_angle - accuracy level in the subroutine delete_equal (angles)
          
real*8 spheres(soln,4),av(2),v,a,vdwr,sa,rw
integer neighbors_number(soln),index_start(soln),neighbors_indices(ki)
integer North_Pole_test

do i=1,soln
  j=solute(i)
  call getat(typ(j),4,at) 
  call vdwradius(at,vdwr)
  spheres(i,1)=rt(1,j)
  spheres(i,2)=rt(2,j)
  spheres(i,3)=rt(3,j)
  spheres(i,4)=vdwr+rw
enddo

ns=soln     

!
!     ns - spheres number
!

!     Study the neighborhood relations
call make_neighbors(1,ns,spheres,neighbors_number,index_start,neighbors_indices,soln,kl,ns,ki)
        
!       If some North Pole is close to the other atom's surface
!     molecule's rotation is necessary
do while (North_Pole_test(1,ns,spheres,neighbors_number,index_start,neighbors_indices,soln,ki).EQ.0)
!  print *,'Rotation after bad North Pole test!'
!         write(10,*)'Rotation after bad North Pole test!'
  sa=0.324d0 ! "Random" sin value
  call spheres_rotation(spheres,soln,ns,sa)! random molecule rotation 
enddo

!     Computation of area and volume as a sum of surface integrals
V=0d0
A=0d0   
do i=1,ns
  call areavolume(i,spheres,neighbors_number,index_start,neighbors_indices,soln,kl,ka,ki,av)
  V=V+av(1)
  A=A+av(2)
enddo
print *,'SOLUTE:   Volume= ',V,' Area= ',A,'  Spheres num= ',ns,' Water radius= ',rw
dnavol=v!*volcorfac
!print *,'DNA Volume with correction factor: ',dnavol 
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

subroutine vdwradius(at,vdwr)
implicit none
integer i,an
real*8 vdw(113),vdwr
character at*2

do i=1,112
  vdw(i)=2.0
enddo

! http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
! Extrapolation is done from Calulated data multiplied by the average of VDW Radii divided by average of Calculated Radii

vdw(1)=1.20!	H	
vdw(2)=1.40!	He	
vdw(3)=1.82!	Li	
vdw(4)=1.53!	Be	
vdw(5)=1.92!	B	
vdw(6)=1.70!	C	
vdw(7)=1.55!	N	
vdw(8)=1.52!	O	
vdw(9)=1.47!	F	
vdw(10)=1.54!	Ne	
vdw(11)=2.27!	Na	
vdw(12)=1.73!	Mg	
vdw(13)=1.84!	Al	
vdw(14)=2.10!	Si	
vdw(15)=1.80!	P	
vdw(16)=1.80!	S	
vdw(17)=1.75!	Cl	
vdw(18)=1.88!	Ar	
vdw(19)=2.75!	K	
vdw(20)=2.31!	Ca	
vdw(21)=2.11!	Sc	
vdw(22)=2.18!	Ti	extrapolated
vdw(23)=2.12!	V	extrapolated
vdw(24)=2.06!	Cr	extrapolated
vdw(25)=1.99!	Mn	extrapolated
vdw(26)=1.93!	Fe	extrapolated
vdw(27)=1.88!	Co	extrapolated
vdw(28)=1.63!	Ni	
vdw(29)=1.40!	Cu	
vdw(30)=1.39!	Zn	
vdw(31)=1.87!	Ga	
vdw(32)=2.11!	Ge	
vdw(33)=1.85!	As	
vdw(34)=1.90!	Se	
vdw(35)=1.85!	Br	
vdw(36)=2.02!	Kr	
vdw(37)=3.03!	Rb	
vdw(38)=2.49!	Sr	
vdw(39)=2.63!	Y	extrapolated
vdw(40)=2.55!	Zr	extrapolated
vdw(41)=2.45!	Nb	extrapolated
vdw(42)=2.35!	Mo	extrapolated
vdw(43)=2.27!	Tc	extrapolated
vdw(44)=2.21!	Ru	extrapolated
vdw(45)=2.14!	Rh	extrapolated
vdw(46)=1.63!	Pd	
vdw(47)=1.72!	Ag	
vdw(48)=1.58!	Cd	
vdw(49)=1.93!	In	
vdw(50)=2.17!	Sn	
vdw(51)=2.06!	Sb	
vdw(52)=2.06!	Te	
vdw(53)=1.98!	I	
vdw(54)=2.16!	Xe	
vdw(55)=3.43!	Cs	
vdw(56)=2.68!	Ba	
vdw(59)=3.06!	Pr	extrapolated
vdw(60)=2.55!	Nd	extrapolated
vdw(61)=2.54!	Pm	extrapolated
vdw(62)=2.95!	Sm	extrapolated
vdw(63)=2.86!	Eu	extrapolated
vdw(64)=2.89!	Gd	extrapolated
vdw(65)=2.79!	Tb	extrapolated
vdw(66)=2.82!	Dy	extrapolated
vdw(68)=2.80!	Er	extrapolated
vdw(69)=2.75!	Tm	extrapolated
vdw(70)=2.75!	Yb	extrapolated
vdw(71)=2.69!	Lu	extrapolated
vdw(72)=2.58!	Hf	extrapolated
vdw(73)=2.48!	Ta	extrapolated
vdw(74)=2.39!	W	extrapolated
vdw(75)=2.33!	Re	extrapolated
vdw(76)=2.29!	Os	extrapolated
vdw(77)=2.23!	Ir	extrapolated
vdw(78)=1.75!	Pt	
vdw(79)=1.66!	Au	
vdw(80)=1.55!	Hg	
vdw(81)=1.96!	Tl	
vdw(82)=2.02!	Pb	
vdw(83)=2.07!	Bi	
vdw(84)=1.97!	Po	
vdw(85)=2.02!	At	
vdw(86)=2.20!	Rn	
vdw(87)=3.48!	Fr	
vdw(88)=2.83!	Ra	
vdw(92)=1.86!	U	

call att2atn(at,an)
if (an.eq.0) an=112
vdwr=vdw(an)
end subroutine

subroutine att2atn(atomtype,atomnumber)
implicit none
integer atomnumber,i,l
character atomtype*2,table*114,at1*2,at2*2
logical loopon
!table='H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBa'
table='h helibeb c n o f nenamgalsip s clark casctiv crmnfeconicuzngageassebrkrrbsry zrnbmotcrurhpdagcdinsnsbtei xecsbael'

loopon=.true.
i=1
l=len_trim(atomtype)
call lcase(atomtype(1:l),at1(1:l),l)
if (l.eq.1) at1(2:2)=' '

do while (loopon.and.i.lt.len(table))
  at2=table(i:i+1)
  if (at1.eq.at2) then
    atomnumber=int((i+1)/2)
    loopon=.false.
  endif
  i=i+2
enddo

if (loopon.and.l.eq.2) then
  i=1
  at1(2:2)=' '
  do while (loopon.and.i.lt.len(table))
    at2=table(i:i+1)
    if (at1.eq.at2) then
      atomnumber=int((i+1)/2)
      loopon=.false.
      write(*,*) 'Warning: ',atomtype,' recognized with the atomic number ',atomnumber
    endif
    i=i+2
  enddo
endif

if (loopon) then
  write(*,*) 'ERROR: Atom type ',atomtype,' not identified as an element.'
  atomnumber=0
endif
end subroutine 

! convert a string to all lower case
subroutine lcase(inchar,outchar,cn)
implicit none
integer i,cn
character s*1
character ( len = cn ) inchar,outchar

do i=1,len_trim(inchar)
  s=inchar(i:i)
  if (iachar(s).ge.65.and.iachar(s).le.90) then
    outchar(i:i)=achar(iachar(s)+32)
  else
    outchar(i:i)=inchar(i:i)
  endif
enddo
outchar=outchar(1:len_trim(inchar))
end subroutine

subroutine make_neighbors(i1,i2,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ns,ki)

!     Determination of neighbors for all atoms. We construct next structure:
!          neighbors_number(i)=neighbors number for ith atom
!              index_start(i)=start of neighbors indices for ith atom in array 
!                               neighbors_indices
!              neighbors_indices - array of neighbors indices for each atom
! neighbors_indices(index_start(i)):neighbors_ind(index_start(i)+neighbors_number(i)-1)

!              For example: 1. atom has neighbors with indices 2, 4, 7
!                                       2. atom has neighbors with indices 1, 3
!                                       3. atom has neighbors with indices 2, 4
!                                       4. atom has neighbors with indices 1, 3
!                       5. atom is subset of some atom 
!                       6. atom has no neighbors 
!                                       7. atom has neighbors with index 1
!              then we have
!               neighbors_number=(3,2,2,2,-1,0,1)
!                        index_start=(1,4,6,8,10,10,10,11)
!                       neighbors_indices(2,4,7,1,3,2,4,1,3,1)
!               
implicit none
integer i,i1,i2,ks,kl,ns,ki,j
real*8 spheres(ks,4)
integer ind(kl),neighbors_number(ks),index_start(ks+1),neighbors_indices(ki),neighbors

index_start(i1)=1
do i=i1,i2
  neighbors_number(i)=neighbors(i,spheres,ind,ks,kl,ns)                
  if (neighbors_number(i).LE.0) then 
!           sphere is subset ot there are no neighbors  
    index_start(i+1)=index_start(i)
  else  ! there are neighbors
    index_start(i+1)=index_start(i)+neighbors_number(i)
    do j=1,neighbors_number(i)
       neighbors_indices(index_start(i)+j-1)=ind(j)
    enddo              
  endif
enddo

return
end subroutine

integer function neighbors(i,spheres,ind,ks,kl,ns)
!     
!     Function 

!     If ith sphere is a subset of other sphere, index_number(i)=-1 and we change its 
!     radius in matrix spheres to -radius !!!
!     If some other sphere is subset of ith sphere, than we change its radius to -radius !!!
implicit none
integer ks,kl,ns,i,k,neighbors_num
real*8 spheres(ks,4),xi,yi,ri,dd,rk,zi
integer ind(kl)

neighbors_num=0
!     i-th sphere data
xi=spheres(i,1)
yi=spheres(i,2)
zi=spheres(i,3)
ri=spheres(i,4)
do k=1,ns
  if (k .NE. i) then
!                 first simple test 
    if(dabs(xi-spheres(k,1)).LT.ri+spheres(k,4)) then
      dd=dsqrt((xi-spheres(k,1))**2+(yi-spheres(k,2))**2+(zi-spheres(k,3))**2)
      rk=spheres(k,4)
      if (dd.LT.ri+rk) then
        if (dd+ri.LE.rk) then
!                        ith sphere is inside of other sphere 
          neighbors_num=-1
          exit
        elseif (dd+rk.GT.ri) then  
!                            kth sphere is neighbor 
          neighbors_num=neighbors_num+1
          ind(neighbors_num)=k
        endif
      endif
    endif
  endif
enddo
neighbors=neighbors_num
return
end function

integer function North_Pole_test(i1,i2,spheres,neighbors_number,index_start,neighbors_indices,ks,ki)
!        Here we check, that North Pole of no sphere lies on other 
!        neighbor sphere

!        dmin - square of minimal distance of the North Pole to neighbor sphere surface


implicit none
integer ks,ki,i1,i2,i,ink,npt,k
real*8 spheres(ks,4),eps_north_pole,dmin,d
integer neighbors_number(ks),index_start(ks),neighbors_indices(ki)

!     Test precision - MAY BE CHANGED
eps_north_pole=1d-4

dmin=10000d0
do i=i1,i2  
  do k=1,neighbors_number(i)
    ink=neighbors_indices(index_start(i)+k-1) ! kth neighbor index
    d=dabs(dsqrt((spheres(i,1)-spheres(ink,1))**2+(spheres(i,2)-spheres(ink,2))**2 &
                +(spheres(i,3)+spheres(i,4)-spheres(ink,3))**2)-spheres(ink,4))
    if (d.LT.dmin) dmin=d 
  enddo
enddo

!print *, dmin

!     minimal distance = dmin        
if (dmin.LT.eps_north_pole) then
  npt=0 ! Bad news!
else
  npt=1 ! O.K.
endif
   
North_Pole_test=npt   
      
return
end function

subroutine spheres_rotation(spheres,ks,ns,sa)
!     Random rotation of molecule about the y-axis
!     after bad North Pole test.
!        Some North Pole is near other spheres surface
implicit none
integer ks,ns,i
real*8 spheres(ks,4),sa,ca,x,z

ca=dsqrt(1d0-sa*sa)         
do i=1,ns
  x=spheres(i,1)
  z=spheres(i,3)
  spheres(i,1)=ca*x-sa*z
  spheres(i,3)=sa*x+ca*z                
enddo
return
end subroutine

subroutine areavolume(i,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ka,ki,av)

!   Function computes i-th part of the whole volume -
!    - the volume of domain inside i-th and outside of
!      all other spheres

implicit none
integer ks,kl,ka,ki,i,j,nls,npos
integer neighbors_number(ks),index_start(ks),neighbors_indices(ki)
integer ind(kl),narcs
integer circles_to_arcs
real*8 spheres(ks,4),av(2)
real*8 circles(kl,4),arcs(ka,3),sphere_local(kl,4),avi(2)
real*8 z1,r1
real*8,parameter :: pi=3.14159265358979323846264d0

!     circles, arcs, sphere_local are described below



!      Determination of i-th sphere's neighbors (row indices in matrix spheres)

if (neighbors_number(i).LT.0) then

!     ith sphere is subset of other sphere
!         sphere(i,4) will be done negative

  av(1)=0d0
  av(2)=0d0
elseif (neighbors_number(i).EQ.0) then 

!     there are no neighbors (nls - number of local spheres = ith sphere + neighbors)

  av(1)=4d0*pi*spheres(i,4)**3/3.d0
  av(2)=4d0*pi*spheres(i,4)**2
else 

!        there are neighbors

  nls=neighbors_number(i)+1
  ind(1)=i
  do j=1,(nls-1)
    ind(j+1)=neighbors_indices(index_start(i)+j-1)
  enddo

!        we will work only with ith and neighbors spheres                            
  call local_spheres(spheres,ind,sphere_local,nls,ks,kl)
  av(1)=0d0
  av(2)=0d0

  call make_ts_circles(sphere_local,circles,kl,nls)
  narcs=circles_to_arcs(circles,arcs,kl,nls,ka)

  npos=0
  do j=1,(nls-1)
    if (circles(j,4).GT.0) npos=npos+1
  enddo

  z1=sphere_local(1,3)
  r1=sphere_local(1,4)
  if (npos.GT.0) then   
!           there exists positive oriented circle 
    call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi) 
    av(1)=av(1)+avi(1)
    av(2)=av(2)+avi(2)
  else 
!           all circles are negative oriented - we compute complement
    call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi) 
    av(1)=av(1)+avi(1)+4d0*pi*sphere_local(1,4)**3/3d0 
    av(2)=av(2)+avi(2)+4d0*pi*sphere_local(1,4)**2 
  endif
endif
return
end subroutine

subroutine local_spheres(spheres,ind,sphere_local,nls,ks,kl)

!     Take sphere_local out of the main array spheres

implicit none
integer ks,nls,kl,i,j
real*8 spheres(ks,4),sphere_local(kl,4)
integer ind(kl)

do i=1,nls
  do j=1,4
    sphere_local(i,j)=spheres(ind(i),j)
  enddo 
enddo 
return
end subroutine

subroutine make_ts_circles(sphere_local,circles,kl,nls)

!     Preparing circles structure for 1st sphere in array        circles
!     according to the paper Hayrjan, Dzurina, Plavka, Busa
!     
!     circles(i,1)=ti 
!     circles(i,2)=si     - ith circle's center coordinates
!     circles(i,3)=ri    - ith circle's radius 
!     circles(i,4)=+1/-1 - circle orientation 

implicit none
integer kl,nls,k
real*8 a,b,c,d,r1,dx,dy,circles(kl,4),sphere_local(kl,4)

r1=sphere_local(1,4)
do k=1,(nls-1)
  dx=sphere_local(1,1)-sphere_local(k+1,1)
  dy=sphere_local(1,2)-sphere_local(k+1,2)        
  a=dx*dx+dy*dy+(sphere_local(1,3)+r1-sphere_local(k+1,3))**2-sphere_local(k+1,4)**2
  b=8d0*r1*r1*dx
  c=8d0*r1*r1*dy
  d=4d0*r1*r1*(dx*dx+dy*dy+(sphere_local(1,3)-r1-sphere_local(k+1,3))**2-sphere_local(k+1,4)**2)
  circles(k,1)=-b/(2d0*a)
  circles(k,2)=-c/(2d0*a)       
  circles(k,3)=dsqrt((b*b+c*c-4d0*a*d)/(4d0*a*a))   
  if (a.GT.0) then
    circles(k,4)=-1d0
  else 
    circles(k,4)=1d0
  endif
enddo
return    
end subroutine

integer function circles_to_arcs(circles,arcs,kl,nls,ka)

!     Computing integration arcs
!     
!     arcs(i,1)=ci            - corresponding circle index
!     arcs(i,2)=sigma     - starting arc angle 
!     arcs(i,3)=delta     - oriented arc angle

!     Arcs (with their orientation) are parts of circles, which
!     bounds are circles intersection points. If the center of
!     arc lies inside all other positive and outside all other
!     negative circles, then we will put it inside arcs structure

implicit none
integer kl,ka,nls,number_arc,nna,i,j,k,new_arcs
real*8 arcs(ka,3),circles(kl,4),arcsnew(ka,3)
real*8, parameter :: pi=3.14159265358979323846264d0

number_arc=0
if (nls.EQ.2) then
!         we have only 1 circle
  number_arc=1
  arcs(1,1)=1d0
  arcs(1,2)=0d0
  arcs(1,3)=2d0*pi*circles(1,4)
else 
!         more than 1 circle
  do i=1,(nls-1)
    nna=new_arcs(i,circles,arcsnew,kl,ka,nls)
    if (nna.GT.0) then
      do j=1,nna
        do k=1,3
          arcs(number_arc+j,k)=arcsnew(j,k)
        enddo                          
      enddo
      number_arc=number_arc+nna
    endif    
  enddo
endif

circles_to_arcs=number_arc
return
end function

integer function new_arcs(i,circles,arcsnew,kl,ka,nls)

!     Function prepares arcs, which are part of i-th circle
!     in circle structure circles.
!     Interesting are these arcs, which are inside other positive
!     circles or outside other negative circles

!     Matrix arcsnew in each row has elements

!     arcsnew(i,1)=ic - ic is the index of arc-circle in circle 
!     arcsnew(i,2)=sigma - sigma is the starting angle of arc
!     arcsnew(i,3)=delta - delta is oriented arc angle

implicit none
integer ka,kl,nls,i,num_arc,num_angle,number_cond,na,jj,j
real*8 circles(kl,4),arcsnew(ka,3),angles(ka),a1,a2,t,ti,ri,b1,b2,d,r,s,si
real*8,parameter :: pi=3.14159265358979323846264d0
integer circle_in_circle,point_in_circle,delete_equal

num_arc=0
num_angle=0

ti=circles(i,1)
si=circles(i,2)
ri=circles(i,3)
do j=1,(nls-1) 
!        composition of angles vector, consisting of intersection points
  if (j .NE. i) then
    t=circles(j,1)
    s=circles(j,2)
    r=circles(j,3)
    d=dsqrt((ti-t)**2+(si-s)**2)
    if ( (d.LT.r+ri) .AND. (dabs(r-ri).LT.d) ) then
!                     2 intersection points exist
      call circles_intersection(i,j,circles,kl,a1,a2,b1,b2)
      angles(num_angle+1)=a1
      angles(num_angle+2)=a2
      num_angle=num_angle+2            
    endif
  endif
enddo
if (num_angle .EQ. 0) then
!         there are no double intersections of i-th circles with others
  number_cond=0 
!             if i-th circle is inside of all other positive and outside of 
!         all other negative circles, it will be new arc
  do j=1,(nls-1)
    if (j.NE.i) number_cond=number_cond+circle_in_circle(i,j,circles,kl)
  enddo
  if (number_cond.EQ.(nls-2)) then
!                 all conditions hold
    num_arc=1
    arcsnew(1,1)=i
    arcsnew(1,2)=0d0
    arcsnew(1,3)=2d0*pi*circles(i,4)
  endif
else 
!         there are double intersection points
  if (circles(i,4).GT.0) then
    call mysort(angles,ka,num_angle)
  else
    call mydsort(angles,ka,num_angle)
  endif
  na=delete_equal(angles,ka,num_angle)
  num_angle=na
  do j=1,(na-1)
    number_cond=0
    do jj=1,(nls-1)
      if (jj.NE.i) then
        t=ti+ri*dcos((angles(j)+angles(j+1))/2d0)
        s=si+ri*dsin((angles(j)+angles(j+1))/2d0)
        number_cond=number_cond+point_in_circle(t,s,jj,circles,kl)
      endif
    enddo
    if (number_cond.EQ.(nls-2)) then
!                         all conditions hold
      num_arc=num_arc+1  
      arcsnew(num_arc,1)=i
      arcsnew(num_arc,2)=angles(j)
      arcsnew(num_arc,3)=angles(j+1)-angles(j)
    endif
  enddo
  number_cond=0
  do j=1,(nls-1)
    if (j.NE.i) then
      t=ti+ri*dcos((angles(1)+2d0*pi+angles(na))/2d0)
      s=si+ri*dsin((angles(1)+2d0*pi+angles(na))/2d0)
      number_cond=number_cond+point_in_circle(t,s,j,circles,kl)
    endif
  enddo
  if (number_cond.EQ.(nls-2)) then
!                     all conditions hold
    num_arc=num_arc+1  
    arcsnew(num_arc,1)=i
    arcsnew(num_arc,2)=angles(na)
    arcsnew(num_arc,3)=angles(1)+circles(i,4)*2d0*pi-angles(na)
  endif
endif

new_arcs=num_arc
return
end function

subroutine circles_intersection(ic1,ic2,circles,kl,a1,a2,b1,b2)

!     Function returns angles of two intersection points
!     of circles with indices ic1 and ic2 in circles structure circles
!     (we will use it ONLY IN CASE, WHEN 2 INTERSECTION POINTS EXIST!!!)

!     a1 and a2 are corresponding angles with respect to the center of 1st circle
!     b1 and b2 are corresponding angles with respect to the center of 2nd circle

implicit none
integer ic1,ic2,kl
real*8 circles(kl,4),eps_deltat,t1,s1,r1,t2,s2,r2,a,b1,b2,b,c,d,a1,a2
real*8,parameter :: pi=3.14159265358979323846264d0

eps_deltat=1d-12
!     (t,s) - circle center, r - circle radius
t1=circles(ic1,1)
s1=circles(ic1,2)
r1=circles(ic1,3) 
t2=circles(ic2,1)
s2=circles(ic2,2)
r2=circles(ic2,3) 
if (dabs(t2-t1).LT.eps_deltat) then
!         t2 .EQ. t1
  B=((r1*r1-r2*r2)/(s2-s1)-(s2-s1))/2d0
  A=dsqrt(r2*r2-B*B)
  if (B.EQ.0) then
    b1=0d0
    b2=pi
  elseif (B.GT.0) then
    b1=datan(dabs(B/A))
    b2=pi-b1
  else
    b1=pi+datan(dabs(B/A))
    b2=3d0*pi-b1        
  endif
  B=B+s2-s1
  if (B.EQ.0) then
    a1=0d0
    a2=pi
  elseif (B.GT.0) then
    a1=datan(dabs(B/A))
    a2=pi-a1
  else
    a1=pi+datan(dabs(B/A))
    a2=3d0*pi-a1        
  endif    
else 
!       t2 .NE. t1
  C=((r1*r1-r2*r2-(s2-s1)**2)/(t2-t1)-(t2-t1))/2d0
  D=(s1-s2)/(t2-t1)
  B=(-C*D+dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
  A=C+D*B
  if (A.EQ.0) then
    if (B.GT.0) then
      b1=pi/2d0
    else 
      b1=-pi/2d0
    endif
  elseif (A.GT.0) then
    b1=datan(B/A)
  else
    b1=pi+datan(B/A)
  endif
  B=B+s2-s1
  A=A+t2-t1
  if (A.EQ.0) then 
    if (B.GT.0) then
      a1=pi/2d0
    else 
      a1=-pi/2d0
    endif
  elseif (A.GT.0) then
    a1=datan(B/A)
  else
    a1=pi+datan(B/A)
  endif
  B=(-C*D-dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
  A=C+D*B
  if (A.EQ.0) then
    if (B.GT.0) then
      b2=pi/2d0
    else 
      b2=-pi/2d0
    endif
  elseif (A.GT.0) then
    b2=datan(B/A)
  else
    b2=pi+datan(B/A)
  endif
  B=B+s2-s1
  A=A+t2-t1
  if (A.EQ.0) then
    if (B.GT.0) then
      a2=pi/2d0
    else 
      a2=-pi/2d0
    endif
  elseif (A.GT.0) then
    a2=datan(B/A)
  else
    a2=pi+datan(B/A)
  endif
endif
if (a1.LT.0) a1=a1+2d0*pi 
if (a2.LT.0) a2=a2+2d0*pi
if (b1.LT.0) b1=b1+2d0*pi
if (b2.LT.0) b2=b2+2d0*pi
return
end subroutine

integer function circle_in_circle(i,k,circles,kl)

!      1  if i-th circle is inside k-th positive circle or 
!                          outside k-th negative circle
!      0  - otherwise

!      WE KNOW, THAT CIRCLES HAVE LESS THAN 2 INTERSECTION POINTS !!!

implicit none
integer i,k,kl
real*8 circles(kl,4),d

d=dsqrt((circles(i,1)+circles(i,3)-circles(k,1))**2+(circles(i,2)-circles(k,2))**2)
if (d.LT.circles(k,3)) then
  if (circles(k,4).GT.0) then
    circle_in_circle=1
  else
    circle_in_circle=0
  endif    
elseif (d.GT.circles(k,3)) then 
  if (circles(k,4).GT.0) then
    circle_in_circle=0
  else
    circle_in_circle=1
  endif
else 
! d=circles(k,3) - right point on k-th circle - touching of circles
  d=dsqrt((circles(i,1)-circles(k,1))**2+(circles(i,2)-circles(k,2))**2)
  if (d.LT.circles(k,3)) then
    if (circles(k,4).GT.0) then
      circle_in_circle=1
    else
      circle_in_circle=0
    endif
  else
    if (circles(k,4).GT.0) then
      circle_in_circle=0
    else
      circle_in_circle=1
    endif
  endif
endif
return
end function

integer function point_in_circle(t,s,k,circles,kl)


!     1  if point (t,s) is inside k-th positive circle 
!                      or outside k-th negative circle
!     0  - otherwise

!     WE KNOW, THAT POINT IS NOT ON THE CIRCLE !!!

implicit none
integer kl,k
real*8 circles(kl,4),t,s,d

d=dsqrt((t-circles(k,1))**2+(s-circles(k,2))**2)
if (d.LT.circles(k,3)) then
  if (circles(k,4).GT.0) then
    point_in_circle=1
  else
    point_in_circle=0
  endif
else
  if (circles(k,4).GT.0) then
    point_in_circle=0
  else
    point_in_circle=1
  endif
endif
return
end function

subroutine mysort(angles,ka,num_angle)

!     Sorting array angles in increasing order
!     num_angle is the angles array length

implicit none
integer ka,num_angle,ii,i,j
real*8 angles(ka),amax

do i=1,(num_angle-1)
  ii=i
  amax=angles(i)
  do j=i+1,num_angle
    if (amax.GT.angles(j)) then
      ii=j
      amax=angles(j)                
    endif
  enddo
  if (ii.NE.i) then
    angles(ii)=angles(i)
    angles(i)=amax
  endif
enddo
return
end subroutine

subroutine mydsort(angles,ka,num_angle)

!     Sorting array angles in decreasing order
!     num_angle is the angles array length

implicit none
integer ka,num_angle,ii,i,j
real*8 angles(ka),amin

do i=1,(num_angle-1)
  ii=i
  amin=angles(i)
  do j=i+1,num_angle
    if (amin.LT.angles(j)) then
      ii=j
      amin=angles(j)                
    endif
  enddo
  if (ii.NE.i) then
    angles(ii)=angles(i)
    angles(i)=amin
  endif
enddo
return
end subroutine

integer function delete_equal(angles,ka,num_angle)

!     Deletion of "equal" (to some precision eps_angle)
!     angles in sorted vector angles

implicit none
integer ka,num_angle,m,i
real*8 angles(ka),anglesnew(ka),eps_angle,angle

eps_angle=1d-12
m=1 
angle=angles(1)
anglesnew(1)=angle
do i=2,num_angle
  if (dabs(angles(i)-angle).GT.eps_angle) then
    angle=angles(i)
    m=m+1
    anglesnew(m)=angle
  endif
enddo
delete_equal=m
do i=1,m
  angles(i)=anglesnew(i)
enddo

return
end function

subroutine avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)

!     Computing integrals over arcs given in arc structure
!     according to paper Hayrian, Dzurina, Plavka, Busa
!     
implicit none
integer kl,ka,narcs,k
real*8 circles(kl,4),arcs(ka,3),avi(2),eps_two_pi,t,s,r,a,b,c,rr,vIone,vItwo,vIthree,vJone,vJtwo,vJthree
real*8 delta_vint,delta_aint,r1,z1,al,be,sb,cb,sa,ca,fract
real*8,parameter :: pi=3.14159265358979323846264d0

eps_two_pi=1d-12
avi(1)=0d0
avi(2)=0d0

do k=1,narcs 
!             cycle over all arcs
  t=circles(int(arcs(k,1)),1)
  s=circles(int(arcs(k,1)),2)
  r=circles(int(arcs(k,1)),3)
  A=(4d0*r1*r1+t*t+s*s+r*r)/2d0
  B=t*r
  C=s*r
  S=dsqrt(A*A-B*B-C*C) 
  rr=r*r-A
  if (dabs(dabs(arcs(k,3))-2d0*pi).LT.eps_two_pi) then
!                 full circle arc
    vIone=2d0*pi/S
    vItwo=2d0*pi*A/(S**3)
    vIthree=pi*(2d0*A*A+B*B+C*C)/(S**5)
    vJone=pi+rr/2d0*vIone
    vJtwo=(vIone+rr*vItwo)/4d0
    vJthree=(vItwo+rr*vIthree)/8d0
    delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
    delta_aint=2d0*vJone*r1**2
    if (arcs(k,3).LT.0) then
      delta_vint=-delta_vint
      delta_aint=-delta_aint
    endif
    avi(1)=avi(1)+delta_vint    
    avi(2)=avi(2)+delta_aint    
  else
!         integration over arcs
    if (arcs(k,3).LT.0) then 
      al=arcs(k,2)+arcs(k,3)
      be=arcs(k,2) 
    else
      be=arcs(k,2)+arcs(k,3)
      al=arcs(k,2) 
    endif 
    vIone=2d0*(pi/2d0-datan((A*dcos((be-al)/2d0)+B*dcos((al+be)/2d0)+C*dsin((al+be)/2d0))/(S*dsin((be-al)/2d0))))/S
    sb=dsin(be)
    cb=dcos(be)
    sa=dsin(al)
    ca=dcos(al)
    vItwo=(fract(A,B,C,sb,cb,1)-fract(A,B,C,sa,ca,1)+A*vIone)/(S*S)
    vIthree=(fract(A,B,C,sb,cb,2)-fract(A,B,C,sa,ca,2)+(fract(A,B,C,sb,cb,1)-fract(A,B,C,sa,ca,1))/A+(2d0*A*A+B*B+C*C)*vItwo/A)/(2d0*S*S)
    vJone=((be-al)+rr*vIone)/2d0
    vJtwo=(vIone+rr*vItwo)/4d0
    vJthree=(vItwo+rr*vIthree)/8d0
    delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
    delta_aint=2d0*vJone*r1**2
    if (arcs(k,3).LT.0) then
      delta_vint=-delta_vint
      delta_aint=-delta_aint
    endif
    avi(1)=avi(1)+delta_vint    
    avi(2)=avi(2)+delta_aint    
  endif
enddo
return
end subroutine

real*8 function fract(A,B,C,sinphi,cosphi,k)

!     Fraction evaluation for integral

implicit none
integer k
real*8 a,b,c,sinphi,cosphi

fract=(-B*sinphi+C*cosphi)/(A+B*cosphi+C*sinphi)**k

return
end function

subroutine readpdbna(input,na)
implicit none
integer*4 j,kode,na
character input*(*)
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
end subroutine

subroutine readpdb(input)
use comun
implicit none
integer*4 j,kode,kk
character input*256
character ln*256,ann*5,rnn*4,frmt*64
real*8 e
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
      read (ln,'(6x,A5,x,A4,x,A4,x,A4,4x,3F8.3,2F6.2,6x,A4)') ann,typ(j),res(j),rnn,rc(1,j),rc(2,j),rc(3,j),e,w(j),segid(j)
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

