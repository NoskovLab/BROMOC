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

subroutine shell_simul
use apfmod
use efpmod
use ioxmod
use stdiomod
use grandmod
use gsbpmod
use constamod
use errormod
use nucleotmod
use charfuncmod   !command parser
use extramod
use splinemod

!Command parser and file name
implicit none
! allocatables
real,allocatable       :: xy(:,:),cg2(:) 
real,allocatable ::    epsLJ(:),sgLJ(:)
real*4,allocatable :: efield(:,:)
integer, allocatable :: nxf(:)  
logical*1, allocatable :: Qefpotread(:)
character*1,allocatable :: secstr(:)
character com*2048,word*1024,fnam*80,wrd5*5
integer unvec(maxopen), numb, totnumb, nions, kode, maxpart, ncl3, in1, in2
character*80 title
character*4 wrd4,ionname
!real*16 dener
real dener
real battery
real r1,r2,r6,z1,v1,v2,y1,y2,x1,x2,x3,xm,ym,zm,z2
integer ix1,iy1,iz1,ix2,iy2,iz2 
real*4 idv
real resol
real vc1(3),vc2(3),vc3(3)
logical*1 endlog, logfinal, Qlsprmf, doions, dodna,Qadj
logical*1 logmemb, logmmij, logphix, logphiv, logsrpmf,logbuff,Qefpott,Qepwrt,logrfpar,Qnohead
logical*1 Qexpl2nd,Qinputpar,Qonlychden
real*8 zero
!for time
integer*8       start,finish,timer
character*8  date
character*10 time
character*5  zone
integer   values(8)
integer   is,cnt,mnp,nnp,nn,maxd
integer   wunit,s1,s2,s3,iseed,wallsi4
integer*1 walls,dnaparams

!for trajectory, fraction of denatued bases files and translocation time
integer    ntras 
real       vfbs,scald,scaldd

!for security outputfile
integer    nsec,ntc

!for short-range ion-ion interactions
integer nindex

!for TEST order
real       tol, delta, maxl,minl,maxlg,minlg

!for MEMBRANE and ION-ION orders
!for effective dielectric constant for DNA
real       conc
!for solvent-induced contribution
real       anumb1, epsn
!for DNA fixed sites
integer      isite
!for DNA matrix rotation
real       sum1, sum2, sum3
real       xold, yold, zold
!forgotten to declare
integer r1i,z1i,itype,jtype,i,j,k,iat,ib,ifirst,ii,ij,ik,ilast
integer iunit,lfnam,ngcmc,nmcm,nbd,ntypold
integer*8 ncycle,nsave,nsfbs
integer,allocatable :: iunitv(:)
real cc0,cc1,cc2,cc3,cc4
real xtras, ytras, ztras, rot(3,3)

!Default parameters and options
!------------------------------
do i = 1, maxopen
  unvec(i) = -1
enddo
Qpar         = .false.
Qsystem      = .false.
Qbuf         = .false.
Qdeby        = .false.
Qdebyhyb     = .false.
Qsolv        = .false.
Qnucl        = .false.
Qdie         = .false.
Qtraj        = .false.
Qtrajcont    = .false.
Qepwrt       = .false.
Qapfor       = .false.
Qwarn        = .false.
Qcontrans    = .false.
Qcontprint   = .false.
Qcountion    = .false.
Qunsplit     = .false.
Qchden       = .false.
frmt         = ''
Qefpott      = .false.
Qenergy      = .true.
Qforces      = .false.
Qnobond      = .true.
Qecyl        = .false.
Qnonbond     = .true.
Qgr          = .false.
Qmemb        = .false.
Qmmij        = .false.
Qphix        = .false.
Qphiv        = .false.
Qsrpmf       = .false.
Qlsprmf      = .false.
Qnmcden      = .false.
Qdiffuse     = .false.
Qprofile     = .false.
Qproxdiff    = .false.
Qrfpar       = .false.
Qljpar       = .false.
Qljsin       = .false.
Qninfo       = .false.
Qpres        = .false.
Qnotrans     = .false.
Qsvdw        = .false.
Qrfpsin      = .false.
doions       = .false.
dodna        = .false.
iseed        = 3141593
ntype        = 0
nold         = 0
ndna         = 0
nbuffer      = 0
nfix         = 0
nsites       = 0
istrs        = 0
inuc         = 0
cgnuc        = -1.0
diffnuc      = 0.1
epsnuc       = 0.184
ionstr       = 0.0
temp         = 300.0
eps          = 0.0
sigma        = 0.0
lx           = 0.0
ly           = 0.0
lz           = 0.0
cx           = 0.0
cy           = 0.0
cz           = 0.0
cdie         = 80.0
Rsphe        = 0.0
afact        = 0.0
kappa        = 0.0
kbtdna       = 0.0
ikbtdna      = 0.0
setframes    = 0
dnaparams    = 0
scalepairing = 1.0
do is = 1, dspline
  xs(is)=0.0
  ys(is)=0.0
  b(is)=0.0
  c(is)=0.0
  d(is)=0.0
enddo

call header(outu)
start = timer()
call date_and_time(date, time, zone, values)
write(outu,'(6x,a,7(i0,a))') 'Started at (YYYY-MM-DD HH:mm:ss.ms): ',values(1),'-',values(2),'-',values(3),' ',values(5),':',values(6),':',values(7),'.',values(8)
write(outu,*)
write(outu,*)

logfinal = .false. 
do while (.not. logfinal)
  call getlin(com,inpu,outu) ! new commands line
  write(outu,'(/a)')'**********************************************'
  write(outu,'(a)') trim(adjustl(com))
  write(outu,'(a/)')'**********************************************'
  call getfirst(com,wrd5) ! new key word
  wrd5=lcase(wrd5)
  !.....MISCELLANEOUS COMMANDS    
  if  (wrd5.eq.'title') then
  !     ---------------
     call getwrd(com,'name',title)
     write(outu,'(6x,a)') 'TITLE: '//trim(title)
  ! **********************************************************************
  elseif (wrd5.eq.'open')then
  !        ---------------
     ! unit [integer,default=1]
     call gtipar(com,'unit',iunit,1)
     if (iunit.le.0) call error ('shell_simul', 'unit is a zero or negative number', faterr)
     if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
     if (unvec(iunit).ne.-1) call error ('shell_simul', 'unit incorrect in OPEN order',faterr)
     call lualloc(unvec(iunit))
     iunit = unvec(iunit)
     ! namefile [character]
     call getwrd(com,'name',fnam)
     lfnam=len_trim(fnam)
     write(outu,105) fnam(1:lfnam),iunit
  105     format(6x,'open file ',a,' as unit ',i3)
     write(outu,*)
     ! write/ read options
     if (check(com,'write')) then
       if (check(com,'file ')) then ! binary file
         open(unit=iunit,file=fnam,form='unformatted')
       else
         open(unit=iunit,file=fnam,form='formatted')
       endif
     elseif (check(com,'read')) then
       if (check(com,'file ')) then
         open(unit=iunit,file=fnam,status='old',form='unformatted')
         rewind(unit=iunit) ! rewinds a file to the beginning
       else
         open(unit=iunit,file=fnam,status='old',form='formatted')
         rewind(unit=iunit) ! rewinds a file to the beginning
       endif
     endif
  ! **********************************************************************
  elseif (wrd5.eq.'close') then
  !        ---------------
     ! unit [integer,default=-1]
     call gtipar(com,'unit',iunit,1)
     if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
     if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
     if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in CLOSE order',faterr)
     ij = unvec(iunit) 
     unvec(iunit) = -1
     iunit = ij
     close(unit=iunit)
     call lunalloc (iunit)
     write(outu,'(6x,a,1x,i3)') 'close unit',iunit
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'nucle') then
  !         ---------------
     if (.not.Qsystem) call error('shell_simul','SYSTEM must be defined previously',faterr)
     if (Qpar) call error ('shell_simul', 'NUCLEOTIDES order is defined after PARTICLE order', faterr)
     ! number of strand [default=0]
     call gtipar (com,'strands',istrs,0)
     if (istrs.gt.2 .or. istrs.lt.0) call error ('shell_simul', 'Incorrect numbers of strands in NUCLEOTIDES order', faterr)
     ! number of nucleotides of each strand [default=0]
     call gtipar (com,'nucleot',inuc,0)
     if (inuc.lt.0) call error ('shell_simul', 'Incorrect numbers of nucleotides in NUCLEOTIDES order', faterr)
! allocate variables
     maxsite=istrs*inuc*3
     mmxsites = maxsite*(maxsite+1)/2
     maxbond = (maxsite-1)*3
     maxang  = (maxsite-2)*4
     maxdihe = (maxsite-3)*4
     maxstack = maxsite*(maxsite-1)/2
     maxbp    = maxsite*(maxsite-1)/2
     maxex    = maxsite*(maxsite-1)/2
     maxqq    = maxsite*(maxsite-1)/2
     maxsolv  = maxsite*(maxsite-1)/2
     allocate (strand(maxsite),typenuc(maxsite),stfx(maxsite),stfree(maxsite),namnucl(maxsite),namsite(maxsite))
     allocate (xnat(maxsite),ynat(maxsite),znat(maxsite),rnat(maxsite),phinat(maxsite))
     allocate (sitebond(maxbond,2),siteangle(maxang,3),sitedihe(maxdihe,4),distbond(maxbond),valangle(maxang))
     allocate (valdihe(maxdihe),bond(mmxsites),angle(mmxsites),typbond(maxbond))
     allocate (sitestack(maxstack,2),sitebp(maxbp,2),siteex(maxex,2),siteqq(maxqq,2),siteslv(maxsolv,2))
     allocate (sgstack(maxstack),sgbp(maxbp),sgex(maxex))
     stfree     = .true.
     ! charge nucleotides [default=-1]  
     call gtdpar(com,'charge',cgnuc,-1.0)
     ! diffusion constant for nucleotides [default=0.1]
     call gtdpar(com,'diffusion',diffnuc,0.01)
     if (diffnuc.lt.0.0) call error ('shell_simul', 'Diffusion coefficient for each nucleotide is negative', faterr)
     ! Parameter to calculate bonded and non-bonded potential
     ! energy terms [default=0.26]        
     call gtdpar(com,'dnatemp',dnatemp,temp)
     if (dnatemp.ne.temp) then
       write(outu,'(6x,a)') 'Different Temperature for DNA activated' 
       write(outu,'(6x,a,f10.5)') '   DNA Temperature: ',dnatemp
       write(outu,'(6x,a,f10.5)') '   Defined DNA Diffusivity: ',diffnuc
       write(outu,'(6x,a)') '   Rescaling DNA Diffusivity proportionally to DNA Temperature' 
       diffnuc=dnatemp/temp*diffnuc
       write(outu,'(6x,a,f10.5)') '   New DNA Diffusivity: ',diffnuc
       kbtdna = kboltz*dnatemp/kcalmol
       ikbtdna = 1.0/kbtdna
     endif
     call gtdpar(com,'eps',epsnuc,epsnuc)
     call gtdpar(com,'scalepairing',scalepairing,scalepairing)
     if (epsnuc.lt.0.0) call error ('shell_simul', 'epsnuc for each nucleotide is negative', faterr)
     Qsolv = check(com,'solv') .and. istrs.eq.2
     if (Qsolv) then
       call gtdpar (com,'conc',conc,0.0)
       if (conc.lt.0.0) call error ('shell_simul', 'Molarity is negative in NUCLEOTIDES order', faterr)
       if (conc.eq.0.0) call error ('shell_simul', 'Conc is zero (needed for solv)', warning)
       epsn = 0.504982*epsnuc*(1.0-1.0/(1.40418-0.268231*inuc))
       anumb1 = 0.474876*(1.0 + 1.0/(0.148378+10.9553*conc))
       epsolv = epsn*anumb1
       write(outu,'(6x,a)') 'Solvent-induced contribution enabled (Recommended to improve double stranded DNA interactions)'
       write(outu,'(6x,a,f10.5)') 'Solvent-induced Parameter= ',epsolv
     endif
     Qinvstr=check(com,'3-5')
     Qexpl2nd=check(com,'explicit2nd')
     QfirstP=check(com,'keepfirstp')
     Qinputpar=check(com,'inputparam')
     if (check(com,'depablo')) then
       dnaparams=0
     elseif (check(com,'charmmvac')) then 
       dnaparams=1
     elseif (check(com,'charmmwat')) then 
       dnaparams=2
     elseif (check(com,'charmmwatbc')) then
       dnaparams=3
     else ! charmmwat
       dnaparams=3
     endif
     if (Qexpl2nd.and.istrs.ne.2) call error('shell_simul','explicit2nd ignored, 2nd strand was not specified',warning)
     Qninfo=check(com,'printinfo')
     if (istrs.gt.0) then 
       if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of date in NUCLEOTIDE order', faterr)
       do j = 1, inuc
         ntype = ntype + 1
         if (ntype.gt.dtype) call error ('shell_simul', 'ntype is greater than dtype', faterr)
         call getfirst(com,word)
         if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of date in NUCLEOTIDE order', faterr)
         ! nucleotide name [character*4]
         atnam(ntype) = ucase(word)
         ! charge nucleotide
         cg(ntype) = cgnuc
         ! diffusion constant
         diffusion(ntype) = diffnuc
         ! Phosphate
         if ((QfirstP.and.j.eq.1).or.j.gt.1) then
           nsites = nsites + 1
           if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
           strand(nsites) = 1
           typenuc(nsites) = ntype
           namnucl(nsites) = atnam(ntype)(1:1)
           namsite(nsites) = 'P'
         endif
         ! Base
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 1                 
         typenuc(nsites) = ntype
         namnucl(nsites) = atnam(ntype)(1:1)
         if (namnucl(nsites).eq.'A') then
           namsite(nsites) = 'Ab' 
         else if (namnucl(nsites).eq.'T') then
           namsite(nsites) = 'Tb'
         else if (namnucl(nsites).eq.'C') then
           namsite(nsites) = 'Cb'
         else if (namnucl(nsites).eq.'G') then
           namsite(nsites) = 'Gb'
         else      
           call error ('shell_simul', 'nucleotide name is not correct', faterr)
         endif  
         ! Sugar
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 1
         typenuc(nsites) = ntype
         namnucl(nsites) = atnam(ntype)(1:1) 
         namsite(nsites) = 'S'
       enddo
     endif 
     ! read explicit second strand
     if (Qexpl2nd) then
       if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of second strand in NUCLEOTIDE order', faterr)
       allocate (secstr(inuc))
       do j=inuc,1,-1
         call getfirst(com,word)
         if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of second strand in NUCLEOTIDE order', faterr)
         ! nucleotide name [character*4]
         secstr(j) = ucase(word)
       enddo
     endif
     if (Qinputpar) then
       if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
       do i=1,6
         do j=1,3
           call getfirst(com,word)
           if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
           !Bases
           cylall(i,j) = chr2real(word)
         enddo
       enddo
       call getfirst(com,word)
       if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
       !Internuclear Distance
       din = chr2real(word)
       call getfirst(com,word)
       if (len_trim(word).eq.0) call error ('shell_simul', 'premature end of inputparam in NUCLEOTIDE order', faterr)
       !Internuclear Angle
       ain = chr2real(word)
       write(outu,'(6x,a)') 'Using User-defined structural parameters for DNA'
     else
       !Cylindrical polar coordinates
       if (dnaparams.eq.0) then 
         !Bases
         cylall(1,1:3) = (/ 0.773, -41.905, -0.051 /)         !Ab
         cylall(2,1:3) = (/ 2.349, -86.119, -0.191 /)         !Tb
         cylall(3,1:3) = (/ 2.296, -85.027, -0.187 /)         !Cb
         cylall(4,1:3) = (/ 0.828, -40.691, -0.053 /)         !Gb
         cylall(5,1:3) = (/ 6.981, -70.197, -1.280 /)         !Sugar
         cylall(6,1:3) = (/ 8.918, -94.038, -2.186 /)         !Phosphate
         din = 3.38
         ain = 36.0
         write(outu,'(6x,a)') 'Using internal De Pablo et al structural parameters for DNA'
       elseif (dnaparams.eq.1) then
         cylall(1,1:3) = (/ 1.5401186532694178, -19.177377752137673, 4.81713596295395771E-002 /)
         cylall(2,1:3) = (/ 2.5874118642077302, -64.556009869361233, 0.22903343233126944      /)
         cylall(3,1:3) = (/ 2.2651275581525985, -69.277179138939957, 0.29825389290205251      /)
         cylall(4,1:3) = (/ 1.6463744611104871, -24.136045962260361, 6.06394873617855691E-003 /)
         cylall(5,1:3) = (/ 7.2988301635144150, -63.884837335190518, 0.45948475145162820      /)
         cylall(6,1:3) = (/ 9.4038525018460799, -89.846061205360229, -0.48592487490127378     /)
         din = 3.3692353762855434 
         ain = 36.790780189800941
         write(outu,'(6x,a)') 'Using structural parameters of CHARMM DNA minimized in vacuo'
       elseif (dnaparams.eq.2) then
         cylall(1,1:3) = (/ 0.83216456560990049,-43.802827117829452, 5.84700711558472225E-002 /)
         cylall(2,1:3) = (/ 2.25434746257693770,-87.586682816786521, 7.58267641147968852E-002 /) 
         cylall(3,1:3) = (/ 2.49227490170237240,-87.699281096559588, 5.24429792951401769E-002 /)
         cylall(4,1:3) = (/ 0.89858765931206686,-21.581894968064670, 0.11155621478344679      /)
         cylall(5,1:3) = (/ 7.01490325321993600,-68.281680932114654, 1.34340627519265107E-002 /)
         cylall(6,1:3) = (/ 9.08267487151476870,-92.072561423506826,-1.0389212165439070       /)
         din = 3.4190096123453726 
         ain = 36.291088852071780
         write(outu,'(6x,a)') 'Using structural parameters of CHARMM DNA minimized in implicit water'
       elseif (dnaparams.eq.3) then
         cylall(1,1:3) = (/ 2.6019796769397963,-81.191725312579493, 0.15848098213497963      /)
         cylall(2,1:3) = (/ 3.3706977830909244,-96.328909516456363, 0.24317884044268287      /)
         cylall(3,1:3) = (/ 3.4919037468548892,-90.384267464738386, 0.24690507752582472      /)
         cylall(4,1:3) = (/ 2.2051192772253185,-68.561436289845915, 0.10838747455031919      /)
         cylall(5,1:3) = (/ 7.0149106362072891,-68.281465812483361, 1.34359749915350605E-002 /)
         cylall(6,1:3) = (/ 9.0826733119583505,-92.072394528330037,-1.0389241599469166       /)
         din = 3.4189906453957799
         ain = 36.290944969867944
         write(outu,'(6x,a)') 'Using structural parameters of CHARMM DNA minimized in implicit water with bases centered at geometric center of bases heavy atoms'
       endif
     endif
     write(outu,'(/6x,a)')       'Particle   xy-dist    xy-ang         z' 
     write(outu,'(6x,a,3f10.4)') 'Ab      ',(cylall(1,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'Tb      ',(cylall(2,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'Cb      ',(cylall(3,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'Gb      ',(cylall(4,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'S       ',(cylall(5,i),i=1,3)
     write(outu,'(6x,a,3f10.4)') 'P       ',(cylall(6,i),i=1,3)
     write(outu,'(/6x,2(a,f10.4)/)') 'Internucleotide: Distance= ',din,' Angle= ',ain

     if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in NUCLEOTIDE order', faterr)
     if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in NUCLEOTIDE order', faterr)
     if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in NUCLEOTIDE order', faterr)
     ! complementary strand
     nsites1st=nsites
     if (istrs.eq.2) then
       ntypold = ntype      
       do j = 0, inuc-1
         ntype = ntype + 1
         if (ntype.gt.dtype) call error ('shell_simul', 'ntype is greater than dtype', faterr)
         if (Qexpl2nd) then
           atnam(ntype)=secstr(j+1)
         else
           ! nucleotide name [character*4]
           if (atnam(ntypold-j)(1:1).eq.'A') then
             atnam(ntype) = 'T'
           else if (atnam(ntypold-j)(1:1).eq.'T') then
             atnam(ntype) = 'A'
           else if (atnam(ntypold-j)(1:1).eq.'C') then
             atnam(ntype) = 'G'
           else
             atnam(ntype) = 'C'
           endif     
         endif 
         ! charge nucleotide
         cg(ntype) = cgnuc
         ! diffusion constant
         diffusion(ntype) = diffnuc          
         ! Phosphate
         if ((QfirstP.and.j.eq.0).or.j.gt.0) then
           nsites = nsites + 1
           if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
           strand(nsites) = 2
           typenuc(nsites) = ntype
           namnucl(nsites) = atnam(ntype)(1:1)
           namsite(nsites) = 'P'
         endif
         ! Base
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 2
         typenuc(nsites) = ntype
         namnucl(nsites) = atnam(ntype)(1:1)           
         if (namnucl(nsites).eq.'A') then
           namsite(nsites) = 'Ab'
         else if (namnucl(nsites).eq.'T') then
           namsite(nsites) = 'Tb'
         else if (namnucl(nsites).eq.'C') then
           namsite(nsites) = 'Cb'
         else if (namnucl(nsites).eq.'G') then
           namsite(nsites) = 'Gb'
         else
           call error ('shell_simul', 'nucleotide name is not correct', faterr)
         endif
         ! Sugar
         nsites = nsites + 1
         if (nsites.gt.maxsite) call error ('shell_simul', 'nsites is greater than maxsite', faterr)
         strand(nsites) = 2
         typenuc(nsites) = ntype
         namnucl(nsites) = atnam(ntype)(1:1)
         namsite(nsites) = 'S'            
       enddo
       if (Qexpl2nd) deallocate (secstr)
     endif
     nold = ntype           
     fctn=celec*cgnuc**2/cdie
  
     write(outu,*)
     write(outu,'(6x,a,i3,a)') 'There are ',istrs*inuc,' nucleotides'
     if (Qinvstr) then
       write(outu,'(6x,a)') 'DNA strands are written in 3-5 direction'
     else
       write(outu,'(6x,a)') 'DNA strands are written in 5-3 direction'
     endif
     if (QfirstP) then 
       write(outu,'(6x,a)') "First 5' Phosphate included"
     else
       write(outu,'(6x,a)') "First 5' Phosphate not included"
     endif
     write(outu,'(6x,a)') 'EPSILON  -> energy scale [Kcal/mol]'
     write(outu,'(6x,a)') 'CHARGE   -> phosphate sites charge [e]'
     write(outu,'(6x,a)') 'DIFFUSION-> diffusion constant [Ang.**2/ps]'
     write(outu,'(6x,a)') 'STRAND  NAME  EPSILON  CHARGE  DIFFUSION' 
     write(outu,'(6x,a)') '----------------------------------------'
     do i = 1, ntype
       if (i.le.inuc) then
         write(outu,'(6x,i3,4x,a1,2x,f8.3,2x,f8.3,2x,e9.2)') 1, atnam(i), epsnuc, cgnuc, diffnuc
       else 
         write(outu,'(6x,i3,4x,a1,2x,f8.3,2x,f8.3,2x,e9.2)') 2, atnam(i), epsnuc, cgnuc, diffnuc        
       endif  
     enddo
     write(outu,*)
     ! Interaction site coordinates for the native structure
     call native_structure
     ! bond streching terms
     call bonds
     ! bond angles terms
     call angles
     ! torsinal angles terms
     call dihedral
     ! nonbonded terms
     call go_qq
     Qnucl = .true.
     call reasig
  
     write(outu,*)  
  ! **********************************************************************
  elseif (wrd5.eq.'parti') then
  !       ---------------
    if (.not.Qsystem) call error('shell_simul','SYSTEM must be defined previously',faterr)
    if (Qdeby) call error('shell_simul','IMPLICIT IONS defined in system, PARTICLE section is senseless',faterr)
    endlog = .false.
    do while (.not.endlog)
      call getlin(com,inpu,outu) ! new commands line 
      endlog = check(com,'end')
      if (.not.endlog) then
        ntype = ntype + 1
        if (ntype.gt.dtype) call error ('shell_simul', 'ntype is greater than dtype', faterr)
        ! name ion type [character*4]
        call getfirst(com,atnam(ntype))
        ! particle charge [real*8,default=0]
        call gtdpar(com,'charge',cg(ntype),0.0)
        ! diffusion constant [real*8,default=0.1]
        call gtdpar(com,'diffusion',diffusion(ntype),0.1)
        if (diffusion(ntype).lt.0.0) call error ('shell_simul', 'diffusion coefficient is negative', faterr)
      endif
    enddo
    Qpar = .true.
    call reasig
    ntc=(nttyp+ndna+1)*nion/2
    allocate (Qchr(ntc),Qefpot(ntc),Qlj(ntc),Qsrpmfi(ntc),cg2(nttyp))
    Qchr=.false.
    Qefpot=.false.
    Qlj=.false.
    Qsrpmfi=.false.
    cg2=0.0
  
  ! Assign charges
    if (ndna.gt.1) cg2(2)=cgnuc
    do i=nold+1,ntype
      cg2(nwtype(i))=cg(i)
    enddo
  
    ! ASSIGN CHARGES USING CDIE
    allocate (fct(ntc))
    fct=0.0
    do i=1,nttyp
      do j=i,nttyp
        if (j.gt.ndna) then
          if (cg2(i).ne.0.0.and.cg2(j).ne.0.0) then
            is=nindex(i,j)
            fct(is)=celec*cg2(i)*cg2(j)/cdie
            Qchr(is)=.true.
          endif
        endif
      enddo
    enddo
  
    write(outu,*)
    write(outu,'(6x,a,i3,a)') 'There are ',ntype-nold,' atom types'
    write(outu,'(6x,a)') 'CHARGE -> ion charge [e]'
    write(outu,'(6x,a)') 'DIFFUSION -> diffusion constant [Ang.**2/ps]'
    write(outu,'(6x,a)') 'NAME---CHARGE---DIFFUSION'
    do i = nold+1, ntype
      write(outu,'(6x,a,2f8.3)') atnam(i),cg(i),diffusion(i)
    enddo
    write(outu,*)
  
  ! **********************************************************************
  elseif (wrd5.eq.'apfor') then
  !       ----------------
    if (.not. Qnucl) call error ('shell_simul', 'APFOR order is defined before NUCLEOTIDE order', faterr)
    ! number of DNA fixed sites [integer, default=0]
    call gtipar(com,'afn',afn,0)  
    if (afn.gt.nsites .or. afn.lt.0) call error ('shell_simul', 'Incorrect number of DNA sites',faterr)
    if (allocated(af)) deallocate (af)
    if (allocated(sn)) deallocate (sn)
    allocate (af(3,afn),sn(afn))
    if (afn.gt.0) then
      do j = 1, afn
        if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in APFOR. New line expected', faterr)
        if (.not.setint(isite,com)) call error ('shell_simul', 'premature end of data in APFOR. Integer Expected.', faterr)
        if (isite.gt.nsites .or. isite.le.0) call error ('shell_simul', 'Incorrect DNA site', faterr)
        call gtdpar(com,'fx',af(1,j),0.0)
        call gtdpar(com,'fy',af(2,j),0.0)
        call gtdpar(com,'fz',af(3,j),0.0)
        sn(j)=isite
      enddo
      Qapfor=.true.
    endif  
    if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in APFOR order', faterr)
    if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in APFOR order', faterr)
    if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in APFOR order', faterr)
    if (Qapfor) then
      write(outu,*)
      write(outu,'(6x,a,i3,a)') 'There are ',afn,' DNA sites to be forced:'
      do i=1,afn
        write(outu,'(6x,2i4,x,a,x,3(3x,f10.5))') i,sn(i),namsite(sn(i)),af(1:3,i)
      enddo
    endif           
  ! **********************************************************************
  elseif (wrd5.eq.'fixsi') then
  !       ----------------
    if (.not. Qnucl) call error ('shell_simul', 'FIXSITE order is defined before NUCLEOTIDE order', faterr)
    ! number of DNA fixed sites [integer, default=0]
    call gtipar(com,'nstfx',nstfx,0) 
    if (nstfx.gt.nsites .or. nstfx.lt.0) call error ('shell_simul', 'Incorrect number of DNA fixed sites', faterr)
    Qstfx = nstfx.gt.0
    if (Qstfx) then
      if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in FIXSITE order', faterr)
      do j = 1, nstfx
        if (.not.setint(isite,com)) call error ('shell_simul', 'premature end of data in FIXSITE order', faterr)
        if (isite.gt.nsites .or. isite.le.0) call error ('shell_simul', 'Incorrect DNA fixed site',faterr)
        stfx(j) = isite
        stfree(isite) = .false.
      enddo
    endif 
    if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in FIXSITE order', faterr)
    if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in FIXSITE order', faterr)
    if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in FIXSITE order', faterr)
    if (Qstfx) then
      write(outu,*)
      write(outu,'(6x,a,i3,a)') 'There are ',nstfx,' DNA fixed sites:'
      write(outu,*) '     ',(stfx(j),'(',namsite(stfx(j)),')  ',j=1,nstfx)
    endif
  ! **********************************************************************
  elseif (wrd5.eq.'notra') then
  !        ---------------
  ! Do not translate, DNA keep the geometric center in a fixed position if off
    if (.not.Qnucl) call error ('shell_simul', 'NOTRANSL must be defined after NUCLEOTIDE', faterr)
    Qnotrans=.not.check(com,'off')
    Qnotrx=check(com,'x')
    Qnotry=check(com,'y')
    Qnotrz=check(com,'z')
    if (.not.(Qnotrx.or.Qnotry.or.Qnotrz)) Qnotrans=.false.
    insites=1.0/nsites
    notrx=sum(x(1:nsites))*insites
    notry=sum(y(1:nsites))*insites
    notrz=sum(z(1:nsites))*insites
    if (Qnotrans) then 
      write(outu,'(6x,a)') 'NOTRANSLATION is on'
      if (Qnotrx) write(outu,'(6x,a,f9.3)') '  DNA centroid will be fixed at x =',notrx
      if (Qnotry) write(outu,'(6x,a,f9.3)') '  DNA centroid will be fixed at y =',notry
      if (Qnotrz) write(outu,'(6x,a,f9.3)') '  DNA centroid will be fixed at z =',notrz
    else
      write(outu,'(6x,a)') 'NOTRANSLATION is off'
    endif
  ! **********************************************************************
  elseif (wrd5.eq.'contr') then
  !        ---------------
  ! Constrain translation with an harmonic potential of DNA by the geometric center
    if (.not.Qnucl) call error ('shell_simul', 'CONTRANSL must be defined after NUCLEOTIDE', faterr)
    Qcontprint=check(com,'print')
    call gtipar(com,'ctn',ctn,0)
    if (allocated(kx)) deallocate (kx)
    if (allocated(ky)) deallocate (ky)
    if (allocated(kz)) deallocate (kz)
    if (allocated(csn)) deallocate (csn)
    if (allocated(contrx)) deallocate (contrx)
    if (allocated(contry)) deallocate (contry)
    if (allocated(contrz)) deallocate (contrz)
    allocate (kx(ctn),ky(ctn),kz(ctn),csn(ctn),contrx(ctn+1),contry(ctn+1),contrz(ctn+1))
    Qcontrans=.not.check(com,'off')
    Qunsplit=check(com,'unsplit').and.istrs.eq.2
    if (.not.Qcontrans) ctn=0
    do i = 1, ctn
      if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in CONTRA. New line expected', faterr)
      if (.not.setint(isite,com)) call error ('shell_simul', 'premature end of data in CONTRA. Integer Expected.', faterr)
      if (isite.gt.nsites .or. isite.lt.0) call error ('shell_simul', 'Incorrect DNA site', faterr)
      call gtdpar(com,'kx',kx(i),0.0)
      call gtdpar(com,'ky',ky(i),0.0)
      call gtdpar(com,'kz',kz(i),0.0)
      csn(i)=isite
      if (csn(i).eq.0) then
        insites=1.0/nsites
        if (Qunsplit) then
          contrx(i)=2.0*sum(x(1:nsites1st))*insites
          contry(i)=2.0*sum(y(1:nsites1st))*insites
          contrz(i)=2.0*sum(z(1:nsites1st))*insites
          contrx(ctn+1)=2.0*sum(x(nsites1st+1:nsites))*insites
          contry(ctn+1)=2.0*sum(y(nsites1st+1:nsites))*insites
          contrz(ctn+1)=2.0*sum(z(nsites1st+1:nsites))*insites
        else
          contrx(i)=sum(x(1:nsites))*insites
          contry(i)=sum(y(1:nsites))*insites
          contrz(i)=sum(z(1:nsites))*insites
        endif
      else
        contrx(i)=x(csn(i))
        contry(i)=y(csn(i))
        contrz(i)=z(csn(i))
      endif
      call gtdpar(com,'x',contrx(i),contrx(i))
      call gtdpar(com,'y',contry(i),contry(i))
      call gtdpar(com,'z',contrz(i),contrz(i))
      Qcontrans=Qcontrans.or.kx(i).ne.0.0.or.ky(i).ne.0.0.or.kz(i).ne.0.0
    enddo
    if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of data in CONTRA order', faterr)
    if (.not.setword(word,com)) call error ('shell_simul', 'premature end of data in CONTRA order', faterr)
    if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in CONTRA order', faterr)
    if (Qcontrans) then
      write(outu,'(/6x,a)') 'CONTRANSLATION is on'
      write(outu,'(6x,a,i3,a)') 'There are ',ctn,' DNA sites/centroid to be constrained:'
      if (Qunsplit) write(outu,'(6x,a)') 'Double stranded DNA will not be splitted'
      do i=1,ctn
        if (csn(i).eq.0) then
          write(outu,'(6x,2i4,2(x,a,x,3(3x,f10.5)))') i,csn(i),'centroid ',kx(i),ky(i),kz(i),'at',contrx(i),contry(i),contrz(i)
        else
          write(outu,'(6x,2i4,2(x,a,x,3(3x,f10.5)))') i,csn(i),namsite(csn(i))//'       ',kx(i),ky(i),kz(i),'at',contrx(i),contry(i),contrz(i)
        endif
      enddo
    else
      write(outu,'(6x,a)') 'CONTRANSLATION is off'
    endif
  
  ! **********************************************************************
  elseif (wrd5.eq.'syste') then
  !        ---------------
    Qwarn = check(com,'showwarn')
    if (Qwarn) write(outu,'(6x,a)') 'Detailed warnings enabled'
    Qdnafree = check(com,'dnafree')
    if (Qdnafree) write(outu,'(6x,a)') 'DNA not forced to be within system boundaries'
    ! random number seed [integer,default=3141593]
    call gtipar(com,'iseed',iseed,iseed)
    call feedseed(iseed) ! if zero or lower will use cpu clock to generate seed
    write(outu,'(6x,a,i0)') 'Using seed: ',iseed
    ! Temperature
    call gtdpar(com,'temp',temp,temp)
    if (temp.lt.0.0) call error ('shell_simul', 'Temperature must be positive',faterr)
    write(outu,'(6x,a,f8.3)') 'Temperature: ',temp
    ! Effective dielectric constant estimated from temperatura and salt concentration
    Qdie = check(com,'calcdie')
    if (Qdie) then
      call gtdpar(com,'conc',conc,0.0)
      if (conc.lt.0.0) then
        call error ('shell_simul', 'Molarity is negative in SIMUL order', faterr)
      elseif (conc.ge.0.0) then
        cdnuc = (1.0-0.2551*conc+5.151e-2*conc**2-6.889e-3*conc**3)*(249.4-0.788*temp+7.20e-4*temp**2)
        cdie=cdnuc
        write(outu,'(6x,a,f8.2,a,f8.5,a,f8.4)')'Dielectric Constant estimated from   Temp=',temp,'   Salt Conc=',conc,'  is ',cdie
      endif
    else 
      ! Dielectric solvent [real*8,default=80 (water dielectric
      ! constant)]
      call gtdpar(com,'cdie',cdie,cdie)
      if (cdie.lt.0.0) then
        call error ('shell_simul', 'Dielectric constant is negative in SRPMF order', faterr)
      endif
      write(outu,'(6x,a,f8.3)') 'Dielectric constant: ',cdie
    endif
  
    ! Logical variable which indicates if Coulomb interactions
    ! are taken into account using the Debye-Hückel approximation
    Qdeby = check(com,'debye')
    Qdebyhyb = check(com,'debyhyb')
    ! Ionic strength [default=0]
    if (Qdeby.and.Qdebyhyb) call error ('shell_simul','Hybrid and Normal Debye-Huckel cannot be combined.',faterr)
    call gtdpar(com,'ionic',ionstr,ionstr)
    if (ionstr.lt.0.0) then
      call error ('shell_simul', 'Ionic strength is negative', faterr)
    elseif (ionstr.gt.0.0) then
       write(outu,'(6x,a,1x,f10.5)') 'Ionic Strength: ',ionstr
    endif
    if ((Qdeby.or.Qdebyhyb).and.ionstr.eq.0.0) then
      if (Qdeby) Qdeby=.false.
      if (Qdebyhyb) Qdebyhyb=.false.
      call error ('shell_simul','Debye-Huckel approximation disabled. Ionic Strength cannot be zero.',warning)
    endif
    if (Qdeby.and.ionstr.gt.0.0) then
       write(outu,'(6x,a)') 'Debye-Huckel approximation enabled'
    endif
    if (Qdebyhyb.and..not.Qdnafree) then
      Qdebyhyb=.false.
      call error ('shell_simul','Hybrid Debye-Huckel approximation disabled. Use it in combination with DNAFREE.',warning)
    endif
    if (.not.Qdebyhyb.and..not.Qdeby.and.Qdnafree) then
      call error ('shell_simul','DNAFREE is recommended to be used with Normal or Hybrid Debye-Huckel (DEBYHYB) approximation.',warning)
    endif
    if (Qdebyhyb.and.Qdnafree.and.ionstr.gt.0.0) then
       write(outu,'(6x,a)') 'Hybrid Debye-Huckel approximation enabled'
    endif

    ! Transmembrane potential [real*8,default=0]
    call gtdpar(com,'volt',voltage,0.0)
    write(outu,'(6x,a,1x,f10.5)') 'Transmembrane Potential (Volts): ',voltage
     
    voltage = voltage*Coulomb/kcalmol
  
    ! Activates calculation of Virial Pressure
    Qpres = check(com,'pres')
    if (Qpres) write(outu,'(6x,a)') 'Virial Pressure enabled'
  
    ! define constant kappa
    if (ionstr.gt.0.0) then 
      kappa = sqrt(const*cdie*temp/ionstr)
      ikappa = 1.0/kappa ! screening factor
    endif
  
    ! define more constants
    kBT      = kboltz*temp/kcalmol
    ikbt = 1.0/kbt
    kbtdna = kbt
    ikbtdna = ikbt
    
    Qecyl=check(com,'ecyl') 
    Qsphere = check(com,'sphere') ! logical*1 variable
    if (Qecyl.and.Qsphere) call error ('shell_simul','Elliptical Cylinder and Spherical system cannot be used together',faterr)
    if (Qsphere) then
       ! radius of spherical system [real*8,default=0]     
       call gtdpar(com,'radi',rsphe,rsphe)
       if (Rsphe.lt.0.0)  call error ('shell_simul', 'radi is lower than zero in SYSTEM order', faterr)
       tvol=4.0/3.0*pi*rsphe**3
       rsphe2=rsphe**2
       lx = 2.0*rsphe
       ly = 2.0*rsphe
       lz = 2.0*rsphe
       maxl=2.0*rsphe+1.0
    elseif (Qecyl) then
      ! Elliptical Cylinder
       ! minor or mayor axis of ellipse
       call gtdpar(com,'lx',lx,lx)
       if (lx.lt.0.0) call error ('shell_simul', 'LX is lower than zero in SYSTEM order', faterr)
       ! minor or mayor axis of ellipse
       call gtdpar(com,'ly',ly,ly)
       if (ly.lt.0.0) call error ('shell_simul', 'LY is lower than zero in SYSTEM order', faterr)
       ! length of cylinder
       call gtdpar(com,'lz',lz,lz)
       if (lz.lt.0.0) call error ('shell_simul', 'LZ is lower than zero in SYSTEM order', faterr)
       tvol=0.25*pi*lx*ly*lz
       iecx=2.0/lx
       iecy=2.0/ly
       maxl=sqrt(max(lx,ly)**2+lz**2)+1.0
    else 
       ! Orthorombic box size along the X-axis [real*8,default=0]
       call gtdpar(com,'lx',lx,lx)
       if (lx.lt.0.0) call error ('shell_simul', 'LX is lower than zero in SYSTEM order', faterr)
       ! Orthorombic box size along the Y-axis [real*8,default=0]
       call gtdpar(com,'ly',ly,ly)
       if (ly.lt.0.0) call error ('shell_simul', 'LY is lower than zero in SYSTEM order', faterr)
       ! Orthorombic box size along the Z-axis [real*8,default=0]
       call gtdpar(com,'lz',lz,lz)
       if (lz.lt.0.0) call error ('shell_simul', 'LZ is lower than zero in SYSTEM order', faterr)
       tvol=lx*ly*lz
       maxl=sqrt(lx**2+ly**2+lz**2)+1.0
    endif
    call gtdpar(com,'cx',cx,cx)        ! Center of System along X-axis
    call gtdpar(com,'cy',cy,cy)        ! Center of System along Y-axis
    call gtdpar(com,'cz',cz,cz)        ! Center of System along Z-axis
    lx2p = cx+0.5*lx
    ly2p = cy+0.5*ly
    lz2p = cz+0.5*lz
    lx2m = cx-0.5*lx
    ly2m = cy-0.5*ly
    lz2m = cz-0.5*lz
    if (Qnucl) fctn=celec*cgnuc**2/cdie
    if (Qpar) then
      fct=0.0
      do i=1,nttyp
        do j=i,nttyp
          if (j.gt.ndna) then
            if (cg2(i).ne.0.0.and.cg2(j).ne.0.0) then
              is=nindex(i,j)
              fct(is)=celec*cg2(i)*cg2(j)/cdie
              Qchr(is)=.true.
            endif
          endif
        enddo
      enddo
    endif
    maxpart=int(avogadro*1e-27*3.0*tvol)+1  ! no more than 3 Molar of particles in the system volume
    if (maxpart.gt.datom) maxpart=datom 
    Qsystem = .true. 
    
    if (Qsphere) then
      write(outu,'(6x,a,f12.3)') 'Spherical system, Radius ',Rsphe
    elseif (Qecyl) then
      write(outu,'(6x,a,2f12.3,a,f12.3)') 'Elliptical Cylinder system. Ellipse diameters x and y:',lx,ly,' Cylinder Length: ',lz
    else
      write(outu,'(6x,3(a,f12.3))') 'Box:  LX ',lx,'  LY ',ly,'  LZ ',lz
    endif
    write(outu,'(6x,3(a,f12.3))') 'Center:  X ',cx,'  Y ',cy,'  Z ',cz
    write(outu,'(6x,a,f20.3)') 'System Total Volume (Ang**3): ',tvol
  ! **********************************************************************
  elseif (wrd5.eq.'buffe') then
  !        ---------------
     if (.not.Qpar) call error ('shell_simul', 'BUFFER order is defined before PARTICLE order', faterr)
     if (.not.Qsystem) call error ('shell_simul', 'BUFFER order is defined before SYSTEM order', faterr)
     endlog = .false.
     do while (.not.endlog)
       call getlin(com,inpu,outu) ! new commands line
       endlog = check(com,'end')
       if (.not.endlog) then
         nbuffer = nbuffer + 1
         if (nbuffer.gt.dbuff) call error ('shell_simul', 'nbuffer greater than dbuff', faterr)
         ! Obtention of ion type atnam(itype)
         call getfirst(com,wrd4)
         call fatnam(atnam,ntype,wrd4,itype)
         ibfftyp(nbuffer)=itype
         ! Intrinsic chemical potential [real*8,default=0] 
         call gtdpar(com,'mu',mu(nbuffer),0.0)
         ! To avoid unexpected large deviation from average number 
         ! of ions in the buffer regions, Kb [default=0] 
         Qbufferbias(nbuffer)= check(com,'bufferbias')
         call gtdpar(com,'kb',kb(nbuffer),0.0)
         call gtdpar(com,'lzmin',LZmin(nbuffer),lz2m) 
         call gtdpar(com,'lzmax',LZmax(nbuffer),lz2p) 
         if (LZmin(nbuffer).gt.0.0) then 
           if (LZmin(nbuffer).ge.LZmax(nbuffer)) call error ('shell_simul', 'LZmin => Lzmax in BUFFER order', faterr) 
         else 
           if (LZmax(nbuffer).lt.0.0.and.abs(LZmin(nbuffer)).le.abs(LZmax(nbuffer))) call error('shell_simul', 'LZmin => Lzmax in BUFFER order',faterr) 
         endif 
         if (Qsphere) then
           ! Minimum position of buffer sphere [real*8,default=0]
           call gtdpar(com,'rmin',Rmin(nbuffer),0.0)
           if (Rmin(nbuffer).lt.0.0) call error ('shell_simul', 'Rmin is lower than zero in BUFFER order', faterr)
           if (Rmin(nbuffer).gt.Rsphe) call error ('shell_simul', 'spherical region is too large in BUFFER order', warning)
           ! Maximum position of buffer sphere [real*8,default=r1]
           call gtdpar(com,'rmax',Rmax(nbuffer),Rsphe)
           if (Rmax(nbuffer).lt.0.0) call error ('shell_simul', 'Rmax is lower than zero in BUFFER order', faterr)
           if (Rmax(nbuffer).gt.Rsphe) call error ('shell_simul', 'spherical region is too large in BUFFER order', faterr)
           if (Rmin(nbuffer).ge.Rmax(nbuffer)) call error ('shell_simul', 'Rmin is greater or equal to Rmax in in BUFFER order', faterr)
           if (LZmin(nbuffer).gt.0.0) then 
             if (LZmin(nbuffer).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmin has an incorrect value', faterr) 
             if (LZmax(nbuffer).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmax has an incorrect value', faterr) 
           else 
             if (abs(LZmin(nbuffer)).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmin has an incorrect value', faterr) 
             if (abs(LZmax(nbuffer)).gt.Rmax(nbuffer)) call error ('shell_simul', 'LZmax has an incorrect value', faterr) 
           endif 
           ! Obtention of buffer volume (spherical system)
           r1 = Rmin(nbuffer)
           r2 = Rmax(nbuffer)
           logbuff = LZmin(nbuffer).eq.lz2m.and.LZmax(nbuffer).eq.lz2p
           if (logbuff) then 
             v1 = 2.0*twopi*r1**3/3.0 
             v2 = 2.0*twopi*r2**3/3.0 
           else 
             if (LZmin(nbuffer).gt.0.0) z1 = LZmin(nbuffer) 
             if (LZmax(nbuffer).lt.0.0) z1 = abs(LZmax(nbuffer))
             if (r1.gt.z1) then 
               v1 = twopi*(r1**3/3.0-(z1*r1**2*0.5-z1**3/6.0)) 
               v2 = twopi*(r2**3/3.0-(z1*r2**2*0.5-z1**3/6.0)) 
             else 
               v1 = 0.0
               v2 = twopi*(r2**3/3.0-(z1*r2**2*0.5-z1**3/6.0)) 
             endif 
           endif 
           volume(nbuffer) = v2-v1
         elseif (Qecyl) then
           ! Obtention of buffer volume (elliptical cylinder system)      
           volume(nbuffer) = 0.25*pi*lx*ly*(LZmax(nbuffer)-LZmin(nbuffer))
         else
           ! Obtention of buffer volume (orthorombic system)      
           volume(nbuffer) = (LX*LY*(LZmax(nbuffer)-LZmin(nbuffer)))
         endif
         ! Concentration for a specific ion in a buffer
         ! [real*8,default=0]
         call gtdpar(com,'conc',density(nbuffer),0.0)
         if (density(nbuffer).lt.0.0) call error ('shell_simul', 'conc is lower than zero in BUFFER order', faterr)
         ! Obtention of average number of ions in the buffer regions
         density(nbuffer) = density(nbuffer)*(avogadro/liter)*(angstrom**3)
         ! avogadro=6.022045D23; liter=0.001; angstrom=1.0E-10 
         avnum(nbuffer) = density(nbuffer)*volume(nbuffer)
         if (int(avnum(nbuffer)).le.0) call error ('shell_simul','Buffer volume or ion density is too low',warning)
         ! Average number of ions in the buffer regions [real*8,default=0]
         call gtdpar(com,'aver',avnum(nbuffer),avnum(nbuffer))
         if (avnum(nbuffer).lt.0.0) call error ('shell_simul', 'aver is lower than zero in BUFFER order', faterr)
         ! Obtention of concentration for a specific ion in a buffer
         density(nbuffer) = avnum(nbuffer)/volume(nbuffer)
         ! Transmembrane potential for a buffer
         call gtdpar(com,'volt',battery,0.0)
         ! OBTENTION OF THE ELECTROCHEMICAL POTENTIAL  
         mu(nbuffer) = mu(nbuffer)+cg(itype)*battery*Coulomb/kcalmol
       endif
     enddo
     Qbuf = .true.
  
     write(outu,*)
     write(outu,'(6x,a,i3,a)') 'There are ',nbuffer,' buffers'
     write(outu,'(6x,a)') 'MU    -> chemical potential [Kcal/mol]'
     if (Qsphere) then
       if (.not.logbuff) then
         write(outu,'(6x,a)') 'LZmin -> minimum position along Z-axis [Ang.]'
         write(outu,'(6x,a)') 'LZmax -> maximum position along Z-axis [Ang.]'
       endif
       write(outu,'(6x,a)') 'Rmin -> minimum radius of buffer sphere [Ang.]'
       write(outu,'(6x,a)') 'Rmax -> maximum radius of buffer sphere [Ang.]'
     else
       write(outu,'(6x,a)') 'LZmin -> minimum position along Z-axis [Ang.]'
       write(outu,'(6x,a)') 'LZmax -> maximum position along Z-axis [Ang.]'
     endif  
     write(outu,'(6x,a)') 'AVER  -> average number for an ion'
     write(outu,'(6x,a)') 'DENS  -> density for an ion [Ang.**(-3)]'
     write(outu,'(6x,a)') 'VOL   -> volume for a buffer [Ang.**3]'
     write(outu,'(6x,a)') 'KB    -> to avoid unexpected large deviation from an average number of ions'
     if (.not.Qsphere) then
       write(outu,'(6x,a)')   'NAME------MU----------LZmin-------LZmax-------AVER------DENS--------VOL----------KB----'
     else
       if (logbuff) then
         write(outu,'(6x,a)') 'NAME------MU----------Rmin--------Rmax----AVER----DENS----VOL----KB----'
       else
         write(outu,'(6x,a)') 'NAME------MU----------LZmin-------LZmax---Rmin----Rmax----AVER---DENS----VOL----KB----'
       endif
     endif
     do ib = 1, nbuffer
        ! Initializations
        nremove(ib)= 0
        ninsert(ib)= 0
        if (Qsphere) then
          if (logbuff) then
            write(outu,'(6x,a,4f12.5,1x,2e12.4,f10.5)') atnam(ibfftyp(ib)),mu(ib),Rmin(ib),Rmax(ib),avnum(ib),density(ib),volume(ib),kb(ib)
          else
            write(outu,'(6x,a,6f12.5,1x,2e12.4,f10.5)') atnam(ibfftyp(ib)),mu(ib),LZmin(ib),LZmax(ib),Rmin(ib),Rmax(ib),avnum(ib),density(ib),volume(ib),kb(ib)
          endif
        else
          write(outu,'(6x,a,4f12.5,1x,2e12.4,f10.5)') atnam(ibfftyp(ib)),mu(ib),LZmin(ib),LZmax(ib),avnum(ib),density(ib),volume(ib),kb(ib)
        endif
     enddo
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'ljsin') then
    if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'LJSIN order is defined before PARTICLE and/or NUCLEOTIDE order', faterr)
    if (Qljpar) call error ('shell_simul', 'LJSIN is defined after LJPAR', faterr)
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type atnam(itype)
        call getfirst(com,wrd4)
        call fatnam(atnam2,nttyp,wrd4,itype)
        call gtdpar(com,'epsilon',eps(itype),0.0)
        call gtdpar(com,'sigma',sigma(itype),0.0)
      endif
    enddo
    Qljsin = .true.
  
    write(outu,*)
    write(outu,'(6x,a)') 'LJ Single Parameters:'
    write(outu,'(6x,a)') '--------------'
    write(outu,'(6x,a)') 'type---epsilon(kcal/mol)---sigma(Ang)'
    write(outu,'(6x,a)') '-------------------------------------'
    do  i = 1, nttyp
      if (eps(i).ne.0.0.and.sigma(i).ne.0.0) then
        write(outu,'(6x,a,2f12.4)') atnam2(i),eps(i),sigma(i)
      endif
    enddo
    ntc=(nttyp+ndna+1)*nion/2
    allocate (epp4(ntc),sgp2(ntc),epsLJ(ntc),sgLJ(ntc))
    epsLJ=0.0
    sgLJ=0.0
    epp4=0.0
    sgp2=0.0
    do i=1,nttyp
      do j=i,nttyp
        if (j.gt.ndna) then
          if (eps(i).gt.0.0.and.sigma(i).gt.0.0.and.eps(j).gt.0.0.and.sigma(j).gt.0.0) then
              is=nindex(i,j)
              epp4(is)=4.0*sqrt(eps(i)*eps(j))
              sgp2(is)=(0.5*(sigma(i)+sigma(j)))**2
              Qlj(is)=.true.
          else
            write(outu,'(6x,a,x,a,x,a)')'Warning: Missing single LJ parameters to compute pairs for:',atnam2(i),atnam2(j)
          endif
        endif
      enddo
    enddo
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'ljpar') then
  !        ---------------     
    if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'LJPAR order is defined before PARTICLE and/or NUCLEOTIDE order', faterr)
  !         if (Qionsite) then
  !           call error ('shell_simul', 'Combination rules for LJ'
  !     & //' parameters are desactivated', warning)
  !           Qionsite = .false.
  !         endif
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type atnam(itype)
        call getfirst(com,wrd4)
        call fatnam(atnam2,nttyp,wrd4,itype)
        ! Obtention of ion type atnam(jtype)
        call getfirst(com,wrd4)
        call fatnam(atnam2,nttyp,wrd4,jtype)
        is=nindex(itype,jtype)
        call gtdpar(com,'epsilon',epsLJ(is),0.0)
        call gtdpar(com,'sigma',sgLJ(is),0.0)
      endif
    enddo
    Qljpar = .true.
  
    write(outu,*)
    write(outu,'(6x,a)') 'LJ Pair Parameters:'
    write(outu,'(6x,a)') '-------------'
    write(outu,'(6x,a)') 'type1---type2---epsilon(kcal/mol)---sigma(Ang)'
    write(outu,'(6x,a)') '-----------------------------------------'
    if (.not.allocated(epp4)) then
      ntc=(nttyp+ndna+1)*nion/2
      allocate (epp4(ntc),sgp2(ntc))
      epp4=0.0
      sgp2=0.0
    endif
    do  i = 1, nttyp
      do j = i, nttyp
        if (j.gt.ndna) then 
          is=nindex(i,j)
          if (epsLJ(is).gt.0.0 .and. sgLJ(is).gt.0.0) then
            write(outu,'(6x,a,a,2f12.4,$)') atnam2(i),atnam2(j),epsLJ(is),sgLJ(is) 
            if (Qlj(is)) then
              write(outu,'(a)') ' (replaced)'
            else
              write(outu,*)
              Qlj(is)=.true.
            endif
            epp4(is)=4.0*epsLJ(is)
            sgp2(is)=sgLJ(is)**2
          else
            if (Qlj(is)) then
              write(outu,'(6x,a,a,2f12.4,a)') atnam2(i),atnam2(j), epp4(is)*0.25,sqrt(sgp2(is)),' (from single LJ)'
            else
              write(outu,'(6x,a,x,a,x,a)') 'Warning: Missing LJ parameters to compute pairs for:',atnam2(i),atnam2(j)
            endif
          endif
        endif
      enddo
    enddo
    deallocate (epsLJ,sgLJ)
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'simul') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'SIMULATION order is defined before PARTICLE and/or NUCLEOTIDES orders', faterr)
     if (Qpres.and..not.Qpar) then
       write(outu,'(6x,a)') 'Warning: Pressure cannot be calculated if free particles are absent, so deactivated'
       Qpres=.false.
     endif
     if (Qpar.and..not.(Qefpott.or.Qljsin.or.Qljpar)) call error ('shell_simul', 'No effective potential and no Lennard Jones parameters defined', warning)
     if (.not.allocated(warn)) allocate (warn(nttyp))
     warn=0
     ! number of steps for BD or MC simulations
     ! [integer,default=0]
     call gti8par(com,'ncycle',ncycle,0)
     if (ncycle.le.0) call error ('shell_simul', 'ncycle is equal or lower than zero in SIMUL order', faterr)
     ! number of steps for GCMC [integer,default=0]
     call gtipar(com,'ngcmc',ngcmc,0)
     if (ngcmc.lt.0) call error ('shell_simul','ngcmc is lower than zero in SIMUL order', faterr)
     if (ngcmc.gt.0 .and. .not.Qbuf) then
       call error ('shell_simul', 'GCMC is turned off because buffers have not been defined', warning)
       ngcmc = 0
     endif 
     ! number of steps for MC [integer,default=0]
     call gtipar(com,'nmcm',nmcm,0)
     if (nmcm.lt.0) call error ('shell_simul', 'nmcm is lower than zero in SIMUL order', faterr)
     if (nmcm.gt.0 .and. .not.Qpar) then
       call error ('shell_simul', 'Metropolis MonteCarlo is turned off because ions have not been defined', warning)
       nmcm = 0
     endif 
     if (nmcm.gt.0) then
     ! maximum displacement (mcmax) for MC
     ! [real*8,default=1]
       call gtdpar(com,'mcmax',mcmax,0.5)
       if (mcmax.lt.0.0) call error ('shell_simul', 'mcmax is lower than zero in SIMUL order', faterr)
       mcmax=2.0*mcmax
     endif
     ! number of steps for BD [integer,default=0]  
     call gtipar(com,'nbd',nbd,0)
     if (nbd.lt.0) call error ('shell_simul', 'nbd is lower than zero in SIMUL order', faterr)
     ! maximum displacement (bdmax) for BD
     if (nbd.gt.0) then
       call gtdpar(com,'bdmax',bdmax,1.0)
       if (bdmax.lt.0.0) call error ('shell_simul', 'bdmax is lower than zero in SIMUL order', faterr)
     endif
     ! Logical variable which indicates if a trajectory file
     Qchdencnt = check(com,'chden')
     if (.not.Qchden) Qchdencnt=.false.
     ! Logical variable which indicates if a trajectory file
     ! will be written
     Qtraj = check(com,'traject')
     if (Qtraj) then
       ! unit number of a trajectory file [integer,default=1]
       Qtrajcont=check(com,'trajcont')
       call gtipar(com,'setframes',setframes,0)
       call gtipar(com,'iuntrj',iuntrj,0)
       if (iuntrj.le.0) call error ('shell_simul', 'iuntrj is zero or a negative number', faterr)
       if (iuntrj.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuntrj).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuntrj = unvec(iuntrj)
     endif
     ! trajectory saving frequency [integer,default=ncycle]
     call gti8par(com,'nsave',nsave,ncycle)
     if (nsave.le.0) call error ('shell_simul', 'nsave is lower than or equal to zero in SIMUL order', faterr)
     if (mod(ncycle,nsave).ne.0) call error ('simul1', 'nsave is not correct. ncycle must be divisible by nsave', faterr)
  
     ! time-step for BD [real*8,default=0.02]
     call gtdpar(com,'dt',dt,0.02)
     if (dt.le.0.0) call error ('shell_simul', 'dt is lower or equal than zero', faterr)
     
     ! print frequency [integer,default=0]
     call gtipar(com,'nprint',nprint,0)
     if (nprint.lt.0) call error ('shell_simul', 'nprint is a negative number',faterr)

     ! Z-position where the number of ion crossing the channel 
     Qcountion=check(com,'countions').and.Qpar
     if (Qcountion) then
       call gtipar(com,'svcntfq',svcntfq,nprint)
       call gtipar(com,'iuncnt',iuncnt,0)
       if (iuncnt.le.0) call error ('shell_simul', 'iuncnt is zero or a negative number', faterr)
       if (iuncnt.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuncnt).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuncnt = unvec(iuncnt)
       call gtcrpar(com,'zcont',cntpts,word)
       if (cntpts.le.0) call error ('shell_simul','zcont cannot be empty or ommited if countions present',faterr)
       allocate (zcont(cntpts),nbackward(1:ntype-nold,cntpts),nforward(1:ntype-nold,cntpts))
       call gtcdpar(word,zcont)
     endif

     ! Logical variable which indicates if radial distribution 
     ! function will be printed
     Qgr      = check(com,'rdf')
     if (Qgr .and. .not.Qpar) then
       call error ('shell_simul', 'radial distribution function is turned off because ions have not been defined',warning)
       Qgr = .false.
     endif
     ! Logical variable which indicates if average density profile 
     ! along the Z-axis will be printed
     Qrho     = check(com,'rho')
     if (Qrho .and. (.not.Qpar.or..not.Qsystem)) then
       call error ('shell_simul', 'average density profile along the Z-axis is turned off because ions and/or system have not been defined', warning)
       Qrho = .false.
     endif
     ! Logical variable which indicates if average density profile 
     ! along the Z-axis will be printed for DNA sites
     Qrdna    = check(com,'rhDNA')
     if (Qrdna .and. (.not.Qnucl.or..not.Qsystem)) then
        call error ('shell_simul', 'average density profile along the Z-axis is turned off because sites and/or system have not been defined', warning)
       Qrdna = .false.
     endif
     ! Logical variable which indicates if various average probability 
     ! distributions will be printed
     Qprob    = check(com,'prob')
     if (Qprob .and. .not.Qpar) then
        call error ('shell_simul', 'average probability distributions are turned off because ions have not been defined', warning)
       Qprob = .false.
     endif 
     ! ion pairing analysis frequency [integer,default=1]
     call gtipar(com,'nanal',nanal,1)
     if (nanal.le.0) call error ('shell_simul', 'nanal is lower or equal than zero in SIMUL order', faterr)
     ! Logical variable which indicates if ion pairing analysis 
     ! (S frequency) will be made
     Qionpair = check(com,'ionpair')
     if (Qionpair .and. (.not.Qpar.or..not.Qsystem)) then
       call error ('shell_simul','ion pairing analysis is turned off because ions and/or system have not been defined', warning)
       Qionpair = .false.
     endif
     ! Logical variable which indicates if energy profile along 
     ! the Z-axis will be made
     Qenerprofile = check(com,'enerprofile')
     if (Qenerprofile .and. (.not.Qpar.or..not.Qsystem)) then
       call error ('shell_simul', 'energy profile along the Z-axis is turned off because ions and/or system have not been defined', warning)
       Qenerprofile = .false.
     endif
     ! Logical variable which indicates if fraction of denatured 
     ! bases will be calculated after simulation of DNA
     Qfbases = check(com,'denatured') 
     ! Logical variable which indicates if translocation duration
     ! for DNA through pore will be calculated
     Qfmemb = check(com,'translocation') .and. Qnucl 
     ntras = 1
     if (Qfmemb) then
       ! unit number of a translocation duration file [integer,default=1]
       call gtipar(com,'unitn',iuntfm,1)
       if (iuntfm.le.0) call error ('shell_simul', 'iuntfm is zero or a negative number', faterr)
       if (iuntfm.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuntfm).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuntfm = unvec(iuntfm)
       call gtipar(com,'ntras',ntras,1)
       ! translocation saving frequency [real*8,deafult=1]
       if (ntras.le.0) call error ('shell_simul', 'ntras < = 0 in SIMUL order',faterr)
       if (mod(ncycle,ntras).ne.0) call error ('simul1', 'ntras is not correct. ncycle must be divisible by ntras', faterr)
     endif
     if (Qprob .or. Qfmemb) then
       ! maximum position of a channel along the 
       ! Z-axis [real*8,default=0]
       call gtdpar(com,'czmax',czmax,0.0)
       ! minimum position of a channel along the 
       ! Z-axis [real*8,default=0]     
       call gtdpar(com,'czmin',czmin,0.0) 
       if (czmin.gt.czmax) call error ('shell_simul', 'minimum position is greater than maximum position of the channel', faterr)
     endif
     if (Qgr) then
       if (nfix.eq.0) then
         call error ('shell_simul', 'RDF cannot be calculated because nfix = 0', warning)
         Qgr = .false.
       endif 
     endif
     if (Qgr) then
       call gtipar(com,'ion',igr,nsites+1)
       if (igr.le.0 .or. igr.gt.(nsites+nfix)) call error ('shell_simul', 'fixed reference ion is not adequate in SIMUL order', faterr)
     endif
     nsfbs = ncycle
     vfbs = 0.0
     if (Qfbases) then 
       if (.not.Qnucl) call error ('shell_simul', 'Fraction of denatured bases cannot calculated because DNA has not been defined', faterr) 
       if (istrs.ne.2) call error ('shell_simul', 'Fraction of denatured bases cannot calculated because DNA has only one strand', faterr) 
       if (inuc.eq.0) call error ('shell_simul', 'Fraction of denatured bases cannot calculated because the number of nucleotides is zero',faterr) 
       ! unit number of a fraction of denatured bases [integer,default=1]
       call gtipar(com,'iunfbs',iunfbs,1) 
       if (iunfbs.le.0) call error ('shell_simul', 'iunfbs is zero or negative number', faterr)
       if (iunfbs.gt.maxopen) call error ('shell_simul', 'iunfbs is greater than maxopen', faterr)
       if (unvec(iunfbs).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order', faterr)
       iunfbs = unvec(iunfbs)
       ! fraction of denatured bases saving frequency
       ! [integer, default=ncycle]
       call gti8par(com,'nsfbs',nsfbs,ncycle)
       if (nsfbs.le.0) call error ('shell_simul', 'nsfbs < = 0 in SIMUL order',faterr)
       if (mod(ncycle,nsfbs).ne.0) call error ('simul1', 'nsfbs is not correct. ncycle must be divisible by nsfbs', faterr)
  
       ! initial time for calculating  fraction of denatured
       ! [real*8,deafult=1]
       call gtdpar(com,'vfbs',vfbs,0.0)
       if (vfbs.lt.0.0) call error ('shell_simul', 'vfbs < 0 n SIMUL order',faterr)
     endif 
     ! Security ouputfile contains coordinates and seed numbers 
     Qsec = check(com,'security') 
     nsec = 1
     if (Qsec) then
       ! unit number of a security file [integer,default=1]
       call gtipar(com,'units',iuntsc,1)
       if (iuntsc.le.0) call error ('shell_simul', 'iuntfm is zero or a negative number', faterr)
       if (iuntsc.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen',faterr)
       if (unvec(iuntsc).eq.-1) call error ('shell_simul', 'unit incorrect in SIMUL order',faterr)
       iuntsc = unvec(iuntsc)
       call gtipar(com,'nsec',nsec,1)
       ! saving frequency security [real*8,deafult=1]
       if (nsec.le.0) call error ('shell_simul', 'nsec < = 0 in SIMUL order',faterr)
       if (mod(ncycle,nsec).ne.0) call error ('simul1', 'nsec is not correct. ncycle must be divisible by nsec', faterr)
     endif
     ! turn off total energy
     Qenergy  = .not.check(com,'noenergy')
     ! turn off nonbond energy
     Qnonbond = .not.check(com,'nononbond')
     ! turn off bond energy
     Qnobond = .not.check(com,'nobond')
  
     call simul1(ncycle, ngcmc, nmcm, nbd, nsave, nsfbs, vfbs, ntras, nsec, iseed)
     
     if (sum(warn).gt.0) then
       write(outu,*) '      Warnings summary:'
       do i=1,nttyp
         write(outu,*) '               ',atnam2(i),warn(i)
       enddo
     endif
  ! **********************************************************************
  elseif (wrd5.eq.'energ') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'ENERGY order is defined before PARTICLE and/or NUCLEOTIDES orders', faterr)
     ! default values 
     ! Qnobond = Qnonbond =.true.
     ! Qmemb = Qmmij = Qphix = Qphiv = Qsrpmf = .false.
     Qnobond  = .not.check(com,'nobond')
     Qnonbond = .not.check(com,'nononbond')
     if (check(com,'membrane')) then
       if (.not.Qmemb) call error ('shell_simul', 'Planar membrane has not be defined before ENERGY order', faterr)
     else
       logmemb = Qmemb
       Qmemb = .false.
     endif
     if (check(com,'mmij')) then
       if (.not.Qmmij) call error ('shell_simul', 'Reaction Field has not be defined before ENERGY order', faterr)
     else
       logmmij = Qmmij
       Qmmij = .false.
     endif
     if (check(com,'phix')) then
       if (.not.Qphix) call error ('shell_simul', 'Static Field has not be defined before ENERGY order', faterr)
     else
      logphix = Qphix
      Qphix = .false.
     endif
     if (check(com,'phiv')) then
       if (.not.Qphiv) call error ('shell_simul', 'Repulsive term has not be defined before ENERGY order', faterr)
     else
       logphiv = Qphiv
       Qphiv = .false.
     endif
     if (check(com,'srpmf')) then
       if (.not.Qsrpmf) call error ('shell_simul', 'Short-range interaction term has not be defined before ENERGY order', faterr)
     else
       logsrpmf = Qsrpmf
       Qsrpmf = .false.
     endif
     if (check(com,'rfpar')) then
       if (.not.Qrfpar) call error ('shell_simul', 'Reaction Field parameter term has not be defined before ENERGY order', faterr)
     else
       logrfpar = Qrfpar
       Qrfpar = .false.
     endif
  
     write(outu,*)
     write(outu,'(6x,a)') 'ENERGY calculation'
     write(outu,'(6x,a)') '------------------'
     if (Qnobond)  write(outu,'(6x,a)') 'Bonding Energy Term'
     if (Qnonbond) write(outu,'(6x,a)') 'Nonbonding Energy Term'
     if (Qmemb)    write(outu,'(6x,a)') 'Planar membrane Term'
     if (Qmmij)    write(outu,'(6x,a)') 'Reaction Field Energy Term'
     if (Qphix)    write(outu,'(6x,a)') 'External Field Energy Term'
     if (Qphiv)    write(outu,'(6x,a)') 'Repulsive Energy Term'
     if (Qsrpmf)   write(outu,'(6x,a)') 'Short-range Interaction Term'
     if (Qrfpar)   write(outu,'(6x,a)') 'Reaction Field Parameter Term'
  
     call ENERGY
     write(outu,'(6x,a,f12.6)') 'Total energy ',ener
     Qmemb = logmemb 
     Qmmij = logmmij 
     Qphix = logphix 
     Qphiv = logphiv 
     Qsrpmf = logsrpmf 
     Qrfpar = logrfpar
  ! **********************************************************************
  elseif (wrd5.eq.'inter') then 
  !        ---------------
    if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'INTERACT order is defined before PARTICLE and/or NUCLEOTIDES orders', faterr)
    ! atom type [integer,default=1]
    call gtipar(com,'atom',iat,1) 
    ! default values
    ! Qnobond = Qnonbond =.true.
    ! Qmemb = Qmmij = Qphix = Qphiv = Qsrpmf = .false.         
    Qnobond  = .not.check(com,'nobond')
    Qnonbond = .not.check(com,'nononbond')
    if (check(com,'membrane')) then
      if (.not.Qmemb) call error ('shell_simul', 'Planar membrane has not be defined before INTERACT order', faterr)
    else
      logmemb = Qmemb
      Qmemb = .false.
    endif
    if (check(com,'mmij')) then
      if (.not.Qmmij) call error ('shell_simul', 'Reaction Field has not be defined before INTERACT order', faterr)
    else
      logmmij = Qmmij
      Qmmij = .false.
    endif
    if (check(com,'phix')) then
      if (.not.Qphix) call error ('shell_simul', 'Static Field has not be defined before INTERACT order', faterr)
    else
     logphix = Qphix
     Qphix = .false.
    endif
    if (check(com,'phiv')) then
      if (.not.Qphiv) call error ('shell_simul', 'Repulsive term has not be defined before INTERACT order', faterr)
    else
      logphiv = Qphiv
      Qphiv = .false.
    endif
    if (check(com,'srpmf')) then
      if (.not.Qsrpmf) call error ('shell_simul', 'Short-range interaction term has not be defined before INTERACT order', faterr)
    else
      logsrpmf = Qsrpmf
      Qsrpmf = .false.
    endif
    if (check(com,'rfpar')) then
      if (.not.Qrfpar) call error ('shell_simul', 'Reaction Field Parameter term has not be defined before INTERACT order', faterr)
    else
      logrfpar = Qrfpar
      Qrfpar = .false.
    endif
  
    write(outu,*)
    write(outu,'(6x,a)') 'INTERACT calculation'
    write(outu,'(6x,a)') '--------------------'         
    if (Qnobond)  write(outu,'(6x,a)') 'Bonding Energy Term' 
    if (Qnonbond) write(outu,'(6x,a)') 'Nonbonding Energy Term'
    if (Qmemb)    write(outu,'(6x,a)') 'Planar membrane Term'
    if (Qmmij)    write(outu,'(6x,a)') 'Reaction Field Energy Term'
    if (Qphix)    write(outu,'(6x,a)') 'External Field Energy Term'
    if (Qphiv)    write(outu,'(6x,a)') 'Repulsive Energy Term'
    if (Qsrpmf)   write(outu,'(6x,a)') 'Short-range Interaction Term'
    if (Qrfpar)   write(outu,'(6x,a)') 'Reaction Field Parameter Term'
  
    ! Calculate the interaction of particle "iat" with the rest of
    ! the system
    if (iat.le.nsites) then
      itype = typenuc(iat)
    else 
      itype = abs(typei(iat))
    endif
    call interact(dener,x(iat),y(iat),z(iat),itype,iat,.true.)
    write(outu,'(6x,a,f12.6)') 'Interaction of particle ',iat,dener
    Qmemb = logmemb 
    Qmmij = logmmij 
    Qphix = logphix 
    Qphiv = logphiv 
    Qsrpmf = logsrpmf
    Qrfpar = logrfpar
  ! **********************************************************************
  elseif (wrd5.eq.'membr') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'MEMBRANE order is defined before PARTICLE and/or NUCLEOTIDE order', faterr)
     if (.not.Qsystem) call error ('shell_simul', 'MEMBRANE order is defined before SYSTEM order', faterr)
     ! Turn on/off repulsive membrane potential in cylindrical pore
     ! region
     Qpore = check(com,'pore')
     ! Thickness of a membrane [real*8,default=0]
     call gtdpar(com,'thick',thick2,0.0)
     if (thick2.le.0.0) call error ('shell_simul', 'Thickness of a membrane is negative or zero', faterr)
     ! Position of center of a membrane along the Z-axis
     ! [real*8,default=0]
     call gtdpar(com,'zmemb',zmemb,0.0)
     if (zmemb.gt.lz2p.or.zmemb.lt.lz2m) call error ('shell_simul', 'Position of center of a membrane along the Z-axis is not correct', faterr)
     tmemb  =  thick2 ! thickness of the membrane
     thick2 =  thick2*0.5
     zmemb1 = -thick2+zmemb ! lower limit of the membrane
     zmemb2 =  thick2+zmemb ! upper limit of the membrane
     if (zmemb2.gt.lz2p) call error ('shell_simul', 'Upper limit of the membrane along the Z-axis is not correct', faterr)
     if (zmemb1.lt.lz2m) call error ('shell_simul', 'Lower limit of the membrane along the Z-axis is not correct', faterr)
     ! Dielectric constant of the membrane [real*8,default=2]
     call gtdpar(com,'epsm',epsm,2.0)
     if (epsm.le.0.0) call error ('shell_simul', 'Dielectric constant of the membrane has a null or negative value', faterr)
     if (ionstr.le.0.0) call error ('shell_simul', 'Ion Strength cannot be lower or equal zero. Define it in SYSTEM.', faterr)
     ! Setup constant
     afact = epsm*voltage/(2.0*epsm+cdie*ikappa*tmemb) ! Eq. (32) paper

     ceps=cdie/epsm
  
     ampl1 = 0.0
     p1 = 1.0
     ampl2 = 0.0
     p2 = 1.0
     rcylinder = 1.0
   
     endlog = .false.
     do while (.not.endlog)
       call getlin(com,inpu,outu) ! new commands line
       endlog = check(com,'end')
       if (.not.endlog) then
         ! Obtention of ion type atnam(itype)      
         call getfirst(com,wrd4)
         call fatnam(atnam2,nttyp,wrd4,itype)
         ! Repulsive membrane potential [real*8,default=0]
         call gtdpar(com,'amplmemb',ampl1(itype),0.0)
         if (ampl1(itype).lt.0.0) call error ('shell_simul', 'Repulsive membrane potential has a negative value', faterr)
         ! Switching region between bulk and membrane regions
         ! [real*8,default=1]
         call gtdpar(com,'pmemb',p1(1,itype),1.0)
         if (p1(1,itype).le.0.0) call error ('shell_simul', 'Switching region between bulk and membrane regions has a negative or zero value', faterr)
         p1(2,itype)=1.0/p1(1,itype)
         if (Qpore) then 
           ! Repulsive potential for a cylindrical pore
           ! [real*8,default=0]
           call gtdpar(com,'amplpore',ampl2(itype),0.0)
           if (ampl2(itype).le.0.0) call error ('shell_simul', 'Repulsive cylindrical pore potential has a negative or zero value', faterr)
           ! Switching region between pore and membrane regions
           ! [real*8,default=1]
           call gtdpar(com,'ppore',p2(1,itype),1.0)
           if (p2(1,itype).le.0.0) call error ('shell_simul', 'Switching region between pore and membrane regions has a negative or zero value', faterr)
           p2(2,itype)=1.0/p2(1,itype)
           ! Radius of a cylindrical pore
           call gtdpar(com,'radi',rcylinder(itype),1.0)
           if (rcylinder(itype).le.0.0) call error ('shell_simul', 'Radius of a cylindrical pore has a negative or zero value', faterr)
         endif
       endif
     enddo
     Qmemb = .true.
     do i=1,nttyp
       if (ampl1(i).lt.0.0) then 
         write(*,'(6x,a)') 'ERROR: Invalid membrane constant for type '//atnam2(i)
         Qmemb=.false.
       endif
       if (ampl2(i).lt.0.0.and.Qpore) then
         write(*,'(6x,a)') 'ERROR: Invalid Pore constant for type '//atnam2(i)
         Qmemb=.false.
       endif
     enddo
     if (.not.Qmemb) stop
  
     write(outu,*)
     write(outu,'(6x,a)') 'MEMBRANE parameters:'
     write(outu,'(6x,a)') '--------------------'
     write(outu,'(6x,a,f12.3,a)') 'Membrane thickness ',tmemb,' [Angs]'
     write(outu,'(6x,a,f12.3,a)') 'Membrane center    ',zmemb,' [Angs]'
     write(outu,'(6x,a,f12.3,a)') 'Voltage            ',voltage*kcalmol/Coulomb,' [Volts]'
     write(outu,'(6x,a,f12.3)') 'Dielectric constant of the membrane ',epsm
     write(outu,'(6x,a,f12.3,a)') 'Ionic strength     ',ionstr,' [mol/L]' 
     write(outu,'(6x,a)') 'AMPLMEMB -> repulsive membrane potential [Kcal/mol]'
     write(outu,'(6x,a)') 'PMEMB -> switching region between bulk and membrane regions [Ang.]'
     if (Qpore) then  
       write(outu,'(6x,a)') 'PORE Enabled'
       write(outu,'(6x,a)') 'AMPLPORE -> repulsive potential for a cylindrical pore [Kcal/mol]'
       write(outu,'(6x,a)') 'PPORE -> switching region between pore and membrane regions [Ang.]'
       write(outu,'(6x,a)') 'RCYL -> radius of a cylindrical pore [Ang.]'
       write(outu,'(6x,a)') 'NAME---AMPLMEMB---PMEMB---AMPLPORE---PPORE---RCYL'   
     else
       write(outu,'(6x,a)') 'NAME---AMPLMEMB---PMEMB'   
     endif
     do i = 1, nttyp
       if (Qpore) then 
         write(outu,'(6x,a,1x,5f12.3)') atnam2(i),ampl1(i),p1(1,i),ampl2(i),p2(1,i),rcylinder(i)
       else
         write(outu,'(6x,a,1x,5f12.3)') atnam2(i),ampl1(i),p1(1,i)
       endif
     enddo
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'srpmf') then ! short-range ion-ion interactions
  !        ---------------
     if (.not.Qpar) call error ('shell_simul', 'SRPMF order is defined before PARTICLE order', faterr)
     if (Qefpott) call error ('shell_simul', 'SRPMF must be defined before EFPOT', faterr)
     ntc=(nttyp+ndna+1)*nion/2
     allocate (c0(ntc),c1(is),c2(ntc),c3(is),c4(ntc))
     c0=0.0
     c1=0.0
     c2=0.0
     c3=0.0
     c4=0.0
     call gtdpar(com,'rth',rth,8.0)
     srpx=rth-0.5
     srpk=6.90775527898/0.5
     srpy=0.001
     rth=rth**2
     endlog = .false.
     do while (.not. endlog)
       call getlin(com,inpu,outu) ! new commands line
       itype = 0 ! counter of ion types
       jtype = 0 ! counter of ion types
       endlog = check(com,'end')
       if (.not.endlog) then
         ! Obtention of ion type atnam(itype)
         call getfirst(com,wrd4)
         call fatnam(atnam2,nttyp,wrd4,itype)
         ! Obtention of ion type atnam(jtype)
         call getfirst(com,wrd4)
         call fatnam(atnam2,nttyp,wrd4,jtype)
         if (itype.eq.0.or.jtype.eq.0) call error ('shell_simul', 'A ion pair is necessary in SRPMF order', faterr)
         is=nindex(itype,jtype)
         ! c_0 parameter [real*8,default=0]
         call gtdpar(com,'c0',c0(is),0.0)
         ! c_1 parameter [real*8,default=0]
         call gtdpar(com,'c1',c1(is),0.0)
         ! c_2 parameter [real*8,default=0]
         call gtdpar(com,'c2',c2(is),0.0)
         ! c_3 parameter [real*8,default=0]
         call gtdpar(com,'c3',c3(is),0.0)
         ! c_4 parameter [real*8,default=0]
         call gtdpar(com,'c4',c4(is),0.0)
         ! Logical variable which indicates if total ion-ion 
         ! interaction energy file will be written  
         Qlsprmf = check(com,'file')
         if (Qlsprmf) then
           ! Unit number to write the total ion-ion interaction energy
           ! [integer,default=1]             
           call gtipar(com,'unit',iunit,1)
           if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
           if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
           if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in SRPMF order', faterr)
           iunit = unvec(iunit)
           is=nindex(itype,jtype)
           ! Parameters of the Lennard-Jones 6-12 potential 
           cc0 = c0(is)
           cc1 = c1(is)
           cc2 = c2(is)
           cc3 = c3(is)
           cc4 = c4(is)
           ! r1 -> Interaction distance. This functions is applied
           ! only up to r1=8.0 Angstrom (0.05 increment)
           do r1i = 2, 160
             r1=real(r1i)*0.05
             r2 = r1**2
             r6 = (sgp2(is)/r2)**3
             ! OBTENTION OF LENNARD-JONES 6-12 POTENTIAL
             evdw  = epp4(is)*r6*(r6-1.0)
             ! OBTENTION OF ELECTROSTATIC INTERACTION BETWEEN IONS
             eelec = fct(is)/r1
             ! OBTENTION OF WATER-MEDIATED SHORT-RANGE ION-ION
             ! INTERACTION  
             esrpmf = cc0*exp((cc1-r1)/cc2)*cos(cc3*pi*(cc1-r1))+cc4*(cc1/r1)**6
             if (evdw+eelec+esrpmf.le.50.0) write(iunit,'(4f12.5)') r1,evdw+eelec+esrpmf,evdw+eelec,esrpmf
           enddo
         endif
       endif
     enddo
     Qsrpmf = .true.
  
     write(outu,'(6x,a)')'Short-range ion-ion interactions are turn on' 
     write(outu,'(6x,a,/,6x,a)') 'Coefficients for the short-range ion-ion interactions','NAME----NAME----C0----C1----C2----C3----C4'  
     ! ndna = nttyp-ntype+nold
     do i = 1,nttyp
       do j=i,nttyp
         if(j.gt.ndna) then
           is=nindex(i,j)
           if (c2(is).ne.0.0) then
             write(outu,'(6x,a,1x,a,5f8.3)') atnam2(i),atnam2(j),c0(is),c1(is),c2(is),c3(is),c4(is)
             c2(is)=1.0/c2(is)
             Qsrpmfi(is)=.true.
           else
              write(outu,'(6x,a,x,a,x,a)') 'Warning: Missing SRPMF parameters for:',atnam2(i),atnam2(j)
           endif
         endif
       enddo
     enddo
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'efpot') then ! effective potential
  !       ---------------
    if (.not.(Qpar.or.Qnucl)) call error ('shell_simul','EFPOT order is defined before PARTICLE and/or NUCLEOTIDE order', faterr)
    cnt=0
    ntc=(nttyp+ndna+1)*nion/2 
    Qepwrt=check(com,'write')
    if (Qepwrt) then 
      call gtipar(com,'unit',wunit,1)
      wunit = unvec(wunit)
    endif
    call gtdpar(com,'scal',scald,kBT)
    call gtdpar(com,'res',res,0.1)
    ires=1.0/res
    if (Qsrpmf) maxl=sqrt(rth)
    call gtdpar(com,'maxdist',maxlg,maxl)
    mnp=int(maxlg*ires)
    call gtdpar(com,'mindist',minlg,2.5)
    nnp=int(minlg*ires)
    minlg=float(nnp)*res
    maxlg=float(mnp)*res
    maxd=int(float(ntc)*maxlg*ires)
    call gtipar(com,'maxdata',maxd,maxd)     
    write(outu,'(6x,a,i0)') 'Max number of data points to be stored (temporary storage array): ',maxd 
    allocate (nxi(ntc),nxf(ntc),xy(2,maxd),dmi(ntc),dm2(2,ntc),sc(3,ntc),Qefpotread(ntc))
    Qefpotread=.false.
    write(outu,'(6x,a)')
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      itype = 0 ! counter of ion types
      jtype = 0 ! counter of ion types
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type atnam(itype)
        call getfirst(com,wrd4)
        call fatnam(atnam2,nttyp,wrd4,itype)
        ! Obtention of ion type atnam(jtype)
        call getfirst(com,wrd4)
        call fatnam(atnam2,nttyp,wrd4,jtype)
        if (itype.eq.0.or.jtype.eq.0) call error ('shell_simul','A ion pair is necessary in EFPOT order', faterr)
        is=nindex(itype,jtype)
        if (check(com,'read')) then
          Qefpot(is) = .true.
          Qefpotread(is) = .true.
          ! Unit number to read the pair effective potential
          ! [integer,default=1]             
          call gtipar(com,'unit',iunit,1)
          if (iunit.le.0) call error('shell_simul','unit is zero or a negative number', faterr)
          if (iunit.gt.maxopen) call error('shell_simul','unit is greater than maxopen', faterr)
          if (unvec(iunit).eq.-1) call error('shell_simul','unit incorrect in EFPOT order', faterr)
          iunit = unvec(iunit)
          call gtdpar(com,'scal',scaldd,scald)
          write(outu,*) '     Scaling ',atnam2(itype),'-',atnam2(jtype),' potential by ',scaldd
          cnt=cnt+1
          if (cnt.gt.maxd) call error ('shell_simul','Number of data is greater than expected. Increase maxdata.',faterr) ! Check if maxdata is correct
          read(iunit,*,IOSTAT=kode) xy(1,cnt),xy(2,cnt)
          xy(2,cnt)=xy(2,cnt)*scaldd
          if (kode.eq.0) nxi(is)=cnt
          do while (kode.eq.0)
            ! Check if resolution is correct {
            if (cnt.ge.1+nxi(is)) then
              if (xy(1,cnt)-xy(1,cnt-1).lt.0.99*res.or.xy(1,cnt)-xy(1,cnt-1).gt.1.01*res) call error ('shell_simul','Unexpected x spacing. Check that all x data spacing is separated by res.',faterr)
            endif
            ! }
            cnt=cnt+1
            if (cnt.gt.maxd) call error ('shell_simul','Number of data is greater than expected. Increase maxdata.',faterr) ! Check if maxdata is correct
            read(iunit,*,IOSTAT=kode) xy(1,cnt),xy(2,cnt)
            xy(2,cnt)=xy(2,cnt)*scaldd
          enddo
          cnt=cnt-1
          nxf(is)=cnt
        elseif (check(com,'build')) then
          if (Qlj(is)) then 
            Qefpot(is)=.true.
            call gtdpar(com,'maxdist',maxl,maxlg)
            mnp=int(maxl*ires)
            call gtdpar(com,'mindist',minl,minlg)
            nnp=int(minl*ires)
            nxi(is)=cnt+1
            cnt=cnt+1+mnp-nnp
            nxf(is)=cnt
            dmi(is)=nnp*res
            dm2(1,is)=(nnp*res)**2
            dm2(2,is)=(mnp*res)**2
          else
            call error ('shell_simul','Cannot build potential, LJ parameters missing.',faterr)
          endif
        endif
      endif
    enddo
    allocate (ep(1:3,cnt))
    write(outu,'(6x,a)') 'Effective potential was activated for the following pairs:'
    do i = 1,nttyp
      do j=i,nttyp
        if(j.gt.ndna) then
          is=nindex(i,j)
          if (Qefpot(is).and.Qefpotread(is)) then
            write(outu,'(6x,a,1x,2a,i5,2(a,f8.3))')atnam2(i),atnam2(j),'  Number of Points:',1-nxi(is)+nxf(is),'  From-To: ',xy(1,nxi(is)),' - ',xy(1,nxf(is))
            call splinepot(is,1-nxi(is)+nxf(is),xy(1,nxi(is):nxf(is)),xy(2,nxi(is):nxf(is)),nxf(is))
          elseif (Qefpot(is).and..not.Qefpotread(is)) then
            nnp=int(dmi(is)*ires)
            mnp=nxf(is)-nxi(is)+nnp
            call discretize(is,nnp,mnp,nxf(is))
            write(outu,'(6x,a,1x,2a,i5,2(a,f8.3),a)')atnam2(i),atnam2(j),'  Number of Points:',1+mnp-nnp,'  From-To: ',dmi(is),' - ',mnp*res,'  (built from current parameters)'
          else
            write(outu,'(6x,a,x,a,x,a,a)')'WARNING: Missing Effective Potential for:',atnam2(i),atnam2(j),' . Using '  &
                      //' Coulombic and Lennard Jones and/or SRPMF parameters on the fly'
          endif
        endif
      enddo
    enddo
    write(outu,*)'     Resolution: ',res
    if (Qepwrt) then
      do i=1,nttyp
        do j=i,nttyp
          if (j.gt.ndna) then
            is=nindex(i,j)
            if (Qefpot(is)) then
              write(wunit,*) atnam2(i),' - ',atnam2(j)
              do k=int((dmi(is)-1.0)*ires*10.0),int((sqrt(dm2(2,is))+10.0)*ires*10.0)
                x1=k*res*0.1
                call getyd(is,x1**2,x2,x3,r1)
                write(wunit,*) x1,x2,x3
              enddo
            endif
          endif
        enddo
      enddo
    endif
    deallocate (xy,nxf,Qefpotread) 
    Qefpott=.true.
    write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'proxd') then ! user defined proximity diffusion
  !        -----------------     
    if (.not.(Qpar.and.Qnucl)) call error ('shell_simul','PROXDIFF need IONS and DNA',faterr)
    if (.not.(Qefpott)) call error ('shell_simul','PROXDIFF must be defined after EFPOT',faterr)
    call gtdpar(com,'beta',beta,2.93)
    call gtdpar(com,'diff0',diff0,0.101)
    call gtdpar(com,'diffcutoff',diffcutoff,maxl)
    ibeta=1.0/beta
    Qproxdiff=.true.
    write(outu,'(6x,a)') 'PROXDIFF ACTIVATED:'
    write(outu,'(6x,a)') '------------------'
    write(outu,'(6x,a,f8.3,a)') 'beta  =      ',beta,' [Angs]'
    write(outu,'(6x,a,f8.3,a)') 'diff0 =      ',diff0,' [Angs^2/ps]'
    write(outu,'(6x,a,f8.3,a)') 'diffcutoff = ',diffcutoff,' [Angs]'

  elseif (wrd5.eq.'chden') then ! user defined proximity diffusion
  !        -----------------    
    if (.not.Qsystem) call error ('shell_simul', 'CHDEN must be defined after SYSTEM', faterr)
    Qchdenorm = check(com,'norm')
    call gtdpar(com,'dcel',dcel4,0.5)
    call gtdpar(com,'xbcen',xbcen4,cx)
    call gtdpar(com,'ybcen',ybcen4,cy)
    call gtdpar(com,'zbcen',zbcen4,cz)
    nclx4=int(lx/dcel4)+1
    ncly4=int(ly/dcel4)+1
    nclz4=int(lz/dcel4)+1
    call gtipar(com,'nclx',nclx4,nclx4)
    call gtipar(com,'ncly',ncly4,ncly4)
    call gtipar(com,'nclz',nclz4,nclz4)
    tranx4 = 0.5*(nclx4-1)*dcel4
    trany4 = 0.5*(ncly4-1)*dcel4
    tranz4 = 0.5*(nclz4-1)*dcel4
    idcel4=1.0/dcel4
    allocate (chden(nclx4*ncly4*nclz4))
    chden=0.0
    Qchden=.true.
    write(outu,'(6x,a)') 'Charge Density (CHDEN) Activated:'
    write(outu,'(6x,a)') '------------------'
    write(outu,'(6x,a,i6)')  'Number of grid point in X   (nclx) = ',nclx4
    write(outu,'(6x,a,i6)')  'Number of grid point in Y   (ncly) = ',ncly4
    write(outu,'(6x,a,i6)')  'Number of grid point in Z   (nclz) = ',nclz4
    write(outu,'(6x,a,f8.3)') 'Grid spacing                (dcel) = ',dcel4
    write(outu,'(6x,a,f8.3)') 'Center of box in X          (xbcen)= ',xbcen4
    write(outu,'(6x,a,f8.3)') 'Center of box in Y          (ybcen)= ',ybcen4
    write(outu,'(6x,a,f8.3)') 'Center of box in Z          (zbcen)= ',zbcen4
    write(outu,'(6x,a,f8.3,a,f8.3)') 'Box in X from ',xbcen4-tranx4,' to ',xbcen4+tranx4
    write(outu,'(6x,a,f8.3,a,f8.3)') 'Box in Y from ',ybcen4-trany4,' to ',ybcen4+trany4
    write(outu,'(6x,a,f8.3,a,f8.3)') 'Box in Z from ',zbcen4-tranz4,' to ',zbcen4+tranz4
  ! **********************************************************************
  elseif (wrd5.eq.'profi') then ! user defined diffusion constant profile
  !        -----------------     
     if (.not.Qpar) call error ('shell_simul', 'PROFILE order is defined before PARTICLE order', faterr)
     ! Unit for file of diffusion constant profile      
     call gtipar(com,'difunit',iunit,0)
     if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
     if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
     if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PROFILE order', faterr)
     iunit = unvec(iunit)
     ! Reading file of diffusion constant profile
     read(iunit,*) nspline
     if (nspline.gt.dspline) call error ('shell_simul', 'nspline is greater than dspline', faterr)
     if (nspline.gt.2) then 
       Qprofile = .true.
       do is = 1, nspline
         read(iunit,*) xs(is),ys(is)
       enddo
       call spline (nspline,xs,ys,b,c,d)
       write(outu,*)'Diffusion coefficient profile activated'
       write(outu,*)'  nspline = ',nspline
       write(outu,*)'  z: defined from ',xs(1),' to ',xs(nspline)
     endif
  ! ****************************************************************
  elseif (wrd5.eq.'diffu') then ! diffusion_constant
  !        ---------------
     if (.not.Qpar) call error ('shell_simul', 'DIFFUSION order is defined before PARTICLE order', faterr)
     if (.not.Qsystem) call error ('shell_simul', 'DIFFUSION order is defined before SYSTEM order', faterr)
     ! pore length [real*8,default=0]
     call gtdpar(com,'porelength',plength2,0.0)
     if (plength2.lt.0.0) call error ('shell_simul', 'Pore lenght is lower than zero in DIFFUSION order', faterr)
     plength2=plength2*0.5  !store half the pore length
     ! membrane center [real*8,default=0]
     call gtdpar(com,'pcenter',pcenter,0.0)
     if (pcenter.gt.lz2p.or. pcenter.lt.lz2m) call error ('shell_simul', 'Position of center of a pore along the Z-axis is not correct', faterr)
  
     endlog = .false.
     do while (.not. endlog) 
       call getlin(com,inpu,outu) ! new commands line
       endlog = check(com,'end')
       if (.not.endlog) then
         ! Obtention of ion type atnam(itype)
         call getfirst(com,wrd4)
         call fatnam(atnam,ntype,wrd4,itype)
         call gtdpar(com,'lowfraction',ampl3(itype),0.0)
         call gtdpar(com,'switchlength',p3(itype),1.0)
         ! Writting file of diffusion constant 
         if (check(com,'file')) then
           call gtipar(com,'unit',iunit,1)
           if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
           if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
           if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in DIFFUSION order', faterr)
           iunit = unvec(iunit)
           do z1i = -90, 90, 1
             z1 = real(z1i)
             call switch3(r1,r2,z1,plength2,p3(itype),pcenter)
             write(iunit,'(2f13.5)') z1,diffusion(itype)*(ampl3(itype)+(1.0-ampl3(itype))*r1)
           enddo
         endif
       endif  
     enddo
     Qdiffuse = .true.
  
     write(outu,*)
     write(outu,'(6x,a)') 'DIFFUSION parameters:'
     write(outu,'(6x,a)') '---------------------'
     write(outu,'(6x,a,f8.3,a)') 'pore length        ',2*plength2,' [Angs]'
     write(outu,'(6x,a,f8.3,a)') 'membrane center    ',   pcenter,' [Angs]'
  
     do i = nold+1, ntype
        write(outu,'(6x,i3,3x,a,1x,5f8.3)') i,atnam(i),ampl3(i),p3(i)
     enddo
     write(outu,*)
  ! *******************************************************************
  elseif (wrd5.eq.'print') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'PRINT order is defined before PARTICLE and/or NUCLEOTIDE orders', warning)
     if (check(com,'system')) then
       if (Qnucl) write (outu,'(6x,a,i4)') 'Total number of sites ',nsites
       if (Qpar) then       
         nions = ntot - nsites
         write(outu,'(6x,a,i4)') 'Total number of ions  ',nions
         write(outu,'(6x,a,i4)') 'Number of fixed ions  ',nfix
         do ib = 1, nbuffer
           write(outu,'(6x,a,2i4)') 'buffer---number of ions ',ib,nat(ib)
         enddo
         write(outu,'(6x,a,i4)') 'Number of moving atoms ',natom
       endif
     elseif (check(com,'coor')) then
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)   
       write(outu,*)    
       write(outu,'(6x,a)') 'Configuration has been written'
       if (check(com,'xyz')) then
         write(iunit,'(i5)') ntot
         write(iunit,*)
         do i=1,nsites
           write(iunit,*) atnam2(typtyp(i)),x(i),y(i),z(i)
         enddo
         do i=nsites+1,ntot       
           write(iunit,*) atnam2(nwtype(abs(typei(i)))),x(i),y(i),z(i)
         enddo
       elseif (check(com,'pdb')) then
         if (.not.Qsphere) write(iunit,'(A6,3f9.3,3f7.2)') 'CRYST1',LX,LY,LZ,90.0,90.0,90.0
         do i=1,nsites
           write (iunit,'(A6,I5,x,A5,A5,I4,4x,3F8.3)') 'ATOM  ',i,atnam2(typtyp(i)),'DNA  ',typtyp(i),x(i),y(i),z(i)
         enddo
         do i=1+nsites,ntot
           write (iunit,'(A6,I5,x,A5,A5,I4,4x,3F8.3)') 'ATOM  ',i,atnam2(nwtype(abs(typei(i)))),'ION  ',nwtype(abs(typei(i))),x(i),y(i),z(i)
         enddo
         write (iunit,'(A)') 'END'
       else
         if (Qnucl) then
           write(iunit,'(i5)') nsites
           do i = 1, nsites
            write(iunit,'(i1,1x,a1,1x,a2,3(1x,f15.8))') strand(i),namnucl(i),namsite(i),x(i),y(i),z(i)
           enddo
         endif
         nions = ntot - nsites
         write(iunit,'(i5)') nions
         if (Qpar) then
           if (check(com,'charmm')) then
             do i = nsites+1, ntot    
               if (Qbuf) then
                 j=ibuffer(i)
               else
                 j=0
               endif
               write(iunit,'(2(i5,1x),2(a4,1x),3(f10.5,1x),i4)') i,i,atnam2(nwtype(abs(typei(i)))),atnam2(nwtype(abs(typei(i)))),x(i),y(i),z(i),j
             enddo
           else
             do i = nsites+1, ntot    
               if (Qbuf) then
                 j=ibuffer(i)
               else
                 j=0
               endif
               write(iunit,'(2(i5,1x),3(f15.8,1x),i4)') i,typei(i),x(i),y(i),z(i),j
             enddo
           endif
         endif
       endif  
     elseif (check(com,'title')) then
       write(outu,'(6x,a)') 'TITLE: '//trim(title)
     elseif (check(com,'iseed')) then
       write(outu,'(a,i10)') ' iseed: ',iseed
     elseif (check(com,'dnacenter')) then
       if (Qnucl) then
         xm=0.0
         ym=0.0
         zm=0.0
         do i=1,nsites
           xm=xm+x(i)
           ym=ym+y(i)
           zm=zm+z(i)
         enddo
         xm=xm/nsites
         ym=ym/nsites
         zm=zm/nsites
         write(outu,'(6x,a,3(f15.8,1x))') 'DNA GEOMETRIC CENTER: ',xm,ym,zm
       else 
         call error ('shell_simul', 'Cannot PRINT dnacenter, DNA is not defined', warning)
       endif
     elseif (check(com,'static')) then
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'res',resol,0.1)
       call gtdpar(com,'x1',x1,0.0) 
       call gtdpar(com,'y1',y1,0.0) 
       call gtdpar(com,'z1',z1,0.0) 
       call gtdpar(com,'x2',x2,0.0) 
       call gtdpar(com,'y2',y2,0.0) 
       call gtdpar(com,'z2',z2,0.0) 
       call staticplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead)
     elseif (check(com,'repul')) then
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'res',resol,0.1)
       call gtdpar(com,'x1',x1,0.0)
       call gtdpar(com,'y1',y1,0.0)
       call gtdpar(com,'z1',z1,0.0)
       call gtdpar(com,'x2',x2,0.0)
       call gtdpar(com,'y2',y2,0.0)
       call gtdpar(com,'z2',z2,0.0)
       call repulplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead)
     elseif (check(com,'rfpar')) then
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'res',resol,0.1)
       call gtdpar(com,'x1',x1,0.0)
       call gtdpar(com,'y1',y1,0.0)
       call gtdpar(com,'z1',z1,0.0)
       call gtdpar(com,'x2',x2,0.0)
       call gtdpar(com,'y2',y2,0.0)
       call gtdpar(com,'z2',z2,0.0)
       if (Qrfpsin) then 
         call rfparplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead,1)
       else
         call rfparplot(x1,y1,z1,x2,y2,z2,resol,iunit,Qnohead,nion)
       endif
     elseif (check(com,'statxd')) then
!  Prints static field values in x,y,z,pot format in an ASCII file. Allows to generate 3D,2D,1D Plots.
!      [PRINT] statxd unit (integer) {noheader} ix1 (integer) ix2 (integer) iy1 (integer) iy2 (integer) iz1 (integer) iz2 (integer)
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtipar(com,'ix1',ix1,0)
       call gtipar(com,'iy1',iy1,0)
       call gtipar(com,'iz1',iz1,0)
       call gtipar(com,'ix2',ix2,nclx1-1)
       call gtipar(com,'iy2',iy2,ncly1-1)
       call gtipar(com,'iz2',iz2,nclz1-1)
       call statxd(ix1,iy1,iz1,ix2,iy2,iz2,iunit,Qnohead)
     elseif (check(com,'repxd')) then
!  Prints repulsion field values in x,y,z,rep format in an ASCII file. Allows to generate 3D,2D,1D Plots.
!      [PRINT] repxd unit (integer) {noheader} ix1 (integer) ix2 (integer) iy1 (integer) iy2 (integer) iz1 (integer) iz2 (integer)
       Qnohead=check(com,'noheader')
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtipar(com,'ix1',ix1,0)
       call gtipar(com,'iy1',iy1,0)
       call gtipar(com,'iz1',iz1,0)
       call gtipar(com,'ix2',ix2,nclx2-1)
       call gtipar(com,'iy2',iy2,ncly2-1)
       call gtipar(com,'iz2',iz2,nclz2-1)
       call repxd(ix1,iy1,iz1,ix2,iy2,iz2,iunit,Qnohead)
     elseif (check(com,'chden')) then
       if (Qchden) then
         call gtipar(com,'unit',iunit,1)
         if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
         if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
         if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
         iunit = unvec(iunit)
         zero=0d0
         idv=sng(1.0/dcel4**3)
         ncl3=nclx4*ncly4*nclz4
         write(iunit) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunit) (zero,i=1,6)
         write(iunit) ((idv*chden(i)),i=1,ncl3)
         write(outu,'(6x,A)') 'Charge density map written'
       else
         call error ('shell_simul', 'nothing to print, CHDEN was not called', warning) 
       endif
     elseif (check(com,'efield')) then
       if (allocated(iunitv)) deallocate (iunitv)
       allocate (iunitv(4))
       Qonlychden=check(com,'onlychden')
       call gtipar(com,'xunit',iunitv(1),0)
       if (iunitv(1).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(1).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(1)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(1) = unvec(iunitv(1))
       call gtipar(com,'yunit',iunitv(2),0)
       if (iunitv(2).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(2).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(2)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(2) = unvec(iunitv(2))
       call gtipar(com,'zunit',iunitv(3),0)
       if (iunitv(3).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(3).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(3)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(3) = unvec(iunitv(3))
       call gtipar(com,'munit',iunitv(4),0)
       if (iunitv(4).le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunitv(4).gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunitv(4)).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunitv(4) = unvec(iunitv(4))
       if (Qchden) then
         ncl3=nclx4*ncly4*nclz4
         allocate(efield(4,ncl3))
         efield=0.0
         zero=0d0
         if (Qchdencnt) then 
           do i=0,nclx4-1
             do j=0,ncly4-1
               do k=0,nclz4-1
                 x1=i*dcel4-trany4+xbcen4 
                 y1=j*dcel4-trany4+ybcen4 
                 z1=k*dcel4-tranz4+zbcen4
                 in1=i*ncly4*nclz4+j*nclz4+k+1
                 do ii=0,nclx4-1
                   do ij=0,ncly4-1
                     do ik=0,nclz4-1
                       x2=ii*dcel4-trany4+xbcen4
                       y2=ij*dcel4-trany4+ybcen4 
                       z2=ik*dcel4-tranz4+zbcen4
                       in2=ii*ncly4*nclz4+ij*nclz4+ik+1
                       if (in1.ne.in2.and.chden(in2).ne.0.0) then
                         r2=1.0/((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                         r1=sqrt(r2)
                         efield(1,in1)=efield(1,in1)+sng(celec2*chden(in2)*(x2-x1)*r1*r2)
                         efield(2,in1)=efield(2,in1)+sng(celec2*chden(in2)*(y2-y1)*r1*r2)
                         efield(3,in1)=efield(3,in1)+sng(celec2*chden(in2)*(z2-z1)*r1*r2)
                         efield(4,in1)=efield(4,in1)+sng(celec2*chden(in2)*r2)
                       endif
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
         endif
         if (.not.Qonlychden) then
           if (Qphix) then
             do i=0,nclx4-1
               do j=0,ncly4-1
                 do k=0,nclz4-1
                   x1=i*dcel4-trany4+xbcen4
                   y1=j*dcel4-trany4+ybcen4
                   z1=k*dcel4-tranz4+zbcen4
                   in1=i*ncly4*nclz4+j*nclz4+k+1
                   call staefield(x1,y1,z1,x2,y2,z2)
                   efield(1,in1)=efield(1,in1)+sng(x2)
                   efield(2,in1)=efield(2,in1)+sng(y2)
                   efield(3,in1)=efield(3,in1)+sng(z2)
                   efield(4,in1)=efield(4,in1)+sng(sqrt(x2**2+y2**2+z2**2))
                 enddo
               enddo
             enddo
           elseif (voltage.ne.0.0.and.Qmemb) then
             x2=afact/Coulomb*kcalmol
             do k=0,nclz4-1
               z1=k*dcel4-tranz4+zbcen4
               if (z1.lt.zmemb1) then ! REGION 1: z=z(iat)-zmemb1 < 0, lim{z->-inf} pot(1) = 0
                 z2 = -x2*ikappa*exp(ikappa*(z1-zmemb1))
               elseif ((z1.ge.zmemb1) .and. (z1.le.zmemb2)) then ! REGION 2
                 z2 = -x2*ceps*ikappa
               elseif (z1.gt.zmemb2) then ! REGION 3: z=z(iat)-zmemb2 > 0, lim{z->inf} pot(2) = voltage
                 z2 = -x2*ikappa*exp(-ikappa*(z1-zmemb2))
               endif
               do i=0,nclx4-1
                 do j=0,ncly4-1
                   x1=i*dcel4-trany4+xbcen4
                   y1=j*dcel4-trany4+ybcen4
                   in1=i*ncly4*nclz4+j*nclz4+k+1
                   efield(3,in1)=efield(3,in1)+sng(z2)
                   efield(4,in1)=efield(4,in1)+sng(z2)
                 enddo
               enddo
             enddo
           endif
         endif
         write(iunitv(1)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(1)) (zero,i=1,6)
         write(iunitv(1)) ((efield(1,i)),i=1,ncl3)
         write(iunitv(2)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(2)) (zero,i=1,6)
         write(iunitv(2)) ((efield(2,i)),i=1,ncl3)
         write(iunitv(3)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(3)) (zero,i=1,6)
         write(iunitv(3)) ((efield(3,i)),i=1,ncl3)
         write(iunitv(4)) nclx4,ncly4,nclz4,dcel4,xbcen4,ybcen4,zbcen4
         write(iunitv(4)) (zero,i=1,6)
         write(iunitv(4)) ((efield(4,i)),i=1,ncl3)
         write(outu,'(6x,A)') 'Electric Field maps written'
         deallocate(efield)
       else
         call error('shell_simul', 'nothing to print, CHDEN was not called', warning)
       endif
       deallocate (iunitv)
     elseif (check(com,'statsup')) then
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in PRINT order', faterr)
       iunit = unvec(iunit)
       call gtdpar(com,'x',x1,0.0)
       call gtdpar(com,'y',y1,0.0)
       call areapot(x1,y1,iunit)
     endif
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'count') then
  !        ---------------
     if (.not.Qbuf) call error ('shell_simul', 'COUNT order is defined before BUFFER order', faterr)
     call COUNT
     write(outu,*)
     write(outu,'(6x,a,i4)') 'Total number of ions ',ntot-nsites
     do ib = 1, nbuffer
       write(outu,'(6x,a,2i4)') 'buffer--number of ions ',ib,nat(ib)
     enddo
     write(outu,'(6x,a)') 'ION----TYPE--BUFFER'
     do i = nsites+1, ntot
       itype = abs(typei(i))
       if (ibuffer(i).ne.0) write(outu,'(6x,i5,1x,a4,1x,i5)') i,atnam(itype),ibuffer(i)
     enddo
  ! **********************************************************************
  elseif (wrd5.eq.'coor') then
  !        ---------------
     if (.not.Qpar .and. .not.Qnucl) call error ('shell_simul', 'COOR order is defined before PARTICLE and/or NUCLEOTIDE orders', faterr)
     if (check(com,'fixion').and.Qpar) then
       call gtipar(com,'nfix',nfix,nfix)
       if (nfix.lt.0) call error ('shell_simul', 'nfix is lower than zero',faterr)
       if (nfix.gt.0) then
         do i = nsites+1, nsites+nfix
           call getlin(com,inpu,outu) ! new commands line
           ! Obtention of ion type atnam(itype)
           call getfirst(com,wrd4)
           call fatnam(atnam,ntype,wrd4,itype)
           typei(i) = itype
           ! X-coordinate for fixed ion
           call gtdpar(com,'xcor',x(i),0.0)
           ! Y-coordinate for fixed ion
           call gtdpar(com,'ycor',y(i),0.0)
           ! Z-coordinate for fixed ion
           call gtdpar(com,'zcor',z(i),0.0)
         enddo
         if (.not.setline(inpu,com)) call error ('shell_simul', 'premature end of date in COOR order', faterr)
         if (.not.setword(word,com)) call error ('shell_simul', 'premature end of date in COOR order', faterr)
         if (lcase(word).ne.'end') call error ('shell_simul', 'END keyword is not included in COOR order', faterr)
       endif
     elseif (check(com,'read')) then
       call gtipar(com,'unit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in COOR order', faterr)
       iunit = unvec(iunit)
       write(outu,'(6x,a,i3)')'Reading coordinates from unit ',iunit
       if (Qnucl) then
         read(iunit,*) ii
         if (ii.ne.nsites) call error ('shell_simul', 'NSITES has a wrong value',faterr)
         write(outu,'(6x,a,i5)') 'nsites ',nsites
         do i = 1, nsites
           read(iunit,'(a)') com
! free format
           read(com,*) strand(i), namnucl(i),namsite(i), x(i), y(i), z(i)
!           read(com,'(i1,1x,a1,1x,a2,1x,f10.5,1x,f10.5,1x,f10.5)') strand(i), namnucl(i),namsite(i), x(i), y(i), z(i)
         enddo
         ntot = nsites
       endif
       if (Qpar) then 
         if (check(com,'charmm')) then
           ! CHARMM format     
           endlog = .true.
           do while (endlog)
             read(iunit,'(a)') com
             endlog = com(1:1).eq.'*'
           enddo 
           read(com,'(i5)') nions
           write(outu,'(6x,a,i5)') 'nions ',nions
           ntot = nsites + nions
           if (ntot.gt.datom) call error ('shell_simul', 'ntot is greater than datom', faterr)
           do i = nsites+1, ntot
             read(iunit,'(a)') com
             read(com,*) ii,ii,ionname,ionname,x(i),y(i),z(i),ibuffer(i)
!             read(com,'(i5,1x,i5,1x,a4,1x,a4,1x,f10.5,1x,f10.5,1x,f10.5,1x,i4)') ii,ii,ionname,ionname,x(i),y(i),z(i),ibuffer(i)
             if (ii.ne.i) call error ('shell_simul', 'incorrect number of ions in COOR order', faterr)
             call fatnam(atnam,ntype,ionname,itype)
             typei(i)=itype
           enddo
         else
           read(iunit,'(i5)') nions
           write(outu,'(6x,a,i5)') 'nions ',nions
           ntot = nsites + nions
           if (ntot.gt.datom) call error ('shell_simul', 'ntot is greater than datom', faterr)
           do i = nsites+1, ntot
             read(iunit,'(a)') com
             read(com,*) ii,typei(i),x(i),y(i),z(i),ibuffer(i)
!             read(com,'(i5,1x,i5,1x,f10.5,1x,f10.5,1x,f10.5,1x,i4)') ii,typei(i),x(i),y(i),z(i),ibuffer(i)
             if (ii.ne.i) call error ('shell_simul', 'incorrect number of ions in COOR order', faterr)
           enddo
         endif
         call gtipar(com,'nfix',nfix,nfix)
         if (nfix.lt.0) call error ('shell_simul', 'nfix is lower than zero',faterr)
         natom = ntot - nfix
       endif  
       write(outu,'(6x,a)') 'coordinates have been read'
       call COUNT
     elseif (check(com,'gener')) then
       dodna=.false.
       doions=.false.
       if (check(com,'all')) then
         doions=Qbuf
         dodna=Qnucl
       else
         if (check(com,'dna')) dodna=Qnucl
         if (check(com,'ions')) doions=Qbuf
       endif
       if (.not.doions.and..not.dodna) call error ('shell_simul', 'all, dna or ions missing after COOR gener', faterr)
       if (dodna) then
         do i = 1, nsites
           x(i) = xnat(i)
           y(i) = ynat(i)
           z(i) = znat(i)
         enddo
         ntot = nsites
         if (ntot.gt.datom) call error ('shell_simul', 'ntot is greater than datom', faterr)
       endif
       if (doions) then
         do ib = 1, nbuffer
           nat(ib) = nint(avnum(ib))
         enddo
         natom = nsites ! initializations
         ifirst = nsites + nfix + 1 ! first position
         ilast  = nsites + nfix + nat(1) ! last position
         if (ilast.gt.datom) call error ('shell_simul', 'numbers of ions/sites is greater than datom', faterr)
         do ib = 1, nbuffer
           itype = ibfftyp(ib)
           natom = natom + nat(ib)
           do i = ifirst, ilast
             typei(i) = itype
             ibuffer(i) = ib
             call INSERT(ib,x(i),y(i),z(i))
             if(z(i).lt.cz) typei(i) = -itype
           enddo
           ifirst = ifirst + nat(ib)
           ilast = ifirst + nat(ib+1) - 1
           if (ilast.gt.datom) then
             call error ('shell_simul', 'numbers of ions/sites is greater than datom', faterr)
           endif
         enddo
         nions = natom - nsites + nfix
         ntot = nsites + nions
         if (ntot.gt.datom) call error ('shell_simul', 'ntot is greater than datom', faterr)
       endif ! Qbuf
       call gtipar(com,'nfix',nfix,nfix)
       if (nfix.lt.0) call error ('shell_simul', 'nfix is lower than zero',faterr)
       call count       
       write(outu,'(6x,a)') 'coordinates have been generated'
     elseif (check(com,'rot').and.Qnucl) then
       ! rotation for DNA sites coordinates
       Qrot = .true.
       if (check(com,'ref')) then
         call gtipar(com,'s1',s1,1)
         call gtipar(com,'s2',s2,nsites)
         call gtipar(com,'s3',s3,2)
         write(outu,'(6x,3(a,x,i0,x))') 'Rotating using references: s1=',s1,'s2=',s2,'s3=',s3
         vc3(1)=x(s1)-x(s2)
         vc3(2)=y(s1)-y(s2)
         vc3(3)=z(s1)-z(s2)
         x3=sqrt(dot_product(vc3,vc3))
         vc3=vc3/x3
         vc1(1)=x(s3)-x(s2)
         vc1(2)=y(s3)-y(s2)
         vc1(3)=z(s3)-z(s2)
         x1=sqrt(dot_product(vc1,vc1))
         vc1=vc1/x1
         call cross_product(vc3,vc1,vc2)
         x2=sqrt(dot_product(vc2,vc2))
         vc2=vc2/x2
         call cross_product(vc2,vc3,vc1)
         x1=sqrt(dot_product(vc1,vc1))
         vc1=vc1/x1
         rot(1,:)=vc1 
         rot(2,:)=vc2
         rot(3,:)=vc3 
       else
       ! rotation matrix elements
         call gtdpar(com,'r11',rot(1,1),1.0)
         call gtdpar(com,'r12',rot(1,2),0.0)
         call gtdpar(com,'r13',rot(1,3),0.0)
         call gtdpar(com,'r21',rot(2,1),0.0)
         call gtdpar(com,'r22',rot(2,2),1.0)
         call gtdpar(com,'r23',rot(2,3),0.0)
         call gtdpar(com,'r31',rot(3,1),0.0)
         call gtdpar(com,'r32',rot(3,2),0.0)
         call gtdpar(com,'r33',rot(3,3),1.0)
       endif
       do i = 1, 3
         sum1 = 0.0
         do j = 1, 3
           sum1 = sum1 + rot(i,j)*rot(i,j)
         enddo
         if (abs(1.0-sum1).ge.1.0e-8) call error ('shel_simul', 'rotation matrix for DNA sites is not orthogonal in COOR order', faterr)
       enddo
       sum1 = 0.0e0
       sum2 = 0.0e0
       sum3 = 0.0e0
       do j = 1, 3
         sum1 = sum1 + rot(1,j)*rot(2,j)
         sum2 = sum2 + rot(2,j)*rot(3,j)
         sum3 = sum3 + rot(1,j)*rot(3,j)
       enddo
       if (abs(sum1).ge.1.0e-8 .or. abs(sum2).ge.1.0e-8 .or.abs(sum3).ge.1.0e-8) call error ('shel_simul', 'rotation matrix for DNA sites is not orthogonal in COOR order', faterr)
       do i = 1, nsites
         xold = x(i)
         yold = y(i)
         zold = z(i)
         x(i) = xold*rot(1,1) + yold*rot(1,2) + zold*rot(1,3)
         y(i) = xold*rot(2,1) + yold*rot(2,2) + zold*rot(2,3)
         z(i) = xold*rot(3,1) + yold*rot(3,2) + zold*rot(3,3)
       enddo
       write(outu,'(6x,a)') 'Coordinates for DNA sites have been rotated. Rotation matrix:'
       do i = 1, 3
         write (outu,'(2x,3(1x,f8.3))') (rot(i,j),j=1,3)
       enddo
       ! translocation for DNA sites coordinates
     elseif (check(com,'tras').and.Qnucl) then
       Qtras=.true.
       ! traslation in x direction of DNA B isoform coordinates
       call gtdpar(com,'xtras',xtras,0.0)
       ! traslation in y direction of DNA B isoform coordinate
       call gtdpar(com,'ytras',ytras,0.0)
       ! traslation in z direction of DNA B isoform coordinate 
       call gtdpar(com,'ztras',ztras,0.0)
       do i = 1, nsites
         x(i) = x(i) + xtras
         y(i) = y(i) + ytras
         z(i) = z(i) + ztras
       enddo
       write(outu,'(6x,a,1x,a,3(f12.6,a))')'Coordinates for DNA sites have been moved by','(',xtras,',',ytras,',',ztras,')'
     elseif (check(com,'setori').and.Qnucl) then
       xm=0.0
       ym=0.0
       zm=0.0
       do i=1,nsites
         xm=xm+x(i)
         ym=ym+y(i)
         zm=zm+z(i)
       enddo
       xm=xm/nsites
       ym=ym/nsites
       zm=zm/nsites
       do i = 1, nsites
         x(i) = x(i) - xm
         y(i) = y(i) - ym
         z(i) = z(i) - zm
       enddo
       call gtdpar(com,'x',xm,xm)
       call gtdpar(com,'y',ym,ym)
       call gtdpar(com,'z',zm,zm)
       do i = 1, nsites
         x(i) = x(i) + xm
         y(i) = y(i) + ym
         z(i) = z(i) + zm
       enddo
       write(outu,'(6x,a,3(f12.6,a))')'DNA geometric center is now at  (',xm,',',ym,',',zm,')'
     endif
     if (Qnucl) then
       write(outu,'(6x,a)') 'Native structure for B isoform of DNA'
       write(outu,*)
       write(outu,'(6x,a,1x,a,1x,a,1x,a)') 'STRAND','NUCLEOT','SITE','CARTESIAN COORDINATES (X,Y,Z)'
       write(outu,'(6x,a,1x,a,1x,a,1x,a)') '------','-------','----','-----------------------------'
       do i = 1, nsites
        write(outu,'(6x,i1,1x,a1,1x,a2,3(1x,f15.8))') strand(i), namnucl(i),namsite(i), x(i), y(i), z(i)
       enddo
     endif
     if (Qpar) then
       if (nfix.gt.0) then
         write(outu,'(6x,a)') 'Fixed Ions'
         write(outu,'(6x,a,1x,a)') 'ION','CARTESIAN COORDINATES (X,Y,Z)'
         write(outu,'(6x,a,1x,a)') '---','-----------------------------'
         do i = nsites+1, nsites+nfix
           write(outu,'(6x,i5,1x,i5,3(1x,f15.8),1x,i4)') i,typei(i),x(i),y(i),z(i),ibuffer(i)
         enddo
       endif
       write(outu,'(6x,a)') 'Free Ions'
       write(outu,'(6x,a,1x,a)') 'ION','CARTESIAN COORDINATES (X,Y,Z)'
       write(outu,'(6x,a,1x,a)') '---','-----------------------------'
       do i = nsites+nfix+1, ntot
         write(outu,'(6x,i5,1x,i5,3(1x,f15.8),1x,i4)') i,typei(i),x(i),y(i),z(i),ibuffer(i)
       enddo
     endif
     write(outu,*)
  ! **********************************************************************
  elseif (wrd5.eq.'test') then
  !        ---------------
     if (Qnucl) call error ('shell_simul', 'DNA sites cannot uses in TEST order', faterr)
     if (.not.Qpar) call error ('shell_simul', 'Particles have to be defined before TEST order', faterr)
     if (nfix.gt.0) call error ('shell_simul', 'nfix > 0 in TEST order',faterr)
     ! tolerance in force test
     call gtdpar(com,'tol',tol,0.001)
     if (tol.le.0.0) call error ('shell_simul', 'Tolerance is negative or null in TEST order', faterr)
     call gtdpar(com,'delta',delta,0.0005)
     if (delta.le.0.0) call error ('shell_simul', 'Time step is negative or null in TEST order', faterr)
     Qnobond  = .not.check(com,'nobond')
     Qnonbond = .not.check(com,'nononbond')
     if (check(com,'membrane')) then
       if (.not.Qmemb) call error ('shell_simul', 'Planar membrane has not be defined before ENERGY order', faterr)
     else
       logmemb = Qmemb
       Qmemb = .false.
     endif
     if (check(com,'mmij')) then
       if (.not.Qmmij) call error ('shell_simul', 'Reaction Field has not be defined before TEST order', faterr)
     else
       logmmij = Qmmij
       Qmmij = .false.
     endif
     if (check(com,'phix')) then
       if (.not.Qphix) call error ('shell_simul', 'Static Field has not be defined before TEST order', faterr)
     else
      logphix = Qphix
      Qphix = .false.
     endif
     if (check(com,'phiv')) then
       if (.not.Qphiv) call error ('shell_simul', 'Repulsive term has not be defined before TEST order', faterr)
     else
       logphiv = Qphiv
       Qphiv = .false.
     endif
     if (check(com,'srpmf')) then
       if (.not.Qsrpmf) call error ('shell_simul', 'Short-range interaction term has not be defined before TEST order', faterr)
     else
       logsrpmf = Qsrpmf
       Qsrpmf = .false.
     endif
     if (check(com,'rfpar')) then
       if (.not.Qrfpar) call error ('shell_simul', 'Reaction Field Parameter term has not be defined before TEST order', faterr)
     else
       logrfpar = Qrfpar
       Qrfpar = .false.
     endif
     write(outu,*)
     write(outu,'(6x,a)') 'TEST calculation'
     write(outu,'(6x,a)') '----------------'
     if (Qnobond)  write(outu,'(6x,a)') 'Bonding Energy Term'
     if (Qnonbond) write(outu,'(6x,a)') 'Nonbonding Energy Term'
     if (Qmemb)    write(outu,'(6x,a)') 'Planar membrane Term'
     if (Qmmij)    write(outu,'(6x,a)') 'Reaction Field Energy Term'
     if (Qphix)    write(outu,'(6x,a)') 'External Field Energy Term'
     if (Qphiv)    write(outu,'(6x,a)') 'Repulsive Energy Term'
     if (Qsrpmf)   write(outu,'(6x,a)') 'Short-range Interaction Term'
     if (Qrfpar)   write(outu,'(6x,a)') 'Reaction Field Parameter Term'
     write(outu,'(6x,a,1x,e11.4)') 'Tolerance in force test',tol
     write(outu,'(6x,a,1x,e11.4)') 'Time step in force test',delta
  
     call energy
     write(outu,'(6x,a,f12.6)') 'Total energy (from ENERGY)      ',ener
     ener = 0.0 ! initialization
     do iat = 1, ntot
       call interact(dener,x(iat),y(iat),z(iat),abs(typei(iat)),iat,.true.)
       ener = ener + dener
     enddo
     ener = ener*0.5
     write(outu,'(6x,a,f12.6)') 'Total energy (from INTERACTION) ',ener
     call testfirst(tol,delta)         
  ! **************************************************************************
  elseif (wrd5.eq.'gsbp') then !  generalized solvent boundary potential
  !        ---------------
  !REPULSIVE POTENTIAL
     Qadj = check(com,'adjust')
     Qphiv = check(com,'phiv')
     if (Qphiv) then
       ! Threshold for 27 and 8 cells
       call gtdpar(com,'thold27',thold27,27.0)
       call gtdpar(com,'thold8',thold8,8.0)
       ! Magnitude of grid-based repulsive potential 'phiv'
       ! [real*8,default=0] 
       call gtdpar(com,'svdw',svdw,50.0)
       ! Trilinear interpolation is used for 'phiv'
       ! [default=3rd-order B-spline interpolation]
!       Qtrln = check(com,'trilinear')
       Qtrln = .not.check(com,'bspline')
       if (Qtrln) then
          write(outu,'(6x,a)') 'Trilinear function will be used for repulsive potential'
       else
          write(outu,'(6x,a)') 'B-spline function will be used for repulsive potential'
       endif
       Qnmcden = check(com,'nmcden')
       if (Qnmcden) then
         write(outu,*)
         write(outu,'(6x,a)') 'Different ion-accessible space is used for different ions and sites'
       endif   
       call gtipar(com,'phivunit',iunit,1)
       if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
       if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
       if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
       iunit = unvec(iunit)
       if (Qnmcden) then
         if (Qnucl .and. Qpar) then
           if (istrs.eq.1) numb = inuc
           if (istrs.eq.2) numb = 2*inuc
           totnumb = ntype - numb + 6    
         else if (Qnucl) then
           totnumb = 6
         else if (Qpar) then
           totnumb = ntype
         else
           call error ('shell_simul', 'GSBP order is defined before PARTICLE and/or NUCLEOTIDES orders', faterr)
         endif
         if (totnumb.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
         do i = 1, totnumb
           call gtipar(com,'munit',iunit,0)
           if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
           if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
           if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
           vecphiv(i) = unvec(iunit)
         enddo
         write(outu,'(6x,a,10i3)') 'Reading grid-based repulsive potential from unit ',(vecphiv(i),i=1,totnumb)
       else
         write(outu,*)
         write(outu,'(6x,a,i3)')'Reading grid-based repulsive potential from unit ',iunit
       endif
       call readphi(iunit,outu,'PHIV',Qadj)
       ! Maximum position of a grid for 'phiv' along the Z-axis
       ! [real*8,default=max map limit]
       call gtdpar(com,'vzmax',vzmax,zbcen2+tranz2)
       if (vzmax.gt.zbcen2+tranz2) call error('shell_simul', 'vzmax cannot be outside boundaries of repulsion map', faterr)
       ! Minimum position of a grid for 'phiv' along the Z-axis
       ! [real*8,default=-min map limit]         
       call gtdpar(com,'vzmin',vzmin,zbcen2-tranz2)
       if (vzmin.lt.zbcen2-tranz2) call error('shell_simul', 'vzmin cannot be outside boundaries of repulsion map', faterr)
       write(outu,'(6x,A,F8.2,A,F8.2,A,F8.2)') 'Uniformed repusive potential will be scaled by ',svdw,' kcal/mol between ',vzmin,' and ',vzmax 
       call gtipar(com,'repwalls',wallsi4,126)
       if (wallsi4.gt.0) then
         if (wallsi4.gt.126) then
           walls=126
         else
           walls=itoi1(wallsi4)
         endif
         call repwalls(walls)
         write(outu,'(6x,a,i0)') 'REPWALLS activated. Number of walls defined: ',walls
       else
         walls=0
         write(outu,'(6x,a,i0)') 'REPWALLS deactivated. Strongly recommended to be used.'
       endif
     endif ! Qhiv
    !REACTION FIELD
    Qmmij = check(com,'mmij')
    if (Qmmij) then 
     !RECTBOX => RECTANGULAR BOX FOR REACTION FIELD
     ! Maximum (minimum) positions in XYZ of a region where Legrende
     ! polynomial basis functions are applied for the reaction field 
     ! calculation [default=0]  
      call gtdpar(com,'xmax',rbxmax,0.0)
      call gtdpar(com,'xmin',rbxmin,0.0)
      call gtdpar(com,'ymax',rbymax,0.0)
      call gtdpar(com,'ymin',rbymin,0.0)
      call gtdpar(com,'zmax',rbzmax,0.0)
      call gtdpar(com,'zmin',rbzmin,0.0)
      if ((rbxmax-rbxmin).gt.0.0 .and. (rbymax-rbymin).gt.0.0.and. (rbzmax-rbzmin).gt.0.0) then    
        xscal = 2.0/(rbxmax-rbxmin)
        yscal = 2.0/(rbymax-rbymin)
        zscal = 2.0/(rbzmax-rbzmin)
        xmin=-1.0/xscal
        xmax= 1.0/xscal
        ymin=-1.0/yscal
        ymax= 1.0/yscal
        zmin=-1.0/zscal
        zmax= 1.0/zscal
        write(outu,'(6x,a)') 'Reaction field will be calculated following region;'
        write(outu,'(6x,2(a,f10.5))') 'X from ',rbxmin,' to ',rbxmax
        write(outu,'(6x,2(a,f10.5))') 'Y from ',rbymin,' to ',rbymax
        write(outu,'(6x,2(a,f10.5))') 'Z from ',rbzmin,' to ',rbzmax
        call gtdpar(com,'rfscal',rfscal,1.0)
        if(rfscal.ne.1.0) write(outu,'(6x,f5.2,a)') rfscal,' scaling factor will be applied for MMIJ rxnfld energy'
      endif
      !SPHERE => SPHERE FOR REACTION FIELD
      call gtdpar(com,'srdist',srdist,0.0)
      if (srdist.lt.0.0) then
        call error ('shell_simul', 'srdist is a negative number',warning)
        srdist = abs(srdist)
      endif
      if (srdist.ne.0.0) then
        write(outu,'(6x,a)') 'Reaction field will be calculated following region;'
        write(outu,'(6x,a,f10.5)') 'A sphere of radius ',srdist
      endif
      call gtipar(com,'mmijunit',iunit,1)
      if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
      if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
      if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
      iunit = unvec(iunit)
      write(outu,*)
      write(outu,'(6x,a,i3)') 'Reading MMIJ matrix from unit ',iunit
      call READMMIJ(iunit,outu)
      if (shapes.EQ.'RECTBOX ') then
        call RECT_NORM(ntpol,xscal,yscal,zscal,lstpx,lstpy,lstpz,bnorm)
      else
        call SPHE_NORM(nmpol,bnorm,srdist)
      endif
    endif ! Qmmij
    !STATIC EXTERNAL FIELD
    Qphix = check(com,'phix')
    if (Qphix) then
      call gtipar(com,'phixunit',iunit,1)
      if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
      if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
      if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
      iunit = unvec(iunit)         
      write(outu,*) 
      write(outu,'(6x,a,i3)') 'Reading static external field from unit ',iunit
      call READPHI(iunit,outu,'PHIX',Qadj)
    endif ! Qphix
    Qrfpar = check(com,'rfpar')
    if (Qrfpar) then
      Qrfpsin=check(com,'rfpsingle')
      call gtcpar(com,'rfparunit',nn,word)
      if (nn.eq.0) call error ('shell_simul','rfparunit cannot be empty or ommited',faterr)
      allocate (iunitv(nn))
      call gtcipar(word,iunitv)
      do i=1,nn
        iunit=iunitv(i)
        if (iunit.le.0) call error ('shell_simul', 'unit is zero or a negative number', faterr)
        if (iunit.gt.maxopen) call error ('shell_simul', 'unit is greater than maxopen', faterr)
        if (unvec(iunit).eq.-1) call error ('shell_simul', 'unit incorrect in GSBP order', faterr)
        iunitv(i) = unvec(iunit) 
      enddo
      call gtdpar(com,'rfparfac',reffac,1.5)
      sqrfac=sqrt(reffac)
      write(outu,*) 
      write(outu,'(6x,a$)') 'Reading Reaction field parameters from units:'
      write(outu,*) iunitv
      write(outu,'(6x,a$)') 'Using RFPAR factor for effective radius parameters of:'
      write(outu,*) reffac
      call readrfpar(iunitv,nn,outu,Qadj)
      deallocate (iunitv)
    endif
    if(Qmmij.and.Qrfpar) call error ('shell_simul','Error: MMIJ and RFPAR are not compatibles',faterr)
    write(outu,*)
  ! **************************************************************************
  elseif (wrd5.eq.'svdw') then
    if (.not.Qpar.and..not.Qnucl) call error ('shell_simul', 'SVDW order is defined before PARTICLE and/or NUCLEOTIDE order', faterr)
    if (.not.Qphiv) call error ('shell_simul', 'Repulsion field has not been read', faterr)
    if (Qnmcden) call error ('shell_simul', 'SVDW not compatible with NMCDEN', faterr)
    allocate (scal(nttyp))
    scal=1.0
    endlog = .false.
    do while (.not. endlog)
      call getlin(com,inpu,outu) ! new commands line
      endlog = check(com,'end')
      if (.not.endlog) then
        ! Obtention of ion type atnam(itype)
        call getfirst(com,wrd4)
        call fatnam(atnam2,nttyp,wrd4,itype)
        call gtdpar(com,'scale',scal(itype),1.0)
      endif
    enddo
    Qsvdw = .true.
    write(outu,'(/6x,a/6x,a)') 'SVDW Type Scaling Factor enabled','Type   Factor'
    do i=1,nttyp
      write(outu,'(6x,a,3x,f10.5)') atnam2(i),scal(i)
    enddo
  elseif (wrd5.eq.'exit') then
  !        ---------------
     write(outu,*)
     finish = timer()-start
     write(outu,'(/6x,a)') 'CPU Time:'
     write(outu,'(8x,a)') 'Format 1:'
     write(com,*) finish,' msec / ',finish/1e3,' sec / ',finish/6e4,' m / ',finish/36e5,' h / ',finish/864e5,' d'
     write(outu,'(8x,a)') trim(com)
     write(outu,'(8x,a)') 'Format 2:'
     values(3)=int(finish/864e5)  ! days
     values(1)=24*values(3) ! hours eq
     values(5)=int(finish/36e5)-values(1) ! hours
     values(1)=60*(values(1)+values(5)) ! min eq
     values(6)=int(finish/6e4)-values(1) ! minutes
     values(1)=60*(values(1)+values(6)) ! sec eq
     values(7)=int(finish/1e3)-values(1) ! seconds
     values(1)=1000*(values(1)+values(7)) ! milisec eq
     values(8)=finish-values(1) ! miliseconds
     write(com,*) values(3),' d, ',values(5),' h, ',values(6),' m, ',values(7),' s, ',values(8),' ms'
     write(outu,'(8x,a)') trim(com)
     call date_and_time(date, time, zone, values)
     write(outu,'(/6x,a,7(i0,a))') 'Finished at (YYYY-MM-DD HH:mm:ss.ms): ',values(1),'-',values(2),'-',values(3),' ',values(5),':',values(6),':',values(7),'.',values(8)
     logfinal = .true.
  106     format(6x,a,f15.2,a)
  else
    write(outu,'(6x,a)') '*ERROR*  Unrecognized command:'
  endif
enddo
return
end
