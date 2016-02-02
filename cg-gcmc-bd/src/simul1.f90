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

subroutine simul1(ncycle, ngcmc, nmcm, nbd, nsave, nsfbs, vfbs, ntras, nsec, iseed)
use ioxmod
use grandmod
use stdiomod
use constamod
use gsbpmod
use nucleotmod
use errormod
use charfuncmod, only: sng,i8toi   !command parser

implicit none
integer*8   ncycle,i,icycle,nsec2,nsave,nsfbs
integer ngcmc, nmcm, nbd, ntras, nsec, j,istep,np1,np2,ii
real    vfbs
integer   iat, jat, itype, ib, n, iz, ir
real    naver
!real*16 dener
real dener
real    current, area, zz, vol1, vol2, pres, vir
integer   ncount(ntype-nold) 
integer   prob(nbuffer,0:datom) 
integer   prob2(ntype,0:datom) 
integer   prob3(ntype,0:datom)      
integer   irmax, nzmax
integer   npoints,iseed,pncross(ntype-nold,cntpts)
parameter (npoints = 10000)
real    dr, idelz,delz, zmini, dist
real    gr(ntype,npoints), rho(ntype,npoints), rDNA(6,npoints), sgr(ntype) 
real    r,flux 
real,parameter :: i3 = 1.0/3.0

!ion pairing S frequency
integer   nion1(ntype,npoints) 
integer   npair1(ntype,npoints),npair2(ntype,npoints) 
integer   npair3(ntype,npoints),npair4(ntype,npoints) 
real    s1,s2,itvol

!energy profiles
integer  nion2(ntype,npoints) 
!real*16   etmp,etot(ntype,npoints),eion(ntype,npoints) 
!real*16   erf(ntype,npoints),esf(ntype,npoints), egvdw(ntype,npoints)
real   etmp,etot(ntype,npoints),eion(ntype,npoints) 
real   erf(ntype,npoints),esf(ntype,npoints), egvdw(ntype,npoints)

!fraction of denatured bases 
real  ibond, ftime
integer i1, i2, ipos, nstf 
logical*1 logvab,ok

!translocation duration
real  zlpore, zgpore
real  radius, fmemb
integer nmemb

! line
character*8192 ln

itvol=1.0/tvol
Qforces = .false.

!OUTPUT
if (Qpar.and.Qnucl) then
  write (outu,'(/6x,a)') 'IONIC & DNA SIMULATIONS'
else if (Qpar.and..not.Qnucl) then
  write (outu,'(/6x,a)') 'IONIC SIMULATIONS'      
else if (Qnucl.and..not.Qpar) then
  write (outu,'(/6x,a)') 'DNA SIMULATIONS'
endif      

if (.not.Qenergy) then
   write(outu,'(6x,a)') 'Energy and forces will be skipped'
   Qforces = .false.
endif

write(outu,131) 'NCYCLE = ',ncycle,'NGCMC  = ',ngcmc
write(outu,131) 'NMCM   = ',nmcm  ,'NBD    = ',nbd
write(outu,131) 'NPRINT = ',nprint,'NSAVE  = ',nsave
write(outu,131) 'SEED   = ',iseed ,'NANAL  = ',nanal
131  format(6x,2(a,i12,5x))

write(outu,*)
write(outu,132) 'TEMPERATURE         = ',temp
write(outu,132) 'DIELECTRIC CONSTANT = ',cdie
132  format(6x,a,f10.3)

if (ncycle*nbd.ne.0) then
  write(outu,133) 'TIME STEP   = ',dt,' [ps]'
  write(outu,133) 'Brownian dynamics simulation of ',(ncycle*nbd)*dt*pico,' [s]'
endif
133    format(6x,a,e10.3,a)

if (Qenergy) then
  write(outu,*)
  call energy
  if (Qpar.and.Qnucl) then
    write(outu,'(10x,a)') 'CYCLE--------Total-PMF--------Nonbonded--------Bonded-----------PHIsf------------PHIrf----------PHIvdW'
    write(outu,'(5x,i10,6f17.4)') 0,ener,eelec+evdw+esrpmf+estack+ebp+eex+eqq+esolv+eqqmx+evdwmx+esrpmfmx+eefpot+ &
     eefpotmx,ebond+eang+edihe,egsbpa,egsbpb,evdwgd   
    if (Qcontrans.and.Qcontprint) then
      write(outu,'(10x,a,f17.10)') 'Constrain Energy=',econ
      do ii=1,ctn
        j=csn(ii)
        if (j.eq.0) then
          xcon=sum(x(1:nsites))*insites-contrx(ii)
          ycon=sum(y(1:nsites))*insites-contry(ii)
          zcon=sum(z(1:nsites))*insites-contrz(ii)
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,'cent',' dx=',xcon,' dy=',ycon,' dz=',zcon
        else
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,namsite(j),' dx=',x(j)-contrx(ii),' dy=',y(j)-contry(ii),' dz=',z(j)-contrz(ii)
        endif
      enddo
    endif 
    write(outu,'(10x,a)') '------------------------------------------------------------------------------------------------------'  
    if (ngcmc.gt.0) then
      ncount(1:ntype-nold) = 0
      do iat = nsites+1, ntot
        itype = abs(typei(iat))-nold
        ncount(itype) = ncount(itype) + 1
      enddo
      write(ln,*) '          ',(atnam(itype),'>',ncount(itype-nold),' | ',itype=nold+1,ntype)
      write(outu,'(a)') trim(ln)
      write(outu,'(10x,a)') '------------------------------------------------------------------------------------------------------'  
    endif         
  else if (Qpar.and..not.Qnucl) then
    write(outu,'(10x,a)') 'CYCLE--------Total-PMF--------Nonbonded------------PHIsf------------PHIrf-----------PHIvdW' 
    write(outu,'(5x,i10,5f17.4)') 0,ener,eelec+evdw+esrpmf+eefpot, egsbpa,egsbpb,evdwgd 
    write(outu,'(10x,a)') '------------------------------------------------------------------------------------------'  
    if (ngcmc.gt.0) then
      ncount(1:ntype-nold) = 0
      do iat = nsites+1, ntot
        itype = abs(typei(iat))-nold
        ncount(itype) = ncount(itype) + 1
      enddo
      write(ln,*) '          ',(atnam(itype),'>',ncount(itype-nold),' | ',itype=nold+1,ntype)
      write(outu,'(a)') trim(ln)
      write(outu,'(10x,a)') '------------------------------------------------------------------------------------------------------'  
    endif         
  else if (Qnucl.and..not.Qpar) then 
    write(outu,'(10x,a)') 'CYCLE----------Total-PMF--------Nonbonded--------Bonded------------PHIsf------------PHIrf------------PHIvdW' 
    write(outu,'(5x,i10,6f17.4)') 0,ener,estack+ebp+eex+eqq+esolv,ebond+eang+edihe,egsbpa,egsbpb,evdwgd 
    if(Qcontrans.and.Qcontprint) then
      write(outu,'(10x,a,f17.10)') 'Constrain Energy=',econ
      do ii=1,ctn
        j=csn(ii)
        if (j.eq.0) then
          xcon=sum(x(1:nsites))*insites-contrx(ii)
          ycon=sum(y(1:nsites))*insites-contry(ii)
          zcon=sum(z(1:nsites))*insites-contrz(ii)
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,'cent',' dx=',xcon,' dy=',ycon,' dz=',zcon
        else
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,namsite(j),' dx=',x(j)-contrx(ii),' dy=',y(j)-contry(ii),' dz=',z(j)-contrz(ii)
        endif
      enddo
    endif 
    write(outu,'(10x,a)') '-----------------------------------------------------------------------------------------------------------'  
  endif
endif ! Qenergy  
   
!INITIALIZE some statistical variables

if (Qrho .or. Qrdna .or. Qionpair .or. Qenerprofile) then
  zmini  = lz2m
  idelz  = 2.0
  delz  = 0.5
  nzmax = int(lz*idelz)
  if (nzmax.gt.npoints) then 
   call error ('simul1', 'nzmax > npoints', warning)
   write(outu,'(3x,a,1x,i10)') 'nzmax',nzmax
   nzmax = npoints
  endif
endif

if (ngcmc.gt.0) then
!i. number probability in the buffer
  do iat = 0, datom
    do ib = 1, nbuffer
      prob(ib,iat)  = 0 
    enddo
  enddo
endif   

if (Qprob) then
!ii. number probability in the system
  do iat = 0, datom
    do itype = nold+1, ntype
      prob3(itype,iat)  = 0   
    enddo
  enddo

!iii. number probability in the channel
  do iat = 0, datom 
    do itype = nold+1, ntype
      prob2(itype,iat)  = 0 
    enddo
  enddo

endif ! Qprob

!iv. number of ion crossing
if (Qpar.and.Qcountion) then
  nforward   = 0
  nbackward  = 0
  pncross = 0
endif

if (Qrho .or. Qgr) then
!v. radial distribution function g(r) and average density profile 
!along the Z-axis rho(z)
  irmax = 0
  dr    = 0.1 ! for g(r): small fixed distance
  do itype = nold+1, ntype 
    sgr(itype) = 0.0
    do iz = 1, npoints
      rho(itype,iz) = 0.0
      gr(itype,iz)  = 0.0
    enddo
  enddo

endif ! Qrho.or.Qgr

!vi. average density profile along the Z-axis for DNA sites rDNA
if (Qrdna) then 
  do iz = 1, npoints
    do itype = 1, nttyp
      rDNA(itype,iz) = 0.0
    enddo
  enddo
endif ! Qrdna

if (Qionpair) then
!vii. ion pairing S frequency
  s1 = 4.0
  s2 = 6.45
  do itype = nold+1, ntype
    do iz = 1, npoints
      nion1(itype,iz)  = 0
      npair1(itype,iz) = 0
      npair2(itype,iz) = 0
      npair3(itype,iz) = 0
      npair4(itype,iz) = 0
    enddo
  enddo
endif ! Qionpair

if (Qenerprofile) then
!viii. energy profiles
  do itype  = nold+1, ntype
    do iz = 1, npoints
      nion2(itype,iz)  = 0
      etot(itype,iz)   = 0.0
      eion(itype,iz)   = 0.0
      erf(itype,iz)    = 0.0
      esf(itype,iz)    = 0.0
      egvdw(itype,iz)  = 0.0
    enddo
  enddo
endif ! Qenerprofile

!ix. trajectory
if (Qtraj.and..not.Qtrajcont) then
  if (setframes.eq.0) then
    nframe = i8toi(ncycle/nsave)
  else
    nframe = setframes
  endif
  write(iuntrj) nframe                 ! number of frames
  write(iuntrj) nttyp                  ! number of ions and nucleotides
  write(iuntrj) (atnam2(ii),ii=1,nttyp)  ! ion and nucleotides types in character*4          
endif

!x. fraction of denatured bases
if (Qfbases) then
  write(iunfbs,'(2x,a)') 'DNA MELTING'
  write(iunfbs,'(2x,a,1x,i5)') 'Total interstrand hydrogen'//' bonding', inuc
  write(iunfbs,'(2x,a,1x,e13.5)') 'Initial time for calculating'//' fraction of denatured', vfbs
  nstf = 0
  do i = 1, ncycle
    ftime = (i*nbd)*dt
    if (mod(i,nsfbs).eq.0 .and. ftime.gt.vfbs) nstf = nstf + 1
  enddo
  write(iunfbs,'(2x,a,1x,i10)') 'Number of steps', nstf
  write(iunfbs,'(2x,a)') 'time [ps]--broken interstrand hydrogen'//' bonding--fraction of denatured bases'
endif

!xi. translocation duration
if (Qfmemb) then
  ok = .not.Qmemb .and. czmin.eq.0.0.and.czmax.eq.0.0
  if (ok) Qfmemb = .false.
endif
if (Qfmemb) then
  write(iuntfm,'(2x,a)') 'time [ps]--fraction of DNA sites inside the pore' 
endif 

!xii. Security ouputfile
if (Qsec) then
  nsec2 = ncycle/nsec
  write(iuntsc,'(2x,a,1x,i15)') 'maximum number of frames',nsec2
endif        

!xii. forces        
do ii = 1, ntot
  fx(ii) = 0.0
  fy(ii) = 0.0
  fz(ii) = 0.0
enddo


!..........................................................................
!     Main loop
!..........................................................................

if (Qpres) pres=0.0

if (Qbuf) then
  do ib = 1, nbuffer
    ntotat(ib) = nint(avnum(ib))
  enddo
endif  

if (Qdiffuse.and.Qprofile) then
  call error('simul1','Cannot use DIFFUSION and PROFILE',faterr)
endif

kbtdt = dt * ikbt
fact0n1=diffnuc*kbtdna*dt
fact0n2=sqrt(2.0*dt*diffnuc)
do ii=1+nold,ntype
  fact1a(ii)=diffusion(ii)*kBTdt
  fact2a(ii)=sqrt(2.0*dt*diffusion(ii))
enddo
if (Qproxdiff) fact2pd=sqrt(2.0*dt*diff0)
 
do icycle = 1, ncycle

  if (ngcmc.gt.0) then
    call grand(ngcmc,prob,icycle)
  endif 

  if (nmcm.gt.0) then
    call metropolis(nmcm)
  endif 

  if (nbd.gt.0) then
    Qforces = .true.
    do istep = 1, nbd
      call energy
      if (Qpres) then
        vir=dot_product(fx(1:ntot),x(1:ntot))+dot_product(fy(1:ntot),y(1:ntot))+dot_product(fz(1:ntot),z(1:ntot))
        pres=pres+(ntot-vir*i3*ikbt)*itvol*kba3bar*temp
      endif
      if (Qpar) then
        if (Qdiffuse) then
          call dynamics1(nsites+nfix+1,ntot) ! using non-uniform diffusion constant
        elseif (Qprofile) then
          call dynamics2(nsites+nfix+1,ntot) ! using diffusion constant profile 
        elseif (Qproxdiff) then
          call dynamics3(nsites+nfix+1,ntot) ! using dna proximity diffusion constant
        else
          call dynamics0(nsites+nfix+1,ntot)
        endif
      endif
      if (Qnucl) then
        call dynamics0nuc(1,nsites)
        if (Qnotrans) then
          if (Qnotrx) x(1:nsites)=x(1:nsites)-sum(x(1:nsites))*insites+notrx
          if (Qnotry) y(1:nsites)=y(1:nsites)-sum(y(1:nsites))*insites+notry
          if (Qnotrz) z(1:nsites)=z(1:nsites)-sum(z(1:nsites))*insites+notrz
        endif
      endif ! Qnucl
    enddo
    Qforces = .false.
  endif  ! nbd   

  ! Calculate charge density
  if (Qchdencnt) call chrden 

 !write a trajectory if required
  if (Qtraj .and. mod(icycle,nsave).eq.0) then
    if (nbd.gt.0) then 
      runtime = (icycle*nbd)*dt !*pico*1.0E9
    else
      runtime = icycle
    endif
    call wrttraj
  endif

  !write a fraction of denatured bases if required
  if (Qfbases) then
    ftime = (icycle*nbd)*dt
    if (mod(icycle,nsfbs).eq.0 .and. ftime.gt.vfbs) then
      ibond = 0.0E0
      do ii = 1, nbp
        i1 = sitebp(ii,1)
        i2 = sitebp(ii,2)
        if (strand(i1).ne.strand(i2)) then
          if (strand(i1).eq.1) then
            ipos = inuc - typenuc(i1) + 1
            logvab = typenuc(i2).eq.(ipos+inuc)
          else
            ipos = inuc - typenuc(i2) + 1
            logvab = typenuc(i1).eq.(ipos+inuc)
          endif
          if (logvab) then
            dist = sqrt((x(i1)-x(i2))**2 + (y(i1)-y(i2))**2 + (z(i1)-z(i2))**2)
            if (dist.ge.(sgbp(i)+2.0)) ibond = ibond + 1.0
          endif
        endif
      enddo   
      fbases = ibond/inuc
      write(iunfbs,'(2x,e17.8,1x,i5,1x,f7.4)') ftime,nint(ibond),fbases
    endif
  endif 

  !write traslocation duration if required
  if (Qfmemb) then
    if (mod(icycle,ntras).eq.0) then
      nmemb = 0     
      if (Qmemb) then
        zlpore = zmemb1
        zgpore = zmemb2
      else
        zlpore = czmin
        zgpore = czmax
      endif 
      do ii = 1, nsites 
        if (z(ii).ge.zlpore .and. z(ii).le.zgpore) then
          nmemb = nmemb + 1
          if (Qmemb) then
            radius = sqrt(x(ii)**2+y(ii)**2)
            itype = typtyp(ii)
            if (radius.gt.rcylinder(itype)) then
              call error ('simul1', 'DNA sites are insides of membrane', warning)
              write(outu,*) '**site--name--x,y,z',ii,namnucl(ii),x(ii),y(ii),z(ii)
            endif 
          endif
        endif 
      enddo  
      fmemb = 1.0*nmemb/nsites
      ftime = (icycle*nbd)*dt
      write(iuntfm,'(2x,e17.8,1x,f7.4)') ftime, fmemb 
    endif
  endif

  call count

  !Security outputfile contains coordinates and seed nnumbers
  if (Qsec) then
    if (mod(icycle,nsec).eq.0) then
      write(iuntsc,'(2x,a,1x,i15)') 'cycle',icycle
      write(iuntsc,'(2x,a,1x,i15,1x,a,1x,i10)') 'ntot',ntot,'seed number',iseed
      write(iuntsc,'(2x,a)') 'Coordinates:atom--x--y--z'
      do ii = 1, ntot
        write(iuntsc,'(2x,i15,1x,f10.5,1x,f10.5,1x,f10.5)') ii,x(ii), y(ii), z(ii)
      enddo
    endif 
  endif

  !average density profile along the Z-axis                  
  if (Qrho) then
    do iat = nsites+nfix+1, ntot
      itype = abs(typei(iat))
      iz = int((z(iat)-zmini)*idelz) + 1
      if (iz.le.nzmax) then
        rho(itype,iz) = rho(itype,iz) + 1.0
      endif
    enddo
  endif ! Qrho
  !average density profile along the Z-axis for DNA sites
  if (Qrdna) then
    do iat = 1, nsites
      itype = typtyp(iat)
      iz = int((z(iat)-zmini)*idelz) + 1
      if (iz.le.nzmax) then
        rDNA(itype,iz) = rDNA(itype,iz) + 1.0
      endif
    enddo
  endif ! Qrdna
  !radial distribution function 
  if (Qgr) then
    do iat = nsites+nfix+1, ntot
      itype = abs(typei(iat))
      dist = (x(iat)-x(igr))**2+(y(iat)-y(igr))**2+(z(iat)-z(igr))**2
      dist = sqrt(dist) ! distance
      ir = int(dist/dr) + 1 ! spherical shell
      if (ir.gt.irmax) irmax = ir
      if (irmax.gt.npoints) irmax = npoints
      if (ir.le.npoints) then
        gr(itype,ir) = gr(itype,ir) + 1.0
      endif
    enddo 
  endif ! Qgr
  !average probability distributions
  if (Qprob) then
  !  channel
    ok = .not.Qmemb .and. czmin.eq.0.0.and.czmax.eq.0.0
    if (.not.ok) then
      ncount(1:ntype-nold) = 0 
      do iat = nsites+nfix+1, ntot
        itype = abs(typei(iat))-nold
        if (Qmemb) then ! cylindrical channel
          r = sqrt(x(iat)**2+y(iat)**2)
          ok = z(iat).ge.zmemb1.and.z(iat).le.zmemb2
          if (ok) ok = rcylinder(nwtype(itype+nold)).gt.0.0.and.r.le.rcylinder(nwtype(itype+nold))
        else ! otherwise
          ok = z(iat).ge.czmin.and.z(iat).le.czmax
        endif
        if (ok) ncount(itype) = ncount(itype) + 1
      enddo
      do itype = nold+1, ntype
        prob2(itype,ncount(itype-nold)) = prob2(itype,ncount(itype-nold)) + 1
      enddo
    endif
    !system
    ncount(1:ntype-nold) = 0 
    do iat = nsites+nfix+1, ntot
      itype = abs(typei(iat))-nold
      ncount(itype) = ncount(itype) + 1
    enddo
    do itype = nold+1, ntype
      prob3(itype,ncount(itype-nold)) = prob3(itype,ncount(itype-nold)) + 1
    enddo
  endif ! Qprob
  !ion pairing analysis (S frequency)
  if (Qionpair) then         
    if (mod(icycle,nanal).eq.0) then
      do itype = nold+1, ntype
        do iat = nsites+nfix+1, ntot
          if (abs(typei(iat)).eq.itype) then 
            iz = int((z(iat)-zmini)*idelz) + 1
            if (iz.le.nzmax) then
              nion1(itype,iz) = nion1(itype,iz) + 1
              np1 = 0
              np2 = 0
              do jat = nsites+nfix+1, ntot
                if (abs(typei(jat)).ne.itype .and.abs(x(iat)-x(jat)).le.s2 .and.abs(y(iat)-y(jat)).le.s2 .and.abs(z(iat)-z(jat)).le.s2) then        
                  r = sqrt((x(iat)-x(jat))**2 + (y(iat)-y(jat))**2 + (z(iat)-z(jat))**2)
                  if (r.le.s1) then
                    np1 = np1 + 1
                  else if (r.gt.s1.and.r.le.s2) then
                    np2 = np2 + 1
                  endif
                endif   
              enddo ! jat
              if (np1.eq.0.and.np2.eq.0) then
                npair1(itype,iz) = npair1(itype,iz) + 1
              else if (np1.eq.1.and.np2.eq.0) then
                npair2(itype,iz) = npair2(itype,iz) + 1
              else if (np1.eq.0.and.np2.eq.1) then
                npair3(itype,iz) = npair3(itype,iz) + 1
              else if (np1.eq.1.and.np2.eq.1) then
                npair4(itype,iz) = npair4(itype,iz) + 1
              endif
            endif  
          endif
        enddo ! iat
      enddo ! itype
    endif 
  endif ! Qionpair
  !enegy profiles along Z-axis
  if (Qenerprofile) then
    if (mod(icycle,nanal).eq.0) then      
      do itype = nold+1, ntype
        do iat= nsites+nfix+1, ntot
          if (abs(typei(iat)).eq.itype) then
            iz = int((z(iat)-zmini)*idelz) + 1
            if (iz.le.nzmax) then
              nion2(itype,iz) = nion2(itype,iz) + 1
              call interact(dener,x(iat),y(iat),z(iat),itype,iat,.true.)
              etmp = dener - egsbpa - egsbpb - evdwgd 
              etot(itype,iz) = etot(itype,iz) + dener    ! total interaction energy
              eion(itype,iz) = eion(itype,iz) + etmp     ! ion-ion interaction energy
              esf(itype,iz) = esf(itype,iz) + egsbpa     ! static field energy
              erf(itype,iz) = erf(itype,iz) + egsbpb     ! reaction field energy
              egvdw(itype,iz) = egvdw(itype,iz) + evdwgd ! repulsive energy
            endif   
          endif
        enddo ! iat 
      enddo ! itype
    endif 
  endif  ! Qenerprofile

 !write counting ions if applies
  if (Qpar.and.Qcountion.and.nbd.gt.0) then
    if (mod(icycle*nbd,svcntfq).eq.0) call ioncountout(icycle*nbd,pncross)
  endif

  !OUTPUT
  if (nprint.gt.0 .and. mod(icycle,nprint).eq.0) then  
    if (ncycle*nbd.eq.0) call energy
    if (Qpar .and. Qnucl) then
      write(outu,'(10x,a)') 'CYCLE----------Total-PMF--------Nonbonded---------Bonded-----------PHIsf------------PHIrf----------PHIvdW'
      write(outu,'(5x,i12,6f17.4)') icycle,ener,eelec+evdw+esrpmf+estack+ebp+eex+eqq+esolv+eqqmx+evdwmx+esrpmfmx+eefpot+ &
           eefpotmx,ebond+eang+edihe,egsbpa,egsbpb,evdwgd
      if(Qcontrans.and.Qcontprint) then
        write(outu,'(10x,a,f17.10)') 'Constrain Energy=',econ
        do ii=1,ctn
          j=csn(ii)
          if (j.eq.0) then
            xcon=sum(x(1:nsites))*insites-contrx(ii)
            ycon=sum(y(1:nsites))*insites-contry(ii)
            zcon=sum(z(1:nsites))*insites-contrz(ii)
            write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,'cent',' dx=',xcon,' dy=',ycon,' dz=',zcon
          else
            write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,namsite(j),' dx=',x(j)-contrx(ii),' dy=',y(j)-contry(ii),' dz=',z(j)-contrz(ii)
          endif
        enddo
      endif 
      write(outu,'(10x,a)') '-------------------------------------------------------------------------------------------------------'
    else if (Qpar.and..not.Qnucl) then        
      write(outu,'(10x,a)') 'CYCLE----------Total-PMF--------Nonbonded-----------PHIsf------------PHIrf----------PHIvdW'
      write(outu,'(5x,i12,5f17.4)') icycle,ener,eelec+evdw+esrpmf+eefpot,egsbpa,egsbpb,evdwgd
      write(outu,'(10x,a)') '------------------------------------------------------------------------------------------'
    else if (Qnucl.and..not.Qpar) then
      write(outu,'(5x,i12,6f17.4)') icycle,ener,estack+ebp+eex+eqq+esolv,ebond+eang+edihe,egsbpa,egsbpb,evdwgd
      if(Qcontrans.and.Qcontprint) then
        write(outu,'(10x,a,f17.10)') 'Constrain Energy=',econ
        do ii=1,ctn
          j=csn(ii)
          if (j.eq.0) then
            xcon=sum(x(1:nsites))*insites-contrx(ii)
            ycon=sum(y(1:nsites))*insites-contry(ii)
            zcon=sum(z(1:nsites))*insites-contrz(ii)
            write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,'cent',' dx=',xcon,' dy=',ycon,' dz=',zcon
          else
            write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,namsite(j),' dx=',x(j)-contrx(ii),' dy=',y(j)-contry(ii),' dz=',z(j)-contrz(ii)
          endif
        enddo
      endif 
    endif        
    if (Qpar.and.ngcmc.gt.0) then
      ncount(1:ntype-nold) = 0
      do iat = nsites+1, ntot
        itype = abs(typei(iat))-nold
        ncount(itype) = ncount(itype) + 1
      enddo
      write(ln,*) '          ',(atnam(itype),'>',ncount(itype-nold),' | ',itype=nold+1,ntype)
      write(outu,'(a)') trim(ln)
      if (Qpres.and.nbd.gt.0) write(outu,'(6x,a,f18.5,a)') 'Pressure: ',pres/(nbd*icycle),' bar'
      write(outu,'(10x,a)') '------------------------------------------------------------------------------------------------------'  
    endif
  endif

enddo ! icycle
!..........................................................................

if (Qchdencnt) chden=sng(1.0/float(ncycle))*chden

call energy
if (nprint.eq.0) then
  if (Qpar .and. Qnucl) then
    write(outu,'(10x,a)') 'CYCLE----------Total-PMF--------Nonbonded--------Bonded-----------PHIsf------------PHIrf----------PHIvdW'
    write(outu,'(5x,i12,6f17.4)') icycle,ener,eelec+evdw+esrpmf+estack+ebp+eex+eqq+esolv+eqqmx+evdwmx+esrpmfmx+eefpot+eefpotmx &
    ,ebond+eang+edihe,egsbpa,egsbpb,evdwgd
    if(Qcontrans.and.Qcontprint) then
      write(outu,'(10x,a,f17.10)') 'Constrain Energy=',econ
      do ii=1,ctn
        j=csn(ii)
        if (j.eq.0) then
          xcon=sum(x(1:nsites))*insites-contrx(ii)
          ycon=sum(y(1:nsites))*insites-contry(ii)
          zcon=sum(z(1:nsites))*insites-contrz(ii)
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,'cent',' dx=',xcon,' dy=',ycon,' dz=',zcon
        else
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,namsite(j),' dx=',x(j)-contrx(ii),' dy=',y(j)-contry(ii),' dz=',z(j)-contrz(ii)
        endif
      enddo
    endif 
    write(outu,'(10x,a)') '--------------------------------------------------------------------------------------------------------'
  else if (Qpar.and..not.Qnucl) then
    write(outu,'(10x,a)') 'CYCLE----------Total-PMF--------Nonbonded-----------PHIsf------------PHIrf----------PHIvdW'
    write(outu,'(5x,i12,5f17.4)') icycle,ener,eelec+evdw+esrpmf+eefpot,egsbpa,egsbpb,evdwgd
    write(outu,'(10x,a)') '----------------------------------------------------------------------------------------'
  else if (Qnucl.and..not.Qpar) then 
    write(outu,'(5x,i12,6f17.4)') icycle,ener,estack+ebp+eex+eqq+esolv,ebond+eang+edihe,egsbpa,egsbpb,evdwgd
    if(Qcontrans.and.Qcontprint) then
      write(outu,'(10x,a,f17.10)') 'Constrain Energy=',econ
      do ii=1,ctn
        j=csn(ii)
        if (j.eq.0) then
          xcon=sum(x(1:nsites))*insites-contrx(ii)
          ycon=sum(y(1:nsites))*insites-contry(ii)
          zcon=sum(z(1:nsites))*insites-contrz(ii)
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,'cent',' dx=',xcon,' dy=',ycon,' dz=',zcon
        else
          write(outu,'(10x,2i4,x,a5,x,3(a,f10.5))') ii,j,namsite(j),' dx=',x(j)-contrx(ii),' dy=',y(j)-contry(ii),' dz=',z(j)-contrz(ii)
        endif
      enddo
    endif 
    write(outu,'(10x,a)') '----------------------------------------------------------------------------------------'
  endif   
  if (Qpar.and.ngcmc.gt.0) then
    ncount(1:ntype-nold) = 0
    do iat = nsites+1, ntot
      itype = abs(typei(iat))-nold
      ncount(itype) = ncount(itype) + 1
    enddo
    write(ln,*) '          ',(atnam(itype),'>',ncount(itype-nold),' | ',itype=nold+1,ntype)
    write(outu,'(a)') trim(ln)
    if (Qpres.and.nbd.gt.0) write(outu,'(6x,a,f18.5,a)') 'Pressure: ',pres/(nbd*icycle),' bar'
      write(outu,'(10x,a)') '------------------------------------------------------------------------------------------------------'  
  endif  
endif

if (Qrdna) then
  write(outu,*)
  write(outu,*) 'DNA average density profile: '
  write(outu,*) 'ntype = ',nttyp,' delz = ',delz,' zmin = ',zmini
  write(outu,*) 'DNA sites types: ',(atnam2(ii),ii=1,nttyp)
  do iz = 1, nzmax
    zz = zmini + (iz*1.0-0.5)*delz
    area = lx*ly
    if (Qsphere) area = pi*(Rmax(1)**2-zz**2)
    do itype = 1, nttyp
      rDNA(itype,iz) = rDNA(itype,iz)/(ncycle*delz*area)
    enddo
    write(outu,'(2f10.3,2x,8(e11.4,2x))') zz, area, (rDNA(ii,iz),ii=1,nttyp)
  enddo 
endif 

if (Qpar) then
  write (outu,*)
  write (outu,'(6x,a)') 'Results from the simulation: '
  write (outu,'(6x,a)') '-----------------------------'

  write (outu,'(6x,a,i4)') 'Total number of ions ',ntot-nsites
  if (Qbuf) then
    call count
    do ib = 1, nbuffer
      write (outu,'(6x,a,2i4)') 'buffer--number of ions ',ib,nat(ib)
    enddo
  endif 

  write (outu,*) 
  write (outu,'(6x,a)') 'Statistics: '
  write (outu,*) 
  ncount = 0
  do ii = nsites+nfix+1, ntot
    itype = abs(typei(ii))-nold
    ncount(itype) = ncount(itype) + 1
  enddo
  write(outu,'(6x,a)') 'Total number of ions per type:' 
  do itype=1,ntype-nold
    write(outu,'(6x,a,a,i0)') atnam(itype+nold),' =  ',ncount(itype)
  enddo
  write (outu,*) 
  
  if (ngcmc.gt.0) then
    if (Qpar.and.Qcountion.and.nbd.gt.0) then
      do ii=1,cntpts 
        current = 0.0
        write (outu,'(6x,a,i0,a,f11.4)') 'Counting Zone= ',ii,'  Counting z-pos= ',zcont(ii)
        do itype = 1, ntype-nold
          write (outu,'(8x,a,a)') 'Particle type ',atnam(itype+nold)
          write (outu,'(8x,2(a,i7))') 'forward   ',nforward(itype,ii),'   backward ',nbackward(itype,ii)
          write (outu,'(8x,a,i7)') 'Net Number of particle going forward:',nforward(itype,ii)-nbackward(itype,ii)
          write(outu,*) 
          flux = (nforward(itype,ii)-nbackward(itype,ii))/(dt*float(nbd*ncycle))
          write (outu,'(8x,a,e11.4)') 'Flux in [particle/ps]   ',flux
          write (outu,'(8x,a,e11.4)') 'Current [in pA]         ',cg(itype+nold)*(flux*coulomb*ipico**2)
          current=current+cg(itype+nold)*(flux*coulomb*ipico**2)
          write (outu,*) 
        enddo ! itype
        write (outu,'(8x,a,e11.4)') 'Total Current [in pA]               ',current
        write (outu,*) 
      enddo
    endif
     
    if (Qbuf) then
      do ib = 1, nbuffer
        write (outu,'(6x,a,i3,4x,2(a,i12))') 'buffer ',ib,' nremove   ',nremove(ib),' ninsert   ',ninsert(ib)
        write (outu,'(6x,a)') 'Distribution of particles in this'//' buffer: '
        naver = 0.0
        do n = 0, datom
          if (prob(ib,n).ne.0) then
            naver = naver + n*prob(ib,n)      
            write (outu,'(6x,i4,f12.3)') n,prob(ib,n)*1.0/(ncycle*ngcmc)
          endif
        enddo ! n
        write (outu,'(6x,a,2f12.3)') 'Average number: ',naver/(ncycle*ngcmc),avnum(ib)
        write (outu,*)
      enddo ! ib
    endif

  endif

  if (Qprob) then       
    ok = .not.Qmemb .and. czmin.eq.0.0.and.czmax.eq.0.0
    if (.not.ok) then   
      write (outu,'(6x,a)') 'Distribution of particles in the'//' channel: '
      if (Qmemb) then ! cylindrical channel
        write (outu,'(6x,a,f7.2,a,f7.2,a)') 'The channel is defined between ',zmemb1,' and ',zmemb2,' along Z-axis'
      else ! otherwise
        write (outu,'(6x,a,f7.2,a,f7.2,a)') 'The channel is defined between ',czmin,' and ',czmax,' along Z-axis' 
      endif
      do itype = nold+1, ntype
        write (outu,'(6x,a,1x,i5)') 'type',itype
        ncount(itype-nold) = 0
        do n = 0, datom
          if (prob2(itype,n).ne.0) then
            ncount(itype-nold) = ncount(itype-nold) + n*prob2(itype,n)      
            write (outu,'(6x,i4,f12.3)') n,prob2(itype,n)*1.0/ncycle
          endif
        enddo
        write (outu,'(6x,a,1x,f12.3)') 'Average number in the'//' channel:',ncount(itype-nold)*1.0/ncycle
        write (outu,*)
      enddo ! itype
    endif

    write (outu,'(6x,a)') 'Distribution of particles in'//' the system: '
    do itype = nold+1, ntype
      write (outu,'(6x,a,1x,i5)') 'type ',itype
      ncount(itype-nold) = 0
      do n = 0, datom
        if (prob3(itype,n).ne.0) then
          ncount(itype-nold) = ncount(itype-nold) + n*prob3(itype,n)
          write (outu,'(6x,i4,2f12.3)') n,prob3(itype,n)*1.0/ncycle
        endif
      enddo
      write(outu,'(6x,a,1x,f12.3)') 'Average number in the'//'  system:',ncount(itype-nold)*1.0/ncycle
      write(outu,*)
    enddo ! itype
  endif ! Qprob

  if (Qrho) then
    write (outu,*)
    write (outu,*) 'Average density profile: '
    write (outu,*) 'ntype = ',ntype-nold,' delz = ',delz,' zmin = ',zmini
    do iz = 1, nzmax
      zz = zmini+(iz*1.0-0.5)*delz
      do itype = nold+1, ntype
        if (Qsphere) then
          logvab = .false.
          ib = 0
          do while (ib.lt.nbuffer .and. .not.logvab)
            ib = ib + 1
            if (itype.eq.ibfftyp(ib)) logvab = .true.
          enddo
          if (.not.logvab) then
            call error ('simul1', 'Rmax is not correct for calculating average density profile along the Z-axis',faterr)
          endif  
          area = pi*(Rmax(ib)**2-zz**2)
        else
          area = LX*LY
        endif
        rho(itype,iz)=rho(itype,iz)/(ncycle*delz*area)
      enddo
      write (outu,'(2f10.3,2x,8(e11.4,2x))') zz,area,(rho(ii,iz),ii=nold+1,ntype)
    enddo
  endif ! Qrho

  if (Qgr) then
    write (outu,*)
    write (outu,*) 'Integration number : '
    do ir = 1, irmax
      do itype = nold+1, ntype
         sgr(itype)=sgr(itype)+gr(itype,ir)/ncycle
      enddo
      write (outu,'(f12.3,8e17.4)') (ir*1.0-0.5)*dr,(sgr(ii),ii=nold+1,ntype)
    enddo
    write(outu,*)
    write(outu,*) 'Radial distribution function : '
    do itype= nold+1, ntype
      ! obtention index for density in the buffer
      ib = 1
      ok = .true.
      do while (ib.le.nbuffer .and. ok)
        if (ibfftyp(ib).eq.itype) then
          ii = ib
          ok = .false.
        endif
        ib = ib + 1
      enddo
      if (ok) call error ('simul1', 'problem to calculate the radial distribution function', faterr)
      do ir = 1, irmax
        ! average number of atoms in 'ir' shell
        gr(itype,ir) = gr(itype,ir)/ncycle 
        ! (vol2 - vol1) -> volume of 'ir' shell
        vol1 = 4.0*pi/3.0*((ir-1)*dr)**3
        vol2 = 4.0*pi/3.0*(ir*dr)**3
        gr(itype,ir)=gr(itype,ir)/((vol2-vol1)*density(ii)) 
      enddo
    enddo
    do ir = 1, irmax
      write (outu,'(f12.3,8e17.4)') (ir*1.0-0.5)*dr,(gr(ii,ir),ii=nold+1,ntype)
    enddo
  endif ! Qgr

  if (Qionpair) then
    write (outu,*)
    write (outu,'(6x,a)') 'Ion Pairing Analysis: '
    do itype = nold+1, ntype
      write (outu,*) atnam(itype)
      do iz = 1, nzmax
        zz = zmini+(iz*1.0-0.5)*delz
        if (nion1(itype,iz).gt.0) then
          write (outu,'(6x,f12.3,5f10.5)') zz, nion1(itype,iz)*1.0/ncycle*nanal, npair1(itype,iz)*1.0/nion1(itype,iz), &
             npair2(itype,iz)*1.0/nion1(itype,iz),npair3(itype,iz)*1.0/nion1(itype,iz),npair4(itype,iz)*1.0/nion1(itype,iz)
        endif
      enddo
      write (outu,*)
    enddo
  endif ! Qionpair

  if (Qenerprofile) then
    write (outu,*)
    write (outu,'(6x,a)') 'Energy profile along the Z-axis: '
    write (outu,'(6x,a)') 'zz--probability--etot--eion--esf--erf'//'--egvdw'
    do itype = nold+1, ntype
      write (outu,'(6x,a,1x,a)') 'atom',atnam(itype)
      do iz = 1, nzmax
        zz = zmini+(iz*1.0-0.5)*delz
        if (nion2(itype,iz).gt.0) then
          etot(itype,iz) = etot(itype,iz)/nion2(itype,iz)
          eion(itype,iz) = eion(itype,iz)/nion2(itype,iz)
          esf(itype,iz)  = esf(itype,iz)/nion2(itype,iz)
          erf(itype,iz)  = erf(itype,iz)/nion2(itype,iz)
          egvdw(itype,iz) = egvdw(itype,iz)/nion2(itype,iz)
          write(outu,'(6x,f12.3,6f10.5)') zz,nion2(itype,iz)*1.0/ncycle*nanal,etot(itype,iz),eion(itype,iz),esf(itype,iz),erf(itype,iz),egvdw(itype,iz)
        endif
      enddo
      write(outu,*)
    enddo
  endif ! Qenerprofile

endif ! Qpar

return
end subroutine

