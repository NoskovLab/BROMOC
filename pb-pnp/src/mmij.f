!    PB-PNP - Poisson-Boltzmann and Poisson-Nernst-Planck Equations Solver
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

      subroutine rxnfld(com)
c------------------------------------------------------------------------
c     calculate the reaction field matrix mmij using legendre polynomials
c
      implicit none
      include 'pnp.fcm'
      include 'mmij.fcm'
      include 'mainio.fcm'

      character*(*) com

c local
      integer*4 lncel
      real*8  ldcel,lxbcen,lybcen,lzbcen,ltran 
      real*8  xscale,yscale,zscale
      real*8  rrxcen,rrycen,rrzcen


c     construct a large box to obtain boundary potentials of the original box.
c     here, we just consider a cubic large box.
      call gtipar(com,'lncel',lncel,33) 
      call gtdpar(com,'ldcel',ldcel,dcel*4.0d0) 
      call gtdpar(com,'lxbcen',lxbcen,0.0d0) 
      call gtdpar(com,'lybcen',lybcen,0.0d0) 
      call gtdpar(com,'lzbcen',lzbcen,0.0d0) 
      if(mod(lncel,2).eq.0) lncel=lncel+1
      ltran=0.5d0*(lncel-1)*ldcel

c use the focussing method to get the boundary potential by default
      qphifocus=.true.

c cutoff value for eigenvalues of overlap integral sij matrix
      call gtdpar(com,'tolsij',tolsij,0.003d0) 

c basis charge scaling factor for the monopole basis function
      call gtdpar(com,'cgscal',cgscal,1.0d0)

c number of legendre polynomials
      call gtipar(com,'xnpol',xnpol,0) 
      call gtipar(com,'ynpol',ynpol,0) 
      call gtipar(com,'znpol',znpol,0) 
      ntpol = xnpol*ynpol*znpol

c orignal orthorhombic region where legendre polynomials are applied
      call gtdpar(com,'xmax',rbxmax,0.0d0) 
      call gtdpar(com,'xmin',rbxmin,0.0d0) 
      call gtdpar(com,'ymax',rbymax,0.0d0) 
      call gtdpar(com,'ymin',rbymin,0.0d0) 
      call gtdpar(com,'zmax',rbzmax,0.0d0) 
      call gtdpar(com,'zmin',rbzmin,0.0d0) 

c adjust the original box size to fit into a grid 
c but this step does not appear to be necessary... (changed by sergey and benoit)
      rbxmin = (int((rbxmin+tranx)/dcel)+0.5d0)*dcel-tranx
      rbxmax = (int((rbxmax+tranx)/dcel)+0.5d0)*dcel-tranx
      rbymin = (int((rbymin+trany)/dcel)+0.5d0)*dcel-trany
      rbymax = (int((rbymax+trany)/dcel)+0.5d0)*dcel-trany
      rbzmin = (int((rbzmin+tranz)/dcel)+0.5d0)*dcel-tranz
      rbzmax = (int((rbzmax+tranz)/dcel)+0.5d0)*dcel-tranz
      rrxcen = (rbxmax+rbxmin)/2.0d0
      rrycen = (rbymax+rbymin)/2.0d0
      rrzcen = (rbzmax+rbzmin)/2.0d0

c if a simulation system is a sphere
      call gtdpar(com,'srcut',srcut,0.0d0) 
      if(srcut.gt.0.0d0) then
         call gtdpar(com,'srxcen',rrxcen,0.0d0) 
         call gtdpar(com,'srycen',rrycen,0.0d0) 
         call gtdpar(com,'srzcen',rrzcen,0.0d0) 
         rbxmin=-srcut+rrxcen
         rbxmax=srcut+rrxcen
         rbymin=-srcut+rrycen
         rbymax=srcut+rrycen
         rbzmin=-srcut+rrzcen
         rbzmax=srcut+rrzcen
      rbxmin = (int((rbxmin+tranx)/dcel)+0.5d0)*dcel-tranx
      rbxmax = (int((rbxmax+tranx)/dcel)+0.5d0)*dcel-tranx
      rbymin = (int((rbymin+trany)/dcel)+0.5d0)*dcel-trany
      rbymax = (int((rbymax+trany)/dcel)+0.5d0)*dcel-trany
      rbzmin = (int((rbzmin+tranz)/dcel)+0.5d0)*dcel-tranz
      rbzmax = (int((rbzmax+tranz)/dcel)+0.5d0)*dcel-tranz
      rrxcen = (rbxmax+rbxmin)/2.0d0
      rrycen = (rbymax+rbymin)/2.0d0
      rrzcen = (rbzmax+rbzmin)/2.0d0
      endif

      xscale=2.0d0/(rbxmax-rbxmin)
      yscale=2.0d0/(rbymax-rbymin)
      zscale=2.0d0/(rbzmax-rbzmin)

c information from inputs

      if(srcut.gt.0.0d0) then
         write(outu,102)
         write(outu,'(6x,a,/,6x,a,f8.3,a,/,6x,a,3f8.3)')
     $     'Spherical inner region : ',
     $     ' radius                   (SRDIST) =',srcut,' [Angs]',
     $     ' center           (XCEN,YCEN,ZCEN) =',rrxcen,rrycen,rrzcen
         write(outu,102)
     $        'Inner region in X from ',rbxmin,' to ',rbxmax
         write(outu,102)
     $        'Inner region in Y from ',rbymin,' to ',rbymax
         write(outu,102)
     $        'Inner region in Z from ',rbzmin,' to ',rbzmax
      else
         write(outu,102)
         write(outu,102)
     $        'Inner region in X from ',rbxmin,' to ',rbxmax
         write(outu,102)
     $        'Inner region in Y from ',rbymin,' to ',rbymax
         write(outu,102)
     $        'Inner region in Z from ',rbzmin,' to ',rbzmax
      endif

      write(outu,'(a)') 
      write(outu,'(6x,a,i5,a)') 
     $     'MIJ matrix will be built using',NTPOL,
     $     ' basis functions from'
      write(outu,101) 
     $     'Number of Legendre Polynomials in X  (XNPOL) = ',xnpol
      write(outu,101) 
     $     'Number of Legendre Polynomials in Y  (YNPOL) = ',ynpol
      write(outu,101)
     $     'Number of Legendre Polynomials in Z  (ZNPOL) = ',znpol

      if(cgscal.gt.1.0d0) then
         write(outu,'(/,6x,a,f7.3,a)') 'Charge scaling factor',
     $        CGSCAL,' will be used for the monopole basis function'
      endif

      if(QphiFOCUS) then
         if(tranx.ge.ltran .or.trany.ge.ltran.or.tranz.ge.ltran)then
            write(outu,'(/,6x,a)')
     $        'GSBP WARNING: Large system should contain original one.'
            stop
         endif
         write(outu,'(/,6x,a,/,6x,a)')
     $      'Following large box will be used to calculate boundary' ,
     $      'potentials of the original box using focussing;'
      else
         write(outu,'(/,6x,a,/,6x,a)')
     $      'Charge density on grid of the following large box will be',
     $      'used to calculate boundary potentials of the original box;'
      endif
      write(outu,102)
     $     'Large Box in X from ',lxbcen-ltran,' to ',lxbcen+ltran
      write(outu,102)
     $     'Large Box in Y from ',lybcen-ltran,' to ',lybcen+ltran
      write(outu,102)
     $     'Large Box in Z from ',lzbcen-ltran,' to ',lzbcen+ltran

 101  format(6x,a,i6,a)
 102  format(6x,a,f8.3,a,f8.3)

c .. now, we are ready to calculate the reaction field matrix mmij

      call rfrect0(lncel,ltran,ldcel,lxbcen,lybcen,lzbcen,
     $     xscale,yscale,zscale,rrxcen,rrycen,rrzcen)

c
c      call rfrect1(natom,x,y,z,cg,ntpol,xnpol,ynpol,znpol,
c     $     xscale,yscale,zscale,mij,coef,
c     $     lstpx,lstpy,lstpz,bnorm)
c
      return
      end

      subroutine rfrect1(natom,x,y,z,cg,ntpol,xnpol,ynpol,znpol,
     $           xscal,yscal,zscal,mij,coef,
     $           lstpx,lstpy,lstpz,bnorm)
c-----------------------------------------------------------------------
c     calculate the reaction field energy and forces on each ions
c
      implicit none
      include 'mainio.fcm'
      include 'consta.fcm'

      real*8  x(*),y(*),z(*),cg(*),coef(*),mij(*),bnorm(*)
      real*8  xscal,yscal,zscal
      integer*4 natom,ntpol,xnpol,ynpol,znpol
      integer*4 lstpx(*),lstpy(*),lstpz(*)
c local
      real*8  xmin,xmax,ymin,ymax,zmin,zmax
      real*8  xg,yg,zg,norm,xs,ys,zs,egsbpb
      real*8  lpolx(15),lpoly(15),lpolz(15)
      real*8  mq(1000)
      integer*4 i,ii,jj,ij,n
      integer*4 xpol,ypol,zpol
c
      xmin=-1.0d0/xscal
      xmax= 1.0d0/xscal
      ymin=-1.0d0/yscal
      ymax= 1.0d0/yscal
      zmin=-1.0d0/zscal
      zmax= 1.0d0/zscal

c calculate q_{lm} coefficients
      do ii=1,ntpol
         coef(ii)=0.d0
      enddo
      do 101 i=1,natom
         xg=x(i)   !-rrxcen
         yg=y(i)   !-rrycen
         zg=z(i)   !-rrzcen
         if(cg(i).eq.0.0d0) goto 101
         if(xg.lt.xmin.or.xg.gt.xmax) goto 101
         if(yg.lt.ymin.or.yg.gt.ymax) goto 101
         if(zg.lt.zmin.or.zg.gt.zmax) goto 101
         xs=xscal*xg
         ys=yscal*yg
         zs=zscal*zg

         lpolx(1)= 1.0d0
         lpoly(1)= 1.0d0
         lpolz(1)= 1.0d0
         lpolx(2)= xs
         lpoly(2)= ys
         lpolz(2)= zs
         lpolx(3)= 0.5d0*(3.0d0*xs*xs-1.0d0)
         lpoly(3)= 0.5d0*(3.0d0*ys*ys-1.0d0)
         lpolz(3)= 0.5d0*(3.0d0*zs*zs-1.0d0)
         do n=3,xnpol-1
            lpolx(n+1)=((2*n-1)*xs*lpolx(n)-(n-1)*lpolx(n-1))/n
         enddo
         do n=3,ynpol-1
            lpoly(n+1)=((2*n-1)*ys*lpoly(n)-(n-1)*lpoly(n-1))/n
         enddo
         do n=3,znpol-1
            lpolz(n+1)=((2*n-1)*zs*lpolz(n)-(n-1)*lpolz(n-1))/n
         enddo

         do ii=1,ntpol
            xpol=lstpx(ii)
            ypol=lstpy(ii)
            zpol=lstpz(ii)
            norm=bnorm(ii)
            coef(ii)=coef(ii)+
     $           cg(i)*norm*
     $           lpolx(xpol+1)*lpoly(ypol+1)*lpolz(zpol+1)
         enddo
 101  enddo

c construct mq array to speed up the calculations
      do ii=1,ntpol
         mq(ii)=0.d0
         do jj=1,ntpol
            ij=(ii-1)*ntpol+jj
            mq(ii)=mq(ii)+mij(ij)*coef(jj)
         enddo
      enddo

c reaction field energy calculation    
      egsbpb=0.d0
      do ii=1,ntpol
         egsbpb=egsbpb+0.5d0*coef(ii)*mq(ii)
      enddo
      egsbpb=egsbpb*celec
c
      write(outu,'(6x,a,/,6x,a,f13.5,a)')
     $     'Reaction field free energy in inner region ',
     $     'calculated from the basis functions     =',
     $     egsbpb,' [kcal/mol]'
c
      return
      end

      subroutine rfrect0(lncel,ltran,ldcel,lxbcen,lybcen,lzbcen,
     $           xscale,yscale,zscale,rrxcen,rrycen,rrzcen)
c-----------------------------------------------------------------------
c     construct the reaction field matrix mmij using the legendre polynomials 
c
      implicit none
      include 'pnp.fcm'
      include 'mmij.fcm'
      include 'mainio.fcm'

      integer*4 lncel
      real*8  ltran,ldcel,lxbcen,lybcen,lzbcen
      real*8  xscale,yscale,zscale
      real*8  rrxcen,rrycen,rrzcen
c local
      integer*4 oldntp,norder
      integer*4 jx1,jx2,jy1,jy2,jz1,jz2,nproi,npgd
      integer*4 lncel3,ncel3,bncel
      logical qmij

c save basis functions
      oldntp=0         ! in this version, we don't consider 
      qmij=.false.     ! any restart of construction of mmij matrix 
      call rect_svpol(xnpol,ynpol,znpol,oldntp,qmij, 
     $     lstpx,lstpy,lstpz,lstpol)

c calculate the normalization constants of the basis functions

      call rect_norm(ntpol,xscale,yscale,zscale,
     $     lstpx,lstpy,lstpz,bnorm)

c obtain all space functions (may,mayx,mayy,mayz)
c array may will be used as the volume exculsion function for ions

      if(qmcden) goto 11

c     initialize all grid functions
      if(kappa2.eq.0.0d0) kappa2=1.0d0
      call initgrid(nclx,ncly,nclz,mcden,kappa2*dcel)
      call initgrid(nclx,ncly,nclz,epsx,epsw)
      call initgrid(nclx,ncly,nclz,epsy,epsw)
      call initgrid(nclx,ncly,nclz,epsz,epsw)

c     membrane
      if(tmemb.gt.0.0d0)then
         call mayer_memb(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        zbcen,tmemb,zmemb,epsm,htmemb,epsh,
     $        mcden,epsx,epsy,epsz,qmcden)
      endif

c     cylinder
      if(rcyln.gt.0.0d0)then
         call mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epsc,xcyln,ycyln,zcyln,rcyln,hcyln,
     $        mcden,epsx,epsy,epsz,qckap,qmcden)
      endif

c here, we just use outer region array (ntra, lstra) to construct
c the dielectric boundary.
      call mayer_mstep(natom,x,y,z,radius,epsp,watr,ionr,
     $     nclx,ncly,nclz,dcel,
     $     tranx,trany,tranz,xbcen,ybcen,zbcen,
     $     mcden,epsx,epsy,epsz,qreen,qmcden)
      if( watr.gt.0.0d0 .and. qreen ) then
         call mayer_reen(natom,x,y,z,radius,epsp,watr,
     $        nclx,ncly,nclz,dcel,
     $        tranx,trany,tranz,xbcen,ybcen,zbcen,
     $        mcden,epsx,epsy,epsz,
     $        pox,poy,poz,mapt,listr,fistr,qmcden)
      endif

c     overlay cylinder

      if(rocyl.gt.0.0d0)then
         call mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epso,xocyl,yocyl,zocyl,rocyl,hocyl,
     $        mcden,epsx,epsy,epsz,qokap,qmcden)
      endif

 11   continue

c find out the grid points inside region of interest
c save them on wa araay
      jx1=int((tranx+rbxmin-xbcen)/dcel)-1
      jx2=int((tranx+rbxmax-xbcen)/dcel)+3
      jy1=int((trany+rbymin-ybcen)/dcel)-1
      jy2=int((trany+rbymax-ybcen)/dcel)+3
      jz1=int((tranz+rbzmin-zbcen)/dcel)-1
      jz2=int((tranz+rbzmax-zbcen)/dcel)+3
      nproi=(jx2-jx1)*(jy2-jy1)*(jz2-jz1)

      call rfrect_gpinrr(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $     xbcen,ybcen,zbcen,mcden,
     $     rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax,srcut,
     $     jx1,jx2,jy1,jy2,jz1,jz2,wa,npgd)

c re-adjust the kappa array
      if(kappa2.eq.1.0d0) then
         kappa2=0.0d0
         call initgrid(nclx,ncly,nclz,mcden,0.0d0)
      else
c     kappa = zero inside the inner region
         if(srcut.gt.0.0d0) then
            call rfsphe_mayer(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,mcden,
     $           rrxcen,rrycen,rrzcen,srcut)
         else
            call rfrect_mayer(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,mcden,
     $           rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax)
         endif
      endif

      if(qmcden) then
         qmcden=.false.
         call initgrid(nclx,ncly,nclz,mcden,0.0d0)
      endif

c setup the overlap matrix (sij) and calculate its inverse
      call rfrect_sij(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $     xbcen,ybcen,zbcen,tolsij,
     $     ntpol,lstpx,lstpy,lstpz,bnorm,
     $     wa,npgd,rrxcen,rrycen,rrzcen,xscale,yscale,zscale,
     $     sij,sijevec,sijeval,tmpwork)

c for charge densities in the large box and boundary potentials in the original box
      lncel3=lncel*lncel*lncel
      ncel3=nclx*ncly*nclz
      bncel=2*(nclx*nclx+ncly*ncly+nclz*nclz)

c ..
c start here for the matrix MMIJ calculations

      do norder=oldntp+1,ntpol

         write(outu,'(a)') 
         write(outu,'(6x,a,i5,a)') 
     $        'Reaction field calculation for',NORDER,' basis function'

c replace FCDEN array with the basis functions (in original box)
         call rect_chc(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $        xbcen,ybcen,zbcen,
     $        norder,lstpx,lstpy,lstpz,bnorm,
     $        fcden,cgscal,wa,npgd,
     $        rrxcen,rrycen,rrzcen,xscale,yscale,zscale)

c here, we will use CION array for charge distribution in the large box 
c interpolated from charge distribution (FCDEN array) in the original box 
         call initgrid(lncel,lncel,lncel,cion,0.0d0)
         call basis_cd_tril(ncel3,cion,kappa2,
     $        lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
     $        lxbcen,lybcen,lzbcen,
     $        fcden,nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen)

c     In solution without membrane ; 
c     ------------------------------
c     KAPPA2 = ZERO; VMEMB  = ZERO; EPSM   = EPSW 
c     EPSH   = EPSW; EPSS   = EPSW; EPSC   = EPSW 

c initialize the (pure solvent) potential array EFFEPS
         call initgrid(nclx,ncly,nclz,effeps,0.0d0) 

c get boundary potentials in the large box
         if(.not.qzerobp) then
            call basis_bp_int(lncel3,epsw,0.0d0,
     $           cion,effeps,lncel,lncel,lncel,ldcel,
     $           ltran,ltran,ltran,lxbcen,lybcen,lzbcen,
     $           lncel,lncel,lncel,ldcel,lxbcen,lybcen,lzbcen,
     $           qphixypbc)
         endif

c solve PBEQ in the large box
         call pbeq4(maxphi,tolphi,epsw,lncel,lncel,lncel,
     $        effeps,cion,.false.,qphixypbc,qphixyzpbc)

c construct boundary potentials from the large box
         call focusphi(lncel,lncel,lncel,ldcel,
     $        ltran,ltran,ltran,lxbcen,lybcen,lzbcen,
     $        nclx,ncly,nclz,dcel,
     $        tranx,trany,tranz,xbcen,ybcen,zbcen,
     $        effeps,phib)

         call mayer_bp_focus(effeps,phib,nclx,ncly,nclz)

c solve PBEQ in the original box
         call pbeq4(maxphi,tolphi,epsw,nclx,ncly,nclz,
     $        effeps,fcden,qphifocus,qphixypbc,qphixyzpbc)

c------------------------------------------------------------------------
c     In solution with membrane (large box - focussing - original box)
c     VMEMB  = ZERO; KAPPA2; EPSM; EPSH; EPSS; EPSC; EPSB
c------------------------------------------------------------------------

c get grid parameters in the large box
         call initgrid(lncel,lncel,lncel,phi,0.0d0)
         call initgrid(lncel,lncel,lncel,mcden,1.0d0)
         call initgrid(lncel,lncel,lncel,epsx,epsw)
         call initgrid(lncel,lncel,lncel,epsy,epsw)
         call initgrid(lncel,lncel,lncel,epsz,epsw)

         if(tmemb.gt.0.0d0)then
            call mayer_memb(lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
     $           lzbcen,tmemb,zmemb,epsm,htmemb,epsh,
     $           mcden,epsx,epsy,epsz,qmcden)
         endif

         if(rsphe.gt.0.0d0)then
            call mayer_sphe(lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
     $           lxbcen,lybcen,lzbcen,
     $           epss,xsphe,ysphe,zsphe,rsphe,
     $           mcden,epsx,epsy,epsz,qskap,qmcden)
         endif

         if(rcyln.gt.0.0d0)then
            call mayer_cyln(lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
     $           lxbcen,lybcen,lzbcen,
     $           epsc,xcyln,ycyln,zcyln,rcyln,hcyln,
     $           mcden,epsx,epsy,epsz,qckap,qmcden)
         endif

         if(bxmax.ne.bxmin)then
            call mayer_box(lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
     $           lxbcen,lybcen,lzbcen,
     $           epsb,bxmax,bymax,bzmax,bxmin,bymin,bzmin,
     $           mcden,epsx,epsy,epsz,qbkap,qmcden)
         endif

         call mayer_mstep(natom,x,y,z,radius,epsp,watr,ionr,
     $        lncel,lncel,lncel,ldcel,
     $        ltran,ltran,ltran,lxbcen,lybcen,lzbcen,
     $        mcden,epsx,epsy,epsz,qreen,qmcden)

         if( watr.gt.0.0d0 .and. qreen ) then
            call mayer_reen(natom,x,y,z,radius,epsp,watr,
     $           lncel,lncel,lncel,ldcel,
     $           ltran,ltran,ltran,lxbcen,lybcen,lzbcen,
     $           mcden,epsx,epsy,epsz,
     $           pox,poy,poz,mapt,listr,fistr,qmcden)
         endif

         if(rocyl.gt.0.0d0)then
            call mayer_cyln(lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
     $           lxbcen,lybcen,lzbcen,
     $           epso,xocyl,yocyl,zocyl,rocyl,hocyl,
     $           mcden,epsx,epsy,epsz,qokap,qmcden)
         endif

c kappa = zero inside the inner region
         if(kappa2.ne.0.0d0) then
            if(srcut.gt.0.0d0) then
               call rfsphe_mayer(lncel,lncel,lncel,ltran,ltran,ltran,
     $              ldcel,lxbcen,lybcen,lzbcen,mcden,
     $              rrxcen,rrycen,rrzcen,srcut)
            else
               call rfrect_mayer(lncel,lncel,lncel,ltran,ltran,ltran,
     $              ldcel,lxbcen,lybcen,lzbcen,mcden,
     $              rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax)
            endif
         endif

c get boundary potentials in the large box
         if(.not.qzerobp) then
            call basis_bp_int(lncel3,epsw,kappa2,
     $           cion,phi,lncel,lncel,lncel,ldcel,
     $           ltran,ltran,ltran,lxbcen,lybcen,lzbcen,
     $           lncel,lncel,lncel,ldcel,lxbcen,lybcen,lzbcen,
     $           qphixypbc)
         endif

c solve pbeq in the large box
         call pbeq1(maxphi,tolphi,kappa2,kap2top,kap2bot,
     $        lncel,lncel,lncel,ldcel,tmemb,zmemb,ltran,lzbcen,
     $        phi,epsx,epsy,epsz,cion,mcden,qphixypbc,qphixyzpbc)

         call initgrid(nclx,ncly,nclz,mcden,1.0d0)
         call initgrid(nclx,ncly,nclz,epsx,epsw)
         call initgrid(nclx,ncly,nclz,epsy,epsw)
         call initgrid(nclx,ncly,nclz,epsz,epsw)

c construct boundary potentials from the large box
         call focusphi(lncel,lncel,lncel,ldcel,
     $        ltran,ltran,ltran,lxbcen,lybcen,lzbcen,
     $        nclx,ncly,nclz,dcel,
     $        tranx,trany,tranz,xbcen,ybcen,zbcen,
     $        phi,phib)

         call initgrid(nclx,ncly,nclz,phi,0.0d0)
         call mayer_bp_focus(phi,phib,nclx,ncly,nclz)

         if(tmemb.gt.0.0d0)then
            call mayer_memb(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           zbcen,tmemb,zmemb,epsm,htmemb,epsh,
     $           mcden,epsx,epsy,epsz,qmcden)
         endif

         if(rsphe.gt.0.0d0)then
            call mayer_sphe(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epss,xsphe,ysphe,zsphe,rsphe,
     $           mcden,epsx,epsy,epsz,qskap,qmcden)
         endif

         if(rcyln.gt.0.0d0)then
            call mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epsc,xcyln,ycyln,zcyln,rcyln,hcyln,
     $           mcden,epsx,epsy,epsz,qckap,qmcden)
         endif

         if(bxmax.ne.bxmin)then
            call mayer_box(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epsb,bxmax,bymax,bzmax,bxmin,bymin,bzmin,
     $           mcden,epsx,epsy,epsz,qbkap,qmcden)
         endif

         call mayer_mstep(natom,x,y,z,radius,epsp,watr,ionr,
     $        nclx,ncly,nclz,dcel,
     $        tranx,trany,tranz,xbcen,ybcen,zbcen,
     $        mcden,epsx,epsy,epsz,qreen,qmcden)

         if( watr.gt.0.0d0 .and. qreen ) then
            call mayer_reen(natom,x,y,z,radius,epsp,watr,
     $           nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           mcden,epsx,epsy,epsz,
     $           pox,poy,poz,mapt,listr,fistr,qmcden)
         endif

         if(rocyl.gt.0.0d0)then
            call mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epso,xocyl,yocyl,zocyl,rocyl,hocyl,
     $           mcden,epsx,epsy,epsz,qokap,qmcden)
         endif

c kappa = zero inside the inner region
         if(kappa2.ne.0.0d0) then
            if(srcut.gt.0.0d0) then
               call rfsphe_mayer(nclx,ncly,nclz,tranx,trany,tranz,
     $              dcel,xbcen,ybcen,zbcen,mcden,
     $              rrxcen,rrycen,rrzcen,srcut)
            else
               call rfrect_mayer(nclx,ncly,nclz,tranx,trany,tranz,
     $              dcel,xbcen,ybcen,zbcen,mcden,
     $              rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax)
            endif
         endif

c solve pbeq in the original box
         call pbeq1(maxphi,tolphi,kappa2,kap2top,kap2bot,
     $        nclx,ncly,nclz,dcel,tmemb,zmemb,tranz,zbcen,
     $        phi,epsx,epsy,epsz,fcden,mcden,.false.,.false.)
c                                          qphixypbc  qphixyzpbc

c         call tmpwrigd(80,nclx,ncly,nclz,dcel,tranx,trany,tranz,
c     $        xbcen,ybcen,zbcen,-200.d0,0.d0,0.d0,200.d0,0.d0,0.d0,
c     $        phi,1.0d0)
c         call tmpwrigd(81,nclx,ncly,nclz,dcel,tranx,trany,tranz,
c     $        xbcen,ybcen,zbcen,0.d0,-200.d0,0.d0,0.d0,200.d0,0.d0,
c     $        phi,1.0d0)
c         call tmpwrigd(82,nclx,ncly,nclz,dcel,tranx,trany,tranz,
c     $        xbcen,ybcen,zbcen,0.d0,0.d0,-200.d0,0.d0,0.d0,200.d0,
c     $        phi,1.0d0)
c         call tmpwrigd(80,lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
c     $        lxbcen,lybcen,lzbcen,-200.d0,0.d0,0.d0,200.d0,0.d0,0.d0,
c     $        phi,1.0d0)
c         call tmpwrigd(81,lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
c     $        lxbcen,lybcen,lzbcen,0.d0,-200.d0,0.d0,0.d0,200.d0,0.d0,
c     $        phi,1.0d0)
c         call tmpwrigd(82,lncel,lncel,lncel,ldcel,ltran,ltran,ltran,
c     $        lxbcen,lybcen,lzbcen,0.d0,0.d0,-200.d0,0.d0,0.d0,200.d0,
c     $        phi,1.0d0)

c calculate the reaction field energy elements m(ij)
         call rect_mij(nclx,ncly,nclz,tranx,trany,tranz,ntpol,norder,
     $        lstpx,lstpy,lstpz,bnorm,dcel,xbcen,ybcen,zbcen,
     $        rrxcen,rrycen,rrzcen,xscale,yscale,zscale,wa,npgd,
     $        effeps,phi,mij,cgscal)

      enddo

      call comp_mmij(ntpol,mmij,mij,sij)

c
      return
      end

      subroutine rect_mij(nclx,ncly,nclz,tranx,trany,tranz,ntpol,norder,
     $           lstpx,lstpy,lstpz,bnorm,dcel,xbcen,ybcen,zbcen,
     $           rrxcen,rrycen,rrzcen,xscale,yscale,zscale,iip0,npgd,
     $           phip,phiw,mij,cgscal)
c-----------------------------------------------------------------------
c     calculate the matrix m(ij)
c
      implicit none

      real*4  phip(*),phiw(*)
      real*8  rrxcen,rrycen,rrzcen,cgscal,xbcen,ybcen,zbcen
      real*8  dcel,tranx,trany,tranz,xscale,yscale,zscale
      real*8  mij(*),bnorm(*)
      integer*4 lstpx(*),lstpy(*),lstpz(*),norder,ntpol
      integer*4 nclx,ncly,nclz,npgd,iip0(*)
c local
      real*8  lpol,xg,yg,zg,dcel3,norm
      real*8  xadd,yadd,zadd
      integer*4 nc3,ncyz,i,ig,jg,kg,ip0
      integer*4 xp,yp,zp,n,ij
c
      ncyz  = ncly*nclz
      nc3   = ncyz*nclx
      dcel3 = dcel*dcel*dcel

      xadd = -tranx+xbcen-rrxcen
      yadd = -trany+ybcen-rrycen
      zadd = -tranz+zbcen-rrzcen
c
      if(norder.eq.1) then
         n=1
         xp=lstpx(n)
         yp=lstpy(n)
         zp=lstpz(n)
         norm=bnorm(n)
         ij=(n-1)*ntpol+norder
         mij(ij)=0.0d0
         do i=1,npgd
            ip0=iip0(i)
            ig=int((ip0-1)/ncyz)+1
            jg=int(mod((ip0-1),ncyz)/nclz)+1
            kg=mod(mod((ip0-1),ncyz),nclz)+1
            xg=dcel*(ig-1)+xadd
            yg=dcel*(jg-1)+yadd
            zg=dcel*(kg-1)+zadd
            mij(ij)=mij(ij)+norm*lpol(xp,xg,xscale)*
     $                      lpol(yp,yg,yscale)*lpol(zp,zg,zscale)*
     $                      dcel3*(phiw(ip0)-phip(ip0))*cgscal
         enddo
      else
         do n=1,norder
            xp=lstpx(n)
            yp=lstpy(n)
            zp=lstpz(n)
            norm=bnorm(n)
            ij=(n-1)*ntpol+norder
            mij(ij)=0.0d0
            do i=1,npgd
               ip0=iip0(i)
               ig=int((ip0-1)/ncyz)+1
               jg=int(mod((ip0-1),ncyz)/nclz)+1
               kg=mod(mod((ip0-1),ncyz),nclz)+1
               xg=dcel*(ig-1)+xadd
               yg=dcel*(jg-1)+yadd
               zg=dcel*(kg-1)+zadd
               mij(ij)=mij(ij)+norm*lpol(xp,xg,xscale)*
     $                         lpol(yp,yg,yscale)*lpol(zp,zg,zscale)*
     $                         dcel3*(phiw(ip0)-phip(ip0))
            enddo
c            write(*,*) n,mij(ij)
         enddo
      endif
c
      return
      end


      subroutine comp_mmij(ntpol,mmij,mij,xxij)
c------------------------------------------------------------------------
c     MMIJ=XXIJ*MIJ*XXIJ
c
      implicit none

      integer*4 ntpol
      real*8  mmij(*),mij(*),xxij(ntpol,ntpol)
c
      integer*4 i,j,ij,ji,m,n,mn
      real*8  sum,xim
c
c      read(90) ((xxij(i,j),j=1,ntpol),i=1,ntpol)
c      read(91) ((mij((i-1)*ntpol+j),j=1,ntpol),i=1,ntpol)

c      write(6,'(10f8.3)') ((xxij(i,j),j=1,10),i=1,10)
c      write(6,*)
c      write(90) ((xxij(i,j),j=1,ntpol),i=1,ntpol)

c
      do i=1,ntpol
         do j=1,ntpol
            if(i.gt.j) then
               ij=(i-1)*ntpol+j
               ji=(j-1)*ntpol+i
               mij(ij)=mij(ji)
            endif
         enddo
      enddo

c      write(6,'(10f10.5)') ((mij((i-1)*ntpol+j),j=1,10),i=1,10)
c      write(6,*)
c      write(91) ((mij((i-1)*ntpol+j),j=1,ntpol),i=1,ntpol)

c
      do m=1,ntpol
         do n=1,ntpol
            mn=(m-1)*ntpol+n
            sum=0.0d0
            do i=1,ntpol
               xim=xxij(i,m)
               do j=1,ntpol
                  ij=(i-1)*ntpol+j
                  sum=sum+xim*mij(ij)*xxij(j,n)
               enddo
            enddo
            mmij(mn)=sum
         enddo
      enddo

c      write(6,'(10f10.5)') ((mmij((i-1)*ntpol+j),j=1,10),i=1,10)
c      write(6,*)
c
      return
      end


      subroutine basis_bp_int(ntprp,epsw,kappa2,
     $           chc,phi,nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           onclx,oncly,onclz,odcel,oxbcen,oybcen,ozbcen,
     $           qphixypbc)
c----------------------------------------------------------------------
c     half number of grid points is used for the debye-huckel approximation
c
      implicit none
      include 'consta.fcm'

      integer*4  ntprp
      integer*4  nclx,ncly,nclz,onclx,oncly,onclz
      real*4   phi(*),chc(*)
      real*8   dcel,odcel,oxbcen,oybcen,ozbcen
      real*8   tranx,trany,tranz,xbcen,ybcen,zbcen,epsw,kappa2
      logical  qphixypbc
c local
      integer*4  il,ipo,ipox,ipoy,ig,jg,kg,ncyz,oncyz
      integer*4  l,m,n
      real*8   otranx,otrany,otranz
      real*8   chi,xi,yi,zi,kappa,xj,yj,zj,dist1,dist2
c
      kappa=sqrt(kappa2/epsw)
      ncyz=ncly*nclz
      oncyz=oncly*onclz
      otranx=0.5d0*(onclx-1)*odcel
      otrany=0.5d0*(oncly-1)*odcel
      otranz=0.5d0*(onclz-1)*odcel

c to consider the boundary conditions of basis set charge distributions
      do 445 il=1,ntprp
         chi=chc(il)*(odcel/2.0d0/twopi)
         if(chi.eq.0.0d0) goto 445
         l=int((il-1)/oncyz)+1
         m=int(mod((il-1),oncyz)/onclz)+1
         n=mod(mod((il-1),oncyz),onclz)+1
         xi=(l-1)*odcel-otranx+oxbcen+tranx-xbcen
         yi=(m-1)*odcel-otrany+oybcen+trany-ybcen
         zi=(n-1)*odcel-otranz+ozbcen+tranz-zbcen

         if(.not.qphixypbc)then

c yz plane
            ipox=(nclx-1)*ncyz
            xj=(nclx-1)*dcel-xi
            do jg=1,ncly,2
               ipoy=(jg-1)*nclz
               yj=(jg-1)*dcel-yi
               do kg=3,nclz-2,2
                  zj=(kg-1)*dcel-zi
                  
                  ipo=ipoy+kg
                  dist1=sqrt(xi*xi+yj*yj+zj*zj)
                  if(kappa.eq.0.0d0) then
                     phi(ipo)=phi(ipo)+chi/epsw/dist1
                  else
                     phi(ipo)=phi(ipo)+chi*exp(-kappa*dist1)/epsw/dist1
                  endif

                  ipo=ipox+ipo
                  dist2=sqrt(xj*xj+yj*yj+zj*zj)
                  if(kappa.eq.0.0d0) then
                     phi(ipo)=phi(ipo)+chi/epsw/dist2
                  else
                     phi(ipo)=phi(ipo)+chi*exp(-kappa*dist2)/epsw/dist2
                  endif

               enddo
            enddo

c xz plane
            yj=(ncly-1)*dcel-yi
            do ig=3,nclx-2,2
               ipox=(ig-1)*ncyz
               xj=(ig-1)*dcel-xi
               do kg=3,nclz-2,2
                  zj=(kg-1)*dcel-zi

                  ipo=ipox+kg
                  dist1=sqrt(xj*xj+yi*yi+zj*zj)
                  if(kappa.eq.0.0d0) then
                     phi(ipo)=phi(ipo)+chi/epsw/dist1
                  else
                     phi(ipo)=phi(ipo)+chi*exp(-kappa*dist1)/epsw/dist1
                  endif

                  ipo=ipo+ncyz-nclz
                  dist2=sqrt(xj*xj+yj*yj+zj*zj)
                  if(kappa.eq.0.0d0) then
                     phi(ipo)=phi(ipo)+chi/epsw/dist2
                  else
                     phi(ipo)=phi(ipo)+chi*exp(-kappa*dist2)/epsw/dist2
                  endif

               enddo
            enddo
         endif

c xy plane
         zj=(nclz-1)*dcel-zi
         do ig=1,nclx,2
            ipox=(ig-1)*ncyz
            xj=(ig-1)*dcel-xi
            do jg=1,ncly,2
               yj=(jg-1)*dcel-yi

               ipo=ipox+(jg-1)*nclz+1
               dist1=sqrt(xj*xj+yj*yj+zi*zi)
               if(kappa.eq.0.0d0) then
                  phi(ipo)=phi(ipo)+chi/epsw/dist1
               else
                  phi(ipo)=phi(ipo)+chi*exp(-kappa*dist1)/epsw/dist1
               endif
               
               ipo=ipo+nclz-1
               dist2=sqrt(xj*xj+yj*yj+zj*zj)
               if(kappa.eq.0.0d0) then
                  phi(ipo)=phi(ipo)+chi/epsw/dist2
               else
                  phi(ipo)=phi(ipo)+chi*exp(-kappa*dist2)/epsw/dist2
               endif

            enddo
         enddo
 445  enddo 

c bilinear interpolation

c xy plane
c elements in odd-numbered rows, interplating horizontally
      do ig=2,nclx-1,2
         ipox=(ig-1)*ncyz
         do jg=1,ncly,2
            ipo=ipox+(jg-1)*nclz+1
            phi(ipo)=(phi(ipo-ncyz)+phi(ipo+ncyz))/2.0d0
            ipo=ipo+nclz-1
            phi(ipo)=(phi(ipo-ncyz)+phi(ipo+ncyz))/2.0d0
         enddo
      enddo
c elements in even-numbered rows, interplating vertically
      do ig=1,nclx
         ipox=(ig-1)*ncyz
         do jg=2,ncly-1,2
            ipo=ipox+(jg-1)*nclz+1
            phi(ipo)=(phi(ipo-nclz)+phi(ipo+nclz))/2.0d0
            ipo=ipo+nclz-1
            phi(ipo)=(phi(ipo-nclz)+phi(ipo+nclz))/2.0d0
         enddo
      enddo

      if(.not.qphixypbc)then
c yz plane
c elements in odd-numbered rows, interplating horizontally
         ipox=(nclx-1)*ncyz
         do jg=2,ncly-1,2
            ipoy=(jg-1)*nclz
            do kg=3,nclz-2,2
               ipo=ipoy+kg
               phi(ipo)=(phi(ipo-nclz)+phi(ipo+nclz))/2.0d0
               ipo=ipox+ipo
               phi(ipo)=(phi(ipo-nclz)+phi(ipo+nclz))/2.0d0
            enddo
         enddo
c elements in even-numbered rows, interplating vertically
         do jg=1,ncly
            ipoy=(jg-1)*nclz
            do kg=2,nclz,2
               ipo=ipoy+kg
               phi(ipo)=(phi(ipo-1)+phi(ipo+1))/2.0d0
               ipo=ipox+ipo
               phi(ipo)=(phi(ipo-1)+phi(ipo+1))/2.0d0
            enddo
         enddo
c xz plane
c elements in odd-numbered rows, interplating horizontally
         do ig=2,nclx-1,2
            ipox=(ig-1)*ncyz
            do kg=3,nclz-2,2
               ipo=ipox+kg
               phi(ipo)=(phi(ipo-ncyz)+phi(ipo+ncyz))/2.0d0
               ipo=ipo+ncyz-nclz
               phi(ipo)=(phi(ipo-ncyz)+phi(ipo+ncyz))/2.0d0
            enddo
         enddo
c elements in even-numbered rows, interplating vertically
         do ig=2,nclx-1
            ipox=(ig-1)*ncyz
            do kg=2,nclz-1,2
               ipo=ipox+kg
               phi(ipo)=(phi(ipo-1)+phi(ipo+1))/2.0d0
               ipo=ipo+ncyz-nclz
               phi(ipo)=(phi(ipo-1)+phi(ipo+1))/2.0d0
            enddo
         enddo
      endif
c
      return
      end


      subroutine basis_cd_tril(ntprp,chc,kappa2,
     $           nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           ochc,onclx,oncly,onclz,odcel,oxbcen,oybcen,ozbcen)
c----------------------------------------------------------------------
c     the trilinear interpolation for charge distribution 
c
      implicit none
      include  'consta.fcm'

      integer*4  ntprp
      integer*4  nclx,ncly,nclz,onclx,oncly,onclz
      real*4   ochc(*),chc(*)
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen,kappa2
      real*8   odcel,oxbcen,oybcen,ozbcen
c local
      integer*4  il,n1,n2,n3,in1,in2,in3,l,m,n
      integer*4  ix,iy,iz,ncyz,oncyz
      real*8   otranx,otrany,otranz
      real*8   chi,xi,yi,zi,ai,bi,ci,fi
c
      ncyz =ncly*nclz
      oncyz=oncly*onclz
      otranx=0.5d0*(onclx-1)*odcel
      otrany=0.5d0*(oncly-1)*odcel
      otranz=0.5d0*(onclz-1)*odcel

c to consider the boundary conditions of basis set charge distributions
      do 460 il=1,ntprp
         chi=ochc(il)*(odcel/2.0d0/twopi)
         if(chi.eq.0.0d0) goto 460
         l=int((il-1)/oncyz)+1
         m=int(mod((il-1),oncyz)/onclz)+1
         n=mod(mod((il-1),oncyz),onclz)+1
         xi=(l-1)*odcel-otranx+oxbcen+tranx-xbcen
         yi=(m-1)*odcel-otrany+oybcen+trany-ybcen
         zi=(n-1)*odcel-otranz+ozbcen+tranz-zbcen

         ix=int(xi/dcel)+1
         iy=int(yi/dcel)+1
         iz=int(zi/dcel)+1
            
c Construct the charge distribution by 8 grid points adjacent to the atom
         do n1=1,2
            in1=ix+n1-1
            ai=abs(xi-(in1-1)*dcel)/dcel
            in1=(in1-1)*ncyz
            ai=1.0d0-ai
c     
            do n2=1,2
               in2=iy+n2-1
               bi=abs(yi-(in2-1)*dcel)/dcel
               in2=(in2-1)*nclz
               bi=ai*(1.0d0-bi)
c     
               do n3=1,2
                  in3=iz+n3-1
                  ci=abs(zi-(in3-1)*dcel)/dcel
                  fi=bi*(1.0d0-ci)
                  in3=in1+in2+in3 

                  chc(in3)=chc(in3)+(fi*chi)*2.0d0*twopi/dcel

               enddo
            enddo
         enddo
 460  enddo
c
      return
      end


      subroutine rect_chc(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,
     $           norder,lstpx,lstpy,lstpz,bnorm,
     $           chc,cgscal,iip0,npgd,
     $           rrxcen,rrycen,rrzcen,xscale,yscale,zscale)
c-----------------------------------------------------------------------
c     store the basis functions (legendre polynomials) in ichc (cden) array
c     note: chc is a charge (not density) arrary and containes 4*pi/h
c
      implicit none
      include 'consta.fcm'

      real*4  chc(*)
      real*8  rrxcen,rrycen,rrzcen,bnorm(*)
      real*8  dcel,tranx,trany,tranz,xscale,yscale,zscale,cgscal
      real*8  xbcen,ybcen,zbcen
      integer*4 nclx,ncly,nclz,npgd,iip0(*)
      integer*4 lstpx(*),lstpy(*),lstpz(*),norder
c local
      real*8  lpol,xg,yg,zg,dcel3,norm
      real*8  xadd,yadd,zadd,factor
      integer*4 nc3,ncyz,i,ig,jg,kg,ip0,xpol,ypol,zpol
c
      xpol=lstpx(norder)
      ypol=lstpy(norder)
      zpol=lstpz(norder)
      norm=bnorm(norder)
c
      ncyz =ncly*nclz
      nc3  =ncyz*nclx
      dcel3=dcel*dcel*dcel

c initialize chc array
      do i=1,nc3
         chc(i)=0.0d0
      enddo

      xadd = -tranx+xbcen-rrxcen
      yadd = -trany+ybcen-rrycen
      zadd = -tranz+zbcen-rrzcen

c the legendre polynomials are normalized inside region of interest
      if(norder.eq.1) then
         factor=dcel3*2.0d0*twopi/dcel/cgscal
         do i=1,npgd
            ip0=iip0(i)
            ig=int((ip0-1)/ncyz)+1
            jg=int(mod((ip0-1),ncyz)/nclz)+1
            kg=mod(mod((ip0-1),ncyz),nclz)+1
            xg=dcel*(ig-1)+xadd
            yg=dcel*(jg-1)+yadd
            zg=dcel*(kg-1)+zadd
            chc(ip0)=norm*lpol(xpol,xg,xscale)*lpol(ypol,yg,yscale)*
     $               lpol(zpol,zg,zscale)*factor
         enddo
      else
         factor=dcel3*2.0d0*twopi/dcel
         do i=1,npgd
            ip0=iip0(i)
            ig=int((ip0-1)/ncyz)+1
            jg=int(mod((ip0-1),ncyz)/nclz)+1
            kg=mod(mod((ip0-1),ncyz),nclz)+1
            xg=dcel*(ig-1)+xadd
            yg=dcel*(jg-1)+yadd
            zg=dcel*(kg-1)+zadd
            chc(ip0)=norm*lpol(xpol,xg,xscale)*lpol(ypol,yg,yscale)*
     $               lpol(zpol,zg,zscale)*factor
         enddo
      endif
c
c      do i=1,npgd
c         ip0=iip0(i)
c         ig=int((ip0-1)/ncyz)+1
c         jg=int(mod((ip0-1),ncyz)/nclz)+1
c         kg=mod(mod((ip0-1),ncyz),nclz)+1
c         xg=(dcel*(ig-1)-tranx-rrxcen)*10.0
c         yg=(dcel*(jg-1)-trany-rrycen)*10.0
c         zg=(dcel*(kg-1)-tranz-rrzcen)*10.0
c         if(chc(ip0).gt.rsmall) then
c            write(50+norder,101)
c     $           'atom',i,' pol arg     1',
c     $           xg,yg,zg,
c     $           chc(ip0),'  1.00      lpol'
c         elseif(chc(ip0).lt.-rsmall) then
c            write(50+norder,101)
c     $           'atom',i,' pol glu     2',
c     $           xg,yg,zg,
c     $           chc(ip0),' 10.00      lpol'
c         else
c            write(50+norder,101)
c     $           'atom',i,' pol cys     3',
c     $           xg,yg,zg,
c     $           chc(ip0),' 10.00      lpol'
c         endif
c      enddo
c 101  format(a,2x,i5,1x,a,4x,3f8.3,f6.2,a)
c
      return
      end


      subroutine rfrect_sij(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,tolsij,
     $           ntpol,lstpx,lstpy,lstpz,bnorm,
     $           iip0,npgd,rrxcen,rrycen,rrzcen,xscale,yscale,zscale,
     $           sij,sijevec,sijeval,tmpwork)
c-----------------------------------------------------------------------
c     setup the overlap matrix (sij)
c
      implicit none
      include 'consta.fcm'

      integer*4 nclx,ncly,nclz,npgd,iip0(*)
      integer*4 lstpx(*),lstpy(*),lstpz(*),ntpol
      real*8  rrxcen,rrycen,rrzcen,bnorm(*)
      real*8  dcel,tranx,trany,tranz,xscale,yscale,zscale
      real*8  xbcen,ybcen,zbcen,tolsij
      real*8  sij(ntpol,ntpol),sijevec(ntpol,ntpol)
      real*8  sijeval(ntpol),tmpwork(ntpol)
c local
      integer*4 i,j,k,ncyz,nc3,ip0,ig,jg,kg,noff
      integer*4 xpol,ypol,zpol
      real*8  dcel3,xg,yg,zg,norm,b(1000),lpol
      real*8  xadd,yadd,zadd,s
c
      if(ntpol.gt.1000) then
         stop 'Increase dimension of B in NRECT_SIJ'
         return
      endif
c
      do i=1,ntpol
         do j=1,ntpol
            sij(j,i)=0.0d0
         enddo
      enddo
      ncyz =ncly*nclz
      nc3  =ncyz*nclx
      dcel3=dcel*dcel*dcel

      xadd = -tranx+xbcen-rrxcen
      yadd = -trany+ybcen-rrycen
      zadd = -tranz+zbcen-rrzcen
c
      do k=1,npgd
         ip0=iip0(k)
         ig=int((ip0-1)/ncyz)+1
         jg=int(mod((ip0-1),ncyz)/nclz)+1
         kg=mod(mod((ip0-1),ncyz),nclz)+1
         xg=dcel*(ig-1) + xadd
         yg=dcel*(jg-1) + yadd
         zg=dcel*(kg-1) + zadd
         do i=1,ntpol
            xpol=lstpx(i)
            ypol=lstpy(i)
            zpol=lstpz(i)
            norm=bnorm(i)
            b(i)=norm*lpol(xpol,xg,xscale)*lpol(ypol,yg,yscale)*
     $                lpol(zpol,zg,zscale)
         enddo
         do i=1,ntpol
            do j=i,ntpol
               sij(i,j)=sij(i,j)+b(i)*b(j)*dcel3
            enddo
         enddo
      enddo

      do i=1,ntpol
         do j=i,ntpol
            s=sij(i,j)
            if(s.lt.rsmall.and.s.gt.-rsmall) then
               sij(i,j)=0.0d0
               s=0.0d0
            endif
            sij(j,i)=s
         enddo
      enddo
c      write(*,*) ((sij(i,j),j=1,ntpol),i=1,ntpol)

c
c     get the transformation matrix x=u*s^(-1/2)
c
      j=0
c get the eigenvectors (sijevec) and eigenvalues (sijeval)
c of the overalp matrix
      call eigrs(sij,ntpol,11,sijeval,sijevec,ntpol,tmpwork,j)
      if (j .gt. 128) then
         j = j - 128
         write (*,'(a,i6,a)')
     *        ' DIAGRS> Failed to converge on root number ', J, '.'
      endif
      write(6,*) 
      write(6,'(6x,a)') 'Sorted eigenvalues of the overlap matrix SIJ;' 
      do i=1,ntpol
         write(6,'(6x,i5,f10.6)') i,sijeval(i)
      enddo

      noff=0
      do i=1,ntpol
         if(sijeval(i).gt.tolsij) then
            noff=i-1
            write(6,*)
            write(6,'(6x,a,i5,a,e10.3)') 'NOFF =',NOFF,' TOLSIJ=',TOLSIJ
            goto 90
         endif
      enddo
 90   continue

c change the order of eigenvalues and column of eigenvectors
      do i=1,ntpol
         s=sqrt(sijeval(ntpol+1-i))
         do j=1,ntpol
            sij(j,i)=sijevec(j,ntpol+1-i)/ s
         enddo
      enddo

c construct x*x^t using sijevec as a temporary array
      do i=1,ntpol
         do j=1,ntpol
            s=0.0d0
            do k=1,ntpol-noff
               s=s+sij(i,k)*sij(j,k)
            enddo
            sijevec(i,j)=s
         enddo
      enddo
c      write(6,'(10f8.3)') ((sijevec(i,j),j=1,10),i=1,10)
c      write(*,*)

c replace sijevec with sij
      do i=1,ntpol
         do j=1,ntpol
            sij(i,j)=sijevec(i,j)
         enddo
      enddo

c      write(6,'(10f8.3)') ((sij(i,j),j=1,10),i=1,10)
c      write(*,*)
c
      return
      end


      subroutine rfrect_mayer(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,fkapa,
     $           rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax)
c-----------------------------------------------------------------------
c     kappa = 0 inside the inner region
c
      implicit none
      include 'consta.fcm'

      real*4  fkapa(*)
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8  rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax
      integer*4 nclx,ncly,nclz
c local
      real*8  xg,yg,zg
      integer*4 ncyz,ig,jg,kg,ip0x,ip0y,ip0
      integer*4 jx1,jx2,jy1,jy2,jz1,jz2
c
      ncyz  = ncly*nclz
c
      jx1=int((tranx+rbxmin-xbcen)/dcel)-1
      jx2=int((tranx+rbxmax-xbcen)/dcel)+3
      jy1=int((trany+rbymin-ybcen)/dcel)-1
      jy2=int((trany+rbymax-ybcen)/dcel)+3
      jz1=int((tranz+rbzmin-zbcen)/dcel)-1
      jz2=int((tranz+rbzmax-zbcen)/dcel)+3

c for debye-huckel screening factors
      do 101 ig=jx1,jx2
         xg=dcel*(ig-1)-tranx+xbcen
         if(xg.lt.rbxmin-rsmall.or.xg.gt.rbxmax+rsmall) goto 101
         ip0x=(ig-1)*ncyz
         do 102 jg=jy1,jy2
            yg=dcel*(jg-1)-trany+ybcen
            if(yg.lt.rbymin-rsmall.or.yg.gt.rbymax+rsmall) goto 102
            ip0y=(jg-1)*nclz+ip0x
            do 103 kg=jz1,jz2
               zg=dcel*(kg-1)-tranz+zbcen
               if(zg.lt.rbzmin-rsmall.or.zg.gt.rbzmax+rsmall) goto 103
               ip0 = ip0y + kg

               fkapa(ip0)=0.0d0

 103        enddo
 102     enddo
 101  enddo
c
      return
      end


      subroutine rfsphe_mayer(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,fkapa,
     $           rrxcen,rrycen,rrzcen,srdist)
c-----------------------------------------------------------------------
c     kappa = 0 inside the inner region
c
      implicit none
      include 'consta.fcm'

      real*4  fkapa(*)
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8  rrxcen,rrycen,rrzcen,srdist
      integer*4 nclx,ncly,nclz
c local
      real*8  xg,yg,zg,xg2,yg2,zg2,dsq1,sr2
      integer*4 IG,JG,KG
      integer*4 IP0,IP0X,IP0Y,NCYZ
      integer*4 NFIL,IXCEN,IYCEN,IZCEN,JX1,JX2,JY1,JY2,JZ1,JZ2

      NFIL=INT(SRDIST/DCEL)+2
      IXCEN=INT((RRXCEN+TRANX-XBCEN)/DCEL)+1
      IYCEN=INT((RRYCEN+TRANY-YBCEN)/DCEL)+1
      IZCEN=INT((RRZCEN+TRANZ-ZBCEN)/DCEL)+1
      JX1=IXCEN-NFIL 
      JX2=IXCEN+NFIL+1
      JY1=IYCEN-NFIL
      JY2=IYCEN+NFIL+1
      JZ1=IZCEN-NFIL
      JZ2=IZCEN+NFIL+1

c For dielectric constants and Debye-Huckel screening factors
c NOTE: dielectric constants are located at the middle point between grids
      ncyz= ncly*nclz
      sr2=srdist*srdist+rsmall
      do ig=jx1,jx2
         xg=dcel*(ig-1)-tranx+xbcen-rrxcen
         xg2=xg*xg
         ip0x=(ig-1)*ncyz
         do jg=jy1,jy2
            yg=dcel*(jg-1)-trany+ybcen-rrycen
            yg2=yg*yg
            ip0y=(jg-1)*nclz+ip0x
            do kg=jz1,jz2
               zg=dcel*(kg-1)-tranz+zbcen-rrzcen
               zg2=zg*zg
               ip0 = ip0y + kg
               if(fkapa(ip0).ne.0.0d0)then
                  dsq1=xg2+yg2+zg2
                  if(dsq1.le.sr2) fkapa(ip0)=0.0d0
               endif
            enddo
         enddo
      enddo
c
      return
      end


      subroutine rfrect_gpinrr(nclx,ncly,nclz,tranx,trany,tranz,dcel,
     $           xbcen,ybcen,zbcen,vexcl,
     $           rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax,srcut,
     $           jx1,jx2,jy1,jy2,jz1,jz2,iip0,npgd)
c-----------------------------------------------------------------------
c     find out the grid points inside region of interest and
c     store in IIP0 array
c
      implicit none
      real*4  vexcl(*)
      real*8  rbxmin,rbxmax,rbymin,rbymax,rbzmin,rbzmax,srcut
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      integer*4 nclx,ncly,nclz,npgd,iip0(*)
      integer*4 jx1,jx2,jy1,jy2,jz1,jz2
c local
      real*8  xg,yg,zg,rrxcen,rrycen,rrzcen,sr2,r2
      integer*4 ncyz,ig,jg,kg,ip0x,ip0y,ip0
c
      rrxcen = (rbxmax+rbxmin)/2.0d0
      rrycen = (rbymax+rbymin)/2.0d0
      rrzcen = (rbzmax+rbzmin)/2.0d0
      sr2=srcut*srcut
c
      npgd=0
      ncyz =ncly*nclz
      do 101 ig=jx1,jx2
         xg=dcel*(ig-1)-tranx+xbcen
         if(xg.lt.rbxmin.or.xg.gt.rbxmax) goto 101
         ip0x=(ig-1)*ncyz
         do 102 jg=jy1,jy2
            yg=dcel*(jg-1)-trany+ybcen
            if(yg.lt.rbymin.or.yg.gt.rbymax) goto 102
            ip0y=(jg-1)*nclz+ip0x
            do 103 kg=jz1,jz2
               zg=dcel*(kg-1)-tranz+zbcen
               if(zg.lt.rbzmin.or.zg.gt.rbzmax) goto 103
               ip0 = ip0y + kg
               
               if(vexcl(ip0).gt.0.0d0) then
                  if(srcut.gt.0.0d0) then
                     xg=xg-rrxcen
                     yg=yg-rrycen
                     zg=zg-rrzcen
                     r2=xg*xg+yg*yg+zg*zg
                     if(r2.le.sr2) then
                        npgd=npgd+1
                        iip0(npgd)=ip0
                     endif
                  else
                     npgd=npgd+1
                     iip0(npgd)=ip0
                  endif
               endif

 103        enddo
 102     enddo
 101  enddo
c
      return
      end


      subroutine initgrid(nclx,ncly,nclz,gridpar,const)
c-----------------------------------------------------------------------
c     initialize a grid paramter
c
      implicit none
      integer*4 nclx,ncly,nclz
      real*4  gridpar(*)
      real*8  const
c local
      integer*4 i,nc3
c
      nc3=nclx*ncly*nclz
      do i=1,nc3
         gridpar(i)=const
      enddo
c
      return
      end


      subroutine rect_norm(ntpol,xscale,yscale,zscale,
     $           lstpx,lstpy,lstpz,bnorm)
c-----------------------------------------------------------------------
c     calculate the normalization constants of the basis functions
c
      implicit none
      include 'consta.fcm'
      integer*4 ntpol,lstpx(*),lstpy(*),lstpz(*)
      real*8  bnorm(*),xscale,yscale,zscale
c local
      real*8  ilxyz
      integer*4 n,lpol,mpol,npol

      ilxyz=xscale*yscale*zscale/8.0d0      ! inverse lxyz

      do n=1,ntpol
         lpol=lstpx(n)*2.0d0
         mpol=lstpy(n)*2.0d0
         npol=lstpz(n)*2.0d0
         bnorm(n)=sqrt((lpol+1.0d0)*(mpol+1.0d0)*(npol+1.0d0)*ilxyz)
      enddo
c
      return
      end


      subroutine rect_svpol(xnpol,ynpol,znpol,oldntp,qmij,
     $           lstpx,lstpy,lstpz,lstpol)
c-----------------------------------------------------------------------
c     store the basis functions
c     in lstpx, lstpy, and lstpz array for legendre polynomials
c
      implicit none
      integer*4 xnpol,ynpol,znpol,oldntp
      integer*4 lstpx(*),lstpy(*),lstpz(*),lstpol(*)
      logical qmij
c local
      integer*4 xpol,ypol,zpol,norder,n,ntpol
c
      ntpol=xnpol*ynpol*znpol
      if(qmij.and.oldntp.eq.ntpol) then
         do norder=1,ntpol
            lstpol(norder)=norder
         enddo
         return
      endif
c
      norder=0
      if(qmij) then
         do norder=1,oldntp
            lstpol(norder)=norder
         enddo
         norder=oldntp
         do xpol=0,xnpol-1
         do ypol=0,ynpol-1
         do 100 zpol=0,znpol-1
            do n=1,oldntp
               if(lstpx(n).eq.xpol.and.lstpy(n).eq.ypol.and.
     $              lstpz(n).eq.zpol) goto 100
            enddo
            norder=norder+1 
            lstpx(norder)=xpol
            lstpy(norder)=ypol
            lstpz(norder)=zpol
            lstpol(norder)=norder
 100     enddo
         enddo
         enddo
      else
         do xpol=0,xnpol-1
         do ypol=0,ynpol-1
         do zpol=0,znpol-1
            norder=norder+1 
            lstpx(norder)=xpol
            lstpy(norder)=ypol
            lstpz(norder)=zpol
            lstpol(norder)=norder
         enddo
         enddo
         enddo
      endif
c
      return
      end


      real*8 function lpol(npol,coor,scal)
c------------------------------------------------------------------------
c     legnedre polynomials
c
      implicit none
      integer*4   npol,n
      real*8    x,coor,scal
      real*8    lpolnp,lpolnpp
      real*8    x2,x3,x4,x5,x6,x7
c
      x=scal*coor
c
      if(npol.le.4) then
         if(npol.eq.0) then
            lpol = 1.d0
         elseif(npol.eq.1) then
            lpol = x
         elseif(npol.eq.2) then
            lpol = 0.5d0*(3.d0*x*x-1.d0)
         elseif(npol.eq.3) then
            lpol = 0.5d0*(5.d0*x*x*x-3.d0*x)
         elseif(npol.eq.4) then
            x2=x*x
            lpol = 0.125d0*(35.d0*x2*x2-30.d0*x2+3.d0)
         endif
      elseif(npol.gt.9) then
         x2=x*x
         x3=x2*x
         x4=x2*x2
         x5=x3*x2
         x6=x4*x2
         x7=x5*x2
         lpolnpp= (1.d0/128.d0)*(6435.d0*x6*x2-12012.d0*x6+6930.d0*x4-
     $                       1260.d0*x2+35.d0)
         lpolnp = (1.d0/128.d0)*(12155.d0*x7*x2-25740.d0*x7+18018.d0*x5-
     $                       4620.d0*x3+315.d0*x)
         do n=10,npol
            lpol=((2.d0*n-1)*x*lpolnp-(n-1.d0)*lpolnpp)/n
            lpolnpp=lpolnp
            lpolnp=lpol
         enddo
      else
         if(npol.eq.5) then
            x2=x*x
            x3=x2*x
            lpol = 0.125d0*(63.d0*x3*x2-70.d0*x3+15.d0*x)
         elseif(npol.eq.6) then
            x2=x*x
            x4=x2*x2
            lpol = (1.d0/16.d0)*(231.d0*x4*x2-315.d0*x4+105.d0*x2-5.d0)
         elseif(npol.eq.7) then
            x2=x*x
            x3=x2*x
            x5=x3*x2
            lpol = (1.d0/16.d0)*
     $             (429.d0*x5*x2-693.d0*x5+315.d0*x3-35.d0*x)
         elseif(npol.eq.8) then
            x2=x*x
            x4=x2*x2
            x6=x4*x2
            lpol = (1.d0/128.d0)*(6435.d0*x6*x2-12012.d0*x6+6930.d0*x4-
     $                        1260.d0*x2+35.d0)
         elseif(npol.eq.9) then
            x2=x*x
            x3=x2*x
            x5=x3*x2
            x7=x5*x2
            lpol = (1.d0/128.d0)*(12155.d0*x7*x2-25740.d0*x7+
     $                            18018.d0*x5-4620.d0*x3+315.d0*x)
         endif
      endif
c
      return
      end

      subroutine tmpwrigd(unit,nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           xfir0,yfir0,zfir0,xlas0,ylas0,zlas0,
     $           smthng,prftr)
c-----------------------------------------------------------------------

      implicit none

      integer*4 unit,nclx,ncly,nclz
      real*8  dcel,tranx,trany,tranz,xbcen,ybcen,zbcen,prftr 
      real*4  smthng(*)

c local variables
      integer*4 ixs,ixf,iys,iyf,izs,izf,i,j,k,l
      real*8 xfir0,xlas0,yfir0,ylas0,zfir0,zlas0
      real*8 xfir,xlas,yfir,ylas,zfir,zlas,xtem,ytem,ztem

      xfir   = xfir0+tranx-xbcen
      ixs=int(xfir/dcel)+1
      if(ixs.le.0)ixs=1
      if(ixs.gt.nclx)ixs=nclx

      yfir   = yfir0+trany-ybcen
      iys=int(yfir/dcel)+1
      if(iys.le.0)iys=1
      if(iys.gt.ncly)iys=ncly

      zfir   = zfir0+tranz-zbcen
      izs=int(zfir/dcel)+1
      if(izs.le.0)izs=1
      if(izs.gt.nclz)izs=nclz

      xlas   = xlas0+tranx-xbcen
      ixf=int(xlas/dcel)+1
      if(ixf.le.0)ixf=1
      if(ixf.gt.nclx)ixf=nclx

      ylas   = ylas0+trany-ybcen
      iyf=int(ylas/dcel)+1
      if(iyf.le.0)iyf=1
      if(iyf.gt.ncly)iyf=ncly

      zlas   = zlas0+tranz-zbcen
      izf=int(zlas/dcel)+1
      if(izf.le.0)izf=1
      if(izf.gt.nclz)izf=nclz

c      write(6,'(3f10.5)') tranx,trany,tranz
c      write(6,'(6f10.5)') xfir0,yfir0,zfir0,xlas0,ylas0,zlas0
c      write(6,'(6f10.5)') xfir,yfir,zfir,xlas,ylas,zlas
c      write(6,'(6i10)')   ixs,iys,izs,ixf,iyf,izf

      do 10 i=ixs,ixf
         xtem=(i-1)*dcel-tranx+xbcen
      do 10 j=iys,iyf
         ytem=(j-1)*dcel-trany+ybcen
      do 10 k=izs,izf
         l=(i-1)*ncly*nclz+(j-1)*nclz+k
         ztem=(k-1)*dcel-tranz+zbcen
         write(unit,20) xtem,ytem,ztem,smthng(l)*prftr
   20    format(3f10.5,2x,e16.7)
   10 continue

      return
      end
