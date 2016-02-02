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

      subroutine mayer(natom,x,y,z,cg,radius,epsw,epsp,kappa2,temp,
     $           watr,ionr,nclx,ncly,nclz,dcel,
     $           phi,fcden,mcden,epsx,epsy,epsz,
     $           ntype,cion,zion,ctop,cbot,kap2top,kap2bot,
     $           vmemb,tmemb,zmemb,epsm,htmemb,epsh,
     $           epso,xocyl,yocyl,zocyl,rocyl,hocyl,qokap,
     $           epsb,bxmax,bxmin,bymax,bymin,bzmax,bzmin,qbkap,
     $           epsc,xcyln,ycyln,zcyln,rcyln,hcyln,qckap,
     $           epsec,xecyln,yecyln,zecyln,ax,by,tmin,hecyln,qeckap,
     $           epss,xsphe,ysphe,zsphe,rsphe,qskap,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           qreen,mapt,pox,poy,poz,listr,fistr,
     $           qpnp,qphixypbc,qphixyzpbc,qphi,qmcden,qnmcden,
     $           qzerobp,qphifocus,phib,qcionfocus,bcion,
     $           qnonlinear,qpartlinear)
c-----------------------------------------------------------------------
c     construction of all the 3d arrays necessary for performing the
c     iterations on both pb and pnp equations.
c
      implicit none
      include 'consta.fcm'
      integer*4  natom,ntype
      integer*4  nclx,ncly,nclz,mapt,listr(*)
      real*4   phi(*),phib(*),fcden(*),mcden(*),cion(*),bcion(*)
      real*4   epsx(*),epsy(*),epsz(*)
      real*8   x(*),y(*),z(*),cg(*),radius(*),ctop(*),cbot(*),zion(*)
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8   epsw,epsp,kappa2,temp,watr,ionr,kap2top,kap2bot
      real*8   pox(*),poy(*),poz(*),fistr(*)
      real*8   vmemb,zmemb,tmemb,epsm,htmemb,epsh
      real*8   epso,xocyl,yocyl,zocyl,rocyl,hocyl
      real*8   epsb,bxmax,bxmin,bymax,bymin,bzmax,bzmin
      real*8   epsc,xcyln,ycyln,zcyln,rcyln,hcyln
      real*8   epsec,xecyln,yecyln,zecyln,ax,by,tmin,hecyln
      real*8   epss,xsphe,ysphe,zsphe,rsphe
      logical  qpnp,qphixypbc,qphixyzpbc,qzerobp,qphifocus,qcionfocus
      logical  qreen,qphi,qmcden,qckap,qskap,qbkap,qokap,qeckap
      logical  qnonlinear,qpartlinear
      logical  qnmcden

c local
      integer*4 nc3,ncyz
      integer*4 i,j,ip0

      ncyz=ncly*nclz
      nc3=nclx*ncyz


c     initialize the grid parameters
c     ==============================
      do i=1,nc3
         fcden(i)=0d0
         epsx(i)=epsw
         epsy(i)=epsw
         epsz(i)=epsw
      enddo

      if(.not.qphi) then   
         do i=1,nc3
            phi(i)=0d0
         enddo
      endif

      if(.not.qmcden) then
         if(qnmcden) then       ! different mcden for ions
            do i=1,ntype
               do j=1,nc3
                  ip0=(i-1)*nc3+j
                  mcden(ip0)=1d0
               enddo
            enddo
         else                   ! same mcden for ions
            do i=1,nc3
               mcden(i)=1d0
            enddo
         endif
      endif

      if(qpnp) then
         do i=1,ntype
            do j=1,nc3
               ip0=(i-1)*nc3+j
               cion(ip0)=0d0
            enddo
         enddo
      endif


c     membrane part (if requested)
c     ============================

      if(tmemb.gt.0.0d0)then
         call mayer_memb(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        zbcen,tmemb,zmemb,epsm,htmemb,epsh,
     $        mcden,epsx,epsy,epsz,qmcden)
      endif


c     sphere part (if requested)
c     ==========================

      if(rsphe.gt.0.0d0)then
         call mayer_sphe(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epss,xsphe,ysphe,zsphe,rsphe,
     $        mcden,epsx,epsy,epsz,qskap,qmcden)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     kyu il lee added elliptic cylinder
c     ==========================
c
      if(ax.gt.0.0d0.and.by.gt.0.0d0)then
c         write(*,*) 'ecynl called'
         call mayer_ecyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epsec,xecyln,yecyln,zecyln,ax,by,tmin,hecyln,
     $        mcden,epsx,epsy,epsz,qeckap,qmcden)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     cylinder part (if requested)
c     ============================

      if(rcyln.gt.0.0d0)then
         call mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epsc,xcyln,ycyln,zcyln,rcyln,hcyln,
     $        mcden,epsx,epsy,epsz,qckap,qmcden)
      endif


c     box part (if requested)
c     =======================
      if(bxmax.ne.bxmin)then
         call mayer_box(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epsb,bxmax,bymax,bzmax,bxmin,bymin,bzmin,
     $        mcden,epsx,epsy,epsz,qbkap,qmcden)
      endif


c only for pbphix
c      call mayer_box(nclx,ncly,nclz,dcel,tranx,trany,tranz,
c     $     xbcen,ybcen,zbcen,
c     $     epsb,100.0d0,100.0d0,5.5d0,-100.0d0,-100.0d0,-5.5d0,
c     $     mcden,epsx,epsy,epsz,qbkap,qmcden)


c     set the potential at the edge of the box (boundary potentials)
c     =============================================================

      if(.not.qphi.and..not.qzerobp.and..not.qphixyzpbc) then
         if(qphifocus) then
            call mayer_bp_focus(phi,phib,nclx,ncly,nclz)
         else
            call mayer_bp_int(natom,x,y,z,cg,radius,epsw,
     $           kappa2,kap2top,kap2bot,
     $           phi,nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,qphixypbc)
         endif
      endif


c     charge distribution
c     ===================

      call mayer_cd_triln(natom,x,y,z,cg,fcden,
     $     nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $     xbcen,ybcen,zbcen,qphixyzpbc,kappa2)


c     solute volume exclusion functions (mayer function)
c     ==================================================
c     when the dielectric boundary is defined by the van der waals surface
c     mcden = 0 : inside solute
c             1 : outside solute

      call mayer_mstep(natom,x,y,z,radius,epsp,watr,ionr,
     $     nclx,ncly,nclz,dcel,
     $     tranx,trany,tranz,xbcen,ybcen,zbcen,
     $     mcden,epsx,epsy,epsz,qreen,qmcden)


c     probe radius dependent part
c     ===========================
c     works when radw >  0. .and. qreen = true
c     then all grid points inside accessible surface
c     are already set as a solute and it is now
c     necessary to extend the dielectric 
c     till the molecular (contact+reentrance) surface

      if( watr.gt.0.0d0 .and. qreen ) then
         call mayer_reen(natom,x,y,z,radius,epsp,watr,
     $        nclx,ncly,nclz,dcel,
     $        tranx,trany,tranz,xbcen,ybcen,zbcen,
     $        mcden,epsx,epsy,epsz,
     $        pox,poy,poz,mapt,listr,fistr,qmcden)
      endif


c     overlay cylinder part (if requested)
c     ============================

      if(rocyl.gt.0.0d0)then
         call mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        xbcen,ybcen,zbcen,
     $        epso,xocyl,yocyl,zocyl,rocyl,hocyl,
     $        mcden,epsx,epsy,epsz,qokap,qmcden)
      endif


c     set the membrane potential
c     ==========================
c     fix the initial membrane potential and the boltzmann distribution
c     for the salt concentration due to presence of offset potential vmemb

      if(vmemb.ne.0.0d0)then
         call mayer_vmemb(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $        zbcen,epsw,kap2top,vmemb,tmemb,zmemb,epsm,
     $        fcden,phi,mcden,qpnp,qphifocus,qphi,
     $        qnonlinear,qpartlinear)
      endif


c     set boundary concentrations
c     =========================== 
c     mayer_bc_lin should be called first to fill cion array

      if(qpnp) then
         if(qnmcden) then
            call mayer_bc_lin2(ntype,nclx,ncly,nclz,dcel,cion,zion,
     $           phi,mcden,ctop,cbot,vmemb,tmemb,zmemb,tranz,zbcen,temp)
            if(qcionfocus) then
               call mayer_bc_focus2(ntype,mcden,cion,bcion,
     $              nclx,ncly,nclz)
            endif
         else
            call mayer_bc_lin1(ntype,nclx,ncly,nclz,dcel,cion,zion,
     $           phi,mcden,ctop,cbot,vmemb,tmemb,zmemb,tranz,zbcen,temp)
            if(qcionfocus) then
               call mayer_bc_focus1(ntype,mcden,cion,bcion,
     $              nclx,ncly,nclz)
            endif
            do i=1,ntype
               do j=1,nc3
                  ip0=(i-1)*nc3+j
                  mcden(ip0)=mcden(j)
               enddo
            enddo
         endif
      endif

c
      return
      end

      subroutine mayer_vmemb(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           zbcen,epsw,kappa2,vmemb,tmemb,zmemb,epsm,
     $           fcden,phi,mcden,qpnp,qphifocus,qphi,
     $           qnonlinear,qpartlinear)
c------------------------------------------------------------------------
c     membrane potential
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz
      real*4   fcden(*),phi(*),mcden(*)
      real*8   dcel,tranx,trany,tranz,zbcen
      real*8   epsw,kappa2,vmemb,tmemb,zmemb,epsm
      logical  qpnp,qphifocus,qphi,qnonlinear,qpartlinear
c local
      integer*4  ncyz,kg,ig,iii,jg,ipo,iiz,iix,iiix,iiy
      real*8   dcel2,zmemb1,zmemb2,kappa,afact,zc,phic
c
      dcel2=dcel*dcel
      ncyz=ncly*nclz
      kappa=sqrt(kappa2/epsw)
      zmemb1=tranz-zbcen-0.5d0*tmemb+zmemb-rsmall
      zmemb2=zmemb1+tmemb+2.0d0*rsmall

c add the analytical transmembrane potential 
      if(.not.qphifocus.and..not.qphi) then
         afact = vmemb/(2.0d0+(epsw/epsm)*kappa*tmemb)
         do kg=1,nclz           !,ncel-1
            zc = (kg-1)*dcel
            if(zc.lt.zmemb1)then
               phic = afact*exp(kappa*(zc-zmemb1))
            elseif((zc.ge.zmemb1).and.(zc.le.zmemb2))then
               phic = afact*((epsw/epsm)*kappa*(zc-zmemb1)+1.0d0)
            elseif(zc.gt.zmemb2)then
               phic = vmemb-afact*exp(-kappa*(zc-zmemb2))
            endif
            do ig=1,nclx
               iii=(ig-1)*ncyz+kg
               do jg=1,ncly
                  ipo=iii+(jg-1)*nclz
                  phi(ipo)=phi(ipo)+phic
               enddo
            enddo
         enddo
      endif

c construct the background charge density for transmembrane potential 
      if(.not.qpnp.and..not.qnonlinear.and..not.qpartlinear) then
         do iiz=1,nclz
            zc=(iiz-1)*dcel
            if(zc.gt.zmemb2)then
               do iix=1,nclx
                  iiix=(iix-1)*ncyz
                  do iiy=1,ncly
                     iii=iiix+(iiy-1)*nclz+iiz
                     fcden(iii)=fcden(iii)+mcden(iii)*kappa2*dcel2*vmemb
                  enddo
               enddo
            endif
         enddo
      endif
c
      return
      end

      subroutine mayer_mstep(natom,x,y,z,radius,epsp,watr,ionr,
     $           nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           mcden,epsx,epsy,epsz,qreen,qmcden)
c----------------------------------------------------------------------
c     when the dielectric boundary is defined by the van der waals surface
c     mcden = 0 : inside solute
c             1 : outside solute
c
      implicit none
      include 'consta.fcm'
      integer*4  natom
      integer*4  nclx,ncly,nclz
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   x(*),y(*),z(*),radius(*),epsp,watr,ionr
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      logical  qreen,qmcden
c local
      integer*4  i,k,l,m,ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2
      integer*4  ipx,ipy,ipz
      integer*4  nfil,ncyz
      real*8   xi,yi,zi,xc,yc,zc,xc1,yc1,zc1,xsq,ysq,zsq,dsq,dsq1
      real*8   sqr,sqr2,sqrw,sqri
c
      ncyz=ncly*nclz
      do 101 i=1,natom
         xi=x(i)+tranx-xbcen
         yi=y(i)+trany-ybcen
         zi=z(i)+tranz-zbcen
         sqr =radius(i)
         if(sqr.lt.0.0d0) goto 101
         sqrw=sqr+watr
         sqri=sqr+ionr
         nfil=int((sqrw+ionr)/dcel+rsmall)+2
         sqr =sqr *sqr  + rsmall
         sqr2=sqr     - 2*rsmall
         sqrw=sqrw*sqrw + rsmall
         sqri=sqri*sqri + rsmall
         ix=int(xi/dcel+rsmall)+1
         iy=int(yi/dcel+rsmall)+1
         iz=int(zi/dcel+rsmall)+1
c
         jx1=ix-nfil+1
         if(jx1.lt.1)jx1=1
         jx2=ix+nfil-1
         if(jx2.gt.nclx)jx2=nclx
         jy1=iy-nfil+1
         if(jy1.lt.1)jy1=1
         jy2=iy+nfil-1
         if(jy2.gt.ncly)jy2=ncly
         jz1=iz-nfil+1
         if(jz1.lt.1)jz1=1
         jz2=iz+nfil-1
         if(jz2.gt.nclz)jz2=nclz

c
         do k=jx1,jx2
            ipx=(k-1)*ncyz
            xc=(k-1)*dcel
            do l=jy1,jy2
               ipy=(l-1)*nclz
               yc=(l-1)*dcel
               do m=jz1,jz2
                  ipz=m+ipy+ipx
                  zc=(m-1)*dcel
                  xsq=(xc-xi)*(xc-xi)
                  ysq=(yc-yi)*(yc-yi)
                  zsq=(zc-zi)*(zc-zi)

c zero the debye-huckel factor (no ions) inside the solute + stern layer
                  if(.not.qmcden)then  
                  if(mcden(ipz).ne.0d0)then
                     dsq=xsq+ysq+zsq
                     if(qreen)then
                        if(ionr.gt.0.0d0) then
                           if(dsq.le.sqri) mcden(ipz)=0d0
                        else
                           if(dsq.le.sqr)then
                              mcden(ipz)=0d0
                           elseif(dsq.gt.sqr2.and.dsq.le.sqrw)then
                             if(mcden(ipz).gt.0d0)mcden(ipz)=-mcden(ipz)
                           endif
                        endif
                     else
                        if(dsq.le.sqri) mcden(ipz)=0.0d0
                     endif
                  endif
                  endif

c protein dielectric constant (epsp) inside the solute
                  if(epsx(ipz).ne.epsp)then
                     xc1=xc+0.5d0*dcel
                     dsq1=(xc1-xi)**2+ysq+zsq
                     if(qreen)then
                        if(dsq1.le.sqr)then
                           epsx(ipz)=epsp 
                        elseif(dsq1.gt.sqr2.and.dsq1.le.sqrw)then
                           if(epsx(ipz).gt.0d0) epsx(ipz)=-epsx(ipz)
                        endif
                     else
                        if(dsq1.le.sqrw) epsx(ipz)=epsp 
                     endif
                  endif

                  if(epsy(ipz).ne.epsp)then
                     yc1=yc+0.5d0*dcel
                     dsq1=(yc1-yi)**2+xsq+zsq
                     if(qreen)then
                        if(dsq1.le.sqr)then
                           epsy(ipz)=epsp 
                        elseif(dsq1.gt.sqr2.and.dsq1.le.sqrw)then
                           if(epsy(ipz).gt.0d0) epsy(ipz)=-epsy(ipz)
                        endif
                     else
                        if(dsq1.le.sqrw) epsy(ipz)=epsp 
                     endif
                  endif

                  if(epsz(ipz).ne.epsp)then
                     zc1=zc+0.5d0*dcel
                     dsq1=(zc1-zi)**2+ysq+xsq
                     if(qreen)then
                        if(dsq1.le.sqr)then
                           epsz(ipz)=epsp 
                        elseif(dsq1.gt.sqr2.and.dsq1.le.sqrw)then
                           if(epsz(ipz).gt.0d0) epsz(ipz)=-epsz(ipz)
                        endif
                     else
                        if(dsq1.le.sqrw) epsz(ipz)=epsp 
                     endif
                  endif
c
               enddo
            enddo
         enddo
 101  enddo
c
      return
      end

      subroutine mayer_cd_triln(natom,x,y,z,cg,fcden,
     $           nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,qphixyzpbc,kappa2)
c----------------------------------------------------------------------
c     the trilinear interpolation for charge distribution 
c
      implicit none
      include 'consta.fcm'
      integer*4  natom,nclx,ncly,nclz
      real*4   fcden(*)
      real*8   x(*),y(*),z(*),cg(*),dcel
      real*8   tranx,trany,tranz,xbcen,ybcen,zbcen,kappa2
      logical  qphixyzpbc
c local
      integer*4  i,n1,n2,n3,in1,in2,in3
      integer*4  ix,iy,iz,ncyz,nc3
      real*8   chi,xi,yi,zi,ai,bi,ci,fi,qtot
c
      ncyz =ncly*nclz
      do 460 i=1,natom
         chi=cg(i)
         xi=x(i)+tranx-xbcen
         yi=y(i)+trany-ybcen
         zi=z(i)+tranz-zbcen

c  check if a charge is outside the calculation box
         if(xi.lt.0.0d0 .or. xi.gt.2*tranx .or.   
     $      yi.lt.0.0d0 .or. yi.gt.2*trany .or.  
     $      zi.lt.0.0d0 .or. zi.gt.2*tranz) goto 460   

         ix=int(xi/dcel+rsmall)+1
         iy=int(yi/dcel+rsmall)+1
         iz=int(zi/dcel+rsmall)+1
            
c construct the charge distribution by 8 grid points adjacent to the atom
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

                  fcden(in3)=fcden(in3)+(fi*chi)*2.0d0*twopi/dcel

               enddo
            enddo
         enddo
 460  enddo

c construct homogeneous neturalizing background charge if kappa=0
c q(i,j,k) = 4 * pi * qtot / (lx*ly*lz) * dcel2 [e/angs] for finite-difference solver

      if(qphixyzpbc.and.kappa2.eq.0.0d0) then
         qtot=0.0d0
         do i=1,natom
            qtot=qtot+cg(i)
         enddo
         qtot=qtot/(nclx*ncly*nclz*dcel*dcel*dcel)
         nc3=ncyz*nclx
         do i=1,nc3
            fcden(i)=fcden(i)-qtot*2.0d0*twopi*dcel*dcel
         enddo
      endif
c
      return
      end


      subroutine mayer_bp_int(natom,x,y,z,cg,radius,epsw,
     $           kappa2,kap2top,kap2bot,
     $           phi,nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,qphixypbc)
c----------------------------------------------------------------------
c     half number of grid points is used for the debye-huckel approximation
c
      implicit none
      include 'consta.fcm'
      integer*4  natom
      integer*4  nclx,ncly,nclz
      real*4   phi(*)
      real*8   x(*),y(*),z(*),cg(*),radius(*),dcel
      real*8   tranx,trany,tranz,xbcen,ybcen,zbcen,epsw
      real*8   kappa2,kap2top,kap2bot
      logical  qphixypbc
c local
      integer*4  i,ipo,ipox,ig,jg,ncyz,ipoy,kg
      real*8   chi,xi,yi,zi,xj,yj,zj,dist1,dist2
      real*8   kappa,rad,kaptop,kapbot
c
      kappa =sqrt(kappa2/epsw)
      kaptop=sqrt(kap2top/epsw)
      kapbot=sqrt(kap2bot/epsw)
      ncyz=ncly*nclz
      do 445 i=1,natom
         chi=cg(i)
         rad=radius(i)
         if(chi.eq.0.0d0) goto 445
         if(rad.lt.0.0d0) rad=0.0d0
         xi=x(i)+tranx-xbcen
         yi=y(i)+trany-ybcen
         zi=z(i)+tranz-zbcen
         zj=(nclz-1)*dcel-zi


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
                     phi(ipo)=phi(ipo)+
     $                        chi*exp(-kappa*(dist1-rad))/epsw/dist1
                  endif

                  ipo=ipox+ipo
                  dist2=sqrt(xj*xj+yj*yj+zj*zj)
                  if(kappa.eq.0.0d0) then
                     phi(ipo)=phi(ipo)+chi/epsw/dist2
                  else
                     phi(ipo)=phi(ipo)+
     $                        chi*exp(-kappa*(dist2-rad))/epsw/dist2
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
                     phi(ipo)=phi(ipo)+
     $                        chi*exp(-kappa*(dist1-rad))/epsw/dist1
                  endif

                  ipo=ipo+ncyz-nclz
                  dist2=sqrt(xj*xj+yj*yj+zj*zj)
                  if(kappa.eq.0.0d0) then
                     phi(ipo)=phi(ipo)+chi/epsw/dist2
                  else
                     phi(ipo)=phi(ipo)+
     $                        chi*exp(-kappa*(dist2-rad))/epsw/dist2
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
               
               ipo=ipox+(jg-1)*nclz+1             ! z=zmin
               dist1=sqrt(xj*xj+yj*yj+zi*zi)
               if(kapbot.eq.0.0d0) then
                  phi(ipo)=phi(ipo)+chi/epsw/dist1
               else
                  phi(ipo)=phi(ipo)+
     $                     chi*exp(-kapbot*(dist1-rad))/epsw/dist1
               endif

               ipo=ipo+nclz-1                     ! z=zmax
               dist2=sqrt(xj*xj+yj*yj+zj*zj)
               if(kaptop.eq.0.0d0) then
                  phi(ipo)=phi(ipo)+chi/epsw/dist2
               else
                  phi(ipo)=phi(ipo)+
     $                     chi*exp(-kaptop*(dist2-rad))/epsw/dist2
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

      subroutine mayer_memb(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           zbcen,tmemb,zmemb,epsm,htmemb,epsh,
     $           mcden,epsx,epsy,epsz,qmcden)
c----------------------------------------------------------------------
c     membrane part 
c
c              <-tmemb->
c    |       hmemb1 hmemb2
c    |           |   |              |
c    |           v   v              |
c    |         |||||||||            |
c    | phi=0   |||||||||  phi=vmemb |
c    |         |||||||||            |
c    |         ^   ^   ^            |
c    |         |   |   |            |
c    |      zmemb1 |   zmemb2       |
c vmemb1         zmemb              vmemb2
c
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   dcel,tranx,trany,tranz,zbcen
      real*8   tmemb,zmemb,epsm,htmemb,epsh
      logical  qmcden
c local
      integer*4  k,l,m,ipx,ipy,ipz,ncyz
      real*8   zmemb1,zmemb2,hmemb1,hmemb2
      real*8   zc,zc1

      ncyz=ncly*nclz
      zmemb1=tranz-zbcen-0.5d0*tmemb+zmemb-rsmall
      zmemb2=zmemb1+tmemb+2.0d0*rsmall
      hmemb1=zmemb1+htmemb
      hmemb2=zmemb2-htmemb

      do k=1,nclx
         ipx=(k-1)*ncyz
         do l=1,ncly
            ipy=(l-1)*nclz
            do m=1,nclz
               ipz=m+ipy+ipx
               zc=(m-1)*dcel

               if(zc.ge.zmemb1.and.zc.le.zmemb2) then
c     mobile ion concentration zero inside the membrane
                  if(.not.qmcden) mcden(ipz)=0.0d0
c     membrane dielectric constant
                  if(zc.ge.hmemb1.and.zc.le.hmemb2) then
                     epsx(ipz)=epsm
                     epsy(ipz)=epsm
                  else
                     epsx(ipz)=epsh
                     epsy(ipz)=epsh
                  endif
               endif
               zc1=zc+0.5d0*dcel
               if(zc1.ge.zmemb1.and.zc1.le.zmemb2) then
                  if(zc1.ge.hmemb1.and.zc1.le.hmemb2) then
                     epsz(ipz)=epsm
                  else
                     epsz(ipz)=epsh
                  endif
               endif

            enddo
         enddo
      enddo
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c kyu il lee added mayer_ecyln (elliptic cone) for oval-shaped pore
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mayer_ecyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epsec,xecyln,yecyln,zecyln,ax,by,tmin,hecyln,
     $           mcden,epsx,epsy,epsz,qeckap,qmcden)
c----------------------------------------------------------------------
c     elliptic cylinder part 
c
      implicit none
      include 'consta.fcm'
      integer  nclx,ncly,nclz
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8   epsec,xecyln,yecyln,zecyln,ax,by,tmin,hecyln
      logical  qeckap,qmcden
c local
      integer  ncyz,nfil,nfilz,ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2
      integer  k,l,m,ipx,ipy,ipz
      real*8   xi,yi,zi,bisq,zisq,xc,yc,zc,xsq,ysq,zsq,dsq,dsq1,dcel2
      real*8   xc1,yc1,zc1,ax2,by2
      real*8   zzzz

      ncyz=ncly*nclz
      dcel2=dcel*dcel
c
      xi=xecyln+tranx-xbcen
      yi=yecyln+trany-ybcen
      zi=zecyln+tranz-zbcen

      ax2 = ax*ax
      by2 = by*by

      bisq=ax
      if(by.gt.ax) bisq=by

      nfil=int(bisq/dcel+rsmall)+2
c      bisq=bisq*bisq+rsmall
      bisq=1.0d0
      ix=int(xi/dcel+rsmall)+1
      iy=int(yi/dcel+rsmall)+1
      iz=int(zi/dcel+rsmall)+1
      nfilz=int(0.5d0*hecyln/dcel+rsmall)+2
      zisq=hecyln*hecyln/4.0d0+rsmall

c      write(*,*) 'qeckap',qeckap

      jx1=ix-nfil+1
      if(jx1.lt.1)jx1=1
      jx2=ix+nfil
      if(jx2.gt.nclx)jx2=nclx
      jy1=iy-nfil+1
      if(jy1.lt.1)jy1=1
      jy2=iy+nfil
      if(jy2.gt.ncly)jy2=ncly
      jz1=iz-nfilz+1
      if(jz1.lt.1)jz1=1
      jz2=iz+nfilz
      if(jz2.gt.nclz)jz2=nclz
c
      do k=jx1,jx2
         ipx=(k-1)*ncyz
         xc=(k-1)*dcel
         do l=jy1,jy2
            ipy=(l-1)*nclz
            yc=(l-1)*dcel
            do m=jz1,jz2
               ipz=m+ipy+ipx
               zc=(m-1)*dcel
               xsq=(xc-xi)*(xc-xi)/ax2
               ysq=(yc-yi)*(yc-yi)/by2
               zsq=(zc-zi)*(zc-zi)
               zzzz=(zc-zi)*(zc-zi)/zisq * (1.0 - tmin) + tmin
c ion-accessibility inside an elliptic cylinder cavity
               dsq=xsq+ysq
               if(.not.qmcden)then
                  if(dsq.le.zzzz.and.zsq.le.zisq)then
                     mcden(ipz)=0.0
                     if(qeckap) mcden(ipz)=1.0
                  endif
               endif

c dielectric constant of a cylinder cavity 
               if(epsec.ne.0.0) then
                  xc1=xc+0.5*dcel
                  dsq1=(xc1-xi)**2/ax2+ysq
                  if(dsq1.le.zzzz.and.zsq.le.zisq)then
                       epsx(ipz)=epsec
                  endif

                  yc1=yc+0.5*dcel
                  dsq1=(yc1-yi)**2/by2+xsq
                  if(dsq1.le.zzzz.and.zsq.le.zisq)then
                        epsy(ipz)=epsec
                  endif

                  zc1=zc+0.5*dcel
                  dsq1=(zc1-zi)**2
                  if(dsq.le.zzzz.and.dsq1.le.zisq)then
                       epsz(ipz)=epsec
                  endif
               endif

            enddo
         enddo
      enddo
c
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c end of mayer_ecyln
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mayer_cyln(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epsc,xcyln,ycyln,zcyln,rcyln,hcyln,
     $           mcden,epsx,epsy,epsz,qckap,qmcden)
c----------------------------------------------------------------------
c     cylinder part 
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8   epsc,xcyln,ycyln,zcyln,rcyln,hcyln
      logical  qckap,qmcden
c local
      integer*4  ncyz,nfil,nfilz,ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2
      integer*4  k,l,m,ipx,ipy,ipz
      real*8   xi,yi,zi,bisq,zisq,xc,yc,zc,xsq,ysq,zsq,dsq,dsq1,dcel2
      real*8   xc1,yc1,zc1

      ncyz=ncly*nclz
      dcel2=dcel*dcel
c
      xi=xcyln+tranx-xbcen
      yi=ycyln+trany-ybcen 
      zi=zcyln+tranz-zbcen 
      bisq=rcyln
      nfil=int(bisq/dcel+rsmall)+2
      bisq=bisq*bisq+rsmall
      ix=int(xi/dcel+rsmall)+1
      iy=int(yi/dcel+rsmall)+1
      iz=int(zi/dcel+rsmall)+1
      nfilz=int(0.5d0*hcyln/dcel+rsmall)+2
      zisq=hcyln*hcyln/4.0d0+rsmall

      jx1=ix-nfil+1
      if(jx1.lt.1)jx1=1
      jx2=ix+nfil
      if(jx2.gt.nclx)jx2=nclx
      jy1=iy-nfil+1
      if(jy1.lt.1)jy1=1
      jy2=iy+nfil
      if(jy2.gt.ncly)jy2=ncly
      jz1=iz-nfilz+1
      if(jz1.lt.1)jz1=1
      jz2=iz+nfilz
      if(jz2.gt.nclz)jz2=nclz
c
      do k=jx1,jx2
         ipx=(k-1)*ncyz
         xc=(k-1)*dcel
         do l=jy1,jy2
            ipy=(l-1)*nclz
            yc=(l-1)*dcel
            do m=jz1,jz2
               ipz=m+ipy+ipx
               zc=(m-1)*dcel
               xsq=(xc-xi)*(xc-xi)
               ysq=(yc-yi)*(yc-yi)
               zsq=(zc-zi)*(zc-zi)

c ion-accessibility inside a cylinder cavity
               dsq=xsq+ysq
               if(.not.qmcden)then
                  if(dsq.le.bisq.and.zsq.le.zisq)then
                     mcden(ipz)=0.0d0 
                     if(qckap) mcden(ipz)=1.0d0 
                  endif
               endif 
               
c dielectric constant of a cylinder cavity 
               if(epsc.ne.0.0d0) then
                  xc1=xc+0.5d0*dcel
                  dsq1=(xc1-xi)**2+ysq
                  if(dsq1.le.bisq.and.zsq.le.zisq) epsx(ipz)=epsc 
 
                  yc1=yc+0.5d0*dcel
                  dsq1=(yc1-yi)**2+xsq
                  if(dsq1.le.bisq.and.zsq.le.zisq) epsy(ipz)=epsc 
 
                  zc1=zc+0.5d0*dcel
                  dsq1=(zc1-zi)**2
                  if(dsq.le.bisq.and.dsq1.le.zisq) epsz(ipz)=epsc 
               endif

            enddo
         enddo
      enddo
c
      return
      end


      subroutine mayer_sphe(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epss,xsphe,ysphe,zsphe,rsphe,
     $           mcden,epsx,epsy,epsz,qskap,qmcden)
c----------------------------------------------------------------------
c     sphere part 
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8   epss,xsphe,ysphe,zsphe,rsphe
      logical  qskap,qmcden
c local
      integer*4  ncyz,nfil,ix,iy,iz,jx1,jx2,jy1,jy2,jz1,jz2
      integer*4  k,l,m,ipx,ipy,ipz
      real*8   xi,yi,zi,bisq,xc,yc,zc,xsq,ysq,zsq,dsq,dsq1,dcel2
      real*8   xc1,yc1,zc1

      ncyz=ncly*nclz
      dcel2=dcel*dcel
c
      xi=xsphe+tranx-xbcen
      yi=ysphe+trany-ybcen 
      zi=zsphe+tranz-zbcen 
      bisq=rsphe
      nfil=int(bisq/dcel+rsmall)+2
      bisq=bisq*bisq+rsmall
      ix=int(xi/dcel+rsmall)+1
      iy=int(yi/dcel+rsmall)+1
      iz=int(zi/dcel+rsmall)+1

      jx1=ix-nfil+1
      if(jx1.lt.1)jx1=1
      jx2=ix+nfil
      if(jx2.gt.nclx)jx2=nclx
      jy1=iy-nfil+1
      if(jy1.lt.1)jy1=1
      jy2=iy+nfil
      if(jy2.gt.ncly)jy2=ncly
      jz1=iz-nfil+1
      if(jz1.lt.1)jz1=1
      jz2=iz+nfil
      if(jz2.gt.nclz)jz2=nclz
c
      do k=jx1,jx2
         ipx=(k-1)*ncyz
         xc=(k-1)*dcel
         do l=jy1,jy2
            ipy=(l-1)*nclz
            yc=(l-1)*dcel
            do m=jz1,jz2
               ipz=m+ipy+ipx
               zc=(m-1)*dcel
               xsq=(xc-xi)*(xc-xi)
               ysq=(yc-yi)*(yc-yi)
               zsq=(zc-zi)*(zc-zi)

c ion-accessibility inside a spherical cavity
               dsq=xsq+ysq+zsq
               if(.not.qmcden)then
                  if(dsq.le.bisq)then
                     mcden(ipz)=0.0d0 
                     if(qskap) mcden(ipz)=1.0d0 
                  endif
               endif 
               
c dielectric constant of a spherical cavity 
               if(epss.ne.0.0d0)then
                  xc1=xc+0.5d0*dcel
                  dsq1=(xc1-xi)**2+ysq+zsq
                  if(dsq1.le.bisq) epsx(ipz)=epss  

                  yc1=yc+0.5d0*dcel
                  dsq1=(yc1-yi)**2+xsq+zsq
                  if(dsq1.le.bisq) epsy(ipz)=epss 
 
                  zc1=zc+0.5d0*dcel
                  dsq1=(zc1-zi)**2+xsq+ysq
                  if(dsq1.le.bisq) epsz(ipz)=epss 
               endif

            enddo
         enddo
      enddo
c
      return
      end


      subroutine mayer_box(nclx,ncly,nclz,dcel,tranx,trany,tranz,
     $           xbcen,ybcen,zbcen,
     $           epsb,bxmax,bymax,bzmax,bxmin,bymin,bzmin,
     $           mcden,epsx,epsy,epsz,qbkap,qmcden)
c----------------------------------------------------------------------
c     sphere part 
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      real*8   epsb,bxmax,bymax,bzmax,bxmin,bymin,bzmin
      logical  qbkap,qmcden
c local
      integer*4 ncyz,ig,jg,kg,ip0x,ip0y,ip0
      integer*4 jx1,jx2,jy1,jy2,jz1,jz2
      real*8  xc,yc,zc,xc1,yc1,zc1,dcel2
      real*8  bxmax1,bymax1,bzmax1,bxmin1,bymin1,bzmin1
c
      ncyz=ncly*nclz
      dcel2=dcel*dcel
      bxmax1=bxmax+tranx-xbcen+rsmall
      bxmin1=bxmin+tranx-xbcen-rsmall
      bymax1=bymax+trany-ybcen+rsmall
      bymin1=bymin+trany-ybcen-rsmall
      bzmax1=bzmax+tranz-zbcen+rsmall
      bzmin1=bzmin+tranz-zbcen-rsmall
c
      jx1=int(bxmin1/dcel)-1
      jx2=int(bxmax1/dcel)+3
      jy1=int(bymin1/dcel)-1
      jy2=int(bymax1/dcel)+3
      jz1=int(bzmin1/dcel)-1
      jz2=int(bzmax1/dcel)+3

      if(jx1.lt.1)jx1=1
      if(jx2.gt.nclx)jx2=nclx
      if(jy1.lt.1)jy1=1
      if(jy2.gt.ncly)jy2=ncly
      if(jz1.lt.1)jz1=1
      if(jz2.gt.nclz)jz2=nclz

c ion-accessibility inside the box
      if(.not.qmcden) then
         do 101 ig=jx1,jx2
            xc=dcel*(ig-1)
            if(xc.lt.bxmin1.or.xc.gt.bxmax1) goto 101
            ip0x=(ig-1)*ncyz
            do 102 jg=jy1,jy2
               yc=dcel*(jg-1)
               if(yc.lt.bymin1.or.yc.gt.bymax1) goto 102
               ip0y=(jg-1)*nclz+ip0x
               do 103 kg=jz1,jz2
                  zc=dcel*(kg-1)
                  if(zc.lt.bzmin1.or.zc.gt.bzmax1) goto 103
                  ip0 = ip0y + kg

                  mcden(ip0)=0.0d0
                  if(qbkap) mcden(ip0)=1.0d0

 103           enddo
 102        enddo
 101     enddo
      endif


c dielectric constants inside the box
c note: dielectric constants are located at the middle point between grids
      if(epsb.ne.0.0d0) then
         do ig=jx1,jx2
            ip0x=(ig-1)*ncyz
            xc=dcel*(ig-1)
            xc1=xc+dcel*0.5d0
            do jg=jy1,jy2
               ip0y=(jg-1)*nclz+ip0x
               yc=dcel*(jg-1)
               yc1=yc+dcel*0.5d0
               do kg=jz1,jz2
                  ip0 = ip0y + kg
                  zc=dcel*(kg-1)
                  zc1=zc+dcel*0.5d0

                  if(xc1.ge.bxmin1.and.xc1.le.bxmax1.and.
     $                yc.ge.bymin1.and. yc.le.bymax1.and.
     $                zc.ge.bzmin1.and. zc.le.bzmax1     )
     $                 epsx(ip0)=epsb
                  if( xc.ge.bxmin1.and. xc.le.bxmax1.and.
     $               yc1.ge.bymin1.and.yc1.le.bymax1.and.
     $                zc.ge.bzmin1.and. zc.le.bzmax1     )
     $              epsy(ip0)=epsb
                  if( xc.ge.bxmin1.and. xc.le.bxmax1.and.
     $                yc.ge.bymin1.and. yc.le.bymax1.and.
     $               zc1.ge.bzmin1.and.zc1.le.bzmax1     )
     $                 epsz(ip0)=epsb

               enddo
            enddo
         enddo
      endif
c
      return
      end


      subroutine mayer_bp_focus(phi,phib,nclx,ncly,nclz)
c----------------------------------------------------------------------
c     boundary potentials are constructed using previous potentials
c     using focussing method
c
      implicit none
      integer*4  nclx,ncly,nclz
      real*4   phi(*),phib(*)
c local
      integer*4  ipo,ig,jg,kg,cnt,ncyz
c
      ncyz=ncly*nclz
      cnt=0

c yz plane
      do jg=1,ncly
         do kg=1,nclz
            ipo=(jg-1)*nclz+kg                 !yz plane at x=xmin
            cnt=cnt+1
            phi(ipo)=phib(cnt)
            ipo=(nclx-1)*ncyz+(jg-1)*nclz+kg   !yz plane at x=xmax
            cnt=cnt+1
            phi(ipo)=phib(cnt)
         enddo
      enddo
         
c xz plane         
      do ig=2,nclx-1
         do kg=1,nclz
            ipo=(ig-1)*ncyz+kg                 !xz plane at y=ymin
            cnt=cnt+1
            phi(ipo)=phib(cnt)
            ipo=(ig-1)*ncyz+ncyz-nclz+kg       !xz plane at y=ymax
            cnt=cnt+1               
            phi(ipo)=phib(cnt)
         enddo
      enddo

c xy plane
      do ig=2,nclx-1
         do jg=2,ncly-1
            ipo=(ig-1)*ncyz+(jg-1)*nclz+1      !xy plane at z=zmin
            cnt=cnt+1
            phi(ipo)=phib(cnt)
            ipo=(ig-1)*ncyz+(jg-1)*nclz+nclz   !xy plane at z=zmax
            cnt=cnt+1
            phi(ipo)=phib(cnt)
         enddo
      enddo
c
      return
      end


      subroutine mayer_bc_lin1(ntype,nclx,ncly,nclz,dcel,cion,zion,
     $           phi,mcden,ctop,cbot,vmemb,tmemb,zmemb,tranz,zbcen,temp)
c------------------------------------------------------------------------
c     construct initial concentrations of each ion inside the box
c     from their bulk concentrations (ctop and cbot) using
c     a linear function (1) or the boltzmann factor (2).
c
c             (ct-cb)        (cb*l-ct) 
c      c(kg)=-------- * kg + ---------     -- (1)
c              (l-1)           (l-1)
c
c     c(i)=c(kg) * exp[-zion/kt*(phi(i)-vmemb*h)],    -- (2)
c     where h=1 when z>zmemb2 and otherwise h=0
c
      implicit none
      include 'consta.fcm' 
      integer*4 ntype,nclx,ncly,nclz
      real*4  cion(*),phi(*),mcden(*)
      real*8  ctop(*),cbot(*),zion(*)
      real*8  dcel,vmemb,tmemb,zmemb,tranz,zbcen,temp
c local
      integer*4 i,ig,jg,kg,ncyz,nc3,ip0x,ip0,ipi
      real*8  ct,cb,cz1,cz2,factor1,factor3
      real*8  zc,zmembcenter,zmemb2,phif

c

      ncyz=ncly*nclz
      nc3=nclx*ncyz
      factor1=celec /(kboltz*temp/kcalmol)       ! 1/(kcal/(mol*e))->1/(e/a)
      zmembcenter=tranz-zbcen+zmemb-rsmall
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall

c initialize concentrations using the linear function

      do i=1,ntype
         ct=ctop(i)
         cb=cbot(i)
         factor3=zion(i)*factor1
         do kg=1,nclz
            zc=(kg-1)*dcel
            cz1=((ct-cb)*kg+(cb*nclz-ct))/(nclz-1)
            if(zc.lt.zmembcenter) then
               cz2=cb
            elseif(zc.ge.zmembcenter) then
               cz2=ct
            endif
            do ig=1,nclx
               ip0x=(ig-1)*ncyz+kg
               do jg=1,ncly
                  ip0=(jg-1)*nclz+ip0x
                  if(mcden(ip0).ne.0.0d0) then
                     ipi=(i-1)*nc3+ip0
                     phif=phi(ip0)*factor3
                     if(phif.eq.0.0d0) then
                        cion(ipi)=cz1
                     else
                        if(vmemb.ne.0.0d0.and.zc.gt.zmemb2) then
                           phif=(phi(ip0)-vmemb)*factor3
                        endif
                        cion(ipi)=cz2*exp(-phif)
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
c
      return
      end

      subroutine mayer_bc_lin2(ntype,nclx,ncly,nclz,dcel,cion,zion,
     $           phi,mcden,ctop,cbot,vmemb,tmemb,zmemb,tranz,zbcen,temp)
c------------------------------------------------------------------------
c     construct initial concentrations of each ion inside the box
c     from their bulk concentrations (ctop and cbot) using
c     a linear function (1) or the boltzmann factor (2).
c
c             (ct-cb)        (cb*l-ct) 
c      c(kg)=-------- * kg + ---------     -- (1)
c              (l-1)           (l-1)
c
c     c(i)=c(kg) * exp[-zion/kt*(phi(i)-vmemb*h)],    -- (2)
c     where h=1 when z>zmemb2 and otherwise h=0
c
      implicit none
      include 'consta.fcm' 
      integer*4 ntype,nclx,ncly,nclz
      real*4  cion(*),phi(*),mcden(*)
      real*8  ctop(*),cbot(*),zion(*)
      real*8  dcel,vmemb,tmemb,zmemb,tranz,zbcen,temp
c local
      integer*4 i,ig,jg,kg,ncyz,nc3,ip0x,ip0,ipi
      real*8  ct,cb,cz1,cz2,factor1,factor3
      real*8  zc,zmembcenter,zmemb2,phif

c

      ncyz=ncly*nclz
      nc3=nclx*ncyz
      factor1=celec /(kboltz*temp/kcalmol)       ! 1/(kcal/(mol*e))->1/(e/a)
      zmembcenter=tranz-zbcen+zmemb-rsmall
      zmemb2=tranz-zbcen+0.5d0*tmemb+zmemb+rsmall

c initialize concentrations using the linear function

      do i=1,ntype
         ct=ctop(i)
         cb=cbot(i)
         factor3=zion(i)*factor1
         do kg=1,nclz
            zc=(kg-1)*dcel
            cz1=((ct-cb)*kg+(cb*nclz-ct))/(nclz-1)
            if(zc.lt.zmembcenter) then
               cz2=cb
            elseif(zc.ge.zmembcenter) then
               cz2=ct
            endif
            do ig=1,nclx
               ip0x=(ig-1)*ncyz+kg
               do jg=1,ncly
                  ip0=(jg-1)*nclz+ip0x
                  ipi=(i-1)*nc3+ip0
                  if(mcden(ipi).ne.0.0d0) then
                     phif=phi(ip0)*factor3
                     if(phif.eq.0.0d0) then
                        cion(ipi)=cz1
                     else
                        if(vmemb.ne.0.0d0.and.zc.gt.zmemb2) then
                           phif=(phi(ip0)-vmemb)*factor3
                        endif
                        cion(ipi)=cz2*exp(-phif)
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
c
      return
      end


      subroutine mayer_bc_focus1(ntype,mcden,cion,bcion,nclx,ncly,nclz)
c----------------------------------------------------------------------
c     boundary concentrations are constructed using previous concentrations 
c     using focussing method.
c
      implicit none
      integer*4  ntype,nclx,ncly,nclz
      real*4   cion(*),bcion(*),mcden(*)
c local
      integer*4  i,ip0,ipi,ig,jg,kg,cnt,cmt,nc3,ncyz,nbound
c
      ncyz=ncly*nclz
      nc3=nclx*ncyz
      nbound=2*(nclx*ncly+ncly*nclz+nclz*nclx)
      cnt=0

c yz plane
      do jg=1,ncly
         do kg=1,nclz

            cnt=cnt+1
            do i=1,ntype

               ip0=(jg-1)*nclz+kg                 !yz plane at x=xmin
               if(mcden(ip0).ne.0.0d0) then         !check ion exclusion region 
                  ipi=(i-1)*nc3+ip0
                  cmt=(i-1)*nbound+cnt
                  cion(ipi)=bcion(cmt)
               endif

               ip0=(nclx-1)*ncyz+(jg-1)*nclz+kg   !yz plane at x=xmax
               if(mcden(ip0).ne.0.0d0) then         !check ion exclusion region 
                  ipi=(i-1)*nc3+ip0
                  cmt=(i-1)*nbound+cnt+1
                  cion(ipi)=bcion(cmt)
               endif

            enddo
            cnt=cnt+1

         enddo
      enddo
         
c xz plane         
      do ig=2,nclx-1
         do kg=1,nclz

            cnt=cnt+1
            do i=1,ntype

               ip0=(ig-1)*ncyz+kg                 !xz plane at y=ymin
               if(mcden(ip0).ne.0.0d0) then         !check ion exclusion region 
                  ipi=(i-1)*nc3+ip0
                  cmt=(i-1)*nbound+cnt
                  cion(ipi)=bcion(cmt)
               endif

               ip0=(ig-1)*ncyz+ncyz-nclz+kg       !xz plane at y=ymax
               if(mcden(ip0).ne.0.0d0) then         !check ion exclusion region 
                  ipi=(i-1)*nc3+ip0
                  cmt=(i-1)*nbound+cnt+1
                  cion(ipi)=bcion(cmt)
               endif

            enddo
            cnt=cnt+1

         enddo
      enddo

c xy plane
      do ig=2,nclx-1
         do jg=2,ncly-1

            cnt=cnt+1
            do i=1,ntype

               ip0=(ig-1)*ncyz+(jg-1)*nclz+1      !xy plane at z=zmin
               if(mcden(ip0).ne.0.0d0) then         !check ion exclusion region
                  ipi=(i-1)*nc3+ip0
                  cmt=(i-1)*nbound+cnt
                  cion(ipi)=bcion(cmt)
               endif

               ip0=(ig-1)*ncyz+(jg-1)*nclz+nclz   !xy plane at z=zmax
               if(mcden(ip0).ne.0.0d0) then         !check ion exclusion region
                  ipi=(i-1)*nc3+ip0
                  cmt=(i-1)*nbound+cnt+1
                  cion(ipi)=bcion(cmt)
               endif

            enddo
            cnt=cnt+1

         enddo
      enddo
c
      return
      end

      subroutine mayer_bc_focus2(ntype,mcden,cion,bcion,nclx,ncly,nclz)
c----------------------------------------------------------------------
c     boundary concentrations are constructed using previous concentrations 
c     using focussing method.
c
      implicit none
      integer*4  ntype,nclx,ncly,nclz
      real*4   cion(*),bcion(*),mcden(*)
c local
      integer*4  i,ip0,ipi,ig,jg,kg,cnt,cmt,nc3,ncyz,nbound
c
      ncyz=ncly*nclz
      nc3=nclx*ncyz
      nbound=2*(nclx*ncly+ncly*nclz+nclz*nclx)
      cnt=0

c yz plane
      do jg=1,ncly
         do kg=1,nclz

            cnt=cnt+1
            do i=1,ntype

               ip0=(jg-1)*nclz+kg                 !yz plane at x=xmin
               ipi=(i-1)*nc3+ip0
               if(mcden(ipi).ne.0.0d0) then         !check ion exclusion region 
                  cmt=(i-1)*nbound+cnt
                  cion(ipi)=bcion(cmt)
               endif

               ip0=(nclx-1)*ncyz+(jg-1)*nclz+kg   !yz plane at x=xmax
               ipi=(i-1)*nc3+ip0
               if(mcden(ipi).ne.0.0d0) then         !check ion exclusion region 
                  cmt=(i-1)*nbound+cnt+1
                  cion(ipi)=bcion(cmt)
               endif

            enddo
            cnt=cnt+1

         enddo
      enddo
         
c xz plane         
      do ig=2,nclx-1
         do kg=1,nclz

            cnt=cnt+1
            do i=1,ntype

               ip0=(ig-1)*ncyz+kg                 !xz plane at y=ymin
               ipi=(i-1)*nc3+ip0
               if(mcden(ipi).ne.0.0d0) then         !check ion exclusion region 
                  cmt=(i-1)*nbound+cnt
                  cion(ipi)=bcion(cmt)
               endif

               ip0=(ig-1)*ncyz+ncyz-nclz+kg       !xz plane at y=ymax
               ipi=(i-1)*nc3+ip0
               if(mcden(ipi).ne.0.0d0) then         !check ion exclusion region 
                  cmt=(i-1)*nbound+cnt+1
                  cion(ipi)=bcion(cmt)
               endif

            enddo
            cnt=cnt+1

         enddo
      enddo

c xy plane
      do ig=2,nclx-1
         do jg=2,ncly-1

            cnt=cnt+1
            do i=1,ntype

               ip0=(ig-1)*ncyz+(jg-1)*nclz+1      !xy plane at z=zmin
               ipi=(i-1)*nc3+ip0
               if(mcden(ipi).ne.0.0d0) then         !check ion exclusion region
                  cmt=(i-1)*nbound+cnt
                  cion(ipi)=bcion(cmt)
               endif

               ip0=(ig-1)*ncyz+(jg-1)*nclz+nclz   !xy plane at z=zmax
               ipi=(i-1)*nc3+ip0
               if(mcden(ipi).ne.0.0d0) then         !check ion exclusion region
                  cmt=(i-1)*nbound+cnt+1
                  cion(ipi)=bcion(cmt)
               endif

            enddo
            cnt=cnt+1

         enddo
      enddo
c
      return
      end


      subroutine focusphi(nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           fnclx,fncly,fnclz,fdcel,
     $           ftranx,ftrany,ftranz,fxbcen,fybcen,fzbcen,
     $           phi,phib)
c----------------------------------------------------------------------
c     boundary potentials of the focussed system are obtained from 
c     the previous potentials using the trilinear interpolation method.
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz,fnclx,fncly,fnclz
      real*4   phi(*),phib(*)
      real*8   dcel,fdcel,tranx,trany,tranz,ftranx,ftrany,ftranz
      real*8   xbcen,ybcen,zbcen,fxbcen,fybcen,fzbcen
c local
      integer*4  ncyz,cnt,n1,in1,n2,in2,n3,in3
      integer*4  ix,iy,iz,ig,jg,kg
      integer*4  jx,jy,jz,jn1,jn2,jn3
      real*8   xi,xj,yi,yj,zi,zj
      real*8   ai,bi,ci,aj,bj,cj,fi,fj
c
      cnt=0
      ncyz=ncly*nclz

c yz plane
      do jg=1,fncly
         do kg=1,fnclz

            cnt=cnt+1           ! yz plane at x=xmin
            phib(cnt)=0d0     ! ipo=(jg-1)*fnclz+kg 
            cnt=cnt+1           ! yz plane at x=xmax
            phib(cnt)=0d0     ! ipx=(nclx-1)*ncyz+(jg-1)*nclz+kg

c find out the previous mapping points
            xi=0.0d0        - ftranx + fxbcen + tranx - xbcen
            xj=2d0*ftranx     - ftranx + fxbcen + tranx - xbcen 
            yi=(jg-1)*fdcel - ftrany + fybcen + trany - ybcen
            zi=(kg-1)*fdcel - ftranz + fzbcen + tranz - zbcen
            ix=int(xi/dcel+rsmall)+1     
            jx=int(xj/dcel+rsmall)+1     
            iy=int(yi/dcel+rsmall)+1
            iz=int(zi/dcel+rsmall)+1

c previous potentials of 8 grid points around current boundary points 
c are interpolated to boundary points.

            do n1=1,2
               in1=ix+n1-1
               jn1=jx+n1-1
               ai=1.0d0-abs(xi-(in1-1)*dcel)/dcel
               aj=1.0d0-abs(xj-(jn1-1)*dcel)/dcel
               in1=(in1-1)*ncyz
               jn1=(jn1-1)*ncyz
               
               do n2=1,2
                  in2=iy+n2-1 
                  bi=1.0d0-abs(yi-(in2-1)*dcel)/dcel
                  in2=(in2-1)*nclz
                  
                  do n3=1,2
                     in3=iz+n3-1
                     ci=1.0d0-abs(zi-(in3-1)*dcel)/dcel
                     fi=ai*bi*ci
                     fj=aj*bi*ci
                     jn3=jn1+in2+in3
                     in3=in1+in2+in3
                     phib(cnt-1)=phib(cnt-1)+fi*phi(in3)
                     phib(cnt)  =phib(cnt)  +fj*phi(jn3)

                  enddo
               enddo
            enddo

         enddo
      enddo
c xz plane         
      do ig=2,fnclx-1
         do kg=1,fnclz

            cnt=cnt+1           ! xz plane at y=ymin  
            phib(cnt)=0.0d0     ! ipo=(ig-1)*ncyz+kg
            cnt=cnt+1           ! xz plane at y=ymax
            phib(cnt)=0.0d0     ! ipy=(ig-1)*ncyz+ncyz-nclz+kg

c find out the previous mapping points
            xi=(ig-1)*fdcel - ftranx + fxbcen + tranx - xbcen
            yi=0.0d0        - ftrany + fybcen + trany - ybcen
            yj=2d0*ftrany     - ftrany + fybcen + trany - ybcen
            zi=(kg-1)*fdcel - ftranz + fzbcen + tranz - zbcen
            ix=int(xi/dcel+rsmall)+1     
            iy=int(yi/dcel+rsmall)+1
            jy=int(yj/dcel+rsmall)+1     
            iz=int(zi/dcel+rsmall)+1

c previous potentials of 8 grid points around current boundary points 
c are interpolated to boundary points.

            do n1=1,2
               in1=ix+n1-1
               ai=1.0d0-abs(xi-(in1-1)*dcel)/dcel
               in1=(in1-1)*ncyz
               
               do n2=1,2
                  in2=iy+n2-1
                  jn2=jy+n2-1
                  bi=1.0d0-abs(yi-(in2-1)*dcel)/dcel
                  bj=1.0d0-abs(yj-(jn2-1)*dcel)/dcel
                  in2=(in2-1)*nclz
                  jn2=(jn2-1)*nclz
                  
                  do n3=1,2
                     in3=iz+n3-1
                     ci=1.0d0-abs(zi-(in3-1)*dcel)/dcel
                     fi=ai*bi*ci
                     fj=ai*bj*ci
                     jn3=in1+jn2+in3
                     in3=in1+in2+in3
                     
                     phib(cnt-1)=phib(cnt-1)+fi*phi(in3)
                     phib(cnt)  =phib(cnt)  +fj*phi(jn3)
                     
                  enddo
               enddo
            enddo
            
         enddo
      enddo

c xy plane
      do ig=2,fnclx-1
         do jg=2,fncly-1

            cnt=cnt+1                ! xy plane at z=zmin
            phib(cnt)=0.0d0          ! ipo=(ig-1)*ncyz+(jg-1)*nclz+1
            cnt=cnt+1                ! xy plane at z=zmax
            phib(cnt)=0.0d0          ! ipz=(ig-1)*ncyz+(jg-1)*nclz+nclz   

c find out the previous mapping points
            xi=(ig-1)*fdcel - ftranx + fxbcen + tranx - xbcen
            yi=(jg-1)*fdcel - ftrany + fybcen + trany - ybcen
            zi=0.0d0        - ftranz + fzbcen + tranz - zbcen
            zj=2d0*ftranz     - ftranz + fzbcen + tranz - zbcen
            ix=int(xi/dcel+rsmall)+1     
            iy=int(yi/dcel+rsmall)+1
            iz=int(zi/dcel+rsmall)+1
            jz=int(zj/dcel+rsmall)+1     

c previous potentials of 8 grid points around current boundary points 
c are interpolated to boundary points.

            do n1=1,2
               in1=ix+n1-1
               ai=1.0d0-abs(xi-(in1-1)*dcel)/dcel
               in1=(in1-1)*ncyz
               
               do n2=1,2
                  in2=iy+n2-1
                  bi=1.0d0-abs(yi-(in2-1)*dcel)/dcel
                  in2=(in2-1)*nclz
                  
                  do n3=1,2
                     in3=iz+n3-1
                     jn3=jz+n3-1
                     ci=1.0d0-abs(zi-(in3-1)*dcel)/dcel
                     cj=1.0d0-abs(zj-(jn3-1)*dcel)/dcel
                     fi=ai*bi*ci
                     fj=ai*bi*cj
                     jn3=in1+in2+jn3
                     in3=in1+in2+in3

                     phib(cnt-1)=phib(cnt-1)+fi*phi(in3)
                     phib(cnt)  =phib(cnt)  +fj*phi(jn3)

                  enddo
               enddo
            enddo

         enddo
      enddo

c
      return
      end


      subroutine focuscion(nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           fnclx,fncly,fnclz,fdcel,
     $           ftranx,ftrany,ftranz,fxbcen,fybcen,fzbcen,
     $           cion,bcion,ntype)
c----------------------------------------------------------------------
c     boundary concentrations of the focussed system are obtained from 
c     the previous potentials using the trilinear interpolation method.
c
      implicit none
      include 'consta.fcm'
      integer*4  nclx,ncly,nclz,fnclx,fncly,fnclz,ntype
      real*4   cion(*),bcion(*)
      real*8   dcel,fdcel,tranx,trany,tranz,ftranx,ftrany,ftranz
      real*8   xbcen,ybcen,zbcen,fxbcen,fybcen,fzbcen
c local
      integer*4  ncyz,nc3,nbound,cnt,cmt,i,n1,in1,n2,in2,n3,in3,ipi,ipj
      integer*4  ix,iy,iz,ig,jg,kg
      integer*4  jx,jy,jz,jn1,jn2,jn3
      real*8   xi,xj,yi,yj,zi,zj
      real*8   ai,bi,ci,aj,bj,cj,fi,fj
c
      ncyz=ncly*nclz
      nc3=nclx*ncyz
      nbound=2*(fnclx*fncly+fncly*fnclz+fnclz*fnclx)
      cnt=0

c yz plane
      do jg=1,fncly
         do kg=1,fnclz

            cnt=cnt+1           
            do i=1,ntype
               cmt=(i-1)*nbound+cnt   ! yz plane at x=xmin
               bcion(cmt)=0.0d0       ! ip0=(jg-1)*fnclz+kg 
               cmt=(i-1)*nbound+cnt+1 ! yz plane at x=xmax
               bcion(cmt)=0.0d0       ! ip0=(jg-1)*fnclz+kg 
            enddo
            cnt=cnt+1

c find out the previous mapping points
            xi=0.0d0        - ftranx + fxbcen + tranx - xbcen
            xj=2*ftranx     - ftranx + fxbcen + tranx - xbcen 
            yi=(jg-1)*fdcel - ftrany + fybcen + trany - ybcen
            zi=(kg-1)*fdcel - ftranz + fzbcen + tranz - zbcen
            ix=int(xi/dcel+rsmall)+1     
            jx=int(xj/dcel+rsmall)+1     
            iy=int(yi/dcel+rsmall)+1
            iz=int(zi/dcel+rsmall)+1

c previous potentials of 8 grid points around current boundary points 
c are interpolated to boundary points.

            do n1=1,2
               in1=ix+n1-1
               jn1=jx+n1-1
               ai=1.0d0-abs(xi-(in1-1)*dcel)/dcel
               aj=1.0d0-abs(xj-(jn1-1)*dcel)/dcel
               in1=(in1-1)*ncyz
               jn1=(jn1-1)*ncyz
               
               do n2=1,2
                  in2=iy+n2-1 
                  bi=1.0d0-abs(yi-(in2-1)*dcel)/dcel
                  in2=(in2-1)*nclz
                  
                  do n3=1,2
                     in3=iz+n3-1
                     ci=1.0d0-abs(zi-(in3-1)*dcel)/dcel
                     fi=ai*bi*ci
                     fj=aj*bi*ci
                     jn3=jn1+in2+in3
                     in3=in1+in2+in3
                     
                     do i=1,ntype
                        ipi=(i-1)*nc3+in3
                        ipj=(i-1)*nc3+jn3
                        cmt=(i-1)*nbound+cnt
                        bcion(cmt-1)=bcion(cmt-1)+fi*cion(ipi)
                        bcion(cmt)  =bcion(cmt)  +fj*cion(ipj)
                     enddo

                  enddo
               enddo
            enddo

         enddo
      enddo

c xz plane         
      do ig=2,fnclx-1
         do kg=1,fnclz

            cnt=cnt+1
            do i=1,ntype
               cmt=(i-1)*nbound+cnt   ! xz plane at y=ymin  
               bcion(cmt)=0.0d0       ! ip0=(ig-1)*ncyz+kg
               cmt=(i-1)*nbound+cnt+1 ! xz plane at y=ymax
               bcion(cmt)=0.0d0       ! ip0=(ig-1)*ncyz+ncyz-nclz+kg
            enddo
            cnt=cnt+1             

c find out the previous mapping points
            xi=(ig-1)*fdcel - ftranx + fxbcen + tranx - xbcen
            yi=0.0d0        - ftrany + fybcen + trany - ybcen
            yj=2*ftrany     - ftrany + fybcen + trany - ybcen
            zi=(kg-1)*fdcel - ftranz + fzbcen + tranz - zbcen
            ix=int(xi/dcel+rsmall)+1     
            iy=int(yi/dcel+rsmall)+1
            jy=int(yj/dcel+rsmall)+1     
            iz=int(zi/dcel+rsmall)+1

c previous potentials of 8 grid points around current boundary points 
c are interpolated to boundary points.

            do n1=1,2
               in1=ix+n1-1
               ai=1.0d0-abs(xi-(in1-1)*dcel)/dcel
               in1=(in1-1)*ncyz
               
               do n2=1,2
                  in2=iy+n2-1
                  jn2=jy+n2-1
                  bi=1.0d0-abs(yi-(in2-1)*dcel)/dcel
                  bj=1.0d0-abs(yj-(jn2-1)*dcel)/dcel
                  in2=(in2-1)*nclz
                  jn2=(jn2-1)*nclz
                  
                  do n3=1,2
                     in3=iz+n3-1
                     ci=1.0d0-abs(zi-(in3-1)*dcel)/dcel
                     fi=ai*bi*ci
                     fj=ai*bj*ci
                     jn3=in1+jn2+in3
                     in3=in1+in2+in3

                     do i=1,ntype
                        ipi=(i-1)*nc3+in3
                        ipj=(i-1)*nc3+jn3
                        cmt=(i-1)*nbound+cnt
                        bcion(cmt-1)=bcion(cmt-1)+fi*cion(ipi)
                        bcion(cmt)  =bcion(cmt)  +fj*cion(ipj)
                     enddo

                     
                  enddo
               enddo
            enddo
            
         enddo
      enddo

c xy plane
      do ig=2,fnclx-1
         do jg=2,fncly-1

            cnt=cnt+1             
            do i=1,ntype
               cmt=(i-1)*nbound+cnt   ! xy plane at z=zmin
               bcion(cmt)=0.0d0       ! ip0=(ig-1)*ncyz+(jg-1)*nclz+1
               cmt=(i-1)*nbound+cnt+1 ! xy plane at z=zmax
               bcion(cmt)=0.0d0       ! ip0=(ig-1)*ncyz+(jg-1)*nclz+nclz   
            enddo
            cnt=cnt+1             

c find out the previous mapping points
            xi=(ig-1)*fdcel - ftranx + fxbcen + tranx - xbcen
            yi=(jg-1)*fdcel - ftrany + fybcen + trany - ybcen
            zi=0.0d0        - ftranz + fzbcen + tranz - zbcen
            zj=2*ftranz     - ftranz + fzbcen + tranz - zbcen
            ix=int(xi/dcel+rsmall)+1     
            iy=int(yi/dcel+rsmall)+1
            iz=int(zi/dcel+rsmall)+1
            jz=int(zj/dcel+rsmall)+1     

c previous potentials of 8 grid points around current boundary points 
c are interpolated to boundary points.

            do n1=1,2
               in1=ix+n1-1
               ai=1.0d0-abs(xi-(in1-1)*dcel)/dcel
               in1=(in1-1)*ncyz
               
               do n2=1,2
                  in2=iy+n2-1
                  bi=1.0d0-abs(yi-(in2-1)*dcel)/dcel
                  in2=(in2-1)*nclz
                  
                  do n3=1,2
                     in3=iz+n3-1
                     jn3=jz+n3-1
                     ci=1.0d0-abs(zi-(in3-1)*dcel)/dcel
                     cj=1.0d0-abs(zj-(jn3-1)*dcel)/dcel
                     fi=ai*bi*ci
                     fj=ai*bi*cj
                     jn3=in1+in2+jn3
                     in3=in1+in2+in3

                     do i=1,ntype
                        ipi=(i-1)*nc3+in3
                        ipj=(i-1)*nc3+jn3
                        cmt=(i-1)*nbound+cnt
                        bcion(cmt-1)=bcion(cmt-1)+fi*cion(ipi)
                        bcion(cmt)  =bcion(cmt)  +fj*cion(ipj)
                     enddo

                  enddo
               enddo
            enddo

         enddo
      enddo
c
      return
      end


      subroutine mayer_reen(natom,x,y,z,radius,epsp,watr,
     $           nclx,ncly,nclz,dcel,
     $           tranx,trany,tranz,xbcen,ybcen,zbcen,
     $           mcden,epsx,epsy,epsz,
     $           pox,poy,poz,mapt,listr,fistr,qmcden)
c------------------------------------------------------------------------
c     molecular (contact+reentrance) surface
c
      implicit none
      include 'consta.fcm'
      integer*4  natom,listr(*)
      integer*4  nclx,ncly,nclz,mapt
      real*4   mcden(*),epsx(*),epsy(*),epsz(*)
      real*8   x(*),y(*),z(*),radius(*),pox(*),poy(*),poz(*),fistr(*)
      real*8   epsp,watr,dcel,tranx,trany,tranz,xbcen,ybcen,zbcen
      logical  qmcden
c local
      integer*4  iseed,il,jl,i,j,l,m,nl,ml,kl,ll,nfil,ix,iy,iz
      integer*4  jx1,jx2,jy1,jy2,jz1,jz2,k,ipx,ipy,ipz,ncyz
      real*8   wrsf,raneti,ranetj,xi,yi,zi,xj,yj,zj,aij,fistr0
      real*8   xp,yp,zp,bisq,xc,yc,zc,xsq,ysq,zsq,dsq,xc1,yc1,zc1,dsq1
      real*8   segpar
c
      iseed=314159265
      ncyz=ncly*nclz

c generate points on the sphere surface and adjust parameter
      wrsf=0d0
      do i=1,natom
         raneti=radius(i)+watr
         if(wrsf.lt.raneti)wrsf=raneti
      enddo
      wrsf=mapt/wrsf**2 + 0.000001d0
      do i=1,mapt
         call genpoint(pox(i),poy(i),poz(i),iseed)
      enddo

c main loop by atoms
      do i=1,natom
         nl=0
         xi=x(i)
         yi=y(i)
         zi=z(i)
         raneti=radius(i)+watr

c create the list of closest neighbors for the atom from the main loop
c according to their screening ability of the neighbors
c 0<segpar<1
c segpar=0  if the neighbor does not reduce the accessible surface of the atom 
c           and it means that it is not a neighbor at all.
c segpar=1  if the neighbor buries the atom totally.
         do j=1,natom
            if(i.eq.j)goto 61
            xj=x(j)
            yj=y(j)
            zj=z(j)
            ranetj=radius(j)+watr
            aij=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))
            fistr0=segpar(raneti,ranetj,aij)
            if(fistr0.eq.0d0)goto 61
            if(fistr0.eq.1d0)goto 60
            nl=nl+1
            fistr(nl)=fistr0
            listr(nl)=j
 61         continue
         enddo

c sort neighbors according to their ability to 
c reduce an accessible surface of the atom
         if(nl.gt.0)call qcksrt(nl,fistr,listr)
         ml=int(wrsf*raneti**2)
         if(ml.gt.mapt)ml=mapt
c loop over points on the sphere surface
         do jl=1,ml
            xp=xi+raneti*pox(jl)
            yp=yi+raneti*poy(jl)
            zp=zi+raneti*poz(jl)
c reject the points which are not on the accessible surface.
            if(nl.gt.0)then
               do kl=1,nl
                  ll=listr(kl)
                  aij=(xp-x(ll))**2+(yp-y(ll))**2+(zp-z(ll))**2
                  if(aij.lt.(radius(ll)+watr)**2)goto 62
               enddo
            endif

c reset grid points inside a probe sphere around any
c of random accessible points as a dielectric media.
            bisq=watr
            nfil=int(bisq/dcel)+2
            xp=xp+tranx-xbcen
            yp=yp+trany-ybcen
            zp=zp+tranz-zbcen
            bisq=bisq*bisq + rsmall

            ix=int(xp/dcel)+1
            iy=int(yp/dcel)+1
            iz=int(zp/dcel)+1

            jx1=ix-nfil+1
            if(jx1.lt.1)jx1=1
            jx2=ix+nfil
            if(jx2.gt.nclx)jx2=nclx
            jy1=iy-nfil+1
            if(jy1.lt.1)jy1=1
            jy2=iy+nfil
            if(jy2.gt.ncly)jy2=ncly
            jz1=iz-nfil+1
            if(jz1.lt.1)jz1=1
            jz2=iz+nfil
            if(jz2.gt.nclz)jz2=nclz
c
            do k=jx1,jx2
               ipx=(k-1)*ncyz
               xc =(k-1)*dcel
               xsq=(xc-xp)*(xc-xp)
               do l=jy1,jy2
                  ipy=(l-1)*nclz
                  yc =(l-1)*dcel
                  ysq=(yc-yp)*(yc-yp)
                  do m=jz1,jz2
                     ipz=m+ipy+ipx
                     zc=(m-1)*dcel
                     zsq=(zc-zp)*(zc-zp)
                     dsq=xsq+ysq+zsq
                     
c if radi > 0 it never happens
                     if(.not.qmcden) then
                     if(mcden(ipz).lt.0d0)then
                        if(dsq.le.bisq)then
                           mcden(ipz)=-mcden(ipz)   ! bulk kappa restored
                        endif
                     endif
                     endif

                     if(epsx(ipz).lt.0d0)then
                        xc1=xc+0.5d0*dcel
                        dsq1=(xc1-xp)**2+ysq+zsq
                        if(dsq1.le.bisq)then
                           epsx(ipz)=-epsx(ipz) ! bulk dielectric constant restored
                        endif
                     endif

                     if(epsy(ipz).lt.0d0)then
                        yc1=yc+0.5d0*dcel
                        dsq1=(yc1-yp)**2+xsq+zsq
                        if(dsq1.le.bisq)then
                           epsy(ipz)=-epsy(ipz) ! bulk dielectric constant restored
                        endif
                     endif

                     if(epsz(ipz).lt.0d0)then
                        zc1=zc+0.5d0*dcel
                        dsq1=(zc1-zp)**2+ysq+xsq
                        if(dsq1.le.bisq)then
                           epsz(ipz)=-epsz(ipz) ! bulk dielectric constant restored
                        endif
                     endif

                  enddo
               enddo
            enddo
 62         continue
         enddo
 60      continue
      enddo
c
      do k=1,nclx
         ipx=(k-1)*ncyz
         do l=1,ncly
            ipy=(l-1)*nclz
            do m=1,nclz
               ipz=m+ipy+ipx
               if(epsz(ipz).lt.0d0)  epsz(ipz)=epsp
               if(epsy(ipz).lt.0d0)  epsy(ipz)=epsp
               if(epsx(ipz).lt.0d0)  epsx(ipz)=epsp
               if(mcden(ipz).lt.0d0) mcden(ipz)=0d0
            enddo
         enddo
      enddo
c
      return
      end


      subroutine genpoint(xp,yp,zp,iseed)
c------------------------------------------------------------------------
c generate random point on the sphere surface (used in mayer)
c
      implicit none
      real*8 xp,yp,zp,rp,random
      integer*4 iseed
c
   10 continue
      xp=2*random(iseed)-1d0
      yp=2*random(iseed)-1d0
      zp=2*random(iseed)-1d0
      rp=xp**2+yp**2+zp**2
      if(rp.gt.1d0)goto 10
      if(rp.eq.0d0)goto 10
      rp=sqrt(rp)
      xp=xp/rp
      yp=yp/rp
      zp=zp/rp
c
      return
      end


      subroutine qcksrt(n,arr,ibrr)
c------------------------------------------------------------------------
c quick sorting of arr (used in mayer)
c
      implicit none
      integer*4 n,ibrr,m,nstack,istack,jstack,l,ir,ib,j,i,iq
      real*8 arr,fm,fa,fc,fmi,a,fx
      parameter (m=7,nstack=50,fm=7875d0,fa=211d0,fc=1663d0,fmi=1d0/fm)
      dimension arr(n),ibrr(n),istack(nstack)
c
      jstack=0
      l=1
      ir=n
      fx=0d0
 10   if(ir-l.lt.m)then
         do 13 j=l+1,ir
            a=arr(j)
            ib=ibrr(j)
            do 11 i=j-1,1,-1
               if(arr(i).le.a)goto 12
               arr(i+1)=arr(i)
               ibrr(i+1)=ibrr(i)
 11         continue
            i=0
 12         arr(i+1)=a
            ibrr(i+1)=ib
 13      continue
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         i=l
         j=ir
         fx=mod(fx*fa+fc,fm)
         iq=l+(ir-l+1)*(fx*fmi)
         a=arr(iq)
         ib=ibrr(iq)
         arr(iq)=arr(l)
         ibrr(iq)=ibrr(l)
 20      continue
 21      if(j.gt.0)then
            if(a.lt.arr(j))then
               j=j-1
               goto 21
            endif
         endif
         if(j.le.i)then
            arr(i)=a
            ibrr(i)=ib
            goto 30
         endif
         arr(i)=arr(j)
         ibrr(i)=ibrr(j)
         i=i+1
 22      if(i.le.n)then
            if(a.gt.arr(i))then
               i=i+1
               goto 22
            endif
         endif
         if(j.le.i)then
            arr(j)=a
            ibrr(j)=ib
            i=j
            goto 30
         endif
         arr(j)=arr(i)
         ibrr(j)=ibrr(i)
         j=j-1
         goto 20
 30      jstack=jstack+2
         if(jstack.gt.nstack) then
            stop 'nstack must be made larger in qcksrt.'
         endif
         if(ir-i.ge.i-l)then
            istack(jstack)=ir
            istack(jstack-1)=i+1
            ir=i-1
         else
            istack(jstack)=i-1
            istack(jstack-1)=l
            l=i+1
         endif
      endif
      goto 10
c
      end


      real*8 function segpar(raneti,ranetj,aij)
c------------------------------------------------------------------------
c segment part formula (used in mayer)
c
      implicit none
      real*8 raneti,ranetj,aij,v1,v2
c
      v1=raneti+ranetj
      v2=abs(raneti-ranetj)
      if(aij.gt.v1)then
         segpar=0d0
      elseif(aij.lt.v2)then
         if(raneti.gt.ranetj)then
            segpar=0d0
         else
            segpar=1d0
         endif
      else
         segpar=(ranetj**2-(raneti-aij)**2)/(4d0*aij*raneti)
      endif
c
      return
      end


      real*8 function random(iseed)
c------------------------------------------------------------------------
c     random number generator: uniform distribution (0,1)
c     iseed: seed for generator. on the first call this has to
c     have a value in the exclusive range (1, 2147483647)
c     and will be replaced by a new value to be used in
c     following call.
c
c     ref: lewis, p.a.w., goodman, a.s. & miller, j.m. (1969)
c     "pseudo-random number generator for the system/360", ibm
c     systems journal 8, 136.
c
c     this is a "high-quality" machine independent generator.
c     integers are supposed to be 32 bits or more.
c     the same algorithm is used as the basic imsl generator.
c
c     author: lennart nilsson
c
      integer*4 iseed
      real*8 dseed,divis,denom,multip
      data  divis/2147483647.d0/
      data  denom /2147483711.d0/
      data  multip/16807.d0/
c
      if(iseed.le.1) iseed=314159
      dseed=multip*iseed
      dseed=mod(dseed,divis)
      random=dseed/denom
      iseed=dseed
c
      return
      end
