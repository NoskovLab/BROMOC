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



c==============================================================================
      subroutine rfpar(ntype,nclx,ncly,nclz,dcel,
     &           xbcen,ybcen,zbcen,
     &           epsx,epsy,epsz,epsw,tranx,trany,tranz,
     &           tmemb,htmemb,zmemb,epsm,epsh,     
     &           rsphe,xsphe,ysphe,zsphe,
     &           bxmax,bymax,bzmax,bxmin,bymin,bzmin,rmax,rion,Qbox,
     &           gnowat,gwater,gsrfen,greff)
c==============================================================================
c
c     2006
c
c     This subroutine prepares the grids for the reaction field parameters
c     and calculates them.
c
c     Input:   ntype (number of ion types)
c              nclx, ncly, nclz (number of grid points in x, y, and z 
c                 direction)
c              dcel (grid spacing)
c              xbcen, ybcen, zbcen (position of grid center)
c              epsx, epsy, epsz (dielectric constants between grid points)
c              epsw (bulk solvent dielectric constant)
c              tranx, trany, tranz (half of the length covered by the
c                 grid in x, y, and z direction, e.g. 
c                 tranx=0.5D0*(nclx-1)*dcel)
c              tmemb (total thickness of membrane along z-axis)
c              htmemb (thickness of headgroup region)
c              zmemb (membrane position along z-axis)
c              epsm (membrane dielectric constant)
c              epsh (membrane headgroup dielectric constant)
c              rsphe (radius of a spherical simulation system)
c              xsphe, ysphe, zsphe (position of the center of a spherical 
c                 simulation system)
c              bxmax, bymax, bzmax, bxmin, bymin, bzmin (boundaries of a
c                 box shaped simulation system)
c              rmax (maximum radius for numerical integration)
c              rion (Born radii of the ion types)
c              Qbox (true means box shaped simulation system,
c                 false means spherical simulation system)              
c     Outpout: gnowat (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are different from epsw)
c              gwater (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are equal to epsw)
c              gsrfen (grids with self reaction field energy parameters)
c              greff (grids with effective radius parameters)
c     Calls:   rfprep, rfcalc              
c              

      implicit none
      include 'consta.fcm'
 
      integer*4 ntype
      integer*4 nclx, ncly, nclz
      real*8 dcel, xbcen, ybcen, zbcen
      real*4 epsx(*), epsy(*), epsz(*) 
      real*8 epsw
      real*8 tranx, trany, tranz
      real*8 tmemb, htmemb, zmemb, epsm, epsh     
      real*8 rsphe, xsphe, ysphe, zsphe
      real*8 bxmax, bymax, bzmax, bxmin, bymin, bzmin 
      real*8 rmax, rion(*) 
      logical Qbox 

      integer*4 gnowat(*), gwater(*)      
      real*4 gsrfen(*), greff(*)
            
      integer*4 jxmin, jxmax, jymin, jymax, jzmin, jzmax      
      
c     Prepare the grids.
      call rfprep(ntype,nclx,ncly,nclz,dcel,
     &     xbcen,ybcen,zbcen,
     &     epsx,epsy,epsz,epsw,tranx,trany,tranz,
     &     rsphe,xsphe,ysphe,zsphe,
     &     bxmax,bymax,bzmax,bxmin,bymin,bzmin,Qbox,
     &     jxmin,jxmax,jymin,jymax,jzmin,jzmax,
     &     gnowat,gwater,gsrfen,greff)                          
                            
c     Calculate reaction field parameters.
      call rfcalc(ntype,
     &     jxmin,jxmax,jymin,jymax,jzmin,jzmax,gnowat,
     &     nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen,
     &     epsx,epsy,epsz,epsw,tranx,trany,tranz,
     &     tmemb,htmemb,zmemb,epsm,epsh,     
     &     rmax,rion,
     &     gsrfen,greff)   
      
      end
      
      
c==============================================================================
      subroutine rfprep(ntype,nclx,ncly,nclz,
     &           dcel,xbcen,ybcen,zbcen,
     &           epsx,epsy,epsz,epsw,tranx,trany,tranz,
     &           rsphe,xsphe,ysphe,zsphe,
     &           bxmax,bymax,bzmax,bxmin,bymin,bzmin,Qbox,
     &           jxmin,jxmax,jymin,jymax,jzmin,jzmax,
     &           gnowat,gwater,gsrfen,greff)
c==============================================================================
c
c     2006
c
c     This subroutine prepares the grids for the reaction field parameters.
c
c     Input:   ntype (number of ion types)
c              nclx, ncly, nclz (number of grid points in x, y, and z 
c                 direction)
c              dcel (grid spacing)
c              xbcen, ybcen, zbcen (position of grid center)
c              epsx, epsy, epsz (dielectric constants between grid points)
c              epsw (bulk solvent dielectric constant)
c              tranx, trany, tranz (half of the length covered by the
c                 grid in x, y, and z direction, e.g. 
c                 tranx=0.5D0*(nclx-1)*dcel)
c              rsphe (radius of a spherical simulation system)
c              xsphe, ysphe, zsphe (position of the center of a spherical 
c                 simulation system)
c              bxmax, bymax, bzmax, bxmin, bymin, bzmin (boundaries of a
c                 box shaped simulation system)
c              Qbox (true means box shaped simulation system,
c                 false means spherical simulation system)              
c     Outpout: jxmin, jxmax, jymin, jymax, jzmin, jzmax (boundary
c                 of a box within the grid which contains only the 
c                 simulation system)
c              gnowat (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are different from epsw)
c              gwater (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are equal to epsw)
c              gsrfen (grids with self reaction field energy parameters,        
c                 gridpoints with negative values are accessible to ions)
c              greff (grids with effective radius parameters)
c     Calls:             
c              

      implicit none
      include 'consta.fcm'
 
      integer*4 ntype, nclx, ncly, nclz
      real*8 dcel, xbcen, ybcen, zbcen
      real*4 epsx(*), epsy(*), epsz(*)
      real*8 epsw
      real*8 tranx, trany, tranz
      real*8 rsphe, xsphe, ysphe, zsphe
      real*8 bxmax, bymax, bzmax, bxmin, bymin, bzmin 
      logical Qbox 
      
      integer*4 jxmin, jxmax, jymin, jymax, jzmin, jzmax
      integer*4 gnowat(*), gwater(*)      
      real*4 gsrfen(*), greff(*)
      
      integer*4 j, n, jx, jy, jz
      integer*4 kx, ky, kz, kaux
      integer*4 kxmin, kxmax, kymin, kymax, kzmin, kzmax
      integer*4 gposx, gposy, gpos
      integer*4 ncyz, nc3, naux
      real*8 xj, yj, zj
      real*8 bxmax1, bxmin1, bymax1, bymin1, bzmax1, bzmin1
      real*8 xsphe1, ysphe1, zsphe1, dsqx, dsqy, dsq
      real*8 rsphsq, maxdif 
      parameter (maxdif=1.0D-8)          
      real*4 srfout, srfin1, srfin2, refout
      parameter (srfout=4.0E0, srfin1=-0.5E0, srfin2=-1.0E0)
      parameter (refout=1.0E-10)
            
      ncyz=ncly*nclz
      nc3=nclx*ncyz
      
      do j=1,nc3
         gwater(j)=0
         gnowat(j)=0
      end do   
      
c     Calculate the grid gwater with the number of next neighbour 
c     epsx, epsy, and epsz values which are equal to epsw and
c     the grid gnowat with the number of next neighbour 
c     epsx, epsy, and epsz values which are different from epsw.
c     Only epsx, epsy, and epsz values within the grid are considered. 
      do jx=1,nclx
         gposx=(jx-1)*ncyz
         do jy=1,ncly
            gposy=gposx+(jy-1)*nclz
            do jz=1,nclz
               gpos=gposy+jz
               if (jx.lt.nclx) then
                  if (dabs(epsx(gpos)-epsw).lt.maxdif) then
                     gwater(gpos)=gwater(gpos)+1
                  else
                     gnowat(gpos)=gnowat(gpos)+1
                  end if
               end if
               if (jy.lt.ncly) then
                  if (dabs(epsy(gpos)-epsw).lt.maxdif) then
                     gwater(gpos)=gwater(gpos)+1
                  else
                     gnowat(gpos)=gnowat(gpos)+1
                  end if
               end if
               if (jz.lt.nclz) then
                  if (dabs(epsz(gpos)-epsw).lt.maxdif) then
                     gwater(gpos)=gwater(gpos)+1
                  else
                     gnowat(gpos)=gnowat(gpos)+1
                  end if
               end if
               if (jx.gt.1) then
                  if (dabs(epsx(gpos-ncyz)-epsw).lt.maxdif) then
                     gwater(gpos)=gwater(gpos)+1
                  else
                     gnowat(gpos)=gnowat(gpos)+1
                  end if  
               end if   
               if (jy.gt.1) then
                  if (dabs(epsy(gpos-nclz)-epsw).lt.maxdif) then
                     gwater(gpos)=gwater(gpos)+1
                  else
                     gnowat(gpos)=gnowat(gpos)+1
                  end if
               end if  
               if (jz.gt.1) then
                  if (dabs(epsz(gpos-1)-epsw).lt.maxdif) then
                     gwater(gpos)=gwater(gpos)+1
                  else
                     gnowat(gpos)=gnowat(gpos)+1
                  end if
               end if  
            end do 
         end do 
      end do   
      
c     Assign the (positive) value srfout to each grid point within 
c     the grid of ion type 1.            
      do j=1,nc3
         gsrfen(j)=srfout 
      end do      
      
      if (Qbox) then           
         bxmax1=bxmax+tranx-xbcen+rsmall
         bxmin1=bxmin+tranx-xbcen-rsmall
         bymax1=bymax+trany-ybcen+rsmall
         bymin1=bymin+trany-ybcen-rsmall
         bzmax1=bzmax+tranz-zbcen+rsmall
         bzmin1=bzmin+tranz-zbcen-rsmall 
      else      
         xsphe1=xsphe+tranx-xbcen
         ysphe1=ysphe+trany-ybcen
         zsphe1=zsphe+tranz-zbcen 
         bxmax1=xsphe1+rsphe+rsmall
         bxmin1=xsphe1-rsphe-rsmall
         bymax1=ysphe1+rsphe+rsmall
         bymin1=ysphe1-rsphe-rsmall
         bzmax1=zsphe1+rsphe+rsmall
         bzmin1=zsphe1-rsphe-rsmall               
         rsphsq=rsphe*rsphe+rsmall            
      end if  

c     Define a box within the grid which contains only the 
c     simulation system.              
      jxmin=INT(bxmin1/dcel)-1
      jxmax=INT(bxmax1/dcel)+3
      jymin=INT(bymin1/dcel)-1
      jymax=INT(bymax1/dcel)+3
      jzmin=INT(bzmin1/dcel)-1
      jzmax=INT(bzmax1/dcel)+3 
         
      if (jxmin.lt.1)    jxmin=1
      if (jxmax.gt.nclx) jxmax=nclx
      if (jymin.lt.1)    jymin=1
      if (jymax.gt.ncly) jymax=ncly
      if (jzmin.lt.1)    jzmin=1
      if (jzmax.gt.nclz) jzmax=nclz
         
c     Assign the (negative) value srfin1 to each grid point within the 
c     grid of ion type 1 which is in contact with the solvent.  
      if (Qbox) then      
c     Grid preparation for a box shaped simulation system.      
         do jx=jxmin,jxmax
            xj=(jx-1)*dcel
            if (xj.ge.bxmin1.and.xj.le.bxmax1) then
               gposx=(jx-1)*ncyz
               do jy=jymin,jymax
                  yj=(jy-1)*dcel
                  if (yj.ge.bymin1.and.yj.le.bymax1) then
                     gposy=gposx+(jy-1)*nclz
                     do jz=jzmin,jzmax
                        zj=(jz-1)*dcel
                        if (zj.ge.bzmin1.and.zj.le.bzmax1) then
                           gpos=gposy+jz
                           if (gwater(gpos).gt.0) then
                              gsrfen(gpos)=srfin1
                           end if
                        end if
                     end do 
                  end if    
               end do
            end if   
         end do         
      else      
c     Grid preparation for a spherical simulation system.      
         do jx=jxmin,jxmax
            xj=(jx-1)*dcel
            gposx=(jx-1)*ncyz
            dsqx=(xj-xsphe1)*(xj-xsphe1)
            do jy=jymin,jymax
               yj=(jy-1)*dcel
               gposy=gposx+(jy-1)*nclz
               dsqy=dsqx+(yj-ysphe1)*(yj-ysphe1)
               do jz=jzmin,jzmax
                  zj=(jz-1)*dcel
                  dsq=dsqy+(zj-zsphe1)*(zj-zsphe1)
                  if (dsq.le.rsphsq) then                  
                     gpos=gposy+jz
                     if (gwater(gpos).gt.0) then
                        gsrfen(gpos)=srfin1
                     end if
                  end if
               end do 
            end do
         end do                 
      end if    
      
c     Extend the number of grid poits where reaction field parameters 
c     shall be calculated. Assign the (negative) value srfin2 to each 
c     grid point which has the value srfout and a neighbour with the 
c     value srfin1.           
      do jx=jxmin,jxmax
         gposx=(jx-1)*ncyz
         if (jx.eq.1) then
            kxmin=jx
         else
            kxmin=jx-1
         end if
         if (jx.eq.nclx) then
            kxmax=jx
         else
            kxmax=jx+1
         end if
         do jy=jymin,jymax
            gposy=gposx+(jy-1)*nclz
            if (jy.eq.1) then
               kymin=jy
            else
               kymin=jy-1
            end if
            if (jy.eq.ncly) then
               kymax=jy
            else
               kymax=jy+1
            end if            
            do jz=jzmin,jzmax
               gpos=gposy+jz
               if (gsrfen(gpos).eq.srfout) then
                  if (jz.eq.1) then
                     kzmin=jz
                  else
                     kzmin=jz-1
                  end if
                  if (jz.eq.nclz) then
                     kzmax=jz
                  else
                     kzmax=jz+1
                  end if
                  do kx=kxmin,kxmax  
                     do ky=kymin,kymax  
                        do kz=kzmin,kzmax 
                           if (.NOT.(kx.eq.jx.and.ky.eq.jy.AND.
     &                        kz.eq.jz)) then
                              kaux=(kx-1)*ncyz+(ky-1)*nclz+kz
                              if (gsrfen(kaux).eq.srfin1) then
                                 gsrfen(gpos)=srfin2
                              end if  
                           end if                                   
                        end do
                     end do
                  end do
               end if 
            end do 
         end do
      end do      
      
c     Set grid values for all ion types. 
      do n=1,ntype 
         naux=(n-1)*nc3 
         do j=1,nc3
            if (n.gt.1) then
c     Copy the gsrfen grid of ion type 1 in order to create the grids for 
c     the other ion types. 
               gsrfen(j+naux)=gsrfen(j)
            end if
            greff(j+naux)=refout
         end do
      end do
                
      end
      
      
c==============================================================================
      subroutine rfcalc(ntype,
     &           jxmin,jxmax,jymin,jymax,jzmin,jzmax,gnowat,
     &           nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen,
     &           epsx,epsy,epsz,epsw,tranx,trany,tranz,
     &           tmemb,htmemb,zmemb,epsm,epsh,     
     &           rmax,rion,
     &           gsrfen,greff)
c==============================================================================
c
c     2006
c
c     This subroutine calculates reaction field parameters for the ion 
c     types.  The results are stored on grids.  The reaction field parameters 
c     of a grid point are calculated only if the corresponding input value
c     of the grid gsrfen is negative.
c
c     Input:   ntype (number of ion types)
c              jxmin, jxmax, jymin, jymax, jzmin, jzmax (boundary
c                 of a box within the grid which contains only the 
c                 simulation system)
c              gnowat (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are different from epsw)
c              nclx, ncly, nclz (number of grid points in x, y, and z 
c                 direction)
c              dcel (grid spacing)
c              xbcen, ybcen, zbcen (position of grid center)
c              epsx, epsy, epsz (dielectric constants between grid points)
c              epsw (bulk solvent dielectric constant)
c              tranx, trany, tranz (half of the length covered by the
c                 grid in x, y, and z direction, e.g. 
c                 tranx=0.5D0*(nclx-1)*dcel)
c              tmemb (total thickness of membrane along z-axis)
c              htmemb (thickness of headgroup region)
c              zmemb (membrane position along z-axis)
c              epsm (membrane dielectric constant)
c              epsh (membrane headgroup dielectric constant)
c              rmax (maximum radius for numerical integration)
c              rion (Born radii of the ion types)
c     Input and Outpout: 
c              gsrfen (grids with self reaction field energy parameters)
c              greff (grids with effective radius parameters)
c     Calls:   intradsph
c              rfint          
c              

      implicit none
      include 'consta.fcm'
      include 'rfintsph.fcm'
 
      integer*4 ntype
      integer*4 jxmin, jxmax, jymin, jymax, jzmin, jzmax
      integer*4 gnowat(*)
      integer*4 nclx, ncly, nclz
      real*8 dcel, xbcen, ybcen, zbcen
      real*4 epsx(*), epsy(*), epsz(*) 
      real*8 epsw
      real*8 tranx, trany, tranz
      real*8 tmemb, htmemb, zmemb, epsm, epsh     
      real*8 rmax, rion(*) 
      
      real*4 gsrfen(*), greff(*) 
      
      integer*4 n, m, j, jx, jy, jz
      integer*4 ncyz, nc3, naux, nsort(ntype)
      integer*4 gposx, gposy, gpos, nrmin(ntype)
c     Order of spherical integration scheme (see rfintsph.fcm). 
      integer*4 nrad, ordsph, nsph
c      parameter (ordsph=17, nsph=n17)
c      parameter (ordsph=29, nsph=n29)
      parameter (ordsph=35, nsph=n35)
c      parameter (ordsph=41, nsph=n41)
      real*8 xj, yj, zj
      real*8 srfen(ntype), reff(ntype)
      real*8 vrad(200), wrad(199) 
      real*8 xsph(nsph), ysph(nsph), zsph(nsph), wsph(nsph)
            
      ncyz=ncly*nclz
      nc3=nclx*ncyz
      
      do n=1,ntype
         nsort(n)=n
      end do
      
c     Sort the ion type numbers into a list nsort(n), so that the ion 
c     radii increase with growing n.       
      do n=1,ntype-1
         do m=n+1,ntype
            If (rion(nsort(n)).gt.rion(nsort(m))) then 
               naux=nsort(n)
               nsort(n)=nsort(m)
               nsort(m)=naux
            end if  
         end do
      end do 
      
c     Calculate number of radial integration points.      
      nrad=INT(rmax)+16
      if (nrad.gt.200) nrad=200
      
c     Generate grid for numerical integration.      
      call intradsph(ntype,nsort,nrad,rmax,rion,ordsph,nsph,
     &   nrmin,vrad,wrad,xsph,ysph,zsph,wsph)
         
c     Calculate reaction field parameters at the grid points with a 
c     negative value. 
      do jx=jxmin,jxmax
         xj=(jx-1)*dcel
         gposx=(jx-1)*ncyz
         do jy=jymin,jymax
            yj=(jy-1)*dcel
            gposy=gposx+(jy-1)*nclz
            do jz=jzmin,jzmax
               zj=(jz-1)*dcel
               gpos=gposy+jz
               if (gsrfen(gpos).lt.0.0E0) then  
c     Calculate reaction field parameters by numerical integration. 
                  call rfint(ntype,nsort,nrad,nsph,
     &               vrad,wrad,xsph,ysph,zsph,wsph,
     &               nrmin,gnowat,
     &               nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen,
     &               epsx,epsy,epsz,epsw,tranz,
     &               tmemb,htmemb,zmemb,epsm,epsh, 
     &               xj,yj,zj,rmax,rion,srfen,reff)
                  do n=1, ntype 
                     naux=(n-1)*nc3     
                     gsrfen(gpos+naux)=SQRT(SNGL(srfen(n)))
                     greff(gpos+naux)=SQRT(SNGL(reff(n)))
                  end do
               end if
            end do 
         end do
      end do

      end
      
 
c==============================================================================
      subroutine rfint(ntype,nsort,nrad,nsph,
     &           vrad,wrad,xsph,ysph,zsph,wsph,
     &           nrmin,gnowat,
     &           nclx,ncly,nclz,dcel,xbcen,ybcen,zbcen,
     &           epsx,epsy,epsz,epsw,tranz,
     &           tmemb,htmemb,zmemb,epsm,epsh, 
     &           xion,yion,zion,
     &           rmax,rion,
     &           srfen,reff)   
c==============================================================================
c
c     2006
c     
c     This subroutine computes itegrals, which are required for 
c     calculating the reaction field parameters of an ion, by 
c     numerical integration using spherical coordinates. 
c  
c     Input:   ntype (number of ion types)
c              nsort (sorted list of ion type numbers with growing Born radii)
c              nrad (number of radial integration points)
c              nsph (number of spherical integration points)
c              vrad (values of radial integration points) 
c              wrad (weights of radial integration points)
c              xsph (x-components of spherical integration points)
c              ysph (y-components of spherical integration points)
c              zsph (z-components of spherical integration points)
c              wsph (weights of spherical integration points)
c              nrmin (grid numbers of the first integration points
c                 for the ion types) 
c              gnowat (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are different from epsw)
c              nclx, ncly, nclz (number of grid points in x, y, and z 
c                 direction)
c              dcel (grid spacing)
c              xbcen, ybcen, zbcen (position of grid center)
c              epsx, epsy, epsz (dielectric constants between grid points)
c              epsw (bulk solvent dielectric constant)
c              tranz (half of the length covered by the grid in z direction, 
c                 tranz=0.5D0*(nclz-1)*dcel)
c              tmemb (total thickness of membrane along z-axis)
c              htmemb (thickness of headgroup region)
c              zmemb (membrane position along z-axis)
c              epsm (membrane dielectric constant)
c              epsh (membrane headgroup dielectric constant)
c              xion, yion, zion (position of the ion)
c              rmax (maximum radius for numerical integration)
c              rion (Born radii of the ions)
c     Outpout: srfen (self reaction field energy parameters)
c              reff (effective radius parameters)
c     Calls:   epssphv 
c
 
      implicit none
      include 'consta.fcm'
 
      integer*4 ntype, nsort(ntype), nrad, nsph
      real*8 vrad(*), wrad(*)
      real*8 xsph(nsph), ysph(nsph), zsph(nsph), wsph(nsph)
      integer*4 nrmin(*), gnowat(*)
      integer*4 nclx, ncly, nclz
      real*8 dcel, xbcen, ybcen, zbcen
      real*4 epsx(*), epsy(*), epsz(*) 
      real*8 epsw
      real*8 tranx, trany, tranz
      real*8 tmemb, htmemb, zmemb, epsm, epsh     
      real*8 xion, yion, zion
      real*8 rmax, rion(*) 
      
      real*8 srfen(ntype), reff(ntype)
       
      logical sphcal(nrad*nsph), case2(ntype)
      integer*4 i, j, jmin, jvol, jmax
      integer*4 m, n, ncyz, ngird, naux
      real*8 epssph(nrad*nsph), rrat(ntype)     
      real*8 sphv(nrad), sphvg1(nrad), sphvg2(nrad)
      real*8 sphvx(nrad), sphvy(nrad), sphvz(nrad) 
      real*8 sphvm(nrad), omsphv(nrad)
      real*8 vol(nrad), volx(nrad), voly(nrad) 
      real*8 volz(nrad), volm(nrad)
      real*8 vola(nrad), volxa(nrad), volya(nrad) 
      real*8 volza(nrad), volma(nrad)
      real*8 intv0(ntype)
      real*8 intv1(nrad), intv1a(nrad)
      real*8 intvg1(nrad)
      real*8 intvg2(nrad), dirg2(nsph), sphg2
      real*8 vrad3(nrad)
      real*8 w, w1, w2
      real*8 rratx, rraty, rratz
      real*8 maxx, maxy, maxz
      real*8 gam1, gam2
      real*8 epsaux
      real*8 wrmin
      real*8 cmemb, zmemb1, zmemb2, hmemb1, hmemb2
      real*8 sphaux 
      real*8 volfac, volaux, volmax
      real*8 volref
      parameter (volref=5000.0D0)
      real*8 maxdif
      parameter (maxdif=1.0D-8)
      real*8 exp1, exp2, exp3
      parameter (exp1=2.8D0, exp2=1.26D0, exp3=0.0D0)
      
      volfac=4.0D0*pi/3.0D0
            
      ncyz=ncly*nclz
      
      maxx=DBle(nclx-1)*dcel
      maxy=DBle(ncly-1)*dcel
      maxz=DBle(nclz-1)*dcel
      
      cmemb=tranz-zbcen+zmemb                  
      zmemb1=cmemb-0.5D0*tmemb
      zmemb2=zmemb1+tmemb
      hmemb1=zmemb1+htmemb
      hmemb2=zmemb2-htmemb

      do 200 n=ntype, 1, -1
         m=nsort(n)
         
         if (n.eq.ntype) then
c     Biggest ion.
            jmax=nrad
            do j=nrmin(m), nrad
c     Spherical integration.        
               call epssphv(j,nsph,xsph,ysph,zsph,wsph,vrad(j),
     &              xion,yion,zion,
     &              epsw,epsm,epsh,epsx,epsy,epsz,
     &              maxx,maxy,maxz,zmemb1,zmemb2,hmemb1,hmemb2,
     &              exp2,exp3,gnowat,
     &              dcel,nclz,ncyz,sphcal,epssph,
     &              sphv(j),sphvx(j),sphvy(j),sphvz(j),sphvm(j),
     &              omsphv(j),sphvg2(j))
               vrad3(j)=vrad(j)**3
            end do
         else 
            if ((rion(m).eq.rion(nsort(n+1)))
     &         .or.case2(nsort(ntype))) then
               case2(m)=case2(nsort(n+1))     
               reff(m)=reff(nsort(n+1))
               srfen(m)=srfen(nsort(n+1))               
               GOTO 200
            end if
            jmax=nrmin(nsort(n+1))+1
            do j=nrmin(m), nrmin(nsort(n+1))
c     Spherical integration.        
               call epssphv(j,nsph,xsph,ysph,zsph,wsph,vrad(j),
     &              xion,yion,zion,
     &              epsw,epsm,epsh,epsx,epsy,epsz,
     &              maxx,maxy,maxz,zmemb1,zmemb2,hmemb1,hmemb2,
     &              exp2,exp3,gnowat,
     &              dcel,nclz,ncyz,sphcal,epssph,
     &              sphv(j),sphvx(j),sphvy(j),sphvz(j),sphvm(j),
     &              omsphv(j),sphvg2(j))
               vrad3(j)=vrad(j)**3
            end do                        
         end if     

c     Compute integral contributions (first integration shell).
         j=nrmin(m)
c     Spherical integration (first integration shell).               
         call epssphv(j,nsph,xsph,ysph,zsph,wsph,rion(m),
     &        xion,yion,zion,
     &        epsw,epsm,epsh,epsx,epsy,epsz,
     &        maxx,maxy,maxz,zmemb1,zmemb2,hmemb1,hmemb2,
     &        exp2,exp3,gnowat,
     &        dcel,nclz,ncyz,sphcal,epssph,
     &        sphv(j),sphvx(j),sphvy(j),sphvz(j),sphvm(j),
     &        omsphv(j),sphvg2(j))
         vrad3(j)=rion(m)**3
         wrmin=0.5D0/rion(m)-0.5D0/vrad(j+1)
         
c     Compute integral contributions.
         do j=nrmin(m)+1, jmax      
            sphaux=sphv(j-1)+sphv(j)
            volaux=0.5D0*volfac*(vrad3(j)-vrad3(j-1))
c     Radial integral contributions.        
            vola(j)=volaux*sphaux     
            volxa(j)=volaux*(sphvx(j-1)+sphvx(j))
            volya(j)=volaux*(sphvy(j-1)+sphvy(j))
            volza(j)=volaux*(sphvz(j-1)+sphvz(j))
            volma(j)=volaux*(sphvm(j-1)+sphvm(j))
            if (j.eq.nrmin(m)+1) then
               intv1a(j)=wrmin*sphaux
            else
               intv1a(j)=wrad(j-1)*sphaux
            end if                         
         end do 
         
c     Compute the values of the integrals which are required for  
c     calculating the reaction field parameters of the ions. 
         do j=nrmin(m), nrad
            naux=(j-1)*nsph
            if (j.eq.nrmin(m)) then
               vol(j)=0.0D0
               volx(j)=0.0D0
               voly(j)=0.0D0
               volz(j)=0.0D0
               volm(j)=0.0D0
               intv1(j)=1.0D0/rion(m)
               intvg1(j)=0.0D0
               intvg2(j)=0.0D0
               gam1=0.0D0
               do i=1, nsph 
                  ngird=naux+i
                  if (sphcal(ngird)) then
                     dirg2(i)=0.0D0
                  else
                     dirg2(i)=1.0D0
                  end if               
               end do
               gam2=2.0D0*omsphv(j)**exp2
            else
               vol(j)=vol(j-1)+vola(j)
               volx(j)=volx(j-1)+volxa(j)
               voly(j)=voly(j-1)+volya(j)
               volz(j)=volz(j-1)+volza(j)
               volm(j)=volm(j-1)+volma(j)
               intv1(j)=intv1(j-1)-intv1a(j)
               gam1=2.0D0*dabs(1.0D0-1.0D0/(vrad(j)*intv1(j)))**exp1
               sphg2=0.0D0
               do i=1, nsph 
                  ngird=naux+i
                  if (sphcal(ngird)) then
                     dirg2(i)=0.0D0
                  else
                     sphg2=sphg2+wsph(i)*dirg2(i)    
                  end if              
               end do
               gam2=2.0D0*sphg2**exp2
            end if
c     Spherical integration.        
            sphvg1(j)=0.0D0 
            sphvg2(j)=0.0D0                       
            do i=1, nsph 
               ngird=naux+i
               if (sphcal(ngird)) then
                  epsaux=epssph(ngird)
                  sphvg1(j)=sphvg1(j)+wsph(i)
     &               *(epsw-epsaux)/(gam1*epsw+epsaux) 
                  sphvg2(j)=sphvg2(j)+wsph(i)
     &               *(epsw-epsaux)/(gam2*epsw+epsaux)
               end if
            end do
            sphvg1(j)=sphvg1(j)*(1.0D0+gam1)/epsw
            sphvg2(j)=sphvg2(j)*(1.0D0+gam2)/epsw
c     Radial integration.
            if (j.eq.nrmin(m)+1) then
               intvg1(j)=intvg1(j-1)
     &            +wrmin*(sphvg1(j-1)+sphvg1(j))            
               intvg2(j)=intvg2(j-1)
     &            +wrmin*(sphvg2(j-1)+sphvg2(j))
            else if (j.gt.nrmin(m)+1) then             
               intvg1(j)=intvg1(j-1)   
     &            +wrad(j-1)*(sphvg1(j-1)+sphvg1(j))
               intvg2(j)=intvg2(j-1)
     &            +wrad(j-1)*(sphvg2(j-1)+sphvg2(j)) 
            end if
         end do
         if ((intvg1(nrad).gt.intvg2(nrad)).or.
     &       (omsphv(nrmin(m)).lt.0.99999d0)) then
            intv0(m)=intvg1(nrad)
            case2(m)=.false.
         else 
            intv0(m)=intvg2(nrad)
            case2(m)=.true.
         end if
         
c     Calculate the effective radii reff.        
         if (vol(nrad).lt.maxdif) then
            rrat(m)=1.0D0
            reff(m)=0.0D0
         else 
            if (vol(nrad).lt.volref) then 
               volmax=vol(nrad) 
            else   
               volmax=volref
            end if 
            jvol=nrad+1
            do j=nrad, nrmin(m), -1 
               if (vol(j).gt.volmax) then
                  jvol=j
               end if               
            end do 
            if (jvol.eq.nrad+1) then
               reff(m)=rmax
               rrat(m)=rmax               
            else if (jvol.eq.nrmin(m)) then
               reff(m)=rion(m)
               rrat(m)=rion(m)               
            else if (jvol.gt.nrmin(m)) then
               w=(vol(jvol)-vol(jvol-1))
               w1=(vol(jvol)-volmax)/w
               w2=(volmax-vol(jvol-1))/w               
               reff(m)=volm(jvol-1)/vol(jvol-1)*w1
     &                +volm(jvol)/vol(jvol)*w2
               rratx=  volx(jvol-1)/vol(jvol-1)*w1
     &                +volx(jvol)/vol(jvol)*w2
               rraty=  voly(jvol-1)/vol(jvol-1)*w1
     &                +voly(jvol)/vol(jvol)*w2
               rratz=  volz(jvol-1)/vol(jvol-1)*w1
     &                +volz(jvol)/vol(jvol)*w2
               rrat(m)=DSQRT(rratx*rratx+rraty*rraty
     &            +rratz*rratz)
            end if 
            rrat(m)=rrat(m)/reff(m) 
         end if
!         reff(m)=reff(m)*dabs(1.0D0-rrat(m))
!     &      *(omsphv(nrmin(m))*intvg2(nrad)/intv0(m))**1.5D0

         if (intv0(m).gt.0.000001D0) then 
           reff(m)=reff(m)*dabs(1.0D0-rrat(m))
     &      *(omsphv(nrmin(m))*intvg2(nrad)/intv0(m))**1.5D0
         else
           reff(m)=reff(m)*dabs(1.0D0-rrat(m))
         end if
     
c     Self reactin field energy integral. 
         srfen(m)=intv0(m)
       
c 40   format(f6.2, 2f20.10) 
c      write(*,40) zion, DSQRT(srfen(m)), DSQRT(reff(m))
c      write(*,*) zion, reff(m), srfen(m)
                    
 200  end do
      
      end
      
            
c==============================================================================
      subroutine epssphv(j,nsph,xsph,ysph,zsph,wsph,rad,
     &           xion,yion,zion,
     &           epsw,epsm,epsh,epsx,epsy,epsz,
     &           maxx,maxy,maxz,zmemb1,zmemb2,hmemb1,hmemb2,
     &           exp2,exp3,gnowat,
     &           dcel,nclz,ncyz,sphcal,epssph,
     &           sphv,sphvx,sphvy,sphvz,sphvm,omsphv,sphvg2)
c==============================================================================
c
c     2006
c     
c     This subroutine calculates epssph and spherical integrals for
c     subroutine rfint. 
c  
c     Input:   j (number of radial integration point)
c              nsph (number of spherical integration points)
c              xsph (x-components of spherical integration points)
c              ysph (y-components of spherical integration points)
c              zsph (z-components of spherical integration points)
c              wsph (weights of spherical integration points)
c              rad (radius of the j-th radial integration point) 
c              xion, yion, zion (position of the ion)
c              epsw (bulk solvent dielectric constant)
c              epsm (membrane dielectric constant)
c              epsh (membrane headgroup dielectric constant)
c              epsx, epsy, epsz (dielectric constants between grid points)
c              maxx, maxy, maxz (size of the grid)
c              zmemb1, zmemb2, hmemb1, hmemb2 (membeane parameters)
c              exp2, exp3 (constants required for calculating sphvg2)
c              gnowat (grid with the number of next neighbour epsx, epsy, 
c                 and epsz values which are different from epsw)
c              dcel (grid spacing)
c              nclz (number of grid points in z direction)
c              ncyz (nclz*nclz)
c     Outpout: sphcal (integrate or do not integrate at the 
c                 gridpoints of the numerical integrator)
c              epssph (values of the local dielectric constant at the 
c                 gridpoints of the numerical integrator)
c              sphv (result of spherical integration)
c              sphvx (result of spherical integration, x-coordinate)
c              sphvy (result of spherical integration, y-coordinate)
c              sphvz (result of spherical integration, z-coordinate)
c              sphvm (result of spherical integration, magnitude)
c              omsphv (1 - sphv)
c              sphvg2 (result of spherical integration, dependent on
c                 the local dielectric constant)
c     Calls:   
c
 
      implicit none
       
      integer*4 j, nsph, nclz, ncyz
      real*8 xsph(nsph), ysph(nsph), zsph(nsph), wsph(nsph)
      real*8 rad
      real*8 xion, yion, zion
      real*8 epsw, epsm, epsh
      real*4 epsx(*), epsy(*), epsz(*)       
      real*8 maxx, maxy, maxz
      real*8 zmemb1, zmemb2, hmemb1, hmemb2
      real*8 exp2, exp3 
      integer*4 gnowat(*)
      real*8 dcel
         
      logical sphcal(*)
      real*8 epssph(*), sphv, sphvx, sphvy, sphvz, sphvm 
      real*8 omsphv, sphvg2
      
      integer*4 i, ngird, naux, jx, jy, jz, gpos    
      real*8 xpos, ypos, zpos 
      real*8 diffx, diffy, diffz
      real*8 maxdif, aux
      real*8 gam2, epsaux
      parameter (maxdif=1.0D-8)

      sphv=0.0D0
      sphvx=0.0D0 
      sphvy=0.0D0 
      sphvz=0.0D0 
      sphvm=0.0D0
      naux=(j-1)*nsph
      do i=1, nsph
         xpos=rad*xsph(i)+xion
         ypos=rad*ysph(i)+yion
         zpos=rad*zsph(i)+zion
         ngird=naux+i
c     Determine the value of the dielectric constant at the spherical
c     integration grid points.            
         if(zpos.ge.zmemb1.and.zpos.le.zmemb2) then
            if(zpos.ge.hmemb1.and.zpos.le.hmemb2) then
               epssph(ngird)=epsm
            else
               epssph(ngird)=epsh
            end if
         else
            epssph(ngird)=epsw
         end if
         if ((xpos.gt.0.0D0).and.(xpos.lt.maxx).AND.
     &      (ypos.gt.0.0D0).and.(ypos.lt.maxy).AND.
     &      (zpos.gt.0.0D0).and.(zpos.lt.maxz)) then
            jx=idnint(xpos/dcel)+1
            jy=idnint(ypos/dcel)+1
            jz=idnint(zpos/dcel)+1
            gpos=(jx-1)*ncyz+(jy-1)*nclz+jz
            if (gnowat(gpos).lt.1) then
               epssph(ngird)=epsw
            else
               diffx=xpos-(jx-1)*dcel
               diffy=ypos-(jy-1)*dcel
               diffz=zpos-(jz-1)*dcel
               if (dabs(diffz).ge.dabs(diffy).and.
     &             dabs(diffz).ge.dabs(diffx)) then
                  if (diffz.ge.0.0D0) then
                     epssph(ngird)=epsz(gpos)
                  else 
                     epssph(ngird)=epsz(gpos-1)
                  end if
               else if (dabs(diffy).gt.dabs(diffz).and.
     &                  dabs(diffy).ge.dabs(diffx)) then
                  if (diffy.ge.0.0D0) then
                     epssph(ngird)=epsy(gpos)
                  else 
                     epssph(ngird)=epsy(gpos-nclz)
                  end if
               else 
                  if (diffx.ge.0.0D0) then
                     epssph(ngird)=epsx(gpos)
                  else 
                     epssph(ngird)=epsx(gpos-ncyz)
                  end if
               end if
            end if
         end if
c     Spherical integration.        
         if (dabs(epssph(ngird)-epsw).gt.maxdif) then
            sphv=sphv+wsph(i)
            aux=wsph(i)*rad
            sphvx=sphvx+aux*xsph(i) 
            sphvy=sphvy+aux*ysph(i) 
            sphvz=sphvz+aux*zsph(i)
            sphvm=sphvm+aux           
            sphcal(ngird)=.true.
         else
            sphcal(ngird)=.false.   
         end if
      end do
      
      omsphv=dabs(1.0D0-sphv)
c      gam2=2.0D0*omsphv**(exp2+exp3*omsphv)              
c      sphvg2=0.0D0                       
c      naux=(j-1)*nsph
c      do i=1, nsph 
c         ngird=naux+i
c         if (sphcal(ngird)) then
c            epsaux=epssph(ngird)
c            sphvg2=sphvg2+wsph(i)
c     &         *(epsw-epsaux)/(gam2*epsw+epsaux)
c         end if
c      end do
c      sphvg2=sphvg2*(1.0D0+gam2)/epsw
           
      end
      
            
c==============================================================================
      subroutine intradsph(ntype,nsort,nrad,rmax,rion,ordsph,nsph,
     &           nrmin,vrad,wrad,xsph,ysph,zsph,wsph)
c==============================================================================
c
c     2006
c
c     This subroutine yields integration points and weights for 
c     numerical integration in spherical coordinates.
c
c     Input:   ntype (number of ion types)
c              nsort (sorted list of ion type numbers with growing Born radii)
c              nrad (number of radial integration points)
c              rmax (maximum radius for numerical integration)
c              rion (Born radii of the ions)
c              ordsph (order of the Gauss quadrature scheme)
c              nsph (number of spherical integration points)
c     Outpout: nrmin (grid numbers of the first integration points
c                 for the ion types) 
c              vrad (values of radial integration points) 
c              wrad (nonstandard weights of radial integration points)
c              xsph (x-components of spherical integration points)
c              ysph (y-components of spherical integration points)
c              zsph (z-components of spherical integration points)
c              wsph (weights of spherical integration points)
c     Calls:   intradaux
c              intsphaux
c

      implicit none
      include 'rfintsph.fcm'
      
      integer*4 ntype, nsort(*), nrad, ordsph, nsph
      real*8 rmax, rion(*)
      
      integer*4 nrmin(*)
      real*8 vrad(*), wrad(*)
      real*8 xsph(nsph), ysph(nsph), zsph(nsph), wsph(nsph)      
      
      call intradaux(ntype,nsort,nrad,rmax,rion,nrmin,vrad,wrad)
      
      if (ordsph.eq.17) then
         call intsphaux(n17,m17,x17,y17,w17,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.23) then 
         call intsphaux(n23,m23,x23,y23,w23,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.29) then 
         call intsphaux(n29,m29,x29,y29,w29,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.35) then 
         call intsphaux(n35,m35,x35,y35,w35,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.41) then 
         call intsphaux(n41,m41,x41,y41,w41,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.47) then 
         call intsphaux(n47,m47,x47,y47,w47,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.53) then 
         call intsphaux(n53,m53,x53,y53,w53,xsph,ysph,zsph,wsph)
      else if (ordsph.eq.59) then 
         call intsphaux(n59,m59,x59,y59,w59,xsph,ysph,zsph,wsph)
      end if  
      
      end   
 
 
c==============================================================================
      subroutine intradaux(ntype,nsort,nrad,rmax,rion,nrmin,vrad,wrad)
c==============================================================================
c
c     2006
c
c     This auxiliary subroutine is called by intradsph. It calculates 
c     integration points and weights for a non-standard numerical 
c     integration scheme applied to perform radial integration. 
c
c     Input:   ntype (number of ion types)
c              nsort (sorted list of ion type numbers with growing Born radii)
c              nrad (number of radial integration points)
c              rmax (maximum radius for numerical integration)
c              rion (Born radii of the ions)
c     Outpout: nrmin (grid numbers of the first integration points
c                 for the ion types) 
c              vrad (values of radial integration points) 
c              wrad (weights of radial integration points)
c             
    
      implicit none
      
      integer*4 ntype, nsort(*), nrad
      real*8 rmax, rion(*)
      
      integer*4 nrmin(*)    
      real*8 vrad(*), wrad(*)
      
      logical test
c     test: Tests checking the correctness of the points and 
c     weights for radial integration (if test=.true.)
      parameter (test=.false.)   
c      parameter (test=.true.)             
      integer*4 j, k, m, n
      integer*4 intaux(nrad), norm
      real*8 eps
      parameter (eps=0.00000000000001D0) 
      
      n=nrad
      
      intaux(1)=1
      do j=2,n
         intaux(j)=intaux(j-1)+j
      end do
      
      norm=intaux(n)
      do j=1,n
         vrad(j)=rmax*DBle(intaux(j))/DBle(norm)  
      end do
      
      do j=1,n-1
         wrad(j)=0.5D0/vrad(j)-0.5D0/vrad(j+1)      
      end do      

      do k=1,ntype
         m=nsort(k)
         nrmin(m)=1
         do j=1,n-1             
            if ((rion(m).lt.vrad(j+1))
     &         .and.(rion(m).ge.vrad(j))) then
               nrmin(m)=j
            end if
         end do         
      end do      
      
      if (test) then
         do j=1,n-1
            write(*,*) j,vrad(j),wrad(j)              
         end do
         write(*,*) n,vrad(n)
         do k=1,ntype
            write(*,*) k,nsort(k),nrmin(nsort(k))
         end do            
      end if
               
      end   
      
         
c==============================================================================
      subroutine intsphaux(nsph,mdata,xdata,ydata,wdata,
     &           xsph,ysph,zsph,wsph)
c==============================================================================
c
c     2006
c
c     This auxiliary subroutine is called by intradsph. It calculates the 
c     integration points and weights of a Gauss quadrature integration 
c     scheme on the unit sphere [B. Delley, J. Comput. Chem. 17(9), 
c     1152-1155 (1996)].
c
c     Input:   nsph (number of spherical integration points)
c              mdata (size of data set) 
c              xdata (x-components of data points)
c              ydata (y-components of data points)
c              wdata (weights of data points)
c     Outpout: xsph (x-components of spherical integration points)
c              ysph (y-components of spherical integration points)
c              zsph (z-components of spherical integration points)
c              wsph (weights of spherical integration points)
c
      
      implicit none
      
      integer*4 nsph, mdata
      real*8 xdata(mdata), ydata(mdata), wdata(mdata)
      
      real*8 xsph(nsph), ysph(nsph), zsph(nsph), wsph(nsph)
           
      logical uneq, error, test
c     test: Tests checking the correctness of the points and 
c     weights for numerical integration (if test=.true.)    
      parameter (test=.false.) 
c      parameter (test=.true.) 
      integer*4 i, j, k, l
      real*8 eps
      parameter (eps=0.00000000000001D0) 
      real*8 xrot, yrot
      parameter (xrot=0.20D0)
      parameter (yrot=0.21D0)
c      parameter (xrot=0.00D0)
c      parameter (yrot=0.00D0)
      real*8 xaux(48), yaux(48), zaux(48)
      real*8 xorig(nsph), yorig(nsph), zorig(nsph)      
      real*8 mag, xsum, ysum, zsum, wsum 
      real*8 wxsum, wysum, wzsum
      
c     Generating the full set of nsph integration points on the 
c     unit sphere.      
      i=0
      do j=1, mdata
         xaux(1)=xdata(j)
         yaux(1)=ydata(j)
         zaux(1)=DSQRT(1.0D0-xaux(1)*xaux(1)-yaux(1)*yaux(1))
         xaux(2)=yaux(1)
         yaux(2)=xaux(1)
         zaux(2)=zaux(1)
         xaux(3)=yaux(1)
         yaux(3)=zaux(1)
         zaux(3)=xaux(1)
         xaux(4)=xaux(1)
         yaux(4)=zaux(1)
         zaux(4)=yaux(1)
         xaux(5)=zaux(1)
         yaux(5)=xaux(1)
         zaux(5)=yaux(1)
         xaux(6)=zaux(1)
         yaux(6)=yaux(1)
         zaux(6)=xaux(1)
         k=1
         do l=1, 6
            xaux(k*6+l)=-xaux(l)
            yaux(k*6+l)=-yaux(l)
            zaux(k*6+l)=-zaux(l)
         end do                  
         k=2
         do l=1, 6
            xaux(k*6+l)=-xaux(l)
            yaux(k*6+l)=+yaux(l)
            zaux(k*6+l)=+zaux(l)
         end do                  
         k=3
         do l=1, 6
            xaux(k*6+l)=+xaux(l)
            yaux(k*6+l)=-yaux(l)
            zaux(k*6+l)=+zaux(l)
         end do                  
         k=4
         do l=1, 6
            xaux(k*6+l)=+xaux(l)
            yaux(k*6+l)=+yaux(l)
            zaux(k*6+l)=-zaux(l)
         end do 
         k=5                 
         do l=1, 6
            xaux(k*6+l)=+xaux(l)
            yaux(k*6+l)=-yaux(l)
            zaux(k*6+l)=-zaux(l)
         end do                  
         k=6                 
         do l=1, 6
            xaux(k*6+l)=-xaux(l)
            yaux(k*6+l)=+yaux(l)
            zaux(k*6+l)=-zaux(l)
         end do                  
         k=7                 
         do l=1, 6
            xaux(k*6+l)=-xaux(l)
            yaux(k*6+l)=-yaux(l)
            zaux(k*6+l)=+zaux(l)
         end do

         i=i+1
         if (i.le.nsph) then       
            xsph(i)=xaux(1)
            ysph(i)=yaux(1)
            zsph(i)=zaux(1)
            wsph(i)=wdata(j) 
         end if
         do k=2, 48
            uneq=.true. 
            do l=1, k-1
               if ((dabs(xaux(l)-xaux(k)).lt.eps).and.
     &             (dabs(yaux(l)-yaux(k)).lt.eps).and.
     &             (dabs(zaux(l)-zaux(k)).lt.eps)) then
                  uneq=.false.                  
                  GOTO 10
               end if 
            end do
 10         continue  
            if (uneq) then
               i=i+1
               if (i.le.nsph) then       
                  xsph(i)=xaux(k)
                  ysph(i)=yaux(k)
                  zsph(i)=zaux(k)
                  wsph(i)=wdata(j) 
               end if
            end if 
         end do
      end do 

c     Rotating the integration points.
c     Rotation 1 (the x-axis is the rotation axis 
c     and xrot is the rotation angle).
      do k=1, nsph 
         yorig(k)=ysph(k)
         zorig(k)=zsph(k)
      end do
      do k=1, nsph 
         ysph(k)=yorig(k)*DCOS(xrot)+zorig(k)*DSIN(xrot)
         zsph(k)=zorig(k)*DCOS(xrot)-yorig(k)*DSIN(xrot)
      end do
c     Rotation 2 (the y-axis is the rotation axis 
c     and yrot is the rotation angle).
      do k=1, nsph 
         xorig(k)=xsph(k)
         zorig(k)=zsph(k)
      end do
      do k=1, nsph 
         xsph(k)=xorig(k)*DCOS(yrot)+zorig(k)*DSIN(yrot)
         zsph(k)=zorig(k)*DCOS(yrot)-xorig(k)*DSIN(yrot)
      end do
                        
c     Tests checking some basic properties of the generated vector
c     components xsph(i), ysph(i), zsph(i) and weights wsph(i).      
      error=.false.
      if (test) then
         if (i.NE.nsph) then
            error=.true.
 30         format('Error: Subroutine intsph generated ',I4,
     &         ' instead of ',I4,' vectors for spherical integration!') 
            write(*,30) i, nsph 
         else
            xsum=0.0D0
            ysum=0.0D0
            zsum=0.0D0
            wsum=0.0D0
            wxsum=0.0D0
            wysum=0.0D0
            wzsum=0.0D0    
            do i=1, nsph
               mag=DSQRT(xsph(i)*xsph(i)+ysph(i)*ysph(i)
     &            +zsph(i)*zsph(i))
               xsum=xsum+xsph(i)
               ysum=ysum+ysph(i)
               zsum=zsum+zsph(i) 
               wsum=wsum+wsph(i) 
               wxsum=wxsum+xsph(i)*wsph(i)
               wysum=wysum+ysph(i)*wsph(i)   
               wzsum=wzsum+zsph(i)*wsph(i)   
               if ((mag.gt.(1.0D0+eps)).or.(mag.lt.(1.0D0-eps))) then
                  error=.true.
 40               format('Error: Subroutine intsph generated a ',
     &               'vector with a wrong magnitude!')
                  write(*,40)
               end if 
            end do           
            if ((xsum.gt.(eps)).or.(xsum.lt.(-eps)).or.
     &          (ysum.gt.(eps)).or.(ysum.lt.(-eps)).or.
     &          (zsum.gt.(eps)).or.(zsum.lt.(-eps))) then
               error=.true.
 50            format('Error: Subroutine intsph generated vectors ',
     &            'with a wrong vector sum (',
     &            D27.20,', ',D27.20,', ',D27.20,')!')
               write(*,50) xsum, ysum, zsum     
            end if
            if ((wsum.gt.(1.0D0+eps)).or.(wsum.lt.(1.0D0-eps))) then 
               error=.true.
 60            format('Error: Subroutine intsph uses weights ',
     &            'yielding a wrong weight sum ',
     &            D27.20,'!')
               write(*,60) wsum   
            end if
            if ((wxsum.gt.(wysum+eps)).or.(wxsum.lt.(wysum-eps)).or.
     &          (wxsum.gt.(wzsum+eps)).or.(wxsum.lt.(wzsum-eps))) then 
               error=.true.
 70            format('Error: Subroutine intsph generated a ',
     &            'direction dependent rusult in an isotropy test (x ',
     &            D27.20,', y ',D27.20,', z ',D27.20,')!')
               write(*,70) wxsum, wysum, wzsum   
            end if
         end if
      end if
            
      end
      
      
      subroutine gridz(itype,nclx,ncly,nclz,grid)
      integer*4 nclx,ncly,nclz
      real*4 grid(*)
      integer*4 j,itype
      write(*,*) itype
      do j=nclz,1,-1  
         write(*,*) j, grid((nclx-1)/2*ncly*nclz 
     &      +(ncly-1)/2*nclz+j+(itype-1)*nclx*ncly*nclz)**2
c         write(*,*) j, grid((j-1)*ncly*nclz 
c     &      +(ncly-1)/2*nclz+nclz/2+(itype-1)*nclx*ncly*nclz)**2
      end do
      end
         
         
   
