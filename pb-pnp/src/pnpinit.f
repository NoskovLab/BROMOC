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

      SUBROUTINE PNPINIT
c------------------------------------------------------------------------
c     Initialize all parameters
c
      implicit none
      include 'pnp.fcm'
      include 'consta.fcm'

c
      QPBEQ=.false.
      QPNP =.false.
      QMMIJ=.false.
      QRFPAR=.false.
c
      NATOM=0                   ! number of explicit atoms
      NTYPE=0                   ! number of ion types


c LINEAR or NON-LINEAR PBEQ solver?
c =================================

      Qnonlinear =.false.       ! non-linear PBEQ solver
      Qpartlinear=.false.       ! partially linearized PBEQ solver
      Qunder     =.false.       ! under-relaxation for non-linear PBEQ solver

c ITERATION PARAMETERS for the SOR (Successive OverRelaxation) method
c ====================
c maxpnp  : maximum number of iterations for both PBEQ and PNP 
c maxphi  : maximum number of iterations for PBEQ
c maxcion : maximum number of iterations for NP
c tolcion : tolerance value for PNP solver
c tolphi  : tolerance value for PBEQ solver
c lambda1 : mixing factor for potential
c lambda2 : mixing factor for concentration

      maxpnp =100
      maxphi =2000
      maxcion=2000
      tolphi =2.0D-6
      tolcion=2.0D-10 
      lambda1=1.0d0
      lambda2=1.0d0

c GRID PARAMETERS
c ===============
c ncel, nclx, ncly, nclz : number of grid points
c xbcen, ybcen, zbcen : center of box
c dcel : grid spacing

      ncel =65
      nclx =ncel
      ncly =ncel
      nclz =ncel
      xbcen=0.0D0
      ybcen=0.0D0
      zbcen=0.0D0
      dcel =0.5D0
      
c PHYSICAL PARAMETERS
c ===================
c watr  : solvent probe radius
c ionr  : ion exclusion radius (Stern layer)
c epsw, epsp : dielectric constants of water and protein
c temp  : temperature in [K]
c conc  : concentration in [mol/L]

      watr=  0.0D0
      ionr=  0.0D0
      epsw= 80.0D0
      epsp=  1.0D0 
      temp=300.0D0 
      conc=  0.0D0 

c MEMBRANE
c ========
c epsm  : dielectric constants of membrane
c epsh  : dielectric constants of head groups
c tmemb : thickness of membrane including head groups
c htmemb: thickness of head groups
c zmemb : membrane position along Z-axis
c vmemb : transmembrane potentail in [volts]

      epsm  =1.0D0 
      epsh  =1.0D0 
      tmemb =0.0D0 
      zmemb =0.0D0 
      vmemb =0.0D0 
      htmemb=0.0D0 

c SPHERE PARAMETERS
c ===================
c rsphe,xyzsphe : radius and positions of a sphere
c epss : dielectric constant of a spherical cavity

      rsphe=0.0D0 
      xsphe=0.0D0 
      ysphe=0.0D0 
      zsphe=0.0D0 
      epss =1.0D0 
      Qskap=.false.

c CYLINDER PARAMETERS
c ===================
c rcyln, hcyln, xyzcyln : radius, height, and positions of a cylinder
c epsc : dielectric constant inside cylinder

      rcyln=0.0D0 
      hcyln=0.0D0 
      xcyln=0.0D0 
      ycyln=0.0D0 
      zcyln=0.0D0 
      epsc =1.0D0 
      Qckap=.false.

c ELLIPTIC CYLINDER PARAMETERS
c ===================
c Ax, By, tmin, hecyln, xyzecyln : radius, height, and positions of a cylinder
c epsec : dielectric constant inside cylinder

      Ax=0.0D0
      By=0.0D0
      tmin=1.0D0
      hecyln=0.0D0
      xecyln=0.0D0
      yecyln=0.0D0
      zecyln=0.0D0
      epsec =1.0D0
      Qeckap=.false.

c BOX PARAMETERS
c ===================
c b{x,y,z}max and b{x,y,z} : X, Y, and Z dimensions of an orthorhombic box
c epsb : dielectric constant of a box

      bxmax=0.0D0 
      bxmin=0.0D0 
      bymax=0.0D0 
      bymin=0.0D0 
      bzmax=0.0D0 
      bzmin=0.0D0 
      epsb =1.0D0 
      Qbkap=.false.
c
      RETURN
      END
