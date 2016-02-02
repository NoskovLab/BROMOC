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

!     common block for generalized solvent boundary potential for program simul
!
!     srdist                                   radius of a sphere  
!     ntpol                                    number of basis functions
!     mmij(*)                                  generalized reaction field matrix
!     bnorm(*)                                 normalization constant for basis functions
!     coef(*)                                  generalized multipole moments
!     lstpl(*)                                 list of l index number in spherical harmonics y_{lm}
!     lstpm(*)                                 list of m index number in spherical harmonics y_{lm}
!     phix(*)                                  static field 
!     egsbpa                                   static field energy
!     egsbpb                                   reaction field energy
!     phiv(*)                                  grid-based uniform repulsive potential 
!                                              (zero inside and one outside protein and membrane)
!     evdwgd                                   repulsive energy
!     svdw                                     scaling factor of phiv
!     nclx, ncly, nclz                         grid parameters
!     dcel                                     mesh size (grid spacing)
!     tranx = half*(nclx-1)*dcel               origin of of corner of the grid in x
!     trany = half*(ncly-1)*dcel               origin of of corner of the grid in y
!     tranz = half*(nclz-1)*dcel               origin of of corner of the grid in z
!     vecphiv(maxopen)                          store different units for grid-based repulsion potential  
!
! gsbp parameters   
module gsbpmod
use ioxmod
implicit none 
!integer   dgrid1,dgrid2,dbasis
!parameter ( dgrid1=201*201*301,dgrid2=201*201*301*2,dbasis=500 )
!parameter ( dgrid1=181*181*181,dgrid2=181*181*181*2,dbasis=500 )
!integer, parameter :: dbasis=500
integer   nclx1,ncly1,nclz1,nclx2,ncly2,nclz2,ntpol
integer,allocatable ::   lstpl(:),lstpm(:)
integer*4,allocatable :: lstpx(:),lstpy(:),lstpz(:)
real,allocatable ::    bnorm(:),coef(:)
real*8,allocatable ::  mmij(:)
real    srdist
integer*4 xnpol,ynpol,znpol,nmpol
integer   vecphiv(maxopen)
real    dcel1,tranx1,trany1,tranz1,dcel2,tranx2,trany2,tranz2
real    xbcen1,ybcen1,zbcen1,xbcen2,ybcen2,zbcen2,idcel2,idcel1
integer    nclx3,ncly3,nclz3                                   ! rfpar
real    dcel3,tranx3,trany3,tranz3,xbcen3,ybcen3,zbcen3,idcel3 ! rfpar
integer nclx4,ncly4,nclz4                                      ! charge density
real*8  dcel4,tranx4,trany4,tranz4,xbcen4,ybcen4,zbcen4,idcel4 ! charge density
real*4,allocatable :: chden(:)                                 ! charge density
logical*1 Qchden,Qchdencnt,Qchdenorm                           ! charge density
real    sqrfac,reffac
real*4,allocatable :: gsrfen(:),greff(:)
real    rbxmax,rbxmin,rbymax,rbymin,rbzmax,rbzmin
real    xscal,yscal,zscal,rfscal
!real*16 egsbpb,egsbpa,evdwgd
real egsbpb,egsbpa,evdwgd
real    svdw,vzmax,vzmin
real*4,allocatable ::  phix(:)
integer*1,allocatable :: phiv(:)
real    xmin,xmax,ymin,ymax,zmin,zmax
real    thold27, thold8
character*8 shapes
logical*1   Qmmij,Qphix,Qphiv,Qtrln,Qnmcden,Qsvdw,Qrfpar,Qrfpsin
real,allocatable ::    scal(:)
end module


