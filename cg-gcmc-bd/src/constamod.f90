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

! Physical constants
module constamod
implicit none
real, parameter :: pi=3.14159265358979323846264338327950288419716939937510
real, parameter :: twopi=2.0*pi
real, parameter :: kboltz=1.380648813E-23 ! [J/K] Uncertainty last 2 digits. 2010 CODATA NIST (kboltz in kcal/mol = 0.00198720416527)
real, parameter :: angstrom=1E-10
real, parameter :: liter=1E-3
real, parameter :: avogadro=6.0221412927E23 ! Uncertainty last 2 digits. 2010 CODATA NIST
real, parameter :: kcal=4184.0, kcalmol=kcal/avogadro
real, parameter :: pico=1E-12,psec=1e-12,fsec=1e-15,ipico=1E12
real, parameter :: Coulomb=1.60217656535E-19 ! Uncertainty last 2 digits. US National Institute of Standards and Technology. June 2011
real, parameter :: epermit=8.854187817620E-12 ! Uncertainty last 2 digits. 2006 CODATA NIST (F/m or A2.s4/kg/m3 or C2/N/m2 or C/V/m)
real, parameter :: celec = Coulomb**2/(4.0*pi*epermit*kcalmol*angstrom) ! 332.06371
real, parameter :: celec2 = Coulomb/(4.0*pi*epermit*angstrom) ! 332.06371
real, parameter :: rsmall = 1E-3
real, parameter :: const = epermit*kboltz*1.0e17/(2.0*avogadro*Coulomb**2)
real, parameter :: kba3bar = kboltz*1e25 ! in Ang**3*bar/K
end module
