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

subroutine bonds 
! Interaction sites are stored in the order (5'-3') 
! B1-S1-P1--B2-S2-P2--B3-S3-P3..., where B=Base, S=Sugar 
! and P=Phosphate. Note that the bond between  a phospate 
! and a sugar belonging to the same nucleotide is S(5')-P, 
! while the bond S(3')-P joins together neighboring 
! residues. 
use stdiomod
use errormod 
use grandmod
use nucleotmod
implicit none
!     Local variables      
integer i, i1, i2, strand1, strand2, strand3 
character*2 namesite1, namesite2, namesite3
character*5 namepr
!     Function to generate the index for a supervector
integer indexi

!     Initializations      
nbond = 0 
do i = 1, mmxsites
  bond(i) = .false.
enddo
      
do i = 1, nsites-1 
  namesite1 = namsite(i)
  namesite2 = namsite(i+1)
  strand1 = strand(i)
  strand2 = strand(i+1)
  if ((i+2).le.nsites) then
    strand3 = strand(i+2)
    namesite3 = namsite(i+2)
  endif  
  if (namesite1.eq.'Ab'.or.namesite1.eq.'Tb'.or.namesite1.eq.'Cb'.or.namesite1.eq.'Gb') then
    nbond = nbond + 1
    if (nbond.gt.maxbond) call error ('bonds', 'The number of bonds exceeds the maximum value', faterr)
    typbond(nbond) = 0 ! intranucleotide bond
    sitebond(nbond,1) = i
    sitebond(nbond,2) = i + 1
    bond(indexi(i,i+1)) = .true.
    if (namesite2.eq.'S ') then 
      if (namesite1.eq.'Ab') then ! S-Ab
        distbond(nbond) = dSAb
      else if (namesite1.eq.'Tb') then ! S-Tb
        distbond(nbond) = dSTb
      else if (namesite1.eq.'Cb') then ! S-Cb
        distbond(nbond) = dSCb
      else ! S-Gb      
        distbond(nbond) = dSGb
      endif  
    else
      call error ('bonds', 'INCORRECT ORDER FOR SITES', faterr)  
    endif  
  else if (namesite1.eq.'S ') then
    if (strand1.eq.strand2) then   
      nbond = nbond + 1
      if (nbond.gt.maxbond) call error ('bonds', 'The number of bonds exceeds the maximum value', faterr)
      sitebond(nbond,1) = i
      sitebond(nbond,2) = i + 1
      bond(indexi(i,i+1)) = .true.      
      if (namesite2.eq.'P ') then ! S(5')-P
        if (Qinvstr) then
          distbond(nbond) = dPS5
          typbond(nbond) = 0 ! intranucleotide bond
        else
          distbond(nbond) = dPS3
          typbond(nbond) = 1 ! internucleotide bond
        endif
      else
        call error ('bonds', 'INCORRECT ORDER FOR SITES', faterr)
      endif 
    endif   
  else ! namesite1.eq.'P '          
    if ((i+2).le.nsites .and. strand1.eq.strand3) then
      nbond = nbond + 1
      if (nbond.gt.maxbond) call error ('bonds', 'The number of bonds exceeds the maximum value', faterr)
      sitebond(nbond,1) = i
      sitebond(nbond,2) = i + 2
      bond(indexi(i,i+2)) = .true. 
      if (namesite3.eq.'S ') then ! S(3')-P    
        if (Qinvstr) then  
          distbond(nbond) = dPS3
          typbond(nbond) = 1 ! internucleotide bond
        else
          distbond(nbond) = dPS5
          typbond(nbond) = 0 ! intranucleotide bond
        endif
      else
        call error ('bonds', 'INCORRECT ORDER FOR SITES', faterr) 
      endif
    endif        
  endif  
enddo

!     OUPUT
if (Qninfo) then
  write (outu,'(/6x,a)') 'DNA: INTRASTRAND BONDS'
  write (outu,'(6x,a)') '----------------------'
  write (outu,'(6x,a)') 'Nbond  strand  namesite1  namesite2  typbond natural bond length'
  do i = 1, nbond
    if (typbond(i).eq.0) then
       namepr = 'intra'
    else
       namepr = 'inter'
    endif 
    i1 = sitebond(i,1)
    i2 = sitebond(i,2)    
    write (outu,'(6x,i3,4x,i2,6x,a2,9x,a2,9x,a5,4x,f6.3)') i, strand(i1),namsite(i1),namsite(i2),namepr,distbond(i)
  enddo  
endif

return
end


