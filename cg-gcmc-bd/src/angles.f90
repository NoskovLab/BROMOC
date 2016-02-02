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

subroutine angles
!Interaction sites are stored in the order (5'-3') 
!B1-S1-P1--B2-S2-P2--B3-S3-P3..., 
!where B=Base, S=Sugar and P=Phosphate. Note that the bond 
!between a phospate and a sugar belonging to the same 
!nucleotide is S(5')-P, while the bond S(3')-P joins together 
!neighboring residues. 
use constamod
use stdiomod
use errormod
use grandmod
use nucleotmod
implicit none 
!Local variables   
real  cte, valueang   
integer i, i1, i2, i3, strand1, strand2, strand3, strand4
character*2 namesite1, namesite2, namesite3, namesite4
!Function to generate the index for a supervector
integer indexi

!Initializations
cte = pi/180.0      
nangle = 0
do i = 1, mmxsites
  angle(i) = .false.
enddo  
 
do i = 1, nsites-2 
  strand1 = strand(i)
  strand2 = strand(i+1)
  strand3 = strand(i+2)
  namesite1 = namsite(i)
  namesite2 = namsite(i+1)
  namesite3 = namsite(i+2)
  if ((i+3).le.nsites) then
    strand4 = strand(i+3)
    namesite4 = namsite(i+3)
  endif       
  if (namesite1.eq.'Ab'.or.namesite1.eq.'Tb'.or.namesite1.eq.'Cb'.or.namesite1.eq.'Gb') then  
    if (strand1.eq.strand3) then 
      nangle = nangle + 1      
      if (nangle.gt.maxang) call error ('angles', 'The number of bond angles exceeds the maximum value', faterr)
      siteangle(nangle,1) = i
      siteangle(nangle,2) = i + 1 ! sugar central site 
      siteangle(nangle,3) = i + 2
      angle(indexi(i,i+1)) = .true.
      angle(indexi(i,i+2)) = .true.
      angle(indexi(i+1,i+2)) = .true.
      if (namesite2.eq.'S ' .and. namesite3.eq.'P ') then  
        if (Qinvstr) then 
          if (namsite(i).eq.'Ab') then ! P-(5')S-Ab
            valangle(nangle) = phPSAb
          elseif (namsite(i).eq.'Tb') then ! P-(5')S-Tb
            valangle(nangle) = phPSTb
          elseif (namsite(i).eq.'Cb') then ! P-(5')S-Cb
            valangle(nangle) = phPSCb
          else ! P-(5')S-Gb      
            valangle(nangle) = phPSGb
          endif
        else
          if (namsite(i).eq.'Ab') then ! P-(3')S-Ab
            valangle(nangle) = phPSAb2
          elseif (namsite(i).eq.'Tb') then ! P-(3')S-Tb
            valangle(nangle) = phPSTb2
          elseif (namsite(i).eq.'Cb') then ! P-(3')S-Cb
            valangle(nangle) = phPSCb2
          else ! P-(5')S-Gb      
            valangle(nangle) = phPSGb2
          endif
        endif 
      else
        call error ('angles', 'INCORRECT ORDER FOR SITES',faterr)
      endif
    endif   
  else if (namesite1.eq.'S ') then  
    if ((i+3).le.nsites .and. strand1.eq.strand4) then
      nangle = nangle + 1     
      if (nangle.gt.maxang) call error ('angles', 'The number of bond angles exceeds the maximum value', faterr)
      siteangle(nangle,1) = i
      siteangle(nangle,2) = i + 1 ! phosphate central site
      siteangle(nangle,3) = i + 3
      angle(indexi(i,i+1)) = .true.
      angle(indexi(i,i+3)) = .true.
      angle(indexi(i+1,i+3)) = .true.
      if (namesite2.eq.'P ' .and. namesite4.eq.'S ') then ! S(5')-P-(3')S
        valangle(nangle) = phSPS
      else  
        call error ('angles', 'INCORRECT ORDER FOR SITES',faterr)
      endif 
    endif         
  else ! namesite1.eq.'P '
    if (namesite2.eq.'Ab'.or.namesite2.eq.'Tb'.or.namesite2.eq.'Cb'.or.namesite2.eq.'Gb'.and.namesite3.eq.'S ') then 
      nangle = nangle + 1
      if (nangle.gt.maxang) call error ('angles', 'The number of bond angles exceeds the maximum value', faterr)
      siteangle(nangle,1) = i 
      siteangle(nangle,2) = i + 2 ! sugar central site
      siteangle(nangle,3) = i + 1     
      angle(indexi(i,i+2)) = .true.
      angle(indexi(i,i+1)) = .true.
      angle(indexi(i+1,i+2)) = .true.
      if (Qinvstr) then
        if (namesite2.eq.'Ab') then !  P-(3')S-Ab
          valangle(nangle) = phPSAb2
        else if (namesite2.eq.'Tb') then !  P-(3')S-Tb 
          valangle(nangle) = phPSTb2
        else if (namesite2.eq.'Cb') then !  P-(3')S-Cb
          valangle(nangle) = phPSCb2 
        else ! P-(3')S-Gb
          valangle(nangle) = phPSGb2   
        endif 
      else
        if (namesite2.eq.'Ab') then !  P-(5')S-Ab
          valangle(nangle) = phPSAb
        else if (namesite2.eq.'Tb') then !  P-(5')S-Tb 
          valangle(nangle) = phPSTb
        else if (namesite2.eq.'Cb') then !  P-(5')S-Cb
          valangle(nangle) = phPSCb
        else ! P-(5')S-Gb
          valangle(nangle) = phPSGb
        endif
      endif 
    else 
      call error ('angles', 'INCORRECT ORDER FOR SITES',faterr)
    endif      
    if ((i+3).le.nsites .and. strand1.eq.strand4) then
      nangle = nangle + 1
      if (nangle.gt.maxang) call error ('angles', 'The number of bond angles exceeds the maximum value', faterr)
      siteangle(nangle,1) = i
      siteangle(nangle,2) = i + 2 ! sugar central site
      siteangle(nangle,3) = i + 3
      angle(indexi(i,i+2)) = .true.
      angle(indexi(i+2,i+3)) = .true.
      angle(indexi(i,i+3)) = .true.
      if (namesite3.eq.'S ' .and. namesite4.eq.'P ') then ! P-(3')S(5')-P
        valangle(nangle) = phPSP    
      else
         call error ('angles', 'INCORRECT ORDER FOR SITES',faterr)
      endif
    endif            
  endif    
enddo

!OUTPUT
if (Qninfo) then
  write (outu,'(/6x,a)') 'DNA: INTRASTRAND BOND ANGLES'
  write (outu,'(6x,a)') '----------------------------'
  write (outu,'(6x,a)')'Nangle  strand  namesite1  namesite2  namesite3  natural bond angle'
  do i = 1, nangle
    i1 = siteangle(i,1)
    i2 = siteangle(i,2)
    i3 = siteangle(i,3)
    valueang = valangle(i)/cte
    write (outu,'(6x,i3,5x,i2,6x,a2,9x,a2,9x,a2,9x,f7.3)') i,strand(i1), namsite(i1), namsite(i2), namsite(i3),valueang
  enddo  
endif

return
end


