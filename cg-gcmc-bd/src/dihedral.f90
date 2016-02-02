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

subroutine dihedral 
! Interaction sites are stored in the order B1-S1-P1--B2-S2-P2--B3-S3-P3..., 
! where P=Phosphate, S=Sugar and B=Base. Note that the bond between a phospate
! and a sugar belonging to the same nucleotide is S(5')-P, while the
! bond S(3')-P joins together neighboring residues.
use constamod
use stdiomod
use errormod
use grandmod
use nucleotmod
implicit none 
!     Local variables     
real  cte, diheang 
integer i, i1, i2, i3, i4, strand1, strand2, strand3, strand4, strand5, strand6
character*2 namesite1, namesite2, namesite3, namesite4, namesite5, namesite6   

cte = pi/180.0
ndihe = 0 ! Initialization

do i = 1, nsites-3 
  namesite1 = namsite(i)
  namesite2 = namsite(i+1) 
  namesite3 = namsite(i+2)
  namesite4 = namsite(i+3)
  strand1 = strand(i)
  strand2 = strand(i+1)
  strand3 = strand(i+2)
  strand4 = strand(i+3)
  if ((i+4).le.nsites) then 
    namesite5 = namsite(i+4)
    strand5 = strand(i+4)
  endif  
  if ((i+5).le.nsites) then
    namesite6 = namsite(i+5)
    strand6 = strand(i+5)
  endif    
  if (namesite1.eq.'Ab'.or.namesite1.eq.'Tb'.or.namesite1.eq.'Cb'.or.namesite1.eq.'Gb') then
    if ((i+4).le.nsites .and. strand1.eq.strand5) then 
      ndihe = ndihe + 1      
      if (ndihe.gt.maxdihe) call error ('dihedral', 'The number of bond angles exceeds the maximum value', faterr)
      sitedihe(ndihe,1) = i
      sitedihe(ndihe,2) = i + 1  
      sitedihe(ndihe,3) = i + 2
      sitedihe(ndihe,4) = i + 4
      if (namesite2.eq.'S '.and.namesite3.eq.'P '.and.namesite5.eq.'S ') then
        if (Qinvstr) then
          if (namesite1.eq.'Ab') then ! Ab-S(5')-P-(3')S
            valdihe(ndihe) = dhAbSPS
          else if (namesite1.eq.'Tb') then ! Tb-S(5')-P-(3')S
            valdihe(ndihe) = dhTbSPS
          else if (namesite1.eq.'Cb') then ! Cb-S(5')-P-(3')S      
            valdihe(ndihe) = dhCbSPS
          else ! Gb-S(5')-P-(3')S
            valdihe(ndihe) = dhGbSPS      
          endif
        else
          if (namesite1.eq.'Ab') then ! Ab-S(3')-P-(5')S
            valdihe(ndihe) = dhSPSAb
          else if (namesite1.eq.'Tb') then ! Tb-S(3')-P-(5')S
            valdihe(ndihe) = dhSPSTb
          else if (namesite1.eq.'Cb') then ! Cb-S(3')-P-(5')S      
            valdihe(ndihe) = dhSPSCb
          else ! Gb-S(3')-P-(5')S
            valdihe(ndihe) = dhSPSGb
          endif
        endif 
      else 
        call error ('dihedral', 'INCORRECT ORDER FOR SITES',faterr) 
      endif 
    endif           
  else if (namesite1.eq.'S ') then
    if (strand1.eq.strand4) then
      if (namesite2.eq.'P ' .and. namesite3.eq.'Ab'.or.namesite3.eq.'Tb'.or.namesite3.eq.'Cb'.or. &
          namesite3.eq.'Gb' .and. namesite4.eq.'S ') then  ! S(5')-P-(3')S-B   
        ndihe = ndihe + 1
        if (ndihe.gt.maxdihe) call error ('dihedral', 'The number of bond angles exceeds the maximum value', faterr)
        sitedihe(ndihe,1) = i 
        sitedihe(ndihe,2) = i + 1 
        sitedihe(ndihe,3) = i + 3
        sitedihe(ndihe,4) = i + 2
        if (Qinvstr) then  
          if (namesite3.eq.'Ab') then ! S(5')-P-(3')S-Ab
            valdihe(ndihe) = dhSPSAb
          else if (namesite3.eq.'Tb') then ! S(5')-P-(3')S-Tb
            valdihe(ndihe) = dhSPSTb
          else if (namesite3.eq.'Cb') then ! S(5')-P-(3')S-Cb    
            valdihe(ndihe) = dhSPSCb
          else ! S(5')-P-(3')S-Gb
            valdihe(ndihe) = dhSPSGb
          endif     
        else
          if (namesite3.eq.'Ab') then ! S(3')-P-(5')S-Ab
            valdihe(ndihe) = dhAbSPS
          else if (namesite3.eq.'Tb') then ! S(3')-P-(5')S-Tb
            valdihe(ndihe) = dhTbSPS
          else if (namesite3.eq.'Cb') then ! S(3')-P-(5')S-Cb
            valdihe(ndihe) = dhCbSPS
          else ! S(3')-P-(5')S-Gb
            valdihe(ndihe) = dhGbSPS
          endif
        endif
      else  
        call error ('dihedral', 'INCORRECT ORDER FOR SITES',faterr)
      endif          
    endif 
    if ((i+4).le.nsites .and. strand1.eq.strand5) then    
      if (namesite2.eq.'P '.and.namesite4.eq.'S '.and.namesite5.eq.'P ') then  ! S(5')-P-(3')S(5')-P
        ndihe = ndihe + 1
        if (ndihe.gt.maxdihe) call error ('dihedral', 'The number of bond angles exceeds the maximum value', faterr)
        sitedihe(ndihe,1) = i
        sitedihe(ndihe,2) = i + 1
        sitedihe(ndihe,3) = i + 3
        sitedihe(ndihe,4) = i + 4
        if (Qinvstr) then
          valdihe(ndihe) = dhSPSP
        else
          valdihe(ndihe) = dhPSPS
        endif
      else
        call error ('dihedral', 'INCORRECT ORDER FOR SITES', faterr)
      endif  
    endif  
  else ! namesite1.eq.'P '
    if ((i+5).le.nsites .and. strand1.eq.strand6) then
      if (namesite3.eq.'S '.and.namesite4.eq.'P '.and.namesite6.eq.'S ') then ! P-(3')S(5')-P-(3')S     
        ndihe = ndihe + 1
        if (ndihe.gt.maxdihe) call error ('dihedral', 'The number of bond angles exceeds the maximum value', faterr)
        sitedihe(ndihe,1) = i
        sitedihe(ndihe,2) = i + 2
        sitedihe(ndihe,3) = i + 3
        sitedihe(ndihe,4) = i + 5
        if (Qinvstr) then
          valdihe(ndihe) = dhPSPS
        else
          valdihe(ndihe) = dhSPSP
        endif 
      else
        call error ('dihedral', 'INCORRECT ORDER FOR SITES',faterr)
      endif 
    endif
  endif  
enddo

!     OUPUT
if (Qninfo) then
  write (outu,'(/6x,a)') 'DNA: INTRASTRAND DIHEDRAL ANGLES'
  write (outu,'(6x,a)') '--------------------------------'
  write (outu,'(6x,a)') 'Ndihe  strand  namesite1  namesite2  namesite3  namesite4  natural dihedral angle'
  do i = 1, ndihe
    i1 = sitedihe(i,1)
    i2 = sitedihe(i,2)
    i3 = sitedihe(i,3)
    i4 = sitedihe(i,4)
    diheang = valdihe(i)/cte
    write (outu,'(6x,i3,4x,i2,6x,a2,9x,a2,9x,a2,9x,a2,9x,f8.3)') i,strand(i1), namsite(i1), namsite(i2), namsite(i3), namsite(i4), diheang  
  enddo
endif

return
end


