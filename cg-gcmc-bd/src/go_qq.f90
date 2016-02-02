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

subroutine go_qq 
! Interaction sites are stored in the order B1-S1-P1--B2-S2-P2--B3-S3-P3..., 
! where P=Phosphate, S=Sugar and B=Base. Note that the bond between a phospate
! and a sugar belonging to the same nucleotide is S(5')-P, while the
! bond S(3')-P joins together neighboring residues. 
use stdiomod
use errormod 
use grandmod
use nucleotmod
implicit none
!Local variables    
real tol
parameter (tol=1.0e-4)
integer i, j, i1, i2 
real  dij
logical*1 logvab1, logvab2
!Function to generate the index for a supervector
integer indexi

!Initializations      
nstack = 0
nbp = 0
nex = 0
nqq = 0
nsolv = 0
 
do i = 2, nsites
  do j = 1, i-1
!    Two sites are excluded form all non-bonded interactions if they
!    constitute a bonb
    if (.not.bond(indexi(i,j))) then
!      Native contact's length        
      dij = sqrt((xnat(j)-xnat(i))**2 + (ynat(j)-ynat(i))**2 + (znat(j)-znat(i))**2)
!      Interaction between nucleotides that form a native contact in
!      the target structure
      if (strand(i).eq.strand(j) .and. (dij-9.0).lt.0.0 .and. abs(dij-9.0).ge.tol) then
        nstack = nstack + 1
        if (nstack.gt.maxstack) then
          call error ('go_qq', 'The number of native contacts exceeds the maximum value', faterr)
        endif   
        sitestack(nstack,1) = i
        sitestack(nstack,2) = j
        sgstack(nstack) = dij
!      Hydrogen bonding between any complementary base pair      
      else if (namsite(i).eq.'Gb'.and.namsite(j).eq.'Cb' .or. &
               namsite(i).eq.'Cb'.and.namsite(j).eq.'Gb' .or. &
               namsite(i).eq.'Ab'.and.namsite(j).eq.'Tb' .or. &
               namsite(i).eq.'Tb'.and.namsite(j).eq.'Ab') then
        nbp = nbp + 1
        if (nbp.gt.maxbp) then
          call error ('go_qq', 'The number of hydrogen bondings exceeds the maximum value', faterr)
        endif         
        sitebp(nbp,1) = i
        sitebp(nbp,2) = j
        if (namsite(i).eq.'Gb' .or. namsite(i).eq.'Cb') then ! Gb-Cb
          sgbp(nbp) = 2.8694
        else ! Ab-Tb
          sgbp(nbp) = 2.9002
        endif 
!     Excluded volume interactions (repulsive core potential)
      else  
        nex = nex + 1
        if (nex.gt.maxex) then   
          call error ('go_qq', 'The number of excluded volume terms exceeds the maximum value', faterr)
        endif 
        siteex(nex,1) = i
        siteex(nex,2) = j
        logvab1 = namsite(i).eq.'P ' .or. namsite(i).eq.'S ' 
        logvab2 = namsite(j).eq.'P ' .or. namsite(j).eq.'S ' 
        if (.not.logvab1 .and. .not.logvab2) then ! mismatched bases
          sgex(nex) = 1.0   
        else ! otherwise          
          sgex(nex) = 6.86
        endif
      endif   
    !  Coulomb interaction between interaction sites
      if (.not.angle(indexi(i,j)) .and. namsite(i).eq.'P '.and.namsite(j).eq.'P ') then
        nqq = nqq + 1    
        if (nqq.gt.maxqq) then
          call error ('go_qq', 'The number of coulomb interactions terms exceeds the maximum value', faterr)
        endif 
        siteqq(nqq,1) = i
        siteqq(nqq,2) = j
      endif     
    endif     
   ! Solvent-induced contribution
    if (Qsolv) then
      if (strand(i).ne.strand(j) .and. namsite(i).eq.'S '.and.namsite(j).eq.'S ') then
        nsolv = nsolv + 1
        if (nsolv.gt.maxsolv) then
          call error ('go_qq', 'The number of Solvent-induced contributions terms exceeds the maximum value', faterr)
        endif
        siteslv(nsolv,1) = i
        siteslv(nsolv,2) = j
      endif 
    endif
  enddo  
enddo 

!OUTPUT
if (Qninfo) then
  write (outu,'(/6x,a)') 'DNA: INTRASTRAND NATIVE CONTACTS'
  write (outu,'(6x,a)') '--------------------------------'
  write (outu,'(6x,a)') 'Nst  strand  namesite1  namesite2  distance'  
  do i = 1, nstack
    i1 = sitestack(i,1)
    i2 = sitestack(i,2)
    write (outu,'(3x,i5,4x,i2,6x,a2,9x,a2,6x,f8.3)') i, strand(i1), namsite(i1), namsite(i2), sgstack(i)
  enddo
  write (outu,'(/6x,a)') 'DNA: HYDROGEN BONDING'
  write (outu,'(6x,a)') '---------------------'      
  write (outu,'(6x,a)') 'Nbp  strand1  strand2  namesite1  namesite2  LJ parameter (distance)'
  do i = 1, nbp
    i1 = sitebp(i,1)
    i2 = sitebp(i,2)
    write (outu,'(3x,i5,4x,i2,7x,i2,7x,a2,9x,a2,9x,f8.3)') i, strand(i1), strand(i2), namsite(i1), namsite(i2), sgbp(i)
  enddo 
  write (outu,'(/6x,a)') 'DNA: EXCLUDE VOLUME'
  write (outu,'(6x,a)') '-------------------'
  write (outu,'(6x,a)') 'Nex  strand1  strand2  namesite1i  namesite2  distance'    
  do i = 1, nex
    i1 = siteex(i,1)
    i2 = siteex(i,2) 
    write (outu,'(3x,i5,4x,i2,7x,i2,7x,a2,9x,a2,6x,f8.3)') i,strand(i1), strand(i2), namsite(i1), namsite(i2), sgex(i)
  enddo        
  write (outu,'(/6x,a)') 'DNA: COULOMB INTERACTIONS'
  write (outu,'(6x,a)') '-------------------------'
  write (outu,'(6x,a)') 'Nqq  strand1  strand2  namenucl1  namenucl2  namesite1 namesite2'   
  do i = 1, nqq
    i1 = siteqq(i,1)
    i2 = siteqq(i,2)      
    write (outu,'(3x,i5,4x,i2,7x,i2,7x,a2,9x,a2,9x,a2,9x,a2)') i, strand(i1), strand(i2), namnucl(i1), namnucl(i2), namsite(i1), namsite(i2)  
  enddo 
  if (Qsolv) then
    write (outu,'(/6x,a)') 'DNA: SOLVENT-INDUCED CONTRIBUTIONS'
    write (outu,'(6x,a)') '----------------------------------'
    write (outu,'(6x,a)') 'Nsolv strand1  strand2  namenucl1  namenucl2  namesite1 namesite2'
    do i = 1, nsolv
      i1 = siteslv(i,1) 
      i2 = siteslv(i,2)
      write (outu,'(3x,i5,4x,i2,7x,i2,7x,a2,9x,a2,9x,a2,9x,a2)') i,strand(i1), strand(i2), namnucl(i1), namnucl(i2), namsite(i1), namsite(i2)
    enddo   
  endif 
endif 

return
end
