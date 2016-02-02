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

subroutine reasig() 
use grandmod
use nucleotmod

implicit none
integer     i, itype
integer     typAb, typTb, typCb, typGb
logical*1     ok(4),okt

itype = 0
if (Qnucl) then
  itype = 2 ! itype = 1 'S'; itype = 2 'P'
  atnam2(1) = 'S'
  atnam2(2) = 'P'
  i = 0
  ok=.false.
  okt=.false.
  do while (.not.okt.and.i.lt.nsites)
    i = i + 1
    if (namsite(i).eq.'Ab') ok(1) = .true.
    if (namsite(i).eq.'Tb') ok(2) = .true.
    if (namsite(i).eq.'Cb') ok(3) = .true.
    if (namsite(i).eq.'Gb') ok(4) = .true.
    okt=ok(1).and.ok(2).and.ok(3).and.ok(4)
  enddo
  if (ok(1)) then
    itype = itype + 1
    atnam2(itype) = 'Ab'
    typAb = itype
  endif
  if (ok(2)) then
    itype = itype + 1
    atnam2(itype) = 'Tb'
    typTb = itype
  endif
  if (ok(3)) then
    itype = itype + 1
    atnam2(itype) = 'Cb'
    typCb = itype
  endif
  if (ok(4)) then
    itype = itype + 1
    atnam2(itype) = 'Gb'
    typGb = itype
  endif
  do i = 1, nsites
!    typat(i) = i
    if (namsite(i).eq.'Ab') then
      typtyp(i) = typAb
    else if (namsite(i).eq.'Tb') then
      typtyp(i) = typTb
    else if (namsite(i).eq.'Cb') then
      typtyp(i) = typCb
    else if (namsite(i).eq.'Gb') then
      typtyp(i) = typGb
    else if (namsite(i).eq.'S ') then
      typtyp(i) = 1
    else
      typtyp(i) = 2
    endif
  enddo
endif ! Qnucl
ndna=itype
if (Qpar) then
  do i = nold+1, ntype
    itype = itype + 1
    nwtype(i) = itype
    atnam2(itype) = atnam(i)
  enddo
endif ! Qpar
nttyp=itype
nion=nttyp-ndna
ndnaxnion=ndna*nion
return
end
