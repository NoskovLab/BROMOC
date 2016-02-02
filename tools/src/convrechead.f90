!    CONVRECHEAD converts 8-bytes record header of unformatted fortran files to 4-bytes'
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

program convert_record_head
implicit none
integer*1 oc
integer*4 kode,i,j,k,record4,head4
integer*8 head8,posi,record8
character*128 input,output
character bs*20,line
integer*1,allocatable :: bigblock(:)

do i=1,len(bs)
  bs(i:i)=achar(8)
enddo
write(*,'(A)') 'CONVRECHEAD converts 8-bytes record header of unformatted fortran files to 4-bytes'
write(*,*)
write(*,'(A$)') 'Input filename: '
read(*,*) input
write(*,*)

open (unit=1,file=input,form='unformatted',access='stream',iostat=kode,action='read')
! Check if header=8
read(1,pos=1,iostat=kode) head8
read(1,pos=1+8+head8,iostat=kode) record8
if (kode.eq.0.and.record8.eq.head8) then
  write(*,'(A)') '8-bytes record header detected'
else 
  write(*,'(A)') '8-bytes record header not detected'
endif

! Check if header=4
read(1,pos=1,iostat=kode) head4
read(1,pos=1+4+head4,iostat=kode) record4
if (kode.eq.0.and.record4.eq.head4) then
  write(*,'(A)') '4-bytes record header detected'
  oc=2
else
  write(*,'(A)') '8-bytes record header not detected'
  oc=1
endif

rewind(1)

write(*,*)
write(*,'(A)') '1- Converts 8-bytes record header of unformatted fortran files to 4-bytes'
write(*,'(A)') '2- Converts 4-bytes record header of unformatted fortran files to 8-bytes'
write(*,'(A,I0,A$)') 'ENTER for detected option (',oc,'): '
read (*,'(A)') line
if (len_trim(line).gt.0) read(line,*) oc
write(*,'(A,I0,A/)') 'Option ',oc,' Selected'

write(*,'(A$)') 'Output filename: '
read(*,*) output
write(*,*)

open (unit=2,file=output,form='unformatted',access='stream',action='write')
posi=0
write(*,'(A,I20$)') 'Reading Offset: ',posi
if (oc.eq.1) then 
  read(1,iostat=kode) head8
  do while (kode.eq.0)
    record4=head8
    write(2) record4
    allocate (bigblock(record4))
    read(1) bigblock
    write(2) bigblock
    deallocate (bigblock)
    read(1) head8
    record4=head8
    write(2) record4
    inquire(unit=1,pos=posi)
    write(*,'(A20,I20$)') bs,posi
    read(1,iostat=kode) head8
  enddo
else
  read(1,iostat=kode) head4
  do while (kode.eq.0)
    record8=head4
    write(2) record8
    allocate (bigblock(record8))
    read(1) bigblock
    write(2) bigblock
    deallocate (bigblock)
    read(1) head4
    record8=head4
    write(2) record8
    inquire(unit=1,pos=posi)
    write(*,'(A20,I20$)') bs,posi
    read(1,iostat=kode) head4
  enddo
endif
close(1)
close(2)
write(*,'(//A/)') 'Normal Termination'
end program

