!    IMC-MACRO - Inverse Monte Carlo for Macromolecules
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

program convertrst
implicit none
integer*4 nav,nmksf,nmks0,imku,npot,imu,iud,nop,i,j,kode
real*4 eners,avexp,ave2,virs,vire
real*4,allocatable :: x(:),y(:),z(:),cors(:),corp(:,:)
integer*4,allocatable :: corsi(:),corpi(:,:)
character fdmp*256,fdmpw*256,label*12,label2*16

write(*,*) 'Convert from IMC v3.1 to IMC v3.2 Restart File'
write(*,'(A$)')'Restart file read from: '
read(*,*) fdmp

write(*,'(A$)')'Number of particles: '
read(*,*) nop

write(*,'(A$)')'Restart file write to: '
read(*,*) fdmpw

allocate (x(nop),y(nop),z(nop))

open(unit=45,file=fdmp,status='old',form='unformatted',IOSTAT=kode)
if (kode.ne.0) then
  write(*,*) 'File ',fdmp,' not found '
  stop
endif
read(45)label
read(45) nav,nmksf,nmks0,imku,npot,imu,iud,eners,avexp,ave2,virs,vire
allocate (cors(npot),corp(npot,npot),corsi(npot),corpi(npot,npot))
read(45) (x(i),y(i),z(i),i=1,nop)
read(45) (cors(i),i=1,npot)
do i=1,npot
  read(45)(corp(i,j),j=i,npot)
enddo
close(45)

do i=1,npot
  corsi(i)=int(cors(i))
  do j=1,npot
    corpi(i,j)=int(corp(i,j))
  enddo
enddo
label2=label
open(unit=45,file=fdmpw,status='unknown',form='unformatted')
write(45)label2
write(45)nav,nmksf,nmks0,nop,npot,imu,iud,eners,avexp,ave2,virs,vire
write(45)(x(i),y(i),z(i),i=1,nop)
write(45)(corsi(i),i=1,npot)
do i=1,npot
  write(45)(corpi(i,j),j=i,npot)
enddo
close(45)

end program
