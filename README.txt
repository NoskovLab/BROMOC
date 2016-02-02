BROMOC Suite  (2011-2014)
************
---------------------------------------------------------------
If you use this program please site us:

BROMOC Suite: Monte Carlo/Brownian dynamics suite for studies of ion permeation and DNA transport in biological and artificial pores with effective potentials. 
Pablo M. De Biase*, Suren Markosyan, and Sergei Noskov* 
J. Comput. Chem., 2014, published online: 15 DEC 2014, DOI: 10.1002/jcc.23799

BROMOC-D: Brownian Dynamics/Monte-Carlo Program Suite to Study Ion and DNA Permeation in Nanopores
Pablo M. De Biase, Carlos J. F. Solano, Suren Markosyan, Luke Czapla, and Sergei Yu. Noskov*
J Chem Theory Comput. Jul 10, 2012; 8(7): 2540-2551.
Published online May 24, 2012. doi:  10.1021/ct3004244

Microsecond simulations of DNA and ion transport in nanopores with novel ion-ion and ion-nucleotides effective potentials
Pablo M. De Biase*, Suren Markosyan and Sergei Noskov*
J Comput Chem. 2014 Apr 5;35(9):711-721.
Article first published online: 12 FEB 2014
DOI: 10.1002/jcc.23544

---------------------------------------------------------------

Author: Pablo M. De Biase
e-mail: pablodebiase@gmail.com

add to .bashrc
export PATH=~/bromoc/bin:~/bromoc/tools/scripts:$PATH

ALL
===
To install all:
inside BROMOC folder

#For GNU Compilers
$ ./install

   --or--

#For INTEL Compilers
$ ./install intel


CG-GCMC-BD
==========

$ cd cg-gcmc-bd/src

#For GNU Compilers
$ make clean -f Makefile_gfortran
$ make -f Makefile_gfortran           
$ make install -f Makefile_gfortran
   
   --or--

#For INTEL Compilers

$ make clean -f Makefile_intel
$ make -f Makefile_intel
$ make install -f Makefile_intel


IMC-MACRO
=========

To compile:

$ cd imc-macro/src

#For GNU Compilers
$ make clean -f Makefile_gfortran
$ make -f Makefile_gfortran          
$ make install -f Makefile_gfortran
$ cd ../scripts
$ ./mklnks

   --or--

#For INTEL Compilers

$ make clean -f Makefile_intel
$ make -f Makefile_intel
$ make install -f Makefile_intel
$ cd ../scripts
$ ./mklnks


PB-PNP
======

To compile:

$ cd pb-pnp/src
#For GNU Compilers
$ make clean -f Makefile_gfortran
$ make -f Makefile_gfortran           
$ make install -f Makefile_gfortran
   
   --or--

#For INTEL Compilers

$ make clean -f Makefile_intel
$ make -f Makefile_intel
$ make install -f Makefile_intel


TOOLS
=====

To compile:

$ cd tools/src

#For GNU Compilers
$ ./makeall
   
   --or--

#For INTEL Compilers
$ ./makeall intel   


