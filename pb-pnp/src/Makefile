#! /usr/bin/make
#
# Main program and default target

obj    = pb-pnp
all : $(obj)

# System dependent compilation commands

FCOMPL = gfortran
FFLAGC = -fPIC -mcmodel=large -O3 -march=native -ffast-math -fno-protect-parens 
#FFLAGC = -C 
CCOMPL = gcc
CFLAGC =
#opt    = -ffast-math -funroll-loops -ftree-vectorize -O3 -ffree-line-length-none 
#opt    = -ffree-line-length-none 
SYSDEP =
UTILS  = 
ldflag = 

# Implicit rules

.SUFFIXES:
.SUFFIXES: .vec .f .c .s .o .a .lib .f90

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.f90.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.s.o:
	as -o $@ $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(opt) -o $@ $<

.lib.a:
	make -f $(<:.lib=.Makefile) $@
	make -f $(<:.lib=.Makefile) clean

# Object files

OBJECTS = charfuncmod.o charman.o cputime.o mayer.o misc.o pbeqsolver.o pnpanal.o pnpinit.o pnpio.o pnpmain-gnu.o pnpsolver.o prep.o matrix.o mmij.o rfpar.o

IFILES = consta.fcm mainio.fcm misc.fcm mmij.fcm pnp.fcm rfintsph.fcm
# Clean objects

clean:
	@rm -f $(OBJECTS) $(obj) *.mod

# Main program target

$(obj): $(SYSDEP) $(VECTOR) $(OBJECTS) $(FFILES) $(UTILS)
	$(FCOMPL) $(FFLAGC) $(ldflag) -o $(obj) $(OBJECTS) $(SYSDEP) $(UTILS)

install:
	@mv $(obj) ../bin/
	@ln -s ../pb-pnp/bin/$(obj) ../../bin/$(obj)

