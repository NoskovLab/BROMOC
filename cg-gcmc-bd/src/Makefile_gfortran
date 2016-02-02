#! /usr/bin/make
#
# Main program and default target

obj    = bromoc
all : $(obj)

# System dependent compilation commands

# -align
#  #   para alinear common y otras cosas, aunque creo que asi no hace nada
# -C -traceback
#  #   para chequear array bounds, formatos, y mas; traceback lo situa con
#   mas precision
#   # -cm
#   elimina warnings de practicas cuestionables ... joder con promolden
#   #   esto es equivalente a -warn nousage
# -w
#  #   disable warning messages
# -warn = -warn all
#   enable all warnings
#   # -debug extended -debug inline_debug_info -debug-parameters all
#   produce informacion para el debug
#   # -g
#   introduce informacion para el debugger
#   # -inline-debug-info
#   mejora la posicion en el fuente de codigo inlined en debug/traceback
#   # -fltconsistency = -mieee-fp = -mp , tambien -mp1
#   codigo consistente para floating point (evita optimizaciones que van
#   #   ligadas a diferente precision); ralentiza, recomienda probar a ver
#   si es necesario, o a lo mejor vale con -mp1, intermedio
#   # -fpe0
#   floating-point exceptions enabled, underflow=0; el default es -fpe3
#   #   en el que estan disabled y que tiene el soporte IEEE completo
# -fp-model strict/precise/fast/fast=2
#  #   permitir optimizaciones agresivas de los calculos fp, default fast
# -ftrapuv
#  #   inicializa variables en el stack a cosas raras para detectar si se
#   usan antes de definir porque daran run-time errors
#   # -I<dir>
#   directorio para buscar INCLUDEs
# -L<dir>
#  #   directorio para buscar librerias (antes de los standard
# -l<string>
#  #   incluir la libreria lib<string>
# -inline-* , -ip* , -ipo*
#  #   varias opciones para controlar que rutinas se compilan "inlining"
# -mcmodel=small/medium/large
#  #   memory model small (2GB code&data, default) medium (2GB code) o
#   large (no restriction) para commons; allocation > 2GB se puede hacer
#   #   con cualquier modelo. small es el mas rapido
# -O=-O2
#  #   standard optimization
# -O3
#  #   como -O, pero con optimizaciones mas agresivas
# -opt-report -vec-report3
#  #   genera un report de la optimizacion y vectorizacion
# -tpp7
#  #   optimiza para Pentium 4 & more; hay muchas opciones de estas, pero
#   esta es la opcion por defecto, asi que la eliminamos.
#   #Opciones promolden
#---Next FFLAGC is to compile checking bounds.
FCOMPL = gfortran
FLAGPAR = -fopenmp
#FFLAGC = -C -fno-align-commons
#Double Precision Real
FFLAGC = -fno-align-commons -fdefault-real-8 
#Single Precision Real
#FFLAGC = -fno-align-commons
CCOMPL = gcc
CFLAGC =
opt    = -O3 -march=native -ffast-math -fno-protect-parens -ffree-line-length-none
#opt    = -Wall -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fdump-core -fcheck=all
#opt = -static -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal -fbacktrace
#opt = -Wall -static -ffpe-trap=invalid,zero,overflow -fbacktrace -g -fcheck=all
#opt  = -Wall -ffree-line-length-none -Wcharacter-truncation -Wconversion-extra -Wintrinsics-std -Wreal-q-constant -Wtabs -Wunderflow -Wunused-dummy-argument -Wunused-parameter -Walign-commons -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal -fbacktrace -g -fdump-core -fcheck=all

# Implicit rules

.SUFFIXES: .vec .f .c .s .o .a .lib .f90

.f90.o: 
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.s.o:
	as -o $@ $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(opt) -o $@ $<

.lib.a:
	make -f $(<:.lib=.Makefile) $@
	make -f $(<:.lib=.Makefile) clean

# Vectorization of selected files where available

#mdc*if vector
#.f.vec:
#	$(FCOMPL) -c $(FFLAGC) $(VFLAG) -o $@ $<
#	rm -f $*.o ; cp $@ $*.o
#
#VECTOR =
#
#$(VECTOR) :
#mdc*endif

# Object files

OBJECTS = extramod.o apfmod.o efpmod.o ioxmod.o constamod.o gsbpmod.o grandmod.o errormod.o\
	nucleotmod.o stdiomod.o strtoolsmod.o charfuncmod.o alpol2.o dalpol2.o splinemod.o sevalmod.o\
	charman.o cputime.o error.o fatnam.o discretize.o splinepot.o cosmphi2.o sinmphi2.o\
	getopt.o ioinit.o lualloc.o lunalloc.o m3.o dm3.o rpowerl2.o rect_norm.o sphe_norm.o\
	sphe_svpol.o cross_product.o indexi.o factori.o angles.o bonds.o count.o dihedral.o countions.o\
	dynamics0nuc.o find.o repwalls.o fixcoor.o go_qq.o grand.o header.o mapsplot.o help.o insert.o\
	interact.o membrane.o metropolis.o move.o native_structure.o rfpar.o feedseed.o random.o readmmij.o\
	chrden.o readphi.o rect_rf0.o rect_rf1.o rgauss.o staefield.o shell_simul.o simul1.o simul.o sphe_rf0.o\
	sphe_rf1.o spline.o staticf0.o staticf1.o switch1.o switch2.o switch3.o testfirst.o vdwgd0spln.o\
	vdwgd0trln.o vdwgd1spln.o vdwgd1trln.o wrttraj.o angbond.o angdih.o reasig.o polynom.o nindex.o

OBJECTSMONO = dynamics0.mono.o dynamics1.mono.o dynamics2.mono.o dynamics3.mono.o energy.mono.o

OBJECTSPARA = dynamics0.para.o dynamics1.para.o dynamics2.para.o dynamics3.para.o energy.para.o

FPARAFILES = dynamics0.para.f90 dynamics1.para.f90 dynamics2.para.f90 dynamics3.para.f90 energy.para.f90

#IFILES = charli.fcm consta.fcm gsbp.fcm grand.fcm error.inc\
#         iox.inc nucleot.inc stdio.inc strtools.inc strtools.ini
 
# Clean objects

clean:
	@rm -f core $(OBJECTS) $(OBJECTSMONO) $(OBJECTSPARA) $(obj) $(obj)-para *.mod

# Main program target

$(obj):	$(OBJECTS) $(OBJECTSMONO)
	$(FCOMPL) $(FFLAGC) $(opt) -o $(obj) $(OBJECTSMONO) $(OBJECTS) $(ldflag)
	$(FCOMPL) $(FFLAGC) $(opt) $(FLAGPAR) -o $(obj)-para $(OBJECTS) $(FPARAFILES)

install:
	@mv $(obj) $(obj)-para ../bin/
	@if [ ! -e ../../bin/$(obj) ] ; then ln -s ../cg-gcmc-bd/bin/$(obj) ../../bin/$(obj); fi
	@if [ ! -e ../../bin/$(obj)-para ] ; then ln -s ../cg-gcmc-bd/bin/$(obj)-para ../../bin/$(obj)-para; fi
 
