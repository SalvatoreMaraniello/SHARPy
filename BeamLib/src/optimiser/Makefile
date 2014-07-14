#
# ---------------- Makefile for opt_xeam: S.M. 08/07/2014 ----------------------
#-------------------------------------------------------------------------------
#
# Reference: 
# a. ../install/Makefile (HH and RP)
# b. ../../wrapper/Makefile (ADR)
# c. http://oreilly.com/catalog/make3/book/
# 
# Prerequisites:
# - check out the path to llapack and lblas libraries (LIBSDIR variable);
# - Change TESTCASE to consider different test problems;
#	- the input_xxx.f90 for the forward code compilation is under the fortran folder
#	- the input_xxx.f90 for the optimiser is under ./main
#
# Usage: 
# - 'make all'    : see 'make exe'
# - 'make exe'    : to create the forward solver exe from sources (good for debugging)
# - 'make opt'    : to compile the optimiser 
# - 'make clean'  : clean up fwd and optimiser
# - 'make run'    : runs the optimisation
# - 'make opttest': cleans all files (excutable included), recompiles them and 
#                   runs the optimiser 
# - 'make rmout'  : removes output from simulation
#
# Unused:
# - 'make xbeamso' : from makefile from HH and RP (not working)
# - 'make xbeamlib': from makefile from HH and RP (not working)
#
# for testing...
# - 'make objlist'
# - 'make dirprint'
# - 'make 'dirlist'
# - 'make libtest': compiles libraries for forward problem
# - 'make xbeamtest': compiles libraries first and xbeam solver and assembly
#
# Developer Notes:
# - The current makefile does not include all the compile options of:
# 	SHARPy/BeamLib/src/fortran/install/Makefile (HH & RP)
# - Commented code has been left at the end of the file for future development
# 
#-------------------------------------------------------------------------------

all: xbeam
	
exe: xbeam

clean:
	rm -f $(OBJOPT) $(OBJFWD) $(OBJWRAP) $(OBJDIR)/*.mod $(OBJDIR)/*.a $(OPTOBJDIR)/*.mod $(OPTOBJDIR)/*.mod

opt: xbeamopt

run: 
	./bin/xbeamopt.exe
	
opttest:
	rm -f ./bin/xbeamopt.exe 
	make clean
	make opt
	make run
	
rmout:
	rm -f ./*_def.txt ./*und.txt ./*.mrb 
	
################################################################################
##  COMPILER OPTIONS.
# 
# F90OPTS: used to compile forward code ref.[../install/Makefile (HH and RP)]
# F90FLAG: to ocmpile object for dynamic library ref.[../../wrapper/Makefile A]
# LIBFLAG used for dynamic library ref.[../../wrapper/Makefile A]
#
# Options List:
# -o: specify output name
# -c: tells to compile an object file
# -J: mod files destination
# -std: 'f95',f2003', f2008',gnu', or legacy' to force compiling with a 
#
# notes;
# gnu: default; works
# f95, f2003: low level error (interface_lapack_v3.f90

F90COMP= gfortran
F90OPTS= -std=gnu -g -m64 -cpp -J $(OBJDIR)
F90FLAG=-fPIC -Wextra -cpp -DNOLAPACK -g
LIBFLAG=-fPIC -shared -ffree-form -g
LIBS= -llapack -lblas
XBEAMLIB= -lxbeam


################################################################################
## I/O SELECTION.
#
# Notes for Developer:
# - when inputting folder names, avoid spaces at the end of any assignment

TESTCASE=sal

# ROOTDIR = SHARPy/BeamLib/src
ROOTDIR=..

# --------------------------------------------------- Forward Solver Directories 
# MAINDIR = SHARPy/BeamLib/src/fortran/main  (input for fwd problem)
# XBEAMDIR  = SHARPy/BeamLib/src/fortran/install/src (fwd problem solver and assembly)
# LIBDIR=SHARPy/BeamLib/src/lib/src (auxiliary libraries)
# OBJDIR  = SHARPy/BeamLib/src/fortran/install/linux/obj (all object for fwd problem)
MAINDIR=$(ROOTDIR)/fortran/main
XBEAMDIR=$(ROOTDIR)/fortran/install/src
LIBDIR=$(ROOTDIR)/lib/src
OBJDIR=$(ROOTDIR)/fortran/install/linux/obj

# -------------------------------------------------------- Optimiser Directories
# OPTDIR = SHARPy/BeamLib/src/optimiser
# OPTMAINDIR (contains input and main files for the optimiser)
# OPTSRCDIR  (contains main source code for optimiser)
# OPTLIBSDIR (contains auxiliary libraries for optimisation)
# OPTOBJDIR  (drops objects for the optimiser)
# OPTBINDIR  (executable)
OPTDIR=.
OPTMAINDIR=$(OPTDIR)/main
OPTSRCDIR=$(OPTDIR)/src
OPTLIBDIR=$(OPTDIR)/lib
OPTOBJDIR=$(OPTDIR)/obj
OPTBINDIR=$(OPTDIR)/bin

# ----------------------------------------------- Wrapper and external Libraries
WRAPDIR=$(ROOTDIR)/wrapper
LIBSDIR=-L/opt/lib


################################################################################
## Source lists.
#
# This is required only to specify the prerequisites for the make targets but
# not to compile each single module

# -------------------------------------------------------------- Forward Problem
OBJLIB=\
     $(OBJDIR)/lib_sparse.o  \
     $(OBJDIR)/interface_lapack_v3.o \
     $(OBJDIR)/lib_out.o     \
     $(OBJDIR)/lib_lu.o      \
     $(OBJDIR)/lib_fem.o     \
     $(OBJDIR)/lib_rot.o     \
     $(OBJDIR)/lib_rotvect.o \
     $(OBJDIR)/lib_bgeom.o   \
     $(OBJDIR)/lib_cbeam3.o  \
     $(OBJDIR)/lib_xbeam.o	 \
     $(OBJDIR)/lib_perturb.o
OBJXBEAM=\
	 $(OBJDIR)/xbeam_shared.o	\
     $(OBJDIR)/xbeam_undef.o	\
     $(OBJDIR)/cbeam3_asbly.o	\
     $(OBJDIR)/cbeam3_solv.o	\
	 $(OBJDIR)/xbeam_asbly.o	\
	 $(OBJDIR)/xbeam_solv.o		\
	 $(OBJDIR)/xbeam_perturb.o	\
	 $(OBJDIR)/xbeam_fdiff.o	 
OBJMAIN=\
     $(OBJDIR)/input_$(TESTCASE).o \
     $(OBJDIR)/main_$(TESTCASE).o

# these have to be sequential, i.e. the order of compilation must follow the 
# code hierarchy 
OBJFWD= $(OBJLIB) $(OBJXBEAM) $(OBJMAIN) 

# ---------------------------------------------------------------------- Wrapper
OBJWRAP=\
     $(OBJDIR)/input_$(TESTCASE).o	\
	 $(OBJDIR)/xbeam_input_wrap.o	\
	 $(OBJDIR)/cbeam3_static_wrap.o

# -------------------------------------------------------------------- Optimiser
OBJOPTLIB=\
    $(OPTOBJDIR)/lib_perf.o

OBJOPTSRC=\
	$(OPTOBJDIR)/opt_fd_perturb.o
    
OBJOPTMAIN=\
    $(OPTOBJDIR)/input_$(TESTCASE).o     \
    $(OPTOBJDIR)/fwd_main.o              \
    $(OPTOBJDIR)/opt_input_$(TESTCASE).o \
    $(OPTOBJDIR)/opt_main_$(TESTCASE).o  


OBJOPT = $(OBJLIB) $(OBJXBEAM) $(OBJOPTLIB) $(OBJOPTSRC) $(OBJOPTMAIN)


################################################################################
## Define make targets
#
# Automatic variables:
# $@: target filename
# $% filename element of an archive member specification
# $< filename of the first prerequisite
# $? names of all prerequisites that are newer than the target, separated by spaces. 

#--------------------------------------------------------------------- Optimiser
xbeamopt: $(OBJOPT)
	$(F90COMP) $(F90OPTS) -o $(OPTBINDIR)/xbeamopt.exe $(OBJOPT) $(LIBS) $(LIBSDIR)	
	
#--------------------------------------------------------------- Forward Problem
xbeam: $(OBJFWD) 
	$(F90COMP) $(F90OPTS) -o $(ROOTDIR)/fortran/bin/xbeam.exe $(OBJFWD) $(LIBS) $(LIBSDIR)

xbeamlib: $(OBJLIB) $(OBJXBEAM) $(OBJWRAP) 
	ar -ru $(ROOTDIR)/fortran/install/linux/obj/libxbeam.a $(OBJLIB) $(OBJXBEAM) $(OBJWRAP)

xbeamso: $(OBJLIB) $(OBJXBEAM) $(OBJWRAP)
	$(F90COMP) -shared -o libadd.so $(OBJLIB) $(OBJXBEAM) $(OBJWRAP) $(LIBSDIR) $(LIBS)

#----------------------------------------------------------------------- testing

optlibtest: $(OBJOPTLIB)

fwdlibtest: $(OBJLIB)
	echo Success!
	
xbeamtest: $(OBJLIB) $(OBJXBEAM)
	echo Success!

dirprint:
	echo Test Directories Path and Content
	echo ROOTDIR: $(ROOTDIR)
	echo Forward Solver
	echo MAINDIR: $(MAINDIR)
	echo XBEAMDIR: $(XBEAMDIR)
	echo LIBDIR: $(LIBDIR)
	echo Object Folder and Wrapper
	echo OBJDIR: $(OBJDIR)
	echo WRAPDIR: $(WRAPDIR)

dirlist:
	ls $(ROOTDIR)
	ls $(MAINDIR)
	ls $(XBEAMDIR)
	ls $(LIBDIR)
	ls $(OBJDIR)
	ls $(OPTDIR)
	ls $(OPTMAINDIR)
	ls $(OPTSRCDIR)
	ls $(OPTLIBSDIR)
	ls $(OPTOBJDIR)
	ls $(OPTBINDIR)	
	ls $(WRAPDIR)

objlist:
	echo $(OBJLIB)
	echo $(OBJXBEAM)
	echo $(OBJMAIN)
	echo $(OBJWRAP)

################################################################################
## COMPILE DEPENDENCIES

# -------------------------------------------------------------------- Optimiser
# Opt Main
$(OPTOBJDIR)/%.o: $(OPTMAINDIR)/%.f90
	$(F90COMP) $(F90OPTS) -o $@ -c $< 	
# Opt SRC
$(OPTOBJDIR)/%.o: $(OPTSRCDIR)/%.f90
	$(F90COMP) $(F90OPTS) -o $@ -c $< 	
# Opt Libraries
$(OPTOBJDIR)/%.o: $(OPTLIBDIR)/%.f90
	$(F90COMP) $(F90OPTS) -o $@ -c $< 
	
# -------------------------------------------------------------- Forward Problem
# Main 
$(OBJDIR)/%_$(TESTCASE).o: $(MAINDIR)/%_$(TESTCASE).f90;
	$(F90COMP) $(F90OPTS) -o $@ -c $< 
# xbeam
$(OBJDIR)/%.o: $(XBEAMDIR)/%.f90;
	$(F90COMP) $(F90OPTS) -o $@ -c $< 
# libraries
$(OBJDIR)/%.o: $(LIBDIR)/%.f90
	$(F90COMP) $(F90OPTS) -o $@ -c $< 
	
# ---------------------------------------------------------------------- Wrapper
$(OBJDIR)/%.o: $(WRAPDIR)/%.c;
	$(F90COMP) $(F90OPTS) -o $@ -c $< 


################################################################################
# Unused Code
#
# F90MAIN=$(MAINDIR)/*.f90       <-- gives the list BUT...
# OBJMAIN=$(F90MAIN:.f90)=.o <--- ...this won't work
# 
# Notes for Developer: failed attempts to list automatically the f90 files
# ls $(MAINDIR)/*.f90                  <-- lists all files in a folder
# ls $(MAINDIR)/*.f90 > test           <-- put them in a write
# grep -i 'f90' test                   <-- take them from a file
# F90MAIN=$(ls $(MAINDIR)/*.f90)       <-- gives nothing...
# comnd="ls $(MAINDIR)/*.f90"
# F90MAIN={eval $(comnd)}
#
#
# Not found out how to add the path to the file name. Object files will be built 
# in the f90 folders and then moved
#F90MAIN=\
input_$(TESTCASE).f90 \
main_$(TESTCASE).f90
#
#OBJFILE=$(F90MAIN:.f90=.o)
#F90LIB=$(LIBSSRCDIR)/lib_rot.f90 \
     $(LIBSSRCDIR)/lib_rotvect.f90
#OBJLIB=$(F90LIB:.f90=.o)
# OBJMAIN=$(OBJDIR)/$(OBJFILE) doesn't work
#
# --------------------------------------------------------------- Create Objects 
# this places the objects in the same folder as f90
#F90_TEST=$(SRC_TEST:.o=.f90)
#SRC_TOOLS=$(F90_TOOLS:.f90=.o)
#F90_XBEAM=$(SRC_XBEAM:.o=.f90)
#F90_TEST=$(SRC_TEST:.o=.f90)
#F90_WRAP=$(SRC_WRAP:.o=.f90)
#F90_MAINLIB=$(SRC_MAINLIB:.o=.f90)
#
#lib: $(OBJLIB)
#	gfortran -o DynLib.so $^ $(LIBFLAG)
#	echo done!!!
#libclean: 
#	rm $(LIBSSRCDIR)/lib*.o
#
# works!	
#%.o: %.f90
#	gfortran -o $@ -c $< $(F90FLAG)
#
# works!
# note that % matches the filename, not the full path, so also this works! 
#lib%.o: lib%.f90
#	gfortran -o $@ -c $<
#    mv $@ $(OBJDIR) ### works, but than dependencies are not found