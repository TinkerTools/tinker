##
###################################################################
##                                                               ##
##  Makefile for Building the Tinker Unit Test Suite             ##
##                                                               ##
###################################################################
##
##  Prerequisite:  build the Tinker library first, so that
##  $(TINKERSRC)/libtinker.a and $(TINKERSRC)/*.mod exist.
##
##  make            build the test executables
##  make test       build and run the tests (pass FILTER=<tag> to select)
##  make clean      remove objects, modules and executables

# test source directory
src := $(PWD)
# Tinker source directory holding libtinker.a and the *.mod files
TINKERSRC := $(src)/../source
# release / debug / profile
opt := release
# gfortran / gfortran8 / ifort / etc.
F77 := gfortran
# FFTW directory
fftw := default__

ifeq ($(fftw), default__)
  FFTWDIR := $(src)/../fftw
else
  FFTWDIR := $(fftw)
endif
FFTW_LIBDIR := -L$(FFTWDIR)/lib
FFTW_LIBS := -lfftw3_threads -lfftw3

# optional test selection filter passed through to the binaries
FILTER :=
# set DETAIL=-v to also print per-check detail lines
DETAIL :=

######################
## Operating System ##
######################
# Darwin (macOS) / Linux / CYGWIN_NT
os__ := $(shell uname -s)

ifeq ($(os__), Darwin)
  SYSLIBS :=
else ifeq ($(os__), Linux)
  SYSLIBS := -Wl,--no-as-needed -ldl
else ifeq ($(shell echo $(os__) | cut -c 1-9), CYGWIN_NT)
  SYSLIBS := -ldl
else
$(error Unknown OS -- $(os__))
endif

##################
## F77 Compiler ##
##################
use_gfortran__ := false
use_ifort__ := false
found__ := false

f77__ := $(shell echo $(F77) | cut -c 1-8)
ifeq ($(f77__), gfortran)
  use_gfortran__ := true
  found__ := true
endif
f77__ := $(shell echo $(F77) | cut -c 1-5)
ifeq ($(f77__), ifort)
  use_ifort__ := true
  found__ := true
endif
ifneq ($(found__), true)
$(error Unknown fortran compiler -- $(F77))
endif

###################
## Compile Flags ##
###################

ifeq ($(use_gfortran__), true)
  ifeq ($(opt), release)
    OPTFLAGS := -O3 -fopenmp
  else
    OPTFLAGS := -g -fopenmp
  endif
  F77FLAGS := -c -I$(TINKERSRC)
  LINKFLAGS := $(OPTFLAGS) -static-libgcc
endif

ifeq ($(use_ifort__), true)
  ifeq ($(opt), release)
    OPTFLAGS := -O3 -no-ipo -no-prec-div -recursive -qopenmp
  else
    OPTFLAGS := -g -qopenmp
  endif
  F77FLAGS := -c -I$(TINKERSRC)
  LINKFLAGS := $(OPTFLAGS) -static-libgcc -static-intel
endif

TINKERLIB := $(TINKERSRC)/libtinker.a

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

COMMON := assert.o test_amoeba.o test_angle.o test_pdbxyz.o test_potential.o \
 test_xyzpdb.o

EXEFILES := tinkertest.x

%.o: $(src)/%.f
	$(F77) $(F77FLAGS) $(OPTFLAGS) $< -o $@

tinkertest.x: $(COMMON) tinkertest.o $(TINKERLIB)
	$(F77) $(LINKFLAGS) -o $@ $(COMMON) tinkertest.o $(TINKERLIB) \
	    $(FFTW_LIBDIR) $(FFTW_LIBS) $(SYSLIBS)

all: $(EXEFILES)

test: $(EXEFILES)
	./tinkertest.x $(FILTER) $(DETAIL)

clean:
	rm -f *.o *.mod *.x

###############################################################
##  Explicit Dependencies on Test Modules                    ##
###############################################################

tinkertest.o: assert.o
