#!/usr/bin/env python

import os
import sys

'''
This script builds a Makefile for the Tinker unit test suite.  It is the
test-directory counterpart of make/buildmake.py: it scans the fixed-form
Fortran test files, categorizes them as modules, programs or subroutines,
and emits a Makefile that compiles them and links each test program
against the prebuilt Tinker library (../source/libtinker.a) and FFTW.

Unlike the main buildmake.py it does NOT rebuild libtinker.a; the library
and its *.mod files are assumed to already exist in the Tinker source
directory.  Module dependencies are emitted only for modules defined among
the test files themselves (e.g. assert); modules that live in the
Tinker library (atoms, energi, ...) are resolved via -I$(TINKERSRC) and at
link time, so they are not listed as make prerequisites.

Usage (run from the test directory):   ./buildmake.py *.f > Makefile
'''

MAKEFILE_CONFIG = '''##
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
'''

UNKNOWN_TYPE = 'UNKNOWN_TYPE'
MODULE_TYPE = 'MODULE_TYPE'
MODULE_FILES = []
PROGRAM_TYPE = 'PROGRAM_TYPE'
PROGRAM_FILES = []
SUBROUTINE_TYPE = 'SUBROUTINE_TYPE'
SUBROUTINE_FILES = []
DEPENDENCY = []


class dependency_record:
    def __init__(self):
        self.target = ''
        self.depend = []

    def __str__(self):
        line = self.target + ':'
        for d in self.depend:
            line = line + ' ' + d + '.o'
        return line


def determine_module_subroutine_program(fortran_filename):
    global MODULE_FILES
    global PROGRAM_FILES
    global SUBROUTINE_FILES

    base = os.path.basename(fortran_filename)  # 'amoeba.f'
    stem, ext = os.path.splitext(base)         # 'amoeba', '.f'
    content = [line.rstrip().lower() for line in open(fortran_filename)]

    use_list = []
    file_type = UNKNOWN_TYPE

    def type_must_be_unknown(t):
        if t != UNKNOWN_TYPE:
            raise BaseException('Cannot parse file: %s' % fortran_filename)

    for line in content:
        if len(line) > 13 and line[0:13] == '      module ':
            type_must_be_unknown(file_type)
            file_type = MODULE_TYPE
            MODULE_FILES.append(stem)
        if len(line) > 14 and line[0:14] == '      program ':
            type_must_be_unknown(file_type)
            file_type = PROGRAM_TYPE
            PROGRAM_FILES.append(stem)
        if len(line) > 10 and line[0:10] == '      use ':
            temp_list = line.split()
            if len(temp_list) > 1:
                word = temp_list[1].rstrip(',')
                if word not in use_list:
                    use_list.append(word)
    if file_type == UNKNOWN_TYPE:
        file_type = SUBROUTINE_TYPE
        SUBROUTINE_FILES.append(stem)

    use_list.sort()
    return stem, file_type, use_list


def categorize_fortran_files(fortran_files):
    global MODULE_FILES
    global PROGRAM_FILES
    global SUBROUTINE_FILES
    global DEPENDENCY

    raw = {}
    for filename in fortran_files:
        stem, filetype, use_list = determine_module_subroutine_program(filename)
        raw[stem] = use_list

    MODULE_FILES.sort()
    PROGRAM_FILES.sort()
    SUBROUTINE_FILES.sort()

    # keep only dependencies on modules defined by the test files; modules
    # from the Tinker library are found via -I$(TINKERSRC) and at link time
    local_modules = set(MODULE_FILES)
    DEPENDENCY = []
    for stem in sorted(raw.keys()):
        depend = [m for m in raw[stem] if m in local_modules and m != stem]
        if depend:
            dr = dependency_record()
            dr.target = stem + '.o'
            dr.depend = depend
            DEPENDENCY.append(dr)


def common_list():
    # shared objects linked into every test program (modules + subroutines)
    temp_list = MODULE_FILES + SUBROUTINE_FILES
    temp_list.sort()
    return [word + '.o' for word in temp_list]


def x_list():
    return [word + '.x' for word in PROGRAM_FILES]


def print_list(init_string, items, end_string=''):
    buff = init_string
    length = len(items)
    if length > 0:
        for i in range(length - 1):
            buff = buff + ' ' + items[i]
            if len(buff) > 60:
                buff = buff + ' \\'
                print(buff)
                buff = ''
        buff = buff + ' ' + items[length - 1]
    print(buff)
    print(end_string)


def target_o():
    print('%.o: $(src)/%.f')
    print('\t$(F77) $(F77FLAGS) $(OPTFLAGS) $< -o $@')
    print('')


def target_x():
    # each test program links the shared objects plus the Tinker library
    for prog in PROGRAM_FILES:
        print('%s.x: $(COMMON) %s.o $(TINKERLIB)' % (prog, prog))
        print('\t$(F77) $(LINKFLAGS) -o $@ $(COMMON) %s.o $(TINKERLIB) \\'
              % prog)
        print('\t    $(FFTW_LIBDIR) $(FFTW_LIBS) $(SYSLIBS)')
        print('')


def all_test_clean():
    print('all: $(EXEFILES)')
    print('')

    print('test: $(EXEFILES)')
    for prog in PROGRAM_FILES:
        print('\t./%s.x $(FILTER) $(DETAIL)' % prog)
    print('')

    print('clean:')
    print('\trm -f *.o *.mod *.x')
    print('')


def print_dependency():
    if not DEPENDENCY:
        return
    print('###############################################################')
    print('##  Explicit Dependencies on Test Modules                    ##')
    print('###############################################################')
    print('')
    for d in DEPENDENCY:
        print(d)


if __name__ == '__main__':
    categorize_fortran_files(sys.argv[1:])

    print(MAKEFILE_CONFIG)

    print_list('COMMON :=', common_list())
    print_list('EXEFILES :=', x_list())

    target_o()
    target_x()
    all_test_clean()
    print_dependency()
