#!/usr/bin/env python

import os
import sys

'''
This script only works with the fixed-from fortran source code in Tinker.
If the fortran file contains either a module or a program, their names are
expected to be the same.
'''

MAKEFILE_CONFIG = '''##
###################################################################
##                                                               ##
##  Makefile for Building the Tinker Molecular Modeling Package  ##
##                                                               ##
###################################################################
##

# source directory
src := $(PWD)
# release / debug / profile / etc.
opt := release
# gfortran / gfortran8 / gfortran-8 / ifort / ifort14 / ifort-14 / etc.
F77 := gfortran
# FFTW directory
fftw := default__

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##  TINKERDIR       Tinker Distribution Directory
##  TINKER_LIBDIR   Libraries needed to build Tinker
##  BINDIR          Directory with Tinker Executables
##  LINKDIR         Linked Copies of Tinker Executables
##  FFTWDIR         FFTW Top-Level Directory
##  FFTW_LIB_DIR    Directory with FFTW Libraries
##  FFTW_LIBS       FFTW Libraries needed to build Tinker
##  APBSDIR         APBS Top-Level Directory
##  APBS_INC_DIR    Directory with APBS Include Files
##  APBS_LIB_DIR    Directory with APBS Libraries
##  APBS_LIBS       APBS Libraries needed to build Tinker

TINKERDIR := $(HOME)/tinker
TINKER_LIBDIR := $(TINKERDIR)/lib
BINDIR := $(TINKERDIR)/bin
LINKDIR := /usr/local/bin

ifeq ($(fftw), default__)
  FFTWDIR := $(TINKERDIR)/fftw
else
  FFTWDIR := $(fftw)
endif
FFTW_LIBDIR := -L$(FFTWDIR)/lib
FFTW_LIBS := -lfftw3_threads -lfftw3

APBSDIR := $(TINKERDIR)/apbs
APBS_INCDIR := -I$(APBSDIR)/include
APBS_LIBDIR := -L$(APBSDIR)/lib
APBS_LIBS := -lapbsmainroutines -lapbs -lmaloc -lapbsblas

######################
## Operating System ##
######################
# Darwin (macOS) / Linux / CYGWIN_NT
os__ := $(shell uname -s)

ifeq ($(os__), Darwin)
else ifeq ($(os__), Linux)
else ifeq ($(shell echo $(os__) | cut -c 1-9), CYGWIN_NT)
else
$(error Unknown OS -- $(os__); Please help with us)
endif

##################
## F77 Compiler ##
##################
use_gfortran__ := false
use_ifort := false
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
$(error Unknown fortran compiler -- $(F77); Please help with us)
endif

###################
## Compile Flags ##
###################

ifeq ($(use_gfortran__), true)
  ifeq ($(os__), Linux)
    ifeq ($(opt), release)
      OPTFLAGS := -Ofast -msse3 -fopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c
    LIBDIR := -L. -L$(TINKER_LIBDIR)/linux
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib
    LINKFLAGS := $(OPTFLAGS) -static-libgcc
    RENAME := rename_bin
  endif

  ifeq ($(os__), Darwin)
    ifeq ($(opt), release)
      OPTFLAGS := -Ofast -mssse3 -fopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c
    LIBDIR := -L. -L$(TINKER_LIBDIR)/macos
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib -c
    LINKFLAGS := $(OPTFLAGS) -static-libgcc
    RENAME := rename_bin
  endif

  ifeq ($(shell echo $(os__) | cut -c 1-9), CYGWIN_NT)
    ifeq ($(opt), release)
      OPTFLAGS := -Ofast -msse3 -fopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c
    LIBDIR := -L. -L$(TINKER_LIBDIR)/windows
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib
    LINKFLAGS := $(OPTFLAGS) -static
    RENAME := rename_exe
  endif
endif

ifeq ($(use_ifort__), true)
  ifeq ($(os__), Linux)
    ifeq ($(opt), release)
      OPTFLAGS := -O3 -no-ipo -no-prec-div -recursive -qopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c -xHost
    LIBDIR := -L. -L$(TINKER_LIBDIR)/linux
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := echo
    LINKFLAGS := $(OPTFLAGS) -static-libgcc -static-intel
    RENAME := rename_bin
  endif

  ifeq ($(os__), Darwin)
    ifeq ($(opt), release)
      OPTFLAGS := -O3 -no-ipo -no-prec-div -mdynamic-no-pic -qopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c -axSSSE3
    LIBDIR := -L. -L$(TINKER_LIBDIR)/macos
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib -c
    LINKFLAGS := $(OPTFLAGS) -static-intel -Wl,-stack_size,0x10000000
    RENAME := rename_bin
  endif
endif

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
        line = self.target+':'
        for d in self.depend:
            line = line+' '+d+'.o'
        return line

def determine_module_subroutine_program(fortran_filename):
    global MODULE_FILES
    global PROGRAM_FILES
    global SUBROUTINE_FILES

    base = os.path.basename(fortran_filename) # 'source/dynamic.f' -> 'dynamic.f'
    stem, ext = os.path.splitext(base) # 'dynamic.f' -> 'dynamic', '.f'
    content = [line.rstrip().lower() for line in open(fortran_filename)]

    use_list = []
    file_type = UNKNOWN_TYPE
    def type_must_be_unknown(t):
        if t != UNKNOWN_TYPE:
            raise BaseException('Cannot parse file: %s' % fortran_filename)

    for line in content:
        # '      module foo'
        if len(line) > 13:
            if line[0:13] == '      module ':
                type_must_be_unknown(file_type)
                file_type = MODULE_TYPE
                MODULE_FILES.append(stem)
        # '      program foo'
        if len(line) > 14:
            if line[0:14] == '      program ':
                type_must_be_unknown(file_type)
                file_type = PROGRAM_TYPE
                PROGRAM_FILES.append(stem)
        # '      use foo'
        if len(line) > 10:
            if line[0:10] == '      use ':
                temp_list = line.split()
                if len(temp_list) > 1:
                    word = temp_list[1]
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

    temp_dict = {}
    for filename in fortran_files:
        stem, filetype, use_list = determine_module_subroutine_program(filename)
        dr = dependency_record()
        dr.target = stem+'.o'
        dr.depend = use_list
        temp_dict[dr.target] = dr
    sorted_keys = [k for k in temp_dict.keys()]
    sorted_keys.sort()

    MODULE_FILES.sort()
    PROGRAM_FILES.sort()
    SUBROUTINE_FILES.sort()
    DEPENDENCY = [temp_dict[k] for k in sorted_keys]

def x_list():
    return [word+'.x' for word in PROGRAM_FILES]

def xobj_list():
    return [word+'.o' for word in PROGRAM_FILES]

def lib_list():
    temp_list = MODULE_FILES+SUBROUTINE_FILES
    temp_list.sort()
    return [word+'.o' for word in temp_list]

def print_list(init_string, list, end_string=''):
    buff = init_string
    length = len(list)
    if length > 0:
        for i in range(length-1):
            buff = buff+' '+list[i]
            if len(buff) > 60:
                buff = buff+' \\'
                print(buff)
                buff = ''
        buff = buff+' '+list[length-1]
    print(buff)
    print(end_string)

def print_dependency():
    print('''###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################
''')
    for d in DEPENDENCY:
        print(d)

def target_o():
    print('%.o: $(src)/%.f')
    print('\t$(F77) $(F77FLAGS) $(OPTFLAGS) $< -o $@')
    print('')

def target_x():
    print('%.x: %.o libtinker.a')
    print('\t$(F77) $(LINKFLAGS) -o $@ $(LIBDIR) $(FFTW_LIBDIR) $^ $(LIBS) $(FFTW_LIBS); strip $@')
    print('')

def all_install_clean_listing():
    print('all: $(EXEFILES)')
    print('')

    print('install: $(RENAME)')
    print('')

    print('clean:')
    print('\trm -f *.o *.mod *.a *.x')
    print('')

    print('listing:')
    print('\tcat $(src)/*.f $(src)/*.c > tinker.txt')
    print('')

def rename_bin_rename_exe_remove_links_create_links():
    print('rename_bin:')
    for x in PROGRAM_FILES:
        print('\tmv %-20s $(BINDIR)/%s' % (x+'.x', x))
    print('')

    print('rename_exe:')
    for x in PROGRAM_FILES:
        print('\tmv %-20s $(BINDIR)/%s.exe' % (x+'.x', x))
    print('')

    print('remove_links:')
    for x in PROGRAM_FILES:
        print('\trm -f $(LINKDIR)/%s' % x)
    print('')

    print('create_links:')
    for x in PROGRAM_FILES:
        print('\tln -s $(BINDIR)/%-15s $(LINKDIR)/%s' % (x, x))
    print('')

if __name__ == '__main__':
    categorize_fortran_files(sys.argv[1:])

    print(MAKEFILE_CONFIG)

    print_list('LIBOBJS :=', lib_list())
    print_list('EXEOBJS :=', xobj_list())
    print('OBJS := $(LIBOBJS) $(EXEOBJS)')
    print('')
    print_list('EXEFILES :=', x_list())

    target_o()
    target_x()
    all_install_clean_listing()
    rename_bin_rename_exe_remove_links_create_links()

    print_list('libtinker.a: $(LIBOBJS)\n\tar $(LIBFLAGS) libtinker.a $(LIBOBJS)', [], '\t$(RANLIB) libtinker.a\n')
    print_dependency()
