#const static char cvsid[] = "$Id: def.makefile,v 1.20 2003/07/22 21:52:55 zhzhu Exp $";

#! /bin/sh

# For debugging:
#GDB_FLAG    = -g
#MORE_GDB_FLAGS    = -fno-for-scope

# For optimizing
OPTIM_FLAGS_FOR_C++ = -O3
OPTIM_FLAGS_FOR_C = -O2

#For performance meter
#PERFORMANCE_METER = -pg

# Skip generation of a temp variable in constructing an object
CONSTRUCTOR = -felide-constructors

# enable PIC (position-independent code) to avoid the
# "relocation truncated to fit: R_X86_64_PC32" error.
# The error is caused by linking a shared library (which requires position-independent code,
# PIC) to a static library (which has not been compiled with PIC). You need to either link
# the shared library against a shared version of the static code (such as is produced
# automatically by libtool), or re-compile the static library with PIC enabled.
#PIC=-fPIC
# enable the right memory model to avoid the "relocation truncated to fit: R_X86_64_PC32"
# issue, due to assumptions about the program and symbols being in the lower 2G
# of address space. With 'medium', symbols with sizes larger than `-mlarge-data-threshold'
# are put into large data or BSS sections and can be located above 2GB.
# With 'large', the compiler makes no assumptions about addresses and sizes of sections,
# at the cost of a slightly larger overhead when calling functions.
MCMODEL=-mcmodel=large

# For compilation warning message:
WARNING_FLAGS  = -ansi
#MORE_WARNING_FLAGS  = -Wall -pedantic 

# For debug 
#PRINT_PROGRESS = -DDEBUG_PROGRESS
#DEBUG_ELEMENT = -DDEBUG_ELEMENT
#DEBUG_GRID = -DDEBUG_GRID
#MEMORY_REPORT = -DMEMORY_REPORT
#PRINT_RUN_TIME = -DTIME
#DEBUG_FLAG= -DDEBUG
#DISABLE_CONSISTENT_STENCIL= -DDISABLE_CONSISTENT_STENCIL
PRINT_ELEMENT = -DPRINT_ELEMENT

# if you installed the newer version gcc and g++, it will be in here
#GCC_ROOT         = /usr/local/bin
# the default place for gcc and g++
GCC_ROOT         = /usr/bin
LD_ROOT         = /usr/bin
BIN_ROOT        = /bin
SHELL		= $(BIN_ROOT)/sh
RM		= $(BIN_ROOT)/rm -f
RMR             = $(BIN_ROOT)/rm -rf
SED        	= $(BIN_ROOT)/sed
MV	        = $(BIN_ROOT)/mv
RANLIB          = /usr/bin/ranlib
AR              = ar
LD              = $(LD_ROOT)/ld
CC              = $(GCC_ROOT)/gcc
C++		 = $(GCC_ROOT)/g++

ROOT = ../..
INC_DIR = $(ROOT)/inc
BIN_DIR = $(ROOT)/bin
UTIL_DIR = $(ROOT)/util
FFTW_INC_DIR = $(ROOT)/src/fftw-2.1.3/include
DEF_MAKEFILE = $(ROOT)/def.makefile

