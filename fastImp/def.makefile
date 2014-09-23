#const static char cvsid[] = "$Id: def.makefile,v 1.35 2003/07/22 21:31:00 zhzhu Exp $";

#! /bin/sh

# For debugging:
#GDB_FLAG    = -g
#MORE_GDB_FLAGS    = -fno-for-scope

# For optimizing
OPTIM_FLAGS_FOR_C++ = -O3
OPTIM_FLAGS_FOR_C = -O2

# For fastimp debug info 
SURF_RUN_TIME = -DPRINT_TIME
SURF_PROGRESS = -DPRINT_PROGRESS
SURF_MEMORY = -DPRINT_MEMORY

# For fastimp system matrix debug
#DEBUG_VERSION = -DDEBUG_VERSION
#DEBUG_PRE_COND = -DDEBUG_PRE_COND
#DEBUG_DIRECT = -DDEBUG_DIRECT
#DEBUG_MEMORY = -DDEBUG_MEMORY
#DIAG_LOCAL_PRECOND = -DDIAG_LOCAL_PRECOND
#USE_FIXED_BUFFER_SIZE = -DUSE_FIXED_BUFFER_SIZE
#DEBUG_MESH_RECT_SPIRAL -DDEBUG_MESH_RECT_SPIRAL
#DEBUG_EMQS = -DDEBUG_EMQS
#DEBUG_RHO = -DDEBUG_RHO
#DEBUG_GMRES = -DDEBUG_GMRES
#NO_GRAY_WINDOW = -DNO_GRAY_WINDOW
#LARGE_STEP_RATIO = -DLARGE_STEP_RATIO
#HIGH_ACCURACY = -DHIGH_ACCURACY
#DEBUG_TRANSLINE = -DDEBUG_TRANSLINE
#DISABLE_SCALING = -DDISABLE_SCALING
#PRINT_GMRES_RESID = -DPRINT_GMRES_RESID
#DEBUG_CURRENT = -DDEBUG_CURRENT
#DEBUG_ITERATIVE = -DDEBUG_ITERATIVE
#DEBUG_BUFFER = -DDEBUG_BUFFER
#DEBUG_CONTACT = -DDEBUG_CONTACT

# Skip generation of a temp variable in constructing an object
CONSTRUCTOR = -felide-constructors

# For profiling
#GPROF_FLAG = -pg

# For compilation warning message:
WARNING_FLAGS  = -ansi
#MORE_WARNING_FLAGS  = -Wall -pedantic 

# For pfft++ debug info 
#PFFT_RUN_TIME = -DTIME
#PFFT_PROGRESS = -DDEBUG_PROGRESS
#DEBUG_ELEMENT = -DDEBUG_ELEMENT
#DEBUG_GRID = -DDEBUG_GRID
#MEMORY_REPORT = -DMEMORY_REPORT
#DEBUG_FLAG= -DDEBUG
#DISABLE_CONSISTENT_STENCIL= -DDISABLE_CONSISTENT_STENCIL
#PRINT_ELEMENT = -DPRINT_ELEMENT
#DISABLE_PIECEWISE_SCHEME = -DDISABLE_PIECEWISE_SCHEME

LD_ROOT         = /usr/bin
# if you installed the newer version gcc and g++, it will be in here
#GCC_ROOT         = /usr/local/bin
# the default place for gcc and g++
GCC_ROOT         = /usr/bin
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
LIB_DIR = $(ROOT)/lib
BIN_DIR = $(ROOT)/bin
UTIL_DIR = $(ROOT)/util
DEF_MAKEFILE = $(ROOT)/def.makefile

PFFT_ROOT =  ../../../pfft++
PFFT_INC_DIR    = $(PFFT_ROOT)/inc
PFFT_SRC_DIR    = $(PFFT_ROOT)/src/pfft
PFFT_LIB_DIR    = $(PFFT_ROOT)/lib

SUPER_LU_INC_DIR    = ../SuperLU/SRC
