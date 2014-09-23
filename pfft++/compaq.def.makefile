#const static char cvsid[] = "$Id: compaq.def.makefile,v 1.2 2002/10/01 16:37:56 zhzhu Exp $";

# For debugging:
#GDB_FLAG    = -g
#MORE_GDB_FLAGS    = -fno-for-scope

# For optimizing
OPTIM_FLAGS_FOR_C++ = -O3
OPTIM_FLAGS_FOR_C = -O6

# to enforce usage of std iostream in Compaq compilor
USE_STD = -D__USE_STD_IOSTREAM

# This is a work around for bugs in Compaq C++ linker
NESTED_ENUMS = -distinguish_nested_enums
USE_CNAME = -pure_cname
ASSUME_GFULLPATH = -assume gfullpath

# Compaq compilor creates a subdirectory under current sub-directory
# and put the instatiation of templates there. This is not the case 
# for g++. To make clean, I have to delete this subdirectory.
COMPAQ_REPOSITORY = cxx_repository

#For performance meter
#PERFORMANCE_METER = -pg

# Skip generation of a temp variable in constructing an object
#CONSTRUCTOR = -felide-constructors

#some compilation options necessary for compaq OS
FOR_COMPAQ_1 = -D_XOPEN_SOURCE

# For compilation warning message:
WARNING_FLAGS  = -ansi
#MORE_WARNING_FLAGS  = -Wall -pedantic 

# For debug info print out
#PRINT_PROGRESS = -DDEBUG_PROGRESS
#DEBUG_ELEMENT = -DDEBUG_ELEMENT
#DEBUG_GRID = -DDEBUG_GRID
#MEMORY_REPORT = -DMEMORY_REPORT
#DEBUG_FLAG= -DDEBUG
#DEBUG_1 = -DDEBUG_LEVEL1
#PRINT_RUN_TIME = -DTIME

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
GCC              = gcc
G++		 = g++
C++		= cxx
CC              = cc

ROOT = ../..
INC_DIR = $(ROOT)/inc
BIN_DIR = $(ROOT)/bin
UTIL_DIR = $(ROOT)/util
FFTW_INC_DIR = $(ROOT)/src/fftw-2.1.3/include
DEF_MAKEFILE = $(ROOT)/def.makefile

