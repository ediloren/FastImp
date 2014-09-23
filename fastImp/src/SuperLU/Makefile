############################################################################
#
#  Program:         SuperLU
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   October 2, 1995
#
#  Modified:        February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################

include make.inc

all: install lib testing

lib: superlulib_complex16 blaslib_complex16
#lib: superlulib tmglib blaslib

clean: cleanlib cleantesting

install:
	( cd INSTALL; $(MAKE) )
#	( cd INSTALL; cp lsame.c ../SRC/; \
#	  cp dlamch.c ../SRC/; cp slamch.c ../SRC/ )

blaslib:
	( cd CBLAS; $(MAKE) )

blaslib_complex16:
	( cd CBLAS; $(MAKE) complex16)

superlulib:
	( cd SRC; $(MAKE) )

superlulib_complex16:
	( cd SRC; $(MAKE) complex16)

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

matlabmex:
	( cd MATLAB; $(MAKE) )

testing:
	( cd TESTING ; $(MAKE) )

cleanlib:
	( cd SRC; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd CBLAS; $(MAKE) clean )

cleantesting:
	( cd INSTALL; $(MAKE) clean )
	( cd TESTING; $(MAKE) clean )
	( cd MATLAB;  $(MAKE) clean )
	( cd EXAMPLE; $(MAKE) clean )

etags:
	( cd SRC; $(MAKE) etags )

