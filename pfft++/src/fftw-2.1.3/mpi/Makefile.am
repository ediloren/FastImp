# The mpi stuff is only compiled if FFTW was configured with
# the --enable-mpi option.  This requires us to use the combination
# of foo_* and EXTRA_* variables below.

lib_LTLIBRARIES = @FFTW_MPI_LIBLIST@
include_HEADERS = @FFTW_MPI_INCLUDELIST@
noinst_PROGRAMS = @FFTW_MPI_PROGLIST@

EXTRA_LTLIBRARIES = libXXX_FFTW_PREFIX_XXXfftw_mpi.la \
                    libXXX_FFTW_PREFIX_XXXrfftw_mpi.la
EXTRA_HEADERS = XXX_FFTW_PREFIX_XXXfftw_mpi.h XXX_FFTW_PREFIX_XXXrfftw_mpi.h 
EXTRA_PROGRAMS = test_sched test_transpose_mpi fftw_mpi_test rfftw_mpi_test
EXTRA_DIST = README.f77

CC=@MPICC@

FFTWDIR=../fftw
RFFTWDIR=../rfftw
INCLUDES = -I$(srcdir)/../tests \
           -I$(srcdir)/$(FFTWDIR) -I$(srcdir)/$(RFFTWDIR) -I$(srcdir)

libXXX_FFTW_PREFIX_XXXfftw_mpi_la_SOURCES = \
                         sched.c sched.h TOMS_transpose.c TOMS_transpose.h \
                         transpose_mpi.c fftwnd_mpi.c fftw_mpi.c \
                         XXX_FFTW_PREFIX_XXXfftw_mpi.h \
			 fftw_f77_mpi.c fftw_f77_mpi.h
libXXX_FFTW_PREFIX_XXXfftw_mpi_la_LDFLAGS = \
                         -version-info @SHARED_VERSION_INFO@ \
                         -rpath $(libdir)

XXX_FFTW_PREFIX1_XXXfftw_mpi.h: fftw_mpi.h
	rm -f XXX_FFTW_PREFIX_XXXfftw_mpi.h
	sed 's/<fftw/<XXX_FFTW_PREFIX_XXXfftw/g;s/<rfftw/<XXX_FFTW_PREFIX_XXXrfftw/g' fftw_mpi.h > XXX_FFTW_PREFIX_XXXfftw_mpi.h

libXXX_FFTW_PREFIX_XXXrfftw_mpi_la_SOURCES = \
                          rfftwnd_mpi.c XXX_FFTW_PREFIX_XXXrfftw_mpi.h \
			  rfftw_f77_mpi.c
libXXX_FFTW_PREFIX_XXXrfftw_mpi_la_LDFLAGS = \
                          -version-info @SHARED_VERSION_INFO@ \
                          -rpath $(libdir)

XXX_FFTW_PREFIX1_XXXrfftw_mpi.h: rfftw_mpi.h
	rm -f XXX_FFTW_PREFIX_XXXrfftw_mpi.h
	sed 's/<fftw/<XXX_FFTW_PREFIX_XXXfftw/g;s/<rfftw/<XXX_FFTW_PREFIX_XXXrfftw/g' rfftw_mpi.h > XXX_FFTW_PREFIX_XXXrfftw_mpi.h

CLEANFILES = XXX_FFTW_PREFIX1_XXXfftw_mpi.h XXX_FFTW_PREFIX1_XXXrfftw_mpi.h

test_sched_SOURCES = test_sched.c
test_sched_LDADD = libXXX_FFTW_PREFIX_XXXfftw_mpi.la \
                   $(FFTWDIR)/libXXX_FFTW_PREFIX_XXXfftw.la @MPILIBS@

test_transpose_mpi_SOURCES = test_transpose_mpi.c
test_transpose_mpi_LDADD = libXXX_FFTW_PREFIX_XXXfftw_mpi.la \
                           $(FFTWDIR)/libXXX_FFTW_PREFIX_XXXfftw.la @MPILIBS@

fftw_mpi_test_SOURCES = fftw_mpi_test.c
fftw_mpi_test_LDADD = ../tests/test_main.o \
                      libXXX_FFTW_PREFIX_XXXfftw_mpi.la \
                      $(FFTWDIR)/libXXX_FFTW_PREFIX_XXXfftw.la @MPILIBS@

rfftw_mpi_test_SOURCES = rfftw_mpi_test.c
rfftw_mpi_test_LDADD = ../tests/test_main.o \
                       libXXX_FFTW_PREFIX_XXXrfftw_mpi.la \
                       libXXX_FFTW_PREFIX_XXXfftw_mpi.la \
		       $(RFFTWDIR)/libXXX_FFTW_PREFIX_XXXrfftw.la \
                       $(FFTWDIR)/libXXX_FFTW_PREFIX_XXXfftw.la \
                       @MPILIBS@

# for some reason, automake tries to use autoheader in order to
# generate config.h.in, and fails because config.h.in is GNU-lly
# incorrect.  Just disable autoheader
AUTOHEADER=echo
