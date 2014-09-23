FFTWDIR = ../fftw
FFTW_INCLUDE = -I$(FFTWDIR)
OTHER_FFTW_OBJECTS = malloc.o twiddle.o planner.o config.o putils.o rader.o
FFTW_LIB_NON_CILK = $(FFTWDIR)/f*.o $(FFTWDIR)/generic.o \
                    $(FFTWDIR)/timer.o $(FFTWDIR)/executor.o \
                    $(FFTWDIR)/wisdom.o $(FFTWDIR)/wisdomio.o 

FFTW_LIB = $(FFTW_LIB_NON_CILK) $(OTHER_FFTW_OBJECTS)

ALL_FFTW_CILK_OBJ = $(OTHER_FFTW_OBJECTS) executor_cilk.o fftwnd_cilk.o  \
                    $(FFTW_LIB_NON_CILK)

LIBFFTW_CILK = libfftw_cilk.a

# Set where you want to install the library for "make install"
prefix = /usr/local
LIBDIR = $(prefix)/lib
INCLUDEDIR = $(prefix)/include

# On systems that have ranlib:
RANLIB = ranlib

# On systems that don't have ranlib
# RANLIB = echo ranlib

all:	$(LIBFFTW_CILK)

tests: test_cilk time_cilk

RM = rm

clean:	
	$(RM) -f *.o core malloc.cilk planner.cilk executor.cilk config.cilk \
              twiddle.cilk test_cilk time_cilk test_cilk.nd test_cilk.ndl \
	      putils.cilk rader.cilk \
              a.out *~ *.s *.bak 

distclean: clean
	$(RM) -f $(LIBFFTW_CILK) test_cilk time_cilk *.out

####################### Cilk Libraries, etc. ####################

# set this to point to the Cilk distribution directory.
CILK_DIR = $(HOME)/sandbox/cilk

# include makefile from distribution directory
-include $(CILK_DIR)/Makefile.common

# set this to a version of gcc that uses the GNU linker:
CC = gcc

# Set this to 1/2/3 to enable critical path measurements (slow)
CILK_CRITICAL_PATH=0

LIBS    = -lm $(CILK_RTS) $(CILK_LIB) $(ARCHLIBS)
HEADERS = $(CILK_HEADERS) $(FFTWDIR)/fftw.h fftw_cilk.cilkh
MY_CFLAGS = $(FFTW_INCLUDE) $(CFLAGS) -O6 -DFFTW_USING_CILK

###################### Building fftw_cilk #######################

AR = ar

config.cilk: $(FFTWDIR)/config.c
	(echo "#include <cilk.h>"; \
	cat $(FFTWDIR)/config.c) >config.cilk

malloc.cilk: $(FFTWDIR)/malloc.c
	cp $(FFTWDIR)/malloc.c malloc.cilk

twiddle.cilk: $(FFTWDIR)/twiddle.c
	cp $(FFTWDIR)/twiddle.c twiddle.cilk

planner.cilk: $(FFTWDIR)/planner.c
	cp $(FFTWDIR)/planner.c planner.cilk

putils.cilk: $(FFTWDIR)/putils.c
	cp $(FFTWDIR)/putils.c putils.cilk

rader.cilk: $(FFTWDIR)/rader.c
	cp $(FFTWDIR)/rader.c rader.cilk


%.o:	%.c $(HEADERS)  Makefile
	$(CC) $(MY_CFLAGS) -c $< -o $@

%.s:	%.c $(HEADERS)  Makefile
	$(CC) $(MY_CFLAGS) -S $< -o $@

%.c:	%.cilk $(HEADERS) $(CILK2C)  Makefile
	$(CILK2C) $(CILK2CFLAGS) -- $(MY_CFLAGS) -- $< -o $@

$(LIBFFTW_CILK): executor_cilk.o fftwnd_cilk.o $(OTHER_FFTW_OBJECTS)
	$(RM) -f $(LIBFFTW_CILK)
	$(AR) rv $(LIBFFTW_CILK) $(ALL_FFTW_CILK_OBJ)
	$(RANLIB) $(LIBFFTW_CILK)

install: $(LIBFFTW_CILK)
	$(CP) $(LIBFFTW_CILK) $(LIBDIR)
	$(RANLIB) $(LIBDIR)/$(LIBFFTW_CILK)
	$(CP) fftw_cilk.cilkh $(INCLUDEDIR)

test_cilk: test_cilk.o $(LIBFFTW_CILK) $(CILK_RTS) $(CILK_LIB)
	$(CC) $(LDFLAGS) test_cilk.o \
              $(LIBS) $(LIBFFTW_CILK) -o test_cilk

time_cilk: time_cilk.o $(LIBFFTW_CILK) $(CILK_RTS) $(CILK_LIB)
	$(CC) $(LDFLAGS) time_cilk.o \
              $(LIBS) $(LIBFFTW_CILK) -o time_cilk

################### Using the Nondeterminator ###################

ND_DIR = $(CILK_DIR)/nd
CILK2CFLAGS_ND = $(CILK2CFLAGS) -nd
CFLAGS_ND  = -I$(ND_DIR) $(MY_CFLAGS)
NDFLAGS =

CILK_RTS_ND = $(ND_DIR)/cilk-rts.nd.a
CILK_RTS_NDL = $(ND_DIR)/cilk-rts.ndl.a
LIBS_ND = -lm $(CILK_RTS_ND) $(CILK_LIB) $(ARCHLIBS)
LIBS_NDL = -lm $(CILK_RTS_NDL) $(CILK_LIB) $(ARCHLIBS)
HEADERS_ND = $(CILK_HEADERS) $(ND_DIR)/cilk-nd.h $(ND_DIR)/cilk-nd-cilk2c.h

# building normal version with the Nondeterminator
test_cilk.nd:   test_cilk.nd.o executor_cilk.nd.o fftwnd_cilk.nd.o \
                $(OTHER_FFTW_OBJECTS) $(CILK_RTS_ND) $(CILK_LIB)
	$(CC) $(LDFLAGS) test_cilk.nd.o executor_cilk.nd.o fftwnd_cilk.nd.o \
              $(LIBS_ND) $(FFTW_LIB) -o $@

%.nd.o: %.nd.c $(HEADERS_ND)
	$(CC) $(CFLAGS_ND) $(NDFLAGS) -DCILK_ND -c $< -o $@

%.nd.c: %.cilk $(HEADERS_ND) $(CILK2C)
	$(CILK2C) $(CILK2CFLAGS_ND) -- $(CFLAGS_ND) $(NDFLAGS)  \
                  -DCILK_ND -- $< -o $@

# building luxurious version with the Nondeterminator
test_cilk.ndl:  test_cilk.ndl.o executor_cilk.ndl.o fftwnd_cilk.ndl.o \
                $(OTHER_FFTW_OBJECTS) $(CILK_RTS_NDL) $(CILK_LIB)
	$(CC) $(LDFLAGS) test_cilk.ndl.o executor_cilk.ndl.o fftwnd_cilk.ndl.o\
              $(LIBS_NDL) $(FFTW_LIB) -o $@

%.ndl.o:        %.ndl.c $(HEADERS_ND)
	$(CC) $(CFLAGS_ND) $(NDFLAGS) -DCILK_NDL -c $< -o $@

%.ndl.c:        %.cilk $(HEADERS_ND)
	$(CILK2C) $(CILK2CFLAGS_ND) -- $(CFLAGS_ND) $(NDFLAGS) \
                  -DCILK_NDL -- $< -o $@
