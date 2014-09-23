include ../../def.makefile

INC_DIR = ../../inc
UTIL_DIR = ../../util
LIB_DIR = ../../lib

#CFLAG = -c -g -ansi -pedantic -DDEBUG 
CFLAGS = -ansi -pedantic $(OPTIM_FLAGS_FOR_C) $(PERFORMANCE_METER)
IFLAG = -I. -I$(INC_DIR)

COMPILER = $(CC)

###############################
# clapack double precision files
LIST1 = dgesvd.o                                                                       \
        dbdsqr.o  dgeqr2.o  dlamch.o  dlarfg.o  dlaset.o  dlasr.o   dorgbr.o  dormbr.o \
        dgebd2.o  dgeqrf.o  dlange.o  dlarft.o  dlasq1.o  dlasrt.o  dorgl2.o  dorml2.o \
        dgebrd.o  dlapy2.o  dlartg.o  dlasq2.o  dlassq.o  dorglq.o  dormlq.o           \
        dgelq2.o  dlabrd.o  dlarf.o   dlas2.o   dlasq3.o  dlasv2.o  dorgqr.o  dormqr.o \
        dgelqf.o  dlacpy.o  dlarfb.o  dlascl.o  dlasq4.o  dorg2r.o  dorm2r.o

# cblas files:
LIST2 = dcopy.o  dgemv.o  dnrm2.o  dscal.o  dtrmm.o  lsame.o  \
        dgemm.o  dger.o   drot.o   dswap.o  dtrmv.o  xerbla.o 

# clapack util files:
LIST3 = ilaenv.o

# f2clibs files:
LIST4 = d_sign.o  i_dnnt.o  pow_dd.o  pow_di.o  s_cat.o  s_cmp.o  s_copy.o

OBJS = $(LIST1) $(LIST2) $(LIST3) $(LIST4) 

ARCH	:= $(shell $(UTIL_DIR)/config.guess)

all: $(LIB_DIR)/$(ARCH)/clapack.a $(DEPEND_FILE)

$(LIB_DIR)/$(ARCH)/clapack.a: $(ARCH) $(OBJS:%=$(ARCH)/%)
	$(AR) ru $@ $(OBJS:%=$(ARCH)/%)
	$(RANLIB) $@

$(ARCH)/%.o: %.c
	$(CC) $(CFLAGS) -c $(IFLAG) $< -o $@

$(ARCH): 
	if [ ! -d $(ARCH) ]; then\
	(mkdir $(ARCH); sleep 1;)\
	fi

DEPEND_FILE = m.depends
$(DEPEND_FILE):
		$(COMPILER) -MM $(IFLAG) $(OBJS:%.o=%.c) > $(DEPEND_FILE)
		$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
		$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)
clean:
	-$(RM) *~; 
	-${RM} $(LIB_DIR)/*.a;
	-$(RM) -rf $(ARCH); 

minorclean:
	$(RM) *~ *.o; \

etags:
	etags *.c *.h;

depend: $(DEPEND_FILE)

include $(DEPEND_FILE)









