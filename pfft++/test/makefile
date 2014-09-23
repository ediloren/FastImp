include ../def.makefile
ROOT_DIR = ..
LIB_DIR = ../lib
INC_DIR = ../inc
FFTW_INC_DIR = ../src/fftw-2.1.3/$(ARCH)/include
FFTW_LIB_DIR = ../src/fftw-2.1.3/$(ARCH)/lib
UTIL_DIR = ../util

Static = driver1
OBJS_Static = driver1.o

Dynamic1 = driver2
OBJS_Dynamic1 = driver2.o

Dynamic2 = driver3
OBJS_Dynamic2 = driver3.o

CFLAG = -c $(GDB_FLAG) $(OPTIM_FLAGS_FOR_C++) $(WARNING_FLAGS) \
	$(PRINT_PROGRESS) $(PRINT_RUN_TIME) \
	$(MEMORY_REPORT) $(DEBUG_ELEMENT) $(DEBUG_GRID) 

IFLAG = -I. -I$(INC_DIR) -I$(FFTW_INC_DIR)

ARCH := $(shell $(UTIL_DIR)/config.guess)

all: Static Dynamic1 Dynamic2 

Static: $(ARCH) $(DEPEND_FILE) $(OBJS_Static:%=$(ARCH)/%)
	${MAKE} -C $(ROOT_DIR);
	$(C++) \
	$(OBJS_Static:%=$(ARCH)/%) $(LIB_DIR)/$(ARCH)/pfft.o -o $(Static) -lm\
	$(FFTW_LIB_DIR)/libfftw.a $(FFTW_LIB_DIR)/librfftw.a \
	$(LIB_DIR)/$(ARCH)/clapack.a;  

Dynamic1: $(ARCH) $(DEPEND_FILE) $(OBJS_Dynamic1:%=$(ARCH)/%)
	${MAKE} -C $(ROOT_DIR);
	$(C++) \
	$(OBJS_Dynamic1:%=$(ARCH)/%) $(LIB_DIR)/$(ARCH)/pfft.o -o $(Dynamic1) -lm\
	$(FFTW_LIB_DIR)/libfftw.a $(FFTW_LIB_DIR)/librfftw.a \
	$(LIB_DIR)/$(ARCH)/clapack.a;  

Dynamic2: $(ARCH) $(DEPEND_FILE) $(OBJS_Dynamic2:%=$(ARCH)/%)
	${MAKE} -C $(ROOT_DIR);
	$(C++) \
	$(OBJS_Dynamic2:%=$(ARCH)/%) $(LIB_DIR)/$(ARCH)/pfft.o -o $(Dynamic2) -lm\
	$(FFTW_LIB_DIR)/libfftw.a $(FFTW_LIB_DIR)/librfftw.a \
	$(LIB_DIR)/$(ARCH)/clapack.a;  

clean:
	rm -rf $(ARCH);
	rm -f *~ core;
	rm -f *.dat
	rm -f $(Static) $(Dynamic1) $(Dynamic2) $(DebugDriver) $(TestAll);

$(ARCH)/%.o: %.cc 
	$(C++) $(CFLAG) $(IFLAG) $< -o $@

$(ARCH):
	if [ ! -d $(ARCH) ]; then\
	(mkdir $(ARCH); sleep 1;)\
	fi

DEPEND_FILE = m.depends

$(DEPEND_FILE):
		$(C++) -MM $(IFLAG) $(OBJS:%.o=./%.cc) $(Static).cc > $(DEPEND_FILE)
		$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
		$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)

depend: $(DEPEND_FILE)
	@echo $(DEPEND_FILE) has been generated;

-include $(DEPEND_FILE)




