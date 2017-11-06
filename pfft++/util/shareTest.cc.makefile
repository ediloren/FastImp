INC_DIR = ../../../inc
FFTW_INC_DIR = ../../../src/fftw-2.1.3/$(ARCH)/include
UTIL_DIR = ../../../util

CFLAG = -c $(GDB_FLAG) $(OPTIM_FLAGS) $(WARNING_FLAGS) \
	$(CONSTRUCTOR) $(MCMODEL) $(PIC) \
        $(PRINT_PROGRESS) $(DEBUG_FLAG) $(PRINT_RUN_TIME) \
	$(MEMORY_REPORT) $(DEBUG_ELEMENT) $(DEBUG_GRID) \
	$(PERFORMANCE_METER)  $(DISABLE_CONSISTENT_STENCIL) \
	$(DISABLE_PIECEWISE_SCHEME) $(PRINT_ELEMENT) \
	$(TEST_FREESPACE)

IFLAG = -I. -I$(INC_DIR) -I$(FFTW_INC_DIR)

COMPILER = $(C++)
ARCH	:= $(shell $(UTIL_DIR)/config.guess)

all: $(ARCH) $(DEPEND_FILE) $(ARCH)/$(MODULE) 
	-rm -f $(MODULE); 
	-ln -s $(ARCH)/$(MODULE) $(MODULE);
	
clean:
	-rm -rf $(ARCH);
	-rm -f $(MODULE) *~ core;
	
$(ARCH)/$(MODULE): $(OBJS:%=$(ARCH)/%) $(ARCH)/$(MODULE).o
	$(COMPILER) -o $(ARCH)/$(MODULE) $(OBJS:%=$(ARCH)/%) $(ARCH)/$(MODULE).o -lm

$(ARCH)/%.o: ../%.cc
	$(COMPILER) $(CFLAG) $(IFLAG) $< -o $@

$(ARCH)/$(MODULE).o: $(MODULE).cc 	
	$(COMPILER) $(CFLAG) $(IFLAG) $< -o $@

$(ARCH):
	if [ ! -d $(ARCH) ]; then\
	(mkdir $(ARCH); sleep 1;)\
	fi

etags:
	etags *.cc *.h;

DEPEND_FILE = m.depends

$(DEPEND_FILE):
		$(COMPILER) -MM $(IFLAG) $(OBJS:%.o=../%.cc) \
		$(MODULE).cc > $(DEPEND_FILE)
		$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
		$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)

depend: $(DEPEND_FILE)
	@echo $(DEPEND_FILE) has been generated;

-include $(DEPEND_FILE)






