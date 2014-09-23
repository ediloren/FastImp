#
# These statements are included in several of the 
# makefiles in the src/* directories. 
#
CFLAG = -c $(GDB_FLAG) $(OPTIM_FLAGS_FOR_C++) $(WARNING_FLAGS) \
	   $(CONSTRUCTOR) \
           $(PRINT_PROGRESS) $(DEBUG_FLAG) $(PRINT_RUN_TIME) \
	   $(MEMORY_REPORT) $(DEBUG_ELEMENT) $(DEBUG_GRID) \
	   $(PERFORMANCE_METER) $(DISABLE_CONSISTENT_STENCIL) \
	   $(DISABLE_PIECEWISE_SCHEME) \
	   $(PRINT_ELEMENT)

IFLAG = -I. -I$(INC_DIR) -I$(FFTW_INC_DIR)

COMPILER = $(C++)
ARCH	:= $(shell $(UTIL_DIR)/config.guess)

all: $(ARCH) $(ARCH)/$(MODULE).o $(DEPEND_FILE)

$(ARCH)/$(MODULE).o: $(OBJS:%=$(ARCH)/%)
	$(RM) $(ARCH)/$(MODULE).o; \
	$(LD) -r $(OBJS:%=$(ARCH)/%) -o $@

$(ARCH)/%.o: %.cc $(DEF_MAKEFILE)
	$(COMPILER) $(CFLAG) $(IFLAG) $< -o $@

$(ARCH): 
	if [ ! -d $(ARCH) ]; then\
	(mkdir $(ARCH); sleep 1;)\
	fi

etags:
	etags *.cc;

DEPEND_FILE = m.depends

$(DEPEND_FILE):
	$(COMPILER) -MM $(IFLAG) $(OBJS:%.o=%.cc) > $(DEPEND_FILE)
	$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
	$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)

depend: $(DEPEND_FILE)
	@echo $(DEPEND_FILE) has been generated;







