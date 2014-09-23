#
# These statements are included in several of the 
# makefiles in the src/* directories. 
#
CFLAG = -c $(DEBUG_FLAGS) $(OPTIM_FLAGS_FOR_C) $(WARNING_FLAGS) \
           $(DEBUG_PRINT_FLAG) $(RUN_TIME_PRINT_FLAG) \
	   $(PERFORMANCE_METER)

IFLAG = -I. -I$(INC_DIR)

COMPILER = $(CC)

ARCH	:= $(shell $(UTIL_DIR)/config.guess)

all:: $(ARCH) $(ARCH)/$(MODULE).o $(DEPEND_FILE)

$(ARCH)/$(MODULE).o: $(OBJS:%=$(ARCH)/%)
	$(RM) $(ARCH)/$(MODULE).o; \
	$(LD) -r $(OBJS:%=$(ARCH)/%) -o $@

$(ARCH)/%.o: %.c
	$(COMPILER) $(CFLAG) $(IFLAG) $< -o $@

$(ARCH): 
	if [ ! -d $(ARCH) ]; then\
	(mkdir $(ARCH); sleep 1;)\
	fi

etags:
	etags *.c *.h;

DEPEND_FILE = m.depends
$(DEPEND_FILE):
	$(COMPILER) -MM $(IFLAG) $(OBJS:%.o=%.c) > $(DEPEND_FILE)
	$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
	$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)

depend: $(DEPEND_FILE)
	@echo $(DEPEND_FILE) has been generated;


clean:
	-$(RM) $(BIN_DIR)/core; 
	-$(RMR) $(ARCH);
	-$(RM) $(DEPEND_FILE);
	-$(RM) *~ ;

minorclean:
	-$(RM) *~;

-include $(DEPEND_FILE)







