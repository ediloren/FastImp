#
# These statements are included in several of the 
# makefiles in the src/* directories. 
#
CFLAG = -c $(GPROF_FLAG) $(GDB_FLAG) $(OPTIM_FLAGS_FOR_C) $(WARNING_FLAGS) $(DEBUG_FLAG) \
	   $(FOR_COMPAQ_1)
IFLAG = -I. -I$(INC_DIR)

COMPILER = $(CC)
DEPEND_COMPILER = $(GCC)

ARCH	:= $(shell $(UTIL_DIR)/config.guess)

#build module.o
all:: $(ARCH) $(ARCH)/$(MODULE).o $(DEPEND_FILE)
$(ARCH)/$(MODULE).o: $(OBJS:%=$(ARCH)/%)
	$(RM) $(ARCH)/$(MODULE).o; \
	$(LD) -r $(OBJS:%=$(ARCH)/%) -o $@
$(ARCH)/%.o: %.c
	$(COMPILER) $(CFLAG) $(IFLAG) $< -o $@

#build *.a
lib:: $(ARCH) $(LIB_DIR)/$(MODULE).a $(DEPEND_FILE)
$(LIB_DIR)/$(MODULE).a: $(OBJS:%=$(ARCH)/%)
	ar ru $@ $(OBJS:%=$(ARCH)/%)
	$(RANLIB) $@


$(ARCH): 
	if [ ! -d $(ARCH) ]; then\
	(mkdir $(ARCH); sleep 1;)\
	fi

etags:
	etags *.c *.h;

DEPEND_FILE = m.depends
$(DEPEND_FILE):
	$(DEPEND_COMPILER) -MM $(IFLAG) $(OBJS:%.o=%.c) > $(DEPEND_FILE)
	$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
	$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)

depend: $(DEPEND_FILE)
	@echo $(DEPEND_FILE) has been generated;

clean:
	-$(RM) $(BIN_DIR)/core; 
	-$(RM) $(LIB_DIR)/$(MODULE).a
	-$(RMR) $(ARCH);
	-$(RM) *~ ;
	-$(RM) $(DEPEND_FILE);

minorclean:
	-$(RM) *~;

-include $(DEPEND_FILE)
