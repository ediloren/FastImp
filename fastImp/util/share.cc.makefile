#const static char cvsid[] = "$Id: share.cc.makefile,v 1.20 2003/07/16 15:32:34 zhzhu Exp $";
#
# These statements are included in several of the 
# makefiles in the src/* directories. 
#
CFLAG = -c $(GPROF_FLAG) $(GDB_FLAG) $(OPTIM_FLAGS_FOR_C++) $(CONSTRUCTOR) \
	   $(WARNING_FLAGS) $(MORE_WARNING_FLAGS) \
           $(DEBUG_FLAG) $(SURF_RUN_TIME) $(SURF_PROGRESS) $(SURF_MEMORY) \
	   $(PFFT_RUN_TIME) $(PFFT_PROGRESS) \
	   $(DEBUG_ELEMENT) $(DEBUG_GRID) $(MEMORY_REPORT) \
           $(DEBUG_PRE_COND) $(DEBUG_ITERATIVE) $(DEBUG_DIRECT) $(DEBUG_CURRENT) \
	   $(DEBUG_GMRES) $(PRINT_GMRES_RESID) $(DEBUG_MEMORY) \
	   $(DEBUG_BUFFER) $(DEBUG_CONTACT) \
	   $(DIAG_LOCAL_PRECOND) $(USE_FIXED_BUFFER_SIZE) \
	   $(DEBUG_MESH_RECT_SPIRAL) $(DEBUG_EMQS) $(DEBUG_RHO) \
	   $(USE_MULTIPLE_SINE) $(HIGH_ACCURACY) $(NO_GRAY_WINDOW) \
	   $(LARGE_STEP_RATIO) $(DEBUG_TRANSLINE) \
	   $(DISABLE_SCALING) $(DISABLE_CONSISTENT_STENCIL) \
	   $(DISABLE_PIECEWISE_SCHEME) $(PRINT_ELEMENT)
 
IFLAG = -I. -I$(INC_DIR) -I$(PFFT_INC_DIR) -I$(SUPER_LU_INC_DIR)

COMPILER = $(C++)
ARCH	:= $(shell $(UTIL_DIR)/config.guess)

all: $(ARCH) $(ARCH)/$(MODULE).o $(DEPEND_FILE)
$(ARCH)/$(MODULE).o: $(OBJS:%=$(ARCH)/%)
	$(RM) $(ARCH)/$(MODULE).o; \
	$(LD) -r $(OBJS:%=$(ARCH)/%) -o $@
$(ARCH)/%.o: %.cc $(DEF_MAKEFILE)
#$(ARCH)/%.o: %.cc 
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
	etags *.cc;

DEPEND_FILE = m.depends

$(DEPEND_FILE):
	$(COMPILER) -MM $(IFLAG) $(OBJS:%.o=%.cc) > $(DEPEND_FILE)
	$(SED) -e "s,\(.*\.o\),$(ARCH:%=%/\1)," $(DEPEND_FILE) > $(DEPEND_FILE).foo
	$(MV) $(DEPEND_FILE).foo $(DEPEND_FILE)

depend: $(DEPEND_FILE)
	@echo $(DEPEND_FILE) has been generated;

clean:
	-$(RM) $(BIN_DIR)/core; 
	-$(RM) $(LIB_DIR)/$(MODULE).a
	-$(RMR) $(ARCH);
	-$(RM) *~ ; 
	-$(RM) $(DEPEND_FILE); \

minorclean:
	-$(RM) *~;

-include $(DEPEND_FILE)
