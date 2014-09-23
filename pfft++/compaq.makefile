include ./def.makefile

PFFT_ROOT =  $(shell pwd)

PFFT_SRC_DIR    = $(PFFT_ROOT)/src/pfft
PFFT_INC_DIR    = $(PFFT_ROOT)/inc
PFFT_BIN_DIR    = $(PFFT_ROOT)/bin
PFFT_TEST_DIR   = $(PFFT_ROOT)/test
PFFT_UTIL       = $(PFFT_ROOT)/util
PFFT_LIB_DIR    = $(PFFT_ROOT)/lib

CLAPACK_SRC_DIR = $(PFFT_ROOT)/src/clapack
CLAPACK_INC_DIR = $(PFFT_ROOT)/inc
CLAPACK_LIB_DIR = $(PFFT_ROOT)/lib

DENSE_SRC_DIR    = $(PFFT_ROOT)/src/dense
DENSE_INC_DIR    = $(PFFT_ROOT)/inc
DENSE_BIN_DIR    = $(PFFT_ROOT)/bin

FFTW_ROOT	= $(PFFT_ROOT)/src/fftw-2.1.3/$(ARCH)
FFTW_SRC_DIR 	= $(FFTW_ROOT)
FFTW_INC_DIR 	= $(FFTW_ROOT)/include
FFTW_LIB_DIR 	= $(FFTW_ROOT)/lib

DOC_DIR = $(PFFT_ROOT)/doc
TEST_DIR = $(PFFT_ROOT)/test

TEST_CALCP_DIR = $(PFFT_SRC_DIR)/test_calcp
TEST_ELEMENT_DIR = $(PFFT_SRC_DIR)/test_element
TEST_SPVEC_DIR = $(PFFT_SRC_DIR)/test_spVec
TEST_SPARSE_DIR = $(PFFT_SRC_DIR)/test_sparse
TEST_VECTOR3D_DIR = $(PFFT_SRC_DIR)/test_vector3D

ARCH	:= $(shell $(PFFT_UTIL)/config.guess)

all:  $(PFFT_LIB_DIR)/$(ARCH) $(FFTW_LIB_DIR)/libfftw.a
	@echo making ...;
	${MAKE} -C $(PFFT_SRC_DIR);
	${MAKE} -C $(CLAPACK_SRC_DIR); 
	${MAKE} -C $(DENSE_SRC_DIR); 
	$(RM) $(PFFT_LIB_DIR)/pfft.o; \
	$(LD) -r \
	$(PFFT_SRC_DIR)/$(ARCH)/pfft.o \
	$(DENSE_SRC_DIR)/$(ARCH)/pfftDense.o \
	$(PFFT_LIB_DIR)/$(ARCH)/clapack.a \
	$(FFTW_LIB_DIR)/libfftw.a \
	$(FFTW_LIB_DIR)/librfftw.a \
	-o $(PFFT_LIB_DIR)/$(ARCH)/pfft.o 

$(FFTW_LIB_DIR)/libfftw.a: $(FFTW_ROOT)
	cd $(FFTW_ROOT); ../configure --prefix=$(FFTW_ROOT);
	${MAKE} -C $(FFTW_SRC_DIR);
	${MAKE} -C $(FFTW_SRC_DIR) install;

$(FFTW_ROOT): 
	if [ ! -d $(FFTW_ROOT) ]; then \
	(mkdir $(FFTW_ROOT); sleep 1;) \
	fi 

$(PFFT_LIB_DIR)/$(ARCH):
	if [ ! -d $(PFFT_LIB_DIR)/$(ARCH) ]; then \
	(mkdir $(PFFT_LIB_DIR)/$(ARCH); sleep 1;) \
	fi

tests: 
	${MAKE} -C $(TEST_CALCP_DIR);
	${MAKE} -C $(TEST_ELEMENT_DIR);
	${MAKE} -C $(TEST_SPARSE_DIR);
	${MAKE} -C $(TEST_SPVEC_DIR);
	${MAKE} -C $(TEST_VECTOR3D_DIR);

$(PFFT_LIB_DIR)/$(ARCH)/pfft.o: all

fftw:  $(FFTW_ROOT) 
	cd $(FFTW_ROOT); ../configure --prefix=$(FFTW_ROOT);
	${MAKE} -C $(FFTW_SRC_DIR);
	${MAKE} -C $(FFTW_SRC_DIR) install;

depend:
	${MAKE} -C $(PFFT_SRC_DIR) depend; 
	${MAKE} -C $(CLAPACK_SRC_DIR) depend; 
	${MAKE} -C $(DENSE_SRC_DIR) depend; 

minorclean:
	${RM} *~; \
	${MAKE} -C $(PFFT_SRC_DIR) minorclean; 
	${MAKE} -C $(CLAPACK_SRC_DIR) minorclean; 
	${MAKE} -C $(DENSE_SRC_DIR) minorclean; 
	${RM} $(PFFT_INC_DIR)/*~; \
	${RM} $(PFFT_TEST_DIR)/core; \
	${RM} $(PFFT_UTIL)/*~; \

clean:	
	${RM} *~; \
	${MAKE} -C $(PFFT_SRC_DIR) clean; 
	${MAKE} -C $(CLAPACK_SRC_DIR) clean; 
	${MAKE} -C $(DENSE_SRC_DIR) clean; 
	-${MAKE} -C $(FFTW_ROOT) distclean; 
	${MAKE} -C $(TEST_CALCP_DIR) clean;
	${MAKE} -C $(TEST_ELEMENT_DIR) clean;
	${MAKE} -C $(TEST_SPARSE_DIR) clean;
	${MAKE} -C $(TEST_SPVEC_DIR) clean;
	${MAKE} -C $(TEST_VECTOR3D_DIR) clean;
	${MAKE} -C $(PFFT_TEST_DIR) clean;
	${RMR} $(PFFT_LIB_DIR)/$(ARCH); \
	${RM} $(PFFT_TEST_DIR)/core; \
	${RM} $(PFFT_INC_DIR)/*~; \
	${RM} $(PFFT_UTIL)/*~; \

realclean:
	${RM} *~; \
	${MAKE} -C $(PFFT_SRC_DIR) clean; 
	${MAKE} -C $(CLAPACK_SRC_DIR) clean; 
	${MAKE} -C $(DENSE_SRC_DIR) clean; 
	-${MAKE} -C $(FFTW_ROOT) distclean;
	${MAKE} -C $(TEST_CALCP_DIR) clean;
	${MAKE} -C $(TEST_ELEMENT_DIR) clean;
	${MAKE} -C $(TEST_SPARSE_DIR) clean;
	${MAKE} -C $(TEST_SPVEC_DIR) clean;
	${MAKE} -C $(TEST_VECTOR3D_DIR) clean;
	${MAKE} -C $(PFFT_TEST_DIR) clean;
	${RMR} $(PFFT_LIB_DIR)/*; \
	${RM} $(PFFT_TEST_DIR)/core; \
	${RM} $(PFFT_INC_DIR)/*~; \
	${RM} $(PFFT_UTIL)/*~; \

tar:
	$(MAKE) realclean; \
	$(RM) ./bak/*.tgz; \
	$(RM) ./*.tgz; \
	$(RM) ./*.tar; \
	tar czvf pfft++.tgz .; \
	mv pfft++.tgz ./bak/.; \
	ls -al ./bak; \
