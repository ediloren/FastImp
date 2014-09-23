include ./def.makefile

SURF_ROOT =  $(shell pwd)

PFFT_ROOT =  ../pfft++
PFFT_LIB_DIR = $(PFFT_ROOT)/lib

SURF_SRC_DIR    = $(SURF_ROOT)/src/surf
SURF_TEST_DIR   = $(SURF_ROOT)/test
SURF_UTIL       = $(SURF_ROOT)/util
SURF_LIB_DIR    = $(SURF_ROOT)/lib
SURF_INC_DIR    = $(SURF_ROOT)/inc
SURF_BIN_DIR    = $(SURF_ROOT)/bin

FILEIO_SRC_DIR    = $(SURF_ROOT)/src/fileIo
FILEIO_INC_DIR    = $(SURF_INC_DIR)

MESH_SRC_DIR    = $(SURF_ROOT)/src/mesh
MESH_INC_DIR    = $(SURF_INC_DIR)

SUPER_LU_SRC_DIR    = $(SURF_ROOT)/src/SuperLU
SUPER_LU_LIB_DIR    = $(SUPER_LU_SRC_DIR)

FFTW_LIB_DIR = $(PFFT_ROOT)/src/fftw-2.1.3/$(ARCH)/lib

ARCH	:= $(shell $(SURF_UTIL)/config.guess)

all:	fastImp

fastImp: $(SURF_BIN_DIR)/$(ARCH) 
	${MAKE} -C $(FILEIO_SRC_DIR); 
	${MAKE} -C $(MESH_SRC_DIR); 
	${MAKE} -C $(SURF_SRC_DIR); 
	${MAKE} -C $(PFFT_ROOT); 
	${MAKE} -C $(SUPER_LU_SRC_DIR) lib;
	$(C++) $(GPROF_FLAG) \
	$(SURF_SRC_DIR)/$(ARCH)/surfModule.o \
	$(FILEIO_SRC_DIR)/$(ARCH)/surfFileIo.o \
	$(MESH_SRC_DIR)/$(ARCH)/surfMesh.o \
	$(PFFT_LIB_DIR)/$(ARCH)/pfft.o \
	-o $(SURF_BIN_DIR)/$(ARCH)/fastimp -lm \
	$(FFTW_LIB_DIR)/libfftw.a \
	$(FFTW_LIB_DIR)/librfftw.a \
	$(SUPER_LU_LIB_DIR)/superlu.a \
	$(PFFT_LIB_DIR)/$(ARCH)/clapack.a \
	$(SUPER_LU_LIB_DIR)/blas.a;

$(SURF_BIN_DIR)/$(ARCH): 
	if [ ! -d $(SURF_BIN_DIR)/$(ARCH) ]; then\
	(mkdir $(SURF_BIN_DIR)/$(ARCH); sleep 1;)\
	fi

depend:
	$(MAKE) clean;
	${MAKE} -C $(SURF_SRC_DIR) depend; 
	${MAKE} -C $(MESH_SRC_DIR) depend; 
	${MAKE} -C $(FILEIO_SRC_DIR) depend; 

minorclean:
	${RM} *~; \
	${MAKE} -C $(SURF_SRC_DIR) minorclean; 
	${MAKE} -C $(MESH_SRC_DIR) minorclean; 
	${MAKE} -C $(FILEIO_SRC_DIR) minorclean; 
	${RM} $(SURF_INC_DIR)/*~; \
	${RM} $(SURF_TEST_DIR)/core; \
	${RM} $(SURF_UTIL)/*~; \

clean:	
	${RM} *~; \
	${MAKE} -C $(SURF_SRC_DIR) clean; 
	${MAKE} -C $(MESH_SRC_DIR) clean; 
	${MAKE} -C $(FILEIO_SRC_DIR) clean; 
	${RMR} $(SURF_BIN_DIR)/$(ARCH); \
	${RM} $(SURF_TEST_DIR)/core; \
	${RM} $(SURF_INC_DIR)/*~; \
	${RM} $(SURF_UTIL)/*~; \

realclean:
	${RM} *~; \
	${RMR} $(SURF_BIN_DIR)/$(ARCH); \
	${RM} $(SURF_TEST_DIR)/core; \
	${RM} $(SURF_INC_DIR)/*~; \
	${RM} $(SURF_UTIL)/*~; \
	${MAKE} -C $(SURF_SRC_DIR) clean; 
	${MAKE} -C $(MESH_SRC_DIR) clean; 
	${MAKE} -C $(FILEIO_SRC_DIR) clean; 
	${MAKE} -C $(SUPER_LU_SRC_DIR) clean; 
	${MAKE} -C $(SUPER_LU_SRC_DIR)/CBLAS clean;

tar:
	$(MAKE) realclean; \
	$(RM) ./bak/*.tgz; \
	$(RM) ./*.tgz; \
	$(RM) ./*.tar; \
	tar czvf fastImp.tgz .; \
	mv fastImp.tgz ./bak/.; \
	ls -al ./bak; \






