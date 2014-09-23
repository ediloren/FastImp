include ../../def.makefile

MODULE = pfft

OBJS = element.o stencil.o gridIndex.o grid.o gridElement.o interpMat.o \
       projectMat.o cubature.o g2gPoint.o \
       kernelIntegration.o utils.o \

include $(UTIL_DIR)/share.cc.makefile

clean:
	$(RM) $(BIN_DIR)/core; 
	$(RM) *~ ; 
	$(RM) $(DEPEND_FILE); 
	$(RMR) $(ARCH); 
	${MAKE} -C ./test_vector3D clean; \
	${MAKE} -C ./test_element clean; \
	${MAKE} -C ./test_sparse clean; \
	${MAKE} -C ./test_spVec clean; \
	${MAKE} -C ./test_calcp clean; \

minorclean:
	$(RM) *~;
	${MAKE} -C ./test_vector3D clean; \
	${MAKE} -C ./test_element clean; \
	${MAKE} -C ./test_sparse clean; \
	${MAKE} -C ./test_spVec clean; \

-include $(DEPEND_FILE)


