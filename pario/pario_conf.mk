# Makefie of parallel IO engine pario

# *.cpp files searching path
vpath %.cpp pario/src pario/src/utils pario/src/xml \
$(DYGILA_DIR)/utils/pioInit $(DYGILA_DIR)/utils/pioInit/action \
$(DYGILA_DIR)/utils/pioStream $(DYGILA_DIR)/utils/pioInit/containerReserve

# Add Ascent include path, linder flags
# of path and binary libs into
# prerequisites & recipes of dyGiLa
# building and linking rule
USE_PARIO := ON
ifeq ($(USE_PARIO), ON)
  include $(ASCENT_DIR)/share/ascent/ascent_config.mk
  APP_OPTS += $(ASCENT_INCLUDE_FLAGS)
  LDFLAGS  += $(ASCENT_LINK_RPATH)
  ifeq ($(ARCH), macos)
    # ascent_config.mk omits the rpath of the zlib that HDF5 links to
    LDFLAGS += -Wl,-rpath,$(wildcard $(ASCENT_DIR)/../zlib-*/lib)
  endif
  ifneq ($(filter $(ARCH),lumi macos),)
    LDLIBS += $(ASCENT_MPI_LIB_FLAGS)
  else
    ifneq ($(filter $(ARCH),lumi-hip-CC roihu-gpu),)
      LDLIBS += $(ASCENT_MPI_CUDA_LIB_FLAGS)
    endif
  endif
endif

# add headers searching directories
APP_OPTS += -I $(DYGILA_DIR)/pario/inc

# Object Path
PARIO_PATH = build/Targets

# pario objects, built by HILA pattern rules
PARIO_OBJECTS = $(PARIO_PATH)/xdmf.o     \
                $(PARIO_PATH)/xml_Amatrix.o \
                $(PARIO_PATH)/xml_pMarker.o \
                $(PARIO_PATH)/xml_massCurrent.o \
                $(PARIO_PATH)/xml_spinCurrent.o \
                $(PARIO_PATH)/pstream.o  \
                $(PARIO_PATH)/init.o     \
                $(PARIO_PATH)/shutdown.o \
                $(PARIO_PATH)/mesh.o     \
                $(PARIO_PATH)/mesh_gapA_FEDensity.o   \
                $(PARIO_PATH)/mesh_insitu_Temperature.o \
                $(PARIO_PATH)/mesh_phaseMarker.o \
                $(PARIO_PATH)/mesh_U1_3phi.o \
                $(PARIO_PATH)/mesh_lVec.o \
                $(PARIO_PATH)/mesh_GradientPhiVec.o \
                $(PARIO_PATH)/mesh_lVecSq.o \
                $(PARIO_PATH)/mesh_massCurrent.o      \
                $(PARIO_PATH)/mesh_spinCurrent.o      \
                $(PARIO_PATH)/mesh_AMatrix.o          \
                $(PARIO_PATH)/mesh_addGhost_verify.o  \
                $(PARIO_PATH)/actions_insitu.o  \
                $(PARIO_PATH)/actions_massCurrent.o \
                $(PARIO_PATH)/actions_spinCurrent.o \
		        $(PARIO_PATH)/actions_spinCurr_exaslice.o \
                $(PARIO_PATH)/actions_AMatrix.o     \
                $(PARIO_PATH)/actions_AMatrix_exaslice.o \
                $(PARIO_PATH)/actions_phaseMarker.o \
                $(PARIO_PATH)/actions_lVec.o \
                $(PARIO_PATH)/actions_GradientPhiVec.o \
                $(PARIO_PATH)/actions_GradientPhiVec_exaslice.o \
                $(PARIO_PATH)/actions_GradientPhiVec_clip.o \
                $(PARIO_PATH)/actions_GradientPhiVec_boxclip.o \
                $(PARIO_PATH)/actions_printTree.o   \
                $(PARIO_PATH)/ghostMask.o           \
                $(PARIO_PATH)/massCurrent.o         \
                $(PARIO_PATH)/spinCurrent.o         \
                $(PARIO_PATH)/Amatrix.o             \
                $(PARIO_PATH)/lVec.o \
                $(PARIO_PATH)/GradientPhiVec.o \
                $(PARIO_PATH)/gapAFETemPMarkerU1lVecSq.o  \
                $(PARIO_PATH)/U1PhaseStreaming.o \
                $(PARIO_PATH)/lVectorStreaming.o \
                $(PARIO_PATH)/GradientPhiVectorStreaming.o

.PHONY: pario
pario: $(PARIO_OBJECTS)
