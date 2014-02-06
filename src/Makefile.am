
AM_CPPFLAGS = @Xyce_INCS@

SUBDIRS = \
  AnalysisPKG \
  CircuitPKG \
  $(DAKOTA_DIR) \
  DeviceModelPKG \
  ErrorHandlingPKG \
  IOInterfacePKG \
  LinearAlgebraServicesPKG \
  LoaderServicesPKG \
  NonlinearSolverPKG \
  ParallelDistPKG \
  TimeIntegrationPKG \
  TopoManagerPKG \
  UtilityPKG \
  MultiTimePDEPKG \
  test

XYCESOURCES = \
  $(srcdir)/CircuitPKG/src/N_CIR_Xyce.C \
  $(srcdir)/CircuitPKG/include/N_CIR_Xyce.h

XYCELINK_BASE = \
  $(top_builddir)/src/CircuitPKG/libCircuit.la \
  $(top_builddir)/src/DeviceModelPKG/libDeviceModel.la \
  $(top_builddir)/src/DeviceModelPKG/Core/libDeviceCore.la \
  $(top_builddir)/src/DeviceModelPKG/OpenModels/libOpenModels.la \
  $(top_builddir)/src/DeviceModelPKG/TCADModels/libTCADModels.la \
  $(top_builddir)/src/DeviceModelPKG/NeuronModels/libNeuronModels.la \
  $(top_builddir)/src/DeviceModelPKG/ADMS/libADMS.la \
  $(top_builddir)/src/ErrorHandlingPKG/libErrorHandling.la \
  $(top_builddir)/src/IOInterfacePKG/libIOInterface.la \
  $(top_builddir)/src/IOInterfacePKG/Output/libOutput.la \
  $(top_builddir)/src/LinearAlgebraServicesPKG/libLinearAlgebraServices.la \
  $(top_builddir)/src/LoaderServicesPKG/libLoaderServices.la \
  $(top_builddir)/src/NonlinearSolverPKG/libNonlinearSolver.la \
  $(top_builddir)/src/ParallelDistPKG/libParallelDist.la \
  $(top_builddir)/src/AnalysisPKG/libAnalysis.la \
  $(top_builddir)/src/TimeIntegrationPKG/libTimeIntegration.la \
  $(top_builddir)/src/TopoManagerPKG/libTopoManager.la \
  $(top_builddir)/src/UtilityPKG/libUtility.la \
  $(top_builddir)/src/MultiTimePDEPKG/libMPDE.la

# conditionally link in external device extensions
if EXTDEV
    EXTLINK = $(XYCELINK_BASE) $(top_builddir)/src/DeviceModelPKG/EXTSC/libEXT.la
else
    EXTLINK = $(XYCELINK_BASE)
endif

# conditionally link in radiation-aware device models.
if RADMODELS
    RADLINK = $(EXTLINK) $(top_builddir)/src/DeviceModelPKG/SandiaModels/libSandiaModels.la
else
    RADLINK = $(EXTLINK)
endif

# conditionally link in radiation-aware device models.
if NONFREEMODELS
    NONFREELINK = $(RADLINK) $(top_builddir)/src/DeviceModelPKG/Xyce_NonFree/libNonFree.la
else
    NONFREELINK = $(RADLINK)
endif

# conditionally link in Dakota link materials
if DAKOTA_LINK
   DAKOTA_DIR = DakotaLinkPKG
   XYCELINK = $(NONFREELINK) $(top_builddir)/src/DakotaLinkPKG/libDakotaLink.la
else
   XYCELINK = $(NONFREELINK)
endif
if DAKOTA_OBJ_NEEDED
  DAKOTA_OBJS = $(DAKOTA_ARCH_DIR)/lib/Cparam.o
endif

# build library (static & shared) for linking (w/libtool) into other codes
lib_LTLIBRARIES = libxyce.la
libxyce_la_SOURCES = $(XYCESOURCES)
libxyce_la_LIBADD = $(XYCELINK)
libxyce_la_LDFLAGS = $(AM_LDFLAGS) @XYCELIBS@

# standalone Xyce executable (monolithic)
if BUILD_XYCE_BINARY
bin_PROGRAMS = Xyce
Xyce_SOURCES = $(srcdir)/CircuitPKG/src/Xyce.C
Xyce_LDADD = libxyce.la
Xyce_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

EXTRA_PROGRAMS = Xygra XygraFile

Xygra_SOURCES = $(srcdir)/CircuitPKG/src/Xygra.C
Xygra_LDADD = libxyce.la
Xygra_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

XygraFile_SOURCES = $(srcdir)/test/XygraTestHarnesses/XygraFile.C
XygraFile_LDADD = libxyce.la
XygraFile_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

endif

# install headers needed to compile other codes with libxyce.*
include_HEADERS = \
  $(srcdir)/CircuitPKG/include/N_CIR_Xyce.h \
  $(srcdir)/CircuitPKG/include/N_CIR_Xygra.h \
  $(srcdir)/UtilityPKG/include/N_UTL_Xyce.h \
  $(srcdir)/IOInterfacePKG/include/N_IO_CmdParse.h \
  $(srcdir)/IOInterfacePKG/include/N_IO_fwd.h \
  $(srcdir)/UtilityPKG/include/N_UTL_MachDepParams.h \
  $(srcdir)/UtilityPKG/include/N_UTL_Misc.h \
  $(srcdir)/UtilityPKG/include/N_UTL_fwd.h \
  $(srcdir)/ParallelDistPKG/include/N_PDS_Manager.h \
  $(srcdir)/DeviceModelPKG/Core/include/N_DEV_ExternalSimulationData.h \
  $(srcdir)/DeviceModelPKG/Core/include/N_DEV_fwd.h \
  $(srcdir)/TimeIntegrationPKG/include/N_TIA_TwoLevelError.h \
  $(top_builddir)/src/Xyce_config.h
