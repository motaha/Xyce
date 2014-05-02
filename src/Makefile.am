
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

XYCELINK = \
  $(top_builddir)/src/CircuitPKG/libCircuit.la \
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
    EXTLINK = $(top_builddir)/src/DeviceModelPKG/EXTSC/libEXT.la
else
    EXTLINK =
endif

# conditionally link in radiation-aware device models.
if RADMODELS
    RADLD = $(top_builddir)/src/DeviceModelPKG/SandiaModels/libSandiaModels.la
else
    RADLD =
endif

# conditionally link in radiation-aware device models.
if NONFREEMODELS
    NONFREELD = $(top_builddir)/src/DeviceModelPKG/Xyce_NonFree/libNonFree.la
else
    NONFREELD =
endif

# conditionally link in Dakota link materials
if DAKOTA_LINK
   DAKOTA_DIR = DakotaLinkPKG
   DAKOTALINK = $(top_builddir)/src/DakotaLinkPKG/libDakotaLink.la
else
   DAKOTALINK =
endif
if DAKOTA_OBJ_NEEDED
  DAKOTA_OBJS =
endif

# build library (static & shared) for linking (w/libtool) into other codes
lib_LTLIBRARIES = libxyce.la
libxyce_la_SOURCES = $(XYCESOURCES)
libxyce_la_LIBADD = $(EXTLINK) $(DAKOTALINK) $(XYCELINK)
libxyce_la_LDFLAGS = $(AM_LDFLAGS) @XYCELIBS@

if BUILD_XYCE_BINARY

if BUILD_XYCE_SHAREABLE

# standalone Xyce executable (shareable object)
bin_PROGRAMS = Xyce
Xyce_SOURCES = $(srcdir)/CircuitPKG/src/Xyce.C
Xyce_LDADD = libxyce.la $(NONFREELD) $(RADLD)
Xyce_LDFLAGS = $(AM_LDFLAGS) $(DAKOTA_OBJS)

EXTRA_PROGRAMS = Xygra XygraFile

Xygra_SOURCES = $(srcdir)/CircuitPKG/src/Xygra.C
Xygra_LDADD = -L.libs -lxyce
Xygra_LDFLAGS = $(AM_LDFLAGS) $(DAKOTA_OBJS)

XygraFile_SOURCES = $(srcdir)/test/XygraTestHarnesses/XygraFile.C
XygraFile_LDADD = -L.libs -lxyce
XygraFile_LDFLAGS = $(AM_LDFLAGS) $(DAKOTA_OBJS)

else

# standalone Xyce executable (monolithic)
bin_PROGRAMS = Xyce
Xyce_SOURCES = $(srcdir)/CircuitPKG/src/Xyce.C
Xyce_LDADD = libxyce.la $(NONFREELD) $(RADLD)
Xyce_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

EXTRA_PROGRAMS = Xygra XygraFile

Xygra_SOURCES = $(srcdir)/CircuitPKG/src/Xygra.C
Xygra_LDADD = libxyce.la $(NONFREELD) $(RADLD)
Xygra_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

XygraFile_SOURCES = $(srcdir)/test/XygraTestHarnesses/XygraFile.C
XygraFile_LDADD = libxyce.la $(NONFREELD) $(RADLD)
XygraFile_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

endif

endif

# install headers needed to compile other codes with libxyce.*
include_HEADERS = \
  $(srcdir)/AnalysisPKG/include/N_ANP_fwd.h \
  $(srcdir)/CircuitPKG/include/N_CIR_Xyce.h \
  $(srcdir)/CircuitPKG/include/N_CIR_Xygra.h \
  $(srcdir)/DeviceModelPKG/Core/include/N_DEV_ExternalSimulationData.h \
  $(srcdir)/DeviceModelPKG/Core/include/N_DEV_fwd.h \
  $(srcdir)/DeviceModelPKG/Core/include/N_DEV_NameLevelKey.h \
  $(srcdir)/ErrorHandlingPKG/include/N_ERH_fwd.h \
  $(srcdir)/IOInterfacePKG/include/N_IO_CmdParse.h \
  $(srcdir)/IOInterfacePKG/include/N_IO_fwd.h \
  $(srcdir)/ParallelDistPKG/include/N_PDS_Manager.h \
  $(srcdir)/ParallelDistPKG/include/N_PDS_fwd.h \
  $(srcdir)/TimeIntegrationPKG/include/N_TIA_TwoLevelError.h \
  $(srcdir)/TopoManagerPKG/include/N_TOP_fwd.h \
  $(srcdir)/UtilityPKG/include/N_UTL_MachDepParams.h \
  $(srcdir)/UtilityPKG/include/N_UTL_Misc.h \
  $(srcdir)/UtilityPKG/include/N_UTL_NameLevelKey.h \
  $(srcdir)/UtilityPKG/include/N_UTL_NoCase.h \
  $(srcdir)/UtilityPKG/include/N_UTL_ReportHandler.h \
  $(srcdir)/UtilityPKG/include/N_UTL_TypeIndex.h \
  $(srcdir)/UtilityPKG/include/N_UTL_Xyce.h \
  $(srcdir)/UtilityPKG/include/N_UTL_fwd.h \
  $(top_builddir)/src/Xyce_config.h
