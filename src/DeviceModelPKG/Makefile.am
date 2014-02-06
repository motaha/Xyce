
AM_CPPFLAGS = @Xyce_INCS@

if EXTDEV
  EXTDIRS = TCADModels EXTSC NeuronModels ADMS
else
  EXTDIRS = TCADModels NeuronModels ADMS
endif

#if RADMODELS
#  SUBDIRS = $(EXTDIRS) SandiaModels
#else
#  SUBDIRS = $(EXTDIRS)
#endif

#if NONFREEMODELS
#  SUBDIRS += Xyce_NonFree
#endif

SUBDIRS=Core OpenModels $(EXTDIRS) $(MAYBE_RAD) $(MAYBE_NONFREE)

noinst_LTLIBRARIES = libDeviceModel.la

libDeviceModel_la_SOURCES =  \
  $(srcdir)/src/N_DEV_CSW.C \
  $(srcdir)/src/N_DEV_Device.C \
  $(srcdir)/src/N_DEV_DeviceBld.C \
  $(srcdir)/src/N_DEV_DeviceLevelKey.C \
  $(srcdir)/src/N_DEV_Factory.C \
  $(srcdir)/src/N_DEV_MOSFET_B1.C \
  $(srcdir)/src/N_DEV_MOSFET_B2.C \
  $(srcdir)/include/N_DEV_CSW.h \
  $(srcdir)/include/N_DEV_DeviceBld.h \
  $(srcdir)/include/N_DEV_DeviceLevelKey.h \
  $(srcdir)/include/N_DEV_Factory.h \
  $(srcdir)/include/N_DEV_MOSFET_B1.h \
  $(srcdir)/include/N_DEV_MOSFET_B2.h \
  $(srcdir)/include/N_DEV_ModelTypes.h \
  $(srcdir)/include/N_DEV_PlaceHolder.h \
  $(srcdir)/include/N_DEV_Property.h


