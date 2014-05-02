
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

