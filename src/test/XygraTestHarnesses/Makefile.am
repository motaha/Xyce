
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = 
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

EXTRA_DIST = XygraFile.C

# standalone XyceLibTest executable
check_PROGRAMS = XygraSource
XygraSource_SOURCES = $(srcdir)/XygraSource.C
XygraSource_LDADD = $(top_builddir)/src/libxyce.la $(RADLD) $(NONFREELD)
XygraSource_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
