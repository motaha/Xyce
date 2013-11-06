
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = $(DAKOTA_ARCH_DIR)/lib/Cparam.o 
endif

EXTRA_DIST = XygraFile.C

# standalone XyceLibTest executable
check_PROGRAMS = XygraSource
XygraSource_SOURCES = $(srcdir)/XygraSource.C
XygraSource_LDADD = $(top_builddir)/src/libxyce.la
XygraSource_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
