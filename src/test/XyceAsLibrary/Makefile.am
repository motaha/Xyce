
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = $(DAKOTA_ARCH_DIR)/lib/Cparam.o 
endif

XYCELIBTESTSOURCES = \
  $(srcdir)/XyceLibTest.C 

# standalone XyceLibTest executable
check_PROGRAMS = XyceLibTest  testHarness1685
XyceLibTest_SOURCES = $(XYCELIBTESTSOURCES)
XyceLibTest_LDADD = $(top_builddir)/src/libxyce.la
XyceLibTest_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
 
BUG1685TESTSOURCES = \
  $(srcdir)/testHarness1685.C 

# standalone testHarness1685  executable
testHarness1685_SOURCES = $(BUG1685TESTSOURCES)
testHarness1685_LDADD = $(top_builddir)/src/libxyce.la
testHarness1685_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)
 
