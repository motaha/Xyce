
AM_CPPFLAGS = @Xyce_INCS@

# needed for Dakota 4.x not 5.0
if DAKOTA_OBJ_NEEDED 
  DAKOTA_OBJS = $(DAKOTA_ARCH_DIR)/lib/Cparam.o 
endif

XYCEUNITTESTSOURCES = \
  $(srcdir)/XyceUnitTest.C \
	$(srcdir)/BlockMatrix.C \
	$(srcdir)/BlockMatrix_Helpers.C \
	$(srcdir)/BlockMatrix_Helpers.h \
	$(srcdir)/HarmonicBalance.C \
	$(srcdir)/FFT.C \
	$(srcdir)/FFT_Helpers.C \
	$(srcdir)/FFT_Helpers.h \
	$(srcdir)/HB_Builder_Helpers.C \
	$(srcdir)/HB_Builder_Helpers.h \
	$(srcdir)/HB_Loader_Helpers.C \
	$(srcdir)/HB_Loader_Helpers.h \
	$(srcdir)/IO_OutputMgr.C \
	$(srcdir)/ModelEvaluator.C

# standalone XyceUnitTest executable
EXTRA_PROGRAMS = XyceUnitTest
XyceUnitTest_SOURCES = $(XYCEUNITTESTSOURCES)
if TRILINOS_DEV_RYTHMOS
XyceUnitTest_LDADD = $(top_builddir)/src/libxyce.la -lrythmos -lstratimikos -lstratimikosifpack -lstratimikosamesos -lstratimikosml -lstratimikosaztecoo -lthyraepetraext -lthyraepetra -lthyra -lrtop -lloca -lnox -lml -laztecoo -lifpack -lamesos -lsuperlu -lumfpack -lepetraext -lepetra -lteuchos -lamd -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread 
else
XyceUnitTest_LDADD = $(top_builddir)/src/libxyce.la 
endif
XyceUnitTest_LDFLAGS = -static $(AM_LDFLAGS) $(DAKOTA_OBJS)

