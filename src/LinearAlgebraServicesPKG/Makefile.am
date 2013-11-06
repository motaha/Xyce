
AM_CPPFLAGS = @Xyce_INCS@

noinst_LTLIBRARIES = libLinearAlgebraServices.la

if BUILD_BELOS_SUPPORT
belos_SOURCES = \
  $(srcdir)/src/N_LAS_BelosSolver.C \
  $(srcdir)/include/N_LAS_BelosSolver.h
endif

if BUILD_ML_SUPPORT
ml_SOURCES = \
  $(srcdir)/include/N_LAS_MLPrecond.h \
  $(srcdir)/src/N_LAS_MLPrecond.C
endif

if BUILD_SHYLU_SUPPORT
shylu_SOURCES = \
  $(srcdir)/src/N_LAS_ShyLUSolver.C \
  $(srcdir)/include/N_LAS_ShyLUSolver.h
endif

if BUILD_KSPARSE
ksparse_SOURCES = \
  $(srcdir)/src/N_LAS_KSparseSolver.C \
  $(srcdir)/ksparse/src/Epetra_CrsKundertSparse.C \
  $(srcdir)/ksparse/src/alloc.c \
  $(srcdir)/ksparse/src/driver.c \
  $(srcdir)/ksparse/src/smplinkrows.c \
  $(srcdir)/ksparse/src/spalloc.c \
  $(srcdir)/ksparse/src/spbuild.c \
  $(srcdir)/ksparse/src/spfactor.c \
  $(srcdir)/ksparse/src/spoutput.c \
  $(srcdir)/ksparse/src/spsmp.c \
  $(srcdir)/ksparse/src/spsolve.c \
  $(srcdir)/ksparse/src/sputils.c \
  $(srcdir)/include/N_LAS_KSparseSolver.h \
  $(srcdir)/ksparse/include/Epetra_CrsKundertSparse.h \
  $(srcdir)/ksparse/include/capabil.h \
  $(srcdir)/ksparse/include/complex.h \
  $(srcdir)/ksparse/include/config.h \
  $(srcdir)/ksparse/include/hw.h \
  $(srcdir)/ksparse/include/iferrmsg.h \
  $(srcdir)/ksparse/include/misc.h \
  $(srcdir)/ksparse/include/smpdefs.h \
  $(srcdir)/ksparse/include/spconfig.h \
  $(srcdir)/ksparse/include/spdefs.h \
  $(srcdir)/ksparse/include/sperror.h \
  $(srcdir)/ksparse/include/spice.h \
  $(srcdir)/ksparse/include/spmat.h \
  $(srcdir)/ksparse/include/spmatrix.h \
  $(srcdir)/ksparse/include/sppars.h \
  $(srcdir)/ksparse/include/strext.h \
  $(srcdir)/ksparse/include/util.h
endif

libLinearAlgebraServices_la_SOURCES = \
  $(belos_SOURCES) \
  $(ml_SOURCES) \
  $(shylu_SOURCES) \
  $(ksparse_SOURCES) \
  $(srcdir)/src/N_LAS_SimpleSolver.C \
  $(srcdir)/src/N_LAS_AmesosSolver.C \
  $(srcdir)/src/N_LAS_AztecOOSolver.C \
  $(srcdir)/src/N_LAS_Problem.C \
  $(srcdir)/src/N_LAS_MultiVector.C \
  $(srcdir)/src/N_LAS_Matrix.C \
  $(srcdir)/src/N_LAS_MOROperators.C \
  $(srcdir)/src/N_LAS_System.C \
  $(srcdir)/src/N_LAS_Builder.C \
  $(srcdir)/src/N_LAS_HBBuilder.C \
  $(srcdir)/src/N_LAS_Vector.C \
  $(srcdir)/src/N_LAS_TransformTool.C \
  $(srcdir)/src/N_LAS_LAFactory.C \
  $(srcdir)/src/N_LAS_SolverFactory.C \
  $(srcdir)/src/N_LAS_BlockVector.C \
  $(srcdir)/src/N_LAS_BlockMatrix.C \
  $(srcdir)/src/N_LAS_BlockSystemHelpers.C \
  $(srcdir)/src/N_LAS_IfpackPrecond.C \
  $(srcdir)/src/N_LAS_HBBlockJacobiPrecond.C \
  $(srcdir)/src/N_LAS_HBBlockJacobiEpetraOperator.C \
  $(srcdir)/src/N_LAS_HBPrecondFactory.C \
  $(srcdir)/src/N_LAS_TrilinosPrecondFactory.C \
  $(srcdir)/include/N_LAS_SimpleSolver.h \
  $(srcdir)/include/N_LAS_AmesosSolver.h \
  $(srcdir)/include/N_LAS_AztecOOSolver.h \
  $(srcdir)/include/N_LAS_Builder.h \
  $(srcdir)/include/N_LAS_HBBuilder.h \
  $(srcdir)/include/N_LAS_Solver.h \
  $(srcdir)/include/N_LAS_Problem.h \
  $(srcdir)/include/N_LAS_Matrix.h \
  $(srcdir)/include/N_LAS_MOROperators.h \
  $(srcdir)/include/N_LAS_MultiVector.h \
  $(srcdir)/include/N_LAS_QueryUtil.h \
  $(srcdir)/include/N_LAS_System.h \
  $(srcdir)/include/N_LAS_Vector.h \
  $(srcdir)/include/N_LAS_TransformTool.h \
  $(srcdir)/include/N_LAS_LAFactory.h \
  $(srcdir)/include/N_LAS_SolverFactory.h \
  $(srcdir)/include/N_LAS_BlockVector.h \
  $(srcdir)/include/N_LAS_BlockMatrix.h \
  $(srcdir)/include/N_LAS_BlockSystemHelpers.h \
  $(srcdir)/include/N_LAS_Preconditioner.h \
  $(srcdir)/include/N_LAS_NoPrecond.h \
  $(srcdir)/include/N_LAS_IfpackPrecond.h \
  $(srcdir)/include/N_LAS_HBBlockJacobiPrecond.h \
  $(srcdir)/include/N_LAS_HBBlockJacobiEpetraOperator.h \
  $(srcdir)/include/N_LAS_PrecondFactory.h \
  $(srcdir)/include/N_LAS_HBPrecondFactory.h \
  $(srcdir)/include/N_LAS_TrilinosPrecondFactory.h \
  $(srcdir)/include/N_LAS_Util.h
