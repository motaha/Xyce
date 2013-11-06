
AM_CPPFLAGS = @Xyce_INCS@ -DXYCEBUILDTIMESTAMP=\"`date +%Y%m%d%H%M`\"

noinst_LTLIBRARIES = libUtility.la

if BUILD_INTEL_FFT_SUPPORT
intelfft_SOURCES = \
  $(srcdir)/src/N_UTL_IntelFFT_Interface.C \
  $(srcdir)/include/N_UTL_IntelFFT_Interface.hpp
else
intelfft_SOURCES =
endif

if BUILD_FFTW_SUPPORT
fftw_SOURCES = \
  $(srcdir)/src/N_UTL_FFTW_Interface.C \
  $(srcdir)/include/N_UTL_FFTW_Interface.hpp
else
fftw_SOURCES =
endif

libUtility_la_SOURCES = \
  $(intelfft_SOURCES) \
  $(fftw_SOURCES) \
  $(srcdir)/src/N_UTL_BreakPoint.C \
  $(srcdir)/src/N_UTL_Demangle.C \
  $(srcdir)/src/N_UTL_NoCase.C \
  $(srcdir)/src/N_UTL_Expression.C \
  $(srcdir)/src/N_UTL_ExpressionData.C \
  $(srcdir)/src/N_UTL_ExpressionInternals.C \
  $(srcdir)/src/N_UTL_FFTInterface.C \
  $(srcdir)/src/N_UTL_LogStream.C \
  $(srcdir)/src/N_UTL_Misc.C \
  $(srcdir)/src/N_UTL_OptionBlock.C \
  $(srcdir)/src/N_UTL_Param.C \
  $(srcdir)/src/N_UTL_Registry.C \
  $(srcdir)/src/N_UTL_Timer.C \
  $(srcdir)/src/N_UTL_Version.C \
  $(srcdir)/include/N_UTL_Interface_Enum_Types.h \
  $(srcdir)/include/N_UTL_Expression.h \
  $(srcdir)/include/N_UTL_ExpressionInternals.h \
  $(srcdir)/include/N_UTL_ExpressionData.h \
  $(srcdir)/include/N_UTL_FFTInterface.hpp \
  $(srcdir)/include/N_UTL_Registry.h \
  $(srcdir)/include/N_UTL_Demangle.h \
  $(srcdir)/include/N_UTL_IndentStreamBuf.h \
  $(srcdir)/include/N_UTL_LogStream.h \
  $(srcdir)/include/N_UTL_FFTInterfaceDecl.hpp \
  $(srcdir)/include/N_UTL_FixedQueue.h \
  $(srcdir)/include/N_UTL_Functors.h \
  $(srcdir)/include/N_UTL_MachDepParams.h \
  $(srcdir)/include/N_UTL_Misc.h \
  $(srcdir)/include/N_UTL_BreakPoint.h \
  $(srcdir)/include/N_UTL_Factory.h \
  $(srcdir)/include/N_UTL_OptionBlock.h \
  $(srcdir)/include/N_UTL_Packable.h \
  $(srcdir)/include/N_UTL_Param.h \
  $(srcdir)/include/N_UTL_Timer.h \
  $(srcdir)/include/N_UTL_Graph.h \
  $(srcdir)/include/N_UTL_Xyce.h \
  $(srcdir)/include/N_UTL_IndentStreamBuf.h \
  $(srcdir)/include/N_UTL_NoCase.h \
  $(srcdir)/include/N_UTL_Version.h

BUILT_SOURCES = timestamp-generator

CLEANFILES = timestamp.h

timestamp-generator:
	echo "/* TIMESTAMP APPLIED */" > timestamp.h


