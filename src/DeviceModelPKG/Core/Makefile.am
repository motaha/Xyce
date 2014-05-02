AM_CPPFLAGS = @Xyce_INCS@

noinst_LTLIBRARIES = libDeviceCore.la

# these are files built by bison and flex
if BUILD_REACTION_PARSER
BUILT_SOURCES= N_DEV_ReactionParser.cxx \
            N_DEV_ReactionParser.hxx \
            N_DEV_ReactionLexer.cxx \
            location.hh \
            position.hh \
	    stack.hh
else
BUILT_SOURCES=
endif

CLEANFILES= $(BUILT_SOURCES)

PARSER_HEADER_DEPS =
if BUILD_REACTION_PARSER
PARSER_HEADER_DEPS += N_DEV_ReactionParser.hxx
endif

#these are extra parser-related files that we want to include conditionally
PARSER_SRCS= $(srcdir)/include/N_DEV_ReactionLexer.h

libDeviceCore_la_SOURCES = \
  $(srcdir)/src/N_DEV_RegisterDevices.C \
  $(srcdir)/src/N_DEV_Algorithm.C \
  $(srcdir)/src/N_DEV_Print.C \
  $(srcdir)/src/N_DEV_CompositeParam.C \
  $(srcdir)/src/N_DEV_Configuration.C \
  $(srcdir)/src/N_DEV_DeviceBlock.C \
  $(srcdir)/src/N_DEV_DeviceEntity.C \
  $(srcdir)/src/N_DEV_DeviceInstance.C \
  $(srcdir)/src/N_DEV_DeviceInterface.C \
  $(srcdir)/src/N_DEV_DeviceMgr.C \
  $(srcdir)/src/N_DEV_DeviceModel.C \
  $(srcdir)/src/N_DEV_DeviceOptions.C \
  $(srcdir)/src/N_DEV_DeviceSensitivities.C \
  $(srcdir)/src/N_DEV_DeviceState.C \
  $(srcdir)/src/N_DEV_DeviceSupport.C \
  $(srcdir)/src/N_DEV_DeviceMaster.C \
  $(srcdir)/src/N_DEV_Dump.C \
  $(srcdir)/src/N_DEV_MatrixLoadData.C \
  $(srcdir)/src/N_DEV_Message.C \
  $(srcdir)/src/N_DEV_NameLevelKey.C \
  $(srcdir)/src/N_DEV_NumericalJacobian.C \
  $(srcdir)/src/N_DEV_LaTexDoc.C \
  $(srcdir)/src/N_DEV_Param.C \
  $(srcdir)/src/N_DEV_Pars.C \
  $(srcdir)/src/N_DEV_RateConstantCalculators.C \
  $(srcdir)/src/N_DEV_Reaction.C \
  $(srcdir)/src/N_DEV_ReactionNetwork.C \
  $(srcdir)/src/N_DEV_Region.C \
  $(srcdir)/src/N_DEV_RegionData.C \
  $(srcdir)/src/N_DEV_RxnSet.C \
  $(srcdir)/src/N_DEV_SolverState.C \
  $(srcdir)/src/N_DEV_Source.C \
  $(srcdir)/src/N_DEV_SourceData.C \
  $(srcdir)/src/N_DEV_SpecieSource.C \
  $(srcdir)/src/N_DEV_Units.C \
  $(srcdir)/include/N_DEV_fwd.h \
  $(srcdir)/include/N_DEV_RegisterDevices.h \
  $(srcdir)/include/N_DEV_Algorithm.h \
  $(srcdir)/include/N_DEV_Print.h \
  $(srcdir)/include/N_DEV_CompositeParam.h \
  $(srcdir)/include/N_DEV_Configuration.h \
  $(srcdir)/include/N_DEV_Const.h \
  $(srcdir)/include/N_DEV_Device.h \
  $(srcdir)/include/N_DEV_DeviceBlock.h \
  $(srcdir)/include/N_DEV_DeviceEntity.h \
  $(srcdir)/include/N_DEV_DeviceInstance.h \
  $(srcdir)/include/N_DEV_DeviceInterface.h \
  $(srcdir)/include/N_DEV_DeviceMgr.h \
  $(srcdir)/include/N_DEV_DeviceModel.h \
  $(srcdir)/include/N_DEV_DeviceOptions.h \
  $(srcdir)/include/N_DEV_DeviceSensitivities.h \
  $(srcdir)/include/N_DEV_DeviceState.h \
  $(srcdir)/include/N_DEV_DeviceSupport.h \
  $(srcdir)/include/N_DEV_DeviceMaster.h \
  $(srcdir)/include/N_DEV_Dump.h \
  $(srcdir)/include/N_DEV_ExternData.h \
  $(srcdir)/include/N_DEV_ExternalSimulationData.h \
  $(srcdir)/include/N_DEV_Interpolators.h \
  $(srcdir)/include/N_DEV_MatrixLoadData.h \
  $(srcdir)/include/N_DEV_Message.h \
  $(srcdir)/include/N_DEV_NameLevelKey.h \
  $(srcdir)/include/N_DEV_NumericalJacobian.h \
  $(srcdir)/include/N_DEV_LaTexDoc.h \
  $(srcdir)/include/N_DEV_Param.h \
  $(srcdir)/include/N_DEV_Pars.h \
  $(srcdir)/include/N_DEV_RateConstantCalculators.h \
  $(srcdir)/include/N_DEV_Reaction.h \
  $(srcdir)/include/N_DEV_ReactionNetwork.h \
  $(srcdir)/include/N_DEV_Region.h \
  $(srcdir)/include/N_DEV_RegionData.h \
  $(srcdir)/include/N_DEV_RxnSet.h \
  $(srcdir)/include/N_DEV_SolverState.h \
  $(srcdir)/include/N_DEV_Source.h \
  $(srcdir)/include/N_DEV_SourceData.h \
  $(srcdir)/include/N_DEV_Specie.h \
  $(srcdir)/include/N_DEV_SpecieSource.h \
  $(srcdir)/include/N_DEV_TransportHelper.h \
  $(srcdir)/include/N_DEV_Units.h

# these can't be in "SOURCES" above or automake screams about YACC.
EXTRA_DIST=$(srcdir)/src/N_DEV_ReactionLexer.l \
  $(srcdir)/src/N_DEV_ReactionParser.yxx

if BUILD_REACTION_PARSER
libDeviceCore_la_SOURCES += $(BUILT_SOURCES) $(PARSER_SRCS)
endif

#We have to do it this way, because automake appears to assume it can use
# yacc and bison interchangably, and we can't.

N_DEV_ReactionParser.cxx: $(srcdir)/src/N_DEV_ReactionParser.yxx
	$(BISON) -o$@ $(srcdir)/src/N_DEV_ReactionParser.yxx

N_DEV_ReactionLexer.cxx: $(srcdir)/src/N_DEV_ReactionLexer.l
	$(LEX) -+ -o$@ $(srcdir)/src/N_DEV_ReactionLexer.l

N_DEV_ReactionParser.hxx: N_DEV_ReactionParser.cxx
location.hh: N_DEV_ReactionParser.cxx
position.hh: N_DEV_ReactionParser.cxx
stack.hh: N_DEV_ReactionParser.cxx

