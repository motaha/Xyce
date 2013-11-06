
AM_CPPFLAGS = @Xyce_INCS@

if EXTDEV
  EXTDIRS = PDESC EXTSC AD ADMS
else
  EXTDIRS = PDESC AD ADMS
endif

#if RADMODELS
#  SUBDIRS = $(EXTDIRS) SandiaModels
#else
#  SUBDIRS = $(EXTDIRS)
#endif

#if NONFREEMODELS
#  SUBDIRS += Xyce_NonFree
#endif

SUBDIRS=$(EXTDIRS) $(MAYBE_RAD) $(MAYBE_NONFREE)

noinst_LTLIBRARIES = libDeviceModel.la

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

libDeviceModel_la_SOURCES = \
  $(srcdir)/src/N_DEV_DeviceInterface.C \
  $(srcdir)/src/N_DEV_ACC.C \
  $(srcdir)/src/N_DEV_BJT.C \
  $(srcdir)/src/N_DEV_CSW.C \
  $(srcdir)/src/N_DEV_Capacitor.C \
  $(srcdir)/src/N_DEV_DAC.C \
  $(srcdir)/src/N_DEV_ADC.C \
  $(srcdir)/src/N_DEV_Device.C \
  $(srcdir)/src/N_DEV_DeviceInstance.C \
  $(srcdir)/src/N_DEV_DeviceBlock.C \
  $(srcdir)/src/N_DEV_DeviceState.C \
  $(srcdir)/src/N_DEV_DeviceModel.C \
  $(srcdir)/src/N_DEV_DeviceMgr.C \
  $(srcdir)/src/N_DEV_DeviceSensitivities.C \
  $(srcdir)/src/N_DEV_DeviceBld.C \
  $(srcdir)/src/N_DEV_MatrixLoadData.C \
  $(srcdir)/src/N_DEV_NumericalJacobian.C \
  $(srcdir)/src/N_DEV_DeviceEntity.C \
  $(srcdir)/src/N_DEV_CompositeParam.C \
  $(srcdir)/src/N_DEV_Bsrc.C \
  $(srcdir)/src/N_DEV_Diode.C  \
  $(srcdir)/src/N_DEV_ISRC.C \
  $(srcdir)/src/N_DEV_Inductor.C \
  $(srcdir)/src/N_DEV_JFET.C \
  $(srcdir)/src/N_DEV_LTRA.C \
  $(srcdir)/src/N_DEV_LTRA_Faddeeva.C \
  $(srcdir)/src/N_DEV_MESFET.C  \
  $(srcdir)/src/N_DEV_MOSFET1.C \
  $(srcdir)/src/N_DEV_MOSFET2.C \
  $(srcdir)/src/N_DEV_MOSFET3.C \
  $(srcdir)/src/N_DEV_MOSFET6.C \
  $(srcdir)/src/N_DEV_MOSFET_B1.C \
  $(srcdir)/src/N_DEV_MOSFET_B2.C \
  $(srcdir)/src/N_DEV_MOSFET_B3.C \
  $(srcdir)/src/N_DEV_MOSFET_B4.C \
  $(srcdir)/src/N_DEV_MOSFET_B3SOI.C \
  $(srcdir)/src/N_DEV_MutIndLin.C \
  $(srcdir)/src/N_DEV_MutIndNonLin.C \
  $(srcdir)/src/N_DEV_MutIndNonLin2.C \
  $(srcdir)/src/N_DEV_NewDevice.C \
  $(srcdir)/src/N_DEV_Xygra.C \
  $(srcdir)/src/N_DEV_OpAmp.C \
  $(srcdir)/src/N_DEV_Digital.C \
  $(srcdir)/src/N_DEV_VDMOS.C \
  $(srcdir)/src/N_DEV_Resistor.C \
  $(srcdir)/src/N_DEV_Resistor3.C \
  $(srcdir)/src/N_DEV_ThermalResistor.C \
  $(srcdir)/src/N_DEV_SW.C \
  $(srcdir)/src/N_DEV_Source.C \
  $(srcdir)/src/N_DEV_SourceData.C \
  $(srcdir)/src/N_DEV_TRA.C \
  $(srcdir)/src/N_DEV_VCCS.C \
  $(srcdir)/src/N_DEV_Vcvs.C \
  $(srcdir)/src/N_DEV_Vsrc.C \
  $(srcdir)/src/N_DEV_ROM.C \
  $(srcdir)/src/N_DEV_RxnSet.C \
  $(srcdir)/src/N_DEV_Region.C \
  $(srcdir)/src/N_DEV_RegionData.C \
  $(srcdir)/src/N_DEV_SpecieSource.C \
  $(srcdir)/src/N_DEV_PlaceHolder.C \
  $(srcdir)/src/N_DEV_DeviceSupport.C \
  $(srcdir)/src/N_DEV_Param.C \
  $(srcdir)/src/N_DEV_DeviceOptions.C \
  $(srcdir)/src/N_DEV_SolverState.C \
  $(srcdir)/src/N_DEV_Reaction.C \
  $(srcdir)/src/N_DEV_ReactionNetwork.C \
  $(srcdir)/src/N_DEV_RateConstantCalculators.C \
  $(srcdir)/src/N_DEV_DeviceLevelKey.C \
  $(srcdir)/src/N_DEV_Pars.C \
  $(srcdir)/src/N_DEV_OutputPars.C \
  $(srcdir)/src/N_DEV_Units.C \
  $(srcdir)/src/N_DEV_Factory.C \
  $(srcdir)/src/N_DEV_RegisterDevices.C \
  $(srcdir)/src/N_DEV_Configuration.C \
  $(srcdir)/include/N_DEV_fwd.h \
  $(srcdir)/include/N_DEV_DeviceLevelKey.h \
  $(srcdir)/include/N_DEV_Pars.h \
  $(srcdir)/include/N_DEV_OutputPars.h \
  $(srcdir)/include/N_DEV_Units.h \
  $(srcdir)/include/N_DEV_Factory.h \
  $(srcdir)/include/N_DEV_Configuration.h \
  $(srcdir)/include/N_DEV_Property.h \
  $(srcdir)/include/N_DEV_DeviceInterface.h \
  $(srcdir)/include/N_DEV_SolverState.h \
  $(srcdir)/include/N_DEV_ExternData.h \
  $(srcdir)/include/N_DEV_ExternalSimulationData.h \
  $(srcdir)/include/N_DEV_MatrixLoadData.h \
  $(srcdir)/include/N_DEV_NumericalJacobian.h \
  $(srcdir)/include/N_DEV_ACC.h \
  $(srcdir)/include/N_DEV_BJT.h \
  $(srcdir)/include/N_DEV_Bsrc.h \
  $(srcdir)/include/N_DEV_CSW.h \
  $(srcdir)/include/N_DEV_Capacitor.h \
  $(srcdir)/include/N_DEV_Const.h \
  $(srcdir)/include/N_DEV_DAC.h \
  $(srcdir)/include/N_DEV_ADC.h \
  $(srcdir)/include/N_DEV_Device.h \
  $(srcdir)/include/N_DEV_DeviceBld.h \
  $(srcdir)/include/N_DEV_DeviceBlock.h \
  $(srcdir)/include/N_DEV_DeviceInstance.h \
  $(srcdir)/include/N_DEV_DeviceMgr.h \
  $(srcdir)/include/N_DEV_DeviceSensitivities.h \
  $(srcdir)/include/N_DEV_DeviceModel.h \
  $(srcdir)/include/N_DEV_DeviceOptions.h \
  $(srcdir)/include/N_DEV_DeviceState.h \
  $(srcdir)/include/N_DEV_DeviceSupport.h \
  $(srcdir)/include/N_DEV_DeviceTemplate.h \
  $(srcdir)/include/N_DEV_DeviceEntity.h \
  $(srcdir)/include/N_DEV_CompositeParam.h \
  $(srcdir)/include/N_DEV_Diode.h \
  $(srcdir)/include/N_DEV_ISRC.h \
  $(srcdir)/include/N_DEV_Inductor.h \
  $(srcdir)/include/N_DEV_Interpolators.h \
  $(srcdir)/include/N_DEV_JFET.h \
  $(srcdir)/include/N_DEV_LTRA.h \
  $(srcdir)/include/N_DEV_LTRA_Faddeeva.h \
  $(srcdir)/include/N_DEV_MESFET.h \
  $(srcdir)/include/N_DEV_MOSFET1.h \
  $(srcdir)/include/N_DEV_MOSFET2.h \
  $(srcdir)/include/N_DEV_MOSFET3.h \
  $(srcdir)/include/N_DEV_MOSFET6.h \
  $(srcdir)/include/N_DEV_MOSFET_B1.h \
  $(srcdir)/include/N_DEV_MOSFET_B2.h \
  $(srcdir)/include/N_DEV_MOSFET_B3.h \
  $(srcdir)/include/N_DEV_MOSFET_B4.h \
  $(srcdir)/include/N_DEV_MOSFET_B3SOI.h \
  $(srcdir)/include/N_DEV_MutIndLin.h \
  $(srcdir)/include/N_DEV_MutIndNonLin.h \
  $(srcdir)/include/N_DEV_MutIndNonLin2.h \
  $(srcdir)/include/N_DEV_NewDevice.h \
  $(srcdir)/include/N_DEV_Xygra.h \
  $(srcdir)/include/N_DEV_OpAmp.h \
  $(srcdir)/include/N_DEV_Digital.h \
  $(srcdir)/include/N_DEV_VDMOS.h \
  $(srcdir)/include/N_DEV_Param.h \
  $(srcdir)/include/N_DEV_PlaceHolder.h \
  $(srcdir)/include/N_DEV_Resistor.h \
  $(srcdir)/include/N_DEV_Resistor3.h \
  $(srcdir)/include/N_DEV_ThermalResistor.h \
  $(srcdir)/include/N_DEV_SW.h \
  $(srcdir)/include/N_DEV_Source.h \
  $(srcdir)/include/N_DEV_SourceData.h \
  $(srcdir)/include/N_DEV_TRA.h \
  $(srcdir)/include/N_DEV_VCCS.h \
  $(srcdir)/include/N_DEV_Vcvs.h \
  $(srcdir)/include/N_DEV_Vsrc.h \
  $(srcdir)/include/N_DEV_ROM.h \
  $(srcdir)/include/N_DEV_RxnSet.h \
  $(srcdir)/include/N_DEV_Region.h \
  $(srcdir)/include/N_DEV_RegionData.h \
  $(srcdir)/include/N_DEV_SpecieSource.h \
  $(srcdir)/include/N_DEV_TransportHelper.h \
  $(srcdir)/include/N_DEV_ModelTypes.h \
  $(srcdir)/include/N_DEV_Specie.h \
  $(srcdir)/include/N_DEV_Reaction.h \
  $(srcdir)/include/N_DEV_ReactionNetwork.h \
  $(srcdir)/include/N_DEV_RateConstantCalculators.h

# these can't be in "SOURCES" above or automake screams about YACC.
EXTRA_DIST=$(srcdir)/src/N_DEV_ReactionLexer.l \
  $(srcdir)/src/N_DEV_ReactionParser.yxx

if BUILD_REACTION_PARSER
libDeviceModel_la_SOURCES += $(BUILT_SOURCES) $(PARSER_SRCS)
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

