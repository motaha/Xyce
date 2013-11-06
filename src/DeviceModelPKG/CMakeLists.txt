

# -- process generated files for this library ---------------------------------

if ( Xyce_REACTION_PARSER )

  # TODO: Use BISON_TARGET and FLEX_TARGET here.

  # create source and header files
  add_custom_command ( OUTPUT ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionParser.cxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionParser.hxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/location.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/position.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/stack.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionLexer.cxx

                      COMMAND ${BISON_EXECUTABLE}  -o ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionParser.cxx
                              ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/src/N_DEV_ReactionParser.yxx

                      COMMAND ${FLEX_EXECUTABLE} -+ -o${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionLexer.cxx
                              ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/src/N_DEV_ReactionLexer.l

                      DEPENDS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/src/N_DEV_ReactionParser.yxx
                              ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/src/N_DEV_ReactionLexer.l ) 

  # force (re-)generation at build time            
  add_custom_target ( reaction_parser_sources 
                      DEPENDS ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionParser.cxx 
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionParser.hxx
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/location.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/position.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/stack.hh
                              ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionLexer.cxx )
  
  # add to list of subpackage target dependencies
  set ( XYCE_TARGET_DEPENDENCIES ${XYCE_TARGET_DEPENDENCIES}
        "reaction_parser_sources" 
        CACHE INTERNAL "X_T_D" )
  
  # save source location and header paths
  set ( REACTION_PARSER_SRC ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionParser.cxx
                            ${Xyce_BINARY_DIR}/src/DeviceModelPKG/N_DEV_ReactionLexer.cxx )

  set ( REACTION_PARSER_INC ${Xyce_BINARY_DIR}/src/DeviceModelPKG )

  # pass list of generated files into global scope
  set ( XYCE_GENERATED_SOURCES ${XYCE_GENERATED_SOURCES}
        ${REACTION_PARSER_SRC}
        CACHE INTERNAL "X_G_S" )

endif ( Xyce_REACTION_PARSER )



# -- process files for this library ------------------------------------------

# class source list
set ( DeviceModel_SOURCES   
      ${REACTION_PARSER_SRC}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceInterface.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ACC.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_BJT.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_CSW.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Capacitor.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DAC.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ADC.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Device.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceInstance.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceState.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceModel.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceSensitivities.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceBld.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MatrixLoadData.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_NumericalJacobian.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceEntity.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_CompositeParam.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Bsrc.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Diode.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ISRC.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Inductor.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_JFET.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_LTRA.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_LTRA_Faddeeva.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MESFET.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET1.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET2.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET3.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET6.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B1.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B2.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B3.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B4.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B3SOI.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MutIndLin.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MutIndNonLin.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MutIndNonLin2.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_NewDevice.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Xygra.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_OpAmp.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Digital.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_VDMOS.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Resistor.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Resistor3.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ThermalResistor.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_SW.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Source.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_SourceData.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_TRA.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_VCCS.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Vcvs.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Vsrc.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ROM.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_RxnSet.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Region.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_RegionData.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_SpecieSource.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_PlaceHolder.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceSupport.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Param.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceOptions.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_SolverState.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Reaction.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ReactionNetwork.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_RateConstantCalculators.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceLevelKey.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Pars.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_OutputPars.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Units.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Factory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_RegisterDevices.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Configuration.C
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${REACTION_PARSER_INC}
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      DeviceModel 
      CACHE INTERNAL "X_L" )


# -- process remaining CMakeLists.txt files ----------------------------------

# optional source subdirectories

if ( Xyce_EXTDEV )
  add_subdirectory( EXTSC )
endif ( Xyce_EXTDEV )

if ( Xyce_RAD_MODELS )
  add_subdirectory( SandiaModels )
endif ( Xyce_RAD_MODELS )

if ( Xyce_NONFREE_MODELS )
  add_subdirectory( Xyce_NonFree )
endif ( Xyce_NONFREE_MODELS )


# required source subdirectories
add_subdirectory( ADMS )
add_subdirectory( AD )
add_subdirectory( PDESC )



