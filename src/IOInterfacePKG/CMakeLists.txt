
# -- process files for this library ------------------------------------------

# class source list
set ( IOInterface_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_NetlistImportTool.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_CircuitBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_CircuitContext.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_CircuitMetadata.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_CmdParse.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_DeviceBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_DistributionTool.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_FunctionBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_mmio.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_Objective.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_OptionBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_OutputMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_Outputter.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_ParameterBlock.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_PkgOptionsMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_RestartMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_RestartNode.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_FourierMgr.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_IO_SpiceSeparatedFieldTool.C
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      IOInterface
      CACHE INTERNAL "X_L" )

# required source subdirectories
add_subdirectory( Output )
