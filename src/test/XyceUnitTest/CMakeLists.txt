
# -- build targets -----------------------------------------------------------

set ( XyceUnitTest_SOURCES 
      ${CMAKE_CURRENT_SOURCE_DIR}/XyceUnitTest.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/BlockMatrix_Helpers.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/HarmonicBalance.C
      ${CMAKE_CURRENT_SOURCE_DIR}/FFT.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/FFT_Helpers.C
      ${CMAKE_CURRENT_SOURCE_DIR}/HB_Builder_Helpers.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/IO_OutputMgr.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/ModelEvaluator.C )

include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )


# create binary
add_executable( XyceUnitTest 
                EXCLUDE_FROM_ALL
                ${XyceUnitTest_SOURCES} )

# link against available Xyce library 
if ( Xyce_ENABLE_SHARED )
  target_link_libraries( XyceUnitTest lib_xyce_shared ${DAKOTA_OBJS} )
else ( Xyce_ENABLE_SHARED )
  target_link_libraries( XyceUnitTest lib_xyce_static ${DAKOTA_OBJS} )
endif ( Xyce_ENABLE_SHARED )
