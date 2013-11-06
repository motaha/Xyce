
# -- process files for this library ------------------------------------------

# class source list  
set ( Analysis_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_AnalysisBase.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_AnalysisInterface.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_AnalysisManager.C  
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_Dakota.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_DCSweep.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_HB.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_AC.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_MOR.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_InnerLevelSolve.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_ModelEvaluator.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_ModelEvaluator_Stateless.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_MPDE.C    
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_OutputMgrAdapter.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_Step.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_SweepParam.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_ANP_Transient.C 
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )




# -- build requested libraries ------------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      Analysis 
      CACHE INTERNAL "X_L" )

