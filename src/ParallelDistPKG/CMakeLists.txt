
# -- process files for this library ------------------------------------------

# class source list
set ( ParallelDist_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_CommFactory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_ParallelMachine.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_MPI.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_GlobalAccessor.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_Manager.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_Migrator.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_MPIComm.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_ParComm.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_ParDir.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_ParMap.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_ParMapFactory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_SerialComm.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_PDS_SerialParComm.C 
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      ParallelDist 
      CACHE INTERNAL "X_L" )
