
# -- process files for this library ------------------------------------------

# class source list
set ( TopoManager_SOURCES 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktGraph.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktGraphCreatorBasic.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktGraphBasic.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktGraphSupport.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktNode.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktNode_Ckt.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktNode_Dev.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktNode_V.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_CktNodeCreator.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_DevInsertionTool.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_Directory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_Indexor.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_Node.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_NodeBlock.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_NodeDevBlock.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_ParNode.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_Topology.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_TopologyMgr.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_TOP_TopoLSUtil.C 
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      TopoManager 
      CACHE INTERNAL "X_L" )

