
# -- process files for this library ------------------------------------------

# class source list
set ( LoaderServices_SOURCES  
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LOA_CktLoader.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LOA_HBLoader.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LOA_Loader.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LOA_LoaderMgr.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LOA_NonlinearEquationLoader.C
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )




# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      LoaderServices
      CACHE INTERNAL "X_L" )

