
# -- process files for this library ------------------------------------------

# class source list
set ( Circuit_SOURCES
      # DO NOT list Xyce.C here.  It is part of the Xyce binary target. 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_CIR_Xyce.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_CIR_Xygra.C
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )

# -- build requested libraries -----------------------------------------------
# libxyce is created at top level.  See ${Xyce_SOURCE_DIR}/src/CMakeLists.txt 
