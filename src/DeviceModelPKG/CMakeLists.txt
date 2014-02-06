

# -- process files for this library ------------------------------------------

# class source list
set ( DeviceModel_SOURCES   
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_CSW.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Device.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceBld.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_DeviceLevelKey.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_Factory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B1.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_MOSFET_B2.C
      CACHE INTERNAL "X_S" )
  

# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
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
add_subdirectory( Core )
add_subdirectory( NeuronModels )
add_subdirectory( OpenModels )
add_subdirectory( TCADModels )



