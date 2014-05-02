
# -- process files for this library ------------------------------------------

# class source list
set ( ADMS_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_RegisterADMSDevices.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ADMSvbic.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ADMSPSP103VA.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ADMSHBT_X.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_DEV_ADMSbsimcmg.C
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      CACHE INTERNAL "X_I" )

# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES}
      ADMS
      CACHE INTERNAL "X_L" )
