# -- process remaining CMakeLists.txt files ----------------------------------

add_subdirectory( Core )
add_subdirectory( OpenModels )

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
add_subdirectory( NeuronModels )
add_subdirectory( TCADModels )



