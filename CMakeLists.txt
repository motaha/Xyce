
# ----------------------------------------------------------------------------
# CMake Setup          
# ----------------------------------------------------------------------------

# -- configure the CMake environment -----------------------------------------

# set CMake minimum version
cmake_minimum_required ( VERSION 2.6 )


# set makefile verbosity level 
set ( CMAKE_VERBOSE_MAKEFILE ON )


# load CMake macros
include ( CheckIncludeFileCXX )
include ( CheckFortranFunctionExists )
include ( CheckFunctionExists )
include ( CheckLibraryExists )
include ( CheckSymbolExists )


# Xyce library macro 
macro ( XYCE_FIND_LIBRARY LIBNAME )
  message ( STATUS "Looking for library ${LIBNAME}" )
  find_library ( LIBRARY_${LIBNAME} 
                 NAMES ${LIBNAME} lib${LIBNAME}
                 PATHS ${CMAKE_LIBRARY_PATH} )

  if ( LIBRARY_${LIBNAME} )
    message ( STATUS "Looking for library ${LIBNAME} - found" )

    # DEBUG:  /link flag broken on win32; use full paths
    if ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBRARY_${LIBNAME}} ${ARCH_LIBRARIES} )
    else ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBNAME} ${ARCH_LIBRARIES} )
    endif ( CMAKE_HOST_WIN32 )
  else ( LIBRARY_${LIBNAME} )
    message ( FATAL_ERROR "Looking for library ${LIBNAME} - not found" )
  endif ( LIBRARY_${LIBNAME} )

endmacro ( XYCE_FIND_LIBRARY LIBNAME )

macro ( XYCE_FIND_OPTIONAL_LIBRARY LIBNAME )
  message ( STATUS "Looking for library ${LIBNAME}" )
  find_library ( LIBRARY_${LIBNAME} 
                 NAMES ${LIBNAME} lib${LIBNAME}
                 PATHS ${CMAKE_LIBRARY_PATH} )

  if ( LIBRARY_${LIBNAME} )
    message ( STATUS "Looking for library ${LIBNAME} - found" )

    # DEBUG:  /link flag broken on win32; use full paths
    if ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBRARY_${LIBNAME}} ${ARCH_LIBRARIES} )
    else ( CMAKE_HOST_WIN32 )
      set ( ARCH_LIBRARIES ${LIBNAME} ${ARCH_LIBRARIES} )
    endif ( CMAKE_HOST_WIN32 )
  endif ( LIBRARY_${LIBNAME} )

endmacro ( XYCE_FIND_OPTIONAL_LIBRARY LIBNAME )

set(Xyce_CMAKE_OPTIONS "" CACHE FILEPATH "List of option files to load")

# load configuration files;  cmd line flags override;  no trap for missing file
foreach (CONFIG_OPTS_FILE ${Xyce_CMAKE_OPTIONS})
  message ( STATUS "Reading in configuration options from ${CONFIG_OPTS_FILE} ..." )
  include ( ${CONFIG_OPTS_FILE} )
endforeach (CONFIG_OPTS_FILE ${Xyce_CMAKE_OPTIONS})

# All of the default options are gathered in one file.
message ( STATUS "Loading ${CMAKE_HOME_DIRECTORY}/defaultConfig.cmake" )
include ( "${CMAKE_HOME_DIRECTORY}/defaultConfig.cmake" )



# -- process build option dependencies; notify user of configuration changes--

# check for consistent library and binary flags                         

if ( Xyce_STATIC_BINARY AND CMAKE_HOST_APPLE )
  message ( STATUS "STATIC_BINARY is not supported on OSX.  Disabling static linkage." )
  set ( Xyce_STATIC_BINARY OFF )
endif ( Xyce_STATIC_BINARY AND CMAKE_HOST_APPLE )

if ( Xyce_STATIC_BINARY AND NOT Xyce_ENABLE_STATIC ) 
  message ( STATUS "STATIC_BINARY requires ENABLE_STATIC.  Enabling static libaries." )
  set ( Xyce_ENABLE_STATIC ON ) 
endif ( Xyce_STATIC_BINARY AND NOT Xyce_ENABLE_STATIC ) 

if ( Xyce_STATIC_BINARY AND Xyce_DISABLE_BINARY )
  message ( STATUS "STATIC_BINARY conflicts with DISABLE_BINARY.  Enabling binary." )
  set ( Xyce_DISABLE_BINARY OFF )
endif ( Xyce_STATIC_BINARY AND Xyce_DISABLE_BINARY )

if ( Xyce_DISABLE_BINARY AND NOT Xyce_ENABLE_SHARED AND NOT Xyce_ENABLE_STATIC )
  message ( STATUS "No targets are enabled.  Resetting defaults." )
  set ( Xyce_DISABLE_BINARY OFF )
  set ( Xyce_ENABLE_STATIC ON )
endif ( Xyce_DISABLE_BINARY AND NOT Xyce_ENABLE_SHARED AND NOT Xyce_ENABLE_STATIC )

set ( Xyce_BELOS ON )
set ( Xyce_USE_ISORROPIA ON )
set ( Xyce_USE_ZOLTAN ON  )

# disable options that conflict with MPI 
if ( Xyce_SUPERLUDIST AND NOT Xyce_PARALLEL_MPI )
  message ( STATUS "SUPERLUDIST requires PARALLEL_MPI.  Disabling SUPERLUDIST" )      
  set ( Xyce_SUPERLUDIST OFF )
endif ( Xyce_SUPERLUDIST AND NOT Xyce_PARALLEL_MPI )


# disable Sandia models if the source directory is missing
if ( Xyce_RAD_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels ) 
  message ( STATUS "SandiaModels directory cannot be found.  Disabling RAD_MODELS" )      
  set ( Xyce_RAD_MODELS OFF )  
endif ( Xyce_RAD_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels ) 

# disable non-free models if the source directory is missing
if ( Xyce_NONFREE_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree ) 
  message ( STATUS "Xyce_NonFree directory cannot be found.  Disabling NONFREE_MODELS" )      
  set ( Xyce_NONFREE_MODELS OFF )  
endif ( Xyce_NONFREE_MODELS AND NOT EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/Xyce_NonFree ) 




# -- prepare CMake environment for system inspection -------------------------

# enable Fortran test   
if ( Xyce_FORTRAN_TEST )
  enable_language ( Fortran )
endif ( Xyce_FORTRAN_TEST )

# Only if Xyce_ARCHDIR is defined...
if ( DEFINED Xyce_ARCHDIR )
  # inspect ARCHDIR
  if (    NOT EXISTS ${Xyce_ARCHDIR}
       OR NOT EXISTS ${Xyce_ARCHDIR}/include
       OR NOT EXISTS ${Xyce_ARCHDIR}/lib )
    message ( FATAL_ERROR "${Xyce_ARCHDIR} failed inspection." )
  endif ()

  # Add it to CMake's prefix paths.
  list ( APPEND CMAKE_PREFIX_PATH ${Xyce_ARCHDIR} )

  # add common include file directories
  set ( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH}
      ${Xyce_ARCHDIR}/include
      ${Xyce_ARCHDIR}/include/trilinos )


  # add common library directories
  set ( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH}
      ${Xyce_ARCHDIR}/lib
      ${Xyce_BINARY_DIR}/src )
endif ( DEFINED Xyce_ARCHDIR )


# set Dakota support header and library paths
if ( Xyce_Dakota AND DEFINED Xyce_ARCHDIR )

  if ( EXISTS ${Xyce_ARCHDIR}/dakota/include )
    set ( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${Xyce_ARCHDIR}/dakota/include )
  endif ( EXISTS ${Xyce_ARCHDIR}/dakota/include )

  if ( EXISTS ${Xyce_ARCHDIR}/reference/lib ) 
    set ( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${Xyce_ARCHDIR}/reference/lib )
  endif ( EXISTS ${Xyce_ARCHDIR}/reference/lib ) 

  if ( EXISTS ${Xyce_ARCHDIR}/dakota/lib )
    set ( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${Xyce_ARCHDIR}/dakota/lib )
  endif ( EXISTS ${Xyce_ARCHDIR}/dakota/lib )

endif ( Xyce_Dakota AND DEFINED Xyce_ARCHDIR )

# set header and linker flags
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} )
include_directories( ${CMAKE_INCLUDE_PATH} )
link_directories( ${CMAKE_LIBRARY_PATH} )


# get compiler vendor
get_property ( COMPILER_VENDOR VARIABLE PROPERTY CMAKE_CXX_COMPILER_ID ) 


# set Windows (native & cygwin) specific flags
if ( CMAKE_HOST_WIN32 )

  set ( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -nologo" )
  set ( WINDOWS_BUILD ON )

  # disable "safe" features of MSVS
  add_definitions ( -D_SECURE_SCL=0 )


  ###### Eric R does this in his local version with a comment "DEBUG ELR: 
  ######  broken static builds!

  add_definitions ( -MD )

  # replace stl components with their MSVS "safe" counterparts
  ## add_definitions ( -D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES=1 -D_CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES_COUNT=1 )

  set ( Xyce_RESTART_NOPACK ON )

  # remove linker conflicts
  ## set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:libcmtd.lib" )
  #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:libmmd.lib /DEFAULTLIB:libmmt.lib")

endif ( CMAKE_HOST_WIN32 )




# -- system inspection --------------------------------------------------------

# -- search for required programs

if ( Xyce_REACTION_PARSER )
  find_package ( BISON 2.3 )
  find_package ( FLEX )
  if ( NOT FLEX_FOUND OR NOT BISON_FOUND )
    message ( STATUS "Could not find bison/flex.  Disabling REACTION_PARSER" )
    set ( Xyce_REACTION_PARSER OFF )  
  endif ( NOT FLEX_FOUND OR NOT BISON_FOUND )
endif ( Xyce_REACTION_PARSER )


# -- search for libraries

# base math library (no usually available under windows)
find_library( LIB_M m )

if ( Xyce_USE_FFTW )
  XYCE_FIND_LIBRARY ( fftw3 )
endif ( Xyce_USE_FFTW )

if ( Xyce_USE_HDF5 )
  XYCE_FIND_LIBRARY ( hdf5 sz z )
  # need to update this so that sz and z are only used if needed 
  # as we don't use those components
  # and if hdf5 is not compiled with parallel support, then don't
  # activate it for mpi parallel builds
  XYCE_FIND_LIBRARY ( sz )
  XYCE_FIND_LIBRARY ( z )
endif ( Xyce_USE_HDF5 )


# add lapack / blas and associated dependencies
if ( COMPILER_VENDOR STREQUAL "Intel" )

  # use the pure linking mesthod from Intel 10+

  find_library ( LIB_MKL_INTERFACE   
                 NAMES mkl_intel mkl_intel_lp64 mkl_intel_c
                 PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )

  find_library ( LIB_MKL_THREADING         
                 NAMES mkl_sequential
                 PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )

  find_library ( LIB_MKL_COMPUTATION       
                 NAMES mkl_core                        
                 PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )

  if ( Xyce_Dakota )
    find_library ( LIB_FORTRAN               
                   NAMES ifcore                           
                   PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )
  endif ( Xyce_Dakota )

  set ( MKL_LIBRARIES 
        ${LIB_MKL_INTERFACE} ${LIB_MKL_THREADING} ${LIB_MKL_COMPUTATION}
        ${LIB_MKL_INTERFACE} ${LIB_MKL_THREADING} ${LIB_MKL_COMPUTATION}
        ${LIB_FORTRAN} )


#if ( CMAKE_HOST_APPLE AND Xyce_PARALLEL_MPI )
#  find_library ( LIBIRC
#                 NAMES irc
#                 PATHS $ENV{INTEL_PATH}/$ENV{INTEL_COM_VER}/lib )
#endif ( CMAKE_HOST_APPLE AND Xyce_PARALLEL_MPI )


  # Windows (native & cygwin) and Open MPI do not require additional threading support
  if ( NOT Xyce_PARALLEL_MPI AND NOT CMAKE_HOST_WIN32 )
    find_library ( LIB_PTHREAD               
                   NAMES pthread                         
                   PATHS ENV "LIBRARY_PATH" ENV "LD_LIBRARY_PATH" ENV "DYLD_LIBRARY_PATH" )
  endif ( NOT Xyce_PARALLEL_MPI AND NOT CMAKE_HOST_WIN32 )


  # collect all libraries
  set ( EXTERNAL_LIBRARIES 
        ${MKL_LIBRARIES} ${LIB_PTHREAD} ${LIBIRC} )



else ( COMPILER_VENDOR STREQUAL "Intel" )

  find_library ( LIB_LAPACK lapack )
  find_library ( LIB_BLAS blas )

  # This stuff is unreliable and should be replaced by proper probing of
  # the fortran link requirements.  
  # Autoconf will run the Fortran compiler with options to display what 
  # libraries it needs, and add those to the link line.  CMake won't do
  # that automatically unless your project contains Fortran files.
  # This mess tries to add -lgfortran for specific system combinations, but
  # it gets it wrong for some (e.g. FreeBSD with clang).  Needed to add
  # a similar wrong block to get mingw builds to work, too.
  # DEBUG:  fix for peculiar CMake behavior for static linked binary and FBSD
  if ( Xyce_STATIC_BINARY OR CMAKE_SYSTEM_NAME STREQUAL "FreeBSD" )
    set ( LIB_FORTRAN "gfortran" )
  endif ( Xyce_STATIC_BINARY OR CMAKE_SYSTEM_NAME STREQUAL "FreeBSD" )

  if ( CMAKE_SYSTEM_NAME STREQUAL "Windows" AND COMPILER_VENDOR STREQUAL "GNU" )
    set ( LIB_FORTRAN "gfortran" )

  # as long as we've got this mingw-specific conditional, make sure to set
  # this, too:
    set ( Xyce_RESTART_NOPACK ON )

  endif ( CMAKE_SYSTEM_NAME STREQUAL "Windows" AND COMPILER_VENDOR STREQUAL "GNU" )

  # collect all libraries
  set ( EXTERNAL_LIBRARIES 
        ${LIB_LAPACK} ${LIB_BLAS} ${LIB_FORTRAN} )

endif ( COMPILER_VENDOR STREQUAL "Intel" )


if ( Xyce_PARALLEL_MPI )

  # expecting openmpi to be present in /usr/local/* or /usr/local/mpi/*
  message ( STATUS "Looking for library mpi" )                 
  file ( GLOB OPENMPI_SEARCH_PATH "/usr/local/openmpi*" "/usr/local/mpi/openmpi*" )
  find_library ( LIB_MPI 
                 NAMES mpi 
                 PATHS ${OPENMPI_SEARCH_PATH}
                 PATH_SUFFIXES lib lib/openmpi 
                 NO_DEFAULT_PATH )
  find_path ( INC_MPI 
                 NAMES mpi mpi.h
                 PATHS ${OPENMPI_SEARCH_PATH}
                 PATH_SUFFIXES include         
                 NO_DEFAULT_PATH )

  if ( NOT LIB_MPI )
    message ( STATUS "Looking for library mpi - not found" )                 
  else ( NOT LIB_MPI )
    message ( STATUS "Looking for library mpi - found" )                 
    set ( ARCH_LIBRARIES ${LIB_MPI} ${ARCH_LIBRARIES} )
    include_directories ( ${INC_MPI} )
  endif ( NOT LIB_MPI )

  if ( Xyce_USE_ZOLTAN OR Xyce_USE_ISORROPIA )
    if ( DEFINED Xyce_ARCHDIR AND EXISTS ${Xyce_ARCHDIR}/include/zoltan )
      include_directories ( ${Xyce_ARCHDIR}/include/zoltan )
    endif ( DEFINED Xyce_ARCHDIR AND EXISTS ${Xyce_ARCHDIR}/include/zoltan )
    XYCE_FIND_LIBRARY ( metis ) 
    XYCE_FIND_LIBRARY ( parmetis )
    XYCE_FIND_LIBRARY ( zoltan )
  endif ( Xyce_USE_ZOLTAN OR Xyce_USE_ISORROPIA )

  add_definitions ( -DXyce_PARALLEL_MPI )
endif ( Xyce_PARALLEL_MPI )

if ( Xyce_AMD OR Xyce_UMFPACK )
  XYCE_FIND_LIBRARY ( amd )
endif ( Xyce_AMD OR Xyce_UMFPACK )

# We no longer support Trilinos 11.0 in CMake. Upgrade to 11.2.X or later.
# XYCE_FIND_LIBRARY ( teuchos )
# The old single teuchos library of <=11.0 has been replaced by these 5:
XYCE_FIND_LIBRARY ( teuchoscore )
XYCE_FIND_LIBRARY ( teuchosremainder )
XYCE_FIND_LIBRARY ( teuchosnumerics )
XYCE_FIND_LIBRARY ( teuchoscomm )
XYCE_FIND_LIBRARY ( teuchosparameterlist )

XYCE_FIND_LIBRARY ( triutils )
XYCE_FIND_LIBRARY ( epetra )
XYCE_FIND_LIBRARY ( epetraext )

if ( Xyce_PARALLEL_MPI AND Xyce_USE_ISORROPIA )
  XYCE_FIND_LIBRARY ( isorropia ) 
endif ( Xyce_PARALLEL_MPI AND Xyce_USE_ISORROPIA )

if ( Xyce_SUPERLUDIST )
  XYCE_FIND_LIBRARY ( superludist )
endif ( Xyce_SUPERLUDIST )

if ( Xyce_UMFPACK )
  XYCE_FIND_OPTIONAL_LIBRARY (suitesparseconfig)
  XYCE_FIND_OPTIONAL_LIBRARY ( colamd )
  XYCE_FIND_OPTIONAL_LIBRARY ( cholmod )
  XYCE_FIND_LIBRARY ( umfpack )
endif ( Xyce_UMFPACK )

if ( Xyce_SHYLU )
  XYCE_FIND_LIBRARY ( shylu )
endif ( Xyce_SHYLU )

if ( Xyce_SUPERLU )
  XYCE_FIND_LIBRARY ( superlu )
endif ( Xyce_SUPERLU )

XYCE_FIND_LIBRARY ( amesos )
XYCE_FIND_LIBRARY ( trilinoscouplings )

if ( Xyce_ML )
  XYCE_FIND_LIBRARY ( ml )
endif ( Xyce_ML )

XYCE_FIND_LIBRARY ( ifpack )
XYCE_FIND_LIBRARY ( aztecoo )

if ( Xyce_BELOS )
  XYCE_FIND_LIBRARY ( belos )
  XYCE_FIND_LIBRARY ( belosepetra )
endif ( Xyce_BELOS )

XYCE_FIND_LIBRARY ( nox )
XYCE_FIND_LIBRARY ( loca )


# check for the dakota library if the user wants to link to it
# this library list is for Dakota 5.2.  We need a more flexible 
# way to get Dakota's dependent libs automatically.  The Dakota
# team is working on this for Dakota 5.2+
if ( Xyce_Dakota  ) 
  XYCE_FIND_LIBRARY ( dakota_src )
  XYCE_FIND_LIBRARY ( nidr )
  XYCE_FIND_LIBRARY ( pecos ) 
  XYCE_FIND_LIBRARY ( pecos_src )
  XYCE_FIND_LIBRARY ( lhs )
  XYCE_FIND_LIBRARY ( mods )
  XYCE_FIND_LIBRARY ( mod )
  XYCE_FIND_LIBRARY ( dfftpack )
  XYCE_FIND_LIBRARY ( sparsegrid )
  XYCE_FIND_LIBRARY ( surfpack )
  XYCE_FIND_LIBRARY ( surfpack )
  XYCE_FIND_LIBRARY ( utilib )
  XYCE_FIND_LIBRARY ( colin )
  XYCE_FIND_LIBRARY ( interfaces )
  XYCE_FIND_LIBRARY ( scolib )
  XYCE_FIND_LIBRARY ( 3po )
  XYCE_FIND_LIBRARY ( pebbl )
  XYCE_FIND_LIBRARY ( tinyxml )
  XYCE_FIND_LIBRARY ( conmin )
  XYCE_FIND_LIBRARY ( dace )
  XYCE_FIND_LIBRARY ( analyzer )
  XYCE_FIND_LIBRARY ( random )
  XYCE_FIND_LIBRARY ( sampling )
  XYCE_FIND_LIBRARY ( bose )
  XYCE_FIND_LIBRARY ( fsudace )
  XYCE_FIND_LIBRARY ( hopspack )
  XYCE_FIND_LIBRARY ( jega )
  XYCE_FIND_LIBRARY ( jega_fe )
  XYCE_FIND_LIBRARY ( moga )
  XYCE_FIND_LIBRARY ( soga )
  XYCE_FIND_LIBRARY ( eutils )
  XYCE_FIND_LIBRARY ( utilities )
  XYCE_FIND_LIBRARY ( ncsuopt )
  XYCE_FIND_LIBRARY ( nlpql )
  XYCE_FIND_LIBRARY ( cport ) 
  XYCE_FIND_LIBRARY ( npsol )
  XYCE_FIND_LIBRARY ( optpp )
  XYCE_FIND_LIBRARY ( psuade )
  XYCE_FIND_LIBRARY ( amplsolver )
  
  XYCE_FIND_LIBRARY ( boost_regex ) 
  XYCE_FIND_LIBRARY ( boost_signals ) 
  XYCE_FIND_LIBRARY ( boost_filesystem )
  XYCE_FIND_LIBRARY ( boost_system )
              
endif ( Xyce_Dakota  )

# -- search for header files

# system headers
CHECK_INCLUDE_FILE_CXX ( cctype HAVE_CCTYPE )
CHECK_INCLUDE_FILE_CXX ( climits HAVE_CLIMITS )
CHECK_INCLUDE_FILE_CXX ( cmath HAVE_CMATH )
CHECK_INCLUDE_FILE_CXX ( cstdio HAVE_CSTDIO )
CHECK_INCLUDE_FILE_CXX ( cstdlib HAVE_CSTDLIB )
CHECK_INCLUDE_FILE_CXX ( cstring HAVE_CSTRING )
CHECK_INCLUDE_FILE_CXX ( fftw3.h HAVE_FFTW3_H )
CHECK_INCLUDE_FILE_CXX ( float.h HAVE_FLOAT_H )
CHECK_INCLUDE_FILE_CXX ( fstream HAVE_FSTREAM )
CHECK_INCLUDE_FILE_CXX ( hdf5.h HAVE_HDF5_H )
CHECK_INCLUDE_FILE_CXX ( iostream HAVE_IOSTREAM )
CHECK_INCLUDE_FILE_CXX ( limits.h HAVE_LIMITS_H )
CHECK_INCLUDE_FILE_CXX ( math.h HAVE_MATH_H )
CHECK_INCLUDE_FILE_CXX ( ostream HAVE_OSTREAM )
CHECK_INCLUDE_FILE_CXX ( ostream.h HAVE_OSTREAM_H )
CHECK_INCLUDE_FILE_CXX ( unistd.h HAVE_UNISTD_H )
CHECK_INCLUDE_FILE_CXX ( values.h HAVE_VALUES_H )
CHECK_INCLUDE_FILE_CXX ( strings.h HAVE_STRINGS_H )
CHECK_INCLUDE_FILE_CXX ( string.h HAVE_STRING_H )
CHECK_INCLUDE_FILE_CXX ( unordered_map HAVE_UNORDERED_MAP)
CHECK_INCLUDE_FILE_CXX ( tr1/unordered_map HAVE_TR1_UNORDERED_MAP)
CHECK_INCLUDE_FILE_CXX ( functional HAVE_FUNCTIONAL)
CHECK_INCLUDE_FILE_CXX ( tr1/functional HAVE_TR1_FUNCTIONAL)


# Trilinos headers
CHECK_INCLUDE_FILE_CXX ( Epetra_MultiVector.h HAVE_EPETRA_MULTIVECTOR_H )
#Don't do :
#CHECK_INCLUDE_FILE_CXX ( Ifpack_CrsRiluk.h HAVE_IFPACK_CRSRILUK_H )
# as it tries to link a test program, which fails due to dependence on teuchos
# libraries.  check_include_file_cxx can't be forced to add libraries to its
# link line!
find_path ( HAVE_IFPACK_CRSRILUK_H Ifpack_CrsRiluk.h )
find_path ( HAVE_AMESOS_KLU_H Amesos_Klu.h )
find_path ( HAVE_LOCA_PARAMETER_VECTOR_H LOCA_Parameter_Vector.H  )
find_path ( HAVE_NOX_ABSTRACT_VECTOR_H NOX_Abstract_Vector.H )

#if ( Xyce_PARALLEL_MPI AND Xyce_USE_ZOLTAN AND NOT Xyce_USE_HMETIS )
  find_path ( HAVE_ZOLTAN_QUERYOBJECT_H Zoltan_QueryObject.h )
#endif ( Xyce_PARALLEL_MPI AND Xyce_USE_ZOLTAN AND NOT Xyce_USE_HMETIS )

#if ( Xyce_UMFPACK )
   find_path(HAVE_AMESOS_UMFPACK_H Amesos_Umfpack.h)
#endif ( Xyce_UMFPACK )

#if ( Xyce_SUPERLU )
   find_path(HAVE_AMESOS_SUPERLU_H Amesos_Superlu.h)
#endif ( Xyce_SUPERLU )

#if ( Xyce_Dakota )
  find_path ( HAVE_DAKOTA_STRATEGY_H DakotaStrategy.H ) 
#endif ( Xyce_Dakota )


# -- search for symbols

CHECK_FUNCTION_EXISTS ( strcasecmp HAVE_STRCASECMP )

#    set ( CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} )
if( HAVE_FLOAT_H )
  CHECK_FUNCTION_EXISTS( _isnan FUNCTION__ISNAN )
  CHECK_FUNCTION_EXISTS( _finite FUNCTION__FINITE )
  if( FUNCTION__ISNAN AND FUNCTION__FINITE ) 
    set( HAVE__ISNAN_AND__FINITE_SUPPORT on )
  endif( FUNCTION__ISNAN AND FUNCTION__FINITE ) 
endif( HAVE_FLOAT_H )


if (NOT (FUNCTION__ISNAN AND FUNCTION__FINITE))
  if ( HAVE_CMATH OR HAVE_MATH_H )
    set ( CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LIB_M} )
    CHECK_FUNCTION_EXISTS ( isnan FUNCTION_ISNAN )
    CHECK_FUNCTION_EXISTS ( isinf FUNCTION_ISINF )
    if ( FUNCTION_ISNAN AND FUNCTION_ISINF )
      set ( HAVE_NAN_INF_SUPPORT ON )
    endif ( FUNCTION_ISNAN AND FUNCTION_ISINF )
  else ( HAVE_CMATH OR HAVE_MATH_H )
  endif ( HAVE_CMATH OR HAVE_MATH_H )
endif(NOT (FUNCTION__ISNAN AND FUNCTION__FINITE))


# message ( STATUS "Warning: Your compiler doesn't support isnan() and isinf().  We will supply a default checker but it is *NOT* guaranteed to work on your platform!" )

if ( HAVE_MATH_H )
  CHECK_SYMBOL_EXISTS ( M_PI math.h MATH_H_HAS_M_PI )
endif( HAVE_MATH_H )

if ( HAVE_FFTW3_H )
  CHECK_SYMBOL_EXISTS ( fftw_plan_dft_1d fftw3.h Xyce_USE_FFTW )
endif( HAVE_FFTW3_H )

if ( HAVE_HDF5_H )
  CHECK_SYMBOL_EXISTS ( H5Fclose hdf5.h Xyce_USE_HDF5 )
endif( HAVE_HDF5_H )

# search for Intel  MKL FFT header and library compatibility
CHECK_INCLUDE_FILE_CXX ( mkl_dfti.h HAVE_MKL_DFTI_H )
if ( HAVE_MKL_DFTI_H )
  if ( CMAKE_HOST_WIN32 AND COMPILER_VENDOR STREQUAL "Intel" )
    set ( Xyce_USE_INTEL_FFT ON )
  else ( CMAKE_HOST_WIN32 AND COMPILER_VENDOR STREQUAL "Intel" )
    CHECK_SYMBOL_EXISTS ( Dfti_CreateDescriptor mkl_dfti.h Xyce_USE_INTEL_FFT )
  endif ( CMAKE_HOST_WIN32 AND COMPILER_VENDOR STREQUAL "Intel" )
endif ( HAVE_MKL_DFTI_H )

if (Xyce_USE_INTEL_FFT OR Xyce_USE_FFTW)
  set (Xyce_USE_FFT ON)
else (Xyce_USE_INTEL_FFT OR Xyce_USE_FFTW)
  set (Xyce_USE_FFT OFF)
endif (Xyce_USE_INTEL_FFT OR Xyce_USE_FFTW)



# -- search for required functions

if ( Xyce_FORTRAN_TEST )
  CHECK_FORTRAN_FUNCTION_EXISTS ( dcopy HAVE_FORTRAN_DCOPY )  
  CHECK_FORTRAN_FUNCTION_EXISTS ( dgetrs HAVE_FORTRAN_DGETRS )  
else ( Xyce_FORTRAN_TEST )
  CHECK_FUNCTION_EXISTS ( dcopy HAVE_DCOPY )  
  CHECK_FUNCTION_EXISTS ( dgetrs HAVE_DGETRS )  
endif ( Xyce_FORTRAN_TEST )

CHECK_FUNCTION_EXISTS ( sqrt HAVE_SQRT )
CHECK_FUNCTION_EXISTS ( getopts HAVE_GETOPTS )

try_compile ( HAVE_IOTA           
              ${CMAKE_BINARY_DIR}
              ${CMAKE_SOURCE_DIR}/config/cmake_check_iota.cxx )


# -- post-processing ----------------------------------------------------------


# use the config.h files for Trilinos
set ( HAVE_CONFIG_H ON )

# flag use of CMake
set ( USING_CMAKE ON )   

# add flag to config files
set ( Xyce_NOX_LOCA_SUPPORT ON )


# set lexer flags
if ( Xyce_REACTION_PARSER ) 
  if ( NOT HAVE_UNISTD_H )
    set ( YY_NO_UNISTD_H ON ) 
  endif ( NOT HAVE_UNISTD_H )
endif ( Xyce_REACTION_PARSER )


# DEBUG:  holdover to old version string approach
# generate the build timestamp and make available from timestamp.h
execute_process ( COMMAND date "+%Y%m%d%H%M" 
                  OUTPUT_VARIABLE XYCEBUILDTIMESTAMP 
                  OUTPUT_STRIP_TRAILING_WHITESPACE )
set ( XYCEBUILDTIMESTAMP "${XYCEBUILDTIMESTAMP}" CACHE INTERNAL "MUTABLE" FORCE )
add_custom_target ( TIMESTAMP_H ALL
                    COMMAND touch ${Xyce_SOURCE_DIR}/src/UtilityPKG/timestamp.h.cmake )
configure_file ( ${Xyce_SOURCE_DIR}/src/UtilityPKG/timestamp.h.cmake
                 ${Xyce_BINARY_DIR}/src/UtilityPKG/timestamp.h )
include_directories ( ${Xyce_BINARY_DIR}/src/UtilityPKG )


# generate the version strings
set ( Xyce_VERSION_STRING_SHORT "${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}" )

if ( DEFINED Xyce_VERSION_PATCH AND NOT (Xyce_VERSION_PATCH STREQUAL "" OR Xyce_VERSION_PATCH STREQUAL "0"))
  set ( Xyce_VERSION_STRING_SHORT "${Xyce_VERSION_STRING_SHORT}.${Xyce_VERSION_PATCH}" )
endif ( DEFINED Xyce_VERSION_PATCH AND NOT (Xyce_VERSION_PATCH STREQUAL "" OR Xyce_VERSION_PATCH STREQUAL "0"))

set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_SHORT}" )

if ( DEFINED Xyce_VERSION_EXTRA )
  set ( Xyce_VERSION_STRING_LONG "${Xyce_VERSION_STRING_LONG} ${Xyce_VERSION_EXTRA}" )
endif ( DEFINED Xyce_VERSION_EXTRA )


if ( Xyce_ATHENA )
  if ( Xyce_RAD_MODELS )
    if( EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA/ )
      find_package(Boost 1.41 REQUIRED COMPONENTS)
      include_directories(${Boost_INCLUDE_DIRS})
      set ( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
      XYCE_FIND_LIBRARY ( boost_filesystem )
      XYCE_FIND_LIBRARY ( boost_system )
    else ( EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA/ )
      message ( STATUS "Xyce_RAD_MODELS is enabled, but the SandiaModels/ATHENA directory was not found." )
      message ( STATUS "Try updating your version of SandiaModels." )
      message ( STATUS "Disabling Xyce_ATHENA." )
      set ( Xyce_ATHENA OFF )
    endif ( EXISTS ${Xyce_SOURCE_DIR}/src/DeviceModelPKG/SandiaModels/ATHENA/ )
  else ( Xyce_RAD_MODELS )
    message ( STATUS "Xyce_RAD_MODELS not enabled.  Disabling Xyce_ATHENA." )
    set ( Xyce_ATHENA OFF )
  endif ( Xyce_RAD_MODELS )
endif ( Xyce_ATHENA )


# -- generate source files and update the include list -----------------------

configure_file ( "${Xyce_SOURCE_DIR}/src/Xyce_config.h.cmake"
                 "${Xyce_BINARY_DIR}/src/Xyce_config.h" )

set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES}
      ${Xyce_BINARY_DIR}/src 
      CACHE INTERNAL "X_I" )




# -- process subdirectories --------------------------------------------------

add_subdirectory ( src )




# -- enable dashboard and CTest ----------------------------------------------

# DEBUG: should this be in Xyce_Regression as an include() instead?
set ( DROP_METHOD "http" )
set ( DROP_SITE "s858869.sandia.gov" )
set ( DROP_LOCATION "/CDash/submit.php?project=RELEASE_6_1" )
include ( CTest )

