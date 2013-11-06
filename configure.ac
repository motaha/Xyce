dnl Process this file with autoconf to produce a configure script.
dnl NEVER edit the configure script itself, as it will be overwritten
dnl when autoconf processes this file.

dnl ********************************************************************
dnl Autoconf boilerplate
dnl Initialize the macros, define the package name, version number and 
dnl bug reporting address
AC_INIT(Xyce,R::6.0,erkeite@sandia.gov)

dnl This is a FILTHY HACK designed to get rid of the troublesome space that
dnl AC_INIT stuffs into PACKAGE_STRING, and which buggers MPILAM's "hcp" 
dnl script.

sed -e 's/PACKAGE_STRING\(.*\) \(.*\)/PACKAGE_STRING\1_\2/g' < confdefs.h > foobie.h && mv foobie.h confdefs.h
PACKAGE_STRING=`echo $PACKAGE_STRING|sed -e 's/ /_/g'`

dnl Identify a single file within the source tree to make sure we're 
dnl configuring the right source
AC_CONFIG_SRCDIR(src/CircuitPKG/src/N_CIR_Xyce.C)

dnl Auxillary Build Tools
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

# Configure should create Xyce_config.h from Xyce_config.h.in
AC_CONFIG_HEADERS(src/Xyce_config.h)

# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

dnl Require version 2.59 or higher of autoconf
AC_PREREQ(2.59)

dnl Make sure the CVS Revision tag that appears in this file gets copied
dnl to the configure script generated from it, but in a way that CVS won't
dnl corrupt in the configure script --- that way you can always tell which
dnl revision of configure.ac was used to generate a given configure

AC_REVISION($Revision: 1.355.2.11 $)

dnl Determine a canonicalized system type (e.g. i386-unknown-freebsd4.4)
AC_CANONICAL_TARGET

dnl Initialize Automake
AM_INIT_AUTOMAKE([1.10.1 foreign])

dnl ********************************************************************
dnl the autoconfing proper 

dnl *****************************************************************
dnl Before we can get to the meat, we need to do some options
dnl *****************************************************************


dnl This is braindamage --- the SGI doesn't use AR to make libraries, it
dnl uses the CXX compiler and not AR.
AM_CONDITIONAL(SGI, test "$Xyce_TARGET" = "sgi10k" -o "$Xyce_TARGET" = "sgin32")

dnl Forte compilers have special limitations.  Also see src/Makefile.am's
AM_CONDITIONAL(SUN, test "$Xyce_TARGET" = "solaris")


dnl ***********************************************************************
dnl The Xyce project has a directory full of pre-built libraries for
dnl each architecture.  Use it if we can.
dnl If the environment variable ARCH_ROOT is set, use it, 
dnl otherwise pick a safe default
dnl Allow user to specify an "ARCH_ROOT", and kludge a -L for it
if test "x$ARCH_ROOT" = x; then
 if test -d /Net/Proj/Xyce/arch; then
   AC_MSG_NOTICE([Using /Net/Proj/Xyce/arch as default ARCH_ROOT])
   ARCH_ROOT=/Net/Proj/Xyce/arch
 fi
fi

dnl But let the user override the default
AC_ARG_WITH(archroot,
[AS_HELP_STRING([--with-archroot=path], [The root of architecture-dependent local libraries and includes (DEPRECATED)])],
[ARCH_ROOT=$withval
 AC_MSG_NOTICE([Using $ARCH_ROOT to find platform libraries and includes by user request.  --with-archroot is a deprecated configure option.  Please use "ARCH_ROOT=" instead.])
]
)

dnl and let the user completely override the ARCH_ROOT/Xyce_ARCH
dnl hierarchy
AC_ARG_WITH(archdir,
[AS_HELP_STRING([--with-archdir=path], [Where to find include and libs for all external libraries instead of ARCH_ROOT/Xyce_ARCH/ (DEPRECATED)])],
[ARCHDIR=$withval
 AC_MSG_NOTICE([Using $ARCHDIR/include and $ARCHDIR/lib for locating external libraries and includes by user request.  --with-archdir is a deprecated configure option.  Please use "ARCHDIR=" instead.])
])

dnl But we can't use ARCH_ROOT to set libraries yet, because we don't
dnl know if we're doing a serial or parallel build yet.
dnl so let's process more options
dnl ********************************************
dnl These are deprecated, unnecessary options from Tammy Kolda that allowed overriding defaults, but there
dnl are already built-in methods to do this through autoconf.  Deprecated!
AC_ARG_WITH(cc,
[AS_HELP_STRING([--with-cc], [Gives the name of the C compiler to use. (DEPRECATED)])],
[CC=$withval
AC_MSG_ERROR([--with-cc is a deprecated configure option.  Please use "CC=" instead.])
]
)

AC_ARG_WITH(cxx,
[AS_HELP_STRING([--with-cxx], [Gives the name of the C++ compiler to use. (DEPRECATED)])],
[CXX=$withval
AC_MSG_ERROR([--with-cxx is a deprecated configure option.  Please use "CXX=" instead.])
]
)

AC_ARG_WITH(f77,
[AS_HELP_STRING([--with-f77], [Gives the name of the fortran compiler to use. (DEPRECATED)])],
[F77=$withval
AC_MSG_ERROR([--with-f77 is a deprecated configure option.  Please use "F77=" instead.])
]
)

AC_ARG_WITH(include,
[AS_HELP_STRING([--with-include], [Specifies include directories ("-Idir1 -Idir2 ...") (DEPRECATED)])],
[
XTRAINC=$withval
AC_MSG_ERROR([Using extra include directories: $withval.   --with-include is a deprecated configure option.  Please use "CPPFLAGS=" instead.])
]
)

AC_ARG_WITH(libs,
[AS_HELP_STRING([--with-libs], [Specifies libraries ("-Ldir -llib ...") (DEPRECATED)])],
[
XTRALIBS=$withval
AC_MSG_ERROR([Using extra libraries: $withval.  --with-libs is a deprecated configure option.  Please use "LIBS=" instead.])
]
)

AC_ARG_WITH(flags,
[AS_HELP_STRING([--with-flags], [Specifies flags ("-Dflag1 -Dflag2 ...") (DEPRECATED)])],
[
XTRAFLAGS=$withval
AC_MSG_ERROR([Using extra flags: $withval.  --with-flags is a deprecated configure option.  Please use "CPPFLAGS" instead.])
]
)

dnl ********************************************


dnl ********************************************
dnl Check to see if the executable is disabled.
dnl RPP: 03/15/2006 - For Xyce/Charon coupling,
dnl we need to turn off building of the executable,
dnl autotools is screwing the link line in Nevada's TPL
AC_ARG_ENABLE(xyce-binary,
[AS_HELP_STRING([--disable-xyce-binary],[Disables building of the xyce executable. (Default is to build the binary)])],
[ENABLE_XYCE_BINARY=$enableval],
[ENABLE_XYCE_BINARY=yes]
)

AC_MSG_CHECKING(whether to build the xyce binary)
if test "X$ENABLE_XYCE_BINARY" != "Xno"; then
  AC_MSG_RESULT(yes)
else
  CPPFLAGS="-DDINT $CPPFLAGS"
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL(BUILD_XYCE_BINARY, test "X$ENABLE_XYCE_BINARY" != "Xno")
dnl ********************************************


dnl ********************************************
dnl Here we determine whether we're doing a parallel build or not
XYCE_DEBUG_OPTION(mpi,no,[parallel build with MPI],Xyce_PARALLEL_MPI,USE_MPI)

if test "$USE_MPI" = yes; then
   if test "$Xyce_TARGET" = "sgi10k" -o "$Xyce_TARGET" = "sgin32"; then
     CXXFLAGS="$CXXFLAGS -DMPI_NO_CPPBIND"
   fi
   Xyce_ARCH="$Xyce_ARCH"_MPI
   Xyce_TARGET="$Xyce_TARGET"_MPI
fi

#Lastly, we allow the user to completely override our choice of Xyce_ARCH
# if need be.

AC_ARG_WITH(xyce_arch,
[AS_HELP_STRING([--with-xyce_arch], [specifies which subdirectory of ARCH_ROOT to use for default locations of many packages (DEPRECATED)])],
[
  Xyce_ARCH=$withval
  AC_MSG_NOTICE([Using Xyce_ARCH: $withval.  --with-xyce_arch is a deprecated configure option.  Please use "Xyce_ARCH=" instead.])
]
)

dnl And now we know precisely what build we're doing, so maybe we can use
dnl ARCH_ROOT
dnl but if user has specified ARCHDIR use it instead of $ARCH_ROOT/$Xyce_ARCH
if test "x$ARCHDIR" = "x"; then
 if test "x$ARCH_ROOT" != "x"; then
  if test "x$Xyce_ARCH" != "x"; then
    if test -d $ARCH_ROOT/$Xyce_ARCH; then
      ARCHDIR=$ARCH_ROOT/$Xyce_ARCH
    fi
  fi  
 fi
fi


dnl see comments at header of XYCE_DEBUG_OPTION in acinclude.m4 for
dnl usage
XYCE_DEBUG_OPTION(verbose_time,no,[verbosity in time integrator])
XYCE_DEBUG_OPTION(verbose_nonlinear,no,[verbosity in nonlinear solver])
XYCE_DEBUG_OPTION(verbose_nox,no,[verbosity in NOX nonlinear solver library])
XYCE_DEBUG_OPTION(verbose_linear,no,[verbosity in linear solver])
XYCE_DEBUG_OPTION(spice_norms,yes,[using SPICE type norms])
XYCE_DEBUG_OPTION(adaptive_solve,no,[adaptive iterative linear solvers])
XYCE_DEBUG_OPTION(bsim3_const,no,[using constants from BSIM3],Xyce_USE_BSIM3_CONST)
XYCE_DEBUG_OPTION(nox_use_vector_copy,yes,[using vector copies in NOX])
XYCE_DEBUG_OPTION(sensitivity,no,[sensitivity capability],Xyce_SENSITIVITY_ENABLE, SENSITIVITY_ENABLE)

dnl optionally compile in support for external device support
XYCE_DEBUG_OPTION(extdev,yes,[extdev extensions],Xyce_EXTDEV,USE_EXTDEV)
AM_CONDITIONAL(EXTDEV,test $USE_EXTDEV = yes)

dnl optionally compile in support for rad models if the directory is there and
dnl the user hasn't specified disable-radmodels
MAYBE_RAD=""
if test -d $srcdir/src/DeviceModelPKG/SandiaModels; then
  XYCE_DEBUG_OPTION(radmodels,yes,[Radiation Model Library],Xyce_RAD_MODELS,RAD_MODELS)

  if test "x$RAD_MODELS" = "xyes"
  then
    MAYBE_RAD="SandiaModels"
  fi
else
  RAD_MODELS=no
fi
AM_CONDITIONAL(RADMODELS,test $RAD_MODELS = yes)
AC_SUBST([MAYBE_RAD])

dnl optionally compile in support for non-free models if the directory is there and
dnl the user hasn't specified disable-nonfree
MAYBE_NONFREE=""    
if test -d $srcdir/src/DeviceModelPKG/Xyce_NonFree; then
  XYCE_DEBUG_OPTION(nonfree,yes,[Non-Free Models Library],Xyce_NONFREE_MODELS,NONFREE_MODELS)
  if test "x$NONFREE_MODELS" = "xyes" 
  then
    MAYBE_NONFREE="Xyce_NonFree"
  fi
else
  NONFREE_MODELS=no
fi
AM_CONDITIONAL(NONFREEMODELS,test $NONFREE_MODELS = yes)
AC_SUBST([MAYBE_NONFREE])

dnl optionally compile in support for Dakota linkage
XYCE_DEBUG_OPTION(dakota,no,[Dakota direct linkage support],Xyce_Dakota, USE_DAKOTA)
AM_CONDITIONAL(DAKOTA_LINK, test $USE_DAKOTA = yes)
dnl For the LHS package in Dakota we need to link in a Dakota object file as well as
dnl the dakota libraries.  
AM_CONDITIONAL(DAKOTA_OBJ_NEEDED, test $USE_DAKOTA = yes )

dnl optionally compile in support for Dakota linkage
XYCE_DEBUG_OPTION(dakota50,no,[Dakota 5.0 direct linkage support],Xyce_Dakota50, USE_DAKOTA50)
AM_CONDITIONAL(DAKOTA_LINK, test $USE_DAKOTA50 = yes -o $USE_DAKOTA = yes )
dnl For the LHS package in Dakota we need to link in a Dakota object file as well as
dnl the dakota libraries. (duplicating this test until I rip out the old dakota support) 
AM_CONDITIONAL(DAKOTA_OBJ_NEEDED, test $USE_DAKOTA50 = yes )

dnl optionally compile in support for Charon linkage
XYCE_DEBUG_OPTION(charon,no,[Charon device support],Xyce_CHARON, USE_CHARON)
AM_CONDITIONAL(BUILD_CHARON_SUPPORT, test $USE_CHARON = yes)

dnl optionally compile in support for Trilinos Development
XYCE_DEBUG_OPTION(trilinos_dev,no,[Trilinos development support],Xyce_TRILINOS_DEV, USE_TRILINOS_DEV)
AM_CONDITIONAL(TRILINOS_DEV, test $USE_TRILINOS_DEV = yes)

dnl optionally compile in support for Trilinos Development with Rythmos
XYCE_DEBUG_OPTION(trilinos_dev_rythmos,no,[Trilinos development support with Rythmos],Xyce_TRILINOS_DEV_RYTHMOS, USE_TRILINOS_DEV_RYTHMOS)
AM_CONDITIONAL(TRILINOS_DEV_RYTHMOS, test $USE_TRILINOS_DEV_RYTHMOS = yes)

dnl compile in support for partitioning
XYCE_DEBUG_OPTION(isorropia,yes,[isorropia],Xyce_USE_ISORROPIA,USE_ISORROPIA)
XYCE_DEBUG_OPTION(zoltan,yes,[zoltan],Xyce_USE_ZOLTAN,USE_ZOLTAN)

# This does nothing anymore --- there are no conditionals in the code that
# need it, nor are there any conditionals in the Makefile.am's
# 
# Philosophical note: it is a very bad thing to port to a platform using
# platform-specific macros.  If something surfaces that requires special 
# treatment, it should be handled by "feature-specific" ifdefs that can be
# probed by configure, that way they'll automatically be applied to any 
# platform on which they're needed.  
# 
#dnl optionally compile in support for Red Storm / PGI compilers 
#XYCE_DEBUG_OPTION(redstorm,no,[Red Storm / PGI support],Xyce_RedStorm, USE_REDSTORM)
#AM_CONDITIONAL(REDSTORM, test $USE_REDSTORM = yes)

dnl optionally compile in support for SuperLU and UmfPack.
XYCE_DEBUG_OPTION(superlu,no,[SuperLU direct solver],Xyce_SUPERLU,SUPERLU)
XYCE_DEBUG_OPTION(umfpack,yes,[UmfPack direct solver],Xyce_UMFPACK,UMFPACK)
XYCE_DEBUG_OPTION(pardiso_mkl,no,[Pardiso direct solver (through MKL)],Xyce_PARDISO_MKL,PARDISO_MKL)

dnl optionally compile in ksparse solver
XYCE_DEBUG_OPTION(ksparse,no,[Kundert's Sparse direct solver],Xyce_KSPARSE,USE_KSPARSE)
AM_CONDITIONAL(BUILD_KSPARSE, test "x$USE_KSPARSE" = "xyes")

dnl optionally compile in support for Belos (as of Xyce 6.0, this is "yes" by default)
XYCE_DEBUG_OPTION(belos,yes,[Belos iterative solver support],Xyce_BELOS, USE_BELOS)
AM_CONDITIONAL(BUILD_BELOS_SUPPORT, test "x$USE_BELOS" = "xyes")

dnl optionally compile in support for ShyLU (we need Trilinos development branch at this time)
XYCE_DEBUG_OPTION(shylu,no,[ShyLU hybrid solver support],Xyce_SHYLU, USE_SHYLU)
AM_CONDITIONAL(BUILD_SHYLU_SUPPORT, test "x$USE_SHYLU" = "xyes")

dnl optionally compile in support for ML (we need the new preconditioning factory)
XYCE_DEBUG_OPTION(ml,no,[ML preconditioner support],Xyce_ML, USE_ML)
AM_CONDITIONAL(BUILD_ML_SUPPORT, test "x$USE_ML" = "xyes")

dnl Must include AMD
XYCE_DEBUG_OPTION(amd,yes,[AMD reordering],Xyce_AMD,AMD)

if test "$USE_MPI" = "yes"; then
  XYCE_DEBUG_OPTION(superludist,no,[SuperLU_Dist linear solver],Xyce_SUPERLUDIST,SUPERLUDIST)
else
  SUPERLUDIST=no
fi

XYCE_DEBUG_OPTION(nox_unstable,no,[NOX/LOCA experimental code - usually requires latest nox and loca development library],Xyce_NOX_UNSTABLE, NOX_UNSTABLE)

XYCE_DEBUG_OPTION(test_soln_var_map,no,[mapping of internal variables to external for debugging])
XYCE_DEBUG_OPTION(debug_device,no,[device package debugging output])
XYCE_DEBUG_OPTION(debug_io,no,[I/O package debugging])
XYCE_DEBUG_OPTION(debug_expression,no,[Expression package debugging output])
XYCE_DEBUG_OPTION(debug_time,no,[time integrator package debugging output])
XYCE_DEBUG_OPTION(debug_analysis,no,[analysis package debugging output])
XYCE_DEBUG_OPTION(debug_parallel,no,[parallel distribution package debugging output])
XYCE_DEBUG_OPTION(debug_distribution,no,[distribution package debugging output])
XYCE_DEBUG_OPTION(debug_topology,no,[topology package debugging output])
XYCE_DEBUG_OPTION(debug_linear,no,[linear solver package debugging output])
XYCE_DEBUG_OPTION(debug_nonlinear,no,[nonlinear solver package debugging output])
XYCE_DEBUG_OPTION(debug_circuit,no,[circuit package debugging output])
XYCE_DEBUG_OPTION(debug_directsolve,no,[direct solver package debugging output])
XYCE_DEBUG_OPTION(debug_restart,no,[restart debugging output])

XYCE_DEBUG_OPTION(op_start,yes,[allow startup from previous operating point])

dnl check if the developer wants to skip the probe of the fortran compiler
XYCE_DEBUG_OPTION(fortran_test,yes,[test for fortran compiler],Xyce_FORTRAN_TEST,FORTRAN_TEST)

# Allow disabling of reaction parser
XYCE_DEBUG_OPTION(reaction_parser,yes,[reaction parser], Xyce_REACTION_PARSER,REACTION_PARSER)

# Allow disabling of FFT libraries
XYCE_DEBUG_OPTION(fft,yes,[support for FFT libraries],Xyce_FFT, USE_FFT)
if test "$USE_FFT" = "yes"; then
  XYCE_DEBUG_OPTION(intel_fft,yes,[support for Intel MKL FFT library],Xyce_IntelFFT, USE_IntelFFT)
  XYCE_DEBUG_OPTION(fftw,yes,[support for FFTW library],Xyce_FFTW, USE_FFTW)
else
  USE_FFTW=no
  USE_IntelFFT=no
fi

XYCE_DEBUG_OPTION([athena], [no], [building the ATHENA device])

# Default to HDF5 support off 
# Note:  a fourth argument *must* be given, because XYCE_DEBUG_OPTION 
# forces a "-D" into CPPFLAGS.  But because this particular option may be
# disabled later when we check if HDF5 libraries and headers are there, the
# fourth argument must not be the same one that Xyce uses to test for HDF5.
# That is, the fourth argument is just "user has asked for this feature"
# not "this feature works"
XYCE_DEBUG_OPTION(hdf5,no,[support for HDF5 output library],Xyce_HDF5, USE_HDF5)

# static exec and lib built by default; use --enable-shared to get shared lib
AC_DISABLE_SHARED

dnl Checks for programs.
# Checks for programs.
AC_PROG_CXX(CC c++ g++ gcc cxx icpc)
AC_PROG_CC(cc gcc icc)
if test "$NO_FORTRAN_TEST" != "yes"; then
  AC_PROG_F77(f77 g77 g95 gfortran ifort)
fi


if test "$REACTION_PARSER" = "yes"; then
 # check for flex --- we use flex-only features, and can't use vanilla lex
 # if don't find flex, don't build parser
 # Same goes for bison, as yacc won't do the job for us, and AC_PROG_YACC
 # will use bison in yacc compatibility mode if we use it.  WE also specifically
 # need bison 2.3, not the bison 1.8x that comes stock with some linuxen
 # Further, bison's LALR1.cc skeleton changed in a way that is incompatible
 # with how we used it with version 2.1 --- Maintaining backward compatibility
 # isn't worth the trouble
 AC_PROG_LEX
 XYCE_CHECK_BISON([ 2.7, 2.6, 2.5, 2.4, 2.3], , [BISON=no])
 if test `basename $LEX` = "flex";
 then
    AC_DEFINE(HAVE_FLEX,1,[define if you have the flex lexical scanner generator])
 else
    AC_MSG_ERROR(["Requested reaction parser, but no FLEX."])
 fi
 if test "$BISON" != "no";
 then
   AC_DEFINE(HAVE_BISON,1,[define if you have the bison parser generator])
 else
   AC_MSG_ERROR(["Requested reaction parser, but no bison of correct version number."])
 fi
else
 LEX=no
 BISON=no
fi

AM_CONDITIONAL(BUILD_REACTION_PARSER, test `basename $LEX` = "flex" -a "$BISON" != "no")

dnl Check to see if valid Fortran compiler was found (HKT; 07/11/06; bug #1091).
if test "$FORTRAN_TEST" != "no" -a "x$F77" = "x" ; then
  AC_MSG_ERROR([Fortran testing requested but no valid Fortran compiler was found in your path.])
fi

AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

dnl *********************************************************************
dnl KLUDGES:
dnl The following kludges add defined constants that are highly machine
dnl and compiler dependent --- so let's do them after we figure out what 
dnl machine, build and compiler we're using

dnl needed to move this thing to AFTER AC_PROG_CXX because putting it before
dnl that test disables the automatic insertion of -g -O2 in the default build.

if test "x$Xyce_TARGET" = "xsgi10k" -o "x$Xyce_TARGET" = "xsgin32" -o "x$Xyce_TARGET" = "xsgi10k_MPI" -o "x$Xyce_TARGET" = "xsgin32_MPI"; then
  LIBS="-lfpe $LIBS"
fi

USING_DEC_CXX=no
XYCE_TEST_USING_DEC_CXX
# add __USE_STD_IOSTREAM on the DEC compilers
if test "$USING_DEC_CXX" != "no"; then
  CPPFLAGS="-D__USE_STD_IOSTREAM $CPPFLAGS"
  # force inclusion of template symbols in cxx_repository; see cxx manpage
  CXXFLAGS="-std gnu $CXXFLAGS"
fi

USING_ICCLT11_CXX=no
XYCE_TEST_USING_ICCLT11_CXX
ADMS_CXXFLAGS="$CXXFLAGS"
# special CXX flags for ADMS directory on the ICC version <=10 compilers
if test "$USING_ICCLT11_CXX" != "no"; then
  ADMS_CXXFLAGS=`echo $CXXFLAGS | sed -e 's/-O3/-O1/' -e 's/-O2/-O1/'`
fi

dnl This is also braindamage --- the DEC leaves miscellaneous droppings
dnl in "cxx_repository" directories and needs to be cleaned out
AM_CONDITIONAL(DEC, test $USING_DEC_CXX = yes)


dnl Include the top-level of the include and lib directories in the search
dnl Add to loader flags so comes before all -l's
if test -d $ARCHDIR/lib; then
  LDFLAGS="$LDFLAGS -L$ARCHDIR/lib"
fi
if test -d $ARCHDIR/include; then
  CPPFLAGS="-I$ARCHDIR/include $CPPFLAGS"
fi

dnl Include the trilinos include directory in search
if test -d $ARCHDIR/include/trilinos; then
  CPPFLAGS="-I$ARCHDIR/include/trilinos $CPPFLAGS"
fi

dnl Include Dakota headers and library paths if asked for Dakota build
if test "x$USE_DAKOTA50" = "xyes"; then
  CPPFLAGS="-DXyce_Dakota $CPPFLAGS"
  if test -d $ARCHDIR/dakota5_0; then
    DAKOTA_ARCH_DIR="$ARCHDIR/dakota5_0"
  fi
  if test -d $ARCHDIR/dakota5_0/include; then
    default_INCS="$default_INCS -I$ARCHDIR/dakota5_0/include"
    CPPFLAGS="-I$ARCHDIR/dakota5_0/include $CPPFLAGS"
  fi
  if test -d $ARCHDIR/dakota5_0/lib; then
    LDFLAGS="$LDFLAGS -L$ARCHDIR/dakota5_0/lib"
  fi
elif test "x$USE_DAKOTA" = "xyes"; then
  if test -d $ARCHDIR/dakota5_0; then
    DAKOTA_ARCH_DIR="$ARCHDIR/dakota"
  fi
  if test -d $ARCHDIR/dakota/include; then
    default_INCS="$default_INCS -I$ARCHDIR/dakota/include"
    CPPFLAGS="-I$ARCHDIR/dakota/include $CPPFLAGS"
  fi
  if test -d $ARCHDIR/reference/lib; then
    LDFLAGS="$LDFLAGS -L$ARCHDIR/reference/lib"
  fi
  if test -d $ARCHDIR/dakota/lib; then
    LDFLAGS="$LDFLAGS -L$ARCHDIR/dakota/lib"
  fi
fi


dnl Add some FLIBS and LDFLAGS for g95 based on Eric Keiter's Intel-Mac (HKT: 1/24/2007; bug# 1098)
if test "x$F77" = "xg95"; then
  FLIBS="-L/sw/lib -lg95 $FLIBS"
  LDFLAGS="-Wl,-framework,Accelerate $LDFLAGS"
fi

dnl *********************************************************************
dnl checks for libraries now
dnl Replace `main' with a function in -lm:
if test "$FORTRAN_TEST" != "no"; then
  AC_F77_LIBRARY_LDFLAGS
fi

AC_SEARCH_LIBS(sqrt, [m])

dnl find a valid blas and lapack library

if test "$FORTRAN_TEST" != "no"; then
  AC_LANG_PUSH(Fortran 77)
  AC_SEARCH_LIBS(dcopy,[essl cxml blas_slu blas])
  AC_SEARCH_LIBS(dgetrs,[essl cxml scs complib.sgimath lapack cxml])
  AC_LANG_POP(Fortran 77)
else
  AC_SEARCH_LIBS(dcopy,[essl cxml blas_slu blas])
  AC_SEARCH_LIBS(dgetrs,[essl cxml scs complib.sgimath lapack cxml])
fi


#It should now be safe to add FLIBS to LIBS
LIBS="${LIBS} ${FLIBS}"

dnl CAREFUL --- if parallel is enabled, epetra needs MPI, so better have
dnl it later in the libs path! (which means detecting it earlier, since
dnl LIBS is assembled right to left)

dnl and if this is a parallel build we need all the Zoltan stuff
if test "x$USE_MPI" = "xyes"; then
  AC_SEARCH_LIBS(MPI_Bcast,[mpi])
  if test "x$USE_ZOLTAN" = "xyes" -o "x$USE_ISORROPIA" = "xyes"; then
    XYCE_LIBRARY_INCLUDE(zoltan)
    AC_CHECK_LIB(metis,METIS_PartGraphKway,[],[AC_MSG_ERROR([FATAL: METIS library not found or not valid])])
    AC_CHECK_LIB(parmetis,ParMETIS_PartKway,[],
       [AC_CHECK_LIB(parmetis,ParMETIS_V3_PartKway,[],
            [AC_MSG_ERROR([FATAL: ParMETIS library not found or not valid])])])
    AC_CHECK_LIB(zoltan,Zoltan_LB_Balance,[],[AC_MSG_ERROR([FATAL: zoltan library not found or not valid])])
  fi
fi

dnl If the user enabled AMD, let's try to add it here:
if test "x$AMD" = "xyes" -o "x$UMFPACK" = "xyes"; then
  AC_CHECK_LIB(amd,amd_valid,[],[AC_MSG_ERROR([FATAL: AMD library does not contain amd_valid and this is required by Trilinos.])])
fi

# Look for teuchos.  First, look for the old one, where teuchos was in a 
# single library.  If not found, look for the new one.
# NOTE: If we find the old one, we're going to accept it, but warn the user
# loudly at the end of the configure run (look at the bottom of this file
# for where the warning happens.
XYCE_TEUCHOS_NOT_FOUND=0
XYCE_LIBRARY_OPTION(teuchos,[],[],[],[XYCE_TEUCHOS_NOT_FOUND=1])

# Now look for the new style if we haven't found the old.  Any one of these
# not being found is a fatal configure error.
if test $XYCE_TEUCHOS_NOT_FOUND = 1; then
  XYCE_LIBRARY_OPTION(teuchoscore)
  XYCE_LIBRARY_OPTION(teuchosremainder)
  XYCE_LIBRARY_OPTION(teuchosnumerics)
  XYCE_LIBRARY_OPTION(teuchoscomm)
  XYCE_LIBRARY_OPTION(teuchosparameterlist)
fi
  
XYCE_LIBRARY_OPTION(triutils)
XYCE_LIBRARY_OPTION(epetra)
XYCE_LIBRARY_OPTION(epetraext)

dnl and if this is a parallel build we (may) need all the Isorropia stuff
if test "x$USE_MPI" = "xyes" -a "x$USE_ISORROPIA" = "xyes"; then
  XYCE_LIBRARY_OPTION(isorropia)
fi

dnl If the user enabled superludist, let's try to add it here:
if test "x$SUPERLUDIST" = "xyes"; then
    XYCE_LIBRARY_OPTION(superludist)
fi

dnl If the user enabled umfpack, let's try to add it here:
if test "x$UMFPACK" = "xyes"; then
  # cholmod may depend on colamd...
  AC_CHECK_LIB(colamd,colamd,[],[])
  # Newer cholmods may depend on suitesparseconfig.  If we find it, add it.
  AC_CHECK_LIB(suitesparseconfig,SuiteSparse_time,[],[])
  # Newer umfpacks depend on cholmod.  If we find it, add it
  AC_CHECK_LIB(cholmod,cholmod_start,[],[])
  AC_CHECK_LIB(umfpack,umfpack_di_solve,[],[AC_MSG_ERROR([FATAL: UmfPack library does not contain umfpack_di_solve and this is required by Trilinos.])])
fi

dnl If the user enabled superlu, let's try to add it here:
if test "x$SUPERLU" = "xyes"; then
  AC_CHECK_LIB(superlu,dgssvx,[],[AC_MSG_ERROR([FATAL: SuperLU library does not contain dgssvx and this is required by Trilinos.])])
fi

XYCE_LIBRARY_OPTION(amesos)
XYCE_LIBRARY_OPTION(trilinoscouplings)

dnl If the user enabled ML, let's add it here:
if test "x$USE_ML" = "xyes"; then 
  XYCE_LIBRARY_OPTION(ml)
fi
XYCE_LIBRARY_OPTION(ifpack)
XYCE_LIBRARY_OPTION(aztecoo)

dnl If the user enabled Belos, let's add it here:
if test "x$USE_BELOS" = "xyes"; then
    XYCE_LIBRARY_OPTION(belos)
    XYCE_LIBRARY_OPTION(belosepetra)
fi

dnl If the user enabled ShyLU, let's add it here:
if test "x$USE_SHYLU" = "xyes"; then
    XYCE_LIBRARY_OPTION(shylu)
fi

XYCE_LIBRARY_OPTION(nox)
XYCE_LIBRARY_OPTION(loca)
CPPFLAGS="-DXyce_NOX_LOCA_SUPPORT $CPPFLAGS"

dnl check for the dakota library if the user wants to link to it
if test "x$USE_DAKOTA50" = "xyes"; then
  XYCE_LIBRARY_OPTION(amplsolver)
  XYCE_LIBRARY_OPTION(cdd)
  XYCE_LIBRARY_OPTION(shared)
  XYCE_LIBRARY_OPTION(conveyor) 
  XYCE_LIBRARY_OPTION(appspack)
  XYCE_LIBRARY_OPTION(nappspack)
  XYCE_LIBRARY_OPTION(3po)
  XYCE_LIBRARY_OPTION(utilib)
  XYCE_LIBRARY_OPTION(pebbl)
  XYCE_LIBRARY_OPTION(colin)
  XYCE_LIBRARY_OPTION(coliny)
  XYCE_LIBRARY_OPTION(quadrature)  
  XYCE_LIBRARY_OPTION(ncsuopt)
  XYCE_LIBRARY_OPTION(newmat)
  XYCE_LIBRARY_OPTION(psuade)
  XYCE_LIBRARY_OPTION(opt) 
  XYCE_LIBRARY_OPTION(npsol)
  XYCE_LIBRARY_OPTION(nlpql)
  XYCE_LIBRARY_OPTION(cport)
  XYCE_LIBRARY_OPTION(jega)
  XYCE_LIBRARY_OPTION(fsudace)
  XYCE_LIBRARY_OPTION(ddace)
  XYCE_LIBRARY_OPTION(conmin)
  XYCE_LIBRARY_OPTION(surfpack)
  XYCE_LIBRARY_OPTION(lhs)
  XYCE_LIBRARY_OPTION(pecos)
  XYCE_LIBRARY_OPTION(dakota)
elif test "x$USE_DAKOTA" = "xyes"; then
  XYCE_LIBRARY_OPTION(fl)
  XYCE_LIBRARY_OPTION(amplsolver)
  XYCE_LIBRARY_OPTION(3po)
  XYCE_LIBRARY_OPTION(appspack)
  XYCE_LIBRARY_OPTION(utilib)
  XYCE_LIBRARY_OPTION(pebbl)
  XYCE_LIBRARY_OPTION(colin)
  XYCE_LIBRARY_OPTION(coliny)
  XYCE_LIBRARY_OPTION(cubature)
  XYCE_LIBRARY_OPTION(gsl)
  XYCE_LIBRARY_OPTION(lhs)
  XYCE_LIBRARY_OPTION(ncsuopt)
  XYCE_LIBRARY_OPTION(newmat)
  XYCE_LIBRARY_OPTION(opt)
  XYCE_LIBRARY_OPTION(npsol)
  XYCE_LIBRARY_OPTION(cport)
  XYCE_LIBRARY_OPTION(jega)
  XYCE_LIBRARY_OPTION(fsudace)
  XYCE_LIBRARY_OPTION(ddace)
  XYCE_LIBRARY_OPTION(conmin)
  XYCE_LIBRARY_OPTION(surfpack)
  XYCE_LIBRARY_OPTION(evidence)
  XYCE_LIBRARY_OPTION(dakota)

  dnl libmpi.a is only needed if we're using the prebuilt binary of Dakota
  dnl not if we've built it from source
  dnl XYCE_LIBRARY_OPTION(mpi)
  dnl XYCE_LIBRARY_OPTION(lam)
  dnl

fi

dnl *****************************************************************
dnl Now include the include directories in CPPFLAGS
dnl We don't need to check directory existence, because we cleared the
dnl variables above if the directory didn't exist.

CPPFLAGS="${CPPFLAGS} ${INCDIRS}"

dnl Lastly, add the current directory to CPPFLAGS
CPPFLAGS="-I. $CPPFLAGS"	

dnl *********************************************************************
dnl Add all the extra things we offered the user the option of including
dnl --------- Set CPPFLAGS and bindir ---------

CPPFLAGS="${CPPFLAGS} ${XTRAFLAGS} ${XTRAINC}"

dnl --------- Set and Check Libs ---------
dnl add the ones we requested on the command line
LIBS="${LIBS} ${XTRALIBS}"

echo $LIBS

AC_MSG_CHECKING(whether libraries (${LIBS}) will link)
AC_LANG_PUSH(C++)
AC_LINK_IFELSE(
             [AC_LANG_PROGRAM([[]],[[]])],
[AC_MSG_RESULT([yes])],
[
AC_MSG_RESULT([no])
AC_MSG_ERROR([FATAL: cannot link with ${LIBS}... check paths, fix with --with-package options.   See config.log for failure message.])
]
)

#dnl put original LDFLAGS back and keep -Ldirs on left
#if test "$Xyce_ARCH" = "solaris"; then
#    LDFLAGS="$TMP_LDFLAGS $LDFLAGS"
#fi

dnl *********************************************************************
dnl Checks for header files.
AC_HEADER_STDC
dnl This is a showstopper, and must be done first!
AC_CXX_HAVE_STL

AC_LANG_PUSH(C++)
AC_CHECK_HEADERS(unistd.h values.h limits.h float.h math.h cmath climits ostream iostream ostream.h fstream cstdio algorithm algo.h cctype cstdlib cstring strings.h string.h unordered_map functional tr1/unordered_map tr1/functional)

dnl *********************************************************************
dnl See if compiler supplies a nan/inf checker.  This must be done after the 
dnl cmath probe, because some platforms have wonky cmath/math.h pairs.
AC_FINITE_NUMBER_CHECK

dnl Need to define HAVE_CONFIG_H because AC_CHECK_HEADER doesn't include the 
dnl CXXFLAGS on the command line when it tries to compile the test program.
dnl Oddly, AC_CHECK_HEADERS does.
AC_CHECK_HEADER(Epetra_MultiVector.h,,[AC_MSG_ERROR([Sorry, Petra headers not found])],)
AC_CHECK_HEADER(Ifpack_CrsRiluk.h,,[AC_MSG_ERROR([Sorry, Ifpack headers not found])],)

if test "x$USE_MPI" = "xyes"; then
if test "x$USE_ZOLTAN" = "xyes"; then
if test "x$USE_HMETIS" = "xno"; then
  AC_CHECK_HEADER(Zoltan_QueryObject.h,,[AC_MSG_ERROR([Sorry, Zoltan headers not found])])
fi 
fi 
fi 
AC_CHECK_HEADER(Amesos_Klu.h,,[AC_MSG_ERROR([Sorry, KLU headers not found])])
if test "x$UMFPACK" = "xyes"; then
  AC_CHECK_HEADER(Amesos_Umfpack.h,,[AC_MSG_ERROR([Sorry, Umfpack headers not found])])
fi
if test "x$SUPERLU" = "xyes"; then
  AC_CHECK_HEADER(Amesos_Superlu.h,,[AC_MSG_ERROR([Sorry, SuperLU headers not found])])
fi

AC_CHECK_HEADER(NOX_Abstract_Vector.H,,[AC_MSG_ERROR([Sorry, NOX headers not found])])
AC_CHECK_HEADER(LOCA_Parameter_Vector.H,,[AC_MSG_ERROR([Sorry, LOCA headers not found])])
AC_CHECK_HEADER([Sacado.hpp], [], [AC_MSG_ERROR([Sorry, Sacado headers not found])])

dnl
dnl check for the dakota headers if the user wants to link to Dakota 
dnl Dakota on the ppc needs "DARWIN" defined
dnl
if test "x$USE_DAKOTA" = "xyes"; then
  if test "x$Xyce_TARGET" = "xppc"; then
    CPPFLAGS="$CPPFLAGS -DDARWIN "
    LDFLAGS="$LDFLAGS  -Wl,-m " 
  fi
  AC_CHECK_HEADER(DakotaStrategy.H,,[AC_MSG_ERROR([Sorry, Dakota headers not found])])
fi

#
# check if the intel math library header is around
# this will define ac_cv_header_mkl_dfti_h = yes or blank depending on the result
# we'll check if we can actually link to a function in the intel math library later
# and then set Xyce_USE_INTEL_FFT as needed.
AC_CHECK_HEADER([mkl_dfti.h])

AC_LANG_POP(C++)

dnl Checks for typedefs, structures, and compiler characteristics.

dnl this one is prereq for other stuff:
dnl AC_CXX_NAMESPACES 
AC_CXX_HAVE_NUMERIC_LIMITS
AC_HEADER_MATH_H_HAS_M_PI
AC_CXX_NEED_PAIR_H

#check if we can use floating point exceptions on linux
##AC_CXX_HAVE_FENV

dnl Checks for library functions.
AC_CHECK_FUNCS(strcasecmp)

dnl Checks for library functions.
AC_CHECK_FUNCS(getopt)

# check if FFT library support is enabled. If so, check for Intel FFT or FFTW.
if test "X$USE_FFT" = "Xyes"; then

AC_DEFINE([Xyce_USE_FFT],[1], [Set to use either FFTW or Intel's FFT library.])
#
# check if the Intel FFT functions are available 
# this should only pass if 
#  1) the right headers are around and (which was checked above in AC_CHECK_HEADER([mkl_dfti.h]))
#  2) the right math libraries have been specified
if test "X$USE_IntelFFT" = "Xyes"; then
  if test "X$ac_cv_header_mkl_dfti_h" = "Xyes"; then
    AC_CHECK_FUNCS( [DftiCreateDescriptor], [AC_DEFINE([Xyce_USE_INTEL_FFT],[1], [Set to use Intel's FFT library])],[AC_MSG_WARN([Could not find Intel FFT Library.  Disabling Intel FFT support])])
  fi
else
  ac_cv_header_mkl_dtfi_h=no
fi

# check if the FFTW functions are available 
#
# check if the FFTW library header is around (version 3+), and if it is,
# check for the library
if test "X$USE_FFTW" = "Xyes"; then
  AC_CHECK_HEADER([fftw3.h], [AC_SEARCH_LIBS( [fftw_plan_dft_1d], [fftw3], [AC_DEFINE([Xyce_USE_FFTW],[1], [Set to use FFTW library])],[AC_MSG_WARN([Could not find FFTW Library.  Disabling FFTW support])])])
else
 ac_cv_header_fftw3_h=no
fi  

# ERROR, neither FFT library was found and FFT was enabled.
if test "X$ac_cv_header_mkl_dfti_h" != "Xyes" -a "X$ac_cv_header_fftw3_h" != "Xyes"; then
  AC_MSG_ERROR(
[Configure has not found an FFT library (FFTW or Intel FFT) in your path.
                  An FFT library is required for frequency-domain analysis capabilities.
                  Reconfigure Xyce with either FFTW or Intel FFT or disable FFT support
                  with --disable-fft.])
fi

else
  ac_cv_header_mkl_dtfi_h=no
  ac_cv_header_fftw3_h=no
fi
AM_CONDITIONAL(BUILD_INTEL_FFT_SUPPORT, test "X$ac_cv_header_mkl_dfti_h" = "Xyes")
AM_CONDITIONAL(BUILD_FFTW_SUPPORT, test "X$ac_cv_header_fftw3_h" = "Xyes")

# check if the HDF5 functions are available 
#
# check if the HDF5 library header is around (version 1.8+), and if it is,
# check for the library (Note. HDF5 can depend on libsz.a and libz.a,
# which are available on their web site.  The libz.a used by hdf5 IS NOT the 
# same as libzc.a used for compressing files.)

if test "X$USE_HDF5" = "Xyes"; then
  AC_SEARCH_LIBS( [inflate], [z] )
  AC_SEARCH_LIBS( [SZ_Decompress], [sz])
  AC_CHECK_HEADER( [hdf5.h],
                   [AC_SEARCH_LIBS( [H5Fclose], [hdf5],
                                    [AC_DEFINE([Xyce_USE_HDF5],[1], [Set to use HDF5 library])],
                                    [AC_MSG_WARN([Could not find HDF5 Library.  Disabling HDF5 support])
                                     ac_cv_header_hdf5_h=no
                                    ])]
                 )
else
 ac_cv_header_hdf5_h=no
fi  
AM_CONDITIONAL(BUILD_HDF5_SUPPORT, test "X$ac_cv_header_hdf5_h" = "Xyes")

#check if the -ansi breaks random_shuffle, as it does on RedHat 7.1.  This
# may be a generic thing for GCCs and GNU libc combinations.
if test "x$GXX" = "xyes"; then
  #this is only a problem if using GCC
  XYCE_CHECK_BRAINDAMAGED_RANDOM_SHUFFLE
  
  if test "$ac_cv_xyce_braindamaged_random_shuffle" = yes; then
    CPPFLAGS="-D__STL_NO_DRAND48 $CPPFLAGS"
  fi
fi

XYCE_CHECK_IOTA
XYCE_CHECK_FLEXIBLE_INSERT

if test "$Xyce_TARGET" = "solaris"; then
  XYCE_CHECK_STL_MEMBER_TEMPLATES
fi

dnl LASTLY:  add -nofor_main to FFLAGS for DEC.  Can't do it earlier than this
if test "$USING_DEC_CXX" != "no"; then
  FFLAGS="$FFLAGS -nofor_main"
fi

dnl necessary Xyce substitution so all the makefiles get it.
dnl This makes @Xyce_INCS@ be replaced with the full string of -I's needed
dnl for all packages to compile without problem.
dnl It is used by setting "INCLUDES = @Xyce_INCS@" in the Makefile.am for
dnl a package.
default_INCS="$default_INCS -I\$(top_srcdir)/src/AnalysisPKG/include -I\$(top_srcdir)/src/CircuitPKG/include -I\$(top_srcdir)/src/IOInterfacePKG/include -I\$(top_srcdir)/src/IOInterfacePKG/Output/include -I\$(top_srcdir)/src/DeviceModelPKG/include -I\$(top_srcdir)/src/DeviceModelPKG/PDESC/include -I\$(top_srcdir)/src/DeviceModelPKG/EXTSC/include -I\$(top_srcdir)/src/DeviceModelPKG/AD/include -I\$(top_srcdir)/src/DeviceModelPKG/ADMS/include -I\$(top_srcdir)/src/ErrorHandlingPKG/include -I\$(top_srcdir)/src/LinearAlgebraServicesPKG/include -I\$(top_srcdir)/src/LoaderServicesPKG/include -I\$(top_srcdir)/src/NonlinearSolverPKG/include -I\$(top_srcdir)/src/ParallelDistPKG/include -I\$(top_srcdir)/src/TimeIntegrationPKG/include -I\$(top_srcdir)/src/TopoManagerPKG/include -I\$(top_srcdir)/src/UtilityPKG/include -I\$(top_srcdir)/src/MultiTimePDEPKG/include -I\$(top_srcdir)/src/DakotaLinkPKG/include -I\$(top_builddir)/src"
if test "x$RAD_MODELS" = "xyes"; then
  default_INCS="$default_INCS -I\$(top_srcdir)/src/DeviceModelPKG/SandiaModels/include"
fi
if test "x$NONFREE_MODELS" = "xyes"; then
  default_INCS="$default_INCS -I\$(top_srcdir)/src/DeviceModelPKG/Xyce_NonFree/include"
fi

if test "x$USE_KSPARSE" = "xyes"; then
  default_INCS="$default_INCS -I\$(top_srcdir)/src/LinearAlgebraServicesPKG/ksparse/include"
fi
AC_SUBST(Xyce_INCS,[$default_INCS])

dnl Set the ATHENA requirements.
AS_IF([test x$enable_athena = xyes],
      [AS_IF([test x$RAD_MODELS = xyes],
             [AS_IF([test -d "${srcdir}/src/DeviceModelPKG/SandiaModels/ATHENA"],
                    [XYCE_LIBRARY_OPTION( boost_filesystem )
                     XYCE_LIBRARY_OPTION( boost_system )
                     AC_DEFINE([Xyce_ATHENA], [1], [Build the ATHENA device.])],
                    [AC_MSG_ERROR([Your SandiaModels directory is too old.  ATHENA is not present.])])],
             [AC_MSG_WARN([ATHENA requested, but RAD_MODELS were not.  Disabling ATHENA.])
              enable_athena=no])
      ])
AM_CONDITIONAL([BUILD_ATHENA], [test x$enable_athena = xyes])

# This is so we can use xyce_verify before doing a make install
XYCE_BUILDDIR=`pwd`
xyce_savedir=`pwd`
cd $srcdir
XYCE_SRCDIR=`pwd`
cd $xyce_savedir
AC_SUBST(XYCE_BUILDDIR)
AC_SUBST(XYCE_SRCDIR)

AC_SUBST(TS_LIBS)

XYCELIBS="${LIBS}"
LIBS=" "
AC_SUBST(XYCELIBS)
AC_SUBST(LIBTOOLFLAGS)

# shouldn't need to export this, but Dakota requires this for one of its object files
AC_SUBST(DAKOTA_ARCH_DIR)

# Ugly stuff:  We must use different CXXFLAGS for Intel compilers of the wrong
# version in the ADMS directory:
AC_SUBST(ADMS_CXXFLAGS,[$ADMS_CXXFLAGS])

#############################################################################
# CPACK STUFF
# FINALLY:  Now that we know everything about this build, let's construct
# a set of variables that we can use to generate a CPack config file, so
# that we can use CPack to make a native installer.
# Sorry, but "make dist" just ain't good enough for us, and we do want
# the ability to create native installers.

# Make sure we at least have a generator defined.  TGZ is just a tarball.  We
# will certainly want to override this on the configure command line instead.
# our "reconfigure.cfg" file has settings for our supported platforms, so 
# this variable is always set on the command line in the build/test cycle.
if test "x$CPACK_GENERATOR" = "x" 
then
   CPACK_GENERATOR="TGZ"
fi
AC_SUBST([CPACK_GENERATOR])

# Now, let's decode the VERSION string that was defined in AC_INIT above.
DEVEL_OR_RELEASE=`echo $VERSION | awk -F: '{print $1}'`
if test $DEVEL_OR_RELEASE = "D" 
then
   DEVEL_OR_RELEASE="DEVELOPMENT"
else
   DEVEL_OR_RELEASE="Release"
fi

DEVEL_PHASE=`echo $VERSION | awk -F: '{print $2}'`
if test "x$DEVEL_PHASE" != "x" 
then
   DEVEL_PHASE="-$DEVEL_PHASE"
fi

VERSION_NUMBER=`echo $VERSION | awk -F: '{print $3}'`
Xyce_VERSION_MAJOR=`echo $VERSION_NUMBER | awk -F. '{print $1}'`
Xyce_VERSION_MINOR=`echo $VERSION_NUMBER | awk -F. '{print $2}'`
if test "x$Xyce_VERSION_MINOR" = "x"
then
  Xyce_VERSION_MINOR=0
fi
Xyce_VERSION_PATCH=`echo $VERSION_NUMBER | awk -F. '{print $3}'`
if test "x$Xyce_VERSION_PATCH" = "x"
then
  Xyce_VERSION_PATCH=0
fi

Xyce_VERSION_LONG="${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}"
if test "x$Xyce_VERSION_EXTRA" != "x"
then
  Xyce_VERSION_LONG="${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}-${Xyce_VERSION_EXTRA}"
fi

Xyce_VERSION="${Xyce_VERSION_MAJOR}.${Xyce_VERSION_MINOR}.${Xyce_VERSION_PATCH}"

AC_SUBST([Xyce_VERSION])
AC_SUBST([Xyce_VERSION_MAJOR])
AC_SUBST([Xyce_VERSION_MINOR])
AC_SUBST([Xyce_VERSION_PATCH])

CPACK_PACKAGE_NAME="Xyce-$DEVEL_OR_RELEASE-$Xyce_VERSION_LONG$DEVEL_PHASE"
AC_SUBST([CPACK_PACKAGE_NAME])
if test "x$CPACK_PACKAGE_FILE_NAME" = "x"
then
  CPACK_PACKAGE_FILE_NAME="Install_Xyce_$Xyce_VERSION_MAJOR.$Xyce_VERSION_MINOR.$Xyce_VERSION_PATCH"
fi
AC_SUBST([CPACK_PACKAGE_FILE_NAME])

CPACK_RPM_PACKAGE_NAME=$CPACK_PACKAGE_NAME
AC_SUBST([CPACK_RPM_PACKAGE_NAME])

# Set CPACK_SYSTEM_NAME and default packaging prefix
# This is a little annoying, but if you package on OSX with anything other than
# / for the install prefix
# you get something like /usr/local/Xyce.../usr/local/bin, which is just
# nuts.
CPACK_PACKAGING_INSTALL_PREFIX="/usr/local"
case "$target" in
  *-freebsd*) 
    CPACK_SYSTEM_NAME="FreeBSD"
    ;;
  *-linux*)
    CPACK_SYSTEM_NAME="Linux"
    ;;
  *-apple-darwin*)
    CPACK_SYSTEM_NAME="Darwin"
    CPACK_PACKAGING_INSTALL_PREFIX="/"

    ;;
  *)
    CPACK_SYSTEM_NAME=$target
    ;;
esac
AC_SUBST([CPACK_SYSTEM_NAME])
AC_SUBST([CPACK_PACKAGING_INSTALL_PREFIX])

# OK, with this one we're going to try to force generated RPMs for OPENMPI
# builds to pull in the openmpi executable package.  Will it work?
if test "x$USE_MPI" = "xyes"  -a "$CPACK_GENERATOR" = "RPM"
then
  CPACK_RPM_DEPENDS_SET_LINE="SET (CPACK_RPM_PACKAGE_REQUIRES \"openmpi >= 1.4\")"
else
  CPACK_RPM_DEPENDS_SET_LINE=""
fi
AC_SUBST([CPACK_RPM_DEPENDS_SET_LINE])

# Finally, figure out what license we're using.

XYCE_LICENSE=ECI
CPACK_RPM_PACKAGE_LICENSE="Export Controlled"
if test $RAD_MODELS = "no"
then
  if test $NONFREE_MODELS = "no"
  then
     XYCE_LICENSE=OS
     CPACK_RPM_PACKAGE_LICENSE="GPLv3"
  else
     XYCE_LICENSE=NonFree
     CPACK_RPM_PACKAGE_LICENSE="NDA Required"
  fi
fi
AC_SUBST([XYCE_LICENSE])
AC_SUBST([CPACK_RPM_PACKAGE_LICENSE])

#
# END CPACK STUFF
############################################################################

AC_CONFIG_FILES([Makefile 
 src/Makefile 
 src/AnalysisPKG/Makefile 
 src/CircuitPKG/Makefile
 src/IOInterfacePKG/Makefile
 src/IOInterfacePKG/Output/Makefile
 src/DakotaLinkPKG/Makefile 
 src/DeviceModelPKG/Makefile 
 src/DeviceModelPKG/PDESC/Makefile 
 src/DeviceModelPKG/AD/Makefile 
 src/DeviceModelPKG/ADMS/Makefile 
 src/DeviceModelPKG/EXTSC/Makefile 
 src/ErrorHandlingPKG/Makefile 
 src/LinearAlgebraServicesPKG/Makefile 
 src/LoaderServicesPKG/Makefile 
 src/NonlinearSolverPKG/Makefile 
 src/ParallelDistPKG/Makefile 
 src/TimeIntegrationPKG/Makefile 
 src/TopoManagerPKG/Makefile 
 src/UtilityPKG/Makefile 
 src/MultiTimePDEPKG/Makefile 
 src/test/Makefile
 src/test/XyceAsLibrary/Makefile
 src/test/FFTInterface/Makefile
 src/test/LinearAlgebraTest/Makefile
 src/test/XygraTestHarnesses/Makefile
 src/IOInterfacePKG/include/N_IO_XMLPath.h
 distribution/CPackConfig.cmake])

#NOTE:
# By putting the makefiles into a shell variable like this we are in fact
# deliberately ignoring the guidance in Automake documentation for a specific
# reason.
# 
# When automake encounters the AC_CONFIG_FILES line here, it cannot tell
# that it needs to create the Makefile.in that AC_CONFIG_FILES will need.
# Thus, it is necessary to tell automake on its command line that it needs
# to do that generation.
#
# Why?  Because these two directories are optional, and are not generally 
# distributed with Xyce --- they are not open source directories, and we can't
# ship them.  But at automake/bootstrap time, automake can't tell that these
# AC_CONFIG_FILES are conditionalized, and if the directories don't exist
# then it pukes with an error that it can't find required files.
#
# So the solution is to confuse automake into thinking it has nothing to do
# automatically, ever, and explicitly tell it what to do in bootstrap if
# the directories actually exist.

if test "x$RAD_MODELS" = "xyes"; then
  rm_makefile=src/DeviceModelPKG/SandiaModels/Makefile
  AC_CONFIG_FILES([$rm_makefile])
fi
if test "x$NONFREE_MODELS" = "xyes"; then
  nf_makefile=src/DeviceModelPKG/Xyce_NonFree/Makefile
  AC_CONFIG_FILES([$nf_makefile])
fi

AC_OUTPUT

# Noisy warning stuff that we shouldn't miss
if test $XYCE_TEUCHOS_NOT_FOUND != 1; then
  AC_MSG_WARN(
[===========================================================])

  AC_MSG_WARN(
[Configure detected a version of Trilinos prior to 11.2.3.  
                   Important bugs were fixed in Trilinos 11.2 that impact 
                   Xyce.  You should update your Trilinos])
  AC_MSG_WARN(
[===========================================================])
fi
