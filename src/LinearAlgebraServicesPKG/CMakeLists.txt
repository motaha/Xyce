
# -- process files for this library ------------------------------------------

# add option class sources
if ( Xyce_BELOS )
  set ( BELOS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_BelosSolver.C )
endif ( Xyce_BELOS )

if ( Xyce_ML )
  set ( ML_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_MLPrecond.C )
endif ( Xyce_ML )

if ( Xyce_SHYLU )
  set ( SHYLU_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_ShyLUSolver.C )
endif ( Xyce_SHYLU )


## DEBUG:  ksparse is in subdir, but is not a lib
if ( Xyce_KSPARSE )
  set ( KSPARSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_KSparseSolver.C 
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_KSparseSolver.C
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/Epetra_CrsKundertSparse.C
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/alloc.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/driver.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/smplinkrows.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/spalloc.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/spbuild.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/spfactor.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/spoutput.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/spsmp.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/spsolve.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/src/sputils.c )
  set ( KSPARSE_INC ${CMAKE_CURRENT_SOURCE_DIR}/ksparse/include )
endif ( Xyce_KSPARSE )


# class source list
set ( LinearAlgebraServices_SOURCES   
      ${BELOS_SRC}
      ${ML_SRC}
      ${SHYLU_SRC}
      ${KSPARSE_SRC}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_SimpleSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_AmesosSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_AztecOOSolver.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_BlockMatrix.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_BlockVector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_BlockSystemHelpers.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_Builder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_HBBlockJacobiPrecond.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_HBBlockJacobiEpetraOperator.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_HBFDJacobianPrecond.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_HBFDJacobianEpetraOperator.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_HBBuilder.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_HBPrecondFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_IfpackPrecond.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_LAFactory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_Matrix.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_MultiVector.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_MOROperators.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_Problem.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_SolverFactory.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_System.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_TransformTool.C 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_TrilinosPrecondFactory.C
      ${CMAKE_CURRENT_SOURCE_DIR}/src/N_LAS_Vector.C
      CACHE INTERNAL "X_S" )


# add header files to project
set ( XYCE_INCLUDE_DIRECTORIES ${XYCE_INCLUDE_DIRECTORIES} 
      ${CMAKE_CURRENT_SOURCE_DIR}/include
      ${KSPARSE_INC}
      CACHE INTERNAL "X_I" )


# -- build requested libraries -----------------------------------------------
set ( XYCE_PACKAGE_LIBRARIES ${XYCE_PACKAGE_LIBRARIES} 
      LinearAlgebraServices 
      CACHE INTERNAL "X_L" )
