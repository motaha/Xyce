//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: HarmonicBalance.C,v $
// Purpose       : This file contains unit tests for Harmonic Balance
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/10/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.4 $
// Revision Date  : $Date: 2011/02/18 21:49:15 $
// Current Owner  : $Author: hkthorn $
//-----------------------------------------------------------------------------

// ---------- Trilinos Includes ----------
#include <Epetra_MultiVector.h>
// ---------- Standard Includes ----------
#include <Teuchos_UnitTestHarness.hpp>
// ---------- HB Includes ----------
#include <HB_Builder_Helpers.h>
#include <HB_Loader_Helpers.h>
#include <N_HB_Builder.h>
#include <N_HB_Loader.h>
#include <N_IO_CmdParse.h>
#include <N_DEV_DeviceInterface.h>
 
// Test Harmonic Balance Builder functions to make sure they're creating the right size vectors.
TEUCHOS_UNIT_TEST( N_HB_Builder, createVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = rcp(hbBuilderRCPtr->createVector(0.0));
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  // This is an expanded real form transpose block vector, so numBlocks = 3, blockLength = 20
  TEST_EQUALITY_CONST( bvec->blockCount(), 3 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 20 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createStateVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = rcp(hbBuilderRCPtr->createStateVector(0.0));
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  // This is a time-domain block vector numBlocks = 10, blockLength = 5
  TEST_EQUALITY_CONST( bvec->blockCount(), 10 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 5 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createTimeDomainBlockVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createTimeDomainBlockVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 10 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 3 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createCompressedRealFormVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createCompressedRealFormVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 1 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 11 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createExpandedRealFormVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createExpandedRealFormVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 1 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 20 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createCompressedRealFormBlockVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createCompressedRealFormBlockVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 11 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 3 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createExpandedRealFormBlockVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createExpandedRealFormBlockVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 20 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 3 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createCompressedRealFormTransposeBlockVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createCompressedRealFormTransposeBlockVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 3 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 11 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, createExpandedRealFormTransposeBlockVector ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<N_LAS_Vector> vec = hbBuilderRCPtr->createExpandedRealFormTransposeBlockVector();
  TEST_EQUALITY_CONST( Teuchos::is_null(vec), false );
  RefCountPtr<N_LAS_BlockVector> bvec = Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(vec,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(bvec), false );
  TEST_EQUALITY_CONST( bvec->blockCount(), 3 ); 
  TEST_EQUALITY_CONST( bvec->blockSize(), 20 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, getSolutionMap ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<const Epetra_Map> solMap = hbBuilderRCPtr->getSolutionMap();
  // This is an expanded real form transpose map for a block vector, so numBlocks = 3, blockLength = 20
  TEST_EQUALITY_CONST( solMap->NumGlobalElements(), 60 );
}

TEUCHOS_UNIT_TEST( N_HB_Builder, getStateMap ) {
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(10,3,5);
  RefCountPtr<const Epetra_Map> stateMap = hbBuilderRCPtr->getStateMap();
  // This is a time domain map for a block vector, so numBlocks = 10, blockLength = 5
  TEST_EQUALITY_CONST( stateMap->NumGlobalElements(), 50 );
}

TEUCHOS_UNIT_TEST( N_HB_Loader, create ) {
  RefCountPtr<N_LOA_Loader> appLoaderRCPtr = createAppLoader();
  RefCountPtr<N_HB_Loader> hbLoaderRCPtr = createHBLoader(appLoaderRCPtr,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(hbLoaderRCPtr), false );
}

TEUCHOS_UNIT_TEST( SinCosLoader, createAppLoader ) {
  RefCountPtr<N_LOA_Loader> appLoader = createAppLoader();
  TEST_EQUALITY_CONST( Teuchos::is_null(appLoader), false );
}

TEUCHOS_UNIT_TEST( SinCosLoader, createAppBuilder ) {
  int numSolVars = 2; // This is the number of unknowns in the SinCosLoader
  RefCountPtr<N_LAS_Builder> appBuilder = createAppBuilder(numSolVars,0);
  TEST_EQUALITY_CONST( Teuchos::is_null(appBuilder), false );
}

/*
TEUCHOS_UNIT_TEST( SinCosLoader, applyDAEMatrices ) {
  int numSolVars = 2; // This is the number of unknowns in the SinCosLoader
  RefCountPtr<N_LOA_Loader> appLoader = createAppLoader();
  RefCountPtr<N_LAS_Builder> appBuilder = createAppBuilder(numSolVars,0);
  // Compute Jacobian through applyDAEMatrices
  RefCountPtr<N_LAS_MultiVector> applyJac = applyJacobian(*appLoader,*appBuilder);
  // Compute Jacobian through loadDAEMatrices
  RefCountPtr<N_LAS_MultiVector> loadJac = loadJacobian(*appLoader,*appBuilder);
  double tol = 1.0e-10;
  for (int i=0 ; i<numSolVars ; ++i) {
    for (int j=0 ; j<numSolVars ; ++j) {
      std::cout << "applyJac["<<i<<"]["<<j<<"] = " << (*applyJac)[i][j] << std::endl;
      std::cout << " loadJac["<<i<<"]["<<j<<"] = " << (*loadJac)[i][j] << std::endl;
      TEST_FLOATING_EQUALITY( (*applyJac)[i][j], (*loadJac)[i][j], tol);
    }
  }
}
*/

TEUCHOS_UNIT_TEST( SinCosLoader, diffJac ) {
  int numSolVars = 2; // This is the number of unknowns in the SinCosLoader
  RefCountPtr<N_LOA_Loader> appLoader = createAppLoader();
  RefCountPtr<N_LAS_Builder> appBuilder = createAppBuilder(numSolVars,0);
  // Compute Jacobian through applyDAEMatrices
  RefCountPtr<N_LAS_MultiVector> applyJac = applyJacobian(*appLoader,*appBuilder);
  // Compute Jacobian through diffJacobian
  RefCountPtr<N_LAS_MultiVector> diffJac = diffJacobian(*appLoader,*appBuilder);
  double tol = 1.0e-8;
  for (int i=0 ; i<numSolVars ; ++i) {
    for (int j=0 ; j<numSolVars ; ++j) {
//      std::cout << "applyJac["<<i<<"]["<<j<<"] = " << (*applyJac)[i][j] << std::endl;
//      std::cout << " diffJac["<<i<<"]["<<j<<"] = " << (*diffJac)[i][j] << std::endl;
      TEST_FLOATING_EQUALITY( (*applyJac)[i][j], (*diffJac)[i][j], tol);
    }
  }
}

/*
TEUCHOS_UNIT_TEST( N_HB_Loader, applyDAEMatrices ) {
  int blocks = 11; // must be odd for the FFT specific code to work or we'll get invalid reads & writes
  int numSolVars = 2; // This is the number of unknowns in the SinCosLoader
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr = createHBBuilder(blocks,numSolVars,0);
  bool matrixFreeFlag = true;
  RefCountPtr<N_LOA_Loader> appLoader = createAppLoader();
  RefCountPtr<N_HB_Loader> hbLoaderRCPtr = createHBLoader(appLoader, matrixFreeFlag);
  RefCountPtr<N_LAS_Builder> appBuilder = createAppBuilder(numSolVars,0);
  registerVectorsOnHBLoader(&*hbLoaderRCPtr,*hbBuilderRCPtr,*appBuilder);
  hbLoaderRCPtr->registerHBBuilder(hbBuilderRCPtr);
  // Set fast time points on the loader, lest we get a segfault when we try to access it.
  std::vector<double> fastTimes(blocks);
  for (int i=0 ; i<blocks ; ++i) {
    fastTimes[i] = i/blocks;
  }
  hbLoaderRCPtr->setFastTimes(fastTimes);
  // We need to set the N_MPDE_DeviceInterface on N_HB_Loader, lest we get a segfault when it tries to access it.
  RefCountPtr<N_MPDE_DeviceInterface> mpdeDevIntPtr;
  {
    mpdeDevIntPtr = rcp(new N_MPDE_DeviceInterface);
  }
  // We need to set an N_DEV_DeviceInterface ptr on N_MPDE_DeviceInterface, lest we get a segfault when it tries to access it.
  RefCountPtr<N_DEV_DeviceInterface> devIntPtr;
  {
    // But first, we need a valid N_IO_CmdParse object to instantiate the N_DEV_DeviceInterface object
    RefCountPtr<N_IO_CmdParse> cp = rcp(new N_IO_CmdParse);
    // This needs to be filled as in the builder case.
    devIntPtr = rcp(N_DEV_DeviceInterface::factory(*cp));
  }
  mpdeDevIntPtr->registerDeviceInterface(devIntPtr);
  hbLoaderRCPtr->registerMPDEDeviceInterface(mpdeDevIntPtr);


  RefCountPtr<N_LAS_MultiVector> applyJac = applyJacobian(*hbLoaderRCPtr,*hbBuilderRCPtr);
  RefCountPtr<N_LAS_MultiVector> diffJac = diffJacobian(*hbLoaderRCPtr,*hbBuilderRCPtr);

  // Compare this with the difference Jacobian:
  double tol = 1.0e-10;
  int N = applyJac->globalLength();
  for (int i=0; i<N ; ++i) {
    for (int j=0; j<N ; ++j) {
      cout << "applyJac["<<i<<"]["<<j<<"] = " << (*applyJac)[i][j] << endl;
      cout << "diffJac["<<i<<"]["<<j<<"] = " << (*diffJac)[i][j] << endl;
      TEST_FLOATING_EQUALITY( (*applyJac)[i][j], (*diffJac)[i][j], tol );
    }
  }
}
*/

