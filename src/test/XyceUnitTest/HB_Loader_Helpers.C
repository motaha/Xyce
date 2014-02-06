//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: HB_Loader_Helpers.C,v $
// Purpose       : This file contains some helper functions for create N_HB_Loaders.
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/10/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.3 $
// Revision Date  : $Date: 2009/12/07 18:06:47 $
// Current Owner  : $Author: elranki $
//-----------------------------------------------------------------------------


#include <HB_Loader_Helpers.h>
#include <HB_Builder_Helpers.h>

#include <N_IO_CmdParse.h>
#include <N_PDS_Manager.h>
#include <N_PDS_SerialComm.h>
#include <N_ERH_ErrorMgr.h>
#include <N_MPDE_Manager.h>
#include <N_MPDE_State.h>
#include <N_LOA_Loader.h>
#include <N_LAS_Builder.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_WarpedPhaseCondition.h>

#include <N_HB_Loader.h>
#include <N_HB_Builder.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>

RefCountPtr<N_HB_Loader> createHBLoader(RefCountPtr<N_LOA_Loader> appLoader, bool matrixFreeFlag) {
  // We need to register a comm object with N_ERH_ErrorMgr or it will segfault on errors
  RefCountPtr<N_PDS_Comm> pdsCommRCPtr = rcp(new N_PDS_SerialComm); 
  N_ERH_ErrorMgr::registerComm(pdsCommRCPtr);

  // Need an N_MPDE_Manager object
  N_IO_CmdParse cmdLine;
  // We need to register the N_PDS_Manager with N_IO_CmdParse or it will segfault
  RefCountPtr<N_PDS_Manager> pdsMgrRCPtr;
  {
    bool isSerial;
    bool procFlag;
    pdsMgrRCPtr = rcp(new N_PDS_Manager(isSerial,procFlag));
  }
  cmdLine.registerParallelMgr(pdsMgrRCPtr);
  // We need a valid cmdLine netlist or the MPDE_Manager will quit
  {
    int iargs = 2;
    char *arg0 = "Xyce";
    char *arg1 = "foo.cir";
    char *cargs[iargs];
    cargs[0] = arg0;
    cargs[1] = arg1;
    cmdLine.parseCommandLine(iargs,cargs);
  }
  RefCountPtr<N_MPDE_Manager> mpdeMgrRCPtr = rcp(new N_MPDE_Manager(cmdLine));
  mpdeMgrRCPtr->setMatrixFreeFlag(matrixFreeFlag);
  // Need an N_MPDE_Discretization object
  N_MPDE_Discretization::Type type = N_MPDE_Discretization::Backward;
  int order = 1;
  RefCountPtr<N_MPDE_Discretization> discRCPtr = rcp(new N_MPDE_Discretization(type,order) );

  // Need an N_MPDE_State object
  RefCountPtr<N_MPDE_State> mpdeState = rcp(new N_MPDE_State);
  // Need an N_MPDE_WarpedPhaseCondition object
  RefCountPtr<N_MPDE_WarpedPhaseCondition> warpPhaseRCPtr = rcp(new N_MPDE_WarpedPhaseCondition(0, 0.0, 0, 0, 0, 0, 0));
  // Now we can create the N_HB_Loader object
  RefCountPtr<N_HB_Loader> hbLoaderRCPtr = rcp(new N_HB_Loader(
        *mpdeState,
        discRCPtr,
        mpdeMgrRCPtr,
        warpPhaseRCPtr
        )
      );
  hbLoaderRCPtr->registerAppLoader(appLoader);
  return(hbLoaderRCPtr);
}

RefCountPtr<N_LOA_Loader> createAppLoader()
{
  RefCountPtr<SinCosLoader> appLoader = rcp(new SinCosLoader);
  return(appLoader);
}

RefCountPtr<N_LAS_Builder> createAppBuilder(int numSolutionVars, int numStateVars)
{
  // We need to register a comm object with N_ERH_ErrorMgr or it will segfault on errors
  RefCountPtr<N_PDS_Comm> pdsCommRCPtr = rcp(new N_PDS_SerialComm); 
  pdsCommRCPtr.release(); // memory leak!
  N_ERH_ErrorMgr::registerComm(pdsCommRCPtr);

  RefCountPtr<N_LAS_Builder> appBuilder = rcp(new N_LAS_Builder);
  // We need an N_PDS_Manager to get maps from.
  RefCountPtr<N_PDS_Manager> pdsMgrRCPtr;
  {
    bool isSerial;
    bool procFlag;
    pdsMgrRCPtr = rcp(new N_PDS_Manager(isSerial,procFlag));
  }
  //pdsMgrRCPtr.release(); // ??? memory leak?
  // We need an N_PDS_ParMap to register with the N_PDS_Manager
  // First, we need to create an Epetra_Map
  RefCountPtr<Epetra_Map> appMap;
  {
    int numGlobalElements = numSolutionVars;
    int IndexBase = 0;
    RefCountPtr<Epetra_Comm> comm = rcp(pdsCommRCPtr->petraComm(),false);
    appMap = rcp(new Epetra_Map(numGlobalElements,IndexBase,*comm));
  }
  appMap.release(); 
  RefCountPtr<N_PDS_ParMap> appPDSParMap = rcp(new N_PDS_ParMap(&*appMap,&*pdsCommRCPtr));
  appPDSParMap.release(); 
  /*
  // Print out the global ids:
  {
    int N = appMap->NumMyPoints();
    int * myGID = appMap->MyGlobalElements();
    cout << endl;
    for (int i=0; i<N; ++i) {
      cout << "myGID["<<i<<"] = "<< myGID[i] << endl;
    }
  }
  */
  RefCountPtr<N_PDS_ParMap> appPDSParMap2 = rcp(new N_PDS_ParMap(&*appMap,&*pdsCommRCPtr));
  appPDSParMap2.release(); 

  RefCountPtr<Epetra_Map> appStateMap;
  {
    int numGlobalElements = numStateVars;
    int IndexBase = 0;
    RefCountPtr<Epetra_Comm> comm = rcp(pdsCommRCPtr->petraComm(),false);
    appStateMap = rcp(new Epetra_Map(numGlobalElements,IndexBase,*comm));
  }
  appStateMap.release(); 
  RefCountPtr<N_PDS_ParMap> appStatePDSParMap = rcp(new N_PDS_ParMap(&*appStateMap,&*pdsCommRCPtr));
  appStatePDSParMap.release(); 

  RefCountPtr<N_PDS_ParMap> appStatePDSParMap2 = rcp(new N_PDS_ParMap(&*appStateMap,&*pdsCommRCPtr));
  appStatePDSParMap2.release(); 

  pdsMgrRCPtr->addParallelMap("SOLUTION",&*appPDSParMap);
  pdsMgrRCPtr->addParallelMap("SOLUTION_OVERLAP_GND",&*appPDSParMap2);
  pdsMgrRCPtr->addParallelMap("STATE",&*appStatePDSParMap);
  pdsMgrRCPtr->addParallelMap("STATE_OVERLAP",&*appStatePDSParMap2);
  // And we need a couple matrix graphs
  RefCountPtr<Epetra_CrsGraph> matrixGraph = rcp(new Epetra_CrsGraph(View, *appMap, numSolutionVars));
  matrixGraph.release(); 

  RefCountPtr<Epetra_CrsGraph> matrixGraph2 = rcp(new Epetra_CrsGraph(View, *appMap, numSolutionVars));
  matrixGraph2.release(); 
  RefCountPtr<Epetra_CrsGraph> matrixGraph3 = rcp(new Epetra_CrsGraph(View, *appMap, numSolutionVars));
  matrixGraph3.release(); 
  RefCountPtr<Epetra_CrsGraph> matrixGraph4 = rcp(new Epetra_CrsGraph(View, *appMap, numSolutionVars));
  matrixGraph4.release(); 

  pdsMgrRCPtr->addMatrixGraph("DAE_DQDX_JAC_OVERLAP_GND",&*matrixGraph);
  pdsMgrRCPtr->addMatrixGraph("DAE_DQDX_JAC",&*matrixGraph2);
  pdsMgrRCPtr->addMatrixGraph("DAE_DFDX_JAC_OVERLAP_GND",&*matrixGraph3);
  pdsMgrRCPtr->addMatrixGraph("DAE_DFDX_JAC",&*matrixGraph4);
  appBuilder->registerPDSManager(pdsMgrRCPtr);
  return(appBuilder);
}

void registerVectorsOnHBLoader(N_HB_Loader* hbLoader, const N_HB_Builder& hbBuilder, const N_LAS_Builder& appBuilder)
{
  hbLoader->registerAppVec ( rcp(appBuilder.createVector()) );
  hbLoader->registerAppNextStaVec ( rcp(appBuilder.createStateVector()) );
  hbLoader->registerAppCurrStaVec ( rcp(appBuilder.createStateVector()) );
  hbLoader->registerAppLastStaVec ( rcp(appBuilder.createStateVector()) );
  hbLoader->registerOmegadQdt2( Teuchos::rcp_dynamic_cast<N_LAS_BlockVector>(rcp(hbBuilder.createVector())) );
  hbLoader->registerAppdQdx( rcp(appBuilder.createDAEdQdxMatrix()) );
  hbLoader->registerAppdFdx( rcp(appBuilder.createDAEdFdxMatrix()) );
  hbLoader->registerMPDEdQdx
    ( Teuchos::rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(hbBuilder.createDAEdQdxMatrix())) );
  hbLoader->registerMPDEdFdx
    ( Teuchos::rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(hbBuilder.createDAEdFdxMatrix())) );
  hbLoader->registerXt(hbBuilder.createTimeDomainBlockVector());
  hbLoader->registerVt(hbBuilder.createTimeDomainBlockVector());
  hbLoader->registerHBdQdx
    ( Teuchos::rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(hbBuilder.createDAEdQdxMatrix())) );
  hbLoader->registerHBdFdx
    ( Teuchos::rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(hbBuilder.createDAEdFdxMatrix())) );
}

RefCountPtr<N_LAS_MultiVector> diffJacobian(N_LOA_Loader& loader, const N_LAS_Builder& builder)
{
  // Call loader.loadDAEVectors to difference the residual to fill a Jacobian 
  RefCountPtr<const Epetra_Map> emap = builder.getSolutionMap();
  int N = emap->NumGlobalElements();
  RefCountPtr<Epetra_MultiVector> ejac = rcp(new Epetra_MultiVector(*emap,N),false);
  RefCountPtr<N_LAS_MultiVector> jac = rcp(new N_LAS_MultiVector(&*ejac,true), true); 
  // Temporary vectors for loads.
  RefCountPtr<N_LAS_Vector> solVector = rcp(builder.createVector());
 RefCountPtr<N_LAS_Vector> currSolVector = rcp(builder.createStateVector());
 RefCountPtr<N_LAS_Vector> lastSolVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Vector> staVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Vector> currStaVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Vector> lastStaVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Vector> staDerivVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Vector> qVector = rcp(builder.createVector());
  RefCountPtr<N_LAS_Vector> fVector = rcp(builder.createVector());
  RefCountPtr<N_LAS_Vector> qVoltLimVector = rcp(builder.createVector());
  RefCountPtr<N_LAS_Vector> fVoltLimVector = rcp(builder.createVector());
  double h = 1.0e-6;
  bool status;
  for (int i=0; i<N; ++i) {
    // Evaluate at identity+h vector
    solVector->putScalar(0.0);
    (*solVector)[i] = 1.0+h;
    status = loader.loadDAEVectors(
        &*solVector,
       &*currSolVector,
       &*lastSolVector,
        &*staVector,
        &*currStaVector,
        &*lastStaVector,
        &*staDerivVector,
        &*qVector,
        &*fVector,
        &*fVoltLimVector,
        &*qVoltLimVector
        ); 
    // Copy value into jac
    RefCountPtr<N_LAS_Vector> jacView = jac->getNonConstVectorView(i);
    *jacView = *fVector;
    // Evaluate at identity vector
    (*solVector)[i] = 1.0;
    status = loader.loadDAEVectors(
        &*solVector,
       &*currSolVector,
       &*lastSolVector,
        &*staVector,
        &*currStaVector,
        &*lastStaVector,
        &*staDerivVector,
        &*qVector,
        &*fVector,
        &*fVoltLimVector,
        &*qVoltLimVector
        ); 
    // Difference to get gradient 
    jacView->update(-1.0/h,*fVector,1.0/h);
  } 
  return jac;
}

RefCountPtr<N_LAS_MultiVector> loadJacobian(N_LOA_Loader& loader, const N_LAS_Builder& builder)
{
  // Call loader.loadDAEVectors to difference the residual to fill a Jacobian 
  RefCountPtr<const Epetra_Map> emap = builder.getSolutionMap();
  int N = emap->NumGlobalElements();
  RefCountPtr<Epetra_MultiVector> ejac = rcp(new Epetra_MultiVector(*emap,N),false);
  RefCountPtr<N_LAS_MultiVector> jac = rcp(new N_LAS_MultiVector(&*ejac,true), true); 
  // Temporary vectors for loads.
  RefCountPtr<N_LAS_Vector> solVector = rcp(builder.createVector());
  RefCountPtr<N_LAS_Vector> staVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Vector> staDerivVector = rcp(builder.createStateVector());
  RefCountPtr<N_LAS_Matrix> dQdx = rcp(builder.createDAEdQdxMatrix());
  RefCountPtr<N_LAS_Matrix> dFdx = rcp(builder.createDAEdFdxMatrix());
  bool status = loader.loadDAEMatrices( &*solVector, &*staVector, &*staDerivVector, &*dQdx, &*dFdx );
  // Copy data out of the N_LAS_Matrix into the N_LAS_MultiVector
  int length = N;
  double coeffs[N];
  int colIndices[N];
  for (int row=0 ; row<N ; ++row) {
    int numEntries = 0;
    dFdx->getRowCopy( row, length, numEntries, coeffs, colIndices);
    for (int col=0 ; col<numEntries ; ++col) {
      (*jac)[row][colIndices[col]] = coeffs[col];
    }
  }
  return jac;
}

RefCountPtr<N_LAS_MultiVector> applyJacobian(N_LOA_Loader& loader, const N_LAS_Builder& builder)
{
  // Call loader.applyDAEMatrices to fill a Jacobian matrix
  RefCountPtr<const Epetra_Map> emap = builder.getSolutionMap();
  int N = emap->NumGlobalElements();
  RefCountPtr<Epetra_MultiVector> ejac = rcp(new Epetra_MultiVector(*emap,N),false);
  RefCountPtr<N_LAS_MultiVector> jac = rcp(new N_LAS_MultiVector(&*ejac,true), true); 
  // Temporary vectors for loads.
  RefCountPtr<N_LAS_Vector> X = rcp(builder.createVector(0.0));
  RefCountPtr<N_LAS_Vector> S = rcp(builder.createStateVector(0.0));
  RefCountPtr<N_LAS_Vector> dSdt = rcp(builder.createStateVector(0.0));
  RefCountPtr<N_LAS_Vector> dQdxV = rcp(builder.createVector(0.0));
  RefCountPtr<N_LAS_Vector> dFdxV = rcp(builder.createVector(0.0));

  // Create a vector to contain columns of the identity matrix as input
  RefCountPtr<N_LAS_Vector> V = rcp(builder.createVector(0.0));

  // loop over NumGlobalElements calling applyDAEMatrices to each column of identity matrix
  for (int i=0 ; i<N ; ++i) {
    // Fill V with i'th column of identity matrix
    V->putScalar(0.0);
    (*V)[i] = 1.0;
    // Call applyDAEMatrices
    loader.applyDAEMatrices(&*X,&*S,&*dSdt,*V,&*dQdxV,&*dFdxV);
    // Get view into ith vector in N_LAS_MultiVector matrix
    RefCountPtr<N_LAS_Vector> matView = jac->getNonConstVectorView(i);
    // Copy dFdxV into matView
    *matView = *dFdxV;
  }
  return jac;
}

