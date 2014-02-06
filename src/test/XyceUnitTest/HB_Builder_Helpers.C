//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: HB_Builder_Helpers.C,v $
// Purpose       : This file contains some helper functions for create N_HB_Builders.
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 9/10/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.2 $
// Revision Date  : $Date: 2008/09/18 17:04:27 $
// Current Owner  : $Author: tscoffe $
//-----------------------------------------------------------------------------

#include <HB_Builder_Helpers.h>

#include <N_IO_CmdParse.h>
#include <N_PDS_Manager.h>
#include <N_PDS_SerialComm.h>
#include <N_ERH_ErrorMgr.h>
#include <N_MPDE_Manager.h>
#include <N_MPDE_Discretization.h>
#include <N_HB_Builder.h>

#include <Epetra_CrsGraph.h>

RefCountPtr<N_HB_Builder> createHBBuilder(int numMPDEBlocks, int numSolutionVars, int numStateVars) {
  // We need to register a comm object with N_ERH_ErrorMgr or it will segfault on errors
  RefCountPtr<N_PDS_Comm> pdsCommRCPtr = rcp(new N_PDS_SerialComm);
  N_ERH_ErrorMgr::registerComm(pdsCommRCPtr);

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
  N_MPDE_Discretization::Type type = N_MPDE_Discretization::Backward;
  int order = 1;
  RefCountPtr<N_MPDE_Discretization> discRCPtr = rcp(new N_MPDE_Discretization(type,order) );
  RefCountPtr<N_HB_Builder> hbBuilderRCPtr;
  {
    int blocks = numMPDEBlocks;
    bool warpMPDEFlag = false;
    hbBuilderRCPtr = rcp(new N_HB_Builder( mpdeMgrRCPtr, blocks, discRCPtr, warpMPDEFlag));
  }
  // Now we have to call generateMaps or it will segfault on createVector
  RefCountPtr<Epetra_Map> baseMap;
  {
    int numGlobalElements = numSolutionVars;
    int IndexBase = 0;
    RefCountPtr<Epetra_Comm> comm = rcp(pdsCommRCPtr->petraComm(),false);
    baseMap = rcp(new Epetra_Map(numGlobalElements,IndexBase,*comm));
  }
  hbBuilderRCPtr->generateMaps(*baseMap);
  // Now we have to call generateHBMaps or it will also segfault on createVector
  hbBuilderRCPtr->generateHBMaps(*baseMap);
  // And we need to call generateStateMaps or it will also segfault on createStateVector
  RefCountPtr<Epetra_Map> stateMap;
  {
    int numGlobalElements = numStateVars;
    int IndexBase = 0;
    RefCountPtr<Epetra_Comm> comm = rcp(pdsCommRCPtr->petraComm(),false);
    stateMap = rcp(new Epetra_Map(numGlobalElements,IndexBase,*comm));
  }
  hbBuilderRCPtr->generateStateMaps(*stateMap);
  // We need to create a Epetra_CrsGraph so we can call generateGraphs
  RefCountPtr<Epetra_CrsGraph> matrixGraph;
  {
    matrixGraph = rcp(new Epetra_CrsGraph(View, *baseMap, numSolutionVars));
  }
  hbBuilderRCPtr->generateGraphs(*matrixGraph,*matrixGraph,*matrixGraph);

  return hbBuilderRCPtr;
}


