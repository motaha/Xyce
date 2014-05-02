//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_Builder.C,v $
//
// Purpose        : Builder class for parallel/serial linear objects
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/22/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.64 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <Epetra_MapColoring.h>

#include <EpetraExt_View_CrsGraph.h>

#include <N_UTL_Functors.h>
#include <N_UTL_fwd.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_GlobalAccessor.h>

#include <N_LAS_Builder.h>
#include <N_LAS_LAFactory.h>
#include <N_LAS_QueryUtil.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <N_IO_CmdParse.h>

#include <N_ERH_ErrorMgr.h>

using std::max;

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::~
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_LAS_Builder::~N_LAS_Builder()
{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createMultiVector
// Purpose       : returns Soln/RHS sized multiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector * N_LAS_Builder::createMultiVector( const int numVectors, const double value ) const
{
    return N_LAS_LAFactory::newMultiVector( value,
                                          *(pdsMgr_->getParallelMap("SOLUTION")),
                                          numVectors,
                                          *(pdsMgr_->getParallelMap("SOLUTION_OVERLAP_GND")) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createStateMultiVector
// Purpose       : returns State sized multiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
N_LAS_MultiVector * N_LAS_Builder::createStateMultiVector( const int numVectors, const double value ) const
{
    return N_LAS_LAFactory::newMultiVector( value,
                                          *(pdsMgr_->getParallelMap("STATE")),
                                          numVectors,
                                          *(pdsMgr_->getParallelMap("STATE_OVERLAP")) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createStoreMultiVector
// Purpose       : returns Store sized multiVector
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
N_LAS_MultiVector * N_LAS_Builder::createStoreMultiVector( const int numVectors, const double value ) const
{
    return N_LAS_LAFactory::newMultiVector( value,
                                          *(pdsMgr_->getParallelMap("STORE")),
                                          numVectors,
                                          *(pdsMgr_->getParallelMap("STORE_OVERLAP")) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createVector
// Purpose       : returns Soln/RHS sized vector
// Special Notes : Takes an initial value argument.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_Builder::createVector( const double value ) const
{
    return N_LAS_LAFactory::newVector( value,
                                     *(pdsMgr_->getParallelMap("SOLUTION")),
                                     *(pdsMgr_->getParallelMap("SOLUTION_OVERLAP_GND")) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createStateVector
// Purpose       : returns State sized vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_Builder::createStateVector( const double value ) const
{
    return N_LAS_LAFactory::newVector( value,
                                     *(pdsMgr_->getParallelMap("STATE")),
                                     *(pdsMgr_->getParallelMap("STATE_OVERLAP")) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createStoreVector
// Purpose       : returns Store sized vector
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_Builder::createStoreVector( const double value ) const
{
    return N_LAS_LAFactory::newVector( value,
                                     *(pdsMgr_->getParallelMap("STORE")),
                                     *(pdsMgr_->getParallelMap("STORE_OVERLAP")) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createMatrix
// Purpose       : returns Matrix initialized based on QueryUtil and ParMap
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_LAS_Builder::createMatrix(const double initialValue) const
{

  N_LAS_Matrix * mat = 0;

  Epetra_CrsGraph * overlapGraph = pdsMgr_->getMatrixGraph( "JACOBIAN_OVERLAP_GND" );
  Epetra_CrsGraph * baseGraph = pdsMgr_->getMatrixGraph( "JACOBIAN" );

  mat = new N_LAS_Matrix( overlapGraph, baseGraph );

  return mat;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createSolnColoring
// Purpose       : Color Map representing variable types in solution vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/08/04
//-----------------------------------------------------------------------------
Epetra_MapColoring * N_LAS_Builder::createSolnColoring() const
{
  const std::vector<char> & charColors = lasQueryUtil_->rowList_VarType();

  int size = charColors.size();
  std::vector<int> colors( size );
  for( int i = 0; i < size; ++i )
  {
    switch( charColors[i] )
    {
      case 'V': colors[i] = 0;
                break;
      case 'I': colors[i] = 1;
                break;
      default : colors[i] = 2;
                break;
    }
  }

  return new Epetra_MapColoring( *(pdsMgr_->getParallelMap("SOLUTION")->petraBlockMap()),
                                 &colors[0], 0 );
}


//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::createInitialConditionColoring
// Purpose       :
// Special Notes : The .IC and .NODESET capabilities will use the variables
//                 which are colored 0.  This will be all voltage nodes not
//                 connected to independent sources.
// Scope         : Public
// Creator       : Eric R. Keiter,  SNL
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
Epetra_MapColoring * N_LAS_Builder::createInitialConditionColoring() const
{
  const std::vector<char> & charColors = lasQueryUtil_->rowList_VarType();
  const std::vector<int> & vsrcGIDColors = lasQueryUtil_->vsrcGIDVec();

  int size = charColors.size();
  int vsrcSize = vsrcGIDColors.size();
  std::vector<int> colors( size );
  for( int i = 0; i < size; ++i )
  {
    switch( charColors[i] )
    {
      case 'V': colors[i] = 0;
                break;
      case 'I': colors[i] = 1;
                break;
      default : colors[i] = 2;
                break;
    }
  }

  N_PDS_ParMap * solnMap = pdsMgr_->getParallelMap("SOLUTION");
  for( int i=0; i < vsrcSize; ++i )
  {
    int vsrcID = vsrcGIDColors[i];
    // Convert the ID from local to global if it is valid and the build is parallel.
    if (vsrcID >= 0)
    {
#ifdef Xyce_PARALLEL_MPI
      vsrcID = solnMap->globalToLocalIndex( vsrcGIDColors[i] );
#endif
      if (vsrcID < size && vsrcID >= 0)
        colors[vsrcID] = 1;
    }
  }

  return new Epetra_MapColoring( *(solnMap->petraBlockMap()), &colors[0], 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::generateParMaps
// Purpose       : Creates parallel maps for SOLN and STATE
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/24/01
//-----------------------------------------------------------------------------
bool N_LAS_Builder::generateParMaps()
{
  int numLocalRows = lasQueryUtil_->numLocalRows();
  std::vector<int> arrayGIDs(numLocalRows);
  arrayGIDs = lasQueryUtil_->rowList_GID();

  int numLocalStateVars = lasQueryUtil_->numLocalStateVars();
  std::vector<int> arrayStateGIDs(numLocalStateVars);
  arrayStateGIDs = lasQueryUtil_->rowList_StateGID();

  int numLocalStoreVars = lasQueryUtil_->numLocalStoreVars();
  std::vector<int> arrayStoreGIDs(numLocalStoreVars);
  arrayStoreGIDs = lasQueryUtil_->rowList_StoreGID();

  int numGlobalRows = lasQueryUtil_->numGlobalRows();
  N_PDS_ParMap * solnMap = pdsMgr_->createParallelMap(
    numGlobalRows, numLocalRows, arrayGIDs);

  int numGlobalStateVars = lasQueryUtil_->numGlobalStateVars();
  N_PDS_ParMap * stateMap = pdsMgr_->createParallelMap(
    numGlobalStateVars, numLocalStateVars, arrayStateGIDs);

  int numGlobalStoreVars = lasQueryUtil_->numGlobalStoreVars();
  N_PDS_ParMap * storeMap = pdsMgr_->createParallelMap(
    numGlobalStoreVars, numLocalStoreVars, arrayStoreGIDs);

  pdsMgr_->addParallelMap("SOLUTION", solnMap);
  pdsMgr_->addParallelMap("STATE", stateMap);
  pdsMgr_->addParallelMap("STORE", storeMap);

  int procCnt = pdsMgr_->getPDSComm()->numProc();
  int numExternRows = lasQueryUtil_->numExternRows();
  int totalRows = numLocalRows + numExternRows;

  int size = lasQueryUtil_->rowList_ExternGID().size();
  for(int i = 0; i < size; ++i)
    arrayGIDs.push_back(lasQueryUtil_->rowList_ExternGID()[i].first);

  int numGlobalExternRows = lasQueryUtil_->numGlobalExternRows();
  int totalGlobalRows = numGlobalRows + numGlobalExternRows;
  N_PDS_ParMap * overlapSolnMap = pdsMgr_->createParallelMap(
	totalGlobalRows,
        totalRows,
        arrayGIDs );

  arrayGIDs.push_back(-1);
  totalGlobalRows += procCnt;
  totalRows++;
  N_PDS_ParMap * overlapGndSolnMap = pdsMgr_->createParallelMap(
	totalGlobalRows,
        totalRows,
        arrayGIDs,
        -1 );

  int numExternStateVars = lasQueryUtil_->numExternStateVars();
  int totalStateVars = numLocalStateVars + numExternStateVars;
  arrayStateGIDs.reserve(totalStateVars);

  size = lasQueryUtil_->rowList_ExternStateGID().size();
  for( int i = 0; i < size; ++i )
    arrayStateGIDs.push_back(lasQueryUtil_->rowList_ExternStateGID()[i].first);

  int numGlobalExternStateVars = lasQueryUtil_->numGlobalExternStateVars();
  int totalGlobalStateVars = numGlobalStateVars + numGlobalExternStateVars;
  N_PDS_ParMap * overlapStateMap = pdsMgr_->createParallelMap(
        totalGlobalStateVars,
        totalStateVars,
        arrayStateGIDs );

  int numExternStoreVars = lasQueryUtil_->numExternStoreVars();
  int totalStoreVars = numLocalStoreVars + numExternStoreVars;
  arrayStoreGIDs.reserve(totalStoreVars);

  size = lasQueryUtil_->rowList_ExternStoreGID().size();
  for( int i = 0; i < size; ++i )
    arrayStoreGIDs.push_back(lasQueryUtil_->rowList_ExternStoreGID()[i].first);

  int numGlobalExternStoreVars = lasQueryUtil_->numGlobalExternStoreVars();
  int totalGlobalStoreVars = numGlobalStoreVars + numGlobalExternStoreVars;
  N_PDS_ParMap * overlapStoreMap = pdsMgr_->createParallelMap(
        totalGlobalStoreVars,
        totalStoreVars,
        arrayStoreGIDs );

  pdsMgr_->addParallelMap( "SOLUTION_OVERLAP", overlapSolnMap );
  pdsMgr_->addParallelMap( "SOLUTION_OVERLAP_GND", overlapGndSolnMap );
  pdsMgr_->addParallelMap( "STATE_OVERLAP", overlapStateMap );
  pdsMgr_->addParallelMap( "STORE_OVERLAP", overlapStoreMap );

  N_PDS_GlobalAccessor * solnGA = pdsMgr_->createGlobalAccessor( "SOLUTION" );
  solnGA->registerExternGIDVector( lasQueryUtil_->rowList_ExternGID() );
  solnGA->generateMigrationPlan();

  N_PDS_GlobalAccessor * stateGA = pdsMgr_->createGlobalAccessor( "STATE" );
  stateGA->registerExternGIDVector( lasQueryUtil_->rowList_ExternStateGID() );
  stateGA->generateMigrationPlan();

  N_PDS_GlobalAccessor * storeGA = pdsMgr_->createGlobalAccessor( "STORE" );
  storeGA->registerExternGIDVector( lasQueryUtil_->rowList_ExternStoreGID() );
  storeGA->generateMigrationPlan();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::generateGraphs
// Purpose       : Generation of Matrix Graphs, stored with ParMgr
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/02
//-----------------------------------------------------------------------------
bool N_LAS_Builder::generateGraphs()
{
  const std::vector< std::list<int> > & rcData = lasQueryUtil_->rowList_ColList();
  int numLocalRows_Overlap = rcData.size();

  int maxNZs = 0;
  std::vector<int> arrayNZs( numLocalRows_Overlap );
  for( int i = 0; i < numLocalRows_Overlap; ++i )
  {
    arrayNZs[i] = rcData[i].size();
    maxNZs = max( maxNZs, arrayNZs[i] );
  }

  Epetra_BlockMap & overlapGndMap = *(pdsMgr_->getParallelMap( "SOLUTION_OVERLAP_GND" )->petraBlockMap());
  Epetra_BlockMap & overlapMap = *(pdsMgr_->getParallelMap( "SOLUTION_OVERLAP" )->petraBlockMap());
  Epetra_BlockMap & localMap = *(pdsMgr_->getParallelMap( "SOLUTION" )->petraBlockMap());

  Epetra_CrsGraph * overlapGndGraph = new Epetra_CrsGraph( Copy, overlapGndMap, &arrayNZs[0] );
  Epetra_CrsGraph * overlapGraph = new Epetra_CrsGraph( Copy, overlapMap, &arrayNZs[0] );

  std::vector<int> indices( maxNZs );
  for( int i = 0; i < numLocalRows_Overlap; ++i )
  {
    copy( rcData[i].begin(), rcData[i].end(), indices.begin() );
    if( arrayNZs[i] )
      overlapGndGraph->InsertGlobalIndices( overlapGndMap.GID(i), arrayNZs[i], &indices[0] );

    if( overlapGndMap.GID(i) != -1 && arrayNZs[i] )
    {
      if( indices[0] == -1 )
        overlapGraph->InsertGlobalIndices( overlapMap.GID(i), arrayNZs[i]-1, &indices[1] );
      else
        overlapGraph->InsertGlobalIndices( overlapMap.GID(i), arrayNZs[i], &indices[0] );
    }
  }
  overlapGndGraph->FillComplete();
  overlapGraph->FillComplete();
  pdsMgr_->addMatrixGraph( "JACOBIAN_OVERLAP_GND", overlapGndGraph );
  pdsMgr_->addMatrixGraph( "JACOBIAN_OVERLAP", overlapGraph );

#ifndef Xyce_PARALLEL_MPI
  //determine what the column map should be for the local view of the graph
  //store a new one if it is different than the overlap map
  const Epetra_BlockMap & oColMap = overlapGraph->ColMap();
  int oColMapSize = oColMap.NumMyElements();
  std::vector<int> nColMapElements;
  for( int i = 0; i < oColMapSize; ++i )
    if( overlapMap.MyGID( oColMap.GID(i) ) ) nColMapElements.push_back( oColMap.GID(i) );

  N_PDS_ParMap * localColMap = 0;

  double lSize = static_cast<double>(nColMapElements.size());
  double tSize;
  pdsMgr_->getPDSComm()->sumAll( &lSize, &tSize, 1 );
  if( tSize < overlapMap.NumGlobalElements() )
  {
    int totalSize = static_cast<int>(tSize);
    localColMap = pdsMgr_->createParallelMap(
        totalSize,
        nColMapElements.size(),
        nColMapElements );
    pdsMgr_->addParallelMap( "SOLUTION_COLUMN", localColMap );
  }

  Epetra_BlockMap * nColMap = &overlapMap;
  if( localColMap ) nColMap = localColMap->petraBlockMap();

  EpetraExt::CrsGraph_View * viewTransform = new EpetraExt::CrsGraph_View( &localMap, nColMap );
  Epetra_CrsGraph * localGraph = &((*viewTransform)( *overlapGraph ));
  pdsMgr_->addMatrixGraph( "JACOBIAN", localGraph, viewTransform );
#else
  if( pdsMgr_->getPDSComm()->isSerial() )
  {
    //determine what the column map should be for the local view of the graph
    //store a new one if it is different than the overlap map
    const Epetra_BlockMap & oColMap = overlapGraph->ColMap();
    int oColMapSize = oColMap.NumMyElements();
    std::vector<int> nColMapElements;
    for( int i = 0; i < oColMapSize; ++i )
      if( overlapMap.MyGID( oColMap.GID(i) ) ) nColMapElements.push_back( oColMap.GID(i) );

    N_PDS_ParMap * localColMap = 0;

    double lSize = static_cast<double>(nColMapElements.size());
    double tSize;
    pdsMgr_->getPDSComm()->sumAll( &lSize, &tSize, 1 );
    if( tSize < overlapMap.NumGlobalElements() )
    {
      int totalSize = static_cast<int>(tSize);
      localColMap = pdsMgr_->createParallelMap(
          totalSize,
          nColMapElements.size(),
          nColMapElements );
      pdsMgr_->addParallelMap( "SOLUTION_COLUMN", localColMap );
    }

    Epetra_BlockMap * nColMap = &overlapMap;
    if( localColMap ) nColMap = localColMap->petraBlockMap();

    EpetraExt::CrsGraph_View * viewTransform = new EpetraExt::CrsGraph_View( &localMap, nColMap );
    Epetra_CrsGraph * localGraph = &((*viewTransform)( *overlapGraph ));
    pdsMgr_->addMatrixGraph( "JACOBIAN", localGraph, viewTransform );
  }
  else
  {
    Epetra_Export exporter( overlapMap, localMap );
    Epetra_CrsGraph * localGraph = new Epetra_CrsGraph( Copy, localMap, 0 );
    localGraph->Export( *overlapGraph, exporter, Add );
#ifdef Xyce_VERBOSE_LINEAR
    Xyce::lout() << "Exported To Local Graph!\n"  << std::endl;
#endif
    localGraph->FillComplete();
#ifdef Xyce_VERBOSE_LINEAR
    Xyce::lout() << "Local Graph Transformed!\n"  << std::endl;
#endif
    pdsMgr_->addMatrixGraph( "JACOBIAN", localGraph );
  }
#endif


  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getSolutionMap
// Purpose       :
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<const Epetra_Map> N_LAS_Builder::getSolutionMap() const
{
  return(Teuchos::rcp(pdsMgr_->getParallelMap("SOLUTION")->petraMap(),false));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getSolutionMaps
// Purpose       :
// Special Notes : This is specifically for ModelEvaluator Interface
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/02/09
//-----------------------------------------------------------------------------
void N_LAS_Builder::getSolutionMaps(
    RCP<N_PDS_ParMap>& soln_map,
    RCP<N_PDS_ParMap>& soln_ognd_map
    ) const
{
  soln_map = Teuchos::rcp(pdsMgr_->getParallelMap("SOLUTION"),false);
  soln_ognd_map = Teuchos::rcp(pdsMgr_->getParallelMap("SOLUTION_OVERLAP_GND"),false);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getStateMap
// Purpose       :
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<const Epetra_Map> N_LAS_Builder::getStateMap() const
{
  return(Teuchos::rcp(pdsMgr_->getParallelMap("STATE")->petraMap(),false));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getStoreMap
// Purpose       :
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
RCP<const Epetra_Map> N_LAS_Builder::getStoreMap() const
{
  return(Teuchos::rcp(pdsMgr_->getParallelMap("STORE")->petraMap(),false));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getStateMap
// Purpose       :
// Special Notes : This is specifically for ModelEvaluator Interface
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/02/09
//-----------------------------------------------------------------------------
void N_LAS_Builder::getStateMap(
    RCP<N_PDS_ParMap> & state_map
    ) const
{
  state_map = Teuchos::rcp(pdsMgr_->getParallelMap("STATE"),false);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getStoreMap
// Purpose       :
// Special Notes : This is specifically for ModelEvaluator Interface
// Scope         : Public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_LAS_Builder::getStoreMap(
    RCP<N_PDS_ParMap> & store_map
    ) const
{
  store_map = Teuchos::rcp(pdsMgr_->getParallelMap("STORE"),false);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getdFdxGraphs
// Purpose       :
// Special Notes : This is specifically for ModelEvaluator Interface
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/02/09
//-----------------------------------------------------------------------------
void N_LAS_Builder::getdFdxGraphs(
    RCP<Epetra_CrsGraph>& dFdx_graph,
    RCP<Epetra_CrsGraph>& dFdx_ognd_graph
    ) const
{
  dFdx_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( "DAE_DFDX_JAC" ));
  dFdx_ognd_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( "DAE_DFDX_JAC_OVERLAP_GND" ));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_Builder::getdQdxGraphs
// Purpose       :
// Special Notes : This is specifically for ModelEvaluator Interface
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/02/09
//-----------------------------------------------------------------------------
void N_LAS_Builder::getdQdxGraphs(
    RCP<Epetra_CrsGraph>& dQdx_graph,
    RCP<Epetra_CrsGraph>& dQdx_ognd_graph
    ) const
{
  dQdx_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( "DAE_DQDX_JAC" ));
  dQdx_ognd_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( "DAE_DQDX_JAC_OVERLAP_GND" ));
}


bool N_LAS_Builder::registerPDSManager(N_PDS_Manager * PDS_Manager)
{
  pdsMgr_ = rcp(PDS_Manager,false); return true;
}

bool N_LAS_Builder::registerPDSManager(RCP<N_PDS_Manager> PDS_Manager)
{
  pdsMgr_ = PDS_Manager; return true;
}

