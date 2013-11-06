//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename      : $RCSfile: N_MPDE_Builder.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.38.6.2 $
// Revision Date  : $Date: 2013/10/03 17:23:47 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_MPDE_Builder.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_Manager.h>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>

#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_PDS_ParMap.h>

#include <N_ERH_ErrorMgr.h>


#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::~N_MPDE_Builder
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_MPDE_Builder::~N_MPDE_Builder()
{
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Vector * N_MPDE_Builder::createVector( double initialValue ) const
{
  if (warpMPDE_)
  {
    // tscoffe/tmei 08/11/05:  Appending an extra row for omega and phi
    return dynamic_cast<N_LAS_Vector*>(
          new N_LAS_BlockVector( Size_, MPDEMap_, BaseMap_, 2 ) );
  }
  else
  {
    return dynamic_cast<N_LAS_Vector*>(
          new N_LAS_BlockVector( Size_, MPDEMap_, BaseMap_ ) );
  }
}



//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
N_LAS_Vector * N_MPDE_Builder::createStateVector( double initialValue ) const
{
  return dynamic_cast<N_LAS_Vector*>(
        new N_LAS_BlockVector( Size_, MPDEStateMap_, BaseStateMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
N_LAS_Vector * N_MPDE_Builder::createStoreVector( double initialValue ) const
{
  return dynamic_cast<N_LAS_Vector*>(
        new N_LAS_BlockVector( Size_, MPDEStoreMap_, BaseStoreMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createDAEdQdxMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_MPDE_Builder::createDAEdQdxMatrix( double initialValue ) const
{
  return createDAEdFdxMatrix( initialValue );
  
  /*  NOTE:  The dQdx matrix is created with the same graph as the dFdx so
             that the construction of the Jacobian is faster
  
  vector< vector<int> > Cols(Size_);
  for( int i = 0; i < Size_; ++i )
  {
    Cols[i].resize(1);
    Cols[i][0] = i;
  }
  if (warpMPDE_)
  {
    // tscoffe/tmei 08/11/05:  Appending an extra row & column for omega and phi
    return dynamic_cast<N_LAS_Matrix*>(
          new N_LAS_BlockMatrix( Size_,
                                 Cols,
                                 *MPDEdQdxGraph_,
                                 *BasedQdxGraph_,
                                 2) );
  }
  else
  {
    return dynamic_cast<N_LAS_Matrix*>(
          new N_LAS_BlockMatrix( Size_,
                                 Cols,
                                 *MPDEdQdxGraph_,
                                 *BasedQdxGraph_ ) );
  }
  */
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createDAEdFdxMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_MPDE_Builder::createDAEdFdxMatrix( double initialValue ) const
{
  vector< vector<int> > Cols(Size_);
  int Start = Disc_->Start();
  int Width = Disc_->Width();

  for( int i = 0; i < Size_; ++i )
  {
    Cols[i].resize(Width);
    for( int j = 0; j < Width; ++j )
    {
      int Loc = i+(j+Start);
      if( Loc < 0 )           Loc += Size_;
      else if( Loc >= Size_ ) Loc -= Size_;
      Cols[i][j] = Loc;
    }
    sort(Cols[i].begin(),Cols[i].end());
  }

  if (warpMPDE_)
  {
    // tscoffe/tmei 08/11/05:  Appending an extra row & column for omega and phi
    return dynamic_cast<N_LAS_Matrix*>(
          new N_LAS_BlockMatrix( Size_,
                                 Cols,
                                 *MPDEdFdxGraph_,
                                 *BasedFdxGraph_,
                                 2) );
  }
  else
  {
    return dynamic_cast<N_LAS_Matrix*>(
          new N_LAS_BlockMatrix( Size_,
                                 Cols,
                                 *MPDEdFdxGraph_,
                                 *BasedFdxGraph_ ) );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createDAEFullMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_MPDE_Builder::createDAEFullMatrix( double initialValue ) const
{
  return createDAEdFdxMatrix( initialValue );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateMaps( const RCP<N_PDS_ParMap>& BaseMap )
{
  //Save copy of base map
  BaseMap_ = BaseMap;
  
  const Epetra_Comm & Comm  = BaseMap_->petraMap()->Comm();

  //determine block offset
  int MaxGID = BaseMap->petraMap()->MaxAllGID();
  int BaseIndex = BaseMap->petraMap()->IndexBase();
  offset_ = 1;
  while( offset_ <= MaxGID ) offset_ *= 10;

  //Setup Block Maps and MPDE Map
  int BlockSize = BaseMap->petraMap()->NumGlobalElements();
  vector<int> BlockGIDs(BlockSize);
  vector<int> MPDEGIDs;
  if (warpMPDE_)
  {
    // tscoffe/tmei 08/02/05:  Added two to size of vector for omega and phi
    MPDEGIDs.resize(Size_*BlockSize+2);
  }
  else
  {
    MPDEGIDs.resize(Size_*BlockSize);
  }
  vector<int> BaseGIDs(BlockSize);

  BaseMap->petraMap()->MyGlobalElements( &BaseGIDs[0] );

  for( int i = 0; i < Size_; ++i )
    for( int j = 0; j < BlockSize; ++j )
    {
      BlockGIDs[j] = BaseGIDs[j] + offset_*i;
      MPDEGIDs[i*BlockSize+j] = BaseGIDs[j] + offset_*i;
    }
  if (warpMPDE_)
  {
    // tscoffe/tmei 08/02/05
    omegaGID_ = offset_*Size_; // next unique GID
    phiGID_ = omegaGID_+1;
    MPDEGIDs[Size_*BlockSize] = omegaGID_;
    MPDEGIDs[Size_*BlockSize+1] = phiGID_;

    // tscoffe/tmei 08/02/05:  Added two to BlockSize*Size_
    //                         This extra size in the Map is not handled in the
    //                         createVector & createMatrix
    epetraMPDEMap_ = rcp(new Epetra_Map( BlockSize*Size_+2,
                                         BlockSize*Size_+2,
                                         &MPDEGIDs[0],
                                         BaseIndex,
                                         Comm ));
  }
  else
  {
    epetraMPDEMap_ = rcp(new Epetra_Map( BlockSize*Size_,
                                         BlockSize*Size_,
                                         &MPDEGIDs[0],
                                         BaseIndex,
                                         Comm ));
  }

  // Create N_PDS_ParMap for the state variables
  MPDEMap_ = rcp(new N_PDS_ParMap( &*epetraMPDEMap_, BaseMap_->pdsComm() ));
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateStateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateStateMaps( const RCP<N_PDS_ParMap>& BaseStateMap )
{
  //Save copy of base map
  BaseStateMap_ = BaseStateMap;

  const Epetra_Comm & Comm  = BaseStateMap_->petraMap()->Comm();

  //determine block offset
  int MaxGID = BaseStateMap->petraMap()->MaxAllGID();
  int BaseIndex = BaseStateMap->petraMap()->IndexBase();
  
  stateOffset_=1;
  while ( stateOffset_ <= MaxGID ) stateOffset_ *= 10;

  //Setup Block Maps and MPDE Map
  int BlockSize = BaseStateMap->petraMap()->NumGlobalElements();
  vector<int> BlockGIDs(BlockSize);
  vector<int> MPDEGIDs;
  MPDEGIDs.resize(Size_*BlockSize);
  vector<int> BaseGIDs(BlockSize);

  BaseStateMap->petraMap()->MyGlobalElements( &BaseGIDs[0] );

  for( int i = 0; i < Size_; ++i )
    for( int j = 0; j < BlockSize; ++j )
    {
      BlockGIDs[j] = BaseGIDs[j] + stateOffset_*i;
      MPDEGIDs[i*BlockSize+j] = BaseGIDs[j] + stateOffset_*i;
    }

  // Create Epetra map for the state variables
  epetraMPDEStateMap_ = rcp(new Epetra_Map( BlockSize*Size_,
                                            BlockSize*Size_,
                                            &MPDEGIDs[0],
                                            BaseIndex,
                                            Comm ));

  // Create N_PDS_ParMap for the state variables
  MPDEStateMap_ = rcp(new N_PDS_ParMap( &*epetraMPDEStateMap_, BaseStateMap_->pdsComm() ));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateStoreMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateStoreMaps( const RCP<N_PDS_ParMap>& BaseStoreMap )
{
  //Save copy of base map
  BaseStoreMap_ = BaseStoreMap;

  const Epetra_Comm & Comm  = BaseStoreMap_->petraMap()->Comm();

  //determine block offset
  int MaxGID = BaseStoreMap->petraMap()->MaxAllGID();
  int BaseIndex = BaseStoreMap->petraMap()->IndexBase();
  
  storeOffset_=1;
  while ( storeOffset_ <= MaxGID ) storeOffset_ *= 10;

  //Setup Block Maps and MPDE Map
  int BlockSize = BaseStoreMap->petraMap()->NumGlobalElements();
  vector<int> BlockGIDs(BlockSize);
  vector<int> MPDEGIDs;
  MPDEGIDs.resize(Size_*BlockSize);
  vector<int> BaseGIDs(BlockSize);

  BaseStoreMap->petraMap()->MyGlobalElements( &BaseGIDs[0] );

  for( int i = 0; i < Size_; ++i )
    for( int j = 0; j < BlockSize; ++j )
    {
      BlockGIDs[j] = BaseGIDs[j] + storeOffset_*i;
      MPDEGIDs[i*BlockSize+j] = BaseGIDs[j] + storeOffset_*i;
    }

  // Create Epetra map for the store variables
  epetraMPDEStoreMap_ = rcp(new Epetra_Map( BlockSize*Size_,
                                            BlockSize*Size_,
                                            &MPDEGIDs[0],
                                            BaseIndex,
                                            Comm ));

  // Create N_PDS_ParMap for the store variables
  MPDEStoreMap_ = rcp(new N_PDS_ParMap( &*epetraMPDEStoreMap_, BaseStoreMap_->pdsComm() ));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateGraphs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateGraphs( const Epetra_CrsGraph & BasedQdxGraph,
                                     const Epetra_CrsGraph & BasedFdxGraph,
                                     const Epetra_CrsGraph & BaseFullGraph )
{
  if( Teuchos::is_null(BaseMap_) ) abort(); //Need to setup Maps first

  //Copies of base graphs
  BasedQdxGraph_ = rcp(new Epetra_CrsGraph( BasedQdxGraph ));
  BasedFdxGraph_ = rcp(new Epetra_CrsGraph( BasedFdxGraph ));
  BaseFullGraph_ = rcp(new Epetra_CrsGraph( BaseFullGraph ));

  int BlockSize = BaseMap_->petraMap()->NumGlobalElements();

  //Construct MPDE dQdX Graph
  //MPDEdQdxGraph_ = new Epetra_CrsGraph( Copy,
  //                                      dynamic_cast<Epetra_BlockMap&>(*MPDEMap_),
  //                                      0 );

  //Construct MPDE dFdX Graph
  MPDEdFdxGraph_ = rcp(new Epetra_CrsGraph( Copy,
                                            *(MPDEMap_->petraBlockMap()),
                                            0 ));
  
  int MaxIndices = BasedQdxGraph_->MaxNumIndices();
  vector<int> Indices(MaxIndices);
  int NumIndices;
  int BaseRow;
  int MPDERow;
  for( int i = 0; i < Size_; ++i )
  {
    for( int j = 0; j < BlockSize; ++j )
    {
      BaseRow = BaseMap_->petraMap()->GID(j);
      BasedQdxGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      for( int k = 0; k < NumIndices; ++k ) Indices[k] += offset_*i;
      //Diagonal Block
      MPDERow = BaseRow + offset_*i;
      //MPDEdQdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &Indices[0] );
      // Insert in F also, as its the default graph for the full Jacobian matrix
      MPDEdFdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &Indices[0] );
    }
  }
#ifdef Xyce_DEBUG_MPDE
  cout << "Q and F graphs before adding (phi,phi) entry:" << endl;
  cout << "MPDEdQdxGraph = [same as MPDEdFdxGraph]" << endl;
  if (mpdeMgr_->debugLevel > 0)
  {
    // MPDEdQdxGraph_->Print(cout);
    cout << "MPDEdFdxGraph = " << endl;
    MPDEdFdxGraph_->Print(cout);
  }
#endif // Xyce_DEBUG_MPDE
  if (warpMPDE_)
  {
    // Add phi equation terms:  \dot{phi(t_1)} = omega(t_1)
    MPDERow = phiGID_;
    NumIndices = 1;
    Indices[0] = phiGID_;
    //MPDEdQdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &Indices[0] );
    // Insert in F also, as its the default graph for the full Jacobian matrix
    MPDEdFdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &Indices[0] );
  }
#ifdef Xyce_DEBUG_MPDE
  cout << "Q and F graphs after adding (phi,phi) entry:" << endl;
  cout << "MPDEdQdxGraph = [same as MPDEdFdxGraph]" << endl;
  if (mpdeMgr_->debugLevel > 0)
  {
    // MPDEdQdxGraph_->Print(cout);
    cout << "MPDEdFdxGraph = " << endl;
    MPDEdFdxGraph_->Print(cout);
  }
#endif // Xyce_DEBUG_MPDE

  MaxIndices = BasedFdxGraph_->MaxNumIndices();
  Indices.resize(MaxIndices);
  vector<int> NewIndices(MaxIndices);
  int DiscStart = Disc_->Start();
  int DiscWidth = Disc_->Width();
  vector<int> Cols(DiscWidth);
  for( int i = 0; i < Size_; ++i )
  {
    for( int j = 0; j < DiscWidth; ++j )
    {
      Cols[j] = i + (j+DiscStart);
      if( Cols[j] < 0 ) Cols[j] += Size_;
      else if( Cols[j] > (Size_-1) ) Cols[j] -= Size_;
    }

    for( int j = 0; j < BlockSize; ++j )
    {
      BaseRow = BaseMap_->petraMap()->GID(j);
      BasedFdxGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );

      MPDERow = BaseRow + offset_*i;
      for( int k = 0; k < DiscWidth; ++k )
      {
        int Shift = Cols[k]*offset_;
        for( int kk = 0; kk < NumIndices; ++kk ) NewIndices[kk] = Indices[kk] + Shift;
        MPDEdFdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &NewIndices[0] );
      }
    }
  }
  if (warpMPDE_)
  {
    // tscoffe 01/15/07 This block adds dependence on omega in the dFdx matrix everywhere that q is nonzero.
    // dqdt1 + omega dqdt2 + f + b
    // q = dqdt1   f = omega dqdt2 + f   b = b
    for( int i = 0; i < Size_; ++i )
    {
      for( int j = 0; j < BlockSize; ++j )
      {
        BaseRow = BaseMap_->petraMap()->GID(j);
        MPDERow = BaseRow + offset_*i;
        //BasedQdxGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
        //if ( NumIndices > 0 )
        //{
        NumIndices = 1;
        NewIndices[0] = omegaGID_; 
        MPDEdFdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &NewIndices[0] );
        //}
      }
    }
    Teuchos::RCP<vector<int> > phaseGraph = warpMPDEPhasePtr_->getPhaseGraph();
    MPDERow = omegaGID_;
    NumIndices = phaseGraph->size();
    MPDEdFdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &((*phaseGraph)[0]) );
    // Enter graph details for phi:  \dot{phi(t_1)} = omega(t_1)
    MPDERow = phiGID_;
    NumIndices = 1;
    NewIndices.clear();
    NewIndices.push_back(omegaGID_);
    MPDEdFdxGraph_->InsertGlobalIndices( MPDERow, NumIndices, &NewIndices[0] );
  }
  MPDEdFdxGraph_->FillComplete();

  //Construct MPDE dQdx Graph
  MPDEdQdxGraph_ = rcp(new Epetra_CrsGraph( *MPDEdFdxGraph_ ));
  MPDEdQdxGraph_->FillComplete();

  //Construct MPDE Full Graph
  MPDEFullGraph_ = rcp(new Epetra_CrsGraph( *MPDEdFdxGraph_ ));
  MPDEFullGraph_->FillComplete();

#ifdef Xyce_DEBUG_MPDE
  if (mpdeMgr_->debugLevel > 0)
  {  
    cout << "Final MPDEdQdxGraph = " << endl;
    MPDEdQdxGraph_->Print(cout);
    cout << "Final MPDEdFdxGraph = " << endl;
    MPDEdFdxGraph_->Print(cout);
  }
#endif // Xyce_DEBUG_MPDE

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::getSolutionMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
Teuchos::RCP<const Epetra_Map> N_MPDE_Builder::getSolutionMap() const
{
  return(epetraMPDEMap_);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::getStateMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
Teuchos::RCP<const Epetra_Map> N_MPDE_Builder::getStateMap() const
{
  return(epetraMPDEStateMap_);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::getStoreMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
Teuchos::RCP<const Epetra_Map> N_MPDE_Builder::getStoreMap() const
{
  return(epetraMPDEStoreMap_);
}

