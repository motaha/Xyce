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
// Filename      : $RCSfile: N_LAS_HBBuilder.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.14 $
// Revision Date  : $Date: 2014/02/24 23:49:23 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_HBBuilder.h>
#include <N_MPDE_Discretization.h>

#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>

#include <N_PDS_ParMap.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_fwd.h>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Utils.hpp>

using Teuchos::rcp;
using Teuchos::RCP;

  // Default Constructor
//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::N_LAS_HBBuilder
// Purpose       : Default constructor 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_HBBuilder::N_LAS_HBBuilder( const int Size,
                            Teuchos::RCP<const N_MPDE_Discretization> Disc )
: numHarmonics_(Size),
  numSolVariables_(0),
  numStateVariables_(0),
  numStoreVariables_(0),
  blockPattern_(Size),
  Disc_(Disc),
  offset_(0),
  stateOffset_(0),
  storeOffset_(0)
{
  // Create the block graph pattern.
  int Start = Disc_->Start();
  int Width = Disc_->Width();

  for( int i = 0; i < numHarmonics_; ++i )
  {
    blockPattern_[i].resize(Width);
    for( int j = 0; j < Width; ++j )
    {
      int Loc = i+(j+Start);
      if( Loc < 0 )           Loc += numHarmonics_;
      else if( Loc >= numHarmonics_ ) Loc -= numHarmonics_;
      blockPattern_[i][j] = Loc;
    }
    std::sort(blockPattern_[i].begin(),blockPattern_[i].end());
  } 
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_HBBuilder::createVector( double initialValue ) const
{
  RCP<N_LAS_Vector> vector =
    createExpandedRealFormTransposeBlockVector(); 
  vector.release(); // Release ownership of the object.
  return(&*vector);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createTimeDomainBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createTimeDomainBlockVector( ) const
{
  RCP<N_LAS_BlockVector> vector = rcp(
        new N_LAS_BlockVector( numHarmonics_, HBMap_, BaseMap_ ) 
        );
  return(vector);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createTimeDomainStateBlockVector( ) 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createTimeDomainStateBlockVector( ) const
{
  RCP<N_LAS_BlockVector> vector = rcp(
        new N_LAS_BlockVector( numHarmonics_, HBStateMap_, BaseStateMap_ ) 
        );
  return(vector);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createTimeDomainStoreBlockVector( ) 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createTimeDomainStoreBlockVector( ) const
{
  RCP<N_LAS_BlockVector> vector = rcp(
        new N_LAS_BlockVector( numHarmonics_, HBStoreMap_, BaseStoreMap_ ) 
        );
  return(vector);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createExpandedRealFormBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Rich Schiek 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createExpandedRealFormBlockVector() const
{
  RCP<N_LAS_BlockVector> vec = rcp(
        new N_LAS_BlockVector( 2*numHarmonics_, HBExpandedRealFormBVMap_, BaseMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createExpandedRealFormTransposeBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Rich Schiek 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createExpandedRealFormTransposeBlockVector() const
{
  RCP<N_LAS_BlockVector> vec = rcp(
        new N_LAS_BlockVector( 2*numHarmonics_, HBExpandedRealFormBVMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createExpandedRealFormTransposeStateBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Rich Schiek 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createExpandedRealFormTransposeStateBlockVector() const
{
  RCP<N_LAS_BlockVector> vec = rcp(
        new N_LAS_BlockVector( 2*numHarmonics_, HBExpandedRealFormStateBVMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createExpandedRealFormTransposeStoreBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockVector> N_LAS_HBBuilder::createExpandedRealFormTransposeStoreBlockVector() const
{
  RCP<N_LAS_BlockVector> vec = rcp(
        new N_LAS_BlockVector( 2*numHarmonics_, HBExpandedRealFormStoreBVMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_HBBuilder::createStateVector( double initialValue ) const
{
  return dynamic_cast<N_LAS_Vector*>(
        new N_LAS_BlockVector( numHarmonics_, HBStateMap_, BaseStateMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
N_LAS_Vector * N_LAS_HBBuilder::createStoreVector( double initialValue ) const
{
  return dynamic_cast<N_LAS_Vector*>(
        new N_LAS_BlockVector( numHarmonics_, HBStoreMap_, BaseStoreMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createDAEdQdxMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_LAS_HBBuilder::createDAEdQdxMatrix( double initialValue ) const
{
  
  /*  NOTE:  The dQdx matrix is created with the same graph as the dFdx so
             that the construction of the Jacobian is faster */
  return createDAEdFdxMatrix( initialValue );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createDAEdFdxMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_LAS_HBBuilder::createDAEdFdxMatrix( double initialValue ) const
{
  return dynamic_cast<N_LAS_Matrix*>(
        new N_LAS_BlockMatrix( numHarmonics_,
                               offset_,
                               blockPattern_,
                               *HBFullGraph_,
                               *BaseFullGraph_ ) );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::createDAEFullMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_LAS_Matrix * N_LAS_HBBuilder::createDAEFullMatrix( double initialValue ) const
{
  return createDAEdFdxMatrix( initialValue );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_LAS_HBBuilder::generateMaps( const RCP<N_PDS_ParMap>& BaseMap, 
                                    const RCP<N_PDS_ParMap>& oBaseMap )
{
  //Save copy of base map
  BaseMap_ = BaseMap;
  //Xyce::dout() << "BaseMap" << std::endl;
  //BaseMap->petraMap()->Print(Xyce::dout());
  oBaseMap_ = oBaseMap;
  //Xyce::dout() << "oBaseMap" << std::endl;
  //oBaseMap->petraMap()->Print(Xyce::dout());

  //Determine block offset
  offset_ = generateOffset( *BaseMap );

  //Check to see if a graph can be generated using 32-bit integers for the 
  //specified number of harmonics.
  int MaxGID = BaseMap->maxGlobalEntity();
  long int hbMaxGID = MaxGID + offset_*(numHarmonics_-1);
  if ( hbMaxGID > Teuchos::OrdinalTraits<int>::max() )
  {
      int allowableHarmonics = (Teuchos::OrdinalTraits<int>::max() - MaxGID + offset_)/offset_ - 1;
      std::string msg = "N_LAS_HBBuilder::generateMaps():  Harmonic Balance map cannot be constructed, reduce number of harmonics to " 
                      + Teuchos::Utils::toString( allowableHarmonics ) + ".";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_ERROR_0, msg);
  }                          

  // Use the block linear system helper to create the block parallel maps
  std::vector<RCP<N_PDS_ParMap> > blockMaps = createBlockParMaps(numHarmonics_, *BaseMap, *oBaseMap);

  HBMap_ = blockMaps[0];
  oHBMap_ = blockMaps[1];
  
  //Xyce::dout() << "HBMap_" << std::endl;
  //Xyce::dout() << *HBMap_ << std::endl; 
  //Xyce::dout() << "oHBMap_" << std::endl;
  //Xyce::dout() << *oHBMap_ << std::endl;
  //Save copy of base map

  // Helpful names for various sizes (subtract 1 for ground node):
  numSolVariables_ = oBaseMap_->numLocalEntities()-1;

  // Expanded Real Form for Block Vector map:
  HBExpandedRealFormBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseMap );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::generateStateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool N_LAS_HBBuilder::generateStateMaps( const RCP<N_PDS_ParMap>& BaseStateMap )
{
  //Save copy of base map
  BaseStateMap_ = BaseStateMap;

  //determine block offset
  stateOffset_ = generateOffset( *BaseStateMap );

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  HBStateMap_ = createBlockParMap(numHarmonics_, *BaseStateMap);

  // Helpful names for various sizes:
  numStateVariables_ = BaseStateMap_->numGlobalEntities();

  // Expanded Real Form for Block Vector map:
  HBExpandedRealFormStateBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseStateMap );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::generateStoreMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_LAS_HBBuilder::generateStoreMaps( const RCP<N_PDS_ParMap>& BaseStoreMap )
{
  //Save copy of base map
  BaseStoreMap_ = BaseStoreMap;

  //determine block offset
  storeOffset_ = generateOffset( *BaseStoreMap );

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  HBStoreMap_ = createBlockParMap(numHarmonics_, *BaseStoreMap);

  // Helpful names for various sizes:
  numStoreVariables_ = BaseStoreMap_->numGlobalEntities();

  // Expanded Real Form for Block Vector map:
  HBExpandedRealFormStoreBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseStoreMap );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::generateGraphs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_LAS_HBBuilder::generateGraphs( const Epetra_CrsGraph & BaseFullGraph )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("N_LAS_HBBuilder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of base graphs
  BaseFullGraph_ = rcp(new Epetra_CrsGraph( BaseFullGraph ));

  int BlockSize = BaseMap_->numLocalEntities();

  //Check to see if a graph can be generated using 32-bit integers for the
  //specified number of harmonics.
  int MaxGID = BaseMap_->maxGlobalEntity();
  long int hbMaxGID = MaxGID + offset_*(numHarmonics_-1);
  if ( hbMaxGID > Teuchos::OrdinalTraits<int>::max() )
  {
      int allowableHarmonics = (Teuchos::OrdinalTraits<int>::max() - MaxGID + offset_)/offset_ - 1;
      std::string msg = "N_LAS_HBBuilder::generateGraphs():  Harmonic Balance graph cannot be constructed, reduce number of harmonics to "
                      + Teuchos::Utils::toString( allowableHarmonics ) + ".";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_ERROR_0, msg);
  }

  //Construct HB Graph [All graphs are the same, so only one needs to be made]
  HBFullGraph_ = rcp(new Epetra_CrsGraph( Copy,
                                          dynamic_cast<Epetra_BlockMap&>(*(HBMap_->petraMap())),
                                          0 ));
  
  int MaxIndices = BaseFullGraph_->MaxNumIndices();
  std::vector<int> Indices(MaxIndices);
  int NumIndices=0, BaseRow=0, HBRow=0;

  for( int i = 0; i < numHarmonics_; ++i )
  {
    for( int j = 0; j < BlockSize; ++j )
    {
      BaseRow = BaseMap_->petraMap()->GID(j);
      BaseFullGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      for( int k = 0; k < NumIndices; ++k ) Indices[k] += offset_*i;
      //Diagonal Block
      HBRow = BaseRow + offset_*i;
      // Insert in F, as its the default graph for the full Jacobian matrix
      HBFullGraph_->InsertGlobalIndices( HBRow, NumIndices, &Indices[0] );
    }
  }
  
#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "N_LAS_HBBuilder::generateGraphs():  HB full graph after inserting diagonal blocks:" << std::endl;
  HBFullGraph_->Print(std::cout);
#endif // Xyce_DEBUG_HB
 
  MaxIndices = BaseFullGraph_->MaxNumIndices();
  Indices.resize(MaxIndices);
  int Shift=0, Index=0;
  int DiscWidth = Disc_->Width();
  std::vector<int> NewIndices(MaxIndices*DiscWidth);
  
  for( int j = 0; j < BlockSize; ++j )
  {
    // Extract the base entries from the base row.
    BaseRow = BaseMap_->petraMap()->GID(j);
    BaseFullGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
  
    for( int i = 0; i < numHarmonics_; ++i )
    {
      // For this harmonic, which row will be inserted.
      HBRow = BaseRow + offset_*i;
      
      // Find all entries from a row before inserting it.
      for( int k = 0; k < DiscWidth; ++k )
      {
        // Find which block column to start at.
        Shift = blockPattern_[i][k]*offset_;  // Actual column index.
        Index = k*NumIndices;  // Pointer to next block of column indices.
        for( int kk = 0; kk < NumIndices; ++kk ) NewIndices[Index+kk] = Indices[kk] + Shift;
      }

      // Insert entire row for all blocks.
      HBFullGraph_->InsertGlobalIndices( HBRow, DiscWidth*NumIndices, &NewIndices[0] );
    }
  }
  HBFullGraph_->FillComplete();

#ifdef Xyce_DEBUG_HB
  Xyce::dout() << "N_LAS_HBBuilder::generateGraphs():  Final HB Graph = " << std::endl;
  HBFullGraph_->Print(std::cout);
#endif // Xyce_DEBUG_HB

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::getSolutionMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<const Epetra_Map> N_LAS_HBBuilder::getSolutionMap() const
{
    return(rcp(HBExpandedRealFormBVMap_->petraMap(),false));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::getStateMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<const Epetra_Map> N_LAS_HBBuilder::getStateMap() const
{
  return(rcp(HBStateMap_->petraMap(),false));
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_HBBuilder::getStoreMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<const Epetra_Map> N_LAS_HBBuilder::getStoreMap() const
{
  return(rcp(HBStoreMap_->petraMap(),false));
}


