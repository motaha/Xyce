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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_LAS_BlockVector.C,v $
//
// Purpose        : Implementation file for Block Vector
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Computational Sciences
//
// Creation Date  : 3/13/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_BlockVector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <N_LAS_BlockSystemHelpers.h>

// ---------  Other Includes  -----------

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::N_LAS_BlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( int numBlocks,
                                      const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                                      const Teuchos::RCP<N_PDS_ParMap> & subBlockMap,
                                      int AugmentRows )
: N_LAS_Vector( *globalMap ),
  blocksViewGlobalVec_(true),
  globalBlockSize_(subBlockMap->numGlobalEntities()),
  localBlockSize_(subBlockMap->numLocalEntities()),
  overlapBlockSize_(subBlockMap->numLocalEntities()),
  numBlocks_(numBlocks),
  augmentCount_(AugmentRows),
  startBlock_(0),
  endBlock_(numBlocks),
  newBlockMap_(subBlockMap),
  newoBlockMap_(subBlockMap),
  blocks_(numBlocks)
{
  //Setup Views of blocks using Block Map
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;

  for( int i = 0; i < numBlocks; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(*newBlockMap_->petraMap()), Loc ) ) );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::N_LAS_BlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( int blockSize,
                                      const Teuchos::RCP<N_PDS_ParMap> & globalMap )
: N_LAS_Vector( *globalMap ),
  blocksViewGlobalVec_( true ),
  globalBlockSize_( blockSize ),
  localBlockSize_( blockSize ),
  overlapBlockSize_( blockSize ),
  numBlocks_( globalMap->numGlobalEntities() / blockSize ),
  augmentCount_( 0 ),
  startBlock_( 0 ),
  endBlock_( globalMap->numGlobalEntities() / blockSize ),
  newBlockMap_( Teuchos::rcp( new N_PDS_ParMap( blockSize, blockSize, globalMap->indexBase(), globalMap->pdsComm() ) ) ),
  newoBlockMap_( Teuchos::rcp( new N_PDS_ParMap( blockSize, blockSize, globalMap->indexBase(), globalMap->pdsComm() ) ) ),
  blocks_( globalMap->numGlobalEntities() / blockSize )
{
  // Determine where these blocks start and end in the grand scheme of things.
  startBlock_ = (int) std::floor( (double)(globalMap->petraMap()->MinMyGID() + 1) / (double)blockSize );
  endBlock_ = (int) std::floor( (double)(globalMap->petraMap()->MaxMyGID() + 1) / (double)blockSize );

  //std::cout << "startBlock_ = " << startBlock_ << ", endBlock_ = " << endBlock_ << std::endl;

  //Setup Views of blocks using Block Map
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc = 0;
  if (globalMap->numLocalEntities() > 0)
  {
    Loc = Ptrs[0];
  }

  for( int i = 0; i < numBlocks_; ++i )
  {
    int myBlockSize = 0;

    // Generate maps where all the entries of the block are owned by one processor.
    if ( (i >= startBlock_) && (i < endBlock_) )
    {
      myBlockSize = blockSize;
    }
    N_PDS_ParMap currBlockMap( blockSize, myBlockSize, globalMap->indexBase(), globalMap->pdsComm() );

    // Create a N_LAS_Vector that views all the block data that is local.
    blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(*(currBlockMap.petraMap())), Loc ) ) );

    if ( (i >= startBlock_) && (i < endBlock_) )
    {
      // Advance the pointer for the local data.
      Loc += blockSize;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::N_LAS_BlockVector
// Purpose       : constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
// Constructor that takes the global map and overlap map (Map / oMap) as well as the
// local block map and overlap map (subBlockMap / osubBlockMap)
N_LAS_BlockVector::N_LAS_BlockVector( int numBlocks,
                                      const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                                      const Teuchos::RCP<N_PDS_ParMap> & subBlockMap,
                                      const Teuchos::RCP<N_PDS_ParMap> & osubBlockMap,
                                      int AugmentRows )
: N_LAS_Vector(*globalMap),
  blocksViewGlobalVec_(false), 
  globalBlockSize_(subBlockMap->numGlobalEntities()),
  localBlockSize_(subBlockMap->numLocalEntities()),
  overlapBlockSize_(osubBlockMap->numLocalEntities()),
  numBlocks_(numBlocks),
  augmentCount_(AugmentRows),
  startBlock_(0),
  endBlock_(numBlocks),
  newBlockMap_(subBlockMap),
  newoBlockMap_(osubBlockMap),
  blocks_(numBlocks)
{
  for( int i = 0; i < numBlocks; ++i )
  {
    blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( *subBlockMap, *osubBlockMap ) );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::N_LAS_BlockVector
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( const N_LAS_BlockVector & rhs )
: N_LAS_Vector( dynamic_cast<const N_LAS_Vector&>(rhs) ),
  blocksViewGlobalVec_( rhs.blocksViewGlobalVec_ ),
  globalBlockSize_( rhs.globalBlockSize_ ),
  localBlockSize_( rhs.localBlockSize_ ),
  overlapBlockSize_( rhs.overlapBlockSize_ ),
  numBlocks_( rhs.numBlocks_ ),
  augmentCount_( rhs.augmentCount_ ),
  startBlock_( rhs.startBlock_ ),
  endBlock_( rhs.endBlock_ ),
  newBlockMap_( rhs.newBlockMap_ ),
  newoBlockMap_( rhs.newoBlockMap_ ),
  blocks_( rhs.blocks_.size() )
{
  if (blocksViewGlobalVec_)
  {
    // If the startBlock_ and endBlock_ cover every block in this vector than this is a time-domain representation
    // or serial simulation, in which case a frequency-domain distinction need not be made.
    if ((startBlock_ == 0) && (endBlock_ == numBlocks_))
    {
      int numBlocks = blocks_.size();

      // Setup Views of blocks using Block Map
      double ** Ptrs;
      aMultiVector_->ExtractView( &Ptrs );
      double * Loc;

      for( int i = 0; i < numBlocks; ++i )
      {
        Loc = Ptrs[0] + overlapBlockSize_*i;
        blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(*(newBlockMap_->petraMap())), Loc ) ) );
      }
    }
    else
    {
      // This is a frequency-domain representation of the block vector, so create views accordingly.
      int blockSize = globalBlockSize_;

      // Setup Views of blocks using Block Map
      double ** Ptrs;
      aMultiVector_->ExtractView( &Ptrs );
      double * Loc = Ptrs[0];

      for( int i = 0; i < numBlocks_; ++i )
      {
        // Create a N_LAS_Vector that views all the block data that is local.
        blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, ((rhs.blocks_[i])->epetraObj()).Map(), Loc ) ) );

        if ( (i >= startBlock_) && (i < endBlock_) )
        {
          // Advance the pointer for the local data.
          Loc += blockSize;
        }
      }
    }
  }
  else
  {
    for( int i = 0; i < numBlocks_; ++i )
    {
      blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( *(rhs.blocks_[i]) ) );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector:::N_LAS_BlockVector
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/10/08
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( const N_LAS_Vector & rhs, const Teuchos::RCP<N_PDS_ParMap> & subBlockMap, int numBlocks )
: N_LAS_Vector( rhs ),
  blocksViewGlobalVec_( true ), 
  globalBlockSize_( subBlockMap->numGlobalEntities() ),
  localBlockSize_( subBlockMap->numLocalEntities() ),
  overlapBlockSize_( subBlockMap->numLocalEntities() ),
  numBlocks_( numBlocks ),
  augmentCount_( 0 ),
  startBlock_( 0 ),
  endBlock_( numBlocks ),
  newBlockMap_( subBlockMap ),
  newoBlockMap_( subBlockMap ),
  blocks_( numBlocks )
{
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;

  for( int i = 0; i < numBlocks_; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    blocks_[i] = Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(*(newBlockMap_->petraMap())), Loc ) ) );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector:::N_LAS_BlockVector
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/10/08
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( const N_LAS_Vector & rhs, int blockSize )
: N_LAS_Vector( rhs ),
  blocksViewGlobalVec_( true ), 
  globalBlockSize_( blockSize ),
  localBlockSize_( blockSize ),
  overlapBlockSize_( blockSize ),
  numBlocks_( rhs.globalLength() / blockSize ),
  augmentCount_( 0 ),
  startBlock_( 0 ),
  endBlock_( rhs.globalLength() / blockSize ),
  blocks_( rhs.globalLength() / blockSize )
{
  // Create the new maps for each block that places all the entries of the block on one processor.
  N_LAS_Vector& rhs_nonconst = const_cast<N_LAS_Vector&>( rhs );
  newBlockMap_ = Teuchos::rcp( new N_PDS_ParMap( blockSize, blockSize, 
                                 ( aMultiVector_->Map() ).IndexBase(),
                                 rhs_nonconst.pdsComm() ) );
  newoBlockMap_ = Teuchos::rcp( new N_PDS_ParMap( blockSize, blockSize, 
                                 ( aMultiVector_->Map() ).IndexBase(),
                                 rhs_nonconst.pdsComm() ) );

  // Determine where these blocks start and end in the grand scheme of things.
  int minMyGID = (aMultiVector_->Map()).MinMyGID();
  int maxMyGID = (aMultiVector_->Map()).MaxMyGID();
  startBlock_ = (int) std::floor( (double)(minMyGID + 1) / (double)blockSize );
  endBlock_ = (int) std::floor( (double)(maxMyGID + 1) / (double)blockSize );

  //Setup Views of blocks using Block Map
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc = Ptrs[0];

  for( int i = 0; i < numBlocks_; ++i )
  {
    // Generate maps where all the entries of the block are owned by one processor.
    int myBlockSize = 0;

    // Generate maps where all the entries of the block are owned by one processor.
    if ( (i >= startBlock_) && (i < endBlock_) )
    {
      myBlockSize = blockSize;
    }
    N_PDS_ParMap currBlockMap( blockSize, myBlockSize, newBlockMap_->indexBase(), newBlockMap_->pdsComm() );

    // Create a N_LAS_Vector that views all the block data that is local.
    blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(*(currBlockMap.petraMap())), Loc ) ) );

    if ( (i >= startBlock_) && (i < endBlock_) )
    {
      // Advance the pointer for the local data.
      Loc += blockSize;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector:::N_LAS_BlockVector
// Purpose       : view constructor
// Special Notes :
// Scope         : Public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( Epetra_Vector * rhs, const Teuchos::RCP<N_PDS_ParMap> & subBlockMap, int numBlocks, bool isOwned )
: N_LAS_Vector( rhs, isOwned ),
  blocksViewGlobalVec_( true ), 
  globalBlockSize_( subBlockMap->numGlobalEntities() ),
  localBlockSize_( subBlockMap->numLocalEntities() ),
  overlapBlockSize_( subBlockMap->numLocalEntities() ),
  numBlocks_( numBlocks ),
  augmentCount_( 0 ),
  startBlock_( 0 ),
  endBlock_( numBlocks ),
  newBlockMap_( subBlockMap ),
  newoBlockMap_( subBlockMap ),
  blocks_( numBlocks )
{
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;

  for( int i = 0; i < numBlocks_; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    blocks_[i] = Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(*(subBlockMap->petraMap())), Loc ) ) );
  }
}
/*
//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::putScalar
// Purpose       : Fills N_LAS_BlockVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_BlockVector::putScalar(const double scalar)
{
  for( int i = 0; i < numBlocks_; ++i )
  {
    blocks_[i]->putScalar( scalar );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::addScalar
// Purpose       : Adds to N_LAS_BlockVector with the value "scalar".
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 01/18/01
//-----------------------------------------------------------------------------
void N_LAS_BlockVector::addScalar(const double scalar)
{
  for( int i = 0; i < numBlocks_; ++i )
  {
    blocks_[i]->addScalar( scalar );
  }
}
*/
//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::assembleGlobalVector
// Purpose       : Fills global N_LAS_MultiVector with the values in each block.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------
void N_LAS_BlockVector::assembleGlobalVector()
{
  if (!blocksViewGlobalVec_)
  {

  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector:::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void N_LAS_BlockVector::printPetraObject(std::ostream &os) const
{
  os << "N_LAS_BlockVector Object (Number of Blocks =" << numBlocks_ << ", View =" << blocksViewGlobalVec_ << ")" << std::endl;

  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlocks_; ++i )
  {
    if (i >= startBlock_ && i < endBlock_)
    {
      os << "Block[" << i << "]\n";
    }
    blocks_[i]->printPetraObject( os );
  }
  os << "Base Object\n";
  os << *aMultiVector_;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}
