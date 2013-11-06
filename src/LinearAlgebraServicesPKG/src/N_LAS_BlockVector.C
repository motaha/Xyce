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
// Revision Date  : $Date: 2013/10/03 17:23:45 $
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
N_LAS_BlockVector::N_LAS_BlockVector( int Size,
                                      Epetra_Map & Map,
                                      Epetra_Map & BlockMap,
                                      int AugmentRows )
: N_LAS_Vector( new Epetra_Vector(dynamic_cast<Epetra_BlockMap&>(Map)) ),
  globalBlockSize_(BlockMap.NumGlobalElements()),
  localBlockSize_(BlockMap.NumMyElements()),
  overlapBlockSize_(BlockMap.NumMyElements()),
  BlockCount_(Size),
  AugmentCount_(AugmentRows),
  AugmentLoc_(0),
  offset_(1),
  BlockMap_(BlockMap),
  oBlockMap_(BlockMap),
  Blocks_(Size)
{
  //Setup Views of blocks using Block Map
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;

  for( int i = 0; i < Size; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(BlockMap_), Loc ) ) );
  }

  if( AugmentCount_ ) AugmentLoc_ = Ptrs[0] + overlapBlockSize_*Size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::N_LAS_BlockVector
// Purpose       : constructor
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( int Size,
                                      const Teuchos::RCP<N_PDS_ParMap> & Map,
                                      const Teuchos::RCP<N_PDS_ParMap> & BlockMap,
                                      int AugmentRows )
: N_LAS_Vector( *Map ),
  globalBlockSize_(BlockMap->petraMap()->NumGlobalElements()),
  localBlockSize_(BlockMap->petraMap()->NumMyElements()),
  overlapBlockSize_(BlockMap->petraMap()->NumMyElements()),
  BlockCount_(Size),
  AugmentCount_(AugmentRows),
  AugmentLoc_(0),
  offset_(1),
  newMap_(Map),
  newBlockMap_(BlockMap),
  BlockMap_(*(BlockMap->petraMap())),
  oBlockMap_(*(BlockMap->petraMap())),
  Blocks_(Size)
{ 
  //Setup Views of blocks using Block Map
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;
  
  for( int i = 0; i < Size; ++i )
  { 
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(BlockMap_), Loc ) ) );
  }
  
  if( AugmentCount_ ) AugmentLoc_ = Ptrs[0] + overlapBlockSize_*Size;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector::N_LAS_BlockVector
// Purpose       : constructor
// Special Notes : 
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( int Size,
                                      Epetra_Map & Map,
                                      Epetra_Map & BlockMap,
                                      Epetra_Map & oBlockMap,
                                      int AugmentRows )
: N_LAS_Vector( new Epetra_Vector(dynamic_cast<Epetra_BlockMap&>(Map)) ),
  globalBlockSize_(BlockMap.NumGlobalElements()),
  localBlockSize_(BlockMap.NumMyElements()),
  overlapBlockSize_(oBlockMap.NumMyElements()),
  BlockCount_(Size),
  AugmentCount_(AugmentRows),
  AugmentLoc_(0),
  offset_(1),
  BlockMap_(BlockMap),
  oBlockMap_(oBlockMap),
  Blocks_(Size)
{
  //Setup Views of blocks using Block Map
  double ** Ptrs;
  oMultiVector_->ExtractView( &Ptrs );
  double * Loc;

  for( int i = 0; i < Size; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(oBlockMap_), Loc ), BlockMap ) );
  }

  if( AugmentCount_ ) AugmentLoc_ = Ptrs[0] + overlapBlockSize_*Size;
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
// local block map and overlap map (BlockMap / oBlockMap)
N_LAS_BlockVector::N_LAS_BlockVector( int Size,
                                      const Teuchos::RCP<N_PDS_ParMap> & Map,
                                      const Teuchos::RCP<N_PDS_ParMap> & oMap,
                                      const Teuchos::RCP<N_PDS_ParMap> & BlockMap,
                                      const Teuchos::RCP<N_PDS_ParMap> & oBlockMap )
: N_LAS_Vector( *Map, *oMap ),
  globalBlockSize_(BlockMap->numGlobalEntities()),
  localBlockSize_(BlockMap->numLocalEntities()),
  overlapBlockSize_(oBlockMap->numLocalEntities()),
  BlockCount_(Size),
  AugmentCount_(0),
  AugmentLoc_(0),
  offset_(1),
  newMap_(Map),
  newoMap_(oMap),
  newBlockMap_(BlockMap),
  newoBlockMap_(oBlockMap),
  BlockMap_(*(BlockMap->petraMap())),
  oBlockMap_(*(oBlockMap->petraMap())),
  Blocks_(Size)
{
  // Compute the offset needed for global indexing.
  int maxGID = Map->maxGlobalEntity();
  int numProcs = Map->pdsComm()->numProc();
  while( offset_ <= maxGID ) offset_ *= 10;

  //Setup Views of blocks using Block Map
  double ** Ptrs;
  oMultiVector_->ExtractView( &Ptrs );
  double * Loc;
  
  for( int i = 0; i < Size; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(oBlockMap_), Loc ), *(newBlockMap_->petraMap()) ) );
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
  globalBlockSize_( rhs.globalBlockSize_ ),
  localBlockSize_( rhs.localBlockSize_ ),
  overlapBlockSize_( rhs.overlapBlockSize_ ),
  BlockCount_( rhs.BlockCount_ ),
  AugmentCount_( rhs.AugmentCount_ ),
  AugmentLoc_( rhs.AugmentLoc_ ),
  offset_( rhs.offset_ ),
  BlockMap_(rhs.BlockMap_),
  oBlockMap_(rhs.oBlockMap_),
  Blocks_( rhs.Blocks_.size() )
{
  //Setup Views of blocks using Block Map
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;

  for( int i = 0; i < BlockCount_; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] =  Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(BlockMap_), Loc ) ) );
  }

  if( AugmentCount_ ) AugmentLoc_ = Ptrs[0] + overlapBlockSize_*BlockCount_;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockVector:::N_LAS_BlockVector
// Purpose       : copy constructor
// Special Notes : 
// Scope         : Public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/10/08
//-----------------------------------------------------------------------------
N_LAS_BlockVector::N_LAS_BlockVector( const N_LAS_Vector & rhs, const Epetra_Map& subBlockMap, int numBlocks )
: N_LAS_Vector( rhs ),
  globalBlockSize_( subBlockMap.NumGlobalElements()),
  localBlockSize_( subBlockMap.NumMyElements()),
  overlapBlockSize_( subBlockMap.NumMyElements()),
  BlockCount_( numBlocks ),
  AugmentCount_( 0 ),
  AugmentLoc_(0),
  offset_(1),
  BlockMap_(subBlockMap),
  oBlockMap_(subBlockMap),
  Blocks_( numBlocks )
{
//  Epetra_MultiVector & mv = this->epetraObj();
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;
 
  for( int i = 0; i < BlockCount_; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] = Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(BlockMap_), Loc ) ) );
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
N_LAS_BlockVector::N_LAS_BlockVector( Epetra_Vector * rhs, const Epetra_Map& subBlockMap, int numBlocks, bool isOwned )
: N_LAS_Vector( rhs, isOwned ),
  globalBlockSize_( subBlockMap.NumGlobalElements()),
  localBlockSize_( subBlockMap.NumMyElements()),
  overlapBlockSize_( subBlockMap.NumMyElements()),
  BlockCount_( numBlocks ),
  AugmentCount_( 0 ),
  AugmentLoc_(0),
  offset_(1),
  BlockMap_(subBlockMap),
  oBlockMap_(subBlockMap),
  Blocks_( numBlocks )
{
//  Epetra_MultiVector & mv = this->epetraObj();
  double ** Ptrs;
  aMultiVector_->ExtractView( &Ptrs );
  double * Loc;
 
  for( int i = 0; i < BlockCount_; ++i )
  {
    Loc = Ptrs[0] + overlapBlockSize_*i;
    Blocks_[i] = Teuchos::rcp( new N_LAS_Vector( new Epetra_Vector( View, dynamic_cast<const Epetra_BlockMap&>(BlockMap_), Loc ) ) );
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
void N_LAS_BlockVector::printPetraObject() const
{
  cout << "N_LAS_BlockVector Object (Size=" << BlockCount_ << ")\n";
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < BlockCount_; ++i )
  {
    cout << "Block[" << i << "]\n";
    Blocks_[i]->printPetraObject();
  }
  cout << "Base Object\n";
  cout << *aMultiVector_;
  cout << *oMultiVector_;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}
