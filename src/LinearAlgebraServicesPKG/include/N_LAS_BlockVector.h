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
// Filename       : $RCSfile: N_LAS_BlockVector.h,v $
//
// Purpose        : Block Vector access
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Computational Sciences
//
// Creation Date  : 3/12/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:44 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_BlockVector_h
#define Xyce_N_LAS_BlockVector_h

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

#include <N_LAS_Vector.h>

// ----------  Other Includes   ----------

#include <Teuchos_RCP.hpp>

// --------  Forward Declarations --------

class Epetra_Map;
class N_PDS_ParMap;

//-----------------------------------------------------------------------------
// Class         : N_LAS_BlockVector
// Purpose       : Provides an abstract interface for block vectors
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 3/12/04
//-----------------------------------------------------------------------------
class N_LAS_BlockVector : public N_LAS_Vector
{
 public:

  // Constructors to map to Petra constructors.
  N_LAS_BlockVector( int BlockCount,
                     Epetra_Map & Map,
                     Epetra_Map & BlockMap,
                     int augmentRows = 0  );

  N_LAS_BlockVector( int BlockCount,
                     const Teuchos::RCP<N_PDS_ParMap> & Map,
                     const Teuchos::RCP<N_PDS_ParMap> & BlockMap,
                     int augmentRows = 0  );

  // Constructors to map to Petra constructors.
  N_LAS_BlockVector( int BlockCount,
                     Epetra_Map & Map,
                     Epetra_Map & BlockMap,
                     Epetra_Map & oBlockMap,
                     int augmentRows = 0  );

  // Constructor that takes the global map and overlap map (Map / oMap) as well as the 
  // local block map and overlap map (BlockMap / oBlockMap)
  N_LAS_BlockVector( int Size,
                     const Teuchos::RCP<N_PDS_ParMap> & Map,
                     const Teuchos::RCP<N_PDS_ParMap> & oMap,
                     const Teuchos::RCP<N_PDS_ParMap> & BlockMap,
                     const Teuchos::RCP<N_PDS_ParMap> & oBlockMap );

  //Copy constructor
  N_LAS_BlockVector( const N_LAS_BlockVector & right );

  //Copy constructor
  N_LAS_BlockVector( const N_LAS_Vector & right, const Epetra_Map& subBlockMap, int numBlocks );

  //View constructor
  N_LAS_BlockVector( Epetra_Vector * vector, const Epetra_Map& subBlockMap, int numBlocks, bool isOwned );

  // Destructor
  virtual ~N_LAS_BlockVector() {};

  // Block accessors
  N_LAS_Vector & block( int Loc ) const
  { return *Blocks_[Loc]; }

  int blockSize() const
  { return globalBlockSize_; }

  int blockCount() const
  { return BlockCount_; }

  double * augmentStart()
  { return AugmentLoc_; }

  void printPetraObject() const;

 private:

  const int globalBlockSize_;
  const int localBlockSize_;
  const int overlapBlockSize_;
  const int BlockCount_;
  const int AugmentCount_;
  double * AugmentLoc_;
  int offset_;

  Teuchos::RCP<N_PDS_ParMap> newMap_, newoMap_;
  Teuchos::RCP<N_PDS_ParMap> newBlockMap_, newoBlockMap_;

  const Epetra_Map & BlockMap_;
  const Epetra_Map & oBlockMap_;

  vector<Teuchos::RCP<N_LAS_Vector> > Blocks_;

};

#endif

