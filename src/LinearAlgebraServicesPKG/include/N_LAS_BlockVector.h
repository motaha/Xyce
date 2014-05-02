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
// Revision Number: $Revision: 1.23 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
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
  N_LAS_BlockVector( int numBlocks,
                     const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                     const Teuchos::RCP<N_PDS_ParMap> & subBlockMap,
                     int augmentRows = 0 );

  // Constructor that uses the block size to divide up the number of elements on
  // each processor into vectors whose values are all "owned" by one processor.
  // NOTE:  This constructor is handy for frequency-domain representations of time-domain vectors.
  N_LAS_BlockVector( int blockSize,
                     const Teuchos::RCP<N_PDS_ParMap> & globalMap );
  
  // Constructors to map to Petra constructors.
  N_LAS_BlockVector( int numBlocks,
                     const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                     const Teuchos::RCP<N_PDS_ParMap> & subBlockMap,
                     const Teuchos::RCP<N_PDS_ParMap> & osubBlockMap,
                     int augmentRows = 0 );

  //Copy constructor
  N_LAS_BlockVector( const N_LAS_BlockVector & right );

  //Copy constructor
  N_LAS_BlockVector( const N_LAS_Vector & right, const Teuchos::RCP<N_PDS_ParMap> & subBlockMap, int numBlocks );

  //Copy constructor
  //NOTE:  This constructor assumes that the N_LAS_Vector is divided up into blockSize subvectors,
  //       whose values are solely owned by one of the processors.
  N_LAS_BlockVector( const N_LAS_Vector & right, int blockSize );

  //View constructor
  N_LAS_BlockVector( Epetra_Vector * vector, const Teuchos::RCP<N_PDS_ParMap> & subBlockMap, int numBlocks, bool isOwned );

  // Destructor
  virtual ~N_LAS_BlockVector() {};

  // Block accessors
  N_LAS_Vector & block( int Loc ) const
  { return *blocks_[Loc]; }

  int blockSize() const
  { return globalBlockSize_; }

  int blockCount() const
  { return numBlocks_; }

  int startBlock() const
  { return startBlock_; }

  int endBlock() const
  { return endBlock_; }

  double * augmentStart()
  { return augmentLoc_; }
 /* 
  // Fill vector with constant value.
  void putScalar(const double scalar);

  // Add to vector with constant value.
  void addScalar(const double scalar);
*/
  // Assemble global vector with blocks
  // NOTE:  The global vector is not always a view of the local vector, so this function ensures
  // that the values are sync'ed up.  Call this before using the global vector for computations.
  void assembleGlobalVector();

  void printPetraObject(std::ostream &os) const;

 private:

  bool blocksViewGlobalVec_;
  const int globalBlockSize_;
  const int localBlockSize_;
  const int overlapBlockSize_;
  const int numBlocks_;
  const int augmentCount_;

  // In frequency domain, whole blocks may be owned by one processor.
  // NOTE:  This assumes they are contiguous.  By default these routines
  //        will return 0 and numBlocks_ (which is sane for the time domain specs).
  int startBlock_, endBlock_;

  double * augmentLoc_;

  Teuchos::RCP<N_PDS_ParMap> newBlockMap_, newoBlockMap_;

  std::vector<Teuchos::RCP<N_LAS_Vector> > blocks_;

};

#endif

