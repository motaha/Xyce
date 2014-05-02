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
// Filename       : $RCSfile: N_LAS_BlockMatrix.h,v $
//
// Purpose        : Specification file for the Abstract interface to sparse
//                  block matrix type.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/12/04
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

#ifndef Xyce_N_LAS_BlockMatrix_h
#define Xyce_N_LAS_BlockMatrix_h

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

#include <N_LAS_Matrix.h>

// ----------  Other Includes   ----------
#include <Teuchos_RCP.hpp>
using Teuchos::RCP;

// ----------  Fwd Declares  -------------

class N_LAS_BlockVector;
class Epetra_CrsGraph;
class Epetra_Map;

//-----------------------------------------------------------------------------
// Class         : N_LAS_BlockMatrix
// Purpose       : Abstract interface to sparse block matrix type.
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
class N_LAS_BlockMatrix : public N_LAS_Matrix
{
 public:

  //Constructors
  N_LAS_BlockMatrix( int size,
                     int offset,
                     const std::vector< std::vector<int> > & blockColumns,
                     const Epetra_CrsGraph & globalGraph,
                     const Epetra_CrsGraph & subBlockGraph,
                     int augmentCount = 0 );

  //Destructor
  ~N_LAS_BlockMatrix() {}

  //Block Access
  N_LAS_Matrix & block( int row, int col );

  int blockSize()
  { return blockSize_; }
  
  int numBlockRows()
  { return numBlockRows_; }

  // Put function for the block sparse-matrix.
  void put(double s);

  // Replace the entries of an augmented row using the row GID.
  void replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices); 

  void replaceAugmentedColumn(int augmentedColumn, const N_LAS_BlockVector & vec);

  // Assemble global matrix with blocks
  // NOTE:  The global matrix is not always a view of the local matrix, so this function ensures
  // that the values are sync'ed up.  Call this before using the global matrix for computations.
  void assembleGlobalMatrix();
 
  void fillComplete();
 
  void printPetraObject(std::ostream &os) const;

 private:

  bool blocksViewGlobalMat_;
  const int blockSize_;
  const int offset_;
  const int numBlockRows_;
  const int augmentCount_;

  std::vector<int> augmentGIDs_, baseNumCols_, baseIndices_;
  const std::vector< std::vector<int> > cols_;
  std::vector< std::vector<Teuchos::RCP<N_LAS_Matrix> > > blocks_;
};

#endif // Xyce_N_LAS_BlockMatrix_h
