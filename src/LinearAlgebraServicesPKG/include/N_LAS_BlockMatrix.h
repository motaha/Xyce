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
// Revision Number: $Revision: 1.14.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:44 $
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
  N_LAS_BlockMatrix( int Size,
                     const vector< vector<int> > & BlockColumns,
                     const Epetra_CrsGraph & Graph,
                     const Epetra_CrsGraph & BaseGraph,
                     int AugmentCount = 0 );

  //Destructor
  ~N_LAS_BlockMatrix();

  //Block Access
  N_LAS_Matrix & block( int row, int col );

  int blockSize()
  { return BlockSize_; }
  
  int numBlockRows()
  { return NumBlockRows_; }

  void replaceAugmentedColumn(int augmentedColumn, const N_LAS_BlockVector & vec);
  
  void printPetraObject() const;

 private:

  const int BlockSize_;
  const int NumBlockRows_;
  const int AugmentCount_;
  vector<int> AugmentGIDs_;
  
  const vector< vector<int> > Cols_;
  vector< vector<N_LAS_Matrix*> > Blocks_;

};

// Nonmember constructor
RCP<N_LAS_BlockMatrix> N_LAS_blockMatrix( 
    const vector< vector<int> > & BlockColumns,
    Epetra_Map & BlockMap,
    const Epetra_CrsGraph & BaseGraph,
    int AugmentCount = 0 
    );



#endif // Xyce_N_LAS_BlockMatrix_h
