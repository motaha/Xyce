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
// Filename       : $RCSfile: N_LAS_BlockMatrix.C,v $
//
// Purpose        : Implementation  for Block Matrix
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/13/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.33 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>

#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>

// ----------   Other Includes -----------

#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_TestForException.hpp>

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::N_LAS_BlockMatrix
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockMatrix::N_LAS_BlockMatrix( int size,
                                      int offset,
                                      const std::vector< std::vector<int> > & blockColumns,
                                      const Epetra_CrsGraph & globalGraph,
                                      const Epetra_CrsGraph & subBlockGraph,
                                      int augmentCount )
: N_LAS_Matrix( &(const_cast<Epetra_CrsGraph&>(globalGraph)),
                &(const_cast<Epetra_CrsGraph&>(globalGraph)) ),
  blocksViewGlobalMat_(true),
  blockSize_( subBlockGraph.NumMyRows() ),
  offset_( offset ),
  numBlockRows_( size ),
  augmentCount_( augmentCount ),
  cols_( blockColumns ),
  blocks_( size )
{
  // Individual blocks cannot be a view of the global matrix because of graph ordering and communication
  if ( globalGraph.Comm().NumProc() > 1 )
  {
    blocksViewGlobalMat_ = false;

    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      blocks_[i].resize( numCols );
      for( int j = 0; j < numCols; ++j )
      {
        Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( Copy, subBlockGraph );
        blocks_[i][j] = Teuchos::rcp( new N_LAS_Matrix( bMat ) );
      }
    }

    // Get the local indices for the sub block so assembling is easier
    baseNumCols_.resize( blockSize_ );
    baseIndices_.resize( subBlockGraph.NumMyNonzeros() );
    int ptr = 0;
    for( int i = 0; i < blockSize_; ++i )
    {
      subBlockGraph.ExtractMyRowCopy( i, subBlockGraph.NumMyNonzeros()-ptr, baseNumCols_[i], &baseIndices_[ptr] );
      ptr += baseNumCols_[i];
    }
  }
  else
  {
    std::vector<int> baseNumCols( blockSize_ );
    std::vector< int* > baseIndices( blockSize_ );
    for( int i = 0; i < blockSize_; ++i )
      subBlockGraph.ExtractMyRowView( i, baseNumCols[i], baseIndices[i] );

    std::vector< double* > Values( blockSize_ );
    int NumEntries;

    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      blocks_[i].resize( numCols );

      for( int j = 0; j < blockSize_; ++j )
        aDCRSMatrix_->ExtractMyRowView( j+blockSize_*i, NumEntries, Values[j] );

      for( int j = 0; j < numCols; ++j )
      {
        Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( View, subBlockGraph );

        for( int k = 0; k < blockSize_; ++k )
          bMat->InsertMyValues( k, baseNumCols[k], Values[k]+j*baseNumCols[k], baseIndices[k] );

        blocks_[i][j] = Teuchos::rcp( new N_LAS_Matrix( bMat ) );
      }
    }
  }

  // Generate the augmented GIDs list.
  if( augmentCount_ )
  {
    augmentGIDs_.resize(augmentCount_);
    int augStart = blockSize_ * size;
    for( int i = 0; i < augmentCount_; ++i )
    {
      augmentGIDs_[i] = globalGraph.RowMap().GID(augStart+i);
    }
  }

  // Communicate the augmented GIDs to all processors.
  // All other processors other than the one that owns the augmented GID will have -1.
  std::vector<int> tmpAugmentGIDs = augmentGIDs_;
  globalGraph.Comm().MaxAll( &tmpAugmentGIDs[0], &augmentGIDs_[0], augmentCount_ );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::block
// Purpose       : Block Accessor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_Matrix & N_LAS_BlockMatrix::block( int row, int col )
{
  for( int i = 0; i < (int)(cols_[row].size()); ++i )
    if( cols_[row][i] == col )
      return *blocks_[row][i];

  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
      "Error!  N_LAS_BlockMatrix::block("<<row<<","<<col<<"):  This block does not exist!"
      );
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::put
// Purpose       : Put function for the sparse-matrix.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/04/00
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::put( double s )
{
  aDCRSMatrix_->PutScalar(s);

  if (!blocksViewGlobalMat_)
  {
    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();
      for ( int j = 0; j < numCols; ++j )
      {
        blocks_[i][j]->put( s );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::fillComplete
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/29/03
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::fillComplete()
{
  // Call fillComplete on all the individual N_LAS_Matrix blocks, 
  // then assemble the global matrix.
  for( int i = 0; i < numBlockRows_; ++i )
  {
    int numCols = cols_[i].size();
    for ( int j = 0; j < numCols; ++j )
    {
      blocks_[i][j]->fillComplete();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::assembleGlobalMatrix
// Purpose       : Fills global N_LAS_Matrix with the values in each block.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 10/03/13
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::assembleGlobalMatrix()
{
  if (!blocksViewGlobalMat_)
  {
    for( int i = 0; i < numBlockRows_; ++i )
    {
      int numCols = cols_[i].size();

      for( int k = 0; k < blockSize_; ++k )
      {
        // Create memory for all the entries for one whole row.
        int length = numCols*baseNumCols_[k];
        std::vector<int> Indices( length );
        std::vector<double> Values( length );

        // For each block column extract the current values from the subblock and correct the column indices.
        // NOTE:  All extractions and insertions are done using global ids, it seems easier that way.
        int ptr = 0;
        for( int j = 0; j < numCols; ++j )
        {
          int numIndices = 0;
          int col = cols_[i][j];
          int globalRow = (blocks_[i][j])->epetraObj().Graph().GRID(k);
          (*blocks_[i][j]).getRowCopy(globalRow, length-ptr, numIndices, &Values[ptr], &Indices[ptr]);
        
          // Correct the column indices for this subblock. 
          for (int idx = 0; idx < numIndices; idx++)
          {
            Indices[ptr+idx] += offset_*col;
          }
          ptr += numIndices;
        }
       
        // Insert the values for all the global columns at the same time.
        aDCRSMatrix_->ReplaceGlobalValues(aDCRSMatrix_->Graph().GRID(k + i*blockSize_), length, &Values[0], &Indices[0]); 
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::printPetraObject
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::printPetraObject(std::ostream &os) const
{
  os << "N_LAS_BlockMatrix Object (Size=" << numBlockRows_ << ", View =" << blocksViewGlobalMat_ << ")" << std::endl;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  for( int i = 0; i < numBlockRows_; ++i )
  {
    int numCols = cols_[i].size();
    for( int j = 0; j < numCols; ++j )
    {
      os << "Block[" << i << "][" << cols_[i][j] << "]\n";
      blocks_[i][j]->printPetraObject(os);
    }
  }
  os << "Base Object\n";
  os << *aDCRSMatrix_;
  os << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::replaceAugmentedColumn
// Purpose       : Replace augmented column of matrix with block vector
// Special Notes : This places values directly in the base object.
//               : This should NOT be used to replace values that are within blocks.
// Scope         : Public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 08/12/05
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::replaceAugmentedColumn(int colGID, const N_LAS_BlockVector & vec)
{
  const Epetra_BlockMap & RowMap = const_cast<N_LAS_BlockVector&>(vec).epetraObj().Map();
  int NumRows = RowMap.NumMyElements();

  int lCol = aDCRSMatrix_->Graph().LCID(colGID);
  for( int i = 0; i < NumRows; ++i )
  {
    int Row = i;
    double Val = vec[i];
    aDCRSMatrix_->ReplaceMyValues( Row, 1, &Val, &lCol );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::replaceAugmentedRow
// Purpose       : Replace augmented row values
// Special Notes : This places values directly in the base object.
//               : This should NOT be used to replace values that are within blocks.
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 11/12/13
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices)
{
  if ( aDCRSMatrix_->Graph().LRID( rowGID ) >= 0 )
  {
    aDCRSMatrix_->ReplaceGlobalValues( rowGID, length, coeffs, colIndices );
  }
}

