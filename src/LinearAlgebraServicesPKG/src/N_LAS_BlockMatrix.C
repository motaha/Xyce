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
// Revision Number: $Revision: 1.22.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>

// ----------   Other Includes -----------

#include <Epetra_Map.h>
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
N_LAS_BlockMatrix::N_LAS_BlockMatrix( int Size,
                                      const vector< vector<int> > & BlockColumns,
                                      const Epetra_CrsGraph & Graph,
                                      const Epetra_CrsGraph & BaseGraph,
                                      int AugmentCount )
: N_LAS_Matrix( &(const_cast<Epetra_CrsGraph&>(Graph)),
                &(const_cast<Epetra_CrsGraph&>(Graph)) ),
  BlockSize_( BaseGraph.NumMyRows() ),
  NumBlockRows_( Size ),
  AugmentCount_( AugmentCount ),
  Cols_( BlockColumns ),
  Blocks_( Size )
{
  vector<int> BaseNumCols( BlockSize_ );
  vector< int* > BaseIndices( BlockSize_ );
  for( int i = 0; i < BlockSize_; ++i )
    BaseGraph.ExtractMyRowView( i, BaseNumCols[i], BaseIndices[i] );

  vector< double* > Values( BlockSize_ );
  int NumEntries;

  Epetra_CrsMatrix & Mat = this->epetraObj();

  for( int i = 0; i < Size; ++i )
  {
    int NumCols = Cols_[i].size();
    Blocks_[i].resize( NumCols );

    for( int j = 0; j < BlockSize_; ++j )
      Mat.ExtractMyRowView( j+BlockSize_*i, NumEntries, Values[j] );

    for( int j = 0; j < NumCols; ++j )
    {
      Epetra_CrsMatrix * bMat = new Epetra_CrsMatrix( View, BaseGraph );

      for( int k = 0; k < BlockSize_; ++k )
        bMat->InsertMyValues( k, BaseNumCols[k], Values[k]+j*BaseNumCols[k], BaseIndices[k] );

      Blocks_[i][j] = new N_LAS_Matrix( bMat );
    }
  }

  if( AugmentCount_ )
  {
    AugmentGIDs_.resize(AugmentCount_);
    int AugStart = BlockSize_ * Size;
    for( int i = 0; i < AugmentCount_; ++i )
    {
      AugmentGIDs_[i] = Graph.RowMap().GID(AugStart+i);
      assert( AugmentGIDs_[i] != -1 );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::~N_LAS_BlockMatrix
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/13/04
//-----------------------------------------------------------------------------
N_LAS_BlockMatrix::~N_LAS_BlockMatrix()
{
  for( int i = 0; i < (int)(Blocks_.size()); ++i )
    for( int j = 0; j < (int)(Blocks_[i].size()); ++j )
      if( Blocks_[i][j] ) delete Blocks_[i][j];
  Blocks_.clear();
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
  for( int i = 0; i < (int)(Cols_[row].size()); ++i )
    if( Cols_[row][i] == col )
      return *Blocks_[row][i];

  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
      "Error!  N_LAS_BlockMatrix::block("<<row<<","<<col<<"):  This block does not exist!"
      );
  //abort(); //Row,Col Not Found
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::print
// Purpose       : Output
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::printPetraObject() const
{
  cout << "N_LAS_BlockMatrix Object (Size=" << NumBlockRows_ << ")\n";
  cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
  for( int i = 0; i < NumBlockRows_; ++i )
  {
    int NumCols = Cols_[i].size();
    for( int j = 0; j < NumCols; ++j )
    {
      cout << "Block[" << i << "][" << Cols_[i][j] << "]\n";
      Blocks_[i][j]->printPetraObject();
    }
  }
  cout << "Base Object\n";
  cout << *oDCRSMatrix_;
  cout << *aDCRSMatrix_;
  cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_BlockMatrix::insertIntoAugmentedColumn
// Purpose       : Augment matrix with block vector
// Special Notes : Zero based indexing for augmentedColumn
// Scope         : Public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 08/12/05
//-----------------------------------------------------------------------------
void N_LAS_BlockMatrix::replaceAugmentedColumn(int augmentedColumn, const N_LAS_BlockVector & vec)
{
  const Epetra_BlockMap & RowMap = const_cast<N_LAS_BlockVector&>(vec).epetraObj().Map();
  int NumRows = RowMap.NumMyElements();

  assert( AugmentGIDs_.size() != 0 );
  int Col = AugmentGIDs_[augmentedColumn];
  int lCol = aDCRSMatrix_->Graph().LCID(Col);
  for( int i = 0; i < NumRows; ++i )
  {
    //int Row = RowMap.LID(i);
    int Row = i;
    double Val = vec[i];
    aDCRSMatrix_->ReplaceMyValues( Row, 1, &Val, &lCol );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_blockMatrix
// Purpose       : Nonmember constructor that does not require a full sized
// block graph since it can be (and indeed must be) constructed from the
// BlockColumns and the BaseGraph.
// Special Notes : 
// Scope         : Public
// Creator       : Todd Coffey, SNL, 1414
// Creation Date : 06/08/09
//-----------------------------------------------------------------------------
RCP<N_LAS_BlockMatrix> N_LAS_blockMatrix( 
    const vector< vector<int> > & BlockColumns,
    Epetra_Map & BlockMap, // Defines the GIDs of the elements in the blocks
    const Epetra_CrsGraph & BaseGraph,
    int AugmentCount
    )
{
  // Create BlockGraph consistent with BaseGraph and BlockColumns
  int NumElements = BaseGraph.NumGlobalRows();
  int NumBlocks = BlockColumns.size();
  int numIndicesPerRow = 0;
  RCP<Epetra_CrsGraph> blockGraph = rcp(new Epetra_CrsGraph(
        Copy,
        dynamic_cast<Epetra_BlockMap&>(BlockMap),
        numIndicesPerRow
        )
      );
  blockGraph.release(); // This will leak memory
  for (int bi=0; bi<NumBlocks ; ++bi) { 
    // Loop over blocks by row
    for (int bj=0 ; bj<(int)(BlockColumns[bi].size()) ; ++bj) { 
      // Loop over blocks by column
      std::vector<int> Indices(NumElements*NumBlocks);
      for (int i=0 ; i<NumElements ; ++i) {
        int NumIndices = 0;
        // Loop over individual block by rows
        int GlobalRow = BlockMap.GID(bi*NumElements+i);
        std::vector<int> localIndices(NumElements);
        int localNumIndices;
        BaseGraph.ExtractGlobalRowCopy(BaseGraph.GRID(i), NumElements, localNumIndices, &localIndices[0]);
        for (int j=0 ; j<localNumIndices ; ++j) {
          // Loop over individual block by column
          Indices[NumIndices] = BlockMap.GID(BlockColumns[bi][bj]*NumElements+localIndices[j]);
          NumIndices++;
        }
        blockGraph->InsertGlobalIndices(GlobalRow,NumIndices,&Indices[0]);
      }
    }
  }
  //std::cout << "-----------------------------------------------------------" << std::endl;
  //std::cout << " N_LAS_blockMatrix:  blockGraph = " << std::endl;
  //blockGraph->Print(std::cout);
  //std::cout << "-----------------------------------------------------------" << std::endl;
  blockGraph->FillComplete();
  RCP<N_LAS_BlockMatrix> matrix = rcp(new N_LAS_BlockMatrix(
        NumBlocks,
        BlockColumns,
        *blockGraph,
        BaseGraph,
        AugmentCount
        )
      );
  return matrix;
}
