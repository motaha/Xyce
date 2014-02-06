//-----------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002, Sandia Corporation, Albuquerque, NM, USA.  Under the
// terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for
// use of this work by or on behalf of the U.S. Government.  Export of this
// program may require a license from the United States Government.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: BlockMatrix.C,v $
// Purpose       : This file contains unit tests for the N_LAS_BlockMatrix
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : 06/04/09
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.2 $
// Revision Date  : $Date: 2009/10/20 22:17:06 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------


#include <Teuchos_UnitTestHarness.hpp>
#include <BlockMatrix_Helpers.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

using Teuchos::RCP;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, createBaseMap ) {
  int numElements = 5;
  RCP<Epetra_Map> map = createBaseMap(numElements);
  TEST_EQUALITY( map->NumGlobalElements(), numElements );
  TEST_EQUALITY( map->NumMyElements(), numElements );
  TEST_EQUALITY_CONST( map->IndexBase(), 0 );
  vector<int> GIDs(numElements);
  map->MyGlobalElements(&GIDs[0]);
  vector<int> correctGIDs(numElements);
  for (int i=0; i<numElements ; ++i) {
    correctGIDs[i] =i;
  }
  TEST_COMPARE_ARRAYS( GIDs, correctGIDs );
}

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, createBlockMap ) {
  int numElements = 3;
  int numBlocks = 2;
  RCP<Epetra_Map> map = createBlockMap(numElements,numBlocks);
  TEST_EQUALITY( map->NumGlobalElements(), numElements*numBlocks );
  TEST_EQUALITY( map->NumMyElements(), numElements*numBlocks );
  TEST_EQUALITY_CONST( map->IndexBase(), 0 );
  vector<int> GIDs(numElements*numBlocks);
  map->MyGlobalElements(&GIDs[0]);
  vector<int> correctGIDs(6);
  correctGIDs[0] = 0;
  correctGIDs[1] = 1;
  correctGIDs[2] = 2;
  correctGIDs[3] = 10;
  correctGIDs[4] = 11;
  correctGIDs[5] = 12;
  TEST_COMPARE_ARRAYS( GIDs, correctGIDs );
}

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, createBaseGraph ) {
  int numElements = 3;
  vector<vector<int> > nonzeros;
  nonzeros.resize(numElements);
  nonzeros[0].resize(2);
  nonzeros[0][0] = 0;
  nonzeros[0][1] = 2;
  nonzeros[1].resize(1);
  nonzeros[1][0] = 2;
  nonzeros[2].resize(3);
  nonzeros[2][0] = 0;
  nonzeros[2][1] = 1;
  nonzeros[2][2] = 2;

  RCP<Epetra_CrsGraph> graph = createBaseGraph(numElements,nonzeros);
  int row;
  int length = numElements;
  int numIndices;
  int indices[length];
  {
    row = 0;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 2 );
    TEST_EQUALITY_CONST( indices[0], 0 );
    TEST_EQUALITY_CONST( indices[1], 2 );
  }
  {
    row = 1;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 1 );
    TEST_EQUALITY_CONST( indices[0], 2 );
  }
  {
    row = 2;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 3 );
    TEST_EQUALITY_CONST( indices[0], 0 );
    TEST_EQUALITY_CONST( indices[1], 1 );
    TEST_EQUALITY_CONST( indices[2], 2 );
  }
}

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, createBlockGraph ) {
  // NumElements = 3
  // NumBlocks = 2
  // nonZeroColumns = [ 0, 1 ], [ 1 ]
  // baseNonZeros = [ 0, 2 ], [ 2 ], [ 0, 1, 2 ]
  // Create the following:
  // blockGraph = 
  // [ 0,    2, | 3,    5 ]
  // [       2, |       5 ]
  // [ 0, 1, 2, | 3, 4, 5 ]
  // [--------------------]
  // [          | 3,    5 ]
  // [          |       5 ]
  // [          | 3, 4, 5 ]
  int numElements = 3;
  int numBlocks = 2;
  vector<vector<int> > nonZeroColumns;
  vector<vector<int> > baseNonZeros;
  // First we'll fill the baseNonZeros:
  baseNonZeros.resize(numElements);
  baseNonZeros[0].resize(2);
  baseNonZeros[0][0] = 0;
  baseNonZeros[0][1] = 2;
  baseNonZeros[1].resize(1);
  baseNonZeros[1][0] = 2;
  baseNonZeros[2].resize(3);
  baseNonZeros[2][0] = 0;
  baseNonZeros[2][1] = 1;
  baseNonZeros[2][2] = 2;
  // Second, we'll fill the nonZeroColumns
  nonZeroColumns.resize(numBlocks);
  nonZeroColumns[0].resize(2);
  nonZeroColumns[0][0] = 0;
  nonZeroColumns[0][1] = 1;
  nonZeroColumns[1].resize(1);
  nonZeroColumns[1][0] = 1;
  // Create the big graph:
  RCP<Epetra_CrsGraph> graph = createBlockGraph(numElements,numBlocks,nonZeroColumns,baseNonZeros);
  int row;
  int length = numElements*numBlocks;
  int numIndices;
  int indices[length];
  {
    row = 0;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 4 );
    TEST_EQUALITY_CONST( indices[0], 0 );
    TEST_EQUALITY_CONST( indices[1], 2 );
    TEST_EQUALITY_CONST( indices[2], 10 );
    TEST_EQUALITY_CONST( indices[3], 12 );
  }
  {
    row = 1;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 2 );
    TEST_EQUALITY_CONST( indices[0], 2 );
    TEST_EQUALITY_CONST( indices[1], 12 );
  }
  {
    row = 2;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 6 );
    TEST_EQUALITY_CONST( indices[0], 0 );
    TEST_EQUALITY_CONST( indices[1], 1 );
    TEST_EQUALITY_CONST( indices[2], 2 );
    TEST_EQUALITY_CONST( indices[3], 10 );
    TEST_EQUALITY_CONST( indices[4], 11 );
    TEST_EQUALITY_CONST( indices[5], 12 );
  }
  // repeat of above
  {
    row = 10;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 2 );
    TEST_EQUALITY_CONST( indices[0], 10 );
    TEST_EQUALITY_CONST( indices[1], 12 );
  }
  {
    row = 11;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 1 );
    TEST_EQUALITY_CONST( indices[0], 12 );
  }
  {
    row = 12;
    graph->ExtractGlobalRowCopy(row,length,numIndices,indices);
    TEST_EQUALITY_CONST( numIndices, 3 );
    TEST_EQUALITY_CONST( indices[0], 10 );
    TEST_EQUALITY_CONST( indices[1], 11 );
    TEST_EQUALITY_CONST( indices[2], 12 );
  }
}

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, createBlockMatrix ) {
  // NumElements = 3
  // NumBlocks = 2
  // nonZeroColumns = [ 0, 1 ], [ 1 ]
  // baseNonZeros = [ 0, 2 ], [ 2 ], [ 0, 1, 2 ]
  // Create the following:
  // blockGraph = 
  // [ 0,    2, | 3,    5 ]
  // [       2, |       5 ]
  // [ 0, 1, 2, | 3, 4, 5 ]
  // [--------------------]
  // [          | 3,    5 ]
  // [          |       5 ]
  // [          | 3, 4, 5 ]
  int numElements = 3;
  int numBlocks = 2;
  vector<vector<int> > nonZeroColumns;
  vector<vector<int> > baseNonZeros;
  // First we'll fill the baseNonZeros:
  baseNonZeros.resize(numElements);
  baseNonZeros[0].resize(2);
  baseNonZeros[0][0] = 0;
  baseNonZeros[0][1] = 2;
  baseNonZeros[1].resize(1);
  baseNonZeros[1][0] = 2;
  baseNonZeros[2].resize(3);
  baseNonZeros[2][0] = 0;
  baseNonZeros[2][1] = 1;
  baseNonZeros[2][2] = 2;
  // Second, we'll fill the nonZeroColumns
  nonZeroColumns.resize(numBlocks);
  nonZeroColumns[0].resize(2);
  nonZeroColumns[0][0] = 0;
  nonZeroColumns[0][1] = 1;
  nonZeroColumns[1].resize(1);
  nonZeroColumns[1][0] = 1;
  // Now we create the matrix:
  RCP<N_LAS_BlockMatrix> matrix = createBlockMatrix(numElements,numBlocks,nonZeroColumns,baseNonZeros);
  // Now we check that the global graph is correct:
  {
    Epetra_CrsMatrix & crsMatrix = matrix->epetraObj();
    const Epetra_CrsGraph& graph = crsMatrix.Graph();
    int row;
    int length = numElements*numBlocks;
    int numIndices;
    int indices[length];
    {
      row = 0;
      graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
      TEST_EQUALITY_CONST( numIndices, 4 );
      TEST_EQUALITY_CONST( indices[0], 0 );
      TEST_EQUALITY_CONST( indices[1], 2 );
      TEST_EQUALITY_CONST( indices[2], 10 );
      TEST_EQUALITY_CONST( indices[3], 12 );
    }
    {
      row = 1;
      graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
      TEST_EQUALITY_CONST( numIndices, 2 );
      TEST_EQUALITY_CONST( indices[0], 2 );
      TEST_EQUALITY_CONST( indices[1], 12 );
    }
    {
      row = 2;
      graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
      TEST_EQUALITY_CONST( numIndices, 6 );
      TEST_EQUALITY_CONST( indices[0], 0 );
      TEST_EQUALITY_CONST( indices[1], 1 );
      TEST_EQUALITY_CONST( indices[2], 2 );
      TEST_EQUALITY_CONST( indices[3], 10 );
      TEST_EQUALITY_CONST( indices[4], 11 );
      TEST_EQUALITY_CONST( indices[5], 12 );
    }
    // repeat of above
    {
      row = 10;
      graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
      TEST_EQUALITY_CONST( numIndices, 2 );
      TEST_EQUALITY_CONST( indices[0], 10 );
      TEST_EQUALITY_CONST( indices[1], 12 );
    }
    {
      row = 11;
      graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
      TEST_EQUALITY_CONST( numIndices, 1 );
      TEST_EQUALITY_CONST( indices[0], 12 );
    }
    {
      row = 12;
      graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
      TEST_EQUALITY_CONST( numIndices, 3 );
      TEST_EQUALITY_CONST( indices[0], 10 );
      TEST_EQUALITY_CONST( indices[1], 11 );
      TEST_EQUALITY_CONST( indices[2], 12 );
    }
  }
  // Pull out each block and make sure the blocks have the right graphs
  for (int bi=0 ; bi<nonZeroColumns.size() ; ++bi) {
    for (int bj=0 ; bj<nonZeroColumns[bi].size() ; ++bj) {
      N_LAS_Matrix & bmat = matrix->block(bi,nonZeroColumns[bi][bj]);
      Epetra_CrsMatrix & crsMatrix = bmat.epetraObj();
      const Epetra_CrsGraph& graph = crsMatrix.Graph();
      int row;
      int length = numElements*numBlocks;
      int numIndices;
      int indices[length];
      {
        row = 0;
        graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
        TEST_EQUALITY_CONST( numIndices, 2 );
        TEST_EQUALITY_CONST( indices[0], 0 );
        TEST_EQUALITY_CONST( indices[1], 2 );
      }
      {
        row = 1;
        graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
        TEST_EQUALITY_CONST( numIndices, 1 );
        TEST_EQUALITY_CONST( indices[0], 2 );
      }
      {
        row = 2;
        graph.ExtractGlobalRowCopy(row,length,numIndices,indices);
        TEST_EQUALITY_CONST( numIndices, 3 );
        TEST_EQUALITY_CONST( indices[0], 0 );
        TEST_EQUALITY_CONST( indices[1], 1 );
        TEST_EQUALITY_CONST( indices[2], 2 );
      }
    }
  }
  TEST_THROW( N_LAS_Matrix & bmat = matrix->block(1,0), std::logic_error );
}

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, createDiagonalBlockMatrix ) {
  // NumElements = 3
  // NumBlocks = 2
  // nonZeroColumns = [ 0, 1 ], [ 1 ]
  // baseNonZeros = [ 0, 2 ], [ 2 ], [ 0, 1, 2 ]
  // Create the following:
  // blockGraph = 
  // [ 0,    2, | 3,    5 ]
  // [       2, |       5 ]
  // [ 0, 1, 2, | 3, 4, 5 ]
  // [--------------------]
  // [          | 3,    5 ]
  // [          |       5 ]
  // [          | 3, 4, 5 ]
  int numElements = 3;
  int numBlocks = 2;
  vector<vector<int> > nonZeroColumns;
  vector<vector<int> > baseNonZeros;
  // First we'll fill the baseNonZeros:
  baseNonZeros.resize(numElements);
  baseNonZeros[0].resize(2);
  baseNonZeros[0][0] = 0;
  baseNonZeros[0][1] = 2;
  baseNonZeros[1].resize(1);
  baseNonZeros[1][0] = 2;
  baseNonZeros[2].resize(3);
  baseNonZeros[2][0] = 0;
  baseNonZeros[2][1] = 1;
  baseNonZeros[2][2] = 2;
  // Second, we'll fill the nonZeroColumns
  nonZeroColumns.resize(numBlocks);
  nonZeroColumns[0].resize(2);
  nonZeroColumns[0][0] = 0;
  nonZeroColumns[0][1] = 1;
  nonZeroColumns[1].resize(1);
  nonZeroColumns[1][0] = 1;
  // Now we create the matrix:
  RCP<N_LAS_BlockMatrix> matrix = createBlockMatrix(numElements,numBlocks,nonZeroColumns,baseNonZeros);
  int length=1;
  std::vector<double> coeffs; coeffs.resize(length);
  std::vector<int> indices; indices.resize(length);
  for (int bi=0 ; bi<nonZeroColumns.size() ; ++bi) {
    for (int bj=0 ; bj<nonZeroColumns[bi].size() ; ++bj) {
      N_LAS_Matrix & bmat = matrix->block(bi,nonZeroColumns[bi][bj]);
      for (int row=0 ; row < 2 ; ++row) {
        coeffs[0] = 1.0;
        indices[0] = row;
        bmat.replaceLocalRow(row,length,&coeffs[0],&indices[0]);
      }
    }
  }
}

TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, matvec ) {
  // Does the matvec work correctly?
  // matrix * vector = [ block(0,0)*vector(0) + block(0,1)*vector(1) ]
  //                   [                        block(0,1)*vector(1) ]
  int numElements = 3;
  int numBlocks = 2;
  vector<vector<int> > nonZeroColumns;
  vector<vector<int> > baseNonZeros;
  // First we'll fill the baseNonZeros:
  baseNonZeros.resize(numElements);
  baseNonZeros[0].resize(1);
  baseNonZeros[0][0] = 0;
  // Second, we'll fill the nonZeroColumns
  nonZeroColumns.resize(numBlocks);
  nonZeroColumns[0].resize(2);
  nonZeroColumns[0][0] = 0;
  nonZeroColumns[0][1] = 1;
  nonZeroColumns[1].resize(1);
  nonZeroColumns[1][0] = 1;
  RCP<Epetra_Map> baseMap = createBaseMap(numElements);
  RCP<Epetra_Map> blockMap = createBlockMap(numElements,numBlocks);
  // First we create the block vector
  RCP<N_LAS_BlockVector> input_vector = rcp(new N_LAS_BlockVector( numBlocks, *blockMap, *baseMap));
  RCP<N_LAS_BlockVector> output_vector = rcp(new N_LAS_BlockVector( numBlocks, *blockMap, *baseMap));
  input_vector->putScalar(0.0);
  output_vector->putScalar(0.0);
  {
    N_LAS_Vector & b0 = input_vector->block(0);
    b0[0] = 2.0;
    b0[1] = 3.0;
    b0[2] = 4.0;
    N_LAS_Vector & b1 = input_vector->block(1);
    b1[0] = 5.0;
    b1[1] = 6.0;
    b1[2] = 7.0;

  }
  // Now we create the block matrix:
  RCP<N_LAS_BlockMatrix> matrix = createBlockMatrix(numElements,numBlocks,nonZeroColumns,baseNonZeros);
  // Now we fill them with data.
  matrix->put(1.0);
  matrix->fillComplete();
  matrix->matvec(false,*input_vector,*output_vector);
  TEST_EQUALITY_CONST( (*output_vector)[0], 7.0 );
  TEST_EQUALITY_CONST( (*output_vector)[1], 0.0 );
  TEST_EQUALITY_CONST( (*output_vector)[2], 0.0 );
  TEST_EQUALITY_CONST( (*output_vector)[3], 5.0 );
  TEST_EQUALITY_CONST( (*output_vector)[4], 0.0 );
  TEST_EQUALITY_CONST( (*output_vector)[5], 0.0 );
}

/*
TEUCHOS_UNIT_TEST( N_LAS_BlockMatrix, create ) {
  int numElements = 2;
  int numBlocks = 2;
  // We need to specify which blocks are non-empty (all)
  std::vector<vector<int> > Cols;
  Cols.resize(numBlocks);
  for (int i=0 ; i<numBlocks ; ++i) {
    Cols[i].resize(numBlocks);
    for (int j=0 ; j<numBlocks ; ++j) {
      Cols[i][j] = j;
    }
  }
  // Then we create a N_LAS_BlockMatrix:
  RCP<N_LAS_BlockMatrix> matrix = createBlockMatrix(numElements,numBlocks,Cols);
  TEST_EQUALITY_CONST( !is_null(matrix), true );
}
*/
