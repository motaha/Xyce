
#include <BlockMatrix_Helpers.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <N_LAS_BlockMatrix.h>

using Teuchos::rcp;

RCP<Epetra_Map> createBaseMap(int NumElements)
{
  Epetra_SerialComm Comm;
  vector<int> baseGIDs;
  baseGIDs.resize(NumElements);
  for (int i=0 ; i<NumElements ; ++i) {
    baseGIDs[i] = i;
  }
  int baseIndex = 0;
  RCP<Epetra_Map> baseMap = rcp(new Epetra_Map(
        NumElements, // Global elements
        NumElements, // Local elements
        &baseGIDs[0], // GIDs
        baseIndex, // 0 or 1
        Comm // communicator
        )
      );
  return baseMap;
}
RCP<Epetra_Map> createBlockMap(int NumElements, int NumBlocks)
{
  Epetra_SerialComm Comm;
  RCP<Epetra_Map> baseMap = createBaseMap(NumElements);
  vector<int> baseGIDs;
  baseGIDs.resize(baseMap->NumGlobalElements()); // Should be NumElements
  baseMap->MyGlobalElements(&baseGIDs[0]);
  int MaxGID = baseMap->MaxAllGID();
  int offset = 1;
  while (offset <= MaxGID ) { 
    offset *= 10;
  }
  vector<int> blockGIDs;
  blockGIDs.resize(NumElements*NumBlocks);
  for (int i=0 ; i<NumBlocks ; ++i) {
    for (int j=0 ; j<NumElements ; ++j) {
      blockGIDs[i*NumElements+j] = baseGIDs[j]+offset*i;
    }
  }
  int blockIndex = 0;
  RCP<Epetra_Map> blockMap = rcp(new Epetra_Map(
        NumElements*NumBlocks, // Global elements
        NumElements*NumBlocks, // Local elements
        &blockGIDs[0], // GIDs
        blockIndex, // 0 or 1
        Comm // communicator
        )
      );
  return blockMap;
}



RCP<Epetra_CrsGraph> createBaseGraph(int NumElements, const vector<vector<int> >& nonzeros)
{
  RCP<Epetra_Map> baseMap = createBaseMap(NumElements);
  int numIndicesPerRow = 0;
  RCP<Epetra_CrsGraph> baseGraph = rcp(new Epetra_CrsGraph(
        Copy,
        dynamic_cast<Epetra_BlockMap&>(*baseMap),
        numIndicesPerRow
        )
      );
  for (int i=0 ; i<NumElements; ++i) {
    int Indices[NumElements];
    int NumIndices = 0;
    for (int j=0 ; j<nonzeros[i].size() ; ++j) {
      NumIndices++;
      Indices[j] = baseMap->GID(nonzeros[i][j]);
    }
    baseGraph->InsertGlobalIndices(baseMap->GID(i),NumIndices,&Indices[0]);
  }
  return baseGraph;
}
RCP<Epetra_CrsGraph> createBlockGraph(
    int NumElements, 
    int NumBlocks, 
    const vector<vector<int> > & nonZeroColumns, // NonZeroColumns.size() == NumBlocks
    const vector<vector<int> > & baseNonZeros // baseNonZeros.size() == NumElements
    )
{
  RCP<Epetra_Map> blockMap = createBlockMap(NumElements,NumBlocks);
  int numIndicesPerRow = 0;
  RCP<Epetra_CrsGraph> blockGraph = rcp(new Epetra_CrsGraph(
        Copy,
        dynamic_cast<Epetra_BlockMap&>(*blockMap),
        numIndicesPerRow
        )
      );
  for (int bi=0; bi<NumBlocks ; ++bi) { 
    // Loop over blocks by row
    for (int bj=0 ; bj<nonZeroColumns[bi].size() ; ++bj) { 
      // Loop over blocks by column
      int Indices[NumElements*NumBlocks];
      for (int i=0 ; i<NumElements ; ++i) {
        int NumIndices = 0;
        // Loop over individual block by rows
        int GlobalRow = blockMap->GID(bi*NumElements+i);
        for (int j=0 ; j<baseNonZeros[i].size() ; ++j) {
          // Loop over individual block by column
          Indices[NumIndices] = blockMap->GID(nonZeroColumns[bi][bj]*NumElements+baseNonZeros[i][j]);
          NumIndices++;
        }
        blockGraph->InsertGlobalIndices(GlobalRow,NumIndices,&Indices[0]);
      }
    }
  }
  return blockGraph;
}

RCP<N_LAS_BlockMatrix> createBlockMatrix(
    int NumElements, 
    int NumBlocks, 
    const vector<vector<int> > & nonZeroColumns,
    const vector<vector<int> > & baseNonZeros
    )
{
  RCP<Epetra_Map> blockMap = createBlockMap(NumElements,NumBlocks);
  RCP<Epetra_CrsGraph> baseGraph = createBaseGraph(NumElements,baseNonZeros);
  baseGraph->FillComplete();
  baseGraph.release(); // This will leak memory
  //RCP<Epetra_CrsGraph> blockGraph = createBlockGraph(NumElements,NumBlocks,nonZeroColumns,baseNonZeros);
  //blockGraph->FillComplete();
  //blockGraph.release(); // This will leak memory
  //RCP<N_LAS_BlockMatrix> matrix = rcp( new N_LAS_BlockMatrix( 
  //      NumBlocks,
  //      nonZeroColumns,
  //      *blockGraph,
  //      *baseGraph
  //      ) 
  //    );
  RCP<N_LAS_BlockMatrix> matrix = N_LAS_blockMatrix( 
        nonZeroColumns,
        *blockMap,
        *baseGraph
        );
  return matrix;
}
