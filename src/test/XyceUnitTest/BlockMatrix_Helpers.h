#ifndef BLOCK_MATRIX_HELPERS_H
#define BLOCK_MATRIX_HELPERS_H

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
#include <vector>
using std::vector;

class Epetra_Map;
class Epetra_CrsGraph;
class N_LAS_BlockMatrix;

RCP<Epetra_Map> createBaseMap(int NumElements);
RCP<Epetra_Map> createBlockMap(int NumElements, int NumBlocks);
RCP<Epetra_CrsGraph> createBaseGraph(int NumElements, const vector<vector<int> > & nonzeros);
RCP<Epetra_CrsGraph> createBlockGraph(
    int NumElements, 
    int NumBlocks, 
    const vector<vector<int> > & nonZeroColumns, // NonZeroColumns.size() == NumBlocks
    const vector<vector<int> > & baseNonZeros // baseNonZeros.size() == NumElements
    );

RCP<N_LAS_BlockMatrix> createBlockMatrix(
    int NumElements, 
    int NumBlocks, 
    const vector<vector<int> > & nonZeroColumns,
    const vector<vector<int> > & baseNonZeros
    );


#endif // BLOCK_MATRIX_HELPERS_H

