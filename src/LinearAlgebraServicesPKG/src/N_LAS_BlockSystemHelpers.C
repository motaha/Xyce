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
// Filename       : $RCSfile: N_LAS_BlockSystemHelpers.C,v $
//
// Purpose        : This is collection of non-member functions that help
//                  in the construction of block linear systems, like those
//                  found in AC analysis.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 06/22/11
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:45 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_UTL_Misc.h>

#include <N_LAS_BlockSystemHelpers.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_BlockMatrix.h>

#include <N_LAS_System.h>
#include <N_LAS_QueryUtil.h>

#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>

#ifdef Xyce_VERBOSE_LINEAR
 #include <N_ERH_ErrorMgr.h>
#endif

//-----------------------------------------------------------------------------
// Function      : createBlockVector
// Purpose       : A helper function for creating a block vector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<N_LAS_BlockVector> createBlockVector( int numBlocks, N_LAS_Vector& subBlockVector )
{
   // Create the parallel block maps based on the distribution of the subBlockVector
   std::vector<Teuchos::RCP<N_PDS_ParMap> > blockMaps
     = createBlockParMaps( numBlocks, *(subBlockVector.pmap()), *(subBlockVector.omap()) );

   // Create the new N_LAS_BlockVector using the parallel maps
   Teuchos::RCP<N_LAS_BlockVector> newvector 
     = Teuchos::rcp( new N_LAS_BlockVector( numBlocks, blockMaps[0], blockMaps[1],
                                            Teuchos::rcp(subBlockVector.pmap(),false),
                                            Teuchos::rcp(subBlockVector.omap(),false)
                                          ) );

   // Return the new block vector
   return newvector;
}

//-----------------------------------------------------------------------------
// Function      : copyToBlockVector
// Purpose       : A helper function that copies the array of N_LAS_Vectors
//               : into an N_LAS_BlockVector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyToBlockVector( std::vector<Teuchos::RCP<N_LAS_Vector> >& inputVectors, N_LAS_BlockVector& blockVector )
{
  int inputVecLength = inputVectors.size();
  int blockVecLength = blockVector.blockCount();

  // If the number of blocks is not the same, throw an error.
  if (inputVecLength != blockVecLength) {}
   
  // Loop over the vectors and copy each one.  
  for (int i=0 ; i<blockVecLength ; ++i)
  {
    blockVector.block(i) = *(inputVectors[i]);
  }
}

//-----------------------------------------------------------------------------
// Function      : copyFromBlockVector
// Purpose       : A helper function that copies a N_LAS_BlockVector to an 
//               : array of N_LAS_Vectors.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyFromBlockVector( N_LAS_BlockVector& blockVector, std::vector<Teuchos::RCP<N_LAS_Vector> >& outputVectors )
{
  int outputVecLength = outputVectors.size();
  int blockVecLength = blockVector.blockCount();

  // If the number of blocks is not the same, throw an error.
  if (outputVecLength != blockVecLength) {}
 
  // Loop over the vectors and copy each one. 
  for (int i=0 ; i<blockVecLength ; ++i)
  {
    *(outputVectors[i]) = blockVector.block(i);
  }
}

//-----------------------------------------------------------------------------
// Function      : createBlockParMaps
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns both the map and overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//----------------------------------------------------------------------------- 
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap )
{
   // The information about the current maps, given by pmap and omap (overlap map)
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int olocalBlockSize = omap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();
   int oglobalBlockSize = omap.numGlobalEntities();

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();
   int oBaseIndex = omap.indexBase();

   // Compute the offset needed for global indexing.
   int maxGID = pmap.petraMap()->MaxAllGID();
   int numProcs = pmap.pdsComm()->numProc();
   int offset = 1;
   while( offset <= maxGID ) offset *= 10;

   // Determine size of block maps
   int numGlobalElements = numBlocks*globalBlockSize;
   int onumGlobalElements = numBlocks*oglobalBlockSize;
   int numLocalElements = numBlocks*localBlockSize;
   int onumLocalElements = numBlocks*olocalBlockSize;

   // Initialize vectors to hold the GIDs for the original and block map.
   vector<int> BaseGIDs(localBlockSize), oBaseGIDs(olocalBlockSize);
   vector<int> GIDs(numLocalElements), oGIDs(onumLocalElements);

   // Extract the global indices.
   pmap.petraMap()->MyGlobalElements( &BaseGIDs[0] );
   omap.petraMap()->MyGlobalElements( &oBaseGIDs[0] );
   
   int gnd_node = 0;  // Will be decremented before first use.

   for( int i = 0; i < numBlocks; ++i )
   {
     // Setting up GIDs for the map without overlap and ground nodes 
     for( int j = 0; j < localBlockSize; ++j )
     {
       GIDs[i*localBlockSize+j] = BaseGIDs[j] + offset*i;
     }

     // Setting up GIDs for the map with overlap and ground nodes
     for( int j = 0; j < olocalBlockSize-1; ++j )
     {
       oGIDs[i*olocalBlockSize+j] = oBaseGIDs[j] + offset*i;
     }
     // Assuming the last GID is the ground node (-1)
     gnd_node--;
     oGIDs[(i+1)*olocalBlockSize-1] = gnd_node;
   }

   // Adapt base index for unique ground node numbering.
   oBaseIndex = Xycemin( oBaseIndex, gnd_node );

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap, oBlockMap;
   blockMap = Teuchos::rcp(new N_PDS_ParMap(numGlobalElements, numLocalElements, GIDs, BaseIndex, pmap.pdsComm()));
   oBlockMap = Teuchos::rcp(new N_PDS_ParMap(onumGlobalElements, onumLocalElements, oGIDs, oBaseIndex, pmap.pdsComm()));

   std::vector<Teuchos::RCP<N_PDS_ParMap> > allMaps;
   allMaps.push_back(blockMap);
   allMaps.push_back(oBlockMap);

   return allMaps;
}

//-----------------------------------------------------------------------------
// Function      : createBlockParMaps
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns both the map and overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//----------------------------------------------------------------------------- 
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps2( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap )
{
   // The information about the current maps, given by pmap and omap (overlap map)
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int olocalBlockSize = omap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();
   int oglobalBlockSize = omap.numGlobalEntities();
   int overlapSize = (olocalBlockSize - 1) - localBlockSize;

   // If overlapSize == -1, then the overlap map is the same size as the parallel map.

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();
   int oBaseIndex = omap.indexBase();

   // Compute the offset needed for global indexing.
   int maxGID = pmap.petraMap()->MaxAllGID();
   int numProcs = pmap.pdsComm()->numProc();
   int offset = 1;
   while( offset <= maxGID ) offset *= 10;

   // Determine size of block maps
   int numGlobalElements = numBlocks*globalBlockSize;
   int onumGlobalElements = numBlocks*(oglobalBlockSize-numProcs) + numProcs;
   int numLocalElements = numBlocks*localBlockSize;
   int onumLocalElements = numBlocks*(olocalBlockSize-1) + 1;

   // Initialize vectors to hold the GIDs for the original and block map.
   vector<int> BaseGIDs(localBlockSize), oBaseGIDs(olocalBlockSize);
   vector<int> GIDs(numLocalElements), oGIDs(onumLocalElements);

   // Extract the global indices.
   pmap.petraMap()->MyGlobalElements( &BaseGIDs[0] );
   omap.petraMap()->MyGlobalElements( &oBaseGIDs[0] );
   
   for( int i = 0; i < numBlocks; ++i )
   {
     // Setting up GIDs for the map without overlap and ground nodes 
     for( int j = 0; j < localBlockSize; ++j )
     {
       GIDs[i*localBlockSize+j] = BaseGIDs[j] + offset*i;
       
       // Load all the local elements first in the block map.
       // All the external nodes should be inserted at the end of the GID list.
       oGIDs[i*localBlockSize+j] = oBaseGIDs[j] + offset*i;
     }

     // Now insert the external (overlap) nodes at the end of the GID list.
     for ( int j = localBlockSize, jj=0; j < olocalBlockSize-1; ++j, ++jj )
     {
       oGIDs[numLocalElements+ i*overlapSize + jj] = oBaseGIDs[j] + offset*i;
     }
   }
   // Insert ground node.
   oGIDs[onumLocalElements-1] = -1;

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap, oBlockMap;
   blockMap = Teuchos::rcp(new N_PDS_ParMap(numGlobalElements, numLocalElements, GIDs, BaseIndex, pmap.pdsComm()));
   oBlockMap = Teuchos::rcp(new N_PDS_ParMap(onumGlobalElements, onumLocalElements, oGIDs, oBaseIndex, pmap.pdsComm()));

   std::vector<Teuchos::RCP<N_PDS_ParMap> > allMaps;
   allMaps.push_back(blockMap);
   allMaps.push_back(oBlockMap);

   return allMaps;
}


//-----------------------------------------------------------------------------
// Function      : createBlockParMap
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns only the map and not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//----------------------------------------------------------------------------- 
Teuchos::RCP<N_PDS_ParMap> createBlockParMap( int numBlocks, N_PDS_ParMap& pmap )
{
   // The information about the current maps, given by pmap 
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();

   // Compute the offset needed for global indexing.
   int maxGID = pmap.petraMap()->MaxAllGID();
   int numProcs = pmap.pdsComm()->numProc();
   int offset = 1;
   while( offset <= maxGID ) offset *= 10;

   // Determine size of block maps
   int numGlobalElements = numBlocks*globalBlockSize;
   int numLocalElements = numBlocks*localBlockSize;

   // Initialize vectors to hold the GIDs for the original and block map.
   vector<int> BaseGIDs(localBlockSize);
   vector<int> GIDs(numLocalElements);

   // Extract the global indices.
   pmap.petraMap()->MyGlobalElements( &BaseGIDs[0] );
   
   for( int i = 0; i < numBlocks; ++i )
   {
     // Setting up GIDs for the map without overlap and ground nodes 
     for( int j = 0; j < localBlockSize; ++j )
     {
       GIDs[i*localBlockSize+j] = BaseGIDs[j] + offset*i;
     }
   }

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap = 
     Teuchos::rcp(new N_PDS_ParMap(numGlobalElements, numLocalElements, GIDs, BaseIndex, pmap.pdsComm()));

   return blockMap;
}

//-----------------------------------------------------------------------------
// Function      : createBlockGraph
// Purpose       : A helper function for creating block parallel graphs.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<Epetra_CrsGraph> createBlockGraph( int offset, std::vector<std::vector<int> >& blockPattern, 
                                                N_PDS_ParMap& blockMap, const Epetra_CrsGraph& baseGraph )
{
  int numBlockRows = blockPattern.size();
  int numMyBaseRows = baseGraph.NumMyRows();
  int maxIndices = baseGraph.MaxNumIndices();
 
  int maxBlockCols = blockPattern[0].size();
  for (int i=1; i<numBlockRows; ++i) { 
    int cols=blockPattern[i].size();
    if (cols > maxBlockCols)
      maxBlockCols = cols;
  }
 
  //Construct block graph based on  [All graphs are the same, so only one needs to be made]
  Teuchos::RCP<Epetra_CrsGraph> newGraph = rcp(new Epetra_CrsGraph( Copy, *(blockMap.petraBlockMap()), 0 ));
  
  std::vector<int> indices(maxIndices);
  int shift=0, index=0, baseRow=0, blockRow=0, numIndices=0;
  int maxNNZs = maxIndices*maxBlockCols;
  std::vector<int> newIndices(maxNNZs);  // Make as large as the combined maximum of indices and column blocks

  for( int j = 0; j < numMyBaseRows; ++j )
  {
    // Extract the base entries from the base row.
    baseRow = baseGraph.GRID(j);
    baseGraph.ExtractGlobalRowCopy( baseRow, maxIndices, numIndices, &indices[0] );

    for( int i = 0; i < numBlockRows; ++i )
    {
      // For this harmonic, which row will be inserted.
      blockRow = baseRow + offset*i;

      int numBlockCols = blockPattern[i].size();

      // Find all entries from a row before inserting it.
      for( int k = 0; k < numBlockCols; ++k )
      {
        // Find which block column to start at.
        shift = blockPattern[i][k]*offset;  // Actual column index.
        index = k*numIndices;  // Pointer to next block of column indices.
        for( int kk = 0; kk < numIndices; ++kk ) newIndices[index+kk] = indices[kk] + shift;
      }

      // Insert entire row for all blocks.
      newGraph->InsertGlobalIndices( blockRow, numBlockCols*numIndices, &newIndices[0] );
    }
  }
  newGraph->TransformToLocal();
 
  return newGraph;
}
