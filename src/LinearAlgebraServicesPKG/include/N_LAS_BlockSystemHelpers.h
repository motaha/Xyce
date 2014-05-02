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
// Filename       : $RCSfile: N_LAS_BlockSystemHelpers.h,v $
//
// Purpose        : This is collection of non-member functions that help
//                  in the construction of block linear systems, like those
//                  found in AC or HB analysis.
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
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_BLOCKSYSTEMHELPERS_H
#define  Xyce_LAS_BLOCKSYSTEMHELPERS_H

// ---------- Standard Includes ----------

#include <Teuchos_RCP.hpp>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_MultiVector;
class N_LAS_BlockVector;
class N_LAS_Matrix;
class N_LAS_BlockMatrix;
class N_LAS_System;
class N_PDS_ParMap;
class Epetra_CrsGraph;

//-----------------------------------------------------------------------------
// Function      : generateOffset 
// Purpose       : A helper function that standardizes how offsets are computed
// Special Notes : Block maps, graphs, vectors, and matrices require a global
//               : numbering scheme.  Xyce uses an offset index to space the 
//               : global ids apart so that they are unique.  The computation
//               : of this offset is performed by this function using the
//               : N_PDS_ParMap from the base block object.
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 10/8/13
//-----------------------------------------------------------------------------
int generateOffset( const N_PDS_ParMap& baseMap );

//-----------------------------------------------------------------------------
// Function      : createBlockVector
// Purpose       : A helper function for creating a block vector.
// Special Notes : This block vector construction is a numBlock replicate of the parallel
//               : distribution defined by the maps in the subBlockVector.
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<N_LAS_BlockVector> createBlockVector( int numBlocks, N_LAS_Vector& subBlockVector, int augmentRows=0 );

//-----------------------------------------------------------------------------
// Function      : createBlockParMaps
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns both the map and overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap );
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps2( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap );

//-----------------------------------------------------------------------------
// Function      : createBlockParMap
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns only the map, not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<N_PDS_ParMap> createBlockParMap( int numBlocks, N_PDS_ParMap& pmap,
                                              int augmentRows=0, std::vector<int>* augmentedGIDs = 0 );

//-----------------------------------------------------------------------------
// Function      : createBlockGraph
// Purpose       : A helper function for creating block parallel graphs.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<Epetra_CrsGraph> createBlockGraph( int offset, std::vector<std::vector<int> >& blockPattern,
                                                N_PDS_ParMap& blockMap, const Epetra_CrsGraph& baseGraph );

//-----------------------------------------------------------------------------
// Function      : createBlockFreqERFParMap
// Purpose       : A helper function for creating block parallel maps for 
//               : the frequency domain.  The map generated here has all the
//               : harmonics for one time point grouped together in expanded
//               : real form.
//               : This function returns only the map, not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<N_PDS_ParMap> createBlockFreqERFParMap( int numHarmonics, N_PDS_ParMap& pmap );

//-----------------------------------------------------------------------------
// Function      : copyToBlockVector
// Purpose       : A helper function that copies the array of N_LAS_Vectors
//               : into a N_LAS_BlockVector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyToBlockVector( std::vector<Teuchos::RCP<N_LAS_Vector> >& inputVectors, N_LAS_BlockVector& blockVector ); 

//-----------------------------------------------------------------------------
// Function      : copyFromBlockVector
// Purpose       : A helper function that copies a N_LAS_BlockVector to an 
//               : array of N_LAS_Vectors.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyFromBlockVector( N_LAS_BlockVector& blockVector, std::vector<Teuchos::RCP<N_LAS_Vector> >& outputVectors ); 

#endif
