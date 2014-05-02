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
// Filename       : $RCSfile: N_TOP_Indexor.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/12/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <map>

#include <N_TOP_Indexor.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>

#include <Epetra_CrsGraph.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : Indexor::globalToLocal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/12/02
//-----------------------------------------------------------------------------
bool Indexor::globalToLocal( const std::string & map_name,
                                   std::vector<int> & ids )
{
  N_PDS_ParMap * map;

  assert( pdsMgr_ != 0 );
  // Never, EVER do work inside an assert argument, or that work will not
  // be done when asserts are disabled.
  map = pdsMgr_->getParallelMap( map_name );
  assert( map != 0 );

  map = pdsMgr_->getParallelMap( map_name );

  for( unsigned int i = 0; i < ids.size(); ++i ) ids[i] = map->globalToLocalIndex( ids[i] );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::localToGlobal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/28/10
//-----------------------------------------------------------------------------
bool Indexor::localToGlobal( const std::string & map_name,
                                   std::vector<int> & ids )
{
  N_PDS_ParMap * map;

  assert( pdsMgr_ != 0 );
  // Never, EVER do work inside an assert argument, or that work will not
  // be done when asserts are disabled.
  map = pdsMgr_->getParallelMap( map_name );
  assert( map != 0 );

  map = pdsMgr_->getParallelMap( map_name );

  for( unsigned int i = 0; i < ids.size(); ++i ) ids[i] = map->localToGlobalIndex( ids[i] );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::setupAcceleratedMatrixIndexing
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool Indexor::setupAcceleratedMatrixIndexing( const std::string & graph_name )
{
  Epetra_CrsGraph * graph = 0;

  assert( pdsMgr_ != 0 );
  // Never, EVER do work inside an assert argument, or that work will not
  // be done when asserts are disabled.
  graph = pdsMgr_->getMatrixGraph( graph_name );
  assert( graph != 0 );

  int NumRows = graph->NumMyRows();
  matrixIndexMap_.clear();
  matrixIndexMap_.resize( NumRows );

  int NumElements;
  int * Elements;
  for( int i = 0; i < NumRows; ++i )
  {
    graph->ExtractMyRowView( i, NumElements, Elements );
    for( int j = 0; j < NumElements; ++j ) matrixIndexMap_[i][ Elements[j] ] = j;
  }

  accelMatrixIndex_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::deleteAcceleratedMatrixIndexing
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool Indexor::deleteAcceleratedMatrixIndexing()
{
  matrixIndexMap_.clear();
  accelMatrixIndex_ = false;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Indexor::matrixGlobalToLocal
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/23/02
//-----------------------------------------------------------------------------
bool Indexor::matrixGlobalToLocal( const std::string & graph_name,
                                         const std::vector<int> & gids,
                                         std::vector< std::vector<int> > & stamp )
{
  Epetra_CrsGraph * graph = 0;

  assert( pdsMgr_ != 0 );
  // Never, EVER do work inside an assert argument, or that work will not
  // be done when asserts are disabled.
  graph = pdsMgr_->getMatrixGraph( graph_name );
  assert( graph != 0 );

  int numRows = stamp.size();

  int numElements;
  int * elements;

  if( accelMatrixIndex_ )
  {
    for( int i = 0; i < numRows; ++i )
    {
      int RowLID = graph->LRID(gids[i]);
      int NumCols = stamp[i].size();
      for( int j = 0; j < NumCols; ++j )
      {
        int lid = graph->LCID(stamp[i][j]);
        stamp[i][j] = matrixIndexMap_[RowLID][lid];
      }
    }
  }
  else
  {
    for( int i = 0; i < numRows; ++i )
    {
      graph->ExtractMyRowView( graph->LRID(gids[i]), numElements, elements );

      std::map<int,int> indexToOffsetMap;
      for( int j = 0; j < numElements; ++j ) indexToOffsetMap[ elements[j] ] = j;

      int numCols = stamp[i].size();
      for( int j = 0; j < numCols; ++j )
      {
        int lid = graph->LCID(stamp[i][j]);
//        assert( indexToOffsetMap.count(lid) );
        stamp[i][j] = indexToOffsetMap[lid];
      }
    }
  }

  return true;
}

} // namespace Topo
} // namespace Xyce
