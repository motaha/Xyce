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
// Filename       : $RCSfile: N_PDS_ParMap.C,v $
//
// Purpose        : Implementation file for abstract base class for the
//                  parallel map data and functions.
//
// Special Notes  : Part of a GoF Abstract Factory.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.27.4.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_CommFactory.h>

#include <N_UTL_Misc.h>

// ----------   Other Includes   ----------

#include <Epetra_Map.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::N_PDS_ParMap
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/2/00
//-----------------------------------------------------------------------------
N_PDS_ParMap::N_PDS_ParMap(int & numGlobalEntities,
                           const int & numLocalEntities,
			   const vector<int> & lbMap,
                           const int index_base,
                           N_PDS_Comm * aComm)
  : petraMap_(0),
    mapOwned_(true),
    comm_(aComm),
    commOwned_(false)
{

#ifdef Xyce_DEBUG_PARALLEL
  static const string msg("N_PDS_ParMap::N_PDS_ParMap - ");
#endif

  if( comm_ == 0 )
  {
    comm_ = N_PDS_CommFactory::create();
    commOwned_ = true;
  }

  int * mArray = lbMap.empty() ? 0 : (int*)(&(lbMap[0]));

  // fix for empty maps
  int nGE = Xycemax( 0, numGlobalEntities );
  int nLE = Xycemax( 0, numLocalEntities );
  // Call the Petra constructor for the true Petra map.
  petraMap_ = new Epetra_Map( nGE,
                              nLE,
                              mArray,
                              index_base,
                              *(comm_->petraComm()) );

#ifdef Xyce_DEBUG_PARALLEL
  cout << "New Petra Map: " << numGlobalEntities << " " <<
	numLocalEntities << endl;
  cout << "  " << petraMap_->NumMyElements() << " " <<
	petraMap_->NumGlobalElements() << endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::N_PDS_ParMap
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/2/00
//-----------------------------------------------------------------------------
N_PDS_ParMap::N_PDS_ParMap( Epetra_Map * map,
                            N_PDS_Comm * aComm )
  : petraMap_(map),
    mapOwned_(false),
    comm_(aComm),
    commOwned_(false)
{
#ifdef Xyce_DEBUG_PARALLEL
  static const string msg("N_PDS_ParMap::N_PDS_ParMap - ");
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::~N_PDS_ParMap
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
N_PDS_ParMap::~N_PDS_ParMap()
{
  if( commOwned_ )
    if( comm_ ) delete comm_;

  if( mapOwned_ )
    if( petraMap_ ) delete petraMap_;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::numGlobalEntities
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::numGlobalEntities() const
{
  return petraMap_->NumGlobalElements();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::numLocalEntities
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::numLocalEntities() const
{
  return petraMap_->NumMyElements();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::indexBase
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::indexBase() const
{
  return petraMap_->IndexBase();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::maxGlobalEntity
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParMap::maxGlobalEntity() const
{
  return petraMap_->MaxAllGID();
}

/*
//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::parMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
const int * N_PDS_ParMap::parMap() const
{
}
*/

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::petraBlockMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
Epetra_BlockMap * N_PDS_ParMap::petraBlockMap()
{
  return dynamic_cast<Epetra_BlockMap*>(petraMap_);
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::globalToLocalIndex
// Purpose       : dereference Global to Local Index
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/6/00
//-----------------------------------------------------------------------------
int N_PDS_ParMap::globalToLocalIndex(const int & global_index)
{
  return petraMap_->LID( global_index );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMap::localToGlobalIndex
// Purpose       : dereference Local to Global Index
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/10/06
//-----------------------------------------------------------------------------
int N_PDS_ParMap::localToGlobalIndex(const int & local_index)
{
  return petraMap_->GID( local_index );
}
