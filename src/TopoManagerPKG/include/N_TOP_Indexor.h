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
// Filename       : $RCSfile: N_TOP_Indexor.h,v $
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
// Revision Number: $Revision: 1.8.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_Indexor_h
#define N_TOP_Indexor_h 1

// ---------- Standard Includes ----------

#include <vector>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

#include <N_TOP_Misc.h>

// ---------- Forward Declarations ----------

class N_PDS_Manager;

//-----------------------------------------------------------------------------
// Class         : N_TOP_Indexor
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
class N_TOP_Indexor
{

public:

  // Constructor
  N_TOP_Indexor( N_PDS_Manager * pds = 0 )
  : pdsMgr_(pds),
    accelMatrixIndex_(false)
  { }

  // Destructor
  ~N_TOP_Indexor() {}

  // Registers the PDS manager.
  bool registerParallelServices(N_PDS_Manager * pds)
  { return (pdsMgr_ = pds); }

  bool globalToLocal( const string & map_name, vector<int> & ids );
  bool localToGlobal( const string & map_name, vector<int> & ids );

  bool setupAcceleratedMatrixIndexing( const string & graph_name );
  bool deleteAcceleratedMatrixIndexing();

  bool matrixGlobalToLocal( const string & graph_name,
                            const vector<int> & gids,
                            vector< vector<int> > & stamp );

private:

  // Copy constructor (private).
  N_TOP_Indexor(const N_TOP_Indexor & right);
  // Assignment operator (private).
  N_TOP_Indexor & operator=(const N_TOP_Indexor & right);
  // Equality operator.
  bool operator==(const N_TOP_Indexor & right) const;
  // Non-equality operator.
  bool operator!=(const N_TOP_Indexor & right) const;

private:

  // Pointer to the PDS manager.
  N_PDS_Manager * pdsMgr_;

  bool accelMatrixIndex_;
  vector< map<int,int> > matrixIndexMap_;

};

#endif
