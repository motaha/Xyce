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
// Filename       : $RCSfile: N_TOP_Directory.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_Directory_h
#define N_TOP_Directory_h 1

// ---------- Standard Includes ----------
#include <string>
#include <vector>
#include <iosfwd>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_TOP_Misc.h>

// ---------- Forward Declarations ----------

class N_TOP_Node;
class N_TOP_Topology;

class N_PDS_Manager;

class N_TOP_DirectoryData;

//-----------------------------------------------------------------------------
// Class         : N_TOP_Directory
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/05/01
//-----------------------------------------------------------------------------
class N_TOP_Directory
{

public:

  // Constructor
  N_TOP_Directory( N_TOP_Topology * topo,
                   N_PDS_Manager * pds = 0 )
  : topMgr_(topo),
    pdsMgr_(pds),
    data_(0)
  { }

  // Destructor
  ~N_TOP_Directory();

  // Registers the PDS manager.
  bool registerParallelServices(N_PDS_Manager * pds)
  { return (pdsMgr_ = pds) != 0; }

  bool generateDirectory();

  bool getGIDs( const vector<NodeID> & idVec,
                vector<int> & gidVec );
  bool getProcs( const vector<NodeID> & idVec,
                 vector<int> & procVec );

  bool getSolnGIDs( const vector<NodeID> & idVec,
                    vector< vector<int> > & gidVec,
                    vector<int> & procVec );
  bool getStateGIDs( const vector<NodeID> & idVec,
                     vector< vector<int> > & gidVec,
                     vector<int> & procVec );
  bool getStoreGIDs( const vector<NodeID> & idVec,
                     vector< vector<int> > & gidVec,
                     vector<int> & procVec );

private:

  // Copy constructor (private).
  N_TOP_Directory(const N_TOP_Directory & right);

  // Assignment operator (private).
  N_TOP_Directory & operator=(const N_TOP_Directory & right);

  // Equality operator.
  bool operator==(const N_TOP_Directory & right) const;

  // Non-equality operator.
  bool operator!=(const N_TOP_Directory & right) const;

private:

  // Pointer to the topology manager.
  N_TOP_Topology * topMgr_;

  // Pointer to the PDS manager.
  N_PDS_Manager * pdsMgr_;

  // Pimpl
  N_TOP_DirectoryData * data_;

};

#endif
