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
// Revision Number: $Revision: 1.15 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_Directory_h
#define N_TOP_Directory_h 1

#include <string>
#include <vector>
#include <iosfwd>

#include <N_UTL_Xyce.h>
#include <N_TOP_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_Misc.h>

namespace Xyce {
namespace Topo {

class DirectoryData;

//-----------------------------------------------------------------------------
// Class         : Directory
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/05/01
//-----------------------------------------------------------------------------
class Directory
{

public:

  // Constructor
  Directory( Topology * topo, N_PDS_Manager * pds = 0 )
  : topMgr_(topo),
    pdsMgr_(pds),
    data_(0)
  { }

  // Destructor
  ~Directory();

  // Registers the PDS manager.
  bool registerParallelServices(N_PDS_Manager * pds)
  { return (pdsMgr_ = pds) != 0; }

  bool generateDirectory();

  bool getGIDs( const std::vector<NodeID> & idVec,
                std::vector<int> & gidVec );
  bool getProcs( const std::vector<NodeID> & idVec,
                 std::vector<int> & procVec );

  bool getSolnGIDs( const std::vector<NodeID> & idVec,
                    std::vector< std::vector<int> > & gidVec,
                    std::vector<int> & procVec );
  bool getStateGIDs( const std::vector<NodeID> & idVec,
                     std::vector< std::vector<int> > & gidVec,
                     std::vector<int> & procVec );
  bool getStoreGIDs( const std::vector<NodeID> & idVec,
                     std::vector< std::vector<int> > & gidVec,
                     std::vector<int> & procVec );

private:

  // Copy constructor (private).
  Directory(const Directory & right);

  // Assignment operator (private).
  Directory & operator=(const Directory & right);

  // Equality operator.
  bool operator==(const Directory & right) const;

  // Non-equality operator.
  bool operator!=(const Directory & right) const;

private:

  // Pointer to the topology manager.
  Topology * topMgr_;

  // Pointer to the PDS manager.
  N_PDS_Manager * pdsMgr_;

  // Pimpl
  DirectoryData * data_;
};

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::Directory N_TOP_Directory;

#endif
