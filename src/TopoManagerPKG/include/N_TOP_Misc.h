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
// Filename       : $RCSfile: N_TOP_Misc.h,v $
//
// Purpose        : Misc. typedefs, defines, and structs for Topology
//
// Special Notes  :
//
// Creator        : Rob Hoekstra
//
// Creation Date  : 7/03/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  _N_TOP_MISC_H
#define  _N_TOP_MISC_H

#include <string>
#include <iostream>
#include <utility>

#include <N_UTL_Xyce.h>

namespace Xyce {

enum NodeTYPE {_VNODE, _DNODE, _CNODE, _PNODE, _NUM_NODE_TYPES};


class NodeID : public std::pair< std::string, int >
{
  public:
 
  // Default constructor
  NodeID()
    :std::pair<std::string, int>()
  {}

  // Basic constructor
  NodeID(const std::string& node, int id)
    :std::pair<std::string, int>( node, id )
  {}
};
inline std::ostream& operator<< (std::ostream &os, const NodeID& n)
{
  return os << "( " << n.first << " , " << n.second << " )";
} 

} // namespace Xyce

#endif
