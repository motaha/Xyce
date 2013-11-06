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
// Revision Number: $Revision: 1.7.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  _N_TOP_MISC_H
#define  _N_TOP_MISC_H

// ---------- Standard Includes ----------

#include <string>
#include <iostream>

#ifdef NEED_PAIR_H
#include <pair.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

// ----------   Fwd Declares   -----------

// ----------   DEFINES  -----------------

// ----------   TYPEDEFS -----------------

// ----------   ENUMS    -----------------

enum NodeTYPE {_VNODE, _DNODE, _CNODE, _PNODE, _NUM_NODE_TYPES};

// ----------   STRUCTS  -----------------

class NodeID : public std::pair< string, int >
{
  public:
 
  // Default constructor
  NodeID()
    :std::pair<string, int>()
  {}

  // Basic constructor
  NodeID(const string& node, int id)
    :std::pair<string, int>( node, id )
  {}
};
inline std::ostream& operator<< (std::ostream &os, const NodeID& n)
{
  return os << "( " << n.first << " , " << n.second << " )";
} 



#endif
