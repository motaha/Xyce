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
// Filename       : $RCSfile: N_TOP_Node.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/11/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_fwd.h>
#include <N_TOP_Node.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : Node::pack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
void Node::pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const
{
  int length = nodeID_.first.length();
  comm->pack( &length, 1, buf, bsize, pos );

  comm->pack( nodeID_.first.c_str(), length, buf, bsize, pos );

  comm->pack( &(nodeID_.second), 1, buf, bsize, pos );

  length = (owned_?1:0);
  comm->pack( &length, 1, buf, bsize, pos );
}

//-----------------------------------------------------------------------------
// Function      : Node::unpack
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
void Node::unpack( char * buf, int bsize, int & pos, N_PDS_Comm * comm )
{
  int length;
  comm->unpack( buf, bsize, pos, &length, 1 );

  nodeID_.first = std::string( (buf+pos), length );
  pos += length;

  comm->unpack( buf, bsize, pos, &(nodeID_.second), 1 );

  comm->unpack( buf, bsize, pos, &length, 1 );
  owned_ = (length!=0);
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & os, Node & node )
{
  return node.put(os);
}


//-----------------------------------------------------------------------------
// Function      : Node::put
// Purpose       :
// Special Notes :
// Scope         : protected
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
std::ostream & Node::put(std::ostream & os) const
{
  os << "NodeID:\t" << nodeID_.first << "\t" << nodeID_.second;
  if (owned_) os << "\tOWNED";
  return os << std::endl;
}

} // namespace Topo
} // namespace Xyce
