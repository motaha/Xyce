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
// Filename       : $RCSfile: N_TOP_CktNode.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.18 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_CktNode.h>
#include <N_TOP_NodeBlock.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : CktNode::CktNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode::CktNode( const NodeBlock & nb,
				GraphNode * gN )
{
  id_ = nb.get_id();
  gID_ = nb.get_gID();
  procNum_ = nb.get_ProcNum();
  isOwned_ = nb.get_IsOwned();
  graphNodePtr_ = gN;
  Offset_ = -1;

  solnVarGIDList_ = nb.get_SolnVarGIDList();
  stateVarGIDList_ = nb.get_StateVarGIDList();
  storeVarGIDList_ = nb.get_StoreVarGIDList();
}

//-----------------------------------------------------------------------------
// Function      : CktNode::extractNodeBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/4/01
//-----------------------------------------------------------------------------
NodeBlock * CktNode::extractNodeBlock()
{
  std::list<tagged_param> emptyList;

  return new NodeBlock( id_, emptyList, emptyList, isOwned_,
	gID_, procNum_, solnVarGIDList_, 
  stateVarGIDList_ , storeVarGIDList_ );
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
std::ostream& operator<< (std::ostream& os, const CktNode& cn)
{
  return cn.put(os);
}

} // namespace Topo
} // namespace Xyce
