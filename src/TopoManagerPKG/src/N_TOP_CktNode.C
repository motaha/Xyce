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
// Revision Number: $Revision: 1.13.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_TOP_CktNode.h>
#include <N_TOP_NodeBlock.h>

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNode::N_TOP_CktNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode::N_TOP_CktNode( const N_TOP_NodeBlock & nb,
				N_TOP_GraphNode * gN )
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
// Function      : N_TOP_CktNode::extractNodeBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/4/01
//-----------------------------------------------------------------------------
N_TOP_NodeBlock * N_TOP_CktNode::extractNodeBlock()
{
  list<tagged_param> emptyList;

  return new N_TOP_NodeBlock( id_, emptyList, emptyList, isOwned_,
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
ostream& operator<< (ostream& os, const N_TOP_CktNode& cn)
{
  return cn.put(os);
}

