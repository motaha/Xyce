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
// Filename       : $RCSfile: N_TOP_DevInsertionTool.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/11/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <string>

// ---------- Xyce Includes --------------

#include <N_TOP_DevInsertionTool.h>
#include <N_TOP_Topology.h>

//#include <N_TOP_System.h>
//#include <N_TOP_Node.h>

#include <N_TOP_NodeBlock.h>

#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace Topology {

//-----------------------------------------------------------------------------
// Function      : Xyce::Topology::DevInsertionTool::insertNode
// Purpose       : Insert Node in System
// Special Notes : 
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/11/03
//-----------------------------------------------------------------------------
void DevInsertionTool::insertNode( const N_TOP_NodeBlock & nBlock,
                                   Teuchos::RefCountPtr<N_DEV_InstanceBlock> iBlockPtr )
{
/*
  Node * new_node = sys_->NodeFactory().createNode( nBlock );

  addGraphNode( new_node );

  if( nBlock.type == "D" ) //device
  {
    list<tagged_param>::const_iterator iterNL = nb.get_NodeList().begin();
    list<tagged_param>::const_iterator endNL = nb.get_NodeList().end();
    for( ; iterNL != endNL; ++iterNL )
    {
      Node * vnode = sys_->NodeFactory().createNode( N_TOP_NodeBlock( iterNL->tag ) );

      addGraphNode( vnode );
      addGraphEdge( new_node, vnode ):
    }
  }
  else if( nBlock.type = "V" )
  {
    //no_op
  }
  else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
    "Topology::DevInsertionTool::insertNode - unrecognized node type\n" );
*/
  topo_.addDevice( nBlock, iBlockPtr );
//  delete iBlock;
}

} //namespace Topology
} //namespace Xyce
