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
// Filename       : $RCSfile: N_TOP_DevInsertionTool.h,v $
//
// Purpose        : Tool for insertion of nodes into a system
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/11/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_Topology_DevInsertionTool_h
#define Xyce_Topology_DevInsertionTool_h 1

// ----------   Std Includes  ------------
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_DEV_fwd.h>
#include <N_UTL_Misc.h>

#include <N_TOP_InsertionTool.h>

// ---------   Fwd Declares ---------------

class N_TOP_Topology;
class N_TOP_NodeBlock;

namespace Xyce {
namespace Topology {

//-----------------------------------------------------------------------------
// Class         : Xyce::Topology::DevInsertionTool
// Purpose       : Insertion of Nodes into a System
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/11/03
//-----------------------------------------------------------------------------
class DevInsertionTool : public InsertionTool
{

 public:

  // Constructor
//  DevInsertionTool( System & sys )
//  : InsertionTool(sys)
//  {}
	 //
  DevInsertionTool( N_TOP_Topology & topo )
  : topo_(topo)
  {}

  // Destructor
  ~DevInsertionTool() {}

  void insertNode( const N_TOP_NodeBlock & nBlock, Teuchos::RefCountPtr<N_DEV_InstanceBlock> iBlock );

 private:

  N_TOP_Topology & topo_;
};

} //namespace Topology
} //namespace Xyce

#endif //Xyce_Topology_DevInsertionTool_h
