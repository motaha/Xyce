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
// Filename       : $RCSfile: N_TOP_CktNodeCreator.C,v $
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

#include <N_TOP_CktNodeCreator.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_V.h>
#include <N_TOP_CktNode_Ckt.h>
#include <N_TOP_CktNode_Dev.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CktNodeCreator
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNodeCreator::CktNodeCreator(const CktNodeCreator &right)
{
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::instance
// Purpose       :
// Special Notes : singleton
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/27/01
//-----------------------------------------------------------------------------
CktNodeCreator * CktNodeCreator::instance()
{
  CktNodeCreator * cPtr = new CktNodeCreator ();
  return cPtr;
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::operator=
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNodeCreator & CktNodeCreator::operator=(const CktNodeCreator &right)
{
  return (*this);
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::operator==
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
bool CktNodeCreator::operator==(const CktNodeCreator &right) const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::operator!=
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
bool CktNodeCreator::operator!=(const CktNodeCreator &right) const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateVoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateVoltageNode( const
		std::string & ID ) const
{
  return new CktNode_V(ID);
}


//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateVoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateVoltageNode( const
		tagged_param & tpID ) const
{
  return new CktNode_V(tpID.tag,static_cast<const int>(tpID.param));
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateVoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateVoltageNode( const
		NodeBlock & nb ) const
{
  return new CktNode_V( nb );
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateDeviceNode( const
		std::string & ID, Device::DeviceInstance *
		instancePtr ) const
{
  return new CktNode_Dev(instancePtr, ID);
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateDeviceNode( const
		std::string & ID, Device::DeviceInstance *
		instancePtr, int gID, std::list<int> varGIDList ) const
{
  return new CktNode_Dev(instancePtr, ID, gID, varGIDList);
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateDeviceNode( const
		NodeBlock & nb, Device::DeviceInstance *
		instancePtr ) const
{
  return new CktNode_Dev(instancePtr, nb);
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/02/03
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateDeviceNode(
	       	const NodeBlock & nb,
		      const Teuchos::RCP<Device::InstanceBlock> ibPtr,
		      Device::DeviceInterface & devIF ) const
{
  return new CktNode_Dev( nb, ibPtr, devIF );
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateSubCktNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateSubCktNode( const
		std::string & ID, CktGraph*
		grphPtr ) const
{
  return new CktNode_Ckt(grphPtr, ID);
}

//-----------------------------------------------------------------------------
// Function      : CktNodeCreator::CreateSubCktNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CktNode* CktNodeCreator::CreateSubCktNode( const
		NodeBlock & nb, CktGraph*
		grphPtr ) const
{
  return new CktNode_Ckt(grphPtr, nb);
}

} // namespace Topo
} // namespace Xyce
