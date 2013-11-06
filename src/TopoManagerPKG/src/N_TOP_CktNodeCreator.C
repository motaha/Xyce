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
// Revision Number: $Revision: 1.13.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_TOP_CktNodeCreator.h>
#include <N_TOP_CktNode.h>
#include <N_TOP_CktNode_V.h>
#include <N_TOP_CktNode_Ckt.h>
#include <N_TOP_CktNode_Dev.h>

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::N_TOP_CktNodeCreator
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNodeCreator::N_TOP_CktNodeCreator(const N_TOP_CktNodeCreator &right)
{
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::instance
// Purpose       :
// Special Notes : singleton
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/27/01
//-----------------------------------------------------------------------------
N_TOP_CktNodeCreator * N_TOP_CktNodeCreator::instance()
{
  N_TOP_CktNodeCreator * cPtr = new N_TOP_CktNodeCreator ();
  return cPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::operator=
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNodeCreator & N_TOP_CktNodeCreator::operator=(const N_TOP_CktNodeCreator &right)
{
  return (*this);
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::operator==
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
bool N_TOP_CktNodeCreator::operator==(const N_TOP_CktNodeCreator &right) const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::operator!=
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
bool N_TOP_CktNodeCreator::operator!=(const N_TOP_CktNodeCreator &right) const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateVoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateVoltageNode( const
		string & ID ) const
{
  return new N_TOP_CktNode_V(ID);
}


//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateVoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateVoltageNode( const
		tagged_param & tpID ) const
{
  return new N_TOP_CktNode_V(tpID.tag,static_cast<const int>(tpID.param));
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateVoltageNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateVoltageNode( const
		N_TOP_NodeBlock & nb ) const
{
  return new N_TOP_CktNode_V( nb );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateDeviceNode( const
		string & ID, N_DEV_DeviceInstance *
		instancePtr ) const
{
  return new N_TOP_CktNode_Dev(instancePtr, ID);
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateDeviceNode( const
		string & ID, N_DEV_DeviceInstance *
		instancePtr, int gID, list<int> varGIDList ) const
{
  return new N_TOP_CktNode_Dev(instancePtr, ID, gID, varGIDList);
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateDeviceNode( const
		N_TOP_NodeBlock & nb, N_DEV_DeviceInstance *
		instancePtr ) const
{
  return new N_TOP_CktNode_Dev(instancePtr, nb);
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateDeviceNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/02/03
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateDeviceNode(
	       	const N_TOP_NodeBlock & nb,
		      const Teuchos::RCP<N_DEV_InstanceBlock> ibPtr,
		      N_DEV_DeviceInterface & devIF ) const
{
  return new N_TOP_CktNode_Dev( nb, ibPtr, devIF );
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateSubCktNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateSubCktNode( const
		string & ID, N_TOP_CktGraph*
		grphPtr ) const
{
  return new N_TOP_CktNode_Ckt(grphPtr, ID);
}

//-----------------------------------------------------------------------------
// Function      : N_TOP_CktNodeCreator::CreateSubCktNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
N_TOP_CktNode* N_TOP_CktNodeCreator::CreateSubCktNode( const
		N_TOP_NodeBlock & nb, N_TOP_CktGraph*
		grphPtr ) const
{
  return new N_TOP_CktNode_Ckt(grphPtr, nb);
}

