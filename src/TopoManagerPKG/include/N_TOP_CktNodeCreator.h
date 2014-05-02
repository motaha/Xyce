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
// Filename       : $RCSfile: N_TOP_CktNodeCreator.h,v $
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
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNodeCreator_h
#define N_TOP_CktNodeCreator_h 1

#include <list>
#include <string>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_DEV_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_Misc.h>

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : CktNodeCreator
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class CktNodeCreator
{

private:

  // Constructor (private).
  CktNodeCreator() { }

public:
  // Singleton
  static CktNodeCreator * instance();

  // Destructor
  ~CktNodeCreator() { }

private:

  // Copy constructor (private).
  CktNodeCreator(const CktNodeCreator & right);

  // Assignment operator (private).
  CktNodeCreator & operator = (const CktNodeCreator & right);

  bool operator == (const CktNodeCreator & right) const;
  bool operator != (const CktNodeCreator & right) const;

public:

  //------- Voltage Node Creators
  CktNode * CreateVoltageNode(const std::string & ID) const;
  CktNode * CreateVoltageNode(const tagged_param & tpID) const;
  CktNode * CreateVoltageNode(const NodeBlock & nb) const;

  //------- Device Node Creators
  CktNode * CreateDeviceNode(const std::string & ID,
                                   Device::DeviceInstance * instancePtr) const;
  CktNode * CreateDeviceNode(const std::string & ID,
                                   Device::DeviceInstance * instancePtr, int gID,
                                   std::list< int > varGIDList) const;
  CktNode * CreateDeviceNode(const NodeBlock & nb,
                                   Device::DeviceInstance * instancePtr) const;
  CktNode * CreateDeviceNode( const NodeBlock & nb,
                                    const Teuchos::RCP<Device::InstanceBlock> ibPtr,
		                                Device::DeviceInterface & devIF ) const;

  //------- Sub Ckt Node Creators
  CktNode * CreateSubCktNode(const std::string & ID,
                                   CktGraph * grphPtr) const;

  CktNode * CreateSubCktNode(const NodeBlock & nb,
                                   CktGraph * grphPtr) const;

};

typedef Xyce::Topo::CktNodeCreator N_TOP_CktNodeCreator;

} // namespace Topo
} // namespace Xyce

#endif
