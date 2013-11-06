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
// Revision Number: $Revision: 1.14.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:51 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktNodeCreator_h
#define N_TOP_CktNodeCreator_h 1

// ---------- Standard Includes ----------
#include <list>
#include <string>
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_DEV_fwd.h>

// ---------- Forward Declarations ----------

class N_TOP_CktNode;
class N_TOP_CktGraph;
class N_TOP_NodeBlock;

//-----------------------------------------------------------------------------
// Class         : N_TOP_CktNodeCreator
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class N_TOP_CktNodeCreator
{

private:

  // Constructor (private).
  N_TOP_CktNodeCreator() { }

public:
  // Singleton
  static N_TOP_CktNodeCreator * instance();

  // Destructor
  ~N_TOP_CktNodeCreator() { }

private:

  // Copy constructor (private).
  N_TOP_CktNodeCreator(const N_TOP_CktNodeCreator & right);

  // Assignment operator (private).
  N_TOP_CktNodeCreator & operator = (const N_TOP_CktNodeCreator & right);

  bool operator == (const N_TOP_CktNodeCreator & right) const;
  bool operator != (const N_TOP_CktNodeCreator & right) const;

public:

  //------- Voltage Node Creators
  N_TOP_CktNode * CreateVoltageNode(const string & ID) const;
  N_TOP_CktNode * CreateVoltageNode(const tagged_param & tpID) const;
  N_TOP_CktNode * CreateVoltageNode(const N_TOP_NodeBlock & nb) const;

  //------- Device Node Creators
  N_TOP_CktNode * CreateDeviceNode(const string & ID,
                                   N_DEV_DeviceInstance * instancePtr) const;
  N_TOP_CktNode * CreateDeviceNode(const string & ID,
                                   N_DEV_DeviceInstance * instancePtr, int gID,
                                   list < int > varGIDList) const;
  N_TOP_CktNode * CreateDeviceNode(const N_TOP_NodeBlock & nb,
                                   N_DEV_DeviceInstance * instancePtr) const;
  N_TOP_CktNode * CreateDeviceNode( const N_TOP_NodeBlock & nb,
                                    const Teuchos::RCP<N_DEV_InstanceBlock> ibPtr,
		                                N_DEV_DeviceInterface & devIF ) const;

  //------- Sub Ckt Node Creators
  N_TOP_CktNode * CreateSubCktNode(const string & ID,
                                   N_TOP_CktGraph * grphPtr) const;

  N_TOP_CktNode * CreateSubCktNode(const N_TOP_NodeBlock & nb,
                                   N_TOP_CktGraph * grphPtr) const;

};

#endif
