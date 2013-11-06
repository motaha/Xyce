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
// Filename       : $RCSfile: N_DEV_ExternDevice_Factory.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/09/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ExternDevice_Factory_h
#define Xyce_N_DEV_ExternDevice_Factory_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>


//-----------------------------------------------------------------------------
// Class         : N_DEV_ExternDevice
// Purpose       : Handles external devices.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
struct N_DEV_ExternDevice
{
  //-----------------------------------------------------------------------------
  // Function      : factory
  // Purpose       : This is the factory function for the external devics
  //                 class.
  //
  // Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
  //                 static pointer was returned) but had to be changed
  //                 so that the library version of Xyce would work
  //                 correctly.
  //
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 3/16/00
  //-----------------------------------------------------------------------------
  static N_DEV_Device * factory(N_DEV_SolverState & ss1,
                                N_DEV_DeviceOptions & do1)
  {
     const string name("External Device");
     const string className("N_DEV_ExternDevice");
     const string defaultModelName ("YEXT level 1 (External Device)");

     N_DEV_InstanceBlock IB(defaultModelName);
     // Set up a "default" instance block:
     N_DEV_Param p ("NX", 11, true);
     p.set("NODE", "DEFAULT");
     IB.params.push_back(p);
     p.set("NODE0.NAME", "VCONNECT0000");
     IB.params.push_back(p);
     p.set("NODE1.NAME", "VCONNECT0001");
     IB.params.push_back(p);

     N_DEV_Device * devptr =
       new Xyce::Device::DeviceTemplate<N_DEV_ExternDeviceModel,
                                N_DEV_ExternDeviceInstance>
                   ( name,
                     className,
                     defaultModelName,
                     IB,
                     false,  // not linear
                     false,   // requires a model
                     ss1,do1);

     return devptr;
  }
};

#endif

