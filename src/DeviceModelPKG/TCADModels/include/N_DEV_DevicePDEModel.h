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
// Filename       : $RCSfile: N_DEV_DevicePDEModel.h,v $
//
// Purpose        : This file contains the PDE device model base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DevicePDEModel_h
#define Xyce_N_DEV_DevicePDEModel_h

// ---------- Standard Includes ----------

// ------------- Xyce Includes ------------
#include <N_DEV_DeviceModel.h>

// ---------- Forward Declarations ----------
namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DevicePDEModel
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
class DevicePDEModel: public DeviceModel
{
public:

  DevicePDEModel(
    const ModelBlock & MB,
    SolverState &      ss1,
    DeviceOptions &    do1)
    : DeviceModel (MB, ss1, do1),
      dopeInfoMap(),
      regionMap(),
      nodeMap()
  {}

  virtual ~DevicePDEModel () {};

private:
  DevicePDEModel(const DevicePDEModel &);

public:
  map<string, DopeInfo *> dopeInfoMap;
  map<string, CompositeParam *> regionMap;
  map<string, CompositeParam *> nodeMap;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DevicePDEModel N_DEV_DevicePDEModel;

#endif

