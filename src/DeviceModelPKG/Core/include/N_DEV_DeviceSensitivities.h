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
// Filename       : $RCSfile: N_DEV_DeviceSensitivities.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/15/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.16.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceSensitivities_h
#define Xyce_N_DEV_DeviceSensitivities_h

// ---------- Standard Includes ----------

#include <vector>
#include <string>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceEntity.h>

// ---------- Forward Declarations ----------
class N_NLS_Manager;

class N_LAS_System;
class N_LAS_Matrix;
class N_LAS_MultiVector;
class N_LAS_Vector;

class N_ERH_ErrorMgr;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceSensitivities
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------

class DeviceSensitivities
{
  // functions:
public:
  DeviceSensitivities ( DeviceOptions & do1,
                              ExternData & extData1,
                              SolverState & solState1,
                              N_LAS_System & lasSys1);

  ~DeviceSensitivities();

  bool registerSensParams (const N_UTL_OptionBlock & OB);

  DeviceEntity* getDeviceEntity(const string & param);

  bool stripParamName
  (const std::string & fullName,
   string & strippedName,
   string & entityName) const;


private:
  DeviceSensitivities (const DeviceSensitivities & right);

  void setUpDeviceEntityMap_ ();

  const DeviceEntity * findDeviceEntity_ (const std::string & param) const;
  DeviceEntity * findDeviceEntity_ (const std::string & param);

public:

private:
  map <string, DeviceEntity*> nameDevEntityMap_;
  std::vector<std::string> sensParamArray_;

  bool entityMapDone_;
  int numSensParams_;

  DeviceMgr     & devMgr_;
  DeviceOptions & devOptions_;
  ExternData    & extData_;
  SolverState   & solState_;
  N_LAS_System        & lasSys_;
  const N_NLS_Manager & nlsMgr_;

  // Petra options block (passed to iterative solver)
  N_UTL_OptionBlock* petraOptionBlockPtr_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceSensitivities N_DEV_DeviceSensitivities;

#endif

