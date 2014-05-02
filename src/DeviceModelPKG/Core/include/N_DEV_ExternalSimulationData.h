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
// Filename       : $RCSfile$
//
// Purpose        : This is a container class for solver information that is
//                  used during coupled simulation runs.  The outer code
//                  will populate this struct and pass it into a xyce
//                  inner solve using the simulateStep method on N_CIR_Xyce.
//
//
// Special Notes  :
//
// Creator        : Roger P. Pawlowski, SNL, Applied Math and Applications
//
// Creation Date  : 07/16/2009
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_ExternalSimulationData_h
#define Xyce_N_DEV_ExternalSimulationData_h

#include <map>
#include <string>
#include <vector>

namespace Xyce {
namespace Device {

struct ExternalSimulationData
{
  bool is_transient;
  double current_time;
  double final_time;
  double current_time_step_size;
  double previous_time_step_size;
  int time_step_number;

  // From N_TIA_TimeIntInfo - copying all for now -- RWH
  int    currentOrder;
  int    numberOfSteps;
  int    usedOrder;
  int    nscsco;

  double pdt;
  double nextTimeStep;
  double currTimeStep;
  double currentTime;
  double nextTime;
  double finalTime;
  double startingTimeStep;
  double bpTol;
  bool   dcopFlag;
  bool   acopFlag;
  bool   inputOPFlag;
  bool   tranopFlag;
  bool   transientFlag;
  bool   dcsweepFlag;
  int    timeStepNumber;
  bool   initTranFlag;
  bool   beginIntegrationFlag;
  int    doubleDCOPStep;
  bool   doubleDCOPEnabled;
  int    stepLoopIter;
  int    timeIntMode;
  bool   sweepSourceResetFlag;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ExternalSimulationData N_DEV_ExternalSimulationData;

#endif

