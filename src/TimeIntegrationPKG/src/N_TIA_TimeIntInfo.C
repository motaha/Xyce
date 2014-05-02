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
// Filename      : $RCSfile: N_TIA_TimeIntInfo.C,v $
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/29/07
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <iostream>

#include <N_TIA_TimeIntInfo.h>
#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for two-level info class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const N_TIA_TimeIntInfo & tiInfo)
{
  os << Xyce::section_divider << std::endl
     << "\tTimeIntInfo:\n"
     << "\tcurrentOrder = " << tiInfo.currentOrder << "\n"
     << "\tnumberOfSteps = " << tiInfo.numberOfSteps << "\n"
     << "\tusedOrder = " << tiInfo.usedOrder << "\n"
     << "\tnscsco = " << tiInfo.nscsco << "\n"
     << "\tpdt = " << tiInfo.pdt << "\n"
     << "\tnextTimeStep = " << tiInfo.nextTimeStep << "\n"
     << "\tcurrTimeStep = " << tiInfo.currTimeStep << "\n"
     << "\tnextTime = " << tiInfo.nextTime << "\n"
     << "\tcurrentTime = " << tiInfo.currentTime << "\n"
     << "\tfinalTime = " << tiInfo.finalTime << "\n"
     << "\tstartingTimeStep = " << tiInfo.startingTimeStep << "\n"
     << "\tbpTol = " << tiInfo.bpTol << "\n"
     << "\tdcopFlag = " << tiInfo.dcopFlag << "\n"
     << "\tacopFlag = " << tiInfo.acopFlag << "\n"
     << "\tinputOPFlag = " << tiInfo.inputOPFlag << "\n"
     << "\ttranopFlag = " << tiInfo.tranopFlag << "\n"
     << "\ttransientFlag = " << tiInfo.transientFlag << "\n"
     << "\tdcsweepFlag = " << tiInfo.dcsweepFlag << "\n"
     << "\tsweepSourceResetFlag = " << tiInfo.sweepSourceResetFlag << "\n"
     << "\ttimeStepNumber = " << tiInfo.timeStepNumber << "\n"
     << "\tinitTranFlag = " << tiInfo.initTranFlag << "\n"
     << "\tbeginIntegrationFlag = " << tiInfo.beginIntegrationFlag << "\n"
     << "\tdoubleDCOPStep = " << tiInfo.doubleDCOPStep << "\n"
     << "\tdoubleDCOPEnabled = " << tiInfo.doubleDCOPEnabled << "\n"
     << "\tstepLoopIter = " << tiInfo.stepLoopIter << "\n"
     << Xyce::section_divider << std::endl
     << std::endl;

  return os;
}
